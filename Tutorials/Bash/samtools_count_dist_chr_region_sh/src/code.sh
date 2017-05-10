#!/bin/bash
#
# This app performs a SAMtools count via the following implementation:
#   - Scatter gather approach is taken
#     - Scatter not in files, but in regions counted with samtools view
#   - Index file is an optional parameter.  If not provided:
#     - Bam will be indexed
#   - SAMtools is provided via execDepends.  Review dxapp.json for usage

main() {

  #
  # Debugging setup
  # ---------------
  # Echo input variable names for easier review in platform logs
  #

  set -e -x -o pipefail
  echo "Value of mappings_sorted_bam: '${mappings_sorted_bam}'"

  #
  # Index bam if no index file provided
  # -----------------------------------
  # Download only the bam, then check for an optional index file.
  # Create an index file if needed.
  # dx upload is used to generate dx link to pass to distributed jobs.
  #

  dx download "${mappings_sorted_bam}"

  if [ -z "${mappings_sorted_bai}" ]; then
    samtools index "${mappings_sorted_bam_name}"
    mappings_sorted_bai=$(dx upload "${mappings_sorted_bam_name}.bai" --brief)
  fi

  echo "Value of index file: ${mappings_sorted_bai}"

  #
  # SECTION: Download and prepare regions for scatter
  # ----------------------------------------
  # To scatter processing into regions of size 10
  #

  # Get region information from BAM headers
  regions=$(samtools view -H "${mappings_sorted_bam_name}" | grep "\@SQ" | sed 's/.*SN:\(\S*\)\s.*/\1/')

  # The following code will "pass" an array of regions in groups of 10 to the count job
  #
  # In bash arrays cannot be passed as a parameter.
  # However since an array approach is helpful in distributed and scatter-gather pattern
  # dx-jobutil-new-job allows an array to be "passed" as a parameter by:
  #  passing multiple -i<parameter> with the same name
  #  The receiving function will combine parameters with the same name to one array
  echo "Segmenting into regions"
  count_jobs=()
  counter=0
  temparray=()
  # Loop through regions, for every 10 create a new job with the correct parameters
  for r in $(echo $regions); do
    if [[ "${counter}" -ge 10 ]]; then
      echo "${temparray[@]}"
      count_jobs+=($(dx-jobutil-new-job -ibam_file="${mappings_sorted_bam}" -ibambai_file="${mappings_sorted_bai}" "${temparray[@]}" count_func))
      temparray=()
      counter=0
    fi
    temparray+=("-iregions=${r}") # Here we add to an array of -i<parameter>'s
    counter=$((counter+1))
  done

  if [[ counter -gt 0 ]]; then # Previous loop will miss last iteration  if its < 10
    echo "${temparray[@]}"
    count_jobs+=($(dx-jobutil-new-job -ibam_file="${mappings_sorted_bam}" -ibambai_file="${mappings_sorted_bai}" "${temparray[@]}" count_func))
  fi

  #
  # SECTION: Merge results
  # -------------
  # We take a similar approach and "pass" an array of readcount files to the merge job
  #   We add -ireadfiles="$count_job":counts_txt so the merge job can combine into an array called readfiles
  #
  # "$count_job":counts_txt" is used as a “job-based object reference” [JBOR]
  #   The merge sum_reads job will be able to run using the outputs of the "count_job"
  #

  echo "Merge count files, jobs:"
  echo "${count_jobs[@]}"
  readfiles=()
  for count_job in "${count_jobs[@]}"; do
    readfiles+=("-ireadfiles=${count_job}:counts_txt")
  done
  # Gather reads and merge results
  echo "file name: ${sorted_bamfile_name}"
  echo "Set file, readfile variables:"
  echo "${readfiles[@]}"
  countsfile_job=$(dx-jobutil-new-job -ifilename="${mappings_sorted_bam_prefix}" "${readfiles[@]}" sum_reads)

  #
  # SECTION: Output results
  # --------------
  # The output of our main() funciton is the result of the "sum_reads" job
  #
  # "sum_reads" job is referenced as the output of the main() funciton as
  # a job-based object reference [JBOR]
  #

  echo "Specifying output file"
  dx-jobutil-add-output counts_txt "${countsfile_job}:read_sum" --class=jobref
}

################################################################################
# SECTION: count_func
# Function that will perform SAMtools count by regions group.
# This function will be run on a new instance, as a result variables from other
# other functions, main() for example, will not be accessible here.
# 
# Arguments:
#   Any needed files must be passed as parameters and downloaded.
#   Parameters are passed as -i<parameter name>:
#     regions: an array of regions to be counted
#     bam_file: Sorted mapping file
#     bambai_file: Index file
#
# Returns:
#   outputs a file with the sum of counts in the region as the first line.
################################################################################
count_func() {

  set -e -x -o pipefail

  echo "Value of bam_file: '${bam_file}'"
  echo "Value of bambai_file: '${bambai_file}'"
  echo "Regions being counted '${regions[@]}'"

  #
  # Download all inputs
  # --------------------
  # dx-download-all-inputs will download all input files into
  # the /home/dnanexus/in directory.  A folder will be created for each
  # input and the file(s) will be download to that directory.
  # path, name, and prefix shell variables will be created for each input.
  #
  # In this example, the following variables are created
  #   bam_file_path: '/home/dnanexus/in/mappings_bam/<bam filename>.bam'
  #   mappings_bam_name: <bam filename>.bam
  #   mappings_bam_prefix: <bam filename>
  #   bambai_file_path: /home/dnanexus/in/index_file/<bambai_file name>
  #   bambai_file_name: <bambai_file name>
  #   bambai_file_prefix: bambai_file ** no extension **
  #   regions_path: '{region file path} {region file path} {region file path} ...'  # array of region paths
  #   regions_name: '{region file name} {region file name} {region file name} ...'  # array of region names
  #   regions_prefix: '{region file prefix} {region file prefix} ...'  # array of region paths
  #

  dx-download-all-inputs

  #
  # SAMtools count on the region group
  # ----------------------------------
  # Move bam and bai to same dir
  # Create output dir
  #

  mkdir workspace
  cd workspace || exit
  mv "${bam_file_path}" .
  mv "${bambai_file_path}" .
  # create output dir
  outputdir="./out/samtool/count"
  mkdir -p "${outputdir}"
  # count
  samtools view -c "${bam_file_name}" "${regions[@]}" >> "${outputdir}/readcounts.txt"

  #
  # Upload result file
  # ------------------
  # The output must be specified as counts_txt for other
  # functions to reference this job via a JBOR
  #

  counts_txt_id=$(dx upload "${outputdir}/readcounts.txt" --brief)
  dx-jobutil-add-output counts_txt "${counts_txt_id}" --class=file
}

#####################################################################
# SECTION: sum_reads
# This function will *gather* all the readcount.txt files generated
# by the count_func.
#
# Arguments:
#   readfiles: Array of count summary files from scattered jobs
#   filename: String filename of the input file.
#             Used to preserve filename in output.
#
# Returns:
#   read_sum: dxlink to output readcounts file to be referenced as a
#             Job-Based Object Reference (JBOR).
#
#####################################################################
sum_reads() {

  set -e -x -o pipefail
  echo "$filename"
  #
  # Download inputs
  # ---------------
  # Read count_func Download all inputs for an idea of what variables are available
  #

  echo "Value of read file array '${readfiles[@]}'"
  dx-download-all-inputs
  echo "Value of read file path array '${readfiles_path[@]}'"

  # Sum files
  echo "Summing values in files"
  readsum=0
  for read_f in "${readfiles_path[@]}"; do
    temp=$(cat "$read_f")
    readsum=$((readsum + temp))
  done

  # Write and upload files
  echo "Total reads: ${readsum}" > "${filename}_counts.txt"

  # Again naming the output read_sum is needed so the main() function's output can
  # reference read_sum as a “job-based object reference” [JBOR]
  read_sum_id=$(dx upload "${filename}_counts.txt" --brief)
  dx-jobutil-add-output read_sum "${read_sum_id}" --class=file
}
