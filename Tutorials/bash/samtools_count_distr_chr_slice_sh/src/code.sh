#!/bin/bash
#
# This app performs a SAMtools count via the following implementation:
#   - Index file is optional and will be created if missing
#   - Bam file will be split by chromosome (1-22, X, Y, and MT ) on a mem1_ssd1_x4 instance
#     - Files will be distributed to mem2_ssd1_x2 instances for SAMtools count processing
#     - Files will be merged in a mem1_ssd1_x4 instance
#     - Review dxapp.json to see how instance type declarations are made for entry points
#       - Specifically "runSpec" and "systemRequirements"
#   - Gather jobs will download result files via dx download -o - pipe
#   - SAMtools is added via execDepends
#

main() {

  #
  # Debugging setup
  # ---------------
  # The -e flag causes bash to exit at any point if there is any error,
  # the -o pipefail flag tells bash to throw an error if it encounters an error within a pipeline,
  # the -x flag causes bash to output each line as it is executed
  #

  set -e -x -o pipefail

  #
  # SECTION: Downloading inputs and get chromosomes list from headers
  # --------------------------------------------------------
  # Examine the header of the bam file, pulling out just the @SQ lines (which contain
  # the chromosomes, grab the column with just the chromosome name, and match it to a
  # regex designed to capture canonical chromosomes
  #

  dx download "${mappings_sorted_bam}"
  chromosomes=$(samtools view -H "${mappings_sorted_bam_name}" | grep "\@SQ" | awk -F '\t' '{print $2}' | awk -F ':' '{if ($2 ~ /^chr[0-9XYM]+$|^[0-9XYM]/) {print $2}}')

  #
  # SECTION: Split bam into multiple chromosomes and send as input to sam_count subjob
  # -------------------------------------------------------------------------
  # Looping through all our chromosomes of interest and seperate to smaller bams.
  #
  # Both .bam and .bai files are needed to split a bam by chromosome
  # First create an indexed file if one was not provided
  #

  if [ -z "${mappings_sorted_bai}" ]; then
    samtools index "${mappings_sorted_bam_name}"
  else
    dx download "${mappings_sorted_bai}" -o "${mappings_sorted_bam_name}".bai
  fi

  count_jobs=()
  for chr in $chromosomes; do
    seg_name="${mappings_sorted_bam_prefix}_${chr}".bam
    samtools view -b "${mappings_sorted_bam_name}" "${chr}" > "${seg_name}"
    bam_seg_file=$(dx upload "${seg_name}" --brief)
    count_jobs+=($(dx-jobutil-new-job -isegmentedbam_file="${bam_seg_file}" -ichr="${chr}" count_func))
  done

  #
  # SECTION: Gather output of count jobs and write to result file
  #------------------------------------------------------
  # Output of single count job is referenced as a job-based object reference
  #

  for job in "${count_jobs[@]}"; do
    readfiles+=("-ireadfiles=${job}:counts_txt")
  done

  sum_reads_job=$(dx-jobutil-new-job "${readfiles[@]}" -ifilename="${mappings_sorted_bam_prefix}" sum_reads)

  #
  # SECTION: Upload chromosome_results.txt from sum_reads subjob as job output
  # -------------------------------------------------
  # Make sure to specify --class=jobref
  #
  echo "Uploading file output"
  dx-jobutil-add-output counts_txt "${sum_reads_job}:read_sum_file" --class=jobref
}

####################################################################
# This function will count the number of reads in the input bam file
#
# Arguments:
#   chr: Specified chromosome for current job
#   segmentedbam_file: Mapping file for specified chromosome (chr)
####################################################################
count_func() {

  echo "Value of segmentedbam_file: '${segmentedbam_file}'"
  echo "Chromosome being counted '${chr}'"

  dx download "${segmentedbam_file}"

  readcount=$(samtools view -c "${segmentedbam_file_name}")

  #
  # Output the readcount to a file. <chromosome name>TAB<readcount result>
  # ----------------------------------------------------------------------
  # Upload readcount.txt as job out [readcount_file]
  # Other functions can reference this jobs output as a "job-based object reference” [JBOR]
  #

  >"${segmentedbam_file_prefix}.txt" printf "${chr}:\t%s\n" "${readcount}"

  readcount_file=$(dx upload "${segmentedbam_file_prefix}".txt --brief)
  dx-jobutil-add-output counts_txt "${readcount_file}" --class=file
}

####################################################################
# This function will *gather* all the readcount.txt files generated
# by the count_func.  It will return a single output file: chromosome_result.txt
#
# Arguments:
#   readfiles: Array of dxlinks to files with counts per chromosome
#   filename: Prefix of initial bam
#
# Returns:
#   read_sum_file: dx link to file containing summary of counts
####################################################################
sum_reads() {

  set -e -x -o pipefail

  printf "Value of read file array %s" "${readfiles[@]}"
  echo "Filename: ${filename}"

  #
  # Sum files
  # ---------
  # Read files are download to stdout and fed to stdin of the cat command
  # A 'dx download -o -' pipe approach can be used with any command that support stdin
  #

  echo "Summing values in files and creating output read file"
  for read_f in "${readfiles[@]}"; do
    echo "${read_f}"
    dx download "${read_f}" -o - >> chromosome_result.txt
  done

  count_file="${filename}_chromosome_count.txt"
  total=$(awk '{s+=$2} END {print s}' chromosome_result.txt)
  echo "Total reads: ${total}" >> "${count_file}"

  #
  # Upload files
  # ------------
  #

  readfile_name=$(dx upload "${count_file}" --brief)
  dx-jobutil-add-output read_sum_file "${readfile_name}" --class=file
}
