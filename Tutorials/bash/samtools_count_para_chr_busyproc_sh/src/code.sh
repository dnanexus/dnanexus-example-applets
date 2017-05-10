#!/bin/bash
#
# This app performs a SAMtools count via the following implementation:
#   - SAMmtools is provided via exec depends. Review dxapp.json for implementation.
#   - samtools view command will be run across multiple processors through control flow.
#     - Parallel by chromosome.
#   - Result will be uploaded via dx-upload-all-outputs
#
# Note: wait -n available in Ubuntu 14.04.  Review dxapp.json "distribution" field
#

main() {

  #
  # SECTION: Debugging boilerplate and input download
  # -------------------------------------------------
  # A workspace directory is created to make files
  # easier to find if ssh to the instance is required
  #

  set -e -x -o pipefail
  echo "Value of mappings_sorted_bam: '${mappings_sorted_bam}'"
  echo "Value of mappings_sorted_bai: '${mappings_sorted_bai}'"

  mkdir workspace
  cd workspace
  dx download "${mappings_sorted_bam}"

  if [ -z "$mappings_sorted_bai" ]; then
    samtools index "$mappings_sorted_bam_name"
  else
    dx download "${mappings_sorted_bai}"
  fi

  #
  # SECTION: Parallel SAMtools count by region
  # ---------------------------------
  # This approach to parallelization in bash requires programatic control of processes
  # to be done from our side.  We will run called to samtools view in the background using
  # "&" bash background process.
  #
  # For the purpose of this tutorial we'll showcase core management
  # While running processes in parallel

  # Get region information from BAM headers
  chromosomes=$(samtools view -H "${mappings_sorted_bam_name}" | grep "\@SQ" | awk -F '\t' '{print $2}' | awk -F ':' '{if ($2 ~ /^chr[0-9XYM]+$|^[0-9XYM]/) {print $2}}')

  for chr in $chromosomes; do
    samtools view -b "${mappings_sorted_bam_name}" "${chr}" -o "bam_${chr}.bam"
    echo "bam_${chr}.bam"
  done > bamfiles.txt

  #$(nproc) are the cores in the instance
  # e.g. mem1_ssd1_x4 instance $(nproc) will equal 4
  busyproc=0
  while read -r b_file; do
    echo "${b_file}"
    if [[ "${busyproc}" -ge "$(nproc)" ]]; then
      echo Processes hit max
      while [[ "${busyproc}" -gt  0 ]]; do
        wait -n # p_id
        busyproc=$((busyproc-1))
      done
    fi
    samtools view -c "${b_file}"> "count_${b_file%.bam}" &
    busyproc=$((busyproc+1))
  done <bamfiles.txt

  #
  # SECTION: Wait for background processes to complete
  # -----------------------------------------
  # When the last samtools view call is made for the last region
  # the loop will exit with active process.  This waits until all
  # processes are done.
  while [[ "${busyproc}" -gt  0 ]]; do
    wait -n # p_id
    busyproc=$((busyproc-1))
  done

  #
  # SECTION: Sum and Upload results
  # --------------
  # Collect counts from 1st line of each file and sum them
  #   Files are aggregated using globs
  # Proper directory structure has to be created for dx-upload-all-outputs to function correctly
  #   For this example, the directory structure is: $HOME/out/<output var specified for this job>/<all files here will be uploaded>
  #
  outputdir="${HOME}/out/counts_txt"
  mkdir -p "${outputdir}"
  cat count* | awk '{sum+=$1} END{print "Total reads = ",sum}' > "${outputdir}/${mappings_sorted_bam_prefix}_count.txt"

  # Upload file
  dx-upload-all-outputs
}
