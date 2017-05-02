---
title: SAMtools count using xargs by chromosome
tutorial_type: parallel
source: samtools_count_par_chr_xargs_sh
language: bash
---
```
#!/bin/bash
# This app performs a SAMtools count via the following implementation:
#   - Bam file is required
#     - If bam is not sorted then will be automatically sorted then indexed
#   - Input bam will be split in smaller bams based on chromosomes
#   - View options will be determined based on user input.
#     - Will be parsed to provided a valid set of inputs
#   - count will be performed in parallel on each bam file using xargs
#   - Samtools is provided as a precompiled binary in resources/usr/bin

main() {

  #
  # Debugging Setup
  # ---------------
  # Set up logging and exit on error.
  #

  echo "Value of mappings_bam: '${mappings_bam}'"
  set -e -x -o pipefail

  #
  # Download and segment bam file
  # -----------------------------
  # If the input bam was not sorted by coordinate, which we can detect by attempting
  # to index it, don't fail, but instead attempt to sort and index again.
  #
  # A file with the names of each bam file is created for use by xargs
  #

  dx download "${mappings_bam}"

  indexsuccess=true
  bam_filename="${mappings_bam_name}"
  samtools index "${mappings_bam_name}" || indexsuccess=false
  if [[ $indexsuccess == false ]]; then
    #sortedbam="${mappings_bam_prefix}_sorted.bam"
    samtools sort -o "${mappings_bam_name}" "${mappings_bam_name}"
    samtools index "${mappings_bam_name}"
    bam_filename="${mappings_bam_name}"
  fi

  # Pull out chromosomes
  # --------------------
  # Examine the header of the bam file, pulling out just the @SQ lines (which contain
  # the chromosomes, grab the column with just the chromosome name, and match it to a
  # regex designed to capture canonical chromosomes

  chromosomes=$(samtools view -H "${bam_filename}" | grep "\@SQ" | awk -F '\t' '{print $2}' | awk -F ':' '{if ($2 ~ /^chr[0-9XYM]+$|^[0-9XYM]/) {print $2}}')

  for chr in $chromosomes; do
    samtools view -b "${bam_filename}" "${chr}" -o "bam_${chr}."bam
    echo "bam_${chr}.bam"
  done > bamfiles.txt

  #
  # xargs SAMtools count
  # --------------------
  # SAMtools count will be performed in parallel and count summary
  # will be written to file.
  #
  # View command will be built up from user input view_options
  #   $view_options is not quoted to allow word splitting
  #

  counts_txt_name="${mappings_bam_prefix}_count.txt"

  sum_reads=$(<bamfiles.txt xargs -I {} samtools view -c $view_options '{}' | awk '{s+=$1} END {print s}')
  echo "Total Count: ${sum_reads}" > "${counts_txt_name}"

  #
  # Upload result file
  # ------------------
  #

  counts_txt_id=$(dx upload "${counts_txt_name}" --brief)
  dx-jobutil-add-output counts_txt "${counts_txt_id}" --class=file
}
```