#!/bin/bash
# samtools_count_git_sh


main() {
  set -e -x -o pipefail

  # Download bam files
  # ------------------
  # Use dx download to download the bam file. $mappings_bam is a dx_link
  # containing the file_id of that file. When a file is specified as an input in
  # the dxapp.json, three convenience variables are created for scripting:
  #    mappings_bam_name   : Basename of the file
  #    mappings_bam_prefix : basename of the file without the file extension
  #    mappings_bam_path   : *[value only valid with dx-download-all]
  #                           Full path up to and including the downloaded file
  dx download "$mappings_bam"

  #
  # Perform SAMtools count
  # -------------------------
  # From the dxapp.json execDepends field:
  #    "destdir": Where the folder SAMtools git repo is cloned to
  #    "build_commands": Commands executed in the cloned directory
  #
  # The compliled SAMtools binary exist in /home/dnanexus/samtools/
  # When we call the samtool view command we must specify the location of
  # the compiled binary, e.g. "samtools/samtool view".
  #
  # To use the command directly without specifying the location of the binary
  # e.g. "samtools view" the compiled binary's location must be in a location
  # specified by the $PATH variable. If desired you can modify the PATH variable:
  #    PATH="${PATH}:\home\dnanexus\samtools"
  #
  count_filename="${mappings_bam_prefix}.txt"
  readcount=$(samtools/samtools view -c "${mappings_bam_name}")
  echo "Total reads: ${readcount}" > "${count_filename}"

  # Upload results
  # --------------
  #
  counts_txt=$(dx upload "${count_filename}" --brief)
  dx-jobutil-add-output counts_txt "${counts_txt}" --class=file
}
