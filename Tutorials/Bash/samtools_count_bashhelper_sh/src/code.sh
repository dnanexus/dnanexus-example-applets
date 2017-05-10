#!/bin/bash
# samtools_count_bashhelper_sh 0.0.1

#
# SECTION: Download bam files
# ------------------
# Use dx-download-all-inputs to download the bam file. This script will go
# through all inputs and download into folders which have the pattern
# /home/dnanexus/in/[INPUT]/.
#

dx-download-all-inputs

#
# SECTION: Create output directory
# -----------------------
# Create an output directory, we will discuss the specifics below
#

mkdir -p out/counts_txt

#
# SECTION: Run samtools view
# -----------------
# Here, we use the bash helper variable mappings_bam_path. This bash variable
# will point to the location of a file after it has been downloaded using
# dx-download-all-inputs. So for [VARIABLE]_path, this wll point to
# /home/dnanexus/in/[VARIABLE]/[VARIABLE_FILENAME].
#
# In the case that the filename of the file mappings_bam is "my_mappings.bam",
# mappings_bam_path will be "/home/dnanexus/in/mappings_bam/my_mappings.bam".
#

samtools view -c "${mappings_bam_path}" > out/counts_txt/"${mappings_bam_prefix}.txt"

#
# SECTION: Upload result
# -------------
# We use dx-upload-all-outputs to upload the data to the platform and associate
# it as the output. dx-upload-all-outputs uses the folder pattern
# /home/dnanexus/out/[VARIABLE]. It will upload the files in that folder
# and then mark is as the output corresponding to [VARIABLE]. In this case,
# the output is called counts_txt. Above, we have made the folder which
# corresponds to this, and put the output there.
#

dx-upload-all-outputs