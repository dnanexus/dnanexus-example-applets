---
categories:
- bash
date: '2018-04-23'
github_link: https://github.com/Damien-Black/dnanexus-example-applets/tree/master/Tutorials/bash/samtools_count_bashhelper_sh
summary: SAMtools count implementation using bash helper functions
title: Bash Helpers
type: Document
---
This applet performs a basic SAMtools count with the aid of bash [helper variables](https://wiki.dnanexus.com/Developer-Tutorials/Sample-Code?bash#Bash-app-helper-variables).

## Download BAM Files
Input files are downloaded using the [`dx-download-all-inputs`](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-download-all-inputs) command. The `dx-download-all-inputs` command will go through all inputs and download into folders with the pattern
`/home/dnanexus/in/[VARIABLE]/[file or subfolder with files]`.
```bash
dx-download-all-inputs
```

## Create Output Directory

We create an output directory in preparation for `dx-upload-all-outputs` DNAnexus command in the (Upload Results)[#Upload Result] section.
```bash
mkdir -p out/counts_txt
```

## Run SAMtools View

After executing the `dx-download-all-inputs` command, there are three helper variables created to aid in scripting. For this applet, the input variable name `mappings_bam` with platform filename `my_mappings.bam` will have a helper variables:
```bash
# [VARIABLE]_path the absolute string path to the file.
$ echo $mappings_bam_path
/home/dnanexus/in/mappings_bam/my_mappings.bam
# [VARIABLE]_prefix the file name minus the longest matching pattern in the dxapp.json file
$ echo $mappings_bam_prefix
my_mappings
# [VARIABLE]_name the file name from the platform
$ echo $mappings_bam_name
my_mappings.bam
```
We use the bash helper variable `mappings_bam_path` to reference the location of a file after it has been downloaded using `dx-download-all-inputs`.
```bash
samtools view -c "${mappings_bam_path}" > out/counts_txt/"${mappings_bam_prefix}.txt"
```

## Upload Result

We use the [`dx-upload-all-outputs`](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-upload-all-outputs) command to upload data to the platform and specify
it as the job's output. The `dx-upload-all-outputs` command expects to find file paths matching the pattern
`/home/dnanexus/out/[VARIABLE]/*`. It will upload matching files and then associate them as the output corresponding to `[VARIABLE]`. In this case,
the output is called `counts_txt`. [Above](#create-output-directory) we created the folders, and we can now place the outputs there.  <!-- TODO: Add multiple dx upload all example -->
```bash
dx-upload-all-outputs
```
