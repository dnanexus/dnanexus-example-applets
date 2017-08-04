---
categories:
- bash
date: '2017-08-04'
title: Bash Helpers
type: Document
---
This applet performs a basic SAMtools count with the aid of bash [helper variables](https://wiki.dnanexus.com/Developer-Tutorials/Sample-Code?bash#Bash-app-helper-variables).

## Download bam files
This script uses the DNAnexus [`dx-download-all-inputs`](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-download-all-inputs) command to download all job inputs. The `dx-download-all-inputs` command will go through all inputs and download into folders which have the pattern
`/home/dnanexus/in/[VARIABLE]/[file or subfolder with files]`.  
```bash
dx-download-all-inputs
```

## Create output directory

We create an output directory in preparation for `dx-upload-all-outputs` DNAnexus command.
```bash
mkdir -p out/counts_txt
```

## Run samtools view

Here, we use the bash helper variable `mappings_bam_path`. `[VARIABLE]_path` bash variable will point to the location of a file after it has been downloaded using `dx-download-all-inputs`. In our tutorial applet, the input variable name `mappings_bam` with filename `my_mappings.bam` will have a helper variable `mappings_bam_path` with value:
    `/home/dnanexus/in/mappings_bam/my_mappings.bam` 
```bash
samtools view -c "${mappings_bam_path}" > out/counts_txt/"${mappings_bam_prefix}.txt"
```

## Upload result

We use the `dx-upload-all-outputs` to upload data to the platform and associate
it as the job's output. `dx-upload-all-outputs` uses the folder pattern
`/home/dnanexus/out/[VARIABLE]`. It will upload the files in that folder
and then associate them as the output corresponding to `[VARIABLE]`. In this case,
the output is called `counts_txt`. [Above](#create-output-directory), we created the folders, and placed the output there.  
```bash
dx-upload-all-outputs
```

## Applet Script
```bash
dx-download-all-inputs

mkdir -p out/counts_txt

samtools view -c "${mappings_bam_path}" > out/counts_txt/"${mappings_bam_prefix}.txt"

dx-upload-all-outputs
```
