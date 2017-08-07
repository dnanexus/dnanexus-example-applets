---
categories:
- parallel
- bash
date: '2017-08-07'
github_link: https://github.com/Damien-Black/dnanexus-example-applets/tree/master/Tutorials/bash/samtools_count_para_chr_busyproc_sh
title: Parallel by Region (sh)
type: Document
---
This applet performs a basic SAMtools count on a series of sliced(by canonical chromosome) bam files in parallel using `wait` (Ubuntu 14.04+).

## How is SAMtools dependency provided?
SAMtools dependency is resolved by declaring an [Apt-Get](https://help.ubuntu.com/14.04/serverguide/apt-get.html) package in the dxapp.json runSpec.execDepends.
```json
  "runSpec": {
    ...
    "execDepends": [
      {"name": "samtools"}
    ]
  }
```
For additional information, please refer to the [execDepends wiki page](https://wiki.dnanexus.com/Execution-Environment-Reference#Software-Packages).

## Debugging boilerplate and input download
With Bash scripts, you can prevent a lot of headaches with the command `set -e -x -o pipefail`:
* `-e` causes the shell to immediately exit if a command returns a non-zero exit code.
* `-x` prints commands as they are executed, very useful for tracking job status or pinpointing exact execution failure.
* `-o pipefail` normally, the return code of pipes is the exit code of the last command, this can create difficult to debug problems. This option makes the return code the first non-zero exit code.
```bash
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
```
The `*.bai` file was an optional job input. We check for a empty or unset `var` using the bash built-in `[[-z ${var}}]]` test. Then download or create a `*.bai` index as needed.

## Parallelized run
Bash's [job control](http://tldp.org/LDP/abs/html/x9644.html) system allows for easy management of multiple processes. In this example, we run bash commands in the background as we control maximum job executions in the foreground.
We place processes in the background using an `&` after a command.
```bash
  chromosomes=$(samtools view -H "${mappings_sorted_bam_name}" | grep "\@SQ" | awk -F '\t' '{print $2}' | awk -F ':' '{if ($2 ~ /^chr[0-9XYM]+$|^[0-9XYM]/) {print $2}}')

  for chr in $chromosomes; do
    samtools view -b "${mappings_sorted_bam_name}" "${chr}" -o "bam_${chr}.bam"
    echo "bam_${chr}.bam"
  done > bamfiles.txt

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
```
```bash
  while [[ "${busyproc}" -gt  0 ]]; do
    wait -n # p_id
    busyproc=$((busyproc-1))
  done
```
## Job output
Once the input bam has been sliced, counted, and summed the output counts_txt is uploaded using [dx-upload-all-outputs](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-upload-all-outputs). The following directory structure required for dx-upload-all-outputs:
```
├── $HOME
│   ├── out
│       ├── < output name in dxapp.json >
│           ├── output file
```

In our applet, we upload all outputs by:
```bash
  outputdir="${HOME}/out/counts_txt"
  mkdir -p "${outputdir}"
  cat count* | awk '{sum+=$1} END{print "Total reads = ",sum}' > "${outputdir}/${mappings_sorted_bam_prefix}_count.txt"

  dx-upload-all-outputs
```
