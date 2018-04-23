---
categories:
- parallel
- bash
date: '2018-04-23'
github_link: https://github.com/Damien-Black/dnanexus-example-applets/tree/master/Tutorials/bash/samtools_count_para_chr_xargs_sh
summary: SAMtools parallel count using xargs
title: Parallel xargs by Chr
type: Document
---
This applet slices a BAM file by canonical chromosome then performs a parallelized `samtools view -c` using xargs. Type `man xargs` for general usage information.

## How is the SAMtools dependency provided?
The SAMtools compiled binary is placed directory in the `<applet dir>/resources` directory. Any files found in the `resources/` directory will be uploaded so that they will be present in the root directory of the worker. In our case:
```
├── Applet dir
│   ├── src
│   ├── dxapp.json
│   ├── resources
│       ├── usr
│           ├── bin
│               ├── < samtools binary >
```
When this applet is run on a worker, the `resources/` folder will be placed in the worker's root directory `/`:
```
/
├── usr
│   ├── bin
│       ├── < samtools binary >
├── home
│   ├── dnanexus
```
`/usr/bin` is part of the `$PATH` variable, so in our script, we can reference the samtools command directly, as in `samtools view -c ...`

## Parallel Run
### Splice BAM
First, we download our BAM file and slice it by canonical chromosome, writing the `*bam` file names to another file.

In order to split a BAM by regions, we need to have a `*.bai` index. You can either create an app(let) which takes the `*.bai` as an input or generate a `*.bai` in the applet. In this tutorial, we generate the `*.bai` in the applet, sorting the BAM if necessary.

```bash
  dx download "${mappings_bam}"

  indexsuccess=true
  bam_filename="${mappings_bam_name}"
  samtools index "${mappings_bam_name}" || indexsuccess=false
  if [[ $indexsuccess == false ]]; then
    samtools sort -o "${mappings_bam_name}" "${mappings_bam_name}"
    samtools index "${mappings_bam_name}"
    bam_filename="${mappings_bam_name}"
  fi


  chromosomes=$(samtools view -H "${bam_filename}" | grep "\@SQ" | awk -F '\t' '{print $2}' | awk -F ':' '{if ($2 ~ /^chr[0-9XYM]+$|^[0-9XYM]/) {print $2}}')

  for chr in $chromosomes; do
    samtools view -b "${bam_filename}" "${chr}" -o "bam_${chr}."bam
    echo "bam_${chr}.bam"
  done > bamfiles.txt
```

### Xargs SAMtools view
In the previous section, we recorded the name of each sliced BAM file into a record file. Now we will perform a `samtools view -c` on each slice using the record file as input.
```bash
  counts_txt_name="${mappings_bam_prefix}_count.txt"

  sum_reads=$(<bamfiles.txt xargs -I {} samtools view -c $view_options '{}' | awk '{s+=$1} END {print s}')
  echo "Total Count: ${sum_reads}" > "${counts_txt_name}"
```

### Upload results
The results file is uploaded using the standard bash process:
1.  Upload a file to the job execution's container.
2.  Provide the DNAnexus link as a job's output using the script `dx-jobutil-add-output <output name>`
```bash
  counts_txt_id=$(dx upload "${counts_txt_name}" --brief)
  dx-jobutil-add-output counts_txt "${counts_txt_id}" --class=file
```
