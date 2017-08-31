---
categories:
- bash
date: '2017-08-31'
github_link: https://github.com/Damien-Black/dnanexus-example-applets/tree/master/Tutorials/bash/samtools_count
title: SAMtools count
type: Document
---
This applet performs a basic `samtools view -c {bam}`, referred to as "SAMtools count", on the DNAnexus platform.

## Download BAM Files
For bash scripts, inputs to a job execution are made environment variables. The inputs from our dxapp.json
```json
{
  "inputSpec": [{
    "name": "mappings_bam",
    "label": "Mapping",
    "class": "file",
    "patterns": ["*.bam"],
    "help": "BAM format file."
  }],
}
```
`mappings_bam`, a [DNAnexus link](https://wiki.dnanexus.com/FAQ#What-are-DNAnexus-links,-and-how-are-they-different-from-using-the-data-object-IDs%3F)
containing the file_id of that file, will be available as an environmental variable. Use `dx download` to download the bam file. By default, when we download a file,
we will keep the filename of the object on the platform.
```bash
dx download "${mappings_bam}"
```

## SAMtools Count
Here, we use the bash helper variable [`mappings_bam_name`](https://wiki.dnanexus.com/Developer-Tutorials/Sample-Code?bash#Bash-app-helper-variables). For file inputs,
the DNAnexus platform creates a bash variable `[VARIABLE]_name` that holds a string representing
the filename of the object on the platform; because we downloaded the file using default parameters, this will be the filename of the object on this
worker as well. We use another helper variable, `[VARIABLE]_prefix`, the filename
of the object minus any suffixes specified in the input field patterns. From the input spec above,
the only pattern present is `'["*.bam"]'`, so the platform will remove the trailing *".bam"* and create the helper variable `[VARIABLE]_prefix`.
```bash
readcount=$(samtools view -c "${mappings_bam_name}")
echo "Total reads: ${readcount}" > "${mappings_bam_prefix}.txt"
```

## Upload Result
Use [`dx upload`](https://wiki.dnanexus.com/Command-Line-Client/Index-of-dx-Commands#upload) to upload data to the platform. This will upload the file into the
job container, a temporary project which holds onto files associated
with the job. when running upload with `--brief`, it will return just the
file-id.
```bash
counts_txt_id=$(dx upload "${mappings_bam_prefix}.txt" --brief)
```

{% include note.html content="While job containers are an integral part of the execution process a deeper discussion goes out of scope of a basic tutorial. Review the [Containers for Execution](https://wiki.dnanexus.com/API-Specification-v1.0.0/Containers-for-Execution) wiki page for more information." %}

## Associate With Output
The output of an applet must be declared before the applet is even built. Looking back to the `dxapp.json` file:
```json
{
  "name": "counts_txt",
  "class": "file",
  "label": "Read count file",
  "patterns": [
    "*.txt"
  ],
  "help": "Output file with Total reads as the first line."
}
```
We declared a **file** type output named `counts_txt`. In the applet script, we must tell the system what file should be associated with the output `counts_txt`. On job completion, usually end of script, this file will be copied from the job container to the Project that launched the job.  
```bash
dx-jobutil-add-output counts_txt "${counts_txt_id}" --class=file
```
