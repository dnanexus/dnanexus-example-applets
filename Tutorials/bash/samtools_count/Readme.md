This applet performs a basic `samtools view -c {bam}` command, referred to as "SAMtools count", on the DNAnexus platform.

## Download BAM Files
For bash scripts, inputs to a job execution become environment variables. The inputs from our `dxapp.json` file are formatted as shown below:
```json
{
  "inputSpec": [
    {
      "name": "mappings_bam",
      "label": "Mapping",
      "class": "file",
      "patterns": ["*.bam"],
      "help": "BAM format file."
    }
  ]
}
```
The object `mappings_bam`, a [DNAnexus link](https://wiki.dnanexus.com/FAQ#What-are-DNAnexus-links,-and-how-are-they-different-from-using-the-data-object-IDs%3F)
containing the file ID of that file, will be available as an environmental variable in the applet's execution. Use the command `dx download` to download the BAM file. By default, when we download a file,
we will keep the filename of the object on the platform.
<!--SECTION: Download bam files -->

## SAMtools Count
Here, we use the bash helper variable [`mappings_bam_name`](https://wiki.dnanexus.com/Developer-Tutorials/Sample-Code?bash#Bash-app-helper-variables). For file inputs,
the DNAnexus platform creates a bash variable `[VARIABLE]_name` that holds a string representing
the filename of the object on the platform; because we downloaded the file using default parameters, this will be the filename of the object on this
worker as well. We use another helper variable, `[VARIABLE]_prefix`, the filename
of the object minus any suffixes specified in the input field patterns. From the input spec above,
the only pattern present is `'["*.bam"]'`, so the platform will remove the trailing *".bam"* and create the helper variable `[VARIABLE]_prefix` for our use.
<!--SECTION: Run samtools view -->

## Upload Result
Use [`dx upload`](https://wiki.dnanexus.com/Command-Line-Client/Index-of-dx-Commands#upload) command to upload data to the platform. This will upload the file into the
job container, a temporary project that holds onto files associated
with the job. When running the command `dx upload` with the flag `--brief`, the command will return just the
file ID.
<!--SECTION: Upload result -->
<!-- INCLUDE: {% include note.html content="While job containers are an integral part of the execution process, a deeper discussion goes out of scope of a basic tutorial. Review the [Containers for Execution](https://wiki.dnanexus.com/API-Specification-v1.0.0/Containers-for-Execution) wiki page for more information." %} -->

## Associate With Output
The output of an applet must be declared before the applet is even built. Looking back to the `dxapp.json` file, we see the following:
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
We declared a **file** type output named `counts_txt`. In the applet script, we must tell the system what file should be associated with the output `counts_txt`. On job completion, usually at the end of the script, this file will be copied from the temporary job container to the project that launched the job.  
<!--SECTION: Associate with output -->
