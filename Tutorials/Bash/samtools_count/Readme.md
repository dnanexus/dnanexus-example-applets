## Download BAM Files
For bash scripts, inputs to a job execution are made environment variables. From our dxapp.json
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
`mappings_bam` will be available as an environmental variable. Use `dx download` to download the bam file. `$mappings_bam` is a [DNAnexus link (dx link)](https://wiki.dnanexus.com/FAQ#What-are-DNAnexus-links,-and-how-are-they-different-from-using-the-data-object-IDs%3F)
containing the file_id of that file. By default, when we download a file,
we will keep the filename of the object on the platform.
<!--SECTION: Download bam files -->

## SAMtools Count
Here, we use the bash helper variable `mappings_bam_name`. For file inputs,
the DNAnexus platform creates a bash variable `[VARIABLE]_name` which holds a string representing
the filename of the object on the platform. Because we have downloaded the
file by default above, this will be the filename of the object on this
worker as well. We further use `[VARIABLE]_prefix`, which will be the filename
of the object, removing any suffixes specified in patterns. In this tutorial,
the only pattern present is `'["*.bam"]'`, so the platform will remove the trailing *".bam"*.
<!--SECTION: Run samtools view -->

## Upload Result
We now upload the data to the platform. This will upload the file into the
job container, a temporary project which holds onto files associated
with the job. when running upload with `--brief`, it will return just the
file-id.  
<!--SECTION: Upload result -->

## Associate With Output
Finally, we tell the system what file should be associated with the output
named `counts_txt` (which was specified in the dxapp.json). The system will then
move this output into whatever folder was specified at runtime in the project
running the job.  
<!--SECTION: Associate with output -->
<!-- INCLUDE: ## Applet Script -->
<!-- FUNCTION: FULL SCRIPT -->