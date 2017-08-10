This applet performs a SAMtools count on an input BAM using Pysam, a python wrapper for SAMtools.

## How is Pysam provided?

Pysam is provided through a `pip install` using the pip package manager in the `dxapp.json`'s `runSpec.execDepends` property:
<!-- Since JSON can't be commented cannot autogenerate below. YAML looking good right now -->

```json
{
 "runSpec": {
    ...
    "execDepends": [
      {"name": "pysam",
         "package_manager": "pip",
         "version": "0.9.1.4"
      }
    ]
    ...
 }
```

The `execDepends` value is a JSON array of dependencies to resolve before the applet source code is run. In this applet, we specify `pip` as our package manager and `pysam version 0.9.1.4` as the dependency to resolve. Pysam is installed to `/usr/local/lib/python2.7/dist-packages` and can be imported by our python script.

## Downloading Inputs   
The fields `mappings_sorted_bam` and `mappings_sorted_bai` are passed to the main function as parameters for our job. These parameters are dictionary objects with key-value pair `{"$dnanexus_link": "<file>-<xxxx>"}`. We handle file objects from the platform through [DXFile](http://autodoc.dnanexus.com/bindings/python/current/dxpy_dxfile.html?highlight=dxfile#module-dxpy.bindings.dxfile) handles. If an index file is not supplied, then a `*.bai` index will be created.
<!--SECTION: Download inputs -->

## Working with Pysam
Pysam provides several methods that mimic SAMtools commands. In our applet example, we want to focus only on canonical chromosomes. Pysam's object representation of a BAM file is `pysam.AlignmentFile`.
<!--SECTION: Get chromosomes regions -->
<!-- INCLUDE: The helper function `get_chr` -->
<!-- FUNCTION: get_chr -->

Once we establish a list of canonical chromosomes, we can iterate over them and perform Pysam's version of `samtools view -c`, `pysam.AlignmentFile.count`.
<!--SECTION: Perform basic pysam count. -->

## Uploading Outputs
Our summarized counts are returned as the job output. We use the `dx-toolkit` python SDK's [`dxpy.upload_local_file`](http://autodoc.dnanexus.com/bindings/python/current/dxpy_dxfile.html?highlight=upload_local_file#dxpy.bindings.dxfile_functions.upload_local_file) function to upload and generate a DXFile corresponding to our tabulated result file.
<!--SECTION: Output -->

Python job outputs have to be a dictionary of key-value pairs, with the keys being job output names as defined in the `dxapp.json` file and the values being the output values for corresponding output classes. For files, the output type is a [DXLink](https://wiki.dnanexus.com/api-specification-v1.0.0/Details-and-Links#Linking). We use the [`dxpy.dxlink`](http://autodoc.dnanexus.com/bindings/python/current/dxpy_functions.html?highlight=dxlink#dxpy.bindings.dxdataobject_functions.dxlink) function to generate the appropriate DXLink value.
