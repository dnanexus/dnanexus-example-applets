## How is the SAMtools dependency provided?
The SAMtools dependency is resolved by declaring an [Apt-Get](http://manpages.ubuntu.com/manpages/xenial/man8/apt-get.8.html) package in the `dxapp.json` `runSpec.execDepends` field.
```json
{
  "runSpec": {
    ...
    "execDepends": [
      {"name": "samtools"}
    ]
  }
```
For additional information, please refer to the [`execDepends` documentation ](https://documentation.dnanexus.com/developer/api/running-analyses/io-and-run-specifications#run-specification).

## Download BAM Files

In order to download the input BAM file, we:
- Establish a DXFile handle by providing a DXLink to dxpy.DXFile()
  - `mappings_bam` is a DXLink containing the file-id of the input BAM.
  - `DXFIle` objects inherit from [`DXDataObjects`](http://autodoc.dnanexus.com/bindings/python/current/dxpy_bindings.html#dxpy.bindings.DXDataObject). A `DXDataObject` has the attribute 'name', the filename for an object on the platform. We will use the name attribute to retain the filename programmatically.
- Download the file to the worker by providing a file ID to `dxpy.download_dxfile`
  - `DXFIle` handle's `get_id()` function returns the file ID on the platform.

## Run `samtools view`
The subprocess takes either a string or an array shell command to be executed. 
  * `shell=True` causes the subprocess to run the command program through shell.
  * `shell=False` causes the subprocess to execute command program directly.
The default shell parameter value is False. If an array command is used, `shell=False` is required. If a string command
is used, `shell=True` is required.
<!-- INCLUDE: {% include note.html content="It's good practice to use `shell=False` when possible, especially when user input is involved. That being said, know your use case, if you need to rely on shell features (even though Python offers modules to mimic these features!) use `shell=True`." %} -->
In this tutorial, we use `shell=False`.
<!-- SECTION: Run samtools view -->

## Upload Result  
The [`dxpy.upload_local_file(filename=<filename>)`](http://autodoc.dnanexus.com/bindings/python/current/dxpy_dxfile.html?highlight=upload_local_file#dxpy.bindings.dxfile_functions.upload_local_file) function
will upload the file to the job container, a temporary project which holds onto files
associated with this job. `dxpy.upload_local_file` returns a DXFile object handle for the
uploaded file.
<!-- SECTION: Upload result -->
We then provide the [DXLink](http://autodoc.dnanexus.com/bindings/python/current/dxpy_functions.html?highlight=dxlink#dxpy.bindings.dxdataobject_functions.dxlink) for our uploaded file as the output of this job using the standard process in Python:
1. Use [`dxpy.dxlink(object_id=<Uploaded DXfile object>)`](http://autodoc.dnanexus.com/bindings/python/current/dxpy_functions.html?highlight=dxlink#dxpy.bindings.dxdataobject_functions.dxlink)
2. Set the value in the job output dictionary to the returned `DXLink`
3. `return` the job output dictionary as the entry point function output

<!-- SECTION: Associate with output -->
