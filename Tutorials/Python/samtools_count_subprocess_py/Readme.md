# SAMtools count with subprocess (python)

## Download bam files

In order to download the input BAM file we:
-   Establish a DXFile handle by providing a DXLink to dxpy.DXFile()
 -   mappings_bam is a DXLink containing the file-id of the input BAM.
 -   DXFIle handles have the attribute 'name', the filename on the platform.
    We will use the name attribute to retain the filename programmatically.
-   Download the file to the worker by providing a file-id to dxpy.download_dxfile
 -   DXFIle handle's `get_id()` function returns the file-id on the platform
<!-- SECTION: Download bam files -->

## Run samtools view
Subprocess takes either a string command or an array of commands that will be joined to a string. 
  * `shell=True` causes subprocess to run the command program through shell.
  * `shell=False` causes subprocess to execute command program directly.
The default shell parameter value is False. If an array command is used, `shell=False` is required. If a string command
is used, `shell=True` is required.
<!-- INCLUDE: {% include note.html content="It's good practice to use `shell=False` when possible, especially when user input is involved. However, know your use case, if you need to rely on shell features (Even though Python offers modules to mimic these features!) use `shell=True`." %} -->
In this tutorial we use `shell=True`
<!-- SECTION: Run samtools view -->

## Upload result
We now upload the counts file to the platform.
  
[`dxpy.upload_local_file(<filename>)`](http://autodoc.dnanexus.com/bindings/python/current/dxpy_dxfile.html?highlight=upload_local_file#dxpy.bindings.dxfile_functions.upload_local_file) function
will upload the file to the job container, a temporary project which holds onto files
associated with this job. `dxpy.upload_local_file` returns a DXFile object handle for the
uploaded file.
<!-- SECTION: Upload result -->
We then provide the [DXLink](http://autodoc.dnanexus.com/bindings/python/current/dxpy_functions.html?highlight=dxlink#dxpy.bindings.dxdataobject_functions.dxlink) for our uploaded file as the output of this job using the standard process in Python.
<!-- SECTION: Associate with output -->

## How is SAMtools dependency provided?
SAMtools dependency is resolved by declaring an [Apt-Get](https://help.ubuntu.com/14.04/serverguide/apt-get.html) dependency in the dxapp.json runSpec.execDepends.
```
"runSpec": {
     ...
     "execDepends": [
       {"name": "samtools"}
     ]
   }
 ```
For additional information, please refer to the [execDepends wiki page](https://wiki.dnanexus.com/Execution-Environment-Reference#Software-Packages).