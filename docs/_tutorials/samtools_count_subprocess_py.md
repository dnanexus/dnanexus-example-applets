---
title: SAMtools count with subprocess
tutorial_type: basic
source: samtools_count_subprocess_py
language: python
---
# SAMtools count with subprocess (python)

## Download bam files

In order to download the input BAM file we:
 - Establish a DXFile handle by providing a DXLink to dxpy.DXFile()
   - mappings_bam is a DXLink containing the file-id of the input BAM.
   - DXFIle handles have the attribute 'name', the filename on the platform.
     We will use the name attribute to retain the filename programmatically.
 - Download the file to the worker by providing a file-id to dxpy.download_dxfile
   - DXFIle handle's `get_id()` function returns the file-id on the platform
```
mappings_bam = dxpy.DXFile(mappings_bam)
mappings_bam_name = mappings_bam.name
dxpy.download_dxfile(mappings_bam.get_id(), mappings_bam_name)
```

## Run samtools view
Subprocess takes either a string command or an array of commands that will be
joined to a string. shell=True causes subprocess to use run a string command using
the program specified by the *SHELL* environmental variable. The default shell
parameter value is False. If an array command is used set `shell=False`. If a string command
is used set `shell=True`.
```
samtools_view_cmd = ["samtools", "view", "-c", mappings_bam_name]
with open(count_txt_name, "w") as f:
    try:
        subprocess.check_call(samtools_view_cmd, stdout=f)
    except subprocess.CalledProcessError as cpe:
        print "Command failed: {cmd}".format(
            cmd=" ".join(samtools_view_cmd))
        raise cpe
```

## Upload result
We now upload the counts file to the platform. [`dxpy.upload_local_file(<filename>)`](http://autodoc.dnanexus.com/bindings/python/current/dxpy_dxfile.html?highlight=upload_local_file#dxpy.bindings.dxfile_functions.upload_local_file) function
will upload the file to the job container, a temporary project which holds onto files
associated with this job. `dxpy.upload_local_file` returns a DXFile object handle for the
uploaded file.
```
counts_txt = dxpy.upload_local_file(count_txt_name)
```

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