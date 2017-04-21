# Pysam tools count (DNAnexus Platform App)

This applet performs a SAMtools count on an input BAM using Pysam, a python wrapper for SAMtools.

## How is Pysam obtained

Pysam is obtained through a `pip install` using the pip package manager in the dxapp.json's `runSpec.execDepends` property:

```
 "runSpec": {
 ...
    "execDepends": [
      {"name": "pysam",
         "package_manager": "pip",
         "version": "0.9.1.4"
      }
    ],
...
```
`execDepends` value is a JSON array of dependencies to resolve before the applet src code is run. In this applet we specify `pip` as our package manager and `pysam version 0.9.1.4` as the dependency to resolve. Pysam is installed to `/usr/local/lib/python2.7/dist-packages` and can be imported by out python script.