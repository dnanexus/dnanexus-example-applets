---
tutorial_type: basic
source: samtools_count_par_chr_subprocess_py
language: python
title: SAMtools count chromosomes in parallel (py)
---
<!-- dx-header -->
SAMtools count chromosomes in parallel (py)

Parallel count of reads in BAM format file. 

## How is SAMtools dependency provided?
SAMtools dependency is resolved by declaring an [Apt-Get](https://help.ubuntu.com/14.04/serverguide/apt-get.html) package in the dxapp.json runSpec.execDepends.
```
  "runSpec": {
    ...
    "execDepends": [
      {"name": "samtools"}
    ]
  }
```
For additional information, please refer to the [execDepends wiki page](https://wiki.dnanexus.com/Execution-Environment-Reference#Software-Packages).