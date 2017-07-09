---
tutorial_type: distributed
source: samtools_count_distr_region_py
language: python
title: SAMtools count scatter gather based on regions
---
# SAMtools count distributed by Chromosome (DNAnexus Platform App)

Distributed count of reads in BAM format file. Documentation to create a distributed applet can be found on the [DNAnexus wiki](https://wiki.dnanexus.com/Developer-Tutorials/Parallelize-Your-App). This readme will focus on the details of this applet.

<hr>## How is SAMtools dependency provided?
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

<hr>## Entry points
Distributed python-interpreter apps use python decorators on functions to [declare entry points](https://wiki.dnanexus.com/Developer-Tutorials/Parallelize-Your-App#Adding-Entry-Points-to-Your-Code). This app has the following entry points as decorated functions:

* *main* 
* *samtoolscount_bam*
* *combine_files*

Entry points are executed on a new worker with their own system requirements. In this example, we *split* and *merge* our files on basic mem1_ssd1_x2 instances and perform our, more intensive, *processing* step on a mem1_ssd1_x4 instance. Instance type can be set in the dxapp.json `runSpec.systemRequirements`:
```
  "runSpec": {
    ...
    "systemRequirements": {
      "main": {
        "instanceType": "mem1_ssd1_x2"
      },
      "samtoolscount_bam": {
        "instanceType": "mem1_ssd1_x4"
      },
      "combine_files": {
        "instanceType": "mem1_ssd1_x2"
      }
    },
    ...
  }
```
<hr>## Overview
### main
The *main* function scatters by region bins, based on user input. If no `*.bai` file is present, The applet generates an index `*.bai`.
