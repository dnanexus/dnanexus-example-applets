# Parallel Cores SAMtools count

This applet tutorial will perform a SAMtools count using parallel threads.

In order to take full advantage of the scalability that cloud computing offers, our scripts have to implement the correct methodologies. This applet tutorial will:
1. Install SAMtools
2. Download BAM file
3. Count Regions in Parallel

This applet tutorial code is extremely simalar to the [_Parallel Threads SAMtools count tutorial_](/python_parallel_tutorial.html#samtools_count_para_chr_subprocess_py), except `multiprocessing` is used instead of `multiprocessing.dummy`.

Parallel count of reads in BAM format file. 

## How is SAMtools dependency provided?
SAMtools dependency is resolved by declaring an [Apt-Get](https://help.ubuntu.com/14.04/serverguide/apt-get.html) package in the dxapp.json runSpec.execDepends.
```json
  "runSpec": {
    ...
    "execDepends": [
      {"name": "samtools"}
    ]
  }
```
For additional information, please refer to the [execDepends wiki page](https://wiki.dnanexus.com/Execution-Environment-Reference#Software-Packages).
<!-- INCLUDE: ## Applet Script -->
<!-- FUNCTION: FULL SCRIPT -->