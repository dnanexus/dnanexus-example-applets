# SAMtools count parallelized by region
This applet performs a basic SAMtools count on a series of sliced(by canonical chromosome) bam files in parallel using wait (Ubuntu 14.04+).

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

## Parallelized run
The bash's built-in [job control](http://tldp.org/LDP/abs/html/x9644.html) system allows for easy management of multiple processes. In this example, we run bash commands in the background as we control maximum job executions in the foreground.

## Job output
Once the input bam has been sliced, counted, and summed the output counts_txt is uploaded using [dx-upload-all-outputs](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-upload-all-outputs). The following directory structure required for dx-upload-all-outputs:
```
├── $HOME
│   ├── out
│       ├── < output name in dxapp.json >
│           ├── output file
```