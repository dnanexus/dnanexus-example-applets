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

## Overview
### Debugging boilerplate and input download
With Bash scripts, you can prevent a lot of headaches with the command `set -e -x -o pipefail`:
* `-e` causes the shell to immediately exit if a command returns a non-zero exit code.
* `-x` prints commands as they are executed, very useful for tracking job status or pinpointing exact execution failure.
* `-o pipefail` normally, the return code of pipes is the exit code of the last command, this can create difficult to debug problems. This option makes the return code the first non-zero exit code.
<!-- SECTION: Debugging boilerplate and input download -->
The `*.bai` file was an optional job input. We check for a empty or unset `var` using the bash built-in `[[-z ${var}}]]` test. Then download or create a `*.bai` index as needed.

### Parallelized run
Bash's built-in [job control](http://tldp.org/LDP/abs/html/x9644.html) system allows for easy management of multiple processes. In this example, we run bash commands in the background as we control maximum job executions in the foreground.
We place processes in the background using an `&` after a command.
<!-- SECTION: Parallel SAMtools count by region -->
<!-- SECTION: Wait for background processes to complete -->
## Job output
Once the input bam has been sliced, counted, and summed the output counts_txt is uploaded using [dx-upload-all-outputs](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-upload-all-outputs). The following directory structure required for dx-upload-all-outputs:
```
├── $HOME
│   ├── out
│       ├── < output name in dxapp.json >
│           ├── output file
```
<!-- INCLUDE: In our applet, we accomplish this by -->
<!-- SECTION: Sum and Upload results -->
<!-- INCLUDE: ## Applet Script -->
<!-- FUNCTION: FULL SCRIPT -->
