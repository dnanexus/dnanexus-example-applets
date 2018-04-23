This applet performs a basic SAMtools count on a series of sliced (by canonical chromosome) BAM files in parallel using `wait` (Ubuntu 14.04+).

## How is the SAMtools dependency provided?
The SAMtools dependency is resolved by declaring an [Apt-Get](https://help.ubuntu.com/14.04/serverguide/apt-get.html) package in the `dxapp.json` `runSpec.execDepends`.
```json
  "runSpec": {
    ...
    "execDepends": [
      {"name": "samtools"}
    ]
  }
```
For additional information, please refer to the [`execDepends` wiki page](https://wiki.dnanexus.com/Execution-Environment-Reference#Software-Packages).

## Debugging
The command `set -e -x -o pipefail` will assist you in debugging this applet:
* `-e` causes the shell to immediately exit if a command returns a non-zero exit code.
* `-x` prints commands as they are executed, which is very useful for tracking the job's status or pinpointing the exact execution failure.
* `-o pipefail` makes the return code the first non-zero exit code. (Typically, the return code of pipes is the exit code of the last command, which can create difficult to debug problems.)
<!-- SECTION: Debugging boilerplate and input download -->
The `*.bai` file was an optional job input. We check for a empty or unset `var` using the bash built-in test `[[ - z ${var}} ]]`. Then, we can download or create a `*.bai` index as needed.

## Parallel Run
Bash's [job control](http://tldp.org/LDP/abs/html/x9644.html) system allows for easy management of multiple processes. In this example, we run bash commands in the background as we control maximum job executions in the foreground.
We place processes in the background using the character `&` after a command.
<!-- SECTION: Parallel SAMtools count by region -->
<!-- SECTION: Wait for background processes to complete -->
## Job Output
Once the input bam has been sliced, counted, and summed, the output `counts_txt` is uploaded using the command [`dx-upload-all-outputs`](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-upload-all-outputs). The following directory structure required for dx-upload-all-outputs is below:
```
├── $HOME
│   ├── out
│       ├── < output name in dxapp.json >
│           ├── output file
```
<!-- INCLUDE: In our applet, we upload all outputs by: -->
<!-- SECTION: Sum and Upload results -->
