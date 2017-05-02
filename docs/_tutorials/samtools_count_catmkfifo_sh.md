---
title: SAMtools count using fifo special files
tutorial_type: basic
source: samtools_count_catmkfifo_sh
language: bash
---

This app performs a SAMtools count via the following implementation:
 * Bam file download is performed via dx cat.
 * Files are handled using FIFO special files.
 * Similar to piping, any storage on the local system is avoided.
 * Files will be uploaded once the job is done via dx-upload-all-outputs.
 * SAMtools is provided as a precompiled binary in resources/usr/bin.

NOTE: 
  When deciding to use FIFO special files to reduce memory requirements on machine
  know that the same limitations that apply to pipes also apply to FIFO files.
   * FIFO special files used with a command/script that requires random access will fail/hang.
   * FIFO special files need both a stdin and stdout to not block execution.

## Debugging setup
---------------
The -e flag causes bash to exit at any point if there is any error,
the -o pipefail flag tells bash to throw an error if it encounters an error within a pipeline,
the -x flag causes bash to output each line as it is executed

Create a working directory so files are easy to find if ssh is needed

In addition we'll echo the JSON inputs to stdout
```
main() {

  set -e -x -o pipefail
  echo "Value of mappings_bam: '${mappings_bam}'"
```

## Downloading as input stream
---------------------------
Using a FIFO special file we create a sequential file stream through cat command

A fifo file operates as a named pipe.
As a result both readinga and writing must happen at the same time.

After a dx download the downloaded BAM filename is stored in the var $sorted_bam_name
```
  mkdir workspace
  mappings_fifo_path="workspace/${mappings_bam_name}"
  mkfifo "${mappings_fifo_path}" # FIFO file is created
  # & runs the command in the background and is needed to ensure the foreground execution isn't blocked
  # Remember fifo files need both a stdin and stdout to not block.
  dx cat "${mappings_bam}" > "${mappings_fifo_path}" &
  input_pid="$!"
```
## Processing input stream to output stream
----------------------------------------------------------
samtools view -c can be performed on the established fifo.

Directory structure created here `~/out/read_count` is required to use dx-upload-all-outputs
All files found in the path `~/out/<output name>` will be uploaded to the corresponding
`<output name>` specified in the dxapp.json
```
  mkdir -p ./out/counts_txt/
```
By providing the samtools view as a reader the named pipe can be completed and used here. we'll establish a stdout to readcount.txt
```
  counts_fifo_path="./out/counts_txt/${mappings_bam_prefix}_counts.txt"

  mkfifo "${counts_fifo_path}" # FIFO file is created, readcount.txt
```
SAMtools count is performed with an `&` to make it a background process. without the `&` our foreground shell execution would be blocked since readcount.txt fifo has a stdin, but no stdout.
```
  samtools view -c "${mappings_fifo_path}" > "${counts_fifo_path}" &
  process_pid="$!"
```
## Provide job output
------------------
By providing the dx-jobutil-add-output as a reader of readcount.txt
all fifo special files in our pipeline can be executed and closed.
```
  dx-upload-all-outputs &
  upload_pid="$!"
```
## Wait for blocked background processes to finish
-----------------------------------------------
Now that all FIFO special files have are open to both
reading and writing background processes can finish.
```
  wait -n  # "$input_pid"
  wait -n  # "$process_pid"
  wait -n  # "$upload_pid"
}
```