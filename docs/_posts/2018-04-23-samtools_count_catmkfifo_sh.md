---
categories:
- bash
date: '2018-04-23'
github_link: https://github.com/Damien-Black/dnanexus-example-applets/tree/master/Tutorials/bash/samtools_count_catmkfifo_sh
summary: SAMtools count using mkfifo files and dx cat download stream (bash)
title: Mkfifo and dx cat
type: Document
---
This applet performs a SAMtools count on an input file while minimizing disk usage. For additional details on using FIFO (named pipes) special files, run the command `man fifo` in your shell.


{% include warning.html content="Named pipes require **BOTH** a *stdin* and *stdout* or they will block a process. In these examples, we place incomplete named pipes in background processes so the foreground script process does not block." %}

To approach this use case, let's focus on what we want our applet to do:
1. Stream the BAM file from the platform to a worker.
2. As the BAM is streamed, count the number of reads present.
3. Output the result into a file.
4. Stream the result file to the platform.

## Stream BAM file from the platform to a worker
First, we establish a named pipe on the worker. Then, we stream to *stdin* of the named pipe and download the file as a stream from the platform using [`dx cat`](https://wiki.dnanexus.com/Command-Line-Client/Index-of-dx-Commands#cat).
```bash
  mkdir workspace
  mappings_fifo_path="workspace/${mappings_bam_name}"
  mkfifo "${mappings_fifo_path}" # FIFO file is created
  dx cat "${mappings_bam}" > "${mappings_fifo_path}" &
  input_pid="$!"
```

| FIFO | *stdin* | *stdout* |
|:--------|:-------:|--------:|
| BAM file   | <span style="color: green">**YES**</span>   | <span style="color: red">**NO**</span>   |

## Output BAM file read count
Now that we have created our FIFO special file representing the streamed BAM, we can just call the `samtools` command as we normally would. The `samtools` command reading the BAM would provide our BAM FIFO file with a *stdout*. However, keep in mind that we want to stream the output back to the platform. We must create a named pipe representing our output file too.
```bash
  mkdir -p ./out/counts_txt/

  counts_fifo_path="./out/counts_txt/${mappings_bam_prefix}_counts.txt"

  mkfifo "${counts_fifo_path}" # FIFO file is created, readcount.txt
  samtools view -c "${mappings_fifo_path}" > "${counts_fifo_path}" &
  process_pid="$!"
```

| FIFO | *stdin* | *stdout* |
|:--------|:-------:|--------:|
| BAM file   | <span style="color: green">**YES**</span>   | <span style="color: green">**YES**</span>   |
| output file   | <span style="color: green">**YES**</span>   | <span style="color: red">**NO**</span>   |

The directory structure created here (`~/out/counts_txt`) is required to use the [`dx-upload-all-outputs`](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-upload-all-outputs) command in the next step.
All files found in the path `~/out/<output name>` will be uploaded to the corresponding `<output name>` specified in the `dxapp.json`.

## Stream the result file to the platform
Currently, we've established a stream from the platform, piped the stream into a `samtools` command, and finally outputting the results to another named pipe. However, our background process is still blocked since we lack a *stdout* for our output file. Luckily, creating an upload stream to the platform will resolve this.

We can upload as a stream to the platform using the commands [`dx-upload-all-outputs`](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-upload-all-outputs) or [`dx upload -`](https://wiki.dnanexus.com/Command-Line-Client/Index-of-dx-Commands?q=dx-upload-all-outputs#upload). Make sure to specify `--buffer-size` if needed.
```bash
  mkdir -p ./out/counts_txt/

  counts_fifo_path="./out/counts_txt/${mappings_bam_prefix}_counts.txt"

  mkfifo "${counts_fifo_path}" # FIFO file is created, readcount.txt
  samtools view -c "${mappings_fifo_path}" > "${counts_fifo_path}" &
  process_pid="$!"
```

| FIFO | *stdin* | *stdout* |
|:--------|:-------:|--------:|
| BAM file   | <span style="color: green">**YES**</span>   | <span style="color: green">**YES**</span>   |
| output file   | <span style="color: green">**YES**</span>   | <span style="color: green">**YES**</span>   |


{% include note.html content="Alternatively, `dx upload -` can upload directly from *stdin*. In this example, we would no longer need to have the directory structure required for `dx-upload-all-outputs`." %}


{% include warning.html content="When uploading a file that exists on disk, `dx upload` is aware of the file size and automatically handles any cloud service provider upload chunk requirements. When uploading as a stream, the file size is not automatically known and `dx upload` uses default parameters. While these parameters are fine for most use cases, you may need to specify upload part size with the `--buffer-size` option." %}

## Wait for background processes
Now that our background processes are no longer blocking the rest of the applet's execution, we simply `wait` in the foreground for those processes to finish.
```bash
  wait -n  # "$input_pid"
  wait -n  # "$process_pid"
  wait -n  # "$upload_pid"
```


{% include note.html content="If we didn't wait the app script would running in the foreground would finish and terminate the job! We wouldn't want that." %}

## How is the SAMtools dependency provided?
The SAMtools compiled binary is placed directly in the `<applet dir>/resources` directory. Any files found in the `resources/` directory will be uploaded so that they will be present in the worker's root directory. In our case:
```
├── Applet dir
│   ├── src
│   ├── dxapp.json
│   ├── resources
│       ├── usr
│           ├── bin
│               ├── < samtools binary >
```
When this applet is run on a worker, the `resources/` folder will be placed in the worker's root directory `/`:

```
/
├── usr
│   ├── bin
│       ├── < samtools binary >
├── home
│   ├── dnanexus
```

`/usr/bin` is part of the `$PATH` variable, so we can reference the samtools command directly in our script as `samtools view -c ...`
