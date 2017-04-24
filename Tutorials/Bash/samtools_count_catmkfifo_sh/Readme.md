# SAMtools count with mkfifo and dx cat (bash)

This applet performs a SAMtools count on an input file while minimizing disk usage. For additional details on using FIFO special files, type `man fifo` into shell.

Downloading as a stream from the platform using [`dx cat`](https://wiki.dnanexus.com/Command-Line-Client/Index-of-dx-Commands#cat)

Upload as a stream to the platform using [dx-upload-all-outputs](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-upload-all-outputs) or [dx upload -](https://wiki.dnanexus.com/Command-Line-Client/Index-of-dx-Commands?q=dx-upload-all-outputs#upload). Make sure to specify --buffer-size if needed.

## How is SAMtools dependency provided?
SAMtools compiled binary is placed directory in the <Applet dir>/resources directory. Any files found in the resources directory will be uploaded so that they will be present in the root directory of workers. In our case:
```
├── Applet dir
│   ├── src
│   ├── dxapp.json
│   ├── resources
│       ├── usr
│           ├── bin
│               ├── < samtools binary >
```
When this applet is run on a worker the resources/ folder will be placed in the worker's root /:
```
.
├── usr
│   ├── bin
│       ├── < samtools binary >
├── home
│   ├── dnanexus
```
/usr/bin is part of the $PATH variable, so in our script, we can reference the samtools command directly, `samtools view -c ...`

