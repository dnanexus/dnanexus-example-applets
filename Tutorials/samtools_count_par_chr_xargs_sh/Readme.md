# SAMtools count using xargs by chromosome
This applet slices a BAM file by canonical chromosome then performs a parallelized samtools view -c using xargs.

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

## Parallelized run
Review src/samtools_xargs_chr_sh.sh to see how xargs can be used to parallelize on a worker.
