This applet slices a BAM file by canonical chromosome then performs a Parallel `samtools view -c` using xargs. Type `man xargs` for general usage information.

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
/
├── usr
│   ├── bin
│       ├── < samtools binary >
├── home
│   ├── dnanexus
```
/usr/bin is part of the `$PATH` variable, so in our script, we can reference the samtools command directly, `samtools view -c ...`

## Parallel run
### Splice BAM
First, we download our BAM file and slice it by canonical chromosome, writing the `*bam` file names to another file.

In order to split a BAM by regions, we need to have a `*.bai` index. You can either create an app(let) which takes the `*.bai` as an input or generate a `*.bai` in the applet. In this tutorial, we generate the `*.bai` on the fly, sorting if necessary.

<!-- SECTION: Download and segment bam file -->

### Xargs SAMtools view
In the previous section, we recorded the name of each sliced BAM file into a record file. Now we will perform a `samtools view -c` on each slice using the record file as input.
<-- SECTION: Xargs SAMtools count -->

### Upload results
The results file is uploaded using the standard bash process.
1.  Upload a file to the job execution's container.
2.  Provide the DNAnexus link as a job's output using `dx-jobutil-add-output <output name>`
<!-- SECTION: Upload result file -->
