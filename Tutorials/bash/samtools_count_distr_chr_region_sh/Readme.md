Documentation to create a distributed applet can be found on the [Developer Tutorials wiki page](https://wiki.dnanexus.com/Developer-Tutorials/Parallelize-Your-App). This readme will focus on the details of this example.

## Entry Points

Distributed bash-interpreter apps use bash functions to [declare entry points](https://wiki.dnanexus.com/Developer-Tutorials/Parallelize-Your-App#Adding-Entry-Points-to-Your-Code). Entry points are executed as subjobs on new workers with their own respective system requirements. This app has the following entry points specified as bash functions:

* *main*
* *count_func*
* *sum_reads*

## main
The *main* function takes the initial `*.bam`, generates an index `*.bai` if needed, and obtains the list of regions from the `*.bam` file. Every 10 regions will be sent, as input, to the *count_func* entry point using [`dx-jobutil-new-job`](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-jobutil-new-job) command.
<!-- SECTION: Download and prepare regions for scatter -->

Job outputs from the *count_func* entry point are referenced as Job Based Object References ([JBOR](https://wiki.dnanexus.com/API-Specification-v1.0.0/Job-Input-and-Output#Job-Dependencies)) and used as inputs for the *sum_reads* entry point.
<!-- SECTION: Merge results -->

Job outputs of the *sum_reads* entry point is used as the output of the *main* entry point via JBOR reference in the [`dx-jobutil-add-output`](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-jobutil-add-output) command.
<!-- SECTION: Output results -->

## count_func
This entry point performs a SAMtools count of the 10 regions passed as input. This execution will be run on a new worker. As a result variables from other functions (e.g. `main()`) will not be accessible here.

Once the output file with counts is created, it is uploaded to the platform and assigned as the entry point's job output `counts_txt` via the command [`dx-jobutil-add-output`](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-jobutil-add-output).
<!-- SECTION: count_func -->

## sum_reads
The *main* entry point triggers this subjob, providing the output of *count_func* as an input JBOR. This entry point gathers all the `readcount.txt` files generated by the *count_func* jobs and sums the totals.

This entry point returns `read_sum` as a JBOR, which is then referenced as job output.
<!-- SECTION: sum_reads -->
<!-- INCLUDE: In the main function, the output is referenced -->
<!-- SECTION: Output results -->
