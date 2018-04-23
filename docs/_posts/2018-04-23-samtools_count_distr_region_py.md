---
categories:
- python
- distributed
date: '2018-04-23'
github_link: https://github.com/Damien-Black/dnanexus-example-applets/tree/master/Tutorials/python/samtools_count_distr_region_py
summary: Count number of reads in SAM or BAM format file using a distributed scatter
  gather approach.  Subjobs perform count in parallel based on contigs
title: Distributed by Region (py)
type: Document
---
The applet will create a count of reads from a BAM format file. Documentation to create a distributed applet can be found on the [DNAnexus wiki](https://wiki.dnanexus.com/Developer-Tutorials/Parallelize-Your-App). This readme will focus on the details of this applet.

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

## Entry Points
Distributed python-interpreter apps use python decorators on functions to [declare entry points](https://wiki.dnanexus.com/Developer-Tutorials/Parallelize-Your-App#Adding-Entry-Points-to-Your-Code). This app has the following entry points as decorated functions:

* *main* 
* *samtoolscount_bam*
* *combine_files*

Entry points are executed on a new worker with their own system requirements. In this example, we *split* and *merge* our files on basic mem1_ssd1_x2 instances and perform our own, more intensive, *processing* step on a mem1_ssd1_x4 instance. Instance type can be set in the dxapp.json `runSpec.systemRequirements`:
```json
  "runSpec": {
    ...
    "systemRequirements": {
      "main": {
        "instanceType": "mem1_ssd1_x2"
      },
      "samtoolscount_bam": {
        "instanceType": "mem1_ssd1_x4"
      },
      "combine_files": {
        "instanceType": "mem1_ssd1_x2"
      }
    },
    ...
  }
```
## main
The *main* function scatters by region bins based on user input. If no `*.bai` file is present, the applet generates an index `*.bai`.
```python
    regions = parseSAM_header_for_region(filename)
    split_regions = [regions[i:i + region_size]
                     for i in xrange(0, len(regions), region_size)]

    if not index_file:
        mappings_bam, index_file = create_index_file(filename, mappings_bam)
```
Regions bins are passed to the *samtoolscount_bam* entry point using the [`dxpy.new_dxjob`](http://autodoc.dnanexus.com/bindings/python/current/dxpy_apps.html?highlight=new_dxjob#dxpy.bindings.dxjob.new_dxjob) function.
```python
    print 'creating subjobs'
    subjobs = [dxpy.new_dxjob(
               fn_input={"region_list": split,
                         "mappings_bam": mappings_bam,
                         "index_file": index_file},
               fn_name="samtoolscount_bam")
               for split in split_regions]

    fileDXLinks = [subjob.get_output_ref("readcount_fileDX")
                   for subjob in subjobs]
```
Outputs from the *samtoolscount_bam* entry points are used as inputs for the *combine_files* entry point. The output of the *combine_files* entry point is used as the output of the main entry point.
```python
    print 'combining outputs'
    postprocess_job = dxpy.new_dxjob(
        fn_input={"countDXlinks": fileDXLinks, "resultfn": filename},
        fn_name="combine_files")

    countDXLink = postprocess_job.get_output_ref("countDXLink")

    output = {}
    output["count_file"] = countDXLink

    return output
```

## samtoolscount_bam
This entry point downloads and creates a `samtools view -c` command for each region in the input bin. The dictionary returned from [`dxpy.download_all_inputs()`]() is used to reference input names and paths.
```python
def samtoolscount_bam(region_list, mappings_bam, index_file):
    """Processing function.

    Arguments:
        region_list (list[str]): Regions to count in BAM
        mappings_bam (dict): dxlink to input BAM
        index_file (dict): dxlink to input BAM

    Returns:
        Dictionary containing dxlinks to the uploaded read counts file
    """
    #
    # Download inputs
    # -------------------------------------------------------------------
    # dxpy.download_all_inputs will download all input files into
    # the /home/dnanexus/in directory.  A folder will be created for each
    # input and the file(s) will be download to that directory.
    #
    # In this example our dictionary inputs has the following key, value pairs
    # Note that the values are all list
    #     mappings_bam_path: [u'/home/dnanexus/in/mappings_bam/<bam filename>.bam']
    #     mappings_bam_name: [u'<bam filename>.bam']
    #     mappings_bam_prefix: [u'<bam filename>']
    #     index_file_path: [u'/home/dnanexus/in/index_file/<bam filename>.bam.bai']
    #     index_file_name: [u'<bam filename>.bam.bai']
    #     index_file_prefix: [u'<bam filename>']
    #

    inputs = dxpy.download_all_inputs()

    # SAMtools view command requires the bam and index file to be in the same
    shutil.move(inputs['mappings_bam_path'][0], os.getcwd())
    shutil.move(inputs['index_file_path'][0], os.getcwd())
    input_bam = inputs['mappings_bam_name'][0]

    #
    # Per region perform SAMtools count.
    # --------------------------------------------------------------
    # Output count for regions and return DXLink as job output to
    # allow other entry points to download job output.
    #

    with open('read_count_regions.txt', 'w') as f:
        for region in region_list:
                view_cmd = create_region_view_cmd(input_bam, region)
                region_proc_result = run_cmd(view_cmd)
                region_count = int(region_proc_result[0])
                f.write("Region {0}: {1}\n".format(region, region_count))
    readcountDXFile = dxpy.upload_local_file("read_count_regions.txt")
    readCountDXlink = dxpy.dxlink(readcountDXFile.get_id())

    return {"readcount_fileDX": readCountDXlink}
```
This entry point returns `{"readcount_fileDX": readCountDXlink}`, a JBOR referencing an uploaded text file. This approach to scatter-gather stores the results in files and uploads/downloads the information as needed. This approach exaggerates a scatter-gather for tutorial purposes. You're able to pass types other than **file** such as **int**.
## combine_files
The *main* entry point triggers this subjob, providing the output of *samtoolscount_bam* as an input. This entry point gathers all the files generated by the *samtoolscount_bam* jobs and sums them.
```python
def combine_files(countDXlinks, resultfn):
    """The 'gather' subjob of the applet.

    Arguments:
        countDXlinks (list[dict]): list of DXlinks to process job output files.
        resultfn (str): Filename to use for job output file.

    Returns:
        DXLink for the main function to return as the job output.

    Note: Only the DXLinks are passed as parameters.
    Subjobs work on a fresh instance so files must be downloaded to the machine
    """
    if resultfn.endswith(".bam"):
        resultfn = resultfn[:-4] + '.txt'

    sum_reads = 0
    with open(resultfn, 'w') as f:
        for i, dxlink in enumerate(countDXlinks):
            dxfile = dxpy.DXFile(dxlink)
            filename = "countfile{0}".format(i)
            dxpy.download_dxfile(dxfile, filename)
            with open(filename, 'r') as fsub:
                for line in fsub:
                    sum_reads += parse_line_for_readcount(line)
                    f.write(line)
        f.write('Total Reads: {0}'.format(sum_reads))

    countDXFile = dxpy.upload_local_file(resultfn)
    countDXlink = dxpy.dxlink(countDXFile.get_id())

    return {"countDXLink": countDXlink}
```

{% include important.html content="While the _main_ entry point triggers the _processing_ and _gathering_ entry points, keep in mind the _main_ entry point **doesn't** do any heavy lifting or _processing_. Notice in the `.runSpec` json [above](#Entry-Points) we start with a lightwieght instance, _scale up_ for the processing entry point, then finally _scale down_ for the _gathering_ step." %}
