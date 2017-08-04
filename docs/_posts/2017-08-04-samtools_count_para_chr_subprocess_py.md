---
categories:
- python
- parallel
date: '2017-08-04'
title: Parallel by Chr (py)
type: Document
---
This applet tutorial will perform a SAMtools count using parallel threads.

In order to take full advantage of the scalability that cloud computing offers, our scripts have to implement the correct methodologies. This applet tutorial will:
1. Install SAMtools
2. Download BAM file
3. Count Regions in Parallel

This applet tutorial code is similar to the [_Parallel Cores SAMtools count tutorial_](/python_parallel_tutorial.html#samtools_count_para_reg_subprocess_py), except `multiprocessing.dummy` is used instead of `multiprocessing`.

## How is SAMtools dependency provided?
SAMtools dependency is resolved by declaring an [Apt-Get](https://help.ubuntu.com/14.04/serverguide/apt-get.html) package in the dxapp.json runSpec.execDepends.
```json
  "runSpec": {
    ...
    "execDepends": [
      {"name": "samtools"}
    ]
  }
```
For additional information, please refer to the [execDepends wiki page](https://wiki.dnanexus.com/Execution-Environment-Reference#Software-Packages).

## Download BAM file

`dxpy.download_all_inputs()` function downloads all input files into the `/home/dnanexus/in` directory. A folder will be created for each input and the file(s) will be download to that directory. For convenience, the `dxpy.download_all_inputs` function returns a dictionary containing the following keys:
* `<var>_path` **string** full absolute path to where the file was downloaded.
* `<var>_name` **string** name of the file, includes extention.
* `<var>_prefix` **string** name of the file minus the longest mattern pattern found in the dxapp.json I/O pattern field.

The path, name, and prefix key-value pattern is repeated for all applet file class inputs specified in the dxapp.json. In this example our dictionary has the following key, value pairs:
```json
{
    mappings_bam_path: [u'/home/dnanexus/in/mappings_bam/SRR504516.bam']
    mappings_bam_name: [u'SRR504516.bam']
    mappings_bam_prefix: [u'SRR504516']
    index_file_path: [u'/home/dnanexus/in/index_file/SRR504516.bam.bai']
    index_file_name: [u'SRR504516.bam.bai']
    index_file_prefix: [u'SRR504516']
}
```
```python
    inputs = dxpy.download_all_inputs()
    shutil.move(inputs['mappings_bam_path'][0], os.getcwd())
```

## Count Regions in Parallel
Before we can perform our parallelized SAMtools count, we must determine the workload for each thread. We arbitrarily set our number of workers to `10` and set the workload per thread to `1` chromosome at a time. There are various ways to achieve multithreaded processing in python. For the sake of simplicity, we use [`multiprocessing.dummy`](https://docs.python.org/2/library/multiprocessing.html#module-multiprocessing.dummy), a convenient wrapper around Python's threading module.
<!-- INCLUDE: {% include note.html content="In addition to Python's `multiprocessing.dummy` 
 module we simplify our multithreaded processing with several helper functions. We won't cover all the helper functions here so feel free to review the applet source code in GitHub to see function implementations." %} -->
```python
    input_bam = inputs['mappings_bam_name'][0]

    bam_to_use = create_index_file(input_bam)
    print "Dir info:"
    print os.listdir(os.getcwd())

    regions = parseSAM_header_for_region(bam_to_use)

    view_cmds = [
        create_region_view_cmd(bam_to_use, region)
        for region
        in regions]

    print 'Parallel counts'
    t_pools = ThreadPool(10)
    results = t_pools.map(run_cmd, view_cmds)
    t_pools.close()
    t_pools.join()

    verify_pool_status(results)
```
Each worker creates a *string* to be called in a `subprocess.Popen` call. We use the `multiprocessing.dummy.Pool.map(<func>, <iterable>)` function to call the helper function `run_cmd` for each *string* in the iterable of *view commands*. Because we perform our multithreaded processing using `subprocess.Popen`, if any process failed we would not be alerted. We verify our closed workers in the `verify_pool_status` helper function.
```python
def verify_pool_status(proc_tuples):
    err_msgs = []
    for proc in proc_tuples:
        if proc[2] != 0:
            err_msgs.append(proc[1])
    if err_msgs:
        raise dxpy.exceptions.AppInternalError("\n".join(err_msgs))
```

{% include important.html content="In this example we use `subprocess.Popen` to process, then, later on, verify our results in `verify_pool_status`. In general, it is considered good practice to use python's built-in subprocess convenience functions. In this case, `subprocess.check_call` would achieve the same goal." %}

## Gather results

Each worker returns a read count of just one region in the BAM file. We sum and output the results as the job output. We use the dx-toolkit python SDK's [`dxpy.upload_local_file`](http://autodoc.dnanexus.com/bindings/python/current/dxpy_dxfile.html?highlight=upload_local_file#dxpy.bindings.dxfile_functions.upload_local_file) function to upload and generate a DXFile corresponding to our result file.
<!-- Gather results and generate applet output -->
For python, job outputs have to be a dictionary of key-value pairs, with the keys being job output names, as defined in the dxapp.json, and the value being the output value for corresponding output class. For files, the output type is a [DXLink](https://wiki.dnanexus.com/api-specification-v1.0.0/Details-and-Links#Linking). We use the [`dxpy.dxlink`](http://autodoc.dnanexus.com/bindings/python/current/dxpy_functions.html?highlight=dxlink#dxpy.bindings.dxdataobject_functions.dxlink) function to generate the appropriate DXLink value.
```python
    resultfn = bam_to_use[:-4] + '_count.txt'
    with open(resultfn, 'w') as f:
        sum_reads = 0
        for res, reg in izip(results, regions):
            read_count = int(res[0])
            sum_reads += read_count
            f.write("Region {0}: {1}\n".format(reg, read_count))
        f.write("Total reads: {0}".format(sum_reads))

    count_file = dxpy.upload_local_file(resultfn)
    output = {}
    output["count_file"] = dxpy.dxlink(count_file)

    return output
```

## Applet Script
```python
import os
import dxpy
import subprocess
import shutil
from itertools import izip
from multiprocessing.dummy import Pool as ThreadPool
import re


class NotIndexedException(Exception):
    pass


def parseSAM_header_for_region(bamfile_path):
    header_cmd = ['samtools', 'view', '-H', bamfile_path]
    print 'parsing SAM headers'
    print " ".join(header_cmd)
    headers_str = subprocess.check_output(header_cmd)
    rgx = re.compile(r'SN:(\S+)\s')
    regions = rgx.findall(headers_str)
    return regions


def run_cmd(cmd_arr):
    proc = subprocess.Popen(
        cmd_arr,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    exit_code = proc.returncode
    if exit_code != 0:
        raise subprocess.CalledProcessError(
            returncode=exit_code,
            cmd=" ".join(cmd_arr),
            output=stdout)
    elif 'is not sorted' in stderr:
        raise NotIndexedException("BAM file is not indexed")
    proc_tuple = (stdout, stderr, exit_code)
    return proc_tuple


def create_index_file(bam_filename):
    print "Creating Index file."
    cmd_index = ['samtools', 'index', bam_filename]
    sorted_filename = bam_filename
    try:
        run_cmd(cmd_index)
    except NotIndexedException:
        print "Sorting BAM"
        sorted_filename = '{prefix}{suffix}'.format(
            prefix=bam_filename[:-4], suffix='.sorted.bam')
        cmd_sort = ['samtools',
                    'sort',
                    bam_filename,
                    bam_filename[:-4] + '.sorted']
        run_cmd(cmd_sort)
        index_cmd = ['samtools', 'index', sorted_filename]
        run_cmd(index_cmd)
    except subprocess.CalledProcessError as cpe:
        raise cpe
    finally:
        return sorted_filename


def verify_pool_status(proc_tuples):
    err_msgs = []
    for proc in proc_tuples:
        if proc[2] != 0:
            err_msgs.append(proc[1])
    if err_msgs:
        raise dxpy.exceptions.AppInternalError("\n".join(err_msgs))


def create_region_view_cmd(input_bam, region):
    view_cmd = ['samtools', 'view', '-c', input_bam, region]
    return view_cmd


@dxpy.entry_point('main')
def main(mappings_bam):
    print u'Creating workspace directory to store downloaded files'
    os.mkdir(u'workspace')
    os.chdir(u'workspace')


    inputs = dxpy.download_all_inputs()
    shutil.move(inputs['mappings_bam_path'][0], os.getcwd())

    input_bam = inputs['mappings_bam_name'][0]

    bam_to_use = create_index_file(input_bam)
    print "Dir info:"
    print os.listdir(os.getcwd())

    regions = parseSAM_header_for_region(bam_to_use)

    view_cmds = [
        create_region_view_cmd(bam_to_use, region)
        for region
        in regions]

    print 'Parallel counts'
    t_pools = ThreadPool(10)
    results = t_pools.map(run_cmd, view_cmds)
    t_pools.close()
    t_pools.join()

    verify_pool_status(results)


    resultfn = bam_to_use[:-4] + '_count.txt'
    with open(resultfn, 'w') as f:
        sum_reads = 0
        for res, reg in izip(results, regions):
            read_count = int(res[0])
            sum_reads += read_count
            f.write("Region {0}: {1}\n".format(reg, read_count))
        f.write("Total reads: {0}".format(sum_reads))

    count_file = dxpy.upload_local_file(resultfn)
    output = {}
    output["count_file"] = dxpy.dxlink(count_file)

    return output
```
