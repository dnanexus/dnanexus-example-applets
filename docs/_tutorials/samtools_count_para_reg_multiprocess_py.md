---
tutorial_type: parallel
source: samtools_count_para_reg_multiprocess_py
language: python
title: SAMtools count regions in parallel
---
# Parallel Cores SAMtools count

This applet tutorial will perform a SAMtools count using parallel threads.

In order to take full advantage of the scalability that cloud computing offers, our scripts have to implement the correct methodologies. This applet tutorial will:
1. Install SAMtools
2. Download BAM file
3. Count Regions in Parallel

This applet tutorial code is extremely simalar to the [_Parallel Threads SAMtools count tutorial_](/python_parallel_tutorial.html#samtools_count_para_chr_subprocess_py), except `multiprocessing` is used instead of `multiprocessing.dummy`.

Parallel count of reads in BAM format file. 

<hr>## How is SAMtools dependency provided?
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
<hr>
## Applet Script

```python
import os
import dxpy
import subprocess
import shutil
from multiprocessing import Pool, cpu_count
from itertools import izip
import re


def create_region_view_cmd(input_bam, region):
    view_cmd = ['samtools', 'view', '-c', input_bam, region]
    return view_cmd


def parse_sam_header_for_region(bamfile_path):
    header_cmd = ['samtools', 'view', '-H', bamfile_path]
    print 'parsing SAM headers:', " ".join(header_cmd)
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
    proc_tuple = (stdout, stderr, exit_code)
    return proc_tuple


def verify_pool_status(proc_tuples):
    all_succeed = True
    err_msgs = []
    for proc in proc_tuples:
        if proc[2] != 0:
            all_succeed = False
            err_msgs.append(proc[1])
    if err_msgs:
        raise dxpy.exceptions.AppInternalError("\n".join(err_msgs))


@dxpy.entry_point('main')
def main(mappings_sorted_bam, mappings_sorted_bai):

    print u'Creating workspace directory to store downloaded files'
    os.mkdir(u'workspace')
    os.chdir(u'workspace')

    inputs = dxpy.download_all_inputs()

    shutil.move(inputs['mappings_sorted_bam_path'][0], os.getcwd())
    shutil.move(inputs['mappings_sorted_bai_path'][0], os.getcwd())

    input_bam = inputs['mappings_sorted_bam_name'][0]
    regions = parse_sam_header_for_region(input_bam)
    view_cmds = [create_region_view_cmd(input_bam, region)
                 for region
                 in regions]

    print "Number of cpus: {0}".format(cpu_count())
    worker_pool = Pool(processes=cpu_count())
    results = worker_pool.map(run_cmd, view_cmds)
    worker_pool.close()
    worker_pool.join()

    verify_pool_status(results)

    resultfn = inputs['mappings_sorted_bam_name'][0]
    resultfn = (
        resultfn[:-4] + '_count.txt'
        if resultfn.endswith(".bam")
        else resultfn + '_count.txt')
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


dxpy.run()
```
