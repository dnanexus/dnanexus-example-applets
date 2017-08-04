#!/usr/bin/env python
#
# This app performs a SAMtools count via the following implementation:
#   - dxpy.download_all_inputs is used.
# SAMtools count in parallel per region (multi core) using multiprocess

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
    """Helper function to match SN regions contained in SAM header

    Returns:
        regions (list[string]): list of regions in bam header
    """
    header_cmd = ['samtools', 'view', '-H', bamfile_path]
    print 'parsing SAM headers:', " ".join(header_cmd)
    headers_str = subprocess.check_output(header_cmd)
    rgx = re.compile(r'SN:(\S+)\s')
    regions = rgx.findall(headers_str)
    return regions


def run_cmd(cmd_arr):
    """Run shell command.
    Helper function to simplify the pool.map() call in our parallelization.
    Raises OSError if command specified (index 0 in cmd_arr) isn't valid
    """
    proc = subprocess.Popen(
        cmd_arr,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    exit_code = proc.returncode
    proc_tuple = (stdout, stderr, exit_code)
    return proc_tuple


def verify_pool_status(proc_tuples):
    """
    Helper to verify worker succeeded.

    As failed commands are detected we write the stderr from that command
    to the job_error.json file. This file will be printed to the platform
    job log on App failure.
    """
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

    # Not required.  Making sure all files generated will be
    # in an easy to find location if ssh is needed.
    print u'Creating workspace directory to store downloaded files'
    os.mkdir(u'workspace')
    os.chdir(u'workspace')

    # dxpy.download_all_inputs will download all input files into
    # the /home/dnanexus/in directory.  A folder will be created for each
    # input and the file(s) will be download to that directory.
    #
    # In this example out dictionary inputs has the following key, value pairs
    #     mappings_sorted_bam_path: [u'/home/dnanexus/in/mappings_sorted_bam/SRR504516.bam']
    #     mappings_sorted_bam_name: u'SRR504516.bam'
    #     mappings_sorted_bam_prefix: u'SRR504516'
    #     mappings_sorted_bai_path: u'/home/dnanexus/in/mappings_sorted_bai/SRR504516.bam.bai'
    #     mappings_sorted_bai_name: u'SRR504516.bam.bai'
    #     mappings_sorted_bai_prefix: u'SRR504516'
    inputs = dxpy.download_all_inputs()

    # SAMtools view command required the bam.bai index file to be in the same
    # directory as the bam when specifying regions.
    #
    # When accessing key, value pairs from inputs dictionary note that the
    # values are stored as a list.
    # If our inputs were specified as array:files
    # our list would contain more than 1 element
    shutil.move(inputs['mappings_sorted_bam_path'][0], os.getcwd())
    shutil.move(inputs['mappings_sorted_bai_path'][0], os.getcwd())

    # Create list of regions to parallelize
    # In parallel create workers for each region and generate view cmds
    input_bam = inputs['mappings_sorted_bam_name'][0]
    regions = parse_sam_header_for_region(input_bam)
    view_cmds = [create_region_view_cmd(input_bam, region)
                 for region
                 in regions]

    # Run in parallel using multiprocessing.  multiprocessing.pool class takes
    # an optional argument 'processes' for how many cores the created workers
    # can work on.  When not specified 'processes' = # of cores on system
    # Assigning processes=cpu_count() explicitly is not needed.
    #
    # We process first then verify the results of processing. Throwing an
    # AppInternalError if any commands failed. Review the verify_pool_status()
    # function for details as to how job error logging works.
    print "Number of cpus: {0}".format(cpu_count())
    worker_pool = Pool(processes=cpu_count())
    results = worker_pool.map(run_cmd, view_cmds)
    worker_pool.close()
    worker_pool.join()

    verify_pool_status(results)

    # Write results to file
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

    # Create dictionary to be returned as output for the job
    # Dictionary must contain keys matching outputs set in dxapp.json
    count_file = dxpy.upload_local_file(resultfn)
    output = {}
    output["count_file"] = dxpy.dxlink(count_file)

    return output


dxpy.run()
