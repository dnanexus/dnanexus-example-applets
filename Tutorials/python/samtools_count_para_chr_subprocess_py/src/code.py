#!/usr/bin/env python

import os
import dxpy
import subprocess
import shutil
from itertools import izip
from multiprocessing.dummy import Pool as ThreadPool
import re


class NotIndexedException(Exception):
    pass


# SECTION: Find Regions
def parseSAM_header_for_region(bamfile_path):
    """
    Helper function to match SN regions contained in SAM header

    Returns:
        list of regions in SAM header
    """
    header_cmd = ['samtools', 'view', '-H', bamfile_path]
    print 'parsing SAM headers'
    print " ".join(header_cmd)
    headers_str = subprocess.check_output(header_cmd)
    rgx = re.compile(r'SN:(\S+)\s')
    regions = rgx.findall(headers_str)
    return regions


# SECTION: Subprocess command
def run_cmd(cmd_arr):
    """
    Run shell command. Raises OSError if command specified (index 0 in cmd_arr) isn't valid.
    The stderr check is due to samtools index returning a 0 exit code when indexing fails

    Raises:
        subprocess.CalledProcessError
        NotIndexedException
    """
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


# SECTION: Generate Index
def create_index_file(bam_filename):
    """Create Index file.
    Sorts BAM if needed.
    """
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


# SECTION: Verify Pool Worker
def verify_pool_status(proc_tuples):
    """
    Helper to verify worker succeeded.

    As failed commands are detected we write the stderr from that command
    to the job_error.json file. This file will be printed to the platform
    job log on App failure.
    """
    err_msgs = []
    for proc in proc_tuples:
        if proc[2] != 0:
            err_msgs.append(proc[1])
    if err_msgs:
        raise dxpy.exceptions.AppInternalError("\n".join(err_msgs))


# SECTION: Region Command
def create_region_view_cmd(input_bam, region):
    view_cmd = ['samtools', 'view', '-c', input_bam, region]
    return view_cmd
# SECTION-END


@dxpy.entry_point('main')
def main(mappings_bam):
    print u'Creating workspace directory to store downloaded files'
    os.mkdir(u'workspace')
    os.chdir(u'workspace')

    # SECTION: Download Inputs
    # -------------------------------------------------------------------
    # dxpy.download_all_inputs will download all input files into
    # the /home/dnanexus/in directory.  A folder will be created for each
    # input and the file(s) will be download to that directory.
    #
    # In this example our dictionary inputs has the following key, value pairs
    #     mappings_bam_path: [u'/home/dnanexus/in/mappings_bam/SRR504516.bam']
    #     mappings_bam_name: [u'SRR504516.bam']
    #     mappings_bam_prefix: [u'SRR504516']
    #     index_file_path: [u'/home/dnanexus/in/index_file/SRR504516.bam.bai']
    #     index_file_name: [u'SRR504516.bam.bai']
    #     index_file_prefix: [u'SRR504516']
    #

    inputs = dxpy.download_all_inputs()
    shutil.move(inputs['mappings_bam_path'][0], os.getcwd())

    # SECTION: Parallel by Chromosome using Subprocess.Popen
    # ----------------------------------------------------------
    # Create view commands based on regions.
    # Verify results
    #
    input_bam = inputs['mappings_bam_name'][0]

    bam_to_use = create_index_file(input_bam)
    print "Dir info:"
    print os.listdir(os.getcwd())

    # Regions to parallize
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

    #
    # SECTION: Gather results and generate applet output
    # -----------------------------------------
    # Format and create output file
    #
    # Create dictionary to be returned as output for the job
    # Dictionary must contain keys matching outputs set in dxapp.json
    # and contain values matching output patterns
    #

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
    # SECTION-END


dxpy.run()
