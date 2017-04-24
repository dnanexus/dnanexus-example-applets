#!/usr/bin/env python
# samtoolscount_chrompara_distr.py

"""SAMtools count distributed by regions

This app performs a SAMtools count via the following implementation:
  - Index file is optional and will be created if missing
  - Sorted bam file is specified as input, but the case where a bam is not sorted is handled
  - Regions will be split based on user input on a mem1_ssd1_x2 instance
      - Files will be distributed to mem1_ssd1_x4 instances for a SAMtools count processing
      - Files will be merged in a mem1_ssd1_x2 instance
      - Review dxapp.json to see how instance type declarations are made for entry points
          - Specifically "runSpec" and "systemRequirements"
  - Gather jobs will download result files via dx download pipe
  - SAMtools is added via execDepends
"""

import os
import dxpy
import subprocess
import shutil
import re


class NotIndexedException(Exception):
    pass


def create_region_view_cmd(input_bam, region):
    view_cmd = ['samtools', 'view', '-c', input_bam, region]
    return view_cmd


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


def create_index_file(bam_filename, bam_dxlink):
    """Create Index file.  Sort BAM if needed"""
    print "Creating Index file."
    index_filename = "{bam}.bai".format(bam=bam_filename)
    cmd_index = ['samtools', 'index', bam_filename]
    sorted_filename = bam_filename
    try:
        run_cmd(cmd_index)
    except NotIndexedException:
        print "Sorting BAM"
        sorted_filename = bam_filename[:-4] + '.sorted.bam'
        cmd_sort = [
            'samtools',
            'sort',
            bam_filename,
            bam_filename[:-4] + '.sorted']
        run_cmd(cmd_sort)
        print "Indexing BAM"
        index_cmd = ['samtools', 'index', sorted_filename]
        index_filename = "{sorted_bam_name}.bai".format(
            sorted_bam_name=sorted_filename)
        run_cmd(index_cmd)
    finally:
        index_file_link = dxpy.dxlink(dxpy.upload_local_file(index_filename))
        aligned_sorted_bam = dxpy.dxlink(dxpy.upload_local_file(sorted_filename))
        return aligned_sorted_bam, index_file_link


def parseSAM_header_for_region(bamfile_path):
    """Helper function to match SN regions contained in SAM header.

    Returns:
        regions (list[string]): List of regions in BAM header
    """
    header_cmd = ['samtools', 'view', '-H', bamfile_path]
    print 'parsing SAM headers'
    headers_str = subprocess.check_output(header_cmd)
    rgx = re.compile(r'SN:(\S+)\s')
    regions = rgx.findall(headers_str)
    return regions


def parse_line_for_readcount(line_read):
    """Helper for getting read count from subjob created count files.

    Arguments:
        line_read (str): header line from BAM file.

    Returns:
        Integer for the amount of reads.
    """
    rgx = re.compile(r':\s(\S+)')
    reads = rgx.search(line_read)
    return int(reads.group(1))


@dxpy.entry_point('samtoolscount_bam')
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


@dxpy.entry_point('combine_files')
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


@dxpy.entry_point('main')
def main(mappings_bam, region_size, index_file=None):
    """The 'scatter' subjob of the applet

    The main function will perform logic to distribute our job
    across multiple workers (instances)

    Returns:
        output (dict): Contains key "count_file" with value DXLink to job output file.
    """
    print 'Creating workspace directory to store downloaded files'
    os.mkdir(u'workspace')
    os.chdir(u'workspace')

    mappings_bam_h = dxpy.DXFile(mappings_bam)
    filename = mappings_bam_h.name
    dxpy.download_dxfile(mappings_bam_h.get_id(), filename)

    #
    # Scatter
    # ------------------------------------------------------
    # Split regions into list of <region size> list
    #
    # Create index file if not provided by user.
    # In order to index bam file needs to be sorted already.
    #   Sort BAM if necessary.
    #   Upload dx file to pass to distributed jobs
    #

    regions = parseSAM_header_for_region(filename)
    split_regions = [regions[i:i + region_size]
                     for i in xrange(0, len(regions), region_size)]

    if not index_file:
        mappings_bam, index_file = create_index_file(filename, mappings_bam)

    #
    # Processing
    # -----------------------------------------------------------------------
    # Run subjob for each distributed region.
    #
    # Note: inputs for subjobs are sent as a dictionary with key value pairs:
    #    key: "region_list"   value: [ [], [], ... ](region sections)
    #    key: "mappings_bam"   value: sorted bam
    #    key: "index_file"    value: bam bai index file
    # The dictionary keys must match the input of the subjob
    #
    # Collect outputs for downstream gather job using dxjob.get_output_ref()
    #
    # Note: Programmatically it's possible to intelligently split workload and
    #    create optimized instance types.  dxpy.new_dxjob takes the optional
    #    parameter: instance_type
    #

    print 'creating subjobs'
    subjobs = [dxpy.new_dxjob(
               fn_input={"region_list": split,
                         "mappings_bam": mappings_bam,
                         "index_file": index_file},
               fn_name="samtoolscount_bam")
               for split in split_regions]

    fileDXLinks = [subjob.get_output_ref("readcount_fileDX")
                   for subjob in subjobs]

    #
    # Gather (Post-processing)
    # -------------------------------------------------------------------------
    # Pass DNAnexus object references to post processing job to combine outputs
    #
    # Create dictionary to be returned as output for the job
    # Dictionary must contain keys matching outputs set in dxapp.json
    #

    print 'combining outputs'
    postprocess_job = dxpy.new_dxjob(
        fn_input={"countDXlinks": fileDXLinks, "resultfn": filename},
        fn_name="combine_files")

    countDXLink = postprocess_job.get_output_ref("countDXLink")

    output = {}
    output["count_file"] = countDXLink

    return output


dxpy.run()
