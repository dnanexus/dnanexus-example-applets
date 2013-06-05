#!/usr/bin/env python
#
# Copyright (C) 2013 DNAnexus, Inc.
#   This file is part of dnanexus-example-applets.
#   You may use this file under the terms of the Apache License, Version 2.0;
#   see the License.md file for more information.


# bwa_aligner 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# See http://wiki.dnanexus.com/Developer-Tutorials/Intro-to-Building-Apps
# for instructions on how to modify this file.
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import os
import dxpy

import subprocess
from multiprocessing import Pool, cpu_count

@dxpy.entry_point('main')
def main(left_reads, indexed_reference, right_reads=None, samse_sampe_params="-r '@RG\tID:1\tPL:ILLUMINA\tPU:None\tLB:1\tSM:1'", aln_params=''):

    # The following line(s) initialize your data object inputs on the platform
    # into dxpy.DXDataObject instances that you can start using immediately.

    left_reads = dxpy.DXFile(left_reads)
    indexed_reference = dxpy.DXFile(indexed_reference)

    if right_reads != None:
        right_reads = dxpy.DXFile(right_reads)

    # The following line(s) download your file inputs to the local file system
    # using variable names for the filenames.

    dxpy.download_dxfile(left_reads.get_id(), "left.fq.gz")
    subprocess.check_call("gzip -d left.fq.gz", shell=True)
    if right_reads != None:
        dxpy.download_dxfile(right_reads.get_id(), "right.fq.gz")
        subprocess.check_call("gzip -d right.fq.gz", shell=True)
    dxpy.download_dxfile(indexed_reference.get_id(), "ref.tar.gz")

    # Fill in your application code here.

    name = left_reads.describe()['name'].rstrip(".gz").rstrip(".fq").rstrip(".fastq")

    subprocess.check_call("tar -xzvf ref.tar.gz", shell=True)

    subprocess.check_call("bwa aln ref.fa left.fq -t %d %s > left.sai" % (cpu_count(), aln_params), shell=True)
    if right_reads == None:
        subprocess.check_call("bwa samse ref.fa left.sai left.fq %s > output.sam" % (samse_sampe_params), shell=True)
    else:
        subprocess.check_call("bwa aln ref.fa right.fq -t %d %s > right.sai" % (cpu_count(), aln_params), shell=True)
        subprocess.check_call("bwa sampe ref.fa left.sai right.sai left.fq right.fq %s > output.sam" % (samse_sampe_params), shell=True)

    subprocess.check_call("samtools view -bS output.sam > output.bam", shell=True)
    subprocess.check_call("samtools sort output.bam %s" % name, shell=True)

    # The following line(s) use the Python bindings to upload your file outputs
    # after you have created them on the local file system.  It assumes that you
    # have used the output field name for the filename for each output, but you
    # can change that behavior to suit your needs.

    BAM = dxpy.upload_local_file("%s.bam" % name);

    # The following line fills in some basic dummy output and assumes
    # that you have created variables to represent your output with
    # the same name as your output fields.

    output = {}
    output["BAM"] = dxpy.dxlink(BAM)

    return output

dxpy.run()
