#!/usr/bin/env python
#
# Copyright (C) 2013 DNAnexus, Inc.
#   This file is part of dnanexus-example-applets.
#   You may use this file under the terms of the Apache License, Version 2.0;
#   see the License.md file for more information.


# picard_collect_insert_size_metrics 0.0.1
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

@dxpy.entry_point('main')
def main(input_SAM, deviations=None, histogram_width=None, min_percent=None, metric_acc_level=None, ref=None, is_sorted=None, stop_after=None):

    # The following line(s) download your file inputs to the local file system
    # using variable names for the filenames.

    dxpy.download_dxfile(input_SAM, "input")
    if ref != None:
        dxpy.download_dxfile(ref, "ref.fa")


    command = "java -Xmx2g -jar /CollectInsertSizeMetrics.jar"
    command += " INPUT=input"
    command += " OUTPUT=insert_distribution.txt"
    command += " HISTOGRAM_FILE=histogram.pdf"
    if deviations != None:
        command += " DEVIATIONS="+str(deviations)
    if histogram_width != None:
        command += " HISTOGRAM_WIDTH="+str(histogram_width)
    if min_percent != None:
        command += " MINIMUM_PCT="+str(histogram_width)
    if metric_acc_level != None:
        for level in metric_acc_level:
            command += " METRIC_ACCUMULATION_LEVEL="+str(level)
    if ref != None:
        command += " REFERENCE_SEQUENCE=ref.fa"
    if is_sorted != None:
        if is_sorted:
            command += " ASSUME_SORTED=true"
        else:
            command += " ASSUME_SORTED=false"
    if stop_after != None:
        command += " STOP_AFTER="+str(stop_after)

    print "Executing:"
    print command

    # CALL the command here:
    subprocess.check_call(command, shell=True)

    # The following line(s) use the Python bindings to upload your file outputs
    # after you have created them on the local file system.  It assumes that you
    # have used the output field name for the filename for each output, but you
    # can change that behavior to suit your needs.

    histogram = dxpy.upload_local_file("histogram.pdf")
    histogram.rename(dxpy.DXFile(input_SAM).describe()['name']+"_histogram.pdf")
    output_dist = dxpy.upload_local_file("insert_distribution.txt")
    output_dist.rename(dxpy.DXFile(input_SAM).describe()['name']+"_insert_dist.txt")

    # The following line fills in some basic dummy output and assumes
    # that you have created variables to represent your output with
    # the same name as your output fields.

    output = {}
    output["histogram"] = dxpy.dxlink(histogram)
    output["output"] = dxpy.dxlink(output_dist)

    return output

dxpy.run()