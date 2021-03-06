#!/usr/bin/env python
#
# Copyright (C) 2013 DNAnexus, Inc.
#   This file is part of dnanexus-example-applets.
#   You may use this file under the terms of the Apache License, Version 2.0;
#   see the License.md file for more information.


# picard_mark_duplicates 0.0.1
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
def main(BAM, params='ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT'):

    # The following line(s) initialize your data object inputs on the platform
    # into dxpy.DXDataObject instances that you can start using immediately.

    BAM = dxpy.DXFile(BAM)

    # The following line(s) download your file inputs to the local file system
    # using variable names for the filenames.

    dxpy.download_dxfile(BAM.get_id(), "input.bam")

    # The following line extracts the name from the file object so that
    # outputs can be named intelligently. It is not automatically generated by
    # the app wizard.

    name = BAM.describe()['name'].rstrip(".bam").rstrip(".dedup")

    # Fill in your application code here.

    subprocess.check_call("java -Xmx4g -jar /opt/jar/MarkDuplicates.jar INPUT=input.bam OUTPUT=%s.dedup.bam METRICS_FILE=%s.dedup.metrics.txt %s" % (name, name, params), shell=True)


    # The following line(s) use the Python bindings to upload your file outputs
    # after you have created them on the local file system.  It assumes that you
    # have used the output field name for the filename for each output, but you
    # can change that behavior to suit your needs.

    BAM = dxpy.upload_local_file("%s.dedup.bam" % name);
    metrics = dxpy.upload_local_file("%s.dedup.metrics.txt" % name);

    # The following line fills in some basic dummy output and assumes
    # that you have created variables to represent your output with
    # the same name as your output fields.

    output = {}
    output["BAM"] = dxpy.dxlink(BAM)
    output["metrics"] = dxpy.dxlink(metrics)

    return output

dxpy.run()
