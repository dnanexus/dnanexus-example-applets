#!/usr/bin/env python
# samtools_count_subprocess_py 0.0.1

"""SAMtools count using subprocess

This app performs a basic SAMtools count on an input bam file:
    - samtools is provided as an apt-get package (review dxapp.json)
    - samtools command will be excecuted using subprocess
"""

import dxpy
import subprocess


@dxpy.entry_point('main')
def main(mappings_bam):

    #
    # SECTION: Download bam files
    # ------------------
    # In order to download the input BAM file we:
    #  - Establish a DXFile handle by providing a DXLink to dxpy.DXFile()
    #    - mappings_bam is a DXLink containing the file-id of the input BAM.
    #    - DXFIle handles have the attribute 'name', the filename on the platform.
    #      We will use the name attribute to retain the filename programmatically.
    #  - Download the file to the worker by providing a file-id to dxpy.download_dxfile
    #    - DXFIle handle's 'get_id()' function returns the file-id on the platform
    #

    mappings_bam = dxpy.DXFile(mappings_bam)
    mappings_bam_name = mappings_bam.name
    dxpy.download_dxfile(mappings_bam.get_id(), mappings_bam_name)

    #
    # SECTION: Run samtools view
    # -----------------
    # Subprocess takes either a string command or an array of commands that will be
    # joined to a string. shell=True causes subprocess to use run a string command using
    # the program specified by the SHELL environmental variable. The default shell
    # parameter value is False. If an array command is used set shell=False. If a string command
    # is used set shell=True.
    #
    # For convenience, we construct our command array before executing the
    # subprocess.check_call function. subprocess.check_call waits until the
    # command has executed then either: returns if the command exits 0 or raises a
    # subprocess.CalledProcessError exception.
    #

    count_txt_name = "{prefix}_count.txt".format(prefix=mappings_bam_name[:-4])

    samtools_view_cmd = ["samtools", "view", "-c", mappings_bam_name]
    with open(count_txt_name, "w") as f:
        try:
            subprocess.check_call(samtools_view_cmd, stdout=f)
        except subprocess.CalledProcessError as cpe:
            print "Command failed: {cmd}".format(
                cmd=cpe.cmd)
            raise cpe

    #
    # SECTION: Upload result
    # -------------
    # We now upload the counts file to the platform. The dxpy.upload_local_file() function
    # will upload the file to the job container, a temporary project which holds onto files
    # associated with this job. dxpy.upload_local_file returns a DXFile object handle for the
    # uploaded file.

    counts_txt = dxpy.upload_local_file(count_txt_name)

    #
    # SECTION: Associate with output
    # ---------------------
    # Finally, we create our output dictionary. The output dictionary should contain
    # keys corresponding to job output names, specified in the dxapp.json. The value
    # is a DXLink corresponding to the output file.
    #
    # output = {}
    # output["counts_txt"] = dxpy.dxlink(counts_txt)#
    #

    output = {}
    output["counts_txt"] = dxpy.dxlink(counts_txt)

    return output
    # SECTION-END


dxpy.run()
