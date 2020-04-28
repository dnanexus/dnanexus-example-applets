#!/usr/bin/env python3
#
# This app performs a Pysam count:
#   - Bam file is required
#   - If the index *.bai is not provided it will be generated
#   - Input bam will be split in smaller bams based on chromosomes
#   - Pysam is provided via pip package manager using execDepends
#   - in the the dxapp.json

import dxpy
import pysam
import unicodedata
import re


def get_chr(bam_alignment, canonical=False):
    """Helper function to return canonical chromosomes from SAM/BAM header

    Arguments:
        bam_alignment (pysam.AlignmentFile): SAM/BAM pysam object
        canonical (boolean): Return only canonical chromosomes
    Returns:
        regions (list[str]): Region strings
    """
    regions = []
    headers = bam_alignment.header
    seq_dict = headers['SQ']

    if canonical:
        re_canonical_chr = re.compile(r'^chr[0-9XYM]+$|^[0-9XYM]')
        for seq_elem in seq_dict:
            if re_canonical_chr.match(seq_elem['SN']):
                regions.append(seq_elem['SN'])
    else:
        regions = [''] * len(seq_dict)
        for i, seq_elem in enumerate(seq_dict):
            regions[i] = seq_elem['SN']

    return regions

@dxpy.entry_point('main')
def main(mappings_sorted_bam, canonical_chr, mappings_sorted_bai=None):
    #
    # SECTION: Download inputs
    # --------------------------------------------------------------------------
    # mappings_sorted_bam and mappings_sorted_bai are passed to the main function
    # as parameters for our job. mappings_sorted_bam and mappings_sorted_bai are
    # dictionary objects with key=dnanexus_link and value=<file-id>.
    #
    # We handle file objects from the platform by first creating a DXFile handler.
    # Then performing dxpy.download_dxfile.
    #
    # If index file is not supplied *.bai index will be created with pysam.index
    #
    # DXFIle.name attribute is converted to ASCII since Pysam does not handle Unicode strings.
    #
    print(mappings_sorted_bai)
    print(mappings_sorted_bam)

    mappings_sorted_bam = dxpy.DXFile(mappings_sorted_bam)
    sorted_bam_name = mappings_sorted_bam.name
    dxpy.download_dxfile(mappings_sorted_bam.get_id(),
                         sorted_bam_name)
    ascii_bam_name = unicodedata.normalize(  # Pysam requires ASCII not Unicode string.
        'NFKD', sorted_bam_name).encode('ascii', 'ignore').decode('ascii')

    if mappings_sorted_bai is not None:
        mappings_sorted_bai = dxpy.DXFile(mappings_sorted_bai)
        dxpy.download_dxfile(mappings_sorted_bai.get_id(),
                             mappings_sorted_bai.name)
    else:
        pysam.index(ascii_bam_name)

    #
    # SECTION: Get chromosomes regions
    # --------------------------------------------------------------
    # Generate Pysam Alignmentfile object.
    #
    # Obtain regions to count.

    mappings_obj = pysam.AlignmentFile(ascii_bam_name, "rb")
    regions = get_chr(mappings_obj, canonical_chr)

    #
    # SECTION: Perform basic pysam count.
    # --------------------------------------------------------------
    # Iterate over regions and sum results of pysam.count().

    total_count = 0
    count_filename = "{bam_prefix}_counts.txt".format(
        bam_prefix=ascii_bam_name[:-4])

    with open(count_filename, "w") as f:
        for region in regions:
            temp_count = mappings_obj.count(region=region)
            f.write("{region_name}: {counts}\n".format(
                region_name=region, counts=temp_count))
            total_count += temp_count

        f.write("Total reads: {sum_counts}".format(sum_counts=total_count))

    #
    # SECTION:Output
    # ----------------------------------------------------------------------------
    # Upload generated count file as counts_txt output specified in the dxapp.json

    counts_txt = dxpy.upload_local_file(count_filename)
    output = {}
    output["counts_txt"] = dxpy.dxlink(counts_txt)

    return output


dxpy.run()
