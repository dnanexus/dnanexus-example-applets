---
categories:
- python
date: '2018-04-23'
github_link: https://github.com/Damien-Black/dnanexus-example-applets/tree/master/Tutorials/python/pysam_count
summary: Counts the number of reads in a BAM file.
title: Pysam
type: Document
---
This applet performs a SAMtools count on an input BAM using Pysam, a python wrapper for SAMtools.

## How is Pysam provided?

Pysam is provided through a `pip install` using the pip package manager in the `dxapp.json`'s `runSpec.execDepends` property:
<!-- Since JSON can't be commented cannot autogenerate below. YAML looking good right now -->

```json
{
 "runSpec": {
    ...
    "execDepends": [
      {"name": "pysam",
         "package_manager": "pip",
         "version": "0.9.1.4"
      }
    ]
    ...
 }
```

The `execDepends` value is a JSON array of dependencies to resolve before the applet source code is run. In this applet, we specify `pip` as our package manager and `pysam version 0.9.1.4` as the dependency to resolve. Pysam is installed to `/usr/local/lib/python2.7/dist-packages` and can be imported by our python script.

## Downloading Inputs   
The fields `mappings_sorted_bam` and `mappings_sorted_bai` are passed to the main function as parameters for our job. These parameters are dictionary objects with key-value pair `{"$dnanexus_link": "<file>-<xxxx>"}`. We handle file objects from the platform through [DXFile](http://autodoc.dnanexus.com/bindings/python/current/dxpy_dxfile.html?highlight=dxfile#module-dxpy.bindings.dxfile) handles. If an index file is not supplied, then a `*.bai` index will be created.
```python
    print mappings_sorted_bai
    print mappings_sorted_bam

    mappings_sorted_bam = dxpy.DXFile(mappings_sorted_bam)
    sorted_bam_name = mappings_sorted_bam.name
    dxpy.download_dxfile(mappings_sorted_bam.get_id(),
                         sorted_bam_name)
    ascii_bam_name = unicodedata.normalize(  # Pysam requires ASCII not Unicode string.
        'NFKD', sorted_bam_name).encode('ascii', 'ignore')

    if mappings_sorted_bai is not None:
        mappings_sorted_bai = dxpy.DXFile(mappings_sorted_bai)
        dxpy.download_dxfile(mappings_sorted_bai.get_id(),
                             mappings_sorted_bai.name)
    else:
        pysam.index(ascii_bam_name)
```

## Working with Pysam
Pysam provides several methods that mimic SAMtools commands. In our applet example, we want to focus only on canonical chromosomes. Pysam's object representation of a BAM file is `pysam.AlignmentFile`.
```python
    mappings_obj = pysam.AlignmentFile(ascii_bam_name, "rb")
    regions = get_chr(mappings_obj, canonical_chr)
```

The helper function `get_chr`
```python
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
```

Once we establish a list of canonical chromosomes, we can iterate over them and perform Pysam's version of `samtools view -c`, `pysam.AlignmentFile.count`.
```python
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
```

## Uploading Outputs
Our summarized counts are returned as the job output. We use the `dx-toolkit` python SDK's [`dxpy.upload_local_file`](http://autodoc.dnanexus.com/bindings/python/current/dxpy_dxfile.html?highlight=upload_local_file#dxpy.bindings.dxfile_functions.upload_local_file) function to upload and generate a DXFile corresponding to our tabulated result file.
```python
    counts_txt = dxpy.upload_local_file(count_filename)
    output = {}
    output["counts_txt"] = dxpy.dxlink(counts_txt)

    return output
```

Python job outputs have to be a dictionary of key-value pairs, with the keys being job output names as defined in the `dxapp.json` file and the values being the output values for corresponding output classes. For files, the output type is a [DXLink](https://wiki.dnanexus.com/api-specification-v1.0.0/Details-and-Links#Linking). We use the [`dxpy.dxlink`](http://autodoc.dnanexus.com/bindings/python/current/dxpy_functions.html?highlight=dxlink#dxpy.bindings.dxdataobject_functions.dxlink) function to generate the appropriate DXLink value.
