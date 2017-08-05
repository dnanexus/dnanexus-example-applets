---
categories:
- Example Applet
date: '2017-08-04'
title: Flexbar Read Trimmer
type: Document
---
# Flexbar Read Trimmer

This is a simple applet that trims reads (paired or unpaired) by quality and/or position.
The applet is a simple wrapper around flexbar. For more information consult the flexbar
manual at:

http://sourceforge.net/p/flexbar/wiki/Manual/

## Inputs

| Input name | Input label | Description | Corresponding flexbar option |
|------------|-------------|-------------|------------------------------|
|`reads_fastqgz` |Reads|A file, in gzipped fastq format, with the reads to be trimmed (or the left reads, for paired pairs).|`-r`
|`reads2_fastqgz`|Reads (right mates)|(Optional) A file, in gzipped fastq format, with the right reads to be trimmed (for paired reads).|`-p`
|`phred64`       |Quality scores encoded in PHRED-64 instead of PHRED-33?|Select this if the quality scores in the input fastq files are encoded using PHRED-64 (Illumina 1.3-1.5) instead of PHRED-33 (Sanger, Illumina 1.8+)|`-f fastq-sanger` or `-f fastq-illumina1.3`
|`max_uncalled`|Maximum number of Ns tolerated per read|The maximum number of uncalled bases (N or .) allowed per read. Reads with more than those will be filtered.|`-u`
|`trim_quality_threshold`|Trim 3' bases lower than quality|Starting from the 3' end of each read, trim bases until a quality score equal or higher than this threshold is encountered (this number should be given without the 33 or 64 offset added).|`-q`
|`trim_left`|Number of bases to trim from 5' end|The number of bases to trim from the 5' end of the reads|`-x`
|`trim_right`|Number of bases to trim from 3' end|The number of bases to trim from the 3' end of the reads|`-y`
|`threads`|Number of threads|The number of threads to use (default is 2)|`-n`
|`min_length`|Minimum read length allowed|Discard reads whose length, after processing, is lower than this number|`-m`

## Outputs

The applet outputs the trimmed reads file (compressed with gzip).

If right reads were given in the input, the trimmed right reads are output as well (compressed with gzip).
