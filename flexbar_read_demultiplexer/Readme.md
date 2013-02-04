# Flexbar Read Demultiplexer

This is a simple applet that demultiplexes barcoded reads using the flexbar executable.

To demultiplex reads, you will need an uncompressed fasta file with the barcode specification,
as well as one, two, or three gzipped fastq files with the reads sequences:

- Unpaired data without separate barcode reads: `left_reads`
- Unpaired data with separate barcode reads: `left_reads`, `barcode_reads`
- Paired data without separate barcode reads: `left_reads`, `right_reads`
- Paired data with separate barcode reads: `left_reads`, `right_reads`, `barcode_reads`

The resulting demultiplexed reads (in gzipped fastq format) will be collected in an array
of files in the "results" output field.

The applet source archive includes some example barcode specifications, suitable for the
Illumina TruSeq DNA and small RNA preparation protocols. These can be found inside the
`example_barcodes` directory.

Please refer to the flexbar documentation for additional information, including possible
parameters that can tweak the demultiplexing behavior:

http://sourceforge.net/p/flexbar/wiki/Manual/

The default parameters are the following:
* `-n 2`: Uses 2 threads.
* `-be LEFT_TAIL`: Expects the barcode to be at the left (5-prime) end of the read. IMPORTANT:
Depending on your protocol you may have to change this to `LEFT`, `RIGHT`, `RIGHT_TAIL`, or `ANY`.
Refer to the flexbar manual for more.
* `-u 10`: Allow up to 10 uncalled bases (`.` or `N`) in reads. Reads with too many uncalled
bases will be filtered out and not considered for demultiplexing.

## Inputs

* **The reads to be demultiplexed (or the left reads, for read pairs)** ``left_reads``: ``file``
* **The right reads to be demultiplexed, for read pairs** ``right_reads``: ``file``
* **The barcode reads, if separate barcode reads are created** ``barcode_reads``: ``file``
* **The barcode specificiation (sequences given as an uncompressed fasta)** ``barcode_fasta``: ``file``
* **Other command-line parameters to pass to the flexbar executable** ``params``: ``string``

## Outputs

* **An array containing all resulting files** ``results``: ``array:file``
