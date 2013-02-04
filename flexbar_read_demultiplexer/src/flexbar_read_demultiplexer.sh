#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x

# Fetch and uncompress inputs
dx download "$left_reads" -o reads.fq.gz
gunzip reads.fq.gz
cmd="-r reads.fq"

if [ "$right_reads" != "" ]
then
  dx download "$right_reads" -o reads2.fq.gz
  gunzip reads2.fq.gz
  cmd="$cmd -p reads2.fq"
fi

if [ "$barcode_reads" != "" ]
then
  dx download "$barcode_reads" -o barcode_reads.fq.gz
  gunzip barcode_reads.fq.gz
  cmd="$cmd -br barcode_reads.fq"
fi

dx download "$barcode_fasta" -o barcodes.fa
cmd="$cmd -b barcodes.fa"

flexbar $cmd -t output $params

for file in output_*
do
  gzip "$file"
  fileid=`dx upload "$file".gz --brief --no-progress`
  dx-jobutil-add-output results "$fileid" --array
done
