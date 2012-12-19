#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x

# Inputs
#
dx download "$fastqgz" -o reads.fq.gz --no-progress

# Processing
#
let lines_per_chunk=reads_per_chunk*4

zcat reads.fq.gz | split -a 10 -l $lines_per_chunk -d - chunk
gzip chunk*

# Outputs
#
# The output prefix is the input reads name with any ".gz", ".fastq" and ".fq" extensions removed
input_name=`dx describe "$fastqgz" --name`
output_prefix="$input_name"
output_prefix="${output_prefix%.gz}"
output_prefix="${output_prefix%.fastq}"
output_prefix="${output_prefix%.fq}"

# Upload each file and call it prefix.1.fq.gz, etc.
index=0
for file in chunk*
do
  let index=index+1
  id=`dx upload $file -o "$output_prefix.${index}.fq.gz" --brief --no-progress`
  dx-jobutil-add-output fastqgz_chunks "$id" --array
done
