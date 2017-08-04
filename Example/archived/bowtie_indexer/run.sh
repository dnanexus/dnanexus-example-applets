#!/bin/bash
#
# Copyright (C) 2013 DNAnexus, Inc.
#
#   This file is part of dnanexus-example-applets.
#   You may use the file under the terms of the Apache License, Version 2.0;
#   see the License.md file for more information.

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x

# Inputs
#
dx download "$fastagz" -o reference.fasta.gz --no-progress

# Processing
#
gunzip reference.fasta.gz
bowtie2-build --quiet reference.fasta reference.fasta

tar zcf reference.bowtie-index.tar.gz reference.fasta.*

# Outputs
#
# Make a nice name for the output, by taking the input filename,
# removing '.fasta.gz' from the end, and adding '.bowtie-index.tar.gz'
input_name=`dx describe "$fastagz" --name`
output_name="${input_name%.fasta.gz}.bowtie-index.tar.gz"
output_file_id=`dx upload reference.bowtie-index.tar.gz -o "$output_name" --brief --no-progress`

dx-jobutil-add-output "index_targz" "$output_file_id"
