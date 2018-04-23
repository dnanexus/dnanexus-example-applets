#!/bin/bash
#
# Copyright (C) 2013 DNAnexus, Inc.
#   This file is part of dnanexus-example-applets.
#   You may use this file under the terms of the Apache License, Version 2.0;
#   see the License.md file for more information.


# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x

#
# Fetch and uncompress inputs
#
dx download "$reads_fastqgz" -o reads.fq.gz
gunzip reads.fq.gz
cmd="-r reads.fq"

if [ "$reads2_fastqgz" != "" ]
then
  dx download "$reads2_fastqgz" -o reads2.fq.gz
  gunzip reads2.fq.gz
  cmd="$cmd -p reads2.fq"
fi

#
# Add "-f" option based on PHRED encoding
#
if [ "$phred64" == "true" ]
then
  cmd="$cmd -f fastq-illumina1.3"
else
  cmd="$cmd -f fastq-sanger"
fi

#
# Run flexbar
#
flexbar $cmd -u $max_uncalled -q $trim_quality_threshold -x $trim_left -y $trim_right -n $threads -m $min_length -t output

#
# Upload outputs (paired or unpaired)
#

function output() {
  #
  # Helper function to upload a result
  #
  # Usage: output <local_file> <input_object> <output_field>
  #
  # Will upload <local_file>, calling it after the name of <input_object>
  # (with .fq/.fastq/.gz removed, and _trimmed.fq.gz added), and set it
  # as the uploaded file id as the value of the <output_field>

  local name=`dx describe --name "$2"`
  name="${name%.gz}"
  name="${name%.fq}"
  name="${name%.fastq}"
  name="${name}_trimmed.fq.gz"

  local fileid=`dx upload "$1" -o "$name" --brief --no-progress`
  dx-jobutil-add-output "$3" "$fileid"
}

if [ "$reads2_fastqgz" != "" ]
then
  gzip output_1.fastq
  gzip output_2.fastq
  output output_1.fastq.gz "$reads_fastqgz" trimmed_reads_fastqgz
  output output_2.fastq.gz "$reads2_fastqgz" trimmed_reads2_fastqgz
else
  gzip output.fastq
  output output.fastq.gz "$reads_fastqgz" trimmed_reads_fastqgz
fi
