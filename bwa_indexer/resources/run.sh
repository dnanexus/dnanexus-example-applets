#!/bin/bash -e

PATH="/:${PATH}"
HOME="`pwd`"

#
# Inputs
#

input_id="`jshon -e fastagz -e '$dnanexus_link' -u < job_input.json`"
input_name="`dx describe $input_id --name`"

dx download "$input_id" -o reference.fasta.gz --no-progress

#
# Processing
#

gunzip reference.fasta.gz
if test `stat -c%s reference.fasta` -gt 2000000000
then
  bwa index -a bwtsw reference.fasta
else
  bwa index -a is reference.fasta
fi

tar zcf reference.bwa062-index.tar.gz reference.fasta*

#
#  Outputs
#
output_name="${input_name%.fasta.gz}.bwa062-index.tar.gz"
output_file=`dx upload reference.bwa062-index.tar.gz -o "$output_name" --brief --no-progress`

echo '{"index_targz": {"$dnanexus_link": "'$output_file'"}}' > "$HOME/job_output.json"
