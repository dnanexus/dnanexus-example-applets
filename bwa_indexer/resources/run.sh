#!/bin/bash -e

PATH="/:${PATH}"
HOME="`pwd`"

#
# Inputs
#

input_file="`jshon -e fastagz -e '$dnanexus_link' -u < job_input.json`"
output_name="`jshon -e output_name -u < job_input.json`"

mkdir inputs
dx download "$input_file" -o inputs/reference.fasta.gz --no-progress

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

tar zcf reference.bwa6.tgz reference.fasta*

#
#  Outputs
#
output_file=`ua reference.bwa6.tgz -n "$output_name" --do-not-compress`

echo '{"file": {"$dnanexus_link": "'$output_file'"}}' > "$HOME/job_output.json"
