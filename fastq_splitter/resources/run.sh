#!/bin/bash -e

PATH="/:${PATH}"
HOME="`pwd`"

#
# Inputs
#

input_id="`jshon -e fastqgz -e '$dnanexus_link' -u < job_input.json`"
input_name="`dx describe $input_id --name`"

dx download "$input_id" -o reads.fq.gz --no-progress

reads_per_chunk="`jshon -e reads_per_chunk < job_input.json`"

#
# Processing
#

let lines_per_chunk=reads_per_chunk*4

zcat reads.fq.gz | split -a 10 -l $lines_per_chunk -d - chunk
gzip chunk*

#
# Outputs
#

output_prefix="$input_name"
output_prefix="${output_prefix%.gz}"
output_prefix="${output_prefix%.fastq}"
output_prefix="${output_prefix%.fq}"

output_files='['
index=0
for file in chunk*
do
  let index=index+1
  id=`dx upload $file -o "$output_prefix.${index}.fq.gz" --brief --no-progress`
  output_files="$output_files"'{"$dnanexus_link": "'$id'"},'
done
output_files="${output_files%?}]"

echo '{"fastqgz_chunks": '"$output_files"'}' > "$HOME/job_output.json"
