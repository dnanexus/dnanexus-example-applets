#!/bin/bash -e

PATH="/:${PATH}"
HOME="`pwd`"

#
# Inputs
#
index=0
for file_id in `jshon -e sorted_bam_files -a -e '$dnanexus_link' -u < job_input.json`
do
  let index=index+1
  file_name="`dx describe $file_id --name`"
  local_name=file.`printf '%04d' $index`.bam
  common_prefix=`printf "%s\n%s\n" "${common_prefix:-$file_name}" "$file_name" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'`

  dx download "$file_id" -o $local_name --no-progress
done
params="`jshon -e params -u < job_input.json 2>/dev/null || echo ''`"

output_name="${common_prefix%.bam}"
output_name="${output_name%.}"
output_name="${output_name%_}"
output_name="${output_name}.bam"

#
# Processing
#

samtools merge $params out.bam file*.bam

#
# Outputs
#

output_id=`dx upload out.bam -o "$output_name" --brief --no-progress`

echo '{"merged_bam": {"$dnanexus_link": "'$output_id'"}}' > "$HOME/job_output.json"
