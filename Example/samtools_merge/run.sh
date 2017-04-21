#!/bin/bash
#
# Copyright (C) 2013 DNAnexus, Inc.
#   This file is part of dnanexus-example-applets.
#   You may use this file under the terms of the Apache License, Version 2.0;
#   see the License.md file for more information.


# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x

# Inputs
#
# Iterate through the "sorted_bam_files" array
index=0
for file_id in "${sorted_bam_files[@]}"
do
  let index=index+1

  # Locally the files will be called file.0001.bam, file.0002.bam, etc.
  formatted_index=`printf "%04d" $index`
  local_name=file.$formatted_index.bam
  dx download "$file_id" -o $local_name --no-progress

  # The following lines find the common prefix of all the file names (code copied from the Internet).
  file_name=`dx describe "$file_id" --name`
  common_prefix=`printf "%s\n%s\n" "${common_prefix:-$file_name}" "$file_name" | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'`
done

# Processing
#
if [ "${#sorted_bam_files[@]}" == "1" ]
then
  # Special case for a single file
  # Warning: in this case we skip running samtools altogether, thus $params is
  # not taken into account (important if "-n" or "-r" was part of $params)
  mv file.0001.bam out.bam
else
  samtools merge $params out.bam file*.bam
fi

# Outputs
#
# The name of the output is taken from the common prefix, by removing "." and "_", and adding ".bam"
output_name="${common_prefix%.bam}"
output_name="${output_name%.}"
output_name="${output_name%_}"
output_name="${output_name}.bam"
output_id=`dx upload out.bam -o "$output_name" --brief --no-progress`
dx-jobutil-add-output merged_bam "$output_id"
