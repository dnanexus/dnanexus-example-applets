#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x

#
# Fetch inputs
#
dx download "$baits" -o baits.txt
dx download "$targets" -o targets.txt
dx download "$bam" -o input.bam
dx download "$reference_fastagz" -o reference.fasta.gz

gunzip reference.fasta.gz

#
# Run CalculateHsMetrics
#

baits_name=`dx describe --name "$baits"`

samtools faidx reference.fasta
java -Xmx2g -jar /CalculateHsMetrics.jar BI=baits.txt N="$baits_name" TI=targets.txt I=input.bam O=hsmetrics.txt R=reference.fasta PER_TARGET_COVERAGE=pertarget.txt

#
# Upload outputs
#

input_name=`dx describe --name "$bam"`
input_name="${input_name%.bam}"

hsmetrics_fileid=`dx upload hsmetrics.txt -o "$input_name.hsmetrics.txt" --brief --no-progress`
dx-jobutil-add-output "hsmetrics_file" "$hsmetrics_fileid"

pertarget_fileid=`dx upload pertarget.txt -o "$input_name.pertarget_hsmetrics.txt" --brief --no-progress`
dx-jobutil-add-output "pertarget_hsmetrics_file" "$pertarget_fileid"
