#!/bin/bash -e

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x

# Inputs
#
dx download "$readsgz" -o reads.fq.gz --no-progress
if [ "$reads2gz" != "" ]
then
  dx download "$reads2gz" -o reads2.fq.gz --no-progress
fi
dx download "$index_targz" -o reference.bowtie-index.tar.gz --no-progress

# Processing
#
tar zxvf reference.bowtie-index.tar.gz
if [ "$reads2gz" != "" ]
then
  bowtie2 -t -p 2 -x reference.fasta -1 reads.fq.gz -2 reads2.fq.gz $params > out.sam
else
  bowtie2 -t -p 2 -x reference.fasta -U reads.fq.gz $params > out.sam
fi
samtools view -bS out.sam > out.bam
samtools sort out.bam out.sorted

#  Outputs
#
reads_name=`dx describe "$readsgz" --name`
output_name="$reads_name"
output_name="${output_name%.gz}"
output_name="${output_name%.fq}"
output_name="${output_name%.fastq}"
output_name="${output_name%_1}"
output_name="${output_name}.bam"
output_file_id=`dx upload out.sorted.bam -o "$output_name" --brief --no-progress`

dx-jobutil-add-output "sorted_bam" "$output_file_id"
