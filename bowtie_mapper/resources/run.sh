#!/bin/bash -e

PATH="/:${PATH}"
HOME="`pwd`"

#
# Inputs
#

reads_id="`jshon -e readsgz -e '$dnanexus_link' -u < job_input.json`"
reads2_id="`jshon -e reads2gz -e '$dnanexus_link' -u < job_input.json 2>/dev/null || echo none`"
genome_id="`jshon -e index_targz -e '$dnanexus_link' -u < job_input.json`"
params="`jshon -e params -u < job_input.json 2>/dev/null || echo ''`"

dx download "$reads_id" -o reads.fq.gz --no-progress
if [ "$reads2_id" != "none" ]; then dx download "$reads2_id" -o reads2.fq.gz --no-progress ; fi
dx download "$genome_id" -o reference.bowtie-index.tar.gz --no-progress

reads_name="`dx describe $reads_id --name`"

#
# Processing
#

tar zxvf reference.bowtie-index.tar.gz
if [ "$reads2_id" != "none" ]
then
  bowtie2 -t -p 2 -x reference.fasta -1 reads.fq.gz -2 reads2.fq.gz $params > out.sam
else
  bowtie2 -t -p 2 -x reference.fasta -U reads.fq.gz $params > out.sam
fi
samtools view -bS out.sam > out.bam
samtools sort out.bam out.sorted

#
#  Outputs
#
output_name="$reads_name"
output_name="${output_name%.gz}"
output_name="${output_name%.fq}"
output_name="${output_name%.fastq}"
output_name="${output_name%_1}"
output_name="${output_name}.bam"
output_file=`dx upload out.sorted.bam -o "$output_name" --brief --no-progress`

echo '{"sorted_bam": {"$dnanexus_link": "'$output_file'"}}' > "$HOME/job_output.json"
