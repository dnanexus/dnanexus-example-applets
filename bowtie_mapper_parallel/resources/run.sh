#!/bin/bash

PATH="/:${PATH}"
HOME="`pwd`"

set -e
set -x

function main() {
  reads_id="`jshon -e readsgz -e '$dnanexus_link' -u < job_input.json`"
  reads2_id="`jshon -e reads2gz -e '$dnanexus_link' -u < job_input.json 2>/dev/null || echo none`"
  genome_id="`jshon -e index_targz -e '$dnanexus_link' -u < job_input.json`"
  bowtie_params="`jshon -e bowtie_params -u < job_input.json 2>/dev/null || echo ' '`"
  merge_params="`jshon -e merge_params -u < job_input.json 2>/dev/null || echo ' '`"

  if [ "$reads2_id" != "none" ]
  then
    split_job_id="`dx run fastq_splitter -ifastqgz=$reads_id -ireads_per_chunk=$reads_per_chunk -y --brief`"
    split2_job_id="`dx run fastq_splitter -ifastqgz=$reads2_id -ireads_per_chunk=$reads_per_chunk -y --brief`"
    dx api job new '{"function": "scatter", "input": {"genome_id": {"$dnanexus_link": "'$genome_id'"}, "bowtie_params": "'"$bowtie_params"'", "merge_params": "'"$merge_params"'", "chunks": {"job": "'$split_job_id'", "field": "fastqgz_chunks"}, "right_chunks": {"job": "'$split2_job_id'", "field": "fastqgz_chunks"}}}' > scatter.json
  else
    split_job_id="`dx run fastq_splitter -ifastqgz=$reads_id -ireads_per_chunk=$reads_per_chunk -y --brief`"
    dx api job new '{"function": "scatter", "input": {"genome_id": {"$dnanexus_link": "'$genome_id'"}, "bowtie_params": "'"$bowtie_params"'", "merge_params": "'"$merge_params"'", "chunks": {"job": "'$split_job_id'", "field": "fastqgz_chunks"}}}' > scatter.json
  fi
  scatter_job_id="`jshon -e id -u < scatter.json`"
  echo '{"sorted_bam": {"job": "'$scatter_job_id'", "field": "sorted_bam"}}' > "$HOME/job_output.json"
}

function scatter() {
  genome_id="`jshon -e genome_id -e '$dnanexus_link' -u < job_input.json`"
  bowtie_params="`jshon -e bowtie_params -u < job_input.json 2>/dev/null || echo ' '`"
  merge_params="`jshon -e merge_params -u < job_input.json 2>/dev/null || echo ' '`"

  chunk_count="`jshon -e chunks -l < job_input.json`"
  sorted_bam_files=""
  for (( i=0; i<$chunk_count; i++ ))
  do
    reads_id="`jshon -e chunks -e $i -e '$dnanexus_link' -u < job_input.json`"
    reads2_id="`jshon -e right_chunks -e $i -e '$dnanexus_link' -u < job_input.json 2>/dev/null || echo 'none'`"
    if [ $reads2_id != 'none' ]
    then
      job_id=`dx run bowtie_mapper -ireadsgz=$reads_id -ireads2gz=$reads2_id -iindex_targz=$genome_id -iparams="$bowtie_params" -y --brief`
    else
      job_id=`dx run bowtie_mapper -ireadsgz=$reads_id -iindex_targz=$genome_id -iparams="$bowtie_params" -y --brief`
    fi
    sorted_bam_files="$sorted_bam_files -isorted_bam_files=${job_id}:sorted_bam"
  done

  samtools_merge_id=`dx run samtools_merge $sorted_bam_files -iparams="$merge_params" -y --brief`
  echo '{"sorted_bam": {"job": "'$samtools_merge_id'", "field": "merged_bam"}}' > "$HOME/job_output.json"
}
