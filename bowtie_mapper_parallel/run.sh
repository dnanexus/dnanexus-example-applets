#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x

# This applet implements the following pipeline:
#
# fastq_splitter -> bowtie_mapper -> samtools_merge
#
# Running this applet launches the main() function, which is in charge of
# setting up the pipeline.

function main() {
  # Check if this is paired reads
  if [ "$reads2gz" != "" ]
  then
    # Launch fastq_splitter for the left reads
    split_job_id=`dx run "Developer Applets:fastq_splitter" -ifastqgz="$readsgz" -ireads_per_chunk="$reads_per_chunk" -y --brief`
    # Launch fastq_splitter for the right reads
    split2_job_id=`dx run "Developer Applets:fastq_splitter" -ifastqgz="$reads2gz" -ireads_per_chunk="$reads_per_chunk" -y --brief`
    # Launch the 'scatter' function of this applet
    scatter_job_id=`dx-jobutil-new-job scatter -iindex_targz="$index_targz" -ibowtie_params="$bowtie_params" -imerge_params="$merge_params" -ichunks="$split_job_id":fastqgz_chunks -iright_chunks="$split2_job_id":fastqgz_chunks`
  else
    # Launch fastq_splitter for the reads
    split_job_id=`dx run "Developer Applets:fastq_splitter" -ifastqgz="$readsgz" -ireads_per_chunk="$reads_per_chunk" -y --brief`
    # Launch the 'scatter' function of this applet
    scatter_job_id=`dx-jobutil-new-job scatter -iindex_targz="$index_targz" -ibowtie_params="$bowtie_params" -imerge_params="$merge_params" -ichunks="$split_job_id":fastqgz_chunks`
  fi
  # Return back to the system the 'sorted_bam' output of the 'scatter' job.
  dx-jobutil-add-output sorted_bam "$scatter_job_id":sorted_bam
}

function scatter() {
  sorted_bam_files=""
  # Go through all the chunks
  for i in "${!chunks[@]}"
  do
    reads="${chunks[$i]}"
    reads2="${right_chunks[$i]}"
    if [ "$reads2" != "" ]
    then
      # Map this chunk of paired reads
      job_id=`dx run "Developer Applets:bowtie_mapper" -ireadsgz="$reads" -ireads2gz="$reads2" -iindex_targz="$index_targz" -iparams="$bowtie_params" -y --brief`
    else
      # Map this chunk of unpaired reads
      job_id=`dx run "Developer Applets:bowtie_mapper" -ireadsgz="$reads" -iindex_targz="$index_targz" -iparams="$bowtie_params" -y --brief`
    fi
    # Construct the line to give to "dx run samtools_merge"
    sorted_bam_files="$sorted_bam_files -isorted_bam_files=${job_id}:sorted_bam"
  done

  # Launch samtools_merge
  samtools_merge_id=`dx run "Developer Applets:samtools_merge" $sorted_bam_files -iparams="$merge_params" -y --brief`
  # Return back to the system the 'merged_bam' output of samtools_merge
  dx-jobutil-add-output sorted_bam "$samtools_merge_id":merged_bam
}
