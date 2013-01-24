#!/bin/bash

dx-build-applet -f bwa_aligner
dx-build-applet -f bwa_indexer
dx-build-applet -f bwa_recalibration_pipeline
dx-build-applet -f contigset_to_fasta_gz
dx-build-applet -f fastq_splitter
dx-build-applet -f gatk_apply_variant_recalibration
dx-build-applet -f gatk_realign_and_recalibrate_applet
dx-build-applet -f gatk_variant_annotator
dx-build-applet -f gatk_variant_caller_applet
dx-build-applet -f gatk_variant_recalibration_pipeline
dx-build-applet -f gatk_variant_recalibrator
dx-build-applet -f import_sam_to_visualize
dx-build-applet -f parallel_bwa
dx-build-applet -f picard_mark_duplicates
dx-build-applet -f picard_merge_sam_files
dx-build-applet -f picard_sam_to_fastq
dx-build-applet -f samtools_view
dx-build-applet -f somatic_sniper
dx-build-applet -f split_bam_interchromosomal_pairs
dx-build-applet -f tumor_normal_snp_pipeline


dx upload bwa_aligner.tar.gz
dx upload bwa_indexer.tar.gz
dx upload bwa_recalibration_pipeline.tar.gz
dx upload contigset_to_fasta_gz.tar.gz
dx upload fastq_splitter.tar.gz
dx upload gatk_apply_variant_recalibration.tar.gz
dx upload gatk_realign_and_recalibrate_applet.tar.gz
dx upload gatk_variant_annotator.tar.gz
dx upload gatk_variant_caller_applet.tar.gz
dx upload gatk_variant_recalibration_pipeline.tar.gz
dx upload gatk_variant_recalibrator.tar.gz
dx upload import_sam_to_visualize.tar.gz
dx upload parallel_bwa.tar.gz
dx upload picard_mark_duplicates.tar.gz
dx upload picard_merge_sam_files.tar.gz
dx upload picard_sam_to_fastq.tar.gz
dx upload samtools_view.tar.gz
dx upload somatic_sniper.tar.gz
dx upload split_bam_interchromosomal_pairs.tar.gz
dx upload tumor_normal_snp_pipeline.tar.gz

