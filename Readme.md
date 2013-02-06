DNAnexus Example Applets
========================

This repository contains some example low-level applets, suitable for
developers to use as inspiration when building their own.

<b><a href="http://wiki.dnanexus.com/Developer-Tutorials/Example-Applets">More information about these applets (wiki.dnanexus.com)</a></b>

Applet Categories
-----------------

The applets below are labeled with brief "tags" that show any significant
features that the applet demonstrates, as follows:

* bash: the applet is written in bash (shell script).
* Invokes other applets: the applet invokes other applets.
* Multiple entry points: the applet makes use of multiple entry points to spawn
  subjobs (frequently used to facilitate parallelism across multiple
  computers).
* Python: the applet is written in Python.

Index of Example Applets
------------------------

### Toolkit

These applets implement simple data transformations that are frequently useful
as parts of larger pipelines.

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/contigset_to_fasta_gz">contigset_to_fasta_gz</a>**: Takes a ContigSet object (one of the DNAnexus type objects which are found in the public project "Reference Genomes"), and produces a gzipped fasta file. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/fastq_splitter">fastq_splitter</a>**: Takes a gzipped fastq file and splits it by number of reads into several smaller files. Useful to parallelize aligning reads. Illustrates: bash

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/flexbar_read_demultiplexer">flexbar_read_demultiplexer</a>**: Demultiplexes indexed (barcoded) reads. Illustrates: bash

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/flexbar_read_trimmer">flexbar_read_trimmer</a>**: Trims reads by quality score and/or position. Illustrates: bash

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/picard_calculate_hs_metrics">picard_calculate_hs_metrics</a>**: Calculates hybrid selection (target enrichment) metrics using the Picard CalculateHsMetrics tool. Illustrates: bash

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/picard_mark_duplicates">picard_mark_duplicates</a>**: Runs MarkDuplicates on a BAM file. Defaults to discarding duplicate reads instead of marking them. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/picard_merge_sam_files">picard_merge_sam_files</a>**: Runs Picard module of same name. Useful as the reduce step of a map-reducte strategy. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/picard_sam_to_fastq">picard_sam_to_fastq</a>**: Wraps Picard module of the same name. Produces gzipped fastq files. Can convert either BAM or SAM files to fastq. Useful for those who receive sequence data in BAM format. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/samtools_merge">samtools_merge</a>**: Runs the merge module of samtools. Useful as the reduce step of a map-reduce strategy. Illustrates: bash

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/samtools_view">samtools_view</a>**: Runs the view module of samtools. Used in the pipeline to split a BAM by region. Could also be used to extract mappings of a certain flag or convert from BAM to SAM. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/split_bam_interchromosomal_pairs">split_bam_interchromosomal_pairs</a>**: This is the only applet which wraps custom written code instead of an open-source program. This splits a BAM file by intra and interchromosomal mappings. This is required by picard_mark_duplicates in a strategy which uses map-reduce to split by genome region. All interchromosmally mapped read pairs must be considered together. Illustrates: Python

### Alignment

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/bowtie_indexer">bowtie_indexer</a>**: Takes a gzipped fasta file and produces a .tar.gz file that bowtie_mapper will take as input. Illustrates: bash

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/bowtie_mapper">bowtie_mapper</a>**: Takes one or two arrays of gzipped fastq files for paired or unpaired reads and an indexed genome from the bowtie_indexer applet and maps them with bowtie. Illustrates: bash

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/bowtie_mapper_parallel">bowtie_mapper_parallel</a>**: Wraps fastq_splitter, bowtie_mapper, and samtools_merge to align reads in parallel. Defaults to provide a worker for every 25,000,000 reads. Illustrates: bash, invokes other applets, multiple entry points

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/bwa_indexer">bwa_indexer</a>**: Takes a gzipped fasta file and produces a .tar.gz file that bwa_aligner will take as input. Defaults to using bwtsw. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/bwa_aligner">bwa_aligner</a>**: Takes one or two arrays of gzipped fastq files for paired or unpaired reads, and runs bwa aln followed by bwa samse/sampe. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/bwa_recalibration_pipeline">bwa_recalibration_pipeline</a>**: Wraps parallel_bwa, split_bam_interchromosomal_pairs, picard_mark_duplicates, gatk_realign_and_recalibrate_applet, and picard_merge_sam_files to align and recalibrate reads in parallel. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/parallel_bwa">parallel_bwa</a>**: Wraps fastq_splitter, bwa_aligner, and picard_merge_sam_files to align reads in parallel. Defaults to providing a worker for every 10,000,000 reads. Illustrates: Python, multiple entry points

### Recalibration and Variant Calling

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/gatk_apply_variant_recalibration">gatk_apply_variant_recalibration</a>**: Runs GATK's ApplyRecalibration on a VCF and a model produced by the gatk_variant_recalibrator applet. This recalibrates the quality of a VCF file and applies a filter to the data. Default filter level is 95% specificity. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/gatk_realign_and_recalibrate_applet">gatk_realign_and_recalibrate_applet</a>**: Performs indel realignment and quality recalibration on a BAM file. Requires dbSNP and known indel files which are available in the datasets directory of the Developer Applets project. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/gatk_variant_annotator_applet">gatk_variant_annotator_applet</a>**: Wraps GATK module of the same name. Can be used to annotate with dbsnp, annotation or comparison VCFs, or a snpEff VCF. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/gatk_variant_caller_applet">gatk_variant_caller_applet</a>**: Runs GATK's UnifiedGenotyper on a BAM file to produce a VCF file. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/gatk_variant_recalibration_pipeline">gatk_variant_recalibration_pipeline</a>**: Wraps the gatk_variant_recalibrator and gatk_apply_variant_recalibration applets to build a model from a VCF(s) and apply it to recalibrate the quality of the VCF. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/gatk_variant_recalibrator">gatk_variant_recalibrator</a>**: Runs GATK's VariantRecalibrator on VCF(s) to produce model files to use in variant recalibration. Illustrates: Python

### Cancer-related

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/tumor_normal_snp_pipeline">tumor_normal_snp_pipeline</a>**: Wraps bwa_recalibration_pipeline and somatic_sniper to align, recalibrate, and call tumor vs. normal snps on two sets of paired reads. Defaults to outputting VCF. Illustrates: Python

**<a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/somatic_sniper">somatic_sniper</a>**: Runs program of same name which takes a tumor BAM file and a normal BAM file and produces variant calls. Defaults to produce output in VCF format. Illustrates: Python
