This repository contains some example low-level applets, suitable for
developers to use as inspiration when building their own.

<a href="http://wiki.dnanexus.com/Developer-Tutorials/Example-Applets">More information about these applets (wiki.dnanexus.com)</a>

Index of applets
----------------

<table>
<tr>
<th>Applet</th><th>Notes</th>
</tr>

<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/bowtie_indexer">bowtie_indexer</a></td><td>Takes a gzipped fasta file and produces a .tar.gz file that bowtie_mapper will take as input</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/bowtie_mapper">bowtie_mapper</a></td><td>Takes one or two arrays of gzipped fastq files for paired or unpaired reads and an indexed genome from the bowtie_indexer applet andmaps them with bowtie</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/bowtie_mapper_parallel">bowtie_mapper_parallel</a></td><td>Wraps fastq_splitter, bowtie_mapper, and samtools_merge to align reads in parallel. Defaults to provide a worker for every 25,000,000 reads.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/bwa_indexer">bwa_indexer</a></td><td>Takes a gzipped fasta file and produces a .tar.gz file that bwa_aligner will take as input. Defaulted to use bwtsw.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/bwa_aligner">bwa_aligner</a></td><td>Takes one or two arrays of gzipped fastq files for paired or unpaired reads, and runs bwa aln followed by bwa samse/sampe.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/bwa_recalibration_pipeline">bwa_recalibration_pipeline</a></td><td>Wraps parallel_bwa, split_bam_interchromosomal_pairs, picard_mark_duplicates, gatk_realign_and_recalibrate_applet, and picard_merge_sam_files to align and recalibrate reads in parallel.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/contigset_to_fasta_gz">contigset_to_fasta_gz</a></td><td>Takes a ContigSet object (one of the DNAnexus type objects which are found in the public project "Reference Genomes"), and produces a gzipped fasta file.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/fastq_splitter">fastq_splitter</a></td><td>Takes a gzipped fastq file and splits it by number of reads into several smaller files. Useful to parallelize aligning reads.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/flexbar_read_demultiplexer">flexbar_read_demultiplexer</a></td><td>Demultiplexes indexed (barcoded) reads.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/flexbar_read_trimmer">flexbar_read_trimmer</a></td><td>Trims reads by quality score and/or position.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/gatk_apply_variant_recalibration">gatk_apply_variant_recalibration</a></td><td>Runs GATK's ApplyRecalibration on a VCF and a model produced by the gatk_variant_recalibrator applet. This recalibrates the quality of a VCF file and applies a filter to the data. Default filter level is 95% specificity.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/gatk_realign_and_recalibrate_applet">gatk_realign_and_recalibrate_applet</a></td><td>Performs indel realignment and quality recalibration on a BAM file. Requires dbSNP and known indel files which are available in the datasets directory of the Developer Applets project.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/gatk_variant_annotator_applet">gatk_variant_annotator</a></td><td>Wraps GATK module of the same name. Can be used to annotate with dbsnp, annotation or comparison VCFs, or a snpEff VCF</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/gatk_variant_caller_applet">gatk_variant_caller_applet</a></td><td>Runs GATK's UnifiedGenotyper on a BAM file to produce a VCF file.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/gatk_variant_recalibration_pipeline">gatk_variant_recalibration_pipeline</a></td><td>Wraps the gatk_variant_recalibrator and gatk_apply_variant_recalibration applets to build a model from a VCF(s) and apply it to recalibrate the quality of the VCF.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/gatk_variant_recalibrator">gatk_variant_recalibrator</a></td><td>Runs GATK's VariantRecalibrator on VCF(s) to produce model files to use in variant recalibration.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/parallel_bwa">parallel_bwa</a></td><td>Wraps fastq_splitter, bwa_aligner, and picard_merge_sam_files to align reads in parallel. Defaulted to provide a worker for every 10,000,000 reads.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/picard_calculate_hs_metrics">picard_calculate_hs_metrics</a></td><td>Calculates hybrid selection (target enrichment) metrics using the Picard CalculateHsMetrics tool.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/picard_mark_duplicates">picard_mark_duplicates</a></td><td>Runs MarkDuplicates on a BAM file. Defaulted to discard duplicate reads instead of marking them.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/picard_merge_sam_files">picard_merge_sam_files</a></td><td>Runs Picard module of same name. Useful as the reduce step of a map-reducte strategy.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/picard_sam_to_fastq">picard_sam_to_fastq</a></td><td>Wraps Picard module of the same name. Produces gzipped fastq files. Can convert either BAM or SAM files to fastq. Useful for those who receive sequence data in BAM format</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/samtools_merge">samtools_merge</a></td><td>Runs the merge module of samtools. Useful as the reduce step of a map-reduce strategy.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/samtools_view">samtools_view</a></td><td>Runs the view module of samtools. Used in the pipeline to split a BAM by region. Could also be used to extract mappings of a certain flag or convert from BAM to SAM.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/somatic_sniper">somatic_sniper</a></td><td>Runs program of same name which takes a tumor BAM file and a normal BAM file and produces variant calls. Defaulted to produce output in VCF format.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/split_bam_interchromosomal_pairs">split_bam_interchromosomal_pairs</a></td><td>This is the only applet which wraps custom written code instead of an open-source program. This splits a BAM file by intra and interchromosomal mappings. This is required by picard_mark_duplicates in a strategy which uses map-reduce to split by genome region. All interchromosmally mapped read pairs must be considered together.</td></tr>
<tr><td><a href="https://github.com/dnanexus/dnanexus-example-applets/tree/master/tumor_normal_snp_pipeline">tumor_normal_snp_pipeline</a></td><td>Wraps bwa_recalibration_pipeline and somatic_sniper to align, recalibrate, and call tumor vs. normal snps on two sets of paired reads. Defaulted to output VCF.</td></tr>
</table>
