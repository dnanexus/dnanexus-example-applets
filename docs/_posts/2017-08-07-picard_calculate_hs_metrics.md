---
categories:
- Example Applet
date: '2017-08-07'
title: Picard CalculateHsMetrics
type: Document
---
# Picard CalculateHsMetrics

This applet is a simple wrapper around Picard CalculateHsMetrics, a tool that calculates hybrid selection
metrics. It is suitable for evaluating the performance of hybrid selection (target enrichment) protocols.

Most hybrid selection (target enrichment) protocols start with some *targets* that need to be selected for;
usually these correspond to exons of genes of interest. The protocols also require generation of some *baits*,
which are usually larger than the targets, and which are used to "capture" the sequences of interest.

Picard CalculateHsMetrics uses a list of bait/target coordinates, the reference genome sequence, and a
BAM file containing alignments, to calculate metrics such how well the baits or targets were covered
as well as any encountered AT/GC biases.

General information about Picard CalculateHsMetrics can be found at this link:

http://picard.sourceforge.net/command-line-overview.shtml#CalculateHsMetrics

Information about the output format of the metrics file can be found at this link:

http://picard.sourceforge.net/picard-metric-definitions.shtml#HsMetrics

## Inputs

* **Baits interval file** (`baits`) and **Targets interval file** (`targets`): Files containing the bait and target
coordinates in Picard Interval format. Both files are required; if you are missing one, you can use the other for both
inputs. The Picard Interval format is a special format which requires a SAM header describing the sequences (names
and lengths) of the reference genome, followed by rows describing the interval coordinates. Each row must be tab-delimited
and must contain five fields: the chromosome name, the low and high coordinates of the interval (inclusive, 1-indexed),
the strand, and the name of the interval. Sample target interval files are provided for the Illumina TruSeq target
enrichment protocol in the "example_targets" directory of the applet source archive.

* **BAM file** (`bam`): The BAM file containing the alignments on which metrics will be calculated.

* **Reference gzipped fasta file** (`reference_fastagz`): A gzipped fasta file with the reference genome sequence.

## Outputs

* **Hybrid selection metrics** (`hsmetrics_file`): A text file containing the hybrid selection metrics, as explained in http://picard.sourceforge.net/picard-metric-definitions.shtml#HsMetrics.

* **Per-target metrics** (`pertarget_hsmetrics_file`): A text file containing metrics calculated per target.

