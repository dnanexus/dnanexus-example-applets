---
categories:
- Example Applet
date: '2017-08-08'
title: GATK Variant Recalibration Pipeline
type: Document
---
"The GATK VariantRecalibrator constructs a model which attempts to estimate the properties of variant calls that are true as opposed to those that are false positives. It produces model files that can then be applied to a VCF files to filter variant calls by desired specificity. The resources files  (known, training, and truth) provide example sites used in model construction. Put all files that should have label "Known" into the known input array, all files that should have label "Training" into the training array, and all "Truth" into the truth array. If a file should be labelled with more than one (e.g. training and truth), put the file into both arrays, it will only be extracted once and with the appropriate combinations of labels."