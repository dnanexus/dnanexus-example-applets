---
categories:
- Example Applet
date: '2017-08-06'
title: TMAP Genome Indexer
type: Document
---
<!-- dx-header -->
# TMAP Genome Indexer (DNAnexus Platform App)

Uses TMAP to index a FASTA file

<!-- /dx-header -->

This applet uses indexes a reference genome for use with the TMAP program.  The only required input is a reference genome in FASTA format that has been compressed with gzip (*.fa.gz or *.fasta.gz).  Optionally the user may provide a parameter string to be used with the TMAP index command.  It will be applied in the form:

tmap index -f reference.fasta $options

where options is the provided string.

The resulting outputs are named "reference.fasta.*" (where * represents the various differently named files) and bundled together in a compressed tar archive.  This archive is uploaded and output.  It can be used directly as the input for the TMAP alignment applet.


## Inputs

* **Reference Genome (gzipped FASTA file)** ``reference``: ``file``
* **Indexing Options** ``options``: ``string``

## Outputs

* **Indexed Genome Archive** ``indexed_ref``: ``file``
