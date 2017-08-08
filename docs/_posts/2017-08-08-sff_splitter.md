---
categories:
- Example Applet
date: '2017-08-08'
title: SFF Splitter
type: Document
---

# SFF Splitter (DNAnexus Platform App)

This applet uses the [Biopython package](https://github.com/biopython/biopython) to split an SFF file into smaller SFF files.  This enables parallel processing of those reads spread across multiple workers.  The default number of reads per output file is 20,000,000.

## Inputs

* **Reads to Split (SFF file)** ``reads``: ``file``
* **Reads Per Output File** ``chunk_size``: ``int``

## Outputs

* **Array of Split Reads Files (SFF)** ``split_reads``: ``array:file``
