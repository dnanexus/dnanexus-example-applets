---
categories:
- Example Applet
date: '2017-08-08'
title: 'TMAP: Torrent Mapping Alignment Program'
type: Document
---
<!-- dx-header -->
# TMAP: Torrent Mapping Alignment Program (DNAnexus Platform App)

A parallelized implementation of TMAP

<!-- /dx-header -->

This applet combines the SFF splitter applet and the TMAP (Torrent Mapping Alignment Program) algorithm to map reads coming from the Ion Torrent series of instruments.  The reads must be in SFF format (and compressed with gzip - *.sff.gz).  A properly indexed reference genome can be generated with the "TMAP Genome Indexer" applet.  

Optionally, the user may specify additional parameters to be used on when executing TMAP.  The by default this string is "map1 -i sff", instructing TMAP to use mapping algorithm 1 and that the input is sff.  More documentation about for TMAP can be found [here](https://github.com/iontorrent/TMAP).

## Inputs

* **Reads (gzipped SFF format)** ``reads``: ``file``
* **Reference Indexed for TMAP** ``indexed_ref``: ``file``
* **Reads Per Mapping Job** ``chunk_size``: ``int``
* **Mapping Parameters** ``map_params``:``string``

## Outputs

* **Mapped Reads (sorted BAM file)** ``mappings``: ``file``
