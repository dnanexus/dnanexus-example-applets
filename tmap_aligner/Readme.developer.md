# tmap_scatter_gather Developer Readme

This app has the following entry points:

* main

  Downloads the SFF file and calls the "sff_splitter" applet to generate smaller files for parallel processing

* map

  For each output of the splitter applet a process job is launched that takes that fraction of the reads and the indexed genome

* process

  TMAP program is called here to align a chunk of the reads.  Output is a BAM file.

* postprocess

  All BAM files from the split jobs are downloaded and merged using Samtools merge.  This is the final output to the applet.

