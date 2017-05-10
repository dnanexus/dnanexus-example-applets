---
title: Basic Applet Structure
sidebar: tutorial_sidebar
permalink: basic_applet_structure.html
summary: "This page will cover just the basics of what an app(let) is. For a deeper understanding refer to the DNAnexus wiki."
toc: false
---
##  Components of an Applet

An app(let) directory skeleton looks like:
```
├── Applet dir
│   ├── src/
│   ├── dxapp.json
│   ├── resources/ (optional)
```

### dxapp.json
The [dxapp.json](https://wiki.dnanexus.com/dxapp.json) file defines the interface of the applet and tells the platform how the app(let) should be run. It will contain all the information necessary to build the app(let) on the platform. An example dxapp.json used in our SAMtools count tutorials:

```json
{
  "name": "samtools_count",
  "title": "SAMtools count",
  "summary": "Basic SAMtools count implementation",
  "dxapi": "1.0.0",
  "version": "0.0.1",
  "inputSpec": [
    {
      "name": "mappings_bam",
      "label": "Mapping",
      "class": "file",
      "patterns": ["*.bam"],
      "help": "BAM format file."
    }
  ],
  "outputSpec": [
    {
      "name": "counts_txt",
      "class": "file",
      "label": "Read count file",
      "patterns": [
        "*.txt"
      ],
      "help": "Output file with Total reads as the first line."
    }
  ],
...
```

### Applet Script

Each applet has a script, written in bash or python2.7, that executes on a worker when a platform job is created. Often, this script just calls other scripts and compiled code.

{% include note.html content="Note that app(let) scripts are executed from the `HOME` directory of a worker." %}

### Resources directory

This directory is optional.

The contents of this directory are packed and included in an applet when it is built. When an applet is executed the contents are unpacked and placed in the workers `\` root directory.
{% include note.html content="Since an app(let)'s' script is executed from the `HOME` directory of a worker, it is common to structure the resources directory: `resource/home/dnanexus/*`, where \* represent the contents to unpack." %}
