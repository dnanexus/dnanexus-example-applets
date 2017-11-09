---
date: 2017-08-02
title: Basic Applet Structure
set: getting-started
set_order: 3
description: Basics of the app(let) directory
type: Document
---
Before you are able to [`dx build`](https://wiki.dnanexus.com/Command-Line-Client/Index-of-dx-Commands#build) an app(let) you must construct a valid app(let) directory. An app(let) directory skeleton looks like:
```bash
Applet directory
├── src
│   └── code.sh  # The applet script
│
├── resources/  # (optional directory)
├── test/
├── dxapp.json
├── Readme.md
└── Readme.developer.md
```

## Dxapp.json
The [dxapp.json](https://wiki.dnanexus.com/dxapp.json) file defines the interface of the applet and tells the platform how the app(let) should be run. As you may have guessed from the extension `*.json` an applet is defined in the [JSON](https://en.wikipedia.org/wiki/JSON) format. It will contain all the information necessary to build the app(let) on the platform. An example `dxapp.json` used in our SAMtools count tutorials:

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

## Applet Script

Each applet has a script, written in bash or python2.7, that executes on a worker when a platform job is created. Often, this script just calls other scripts and compiled binaries. These scripts will play a fundamental role as we discuss parallel and distributed applications.

**Note:** App(let) scripts are executed from the `$HOME` directory of a worker. Keep this in mind when performing file operations and relative directory changes.

## Resources Directory

This directory is optional.

The contents of this directory are packed and included in an applet when it is built. When an applet is executed the contents are unpacked and placed in the workers `\` root directory. Since an app(let)s' script is executed from the `HOME` directory of a worker, it is common to structure the resources directory: `resource/home/dnanexus/*`, where \* represent the contents to unpack.

## Test Folder

This folder is not packaged, nor does it have an impact on app(let) creation. It does, however, touch on the very important topics of Continuous Integration, App development, and App maintenance. These are advanced development topics outside the scope of basic DNAnexus tutorials and examples. I mention this folder only because the `dx-app-wizard` generates a template to help users write platform compatible Python functional test. Again, these topics are outside the scope of a basic tutorial or app(let) example.

## Readme Files

These are the markdown files that will be rendered and displayed to users of your app(let). In the `Readme.md`, it is common to include basic usage information such as inputs details, output details, or extended descriptions of parameters. In the `Readme.developer.md`, it is common to include advanced usage information such as instance type selection, entry point descriptions, or known bugs. You can find more information about Readmes on the wiki [App Build Process](https://wiki.dnanexus.com/Developer-Tutorials/App-Build-Process#Readme-files) page.