---
date: 2017-08-02
title: Development Basics
description: What it means to build an App or Applet
set: getting-started
set_order: 4
categories:
  - getting-started
type: Document
---

This page covers essential principles for developing applets. While we won't cover everything there is, we will provide a good jump-off point. To see these principles in action review the specific tutorials that cover each concept.

If you'd like a more in-depth tutorial, see the [Advanced App Tutorials](https://wiki.dnanexus.com/Developer-Tutorials/Advanced-App-Tutorial) page in the wiki.

{% include app-applet-forward.md %}

<!-- Cover Input/Output (I/O), Access Network, Permission and project Brief overview -->
## Inputs and Outputs

You can envision an app(let) as a [black box](https://en.wikipedia.org/wiki/Black_box) that takes an input and gives you an output. Throughout our documentation, you may see the term `I/O`, this just means "Input/Output". We describe the input and output specification in the `dxapp.json` file's `inputSpec` and `outputSpec` respectively:
```json
{
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
      "help": "Output file with total counts as the first line."
    }
  ]
}
```
We offer a great deal of customization with I/O fields. The full details of what properties and values are available in I/O specifications can be found on our [Job Input and Output](https://wiki.dnanexus.com/API-Specification-v1.0.0/Job-Input-and-Output) wiki page. When viewing tutorials and examples, build the app(let) for yourself and see how the I/O specification is presented on the platform.

## Dependencies

Before your app(let) executes there will certain preconditions you must resolve. One common dependency you'll see in our tutorials is [SAMtools](http://www.htslib.org/doc/samtools.html). To resolve this dependency we often use the `execDepends` property in the `dxapp.json`.
```json
{
    "execDepends": [
      {"name": "samtools"}
    ]
}
```
As you read through each tutorial pay close attention to how we resolve the SAMtools dependency. There are various approaches you can take; each better suited for specific use cases.

## Entry Points

Entry points are how we distribute work across many workers.

You can envision entry points as standalone functions that are executed on their own worker, these sub-executions are what we call sub-jobs. They are defined in the `dxapp.json` file and referenced in your applet script:
* **dxapp.json** we can set specs for each entry point
```json
{
  "runSpec": {
    "systemRequirements": {
      "main": {
        "instanceType": "mem1_ssd1_x4"
      },
      "count_func": {
        "instanceType": "mem1_ssd1_x2"
      },
      "sum_reads": {
        "instanceType": "mem1_ssd1_x4"
      }
    }
}
}
```
* **Python** scripts refer to them as decorators on functions
```python
@dxpy.entry_point('main')
def main(**job_inputs):
    # do work
```
* **bash** scripts refer to them as the function itself
```bash
main() {
    # do work
}
```
Distributed computing, the primary use case of entry points, takes advantage of the logical separation of work in app(let) scripts. Because each entry point is executed on its own worker you are able to alter instance type and dependencies in the `systemRequirements` field of your `dxapp.json`. Our distributed tutorials will showcase how to take advantage of the options entry points offer.

A detailed explanation of entry points can be found on the [Applets and Entry](https://wiki.dnanexus.com/API-Specification-v1.0.0/Applets-and-Entry-Points) Points wiki page.