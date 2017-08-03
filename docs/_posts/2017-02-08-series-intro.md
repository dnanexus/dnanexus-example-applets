---
date: 2017-08-02
title: Specification Overview
description: What it means to build an App or Applet
set: applet-tutorials
set_order: 1
categories:
  - development-basics
type: Document
---

This page covers specific principles for developing applets. While we won't cover everything there is, we will cover enough to start a successful developer process. To see these principles in action review the specific tutorials that cover each concept.

{% include app-applet-forward.md %}

<!-- Cover Input/Output (I/O), Access Network, Permission and project Brief overview -->
## Entry Points

You can envision entry points as standalone functions that are executed on their own worker. They are defined in the `dxapp.json` file and referenced in your applet script:
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

Distributed computing, the primary use case of entry points, take advantage of the logical separation of work in applet scripts. Because each entry point is executed on its own worker you can alter instance type and dependencies in the `systemRequirements` field of your `dxapp.json`. Our distributed tutorials will showcase how to take advantage of the options entry points offer.

A detailed explanation of entry points can be found on the [Applets and Entry](https://wiki.dnanexus.com/API-Specification-v1.0.0/Applets-and-Entry-Points) Points wiki page.
