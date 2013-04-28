# BEDTools applets

This directory contains a script, `generate-bedtools-applets`, which automatically generates a DNAnexus file-based applet for each [BEDTools](https://code.google.com/p/bedtools/) subcommand. It does so by parsing the help printed by each tool into an applet spec.

#### Generating the applets

The script has one dependency other than the [SDK](http://wiki.dnanexus.com/Downloads#DNAnexus-Platform-SDK): [`jinja2`](http://jinja.pocoo.org/docs/). Install it with `sudo pip install jinja2` or equivalent. Then, generate the applets by running `make` in this directory.

#### Building the applets
For each applet that you want to build, run `dx-build-applet` on its directory, e.g. `dx-build-applet intersect`.
