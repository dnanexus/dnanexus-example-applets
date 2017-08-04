<!-- dx-header -->
# Circos (DNAnexus Platform App)

Generates Circos (circular layout) plots

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
http://wiki.dnanexus.com/.
<!-- /dx-header -->

Creates visualizations based on genomic (or other) data arranged in a circular layout. See http://circos.ca

## How does this differ from the Circos command-line tool?

The DNAnexus Circos app is largely the same as the Circos command-line 
tool. However, because of the way that apps get executed on DNAnexus, 
some tailoring of your configuration file will probably be required.

To start with, all input data and configuration files get copied into the same 
directory. This means that whenever you reference data or another configuration
file (that isn't part of the basic Circos distribution--those remain 
unchanged in location), you'll just enter the filename without any subfolders.
This also means that if you're executing any tutorial files from the 
[Circos tutorials pack](http://circos.ca/software/download/tutorials/),
you'll need to edit the main configuration file to remove paths to data files.

Also, you won't have access to many of the command-line flags detailed in the 
[Circos Parameters](http://circos.ca/documentation/tutorials/configuration/runtime_parameters/),
so you'll have to make sure that those are all included in your
configuration file when you upload it to DNAnexus.

Finally, the DNAnexus Circos app will by default output both SVG and PNG files
for you to download. Additionally, it will also output a report with the 
plot included for you to view on DNAnexus.

## Inputs

* **Circos configuration** ``circos``: ``file``
* **Extra files** ``extra``: ``array:file`` (optional)

### Advanced

* **Output Report Name** ``report_name``: ``string`` (optional, default "Circos Report")
* **Output Base Name** ``base_name``: ``string`` (optional, default "circos")
* **SVG Report** ``use_svg``: ``boolean`` (optional, default true)
* **Output PNG** ``png``: ``boolean`` (optional, default true)
* **Output SVG** ``svg``: ``boolean`` (optional, default true)

### Circos options

* **Show ticks** ``ticks``: ``boolean`` (optional, default true)
* **Show tick labels** ``ticklabels``: ``boolean`` (optional, default true)

## Outputs

* **Report** ``report``: ``record``
* **Image (SVG)** ``svg``: ``file`` (optional)
* **Image (PNG)** ``png``: ``file`` (optional)
* **HTML Map** ``html``: ``file`` (optional)

