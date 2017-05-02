---
title: Basic Applet Structure
sidebar: tutorial_sidebar
permalink: basic_applet_structure.html
layout: default
---

##  Components of an Applet

There are three main parts of an applet on DNAnexus:

- **[dxapp.json file](https://wiki.dnanexus.com/dxapp.json)**

The dxapp.json file defines the interface of the applet, and it tells the platform
some details about how the applet should be run.

- **Applet script**

Each applet has a script, written in bash or python2.7, that the platform runs
within the worker when the applet it run. Often, this script just calls other scripts
and compiled code.

- **Applet resources**

Optionally, each applet can have resource files packaged with it. These will be
placed on the worker before the applet script is run. There are a few ways to package
applet resources, but for now we will just work with the `./resources` directory
of the applet.

### dxapp.json

Based on the command above, our `dxapp.json` file is straighforward. The applet will
have one input for the HISAT2 reference tarball and an input for each FASTQ file.
Its output will be a SAM file produced by histat2. The rest of the `dxapp.json`
file is mostly some names and descriptions of the applet:

```json
{
  "name": "hisat2",
  "title": "HISAT2",
  "summary": "Runs a simple HISAT2 command.",
  "inputSpec": [
    {
      "name": "hisat2_index_targz",
      "label": "HISAT2 Index Tarball",
      "class": "file",
      "help": "Tar.gz'd HISAT2 index for a reference genome, produced using hisat2-build."
    },
    {
      "name": "mate1_fastq",
      "label": "Mate 1 FASTQ",
      "class": "file",
      "help": "FASTQ file containing mate 1 reads."
    },
    {
      "name": "mate2_fastq",
      "label": "Mate 2 FASTQ",
      "class": "file",
      "help": "FASTQ file containing mate 2 reads."
    }
  ],
  "outputSpec": [
    {
      "name": "aligned_sam",
      "label": "Aligned SAM",
      "class": "file",
      "help": "SAM file with alignments reported by HISAT2"
    }
  ],
```

The remainder of the `dxapp.json` depends on whether you want to write
your applet script in python or bash:

<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#python" data-toggle="tab">python</a></li>
    <li><a href="#bash" data-toggle="tab">bash</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="python">

<pre class="highlight"><code><span class="w">  </span><span class="s2">"runSpec"</span><span class="err">:</span><span class="w"> </span><span class="p">{</span><span class="w">
    </span><span class="nt">"file"</span><span class="p">:</span><span class="w"> </span><span class="s2">"src/script.py"</span><span class="p">,</span><span class="w">
    </span><span class="nt">"interpreter"</span><span class="p">:</span><span class="w"> </span><span class="s2">"python2.7"</span><span class="w">
  </span><span class="p">}</span><span class="w">
</span><span class="err">}</span><span class="w">
</span></code></pre>

  </div>
  <div role="tabpanel" class="tab-pane" id="bash">

 <pre class="highlight"><code><span class="w">  </span><span class="s2">"runSpec"</span><span class="err">:</span><span class="w"> </span><span class="p">{</span><span class="w">
    </span><span class="nt">"file"</span><span class="p">:</span><span class="w"> </span><span class="s2">"src/script.sh"</span><span class="p">,</span><span class="w">
    </span><span class="nt">"interpreter"</span><span class="p">:</span><span class="w"> </span><span class="s2">"bash"</span><span class="w">
  </span><span class="p">}</span><span class="w">
</span><span class="err">}</span><span class="w">
</span></code></pre> 
  
  </div>
</div>

### Applet Script

We already know the command that we need to run for HISAT2, but in order to run that
command on DNAnexus, we have to take care of a few other things. First, we need to
download the input files to the local storage of the worker running the applet. Then,
we run the hisat2 command. Then, we upload the output SAM file and indicate that it
is the output associated with the `aligned_sam` entry in the applet's outputSpec.

We will start with the command we want to run, and then add in the DNAnexus-specific
parts:

```shell
#!/bin/bash

main() {

    hisat2 --dta -x $index_basename -1 mate1.fastq -2 mate2.fastq -S hisat2_output.sam

}
```

Now, we will add some commands to download the input files from the platform to the
local storage of the worker:


```shell
#!/bin/bash

main() {
    # Download the inputs to the worker's local storage
    dx download "$hisat2_index_targz" -o hisat2_index.tar.gz
    dx download "$mate1_fastq" -o mate1.fastq
    dx download "$mate2_fastq" -o mate2.fastq

    hisat2 --dta -x $index_basename -1 mate1.fastq -2 mate2.fastq -S hisat2_output.sam
}
```

Next, we need to extract the HISAT2 reference tarball:

```shell
#!/bin/bash

main() {
    # Download the inputs to the worker's local storage
    dx download "$hisat2_index_targz" -o hisat2_index.tar.gz
    dx download "$mate1_fastq" -o mate1.fastq
    dx download "$mate2_fastq" -o mate2.fastq

    # Extract the tarball containing the HISAT2 reference
    tar xf hisat2_index.tar.gz

    hisat2 --dta -x $index_basename -1 mate1.fastq -2 mate2.fastq -S hisat2_output.sam
}
```

Then, we will upload the output file, and associate with the correct field
from the outputSpec of the `dxapp.json`.

```shell
#!/bin/bash

main() {
    # Download the inputs to the worker's local storage
    dx download "$hisat2_index_targz" -o hisat2_index.tar.gz
    dx download "$mate1_fastq" -o mate1.fastq
    dx download "$mate2_fastq" -o mate2.fastq

    # Extract the tarball containing the HISAT2 reference
    tar xf hisat2_index.tar.gz

    hisat2 --dta -x $index_basename -1 mate1.fastq -2 mate2.fastq -S hisat2_output.sam

    # Upload the resulting file and associate it with the aligned_sam output
    uploaded_id=$(dx upload hisat2_output.sam --brief)
    dx-jobutil-add-output aligned_sam $uploaded_id
}
```

Finally, we need to figure out what string to give to the `hisat2` executable for
"index basename". The HISAT2 website has reference tarballs ready for download, so
we can look at one to see what is inside:

```shell
$ tar tf grch38.tar.gz 
grch38/
grch38/genome.5.ht2
grch38/genome.2.ht2
grch38/make_grch38.sh
grch38/genome.3.ht2
grch38/genome.4.ht2
grch38/genome.7.ht2
grch38/genome.1.ht2
grch38/genome.6.ht2
grch38/genome.8.ht2
```

In this case, the index basename should be "grch38/genome", so we add some string manipulation
logic to our applet to set that variable correctly.

```shell
#!/bin/bash

main() {

    # Download the inputs to the worker's local storage
    dx download "$hisat2_index_targz" -o hisat2_index.tar.gz
    dx download "$mate1_fastq" -o mate1.fastq
    dx download "$mate2_fastq" -o mate2.fastq

    # Extract the tarball containing the HISAT2 reference
    tar xf hisat2_index.tar.gz

    # Get the index basename to use when calling hisat2
    index_filename=$(dx describe --name "$hisat2_index_targz")
    index_basename=${index_filename%.tar.gz}/genome

    # Run hisat2
    hisat2 --dta -x $index_basename -1 mate1.fastq -2 mate2.fastq -S hisat2_output.sam

    # Upload the resulting file and associate it with the aligned_sam output
    uploaded_id=$(dx upload hisat2_output.sam --brief)
    dx-jobutil-add-output aligned_sam $uploaded_id
}
```

And that is it. This full script is at `applets/hisat2/bash/41_script.sh` 

The python version is similar, but obviously with python versions of each command
used:

```python
import os
import subprocess

import dxpy

@dxpy.entry_point("main")
def main(hisat2_index_targz, mate1_fastq, mate2_fastq):

    # First, download all the input files to local storage
    dxpy.download_dxfile(hisat2_index_targz, "hisat2_index.tar.gz")
    dxpy.download_dxfile(mate1_fastq, "mate1.fastq")
    dxpy.download_dxfile(mate2_fastq, "mate2.fastq")


    # Second, extract the index tarball
    proc = subprocess.Popen(["tar", "xf", "hisat2_index.tar.gz"])
    proc.wait()

    # Third, figure out what the basename of the hisat2 reference index is
    # This depends on the basename following the pattern used in the indexes
    # distributed by the authors of HISAT2, that is the index in grch37.tar.gz
    # will extract to grch37/genome*
    index_basename = os.path.join(dxpy.DXFile(hisat2_index_targz).name[:-len(".tar.gz")],
                                  "genome")

    # Prepare the hisat2 command and run it.
    hisat2_cmd_template = ("hisat2 --dta -x {index_basename} -1 {mate1_fastq} "
                           "-2 {mate2_fastq} -S {hisat2_output_sam}")
    hisat2_cmd = hisat2_cmd_template.format(
        index_basename=index_basename,
        mate1_fastq="mate1.fastq",
        mate2_fastq="mate2.fastq",
        hisat2_output_sam="hisat2_output.sam")
    subprocess.check_call(hisat2_cmd, shell=True)

    # Upload the output SAM file.
    uploaded_dxfile = dxpy.upload_local_file("hisat2_output.sam")

    # Return the ID of the uploaded SAM file associated with the "aligned_sam"
    # field in the outputSpec in dxapp.json.
    return {"aligned_sam": dxpy.dxlink(uploaded_dxfile.get_id())}
```

##  Building the Applet

Now we can actually build this applet to the DNAnexus platform. This requires creating
a directory structure for the applet and running `dx build`. The directory should look
like this:

```shell
hisat2/dxapp.json
hisat2/src/script.py
hisat2/resources/usr/bin/hisat2
hisat2/resources/usr/bin/hisat2-build-s
hisat2/resources/usr/bin/hisat2-build-l-debug
...and more hisat2 executables
```

We can create this with the following commands:

```shell
$ mkdir hisat2
$ mkdir hisat2/src
$ mkdir -p hisat2/resources/usr/bin
$ cp applets/hisat2/bash/41_script.sh src/script.sh
$ cp applets/hisat2/bash/41_dxapp.json dxapp.json
$ cp resources/hisat2/hisat2-2.0.4/hisat2* resources/usr/bin
$ cp resources/hisat2/hisat2-2.0.4/extract* resources/usr/bin
$ dx build hisat2
```

```shell
$ mkdir hisat2
$ mkdir hisat2/src
$ mkdir -p hisat2/resources/usr/bin
$ cp applets/hisat2/python/41_script.py src/script.py
$ cp applets/hisat2/python/41_dxapp.json dxapp.json
$ cp resources/hisat2/hisat2-2.0.4/hisat2* resources/usr/bin
$ cp resources/hisat2/hisat2-2.0.4/extract* resources/usr/bin
$ dx build hisat2
```

<ul id="profileTabs" class="nav nav-tabs">
    <li class="active"><a href="#python2" data-toggle="tab">python</a></li>
    <li><a href="#bash2" data-toggle="tab">bash</a></li>
</ul>
<div class="tab-content">
  <div role="tabpanel" class="tab-pane active" id="python2">

<pre class="highlight"><code><span class="gp">$ </span>mkdir hisat2
<span class="gp">$ </span>mkdir hisat2/src
<span class="gp">$ </span>mkdir -p hisat2/resources/usr/bin
<span class="gp">$ </span>cp applets/hisat2/python/41_script.py src/script.py
<span class="gp">$ </span>cp applets/hisat2/python/41_dxapp.json dxapp.json
<span class="gp">$ </span>cp resources/hisat2/hisat2-2.0.4/hisat2<span class="k">*</span> resources/usr/bin
<span class="gp">$ </span>cp resources/hisat2/hisat2-2.0.4/extract<span class="k">*</span> resources/usr/bin
<span class="gp">$ </span>dx build hisat2
</code></pre>

  </div>
  <div role="tabpanel" class="tab-pane" id="bash2">

<pre class="highlight"><code><span class="gp">$ </span>mkdir hisat2
<span class="gp">$ </span>mkdir hisat2/src
<span class="gp">$ </span>mkdir -p hisat2/resources/usr/bin
<span class="gp">$ </span>cp applets/hisat2/bash/41_script.sh src/script.sh
<span class="gp">$ </span>cp applets/hisat2/bash/41_dxapp.json dxapp.json
<span class="gp">$ </span>cp resources/hisat2/hisat2-2.0.4/hisat2<span class="k">*</span> resources/usr/bin
<span class="gp">$ </span>cp resources/hisat2/hisat2-2.0.4/extract<span class="k">*</span> resources/usr/bin
<span class="gp">$ </span>dx build hisat2
</code></pre>
  
  </div>
</div>


And we can run it through the command line or the web UI.
