---
categories:
- bash
date: '2017-08-08'
github_link: https://github.com/Damien-Black/dnanexus-example-applets/tree/master/Tutorials/bash/samtools_count_binary_sh
title: Precompiled Binary
type: Document
---
This tutorial showcases packaging a precompiled binary in the `resources/` directory of an app(let).

## Precompiling a Binary
In this applet, the SAMtools binary was precompiled on an Ubuntu 14.04 machine. A user can do this compilation on an Ubuntu 14.04 machine of their own, or can utilize the Cloud Workstation app to build and compile a binary. On the Cloud Workstation, the user has full access to download the SAMtools source code and compile in the worker environment, ensuring the binary will run on future workers.

See [Cloud Workstation](https://wiki.dnanexus.com/Developer-Tutorials/Cloud-Workstations) in the App library for more information.
## Resources Directory
The SAMtools precompiled binary is placed in the `<Applet dir>/resources/` directory. Any files found in the `resources/` directory will be packaged, uploaded to the platform, and then unpackaged in the root directory `\` of the worker. In our case the `resources/` dir is:
```
├── Applet dir
│   ├── src
│   ├── dxapp.json
│   ├── resources
│       ├── usr
│           ├── bin
│               ├── < samtools binary >
```
When this applet is run on a worker the `resources/` directory will be placed in the worker's root `/`:
```
/
├── usr
│   ├── bin
│       ├── < samtools binary >
├── home
│   ├── dnanexus
│   	├── applet script
```
We are able to access the SAMtools command because the respective binary is visible from the default `$PATH` variable. The directory`/usr/bin/` is part of the `$PATH` variable, so in our script we reference the samtools command directly:
```bash
samtools view -c "${mappings_bam_name}" > "${mappings_bam_prefix}.txt"
```


See [The Resources/ Directory and its Use](https://wiki.dnanexus.com/Developer-Tutorials/App-Build-Process#The-resources/-directory-and-its-use) wiki page for more information.
