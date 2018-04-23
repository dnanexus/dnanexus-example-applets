---
categories:
- bash
date: '2018-04-23'
github_link: https://github.com/Damien-Black/dnanexus-example-applets/tree/master/Tutorials/bash/samtools_count_git_sh
summary: Fetches and makes samtools via the GitHub repo
title: Git Dependency
type: Document
---
## What does this applet do?

This applet performs a basic SAMtools count of alignments present in an input BAM.

## Prerequisites

The app must have network access to the hostname where the git repository is located. In this example, `access.network` is set to:
```json
"access": {
  "network": ["github.com"]
}
```
Additional documentation on `access` and `network` fields can be found on the [Execution Environment Reference](https://wiki.dnanexus.com/Execution-Environment-Reference#Network-Access) wiki page.

## How is the SAMtools dependency added?

SAMtools is cloned and built from the [SAMtools GitHub](https://github.com/samtools/samtools) repository. Let's take a closer look at the `dxapp.json` file's `runSpec.execDepends` property:
```json
  "runSpec": {
 ...
    "execDepends": [
        {
        "name": "htslib",
        "package_manager": "git",
        "url": "https://github.com/samtools/htslib.git",
        "tag": "1.3.1",
        "destdir": "/home/dnanexus"
        },
        {"name": "samtools",
        "package_manager": "git",
        "url": "https://github.com/samtools/samtools.git",
        "tag": "1.3.1",
        "destdir": "/home/dnanexus",
        "build_commands": "make samtools"
        }
    ],
...
  }
```
The [`execDepends`](https://wiki.dnanexus.com/Execution-Environment-Reference?q=execDepends#Software-Packages) value is a JSON array of dependencies to resolve before the applet source code is run. In this applet, we specify the following git fetch dependency for htslib and SAMtools. Dependencies are resolved in the order they're specified. Here we **must** specify htslib first, before samtools `build_commands`, due to newer versions of SAMtools depending on htslib. An overview of the each property in the git dependency:

* `package_manager` - Details the type of dependency and how to resolve.  [supplementary details](https://wiki.dnanexus.com/Execution-Environment-Reference#Software-Packages).
* `url` - Must point to the server containing the repository. In this case, a github url.
* `tag`/`branch` - Git tag/branch to fetch.
* `destdir` - Directory on worker to which the git repo is cloned.
* `build_commands` - If needed, build commands to execute. We know our first dependency, htslib, is built when we build SAMtools; as a result, we only specify "build_commands" for the SAMtools dependency.


{% include note.html content="`build_commands` are executed from the `destdir`; use `cd` when appropriate." %}

## How is SAMtools called in our src script?

Because we set `"destdir": "/home/dnanexus"` in our `dxapp.json`, we know the git repo is cloned to the same directory from which our script will execute. Our example directory's structure:
```
├── home
│   ├── dnanexus
│       ├── < app script >
│       ├── htslib
│       ├── samtools
│           ├── < samtools binary >
```
Our samtools command from the app script is `samtools/samtools`.

## Applet Script
```bash
main() {
  set -e -x -o pipefail

  dx download "$mappings_bam"

  count_filename="${mappings_bam_prefix}.txt"
  readcount=$(samtools/samtools view -c "${mappings_bam_name}")
  echo "Total reads: ${readcount}" > "${count_filename}"

  counts_txt=$(dx upload "${count_filename}" --brief)
  dx-jobutil-add-output counts_txt "${counts_txt}" --class=file
}
```

{% include note.html content="We could've built samtools in a destination within our `$PATH` or added the binary directory to our `$PATH`. Keep this in mind for your app(let) development" %}
