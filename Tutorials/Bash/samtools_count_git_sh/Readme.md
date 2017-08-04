## What does this applet do?

This applet performs a basic SAMtools count of alignments present in an input BAM.

## Prerequisites

The app must have network access to the hostname where the git repository is located. In this example `access.network` is set to:
```json
"access": {
  "network": ["github.com"]
}
```
Additional documentation can be found [here](https://wiki.dnanexus.com/Execution-Environment-Reference#Network-Access)

## How is the SAMtools dependency added?

SAMtools master branch is fetched from the [SAMtools GitHub repository](https://github.com/samtools/samtools). Let's take a closer look at the dxapp.json's `runSpec.execDepends` property:
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
        "tag": "0.1.19",
        "destdir": "/home/dnanexus",
        "build_commands": "make samtools"
        }
    ],
...
  }
```
`execDepends` value is a JSON array of dependencies to resolve before the applet src code is run. In this applet we specify the following git fetch dependency for htslib and SAMtools. Dependencies are resolved in the order they're specified. Here we **must** specify htslib first, before samtools `build_commands`, due to newer versions of SAMtools being dependent on htslib. An overview of the each property in a dependency:

* `package_manager` - Details the type of dependency and how to resolve.  [supplementary details](https://wiki.dnanexus.com/Execution-Environment-Reference#Software-Packages).
* `url` - Must point to the server containing the repository. In this case, a github url.
* `tag`/`branch` - Git tag/branch to fetch.
* `destdir` - Directory on worker to clone git repo to
* `build_commands` - If needed, build commands to execute. We know our first dependency, htslib, is built when we build SAMtools, as a result, we only specify "build_commands" for the SAMtools dependency.

<!-- INCLUDE: {% include note.html content="`build_commands` are executed from the `destdir`, use `cd` when appropriate." %} -->

## How is SAMtools called in our src script?

Because we set `"destdir": "/home/dnanexus"` in our dxapp.json we know the git repo is cloned to the same directory our script will execute from. Our examples directory structure:
```
├── home
│   ├── dnanexus
│       ├── < app script >
│       ├── htslib
│       ├── samtools
│           ├── < samtools binary >
```
Our samtools command from out app script is `samtools/samtools`.
<!-- INCLUDE: ## Applet Script -->
<!-- FUNCTION: FULL SCRIPT -->
