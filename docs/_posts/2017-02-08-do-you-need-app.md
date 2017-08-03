---
date: 2017-08-02
title: Assessing Use Cases
categories:
  - getting-started
description: What problem are you solving, is app development the right approach?
type: Document
---

Before embarking on a new journey it's important to ask, *"Where am I going?"* and *"What am I looking for?"*. This page takes a step back and assesses alternative solutions and approaches that may get you to the real answer you want.

Look before you leap.

##  Alternatives to Applets
Before you start implementing your own app(let), you should explore the [DNAnexus App Library](https://platform.dnanexus.com/apps). Our App Library, which we continuously update, contains common bioinformatics tools implemented on the platform. The Apps are robust and generalized in their implementation and often expose many options to the user. There may be cases where you want to make a slight tweak to an existing App, for this reason, we open source our Apps. Below is an example of getting the source code for the [`BWA-MEM FASTQ Read Mapper`](https://platform.dnanexus.com/app/bwa_mem_fastq_read_mapper) using [`dx find apps`](https://wiki.dnanexus.com/Command-Line-Client/Index-of-dx-Commands#find-apps), [`dx describe`](https://wiki.dnanexus.com/Command-Line-Client/Index-of-dx-Commands#describe), and [`dx get`](https://wiki.dnanexus.com/Command-Line-Client/Index-of-dx-Commands#get):

```bash
~ $ dx find apps | grep "BWA-MEM FASTQ Read Mapper"
x BWA-MEM FASTQ Read Mapper (bwa_mem_fastq_read_mapper), v1.5.4
  # We will use the unique App name "bwa_mem_fastq_read_mapper" to get the app-id
~ $ dx describe bwa_mem_fastq_read_mapper
Result 1:
ID                  app-xxxx  # get this ID
Class               app
Billed to           org-dnanexus
Name                bwa_mem_fastq_read_mapper
...
~ $ dx get app-xxxx
Creating "./bwa_mem_fastq_read_mapper" output directory
Downloading application data
Unpacking resources
# Now we have the populated app(let) directory structure for BWA-MEM
# In the folder bwa_mem_fastq_read_mapper/
# We make necessary minor changes and rebuild
dx build bwa_mem_fastq_read_mapper
# app-xxxx
```

##  Swiss Army Knife

For running quick one-liners on a worker we recommend the swiss army knife.

For example, you might want to merge overlapping entries in a BED file
using a `bedtools` command. You could create an applet that has a BED file input
and a merged BED file output, or you could just type the command and run it using
the Swiss Army Knife app.

The Swiss Army Knife app takes an arbitrary list of files and a command string as
inputs. In commands, it’s often useful to refer to input files using wildcards, like
`*.bed`, so the command can be independent of the specific input files. It downloads all the files, runs the command, and uploads any new files. So, we could just add our BED file as an input and type in the appropriate command. From command line using [`dx run`](https://wiki.dnanexus.com/Command-Line-Client/Index-of-dx-Commands#run):

```bash
~ $ dx run swiss-army-knife
Entering interactive mode for input selection.

Input:   Command line (cmd)
Class:   string

Enter string value ('?' for more options)
cmd: bedtools merge -i *.bed > merged.bed

Select an optional parameter to set by its # (^D or <ENTER> to finish):

 [0] Input files (in)
 [1] Optional Docker image identifier (image)

Optional param #: 0

Input:   Input files (in)
Class:   array:file

Enter file values, one at a time (^D or <ENTER> to finish, <TAB> twice for compatible files in
current directory, '?' for more options)
in[0]: /Bed-Files/target.bed
# Add more files, confirm inputs, then run.
```

The Swiss Army Knife app has a number of tools built-in that you can refer to in commands. The list of tools and more information are available on the [Swiss Army Knife apps' page](https://platform.dnanexus.com/app/swiss-army-knife).

##  Cloud Workstation
For exploratory development or analysis, we recommend Cloud Workstation. The [Cloud Workstation](https://wiki.dnanexus.com/Developer-Tutorials/Cloud-Workstations) app creates a DNAnexus worker for a specified period of time that
, if you have [configured and enabled SSH](https://wiki.dnanexus.com/developer-tutorials/connecting-to-jobs), you can SSH into. Useful when you want to:
* Work interactively on a remote computer that has fast
access to files stored on the DNAnexus platform
* Work in an environment that has fast access to the internet
* Work on a machine with high processing power or high storage.

Running the Cloud Workstation app from the terminal:

**Set up SSH access with `dx ssh_config`:**

```shell
$ dx ssh_config
Select an SSH key pair to use when connecting to DNAnexus jobs. The public key will be saved to your
DNAnexus account (readable only by you). The private key will remain on this computer.                              
dx is already configured to use the SSH key pair at:
    /home/me/.dnanexus_config/ssh_id
    /home/me/.dnanexus_config/ssh_id.pub
0) Use this SSH key pair
1) Select or create another SSH key pair...                                                                                             
Pick a numbered choice: 0
Updated public key for user user-xxxx
Your account has been configured for use with SSH. Use dx run with the --allow-ssh, --ssh,
or --debug-on options to launch jobs and connect to them.
```
    
**Run the app:**
    
```shell
$ dx run app-cloud_workstation --allow-ssh                                                
Select an optional parameter to set< by its # (^D or <ENTER> to finish):        
[0] Maximum Session Length (suffixes allowed: s, m, h, d, w, M, y) (max_session_length) [default="1h"]                                 
[1] Files (fids)                                         
              
Optional param #:
Using input JSON:
{}                                                                                                                                      
Confirm running the executable with this input [Y/n]: y                                                                                 
Calling app-F40jZqQ9fZF4ggVZ1qGYQz5z with output destination project-xxxx:/                                         
Job ID: job-xxxx                                                                                                    
Watch launched job now? [Y/n] n
```
    
**SSH into the job running the app**
    
```shell
$ dx ssh job-xxxx
```

The Cloud Workstation app takes a rather useful approach for exploratory development on the platform. Feel free to look at the app code for Cloud Workstation v1.0.3, which you can obtain with `dx get app-F40jZqQ9fZF4ggVZ1qGYQz5z`.
