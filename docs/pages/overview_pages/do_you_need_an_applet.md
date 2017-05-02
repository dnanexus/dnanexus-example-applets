---
title: Do You Need an Applet?
keywords: getting_started
sidebar: tutorial_sidebar
permalink: do_you_need_an_applet.html
layout: default
---

##  Alternatives to Applets

Before we start implementing our HISAT2 applet, we should recall that you can
run tools on DNAnexus without creating an applet at all. This is helpful when
you are performing a one-off task or just want to do some exploration. These
approaches give up some reproducibility and auditability, but they are quick
and easy to use.

##  Swiss Army Knife

In some cases, you want to run a command or two on some files and get the
output. For example, you might want to merge overlapping entries in a BED file
using a `bedtools` command. You could create an applet that has a BED file input
and a merged BED file output, or you could just type the command and run it using
the Swiss Army Knife app.

The Swiss Army Knife app takes an arbitrary list of files and a command string as
inputs. It downloads all the files, runs the command, and uploads any new files.
So, we could just add our BED file as an input and type in the appropriate command:

When using Swiss Army Knife, it’s often useful to refer to input files using wildcards, like
`*.bed`, so the command can be independent of the specific input files. It also
allows you to refer to multiple files at once.

The Swiss Army Knife app has a number of tools built-in, so you can refer to them
in commands. The list of tools and more information are available on the
[apps page](https://platform.dnanexus.com/app/swiss-army-knife).

##  Cloud Workstation

In other cases, you want to work interactively on a remote computer that has fast
access to files stored on the DNAnexus platform as well as the rest of the internet.
You also may want to work on a computer with specific resources, like 32 cores or
200 GiB of memory. This can be accomplished with the Cloud Workstation app. The
Cloud Workstation app simply runs a DNAnexus worker for a specified period of time
that you can SSH into.

The Cloud Workstation app is described in the
[DNAnexus Wiki](https://wiki.dnanexus.com/Developer-Tutorials/Cloud-Workstations).

You can run the Cloud Workstation app from the terminal. First you need to set up
SSH access with `dx ssh_config`.

```shell
$ dx ssh_config
Select an SSH key pair to use when connecting to DNAnexus jobs. The public key will be saved to your
DNAnexus account (readable only by you). The private key will remain on this computer.                          
dx is already configured to use the SSH key pair at:
    /home/vagrant/.dnanexus_config/ssh_id
    /home/vagrant/.dnanexus_config/ssh_id.pub
0) Use this SSH key pair
1) Select or create another SSH key pair...                                                                                         
Pick a numbered choice: 0
Updated public key for user user-mkinsella_ugm
Your account has been configured for use with SSH. Use dx run with the --allow-ssh, --ssh,
or --debug-on options to launch jobs and connect to them.
```

Then run the app:

```shell
$ dx run app-cloud_workstation --allow-ssh                                                
Select an optional parameter to set< by its # (^D or <ENTER> to finish):        
 [0] Maximum Session Length (suffixes allowed: s, m, h, d, w, M, y) (max_session_length) [default="1h"]                             
 [1] Files (fids)                                         
              
Optional param #:
Using input JSON:
{}                                                                                                                                  
Confirm running the executable with this input [Y/n]: y                                                                             
Calling app-Bpx83Y00zV01pZ0zGGpyk0y0 with output destination project-BzYQ2Y005vP8j0j081y1pjqQ:/                                     
Job ID: job-BzYVf3j05vP5bQxz15q4g0pg                                                                                                
Watch launched job now? [Y/n] n
```

And ssh into the job running the app

```shell
$ dx ssh job-BzYVf3j05vP5bQxz15q4g0pg
```
