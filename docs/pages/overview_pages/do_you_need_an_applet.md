---
title: Do You Need an Applet?
keywords: getting_started
sidebar: tutorial_sidebar
permalink: do_you_need_an_applet.html
toc: false
---

##  Alternatives to Applets
Before you start implementing your own app(let), you should explore the [DNAnexus App Library](https://platform.dnanexus.com/apps). Our App Library, which we continuously update, contains common bioinformatics tools implemented on the platform. For more exploratory or speedy processing we recommend the swiss army knife or cloud workstation apps.

##  Swiss Army Knife

Useful for quickly processing files with a couple commands and getting the
output. For example, you might want to merge overlapping entries in a BED file
using a `bedtools` command. You could create an applet that has a BED file input
and a merged BED file output, or you could just type the command and run it using
the Swiss Army Knife app.

The Swiss Army Knife app takes an arbitrary list of files and a command string as
inputs. It downloads all the files, runs the command, and uploads any new files.
So, we could just add our BED file as an input and type in the appropriate command:

When using Swiss Army Knife, it’s often useful to refer to input files using wildcards, like
`*.bed`, so the command can be independent of the specific input files.

The Swiss Army Knife app has a number of tools built-in that you can refer to them
in commands. The list of tools and more information are available on the
[Swiss Army Knife's apps page](https://platform.dnanexus.com/app/swiss-army-knife).

##  Cloud Workstation

Useful when you want to work interactively on a remote computer that has fast
access to files stored on the DNAnexus platform as well as the rest of the internet
 or when you want to work on a computer with specific resources, like 32 cores or
200 GiB of memory. The [Cloud Workstation](https://wiki.dnanexus.com/Developer-Tutorials/Cloud-Workstations) app creates a DNAnexus worker for a specified period of time that
, if you have [configured and enabled SSH](https://wiki.dnanexus.com/developer-tutorials/connecting-to-jobs), you can SSH into.

Running the Cloud Workstation app from the terminal:

1.  Set up SSH access with `dx ssh_config`

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
    
2.  Then run the app:
    
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
    
3.  SSH into the job running the app
    
    ```shell
    $ dx ssh job-xxxx
    ```

The Cloud Workstation app takes a rather useful approach for exploratory development on the platform. Feel free to look at the app code for Cloud Workstation v1.0.3, which you can obtain with a `dx get app-F40jZqQ9fZF4ggVZ1qGYQz5z`.

