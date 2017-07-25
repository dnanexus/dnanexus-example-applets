---
title: DNAnexus Example Applets
keywords: getting_started
sidebar: tutorial_sidebar
permalink: index.html
toc: false
---
## Welcome

Developers new to the [DNAnexus platform](https://platform.dnanexus.com/login) may find it easier to learn by doing. This page contains a collection of simple tutorials and examples intended to showcase common tasks and methodologies when creating an app(let) on the DNAnexus platform.

The [tutorials][index] pages are not meant to show realistic everyday examples, but rather provide a strong starting point for app(let) developers. A user seeking realistic applications can find DNAnexus developed app(let)s in the [example](TODO internal link) pages. After reading through the tutorials and examples you should be able to develop app(let)s that:

- **Run efficiently**, make use of cloud computing methologies.
- **Are easy to debug**, let developers understand and resolve issues.
- **Use the scale of the cloud**, take advantage of the DNAnexus platform's.
  flexibility
- **Are easy to use**, reduce support and enabling collaboration.

## Prerequisites

We assume you already have a DNAnexus account and are looking to create app(let)s on the platform. If you need to create an account please visit our [sign up page](https://platform.dnanexus.com/register) and consider giving the [website quickstart](https://wiki.dnanexus.com/UI/Quickstart) a read. Additionally, the DNAnexus SDK, dx-toolkit, must be installed to run several of the commands used in the tutorials and examples. Follow these instructions to install dx-toolkit and begin developing app(let)s.

### Install the DNAnexus SDK

The DNAnexus SDK, dx-toolkit, provides a collection of tools to aid in the app(let) development process. While it is possible to develop and publish an app(let) using our API directly, we recommend using our SDK:
1.  Download the dx-toolkit from [our wiki downloads page](https://wiki.dnanexus.com/Downloads#DNAnexus-Platform-SDK).
2.  Set up and source the dx-toolkit environment
  * On untar/execute the [SDK version](https://wiki.dnanexus.com/Downloads#DNAnexus-Platform-SDK) installed
  *  `dx login` from the CLI after successfully sourcing
  * Enter your platform username and password.
3.  Use the `dx-app-wizard` command to create a test applet.
  *  The dx-app-wizard should walk you through a series of step that, upon completion, generate a full applet directory.
4.  Verify you're logged in to the platform from the CLI.
  *  From command line execute the `dx whoami` command. You should see your platform username.
5.  Build the sample app(let)
  * `dx build <directory dx-app-wizard created>`

If you are able to perform the steps above, you're ready to start creating app(let)s!

{% include links.html %}
