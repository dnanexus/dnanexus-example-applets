Developers new to the [DNAnexus platform](https://platform.dnanexus.com/login) may find it easier to learn by doing. This page contains a collection of simple tutorials and examples intended to showcase common tasks and methodologies when creating an app(let) on the DNAnexus platform.

The tutorial pages are not meant to show everyday examples, but rather provide a strong starting point for app(let) developers. A user seeking realistic examples can examine DNAnexus-developed app(let)s in the Examples pages or the [Developer Applets](https://platform.dnanexus.com/projects/B406G0x2fz2B3GVk65200003/data/) resource project (platform login required). After reading through the tutorials and examples, you should be able to develop app(let)s that:

- **Run efficiently**: make use of cloud computing methodologies
- **Are easy to debug**: let developers understand and resolve issues
- **Use the scale of the cloud**: take advantage of the DNAnexus platform's flexibility
- **Are easy to use**: reduce support and enable collaboration

## Prerequisites

We assume that you already have a DNAnexus account and are looking to create app(let)s on the platform. If you need to create an account, please visit our [sign up page](https://platform.dnanexus.com/register) and look over the [website quickstart](https://wiki.dnanexus.com/UI/Quickstart). Additionally, the DNAnexus SDK, `dx-toolkit`, must be installed in order to run commands used in the tutorials and examples.

## Install the DNAnexus SDK

The DNAnexus SDK, `dx-toolkit`, provides a collection of tools to aid in the app(let) development process. While it is possible to develop and publish an app(let) using our API directly, we recommend using our SDK:

### Download `dx-toolkit`

`dx-toolkit` can be found on our [Downloads](https://wiki.dnanexus.com/Downloads#DNAnexus-Platform-SDK) page. Alternatively, it can be built from the [github source code](https://github.com/dnanexus/dx-toolkit).

```bash
$ curl -O https://wiki.dnanexus.com/images/files/dx-toolkit-v0.228.0-osx.tar.gz
$ tar -xvf dx-toolkit-v0.228.0-osx.tar.gz
  # A folder dx-toolkit/ will be created in the current directory
```

## Source the `dx-toolkit` Environment

```bash
$ source dx-toolkit/environment
  # You are now ready to use dx commands
$ dx --version
  # Downloaded dx-toolkit version outputs
```

## Log In to the Platform

From the CLI, use the command [`dx login`](https://wiki.dnanexus.com/Command-Line-Client/Index-of-dx-Commands#login). You may wish to log in using a [token](https://wiki.dnanexus.com/Command-Line-Client/Login-and-Logout#Authentication-Tokens).

```bash
$ dx login
  Username:  # Platform username
  Password:  # Platform password
  # If this is your first time logging in, you will be prompted to select a project
  Note: Use dx select --level VIEW or dx select --public to select from projects for
  which you only have VIEW permissions.

  Available projects (CONTRIBUTE or higher):
  0) DNAnexus Example (CONTRIBUTE)
  1) My Research Project (ADMINISTER)
  
  Pick a numbered choice or "m" for more options: 1
  Setting current project to: My Research Project
```

## Create a Test Applet

The [`dx-app-wizard`](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-app-wizard) generates the directory structure required to build an app. For now, there will be no inputs or outputs; we'll just use the default options when prompted.

```
$ dx-app-wizard
DNAnexus App Wizard, API v1.0.0

Basic Metadata

Please enter basic metadata fields that will be used to describe your app.  Optional fields are denoted by options with square brackets.  At the end of this wizard, the files necessary for building your app will be generated from the answers you provide.

  # We'll name our applet "My First Applet"

  # When the dx-app-wizard is done it'll construct a full, valid applet directory.
  # In this case, we'll have the following directory structure:
  my_first_applet
  ├── src
  │   └── my_first_applet.sh
  │
  ├── resources/
  ├── test/
  ├── dxapp.json
  ├── Readme.md
  └── Readme.developer.md
# As you read through the wiki and go through these tutorials, you'll understand the role that each file and directory plays in a developed app(let)
```

## Build the Sample App(let)

```bash
$ dx build my_first_applet
{"id": "applet-xxxx"} # This is the entity-ID of the applet created in the project
  # The next time that you or a colleague browse your project, they'll see this applet.
```

If you are able to perform the steps above, you're ready to start creating app(let)s!