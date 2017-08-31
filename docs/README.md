# DNAnexus Example Applets GitHub page

This is the auto-generated static site for [DNAnexus Example Applets](https://github.com/dnanexus/dnanexus-example-applets) GitHub repository.

## Background

The inspiration for this site comes from the Cloud Cannon Knowledge *base* template for Jekyll. Cloud Cannon [Source](https://github.com/CloudCannon/base-jekyll-template) and [live demo](https://orange-ape.cloudvent.net/).

In order to substantially alter this static site, you'll need to learn Jekyll and Liquid. Cloud Cannon has a pretty nice step-by-step tutorial with videos at [CloudCannon Academy](https://learn.cloudcannon.com/).

## Getting Set Up for Local Dev

### Ruby

Jekyll and this generated website work primarily with [Ruby Gems](http://guides.rubygems.org/what-is-a-gem/). Before you can use gems, you need to have Ruby. MacOSX does come with Ruby pre-installed but you'll quickly realize MacOSX has Ruby 2.0.0 installed while Jekyll and some gems need newer versions of Ruby. I resolved this dependency using [Ruby Version Manager](https://rvm.io/) and I'll show that here, along with a caveat:
```bash
$ \curl -sSL https://get.rvm.io | bash -s stable
$ rvm install ruby-2.x.x  # Make x whatever version you need, I used 2.4.0 when creating this site
$ rvm use ruby-2.x.x --default  # This sets 2.x.x as your default Ruby version 
```

If you ran those commands you have rvm setup and working, *BUT* dx-toolkit will no longer be able to source correctly (last checked Sept 30, 2017). The nature of the conflict with RVM is unknown to me but I do know a workaround, disable rvm on terminal startup:

```bash
$ vi ~/.bash_profile  # or whatever script is sourced on terminal startup. On my OSX El Capitan it was this file

# While in vim search for this line:
# [[ -s "$HOME/.rvm/scripts/rvm" ]] && source "$HOME/.rvm/scripts/rvm" # Load RVM into a shell session *as a function*
# And comment it out. Now RVM will no longer be loaded into shells on startup
```

Feel free to install Ruby anyway you see fit, you don't have to use RVM.

### Bundler

I *HIGHLY* recommend using [Bundler](http://bundler.io/) to manage Ruby gems. It effectively a package manager (i.e. pip) for Ruby gems. Whenever installing or executing a gem use the `bundle` command. Here's how you install:
```bash
$ gem install bundler
```

That's it. You can read their [GemFiles page](http://bundler.io/gemfile.html) to get an understanding of how this site uses bundler, and why it's *SO* helpful.

### Jekyll

Basic installation instructions for Jekyll can be found on their [main site](https://jekyllrb.com/docs/installation/). You will need to install XCode, which *should* be already installed as it is a prerequisite for dx-toolkit. Again I highly recommend using bundler to install.
```bash
$ cd "$repository_root"/docs
$ bundle install
# Bundler will read the Gemfile located here and install all the dependencies to get you up and running.
# This installs a bit more than Jekyll but it gets your environment in the correct state.
```

## Basic Build Process

### Overview

This site is built statically generated from the [DNAnexus Example Applets](https://github.com/dnanexus/dnanexus-example-applets) GitHub repository. The two branches you'll have to be aware of:
* `master` - Where the documentation and RakeFile lives. The `docs/` directory contains the RakeFile the performs a majority of the build task.
* `gh-pages` - The generated static site. You will rarely (hopefully never) have to commit directly to this branch.

### Update Content

You'll rarely (read never) have to touch the `{Repository root}/docs` folder or any of its contents. You should add content normally to `{Repository root}/Example` and `{Repository root}/Tutorials` folder. The rest will be handled by Scripts and Jenkins.

### Generate Pages locally

On `master` once documentation and other site changes are made in `{Repository root}/Example` and `{Repository root}/Tutorials`, you must update the `{Repository root}/docs/_post` with changed content.

The `site_rehydrate.py` script found in `{Repository root}/docs/scripts` takes care of parsing through all the applet `dxapp.json`, `Readme.md`, `src/code.*` files in `Tutorials` and `Examples`. The simplest command to auto-generate and overwrite old pages is:
```bash
python scripts/site_rehydrate.py --overwrite
# For help use the -h argument
```

Once this is done you'll have to **manually** commit the changes to Master branch. Failure to commit the changes and moving on to the next gh-page push will correctly update the static site with the new documentation, but not update master branch. On the website, everyone will see new content but on GitHub no one will see it.
More of an annoyance, please commit changes to master.

### Rakefile

Next the Rakefile, ruby version of a Makefile, found in `docs/RakeFile`:
```bash
cd {Repository root}/docs/
rake site:publish
``` 
This Rake file will build `docs/` directory using Jekyll then  *--force* commit to the `gh-pages`. While this process is in development the commits will be forced. Once we lock down the `gh-pages` branch to only RakeFile commits we can remove *--force*

There you have it! Once the RakeFile is done you should see an updated commit on the `gh-pages` branch and the static site should be updated.

## Caveats

### POST directory

Not all files in `{Repository root}/docs/_posts` are autogenerated! Some site only documents are kept there. So be careful when deleting files there. In the future, site specific files will be moved out and copied into `_posts`... in the future.
TODO: About CSS nuance (rouge) and any other caveats.

## TODO

- Make sure links are for DNAnexus repo. Add links to direct files when possible, once merged to our Git-Org
- In Rakefile, before site_rehydrate script, check for working directory changes and exit if present.

