# DNAnexus Example Applets GitHub page

This is the auto-generated static site for [DNAnexus Example Applets](https://github.com/dnanexus/dnanexus-example-applets) GitHub repository.

## Background

The inspiration for this site comes from the Cloud Cannon Knowledge *base* template for Jekyll. Cloud Cannon [Source](https://github.com/CloudCannon/base-jekyll-template) and [live demo](https://orange-ape.cloudvent.net/).

In order to substantially alter this static site, you'll need to learn some Jekyll and Liquid. Cloud Cannon has a pretty nice step-by-step tutorial with videos at [CloudCannon Academy](https://learn.cloudcannon.com/).

## Getting Set Up for Local Dev

### Ruby

Jekyll and this generated website work primarily with [Ruby Gems](http://guides.rubygems.org/what-is-a-gem/). Before you can use gems, you need to have Ruby. MacOSX does come with Ruby pre-installed but you'll quickly realize MacOSX has Ruby 2.0.0 installed while Jekyll and some gems need newer versions of Ruby. I resolved this dependency using [Ruby Version Manager](https://rvm.io/) and I'll show that here, along with a caveat:
```bash
$ \curl -sSL https://get.rvm.io | bash -s stable
$ rvm install ruby-2.x.x  # Make x whatever version you need, I used 2.4.0 when creating this site
$ rvm use ruby-2.x.x --default  # This sets 2.x.x as your default Ruby version 
```

If you ran those commands you have rvm setup and working, **BUT** dx-toolkit will no longer be able to source correctly (last checked Nov 10, 2017). The nature of the conflict with RVM is unknown to me but I do know a workaround, prevent rvm sourcing on shell startup:

```bash
$ vi ~/.bash_profile  # or whatever script is sourced on terminal startup. On my OSX El Capitan it was this file

# While in vim search for this line:
# [[ -s "$HOME/.rvm/scripts/rvm" ]] && source "$HOME/.rvm/scripts/rvm" # Load RVM into a shell session *as a function*
# And comment it out. Now RVM will no longer be sourced into shells on startup
```

Feel free to install Ruby anyway you see fit, you don't have to use RVM.

### Bundler

I *HIGHLY* recommend using [Bundler](http://bundler.io/) to manage Ruby gems. It effectively a package manager (i.e. pip) for Ruby gems:
```bash
$ gem install bundler
```

That's it. You can read their [GemFiles page](http://bundler.io/gemfile.html) to get an understanding of how this site uses bundler, and why it's *SO* helpful. Whenever installing or executing a gem use the `bundle` command.

### Jekyll

Basic installation instructions for Jekyll can be found on their [main site](https://jekyllrb.com/docs/installation/). You will need to install XCode, which *should* already be installed as it is a prerequisite for dx-toolkit. Again I highly recommend using bundler to install.
```bash
$ cd "$repository_root"/docs
$ bundle install
# Bundler will read the Gemfile located here and install all the dependencies to get you up and running.
# This installs a bit more than Jekyll but it gets your environment in the correct state.
```

## Basic Build Process

### Overview

This site is statically generated from the [DNAnexus Example Applets](https://github.com/dnanexus/dnanexus-example-applets) GitHub repository. The two branches you'll have to be aware of:
* `master` - Where the documentation and RakeFile lives. The `docs/` directory contains the RakeFile that performs a majority of the build task.
* `gh-pages` - The generated static site. You will rarely (hopefully never) have to commit directly to this branch. The RakeFile will do this work.

### Update Content

You'll rarely (read never) have to touch the `{Repository root}/docs` folder or any of its contents. You should add content normally to `{Repository root}/Example` and `{Repository root}/Tutorials` folder. The rest will be handled by Scripts and Jenkins.

### Generate Pages locally

On `master` once documentation and other site changes are made in `{Repository root}/Example` and `{Repository root}/Tutorials`, you must execute the `site_rehydrate.py` script which updates the `{Repository root}/docs/_post` with changed content.

The `site_rehydrate.py` script found in `{Repository root}/docs/scripts` takes care of parsing through all the applet `dxapp.json`, `Readme.md`, `src/code.*` files in `Tutorials` and `Examples`. The simplest command to auto-generate and overwrite old pages is:
```bash
python scripts/site_rehydrate.py --overwrite
# For help use the -h argument
```

Once pages are automatically generated preview them using Jekyll.
```bash
$ bundle exec jekyll serve
Configuration file: {Repository root}/docs/_config.yml
Configuration file: {Repository root}/docs/_config.yml
            Source: {Repository root}/docs
       Destination: {Repository root}/docs/_site
 Incremental build: disabled. Enable with --incremental
      Generating... 
                    done in 2.6 seconds.
 Auto-regeneration: enabled for '{Repository root}/docs'
Configuration file: {Repository root}/docs/_config.yml
    Server address: http://127.0.0.1:4000/dnanexus-example-applets/
  Server running... press ctrl-c to stop.
```

Copy and paste the `Server address` into a browser. Now you see your dev version of dnanexus-example-apps.

### Rakefile

We use a Rakefile, ruby version of a Makefile, found in `docs/RakeFile` to publish to `gh-pages` branch:
```bash
cd {Repository root}/docs/
rake site:publish
``` 
This Rake file will build `docs/` directory using Jekyll then  *--force* commit to the `gh-pages`. You CANNOT publish if you have working/staged/untracked changes in your current git working tree. The goal is to keep some level of version control between the `gh-pages` branch and `master`.

While this process is in development the commits will be forced. Once we lock down the `gh-pages` branch to only RakeFile commits we can remove *--force*.

There you have it! Once the RakeFile is done you should see an updated commit on the `gh-pages` branch and the static site should be updated.

## Caveats

### \_post/ directory

Not all files in `{Repository root}/docs/_posts` are auto-generated! Some site only documents are kept there. So be careful when deleting files there. In the future, site specific files will be moved out ... in the future.

### site.baseurl

GitHub pages serves our site. It does not serve the site using Jekyll. Our build process builds our \_site directory then copies/commits it to gh-pages branch.

way to write liquid links:

for internal links: `{{ site.url }}{{ site.baseurl }}{{ post.url }}`
for css/js/img and such: `{{ site.url }}{{ site.baseurl }}/path/to/css.css`

Depending on where Jekyll serves, development (local port) versus production (on github), the variable baseurl resolves differently. As long as you use the full path to links you'll be able to follow links on your local machine and on a production webpage.

Our rakefile applies [a workaround](https://github.com/jekyll/jekyll/issues/332) to keep the production and dev process the same:
```ruby
Jekyll::Site.new(Jekyll.configuration({
  "source"      => ".",
  "destination" => "_site",
  "baseurl"     => ""  # we set baseurl when we build our _site
})).process
```

### CSS and SASS

We make use of [Sass](http://sass-lang.com/) to compile our css sheet. You can look at the main site to figure out how their compiler works, but one piece of advice right now.

The compiler tries to put comments at the top of the compiled css file (last checked Nov 10, 2017). Meaning the order of your `@import` statements can be ignored.
