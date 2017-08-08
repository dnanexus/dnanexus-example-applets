# DNAnexus Example Applets GitHub page

This is the auto-generated static site for [DNAnexus Example Applets](https://github.com/dnanexus/dnanexus-example-applets) GitHub page.

## Overview

Inspiration for this site comes from the Cloud Cannon Knowledge base template for Jekyll. Cloud Cannon [Source](https://github.com/CloudCannon/base-jekyll-template) and [live demo](https://orange-ape.cloudvent.net/).

In order to substantially alter this static site, you'll need to learn Jekyll and Liquid. Cloud Cannon has a pretty nice step-by-step tutorial with videos at [CloudCannon Academy](https://learn.cloudcannon.com/).

## Basic Build Process

This site is built statically from the Example applets repo. The two branches you'll have to be aware off:
`master` - Where the documentation and RakeFile lives. The `docs/` directory contains the RakeFile the performs a majority of the build task.
`gh-pages` - The generated static site. You will rarely have to commit directly to this branch.

On `master` once documentation and other site changes are made, you must first update the `docs/_post` with changed content. The `site_rehydrate.py` script found in `docs/scripts` takes care of parsing through all the applet `dxapp.json`, `Readme.md`, `src/code.*` files in `Tutorials` and `Examples`. A letter section will cover its usage. For now, a simple command:
```python
python scripts/site_rehydrate.py --overwrite
```

Once this is done you'll have to **manually** commit the changes to Master branch. Failure to commit the changes and moving on to the next will update the static site with the new documentation without updating master. More of an annoyance, please commit changes to master.

Next the Rakefile, ruby version of a Makefile, found in `docs/RakeFile` must be executed:
```bash
cd docs/
rake site:publish
``` 
This Rake file will build `docs/` directory using Jekyll then  *--force* commit to the `gh-pages`. While this process is in development the commits will be forced. Once we lock down the `gh=pages` to only RakeFile commits we can remove *--force*

There you have it! Once the RakeFile is done you should see an updated commit on the `gh-pages` branch and the static site should be updated.

## Updating

TODO: how to use `site_rehydrate.py` script.

## Details

TODO: About CSS nuance (rouge) and any other caveats