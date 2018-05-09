# DX App(let) Tutorials

Developers new to the [DNAnexus platform](https://platform.dnanexus.com/login) may find it easier to learn by doing. This repository contains a collection of simple tutorials intended to showcase common tasks and methodologies used to create an app(let) on the DNAnexus platform. The tutorials are not meant to show realistic everyday examples; rather, they will provide a strong starting point for app(let) developers.

## Contributing to Tutorials

Follow the standard [app(let) building process](https://wiki.dnanexus.com/Developer-Tutorials/Intro-to-Building-Apps). Make sure to verify your app(let) compiles and runs. Be sure to provide helpful contextual comments in the app(let) source code.

The Readme.md for each app(let) is parsed and added to the statically generated site. When writing a readme you can reference sections or functions in the app(let) source code and include additional comments/HTML you want to only appear on the generated site using markdown comments in the following style:

```
<!-- SECTION: Section name to insert -->
<!-- INCLUDE: Text to insert -->
<!-- FUNCTION: func_name to insert -->
```
Note the file **MUST** be named Readme.md and be in the root of the app directory.

### Adding Sections

This applies to both Python and bash. In the app(let) source code you can specify the start of a section by making a single line comment
```
    # SECTION: Name Of Sections
```

You can end a section by either starting a new section or by explicitly making a single line comment 
```
    # SECTION-END
```

In the Readme.md for the app(let) you can specify a comment with the section name.
```
<!-- SECTION: Name Of Sections -->
```

### Inserting text

To insert test in the Readme.md include markdown comments with the following comment:
```
<!-- INCLUDE: Text to insert -->
```

This inserts the text directly into the Markdown in the `${repo_root}/docs/_posts/` directory right before rendering the static site. This means you can insert any special formatting, HTML, or Jekyll [includes](https://jekyllrb.com/docs/includes/) you want. Look to the [mkfifo tutorial Readme.md](bash/samtools_count_catmkfifo_sh/Readme.md) for examples.

### Inserting functions

Functions declared in python, with `def func_name:`, or in bash, with `func_name() {...}` or `function func_name {...}` syntax, you can reference the function in the Readme.md.

To insert a function in the Readme.md include a markdown comment with the following style referencing the function name:
```
<!-- FUNCTION: func_name -->
```
Look to the [SAMtools distribution by region tutorial Readme.md](python/samtools_count_distr_region_py/Readme.md) for examples


## Resources

* [DNAnexus Wiki](https://wiki.dnanexus.com/Home)
* [Developer app(let)s resource project](https://platform.dnanexus.com/projects/B406G0x2fz2B3GVk65200003/data/) (platform login required)