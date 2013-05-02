# circos Developer Readme

This app is primarily a wrapper for the Circos Perl application available 
from the [Circos website](http://circos.ca/). Please refer to that web 
site for information on how to write configuration files to generate
Circos plots, as well as an overview of Circos capabilities.

## Running this app with additional computational resources

This app has the following entry points:

* main

When running this app, you can override the instance type to be used by
providing the ``systemRequirements`` field to ```/applet-XXXX/run``` or
```/app-XXXX/run```, as follows:

    {
      systemRequirements: {
        "main": {"instanceType": "dx_m1.large"}
      },
      [...]
    }

See <a
href="http://wiki.dnanexus.com/API-Specification-v1.0.0/IO-and-Run-Specifications#Run-Specification">Run
Specification</a> in the API documentation for more information about the
available instance types.
