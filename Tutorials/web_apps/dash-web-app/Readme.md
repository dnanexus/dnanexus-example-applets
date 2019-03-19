<!-- dx-header -->
# Dash Example Web App (DNAnexus Platform App)

An example of building a web applet with Dash

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://wiki.dnanexus.com/.
<!-- /dx-header -->

## Building your own web app
To create a web app, you need to:
1. Make sure you have Titan or Apollo access on DNAnexus, otherwise the ability to make web apps is not enabled. 
2. Launch a local server within the app on port 443. Dash, R Shiny, Flask, and many other programs let you do this, just make sure to set the port to 443.
3. Add this to the `dxapp.json`: `"httpsApp": {"ports":[443], "shared_access": "VIEW"}` at the top level of the JSON.
4. Then go to the job's URL https://job-xxxx.dnanexus.cloud to see it once it has started running.


The main code is in `dash-web-app/resources/home/dnanexus/my_app.py` with a local launcher script called `local_test.py` in the same folder. This is a web app made with Dash, which in turn uses Flask underneath, and it can be launched locally and as an applet on DNAnexus.

## Local use
Install libraries:
```
pip install dash==0.39.0  # The core dash backend
pip install dash-html-components==0.14.0  # HTML components
pip install dash-core-components==0.44.0  # Supercharged components
pip install dash-table==3.6.0  # Interactive DataTable component (new!)
pip install dash-daq==0.1.0  # DAQ components (newly open-sourced!)
pip install pandas # for leading CSV data into python
```

To launch the web app locally:
```
cd dash-web-app/resources/home/dnanexus/
python local_test.py
```

## Creating an applet on DNAnexus
Source dx-toolkit and log in, then run dx-app-wizard with default options.

### Creating the asset
dash-asset specifies all the packages and versions we need.
We take these from the Dash installation guide (https://dash.plot.ly/installation)
```
pip install dash==0.39.0  # The core dash backend
pip install dash-html-components==0.14.0  # HTML components
pip install dash-core-components==0.44.0  # Supercharged components
pip install dash-table==3.6.0  # Interactive DataTable component (new!)
pip install dash-daq==0.1.0  # DAQ components (newly open-sourced!)
```

We put these into dash-asset/dxasset.json:
```
{
  ...
  "execDepends": [
    {"name": "dash", "version":"0.39.0", "package_manager": "pip"},
		{"name": "dash-html-components", "version":"0.14.0", "package_manager": "pip"},
		{"name": "dash-core-components", "version":"0.44.0", "package_manager": "pip"},
		{"name": "dash-table", "version":"3.6.0", "package_manager": "pip"},
		{"name": "dash-daq", "version":"0.1.0", "package_manager": "pip"}
  ],
	...
}

Build the asset:
```
dx build_asset dash-asset
```

## Use the asset from the applet
Add this asset to the applet's dxapp.json:
```
"runSpec": {
	...
	"assetDepends": [
    {
      "id": "record-xxxx
    }
  ]
	...
}
```

## Build the applet
Now build and run the applet itself:
```
dx build -f dash-web-app
dx run dash-web-app

# If you need to ssh into the job to check on things when/if it fails:
# dx ssh job-xxxx
```

Then go to that job's special URL https://job-xxxx.dnanexus.cloud/ and see the result!
