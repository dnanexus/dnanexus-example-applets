This is an example web app made with Dash, which in turn uses Flask underneath.

## Creating the web application
After configuring an `app` with Dash, we start the server on port 443.
```
app.run_server(host='0.0.0.0', port=443)
```
Inside the `dxapp.json`, you would add `"httpsApp": {"ports":[443], "shared_access": "VIEW"}` to tell the worker to expose this port.

Note that for all web apps, if everything is running smoothly and no errors are encountered (the ideal case), the line of code that starts the server will keep it running forever. The applet stops only when it is terminated. This also means that any lines of code after the server starts will not be executed.

The rest of these instructions apply to building any applet with dependencies stored in an asset.

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
```

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
```
You can always use `dx ssh job-xxxx` to ssh into the worker and inspect what's going on or experiment with quick changes 
Then go to that job's special URL https://job-xxxx.dnanexus.cloud/ and see the result!


## Optional local testing
The main code is in `dash-web-app/resources/home/dnanexus/my_app.py` with a local launcher script called `local_test.py` in the same folder. This allows us to launch the same core code in the applet locally to quickly iterate. This is optional because you can also do all testing on the platform itself.

Install locally the same libraries listed above.

To launch the web app locally:
```
cd dash-web-app/resources/home/dnanexus/
python local_test.py
```
Once it spins up, you can go to that job's designated URL based on its job ID, https://job-xxxx.dnanexus.cloud/, to see the result.
