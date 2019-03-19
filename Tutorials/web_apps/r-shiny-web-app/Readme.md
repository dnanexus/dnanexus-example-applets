<!-- dx-header -->
# R Shiny Example Web App (DNAnexus Platform App)

Demonstrates how to use R Shiny on DNAnexus

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

## Creating an applet on DNAnexus
Build the asset with the libraries first:
```
dx build_asset shiny-asset
```

Take the record ID it outputs and add it to the dxapp.json for the applet.
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

Then build the applet
```
dx build -f r-shiny-web-app
dx run --debug-on All r-shiny-web-app
```
Then go to that job's special URL https://job-xxxx.dnanexus.cloud/ and see the result!
