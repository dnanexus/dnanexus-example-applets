This is an example web applet that demonstrates how to build and run an R Shiny application on DNAnexus.

## Creating the web application
Inside the `dxapp.json`, you would add `"httpsApp": {"ports":[443], "shared_access": "VIEW"}` to tell the worker to expose this port.

R Shiny needs two scripts, `server.R` and `ui.R`, which should be under `resources/home/dnanexus/my_app/`. When a job starts based on this applet, the `resources` directory is copied onto the worker, and since the `~/` path on the worker is `/home/dnanexus`, that means you now have `~/my_app` with those two scripts inside.

From the main applet script `code.sh`, we simply start shiny pointing to `~/my_app` repo, serving its mini-application on port 443.

```
main() {
  R -e "shiny::runApp('~/my_app', host='0.0.0.0', port=443)"
}
```

Note that for all web apps, if everything is running smoothly and no errors are encountered (the ideal case), the line of code that starts the server will keep it running forever. The applet stops only when it is terminated. This also means that any lines of code after the server starts will not be executed.

## Modifying this example for your own applet
To make your own applet with R Shiny, simply copy the source code from this example and modify `server.R` and `ui.R` inside `resources/home/dnanexus/my_app`.

### How to rebuild the shiny asset
```
dx build_asset shiny-asset
```
This will output a record ID `record-xxxx` that you can then put into the applet's `dxapp.json` in place of the existing one:
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

### Build the applet
Now build and run the applet itself:
```
dx build -f r-shiny-web-app
dx run r-shiny-web-app
```
Once it spins up, you can go to that job's designated URL based on its job ID, https://job-xxxx.dnanexus.cloud/, to see the result.
