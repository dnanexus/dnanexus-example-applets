# Example web applications

See README files inside each example directory.

In short, to create a web app, you need to:
1. Launch a local server within the app on port 443. Dash, R Shiny, Flask, and many other programs let you do this, just make sure to set the port to 443.
2. Add this to the `dxapp.json`: `"httpsApp": {"ports":[443], "shared_access": "VIEW"}` at the top level of the JSON
3. Then go to the job's URL https://job-xxxx.dnanexus.cloud to see it once it has started running.
