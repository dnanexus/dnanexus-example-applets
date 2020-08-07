<!-- dx-header -->
# RStudio Server (DNAnexus Platform App)

RStudio Server on DNAnexus

This is an example code of wrapping the RStudio Server software into a DNAnexus platform applet
 that can then be run as a web app.

The RStudio image and setup instructions: https://hub.docker.com/r/rocker/rstudio/

<!-- /dx-header -->

<!-- Insert a description of your app here -->
## Procedure of creating the applet


##### 1. Pull ready-to-use RStudio docker image from Docker Hub to your local computer by running in terminal:
```bash
docker pull rocker/rstudio
```


##### 2. Save it to a file `rstudio.docker.gz`:
```bash
docker save -o rstudio.docker.gz rocker/rstudio
```


##### 3. Log in
Log into the DNAnexus platform and navigate to the project where you'd like to upload the RStudio docker image and build the applet.


##### 4. Upload the file to your DNAnexus project e.g. to `/rstudio/rstudio.docker.gz`
```bash
mkdir /rstudio
dx cd /rstudio
dx upload rstudio.docker.gz
```


##### 5. Build and deploy RStudio web applet
```bash
# Navigate to the directory with the source code of this applet
# If you cloned the whole repository then it will be
# /dnanexus-example-applets/Tutorials/webapp/rstudio_server
cd .../rstudio_server

# Build the applet
# Please note that your current directory must be the directory with this applet's code.
dx build .
```

If this is the first time of running this command, the applet will be built and placed in the current directory of your DNAnexus project.

If there is a previously built applet with the same name in the current directory, the build command will throw error that the applet already exists at this path. In this case, run `dx build -f .` instead to overwrite the old version of the applet. If you prefer to archive old versions of built applets, add `-a` argument instead of `-f`. This will move the existing applet to `/.Applet_archive` directory of your DNAnexus project, add a timestamp to it, and replace it with the new applet in the current directory. For help about the build process, `run dx build -h` or go to our Documentation page.


##### 6. Start the RStudio DNAnexus web applet
To start the applet, click on it. In the app's window, click on "Run as Analysis…" button. The app will start launching that can be seen in the platform's Monitor tab.

In a few minutes, the app's status will change to Running. After a moment, a web URL will appear in the "Worker URL" column of the table. This means that the app is ready to accept connections.

NOTE 1: If you don't see the URL for a while, even though the job is running, just reload the page.


##### 7. Use the app
Click on the app's URL that appears in the "Worker URL" column of the Monitor tab. A new  web browser tab will open and display the login page of the RStudio Server. Enter `rstudio` as username and `yourpasswordhere` as password, log into the server, and enjoy! :-)

NOTE 1: If you prefer to set a different password for your RStudio Server, change it in the file `rstudio_server.sh` of this source code.

NOTE: If you got "502 Bad Gateway" browser error when you loaded the page, please wait about a minute or so and reload the page. This happens when the web server part is ready and the URL is available but internally the RStudio Server's docker container has not finished loading and is not ready to listen to web connections yet.
