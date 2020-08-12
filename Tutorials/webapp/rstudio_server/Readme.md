<!-- dx-header -->
# RStudio Server (DNAnexus Platform App)

This is an example code of wrapping the RStudio Server software into a DNAnexus platform's web applet.

The main source file is `src/rstudio_server.sh`. So, modify this file to change the behavior of the applet.

## Used technologies
The app relies on the following technologies:
##### 1. RStudio Server's docker image 
We use the image from the Rocker project and its default setup instructions. If you would like to tweak RStudio settings , please refer to https://hub.docker.com/r/rocker/rstudio/ .

##### 2. dxFUSE
To mount the parent DNAnexus project as a local directory, we used dxFUSE ( https://github.com/dnanexus/dxfuse/ ). The project files will appear in the home directory. 

<!-- /dx-header -->

<!-- Insert a description of your app here -->
## Procedure of creating the applet


##### 1. Pull ready-to-use RStudio docker image from Docker Hub to your local computer by running in terminal:
We are going to pull the image locally and upload it to the DNAnexus project instead of pulling the image directly from the Docker Hub. Caching the image in the project causes applet to load faster, as well as eliminates possibility that the applet will fail to load if the internet connection between DNAnexus and Docker Hub has problems.
 
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
dx mkdir /rstudio
dx cd /rstudio
dx upload rstudio.docker.gz
```


##### 5. Build and deploy RStudio web applet

Navigate to the directory with the source code of this applet. If you cloned the whole repository then the path will be `/dnanexus-example-applets/Tutorials/webapp/rstudio_server`.

In the following command, replace `.../rstudio_server` with the correct path to the applet's directory:
```bash
cd .../rstudio_server
```

The contents of the underlying directory must look like this:
```
└── rstudio_server
    ├── Readme.developer.md
    ├── Readme.md
    ├── dxapp.json
    ├── resources
    ├── src
    │   └── rstudio_server.sh
    └── test
```

Now build and deploy the applet:
```bash
dx build .
```

If this is the **first time** of running this command, the applet will be built and placed in the current directory of your DNAnexus project.

**If there is a previously built applet with the same name in the current directory, the build command will throw error** that the applet already exists at this path. In this case, run `dx build -f .` instead to overwrite the old version of the applet. If you prefer to archive old versions of built applets, add `-a` argument instead of `-f`. This will move the existing applet to `/.Applet_archive` directory of your DNAnexus project, add a timestamp to it, and replace it with the new applet in the current directory. For help about the build process, `run dx build -h` or go to our Documentation page at https://documentation.dnanexus.com/developer/apps/app-build-process .


##### 6. Start the RStudio DNAnexus web applet
To start the applet, click on it. In the app's window, click on "Run as Analysis…" button. The app will start launching that can be seen in the platform's Monitor tab.

In a few minutes, the app's status will change to Running. After a moment, a web URL will appear in the "Worker URL" column of the table. This means that the app is ready to accept connections.

**_NOTE:_**: If you don't see the URL for a while, even though the job is running, just reload the page.


##### 7. Use the app
Click on the app's URL that appears in the "Worker URL" column of the Monitor tab. A new  web browser tab will open and display the login page of the RStudio Server. Enter `rstudio` as username and `pass` as password, and log into the server. You will see familiar RStudio interface. Your parent DNAnexus project will be mounted to `projects` under your home directory. You can read files from it, work on them, and write back.

Enjoy! :-)

**_NOTE 1:_** If you prefer to set a different password for your RStudio Server, change it in the file `rstudio_server.sh` of this source code.

**_NOTE 2:_** If you got "502 Bad Gateway" browser error when you loaded the page, please wait about a minute or so and reload the page. This happens when the web server part is ready and the URL is available but internally the RStudio Server's docker container has not finished loading and is not ready to listen to web connections yet.




# Working with files
This RStudio server setup allows you to read your files quickly from your parent DNAnexus. Your project's files are available in the  `project` folder under your home directory. Just navigate to that folder and open your R files. If you would like to work on these scripts for a while, it is better to temporarily save them locally in your RStudio server. To do so, click File > Save As... and then save the file in your home directory anywhere outside the `~/project` folder. Once you are done working with the file, click Save As again and save back to the `~/project` folder. This way, you will avoid "bothering" the parent DNAnexus project for every file save. 
**_NOTE:_** The dxFUSE does not immediately upload files to the project. Therefore, give it a couple of minutes and then see that your file will appear in the project.


# FAQ
##### How do I install Tidyverse?
If you try to install Tidyverse library and load it, you may get an error about a missing libxml2 library, which apparently is not part of the Rocker image. To solve this, run in Terminal:
```bash
sudo apt install libxml2
```
This will install the essential library.

Now, run in the R Console:
```R
install.packages("tidyverse")
library(tidyverse)
```

There is also another Docker image that you can use. In the Step 1 above, run
```bash
docker pull rocker/tidyverse
```
instead of
```bash
docker pull rocker/rstudio
```
This image contains pre-installed Tidyverse . Read more at:
https://github.com/rocker-org/rocker/wiki and https://hub.docker.com/r/rocker/rstudio/

##### How do I turn on log in window?
If you want to see the Log In window every time you start a new RStudio server (which I definitely don't), modify the docker run command in the `rstudio_server.sh` file by removing `-e DISABLE_AUTH=true`.

Now it would look like this:
```bash
docker run --rm \
  -p 443:8787 \
  -e ROOT=TRUE \
  -e PASSWORD=pass \
  -v ${PROJ_PATH}:/home/rstudio/project \
  rocker/rstudio
```
Please note that the username will be `rsutio` and `pass` will be your password, as set in the command line. 

