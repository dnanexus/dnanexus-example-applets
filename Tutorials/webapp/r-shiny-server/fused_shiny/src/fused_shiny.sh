#!/bin/bash
# dxshiny 0.0.1

set -eux

main() {

    # Download the input: your Shiny app's source code archive from the project
    echo "Downloading Shiny app's source code: '$app_gz'"
    dx download "$app_gz" -o app.gz

    # Unpack the app code. Here, `--strip 1` allows to skip the top folder
    # of the archive and unpack the code directly underneath ./app/
    # This allows you to go to the job's URL and launch your Shiny app directly.
    # Otherwise, the code would be one folder level deeper and one would have to
    # click on the app folder's link in the browser.
    mkdir app
    tar -zxvf app.gz -C app --strip 1

    # Mount the parent project using dxFUSE
    wget https://github.com/dnanexus/dxfuse/releases/download/v0.21/dxfuse-linux
    chmod +x dxfuse-linux
    source environment >& /dev/null
    FUSE_MOUNT=$HOME/projects
    mkdir -p $FUSE_MOUNT
    sudo -E ./dxfuse-linux -uid $(id -u) -gid $(id -g) -verbose 2 $FUSE_MOUNT $DX_PROJECT_CONTEXT_ID
    PROJ_NAME=$(dx describe $DX_PROJECT_CONTEXT_ID --name)
    PROJ_PATH=${FUSE_MOUNT}/${PROJ_NAME}
    echo "Project ${PROJ_NAME} mounted as ${PROJ_PATH}"

    # Set up the Shiny Server
    # Assuming that the archive of the docker image of the Shiny Server is
    # in your project at /rshiny/rshiny.docker.gz
    dx download $DX_PROJECT_CONTEXT_ID:/rshiny/rshiny.docker.gz
    # ($DX_PROJECT_CONTEXT_ID contains your project ID)
    # Now load the image from the archive, create docker image
    docker load -i rshiny.docker.gz
    # Run the docker image. Attach your app's folder as a volume.
    #docker run --rm -p 443:3838 -v $PWD/app:/srv/shiny-server/ rocker/shiny
    docker run --rm -p 443:3838 -v $PWD/app:/srv/shiny-server/ \
      -v $PROJ_PATH:/srv/project/ \
      rocker/shiny
}
