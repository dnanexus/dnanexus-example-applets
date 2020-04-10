#!/bin/bash
# dxshiny 0.0.1

set -eux

main() {

    # Download the input: your Shiny app's source code archive from the project
    echo "Value of app_gz: '$app_gz'"
    dx download "$app_gz" -o app.gz

    # Unpack the app code. Here, `--strip 1` allows to skip the top folder
    # of the archive and unpack the code directly underneath ./app/
    # This allows you to go to the job's URL and launch your Shiny app directly.
    # Otherwise, the code would be one folder level deeper and one would have to
    # click on the app folder's link in the browser.
    mkdir app
    tar -zxvf app.gz -C app --strip 1

    # Set up the Shiny Server
    # Assuming that the archive of the docker image of the Shiny Server is
    # in your project at /rshiny/rshiny.docker.gz
    dx download $DX_PROJECT_CONTEXT_ID:/rshiny/rshiny.docker.gz
    # ($DX_PROJECT_CONTEXT_ID contains your project ID)
    # Now load the image from the archive, create docker image
    docker load -i rshiny.docker.gz
    # Run the docker image. Attach your app's folder as a volume.
    docker run --rm -p 443:3838 -v $PWD/app:/srv/shiny-server/ rocker/shiny
}
