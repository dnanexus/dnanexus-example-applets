#!/bin/bash
# rstudio_server 0.0.1

set -eux

main() {
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

    # Set up the RStudio Server
    # Assuming that the archive of the docker image of the RStudio Server is
    # in your project at /rstudio/rstudio.docker.gz
    dx download $DX_PROJECT_CONTEXT_ID:/rstudio/rstudio.docker.gz
    # ($DX_PROJECT_CONTEXT_ID contains your project ID)
    # Alternatively, the file can be downloaded using its file ID.
    # This is more flexible because the file can be in any folder path.

    # Now load the image from the archive, create docker image
    docker load -i rstudio.docker.gz

    # Run the docker image. Attach your app's folder as a volume.
    docker run --rm \
      -p 443:8787 \
      -e ROOT=TRUE \
      -e PASSWORD=pass \
      -e DISABLE_AUTH=true \
      -v ${PROJ_PATH}:/home/rstudio/project \
      rocker/rstudio
    # Here:
    # `-p 443:8787` maps RStudio's internal port 8787 to the external HTTPS
    # `-e ROOT=TRUE` allows you to run commands with `sudo`
    # `-e PASSWORD=pass` defines password for the default user `rstudio`
    # `-e DISABLE_AUTH=true` turns off the Log In window and logs into server directly
    # `-v ${PROJ_PATH}:/home/rstudio/project` mounts the parent project to the `project` under your home directory
    # `rocker/rstudio` is loaded docker image name; the Rocker project has derivative images too e.g. rocker/tidyverse
}
