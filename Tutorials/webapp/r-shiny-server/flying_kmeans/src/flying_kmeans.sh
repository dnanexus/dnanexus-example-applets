#!/bin/bash
# flying_kmeans 0.0.1
set -eux
main() {
    # get K-means app code
    mkdir kmeans_app
    url=https://raw.githubusercontent.com/rstudio/shiny-examples/master/050-kmeans-example
    wget -P kmeans_app/ $url/DESCRIPTION $url/server.R $url/ui.R
    # pull and run Shiny Server docker image
    # attach our K-means app's folder as a volume
    docker run --rm -p 443:3838 -v $PWD/kmeans_app:/srv/shiny-server/ rocker/shiny
}
