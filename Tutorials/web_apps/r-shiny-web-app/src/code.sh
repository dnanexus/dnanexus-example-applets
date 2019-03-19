#!/bin/bash
set -e -x -o pipefail

main() {
    echo "----Before starting R Shiny server----"

    R -e "shiny::runApp('~/my_app', host='0.0.0.0', port=443)"

    echo "----After running R Shiny server (this should probably not appear)----"
}
