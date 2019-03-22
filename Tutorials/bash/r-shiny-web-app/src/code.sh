#!/bin/bash
set -e -x -o pipefail

main() {
    R -e "shiny::runApp('~/my_app', host='0.0.0.0', port=443)"
}
