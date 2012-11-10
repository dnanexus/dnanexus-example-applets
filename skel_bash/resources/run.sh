#!/bin/bash -e

PATH="/:${PATH}"
HOME="`pwd`"

# Example for retrieving file-123 from input: {"alignment": {"$dnanexus_link": "file-123"}}
#jshon -e alignment -e '$dnanexus_link' -u < job_input.json

echo '{}' > "$HOME/job_output.json"
