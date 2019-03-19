#!/bin/bash
set -e -x -o pipefail

main() {
    mkdir -p LOGS_FOR_TENSORBOARD

    echo "----Starting training script----"
    python mnist_tensorboard_example.py --log_dir LOGS_FOR_TENSORBOARD &

    echo "----Set httpsAppState=running----"
    dx set_properties $DX_JOB_ID httpsAppState=running

    echo "----Starting Tensorboard server----"
    tensorboard  --logdir LOGS_FOR_TENSORBOARD --host 0.0.0.0 --port 443

}
