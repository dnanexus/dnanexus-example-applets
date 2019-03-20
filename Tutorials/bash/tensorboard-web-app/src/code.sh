#!/bin/bash
set -e -x -o pipefail

main() {
    mkdir -p LOGS_FOR_TENSORBOARD

    echo "----Starting training script----"
    # Start the training script and put it into the background,
    # so the next lines of code will run immediately
    python mnist_tensorboard_example.py --log_dir LOGS_FOR_TENSORBOARD &

    echo "----Set httpsAppState=running----"
    dx set_properties $DX_JOB_ID httpsAppState=running

    echo "----Starting Tensorboard server----"
    # Run TensorBoard
    tensorboard  --logdir LOGS_FOR_TENSORBOARD --host 0.0.0.0 --port 443

}
