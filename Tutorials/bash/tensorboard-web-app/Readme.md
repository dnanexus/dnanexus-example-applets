This example demonstrates how to run TensorBoard inside a DNAnexus applet.

TensorBoard is a web application used to visualize and inspect what is going on inside TensorFlow training. To use TensorBoard, our training script in TensorFlow needs to include code that saves various data to a log directory where TensorBoard can then find the data to display it.

This example uses an example script from the TensorBoard authors. For more guidance on how to use TensorBoard, check out the tensorflow website ([external link](https://www.tensorflow.org/guide/summaries_and_tensorboard)).

## Creating the web application
The applet code runs a training script, which is placed in `resources/home/dnanexus/` to make it available in the current working directory of the worker, and then it starts tensorboard on port 443 (HTTPS).
```
# Start the training script and put it into the background,
# so the next line of code will run immediately
python mnist_tensorboard_example.py --log_dir LOGS_FOR_TENSORBOARD &

# Run TensorBoard
tensorboard  --logdir LOGS_FOR_TENSORBOARD --host 0.0.0.0 --port 443
```
We run the training script in the background to start TensorBoard immediately, which will let us see the results while training is still running. This is particularly important for long-running training scripts.

Note that for all web apps, if everything is running smoothly and no errors are encountered (the ideal case), the line of code that starts the server will keep it running forever. The applet stops only when it is terminated. This also means that any lines of code after the server starts will not be executed.

As with all web apps, the `dxapp.json` must include `"httpsApp": {"ports":[443], "shared_access": "VIEW"}` to tell the worker to expose port 443.

## Creating an applet on DNAnexus
Build the asset with the libraries first:
```
dx build_asset tensorflow_asset
```

Take the record ID it outputs and add it to the dxapp.json for the applet.
```
"runSpec": {
	...
	"assetDepends": [
    {
      "id": "record-xxxx
    }
  ]
	...
}
```

Then build the applet
```
dx build -f tensorboard-web-app
dx run tensorboard-web-app
```
Once it spins up, you can go to that job's designated URL based on its job ID, https://job-xxxx.dnanexus.cloud/, to see the result.
