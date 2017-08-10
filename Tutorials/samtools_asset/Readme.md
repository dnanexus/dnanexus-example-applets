This tutorial details how to build a SAMtools asset. For detailed information on building assets, refer to the [Asset Build Process wiki page](https://wiki.dnanexus.com/Developer-Tutorials/Asset-Build-Process).

## Building an Asset

The [`dx build_asset`](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-build-asset) command takes a source directory and builds an asset from the `dxasset.json` file and the optional Makefile and resources folder. In this example, we build an asset based on SAMtools version 0.1.19-1 available through the [Apt-Get](https://help.ubuntu.com/community/AptGet/Howto) package manager. An overview of the properties present in the `dxasset.json` file can be found on the [asset metadata wiki page](https://wiki.dnanexus.com/Developer-Tutorials/Asset-Build-Process#dxasset.json-(Asset-metadata)).

Here, once we've created the directory `samtools_asset` with the `dxasset.json` file, we run the following command from the terminal:
```
dx build_asset samtools_asset
```
This will trigger the app executable [Create asset bundle for Ubuntu 14.04](https://platform.dnanexus.com/app/app-ByQGZQ007G4yg1qP2kvGzBY3/info) to generate and upload an asset named *samtools_aptget_asset* to the platform.

## Using the Asset in an App(let)

The generated asset will exist as an [asset bundle](https://wiki.dnanexus.com/Asset-Bundle) on the platform. In the `runSpec.assetDepends` portion of the `dxapp.json` file, we can reference the created asset bundle by its [record ID](https://wiki.dnanexus.com/API-Specification-v1.0.0/Records).
```
"runSpec": {
  ...
  "assetDepends": [{"id": "record-xxxx"}],
  ...
 }
```
Now, when a worker runs an app(let) referencing this record ID, the asset bundle *samtools_aptget_asset* will be unarchived before the app(let) script is excecuted, and the command `samtools` will be usable within the app(let) script.

[DX BUILD ASSET SOURCE CODE LINK?](https://github.com/dnanexus/dx-toolkit/blob/74542f596747ee9e3c054825138f7b8082e132b0/src/python/dxpy/scripts/dx_build_asset.py)