# samtools_count_bashhelper_sh (DNAnexus Platform App)

This applet performs a basic SAMtools count with the aid of bash [helper variables](https://wiki.dnanexus.com/Developer-Tutorials/Advanced-App-Tutorial#Writing-the-applet-script).

## Download bam files
This script uses the DNAnexus [`dx-download-all-inputs`](https://wiki.dnanexus.com/Helpstrings-of-SDK-Command-Line-Utilities#dx-download-all-inputs) script to download all job inputs. For file type inputs, this script will go through all inputs and download into folders which have the pattern
_*/home/dnanexus/in/[INPUT]/*_.  
<!--SECTION: Download bam files -->

## Create output directory

Create an output directory in preparation for `dx-upload-all-outputs` DNAnexus script.
<!--SECTION: Create output directory -->

## Run samtools view

Here, we use the bash helper variable `mappings_bam_path`. This bash variable
will point to the location of a file after it has been downloaded using
`dx-download-all-inputs`. So for *[VARIABLE]_path*, this wll point to
*/home/dnanexus/in/[VARIABLE]/[VARIABLE_FILENAME]*.
In the case that the filename of the file mappings_bam is *my_mappings.bam*,
mappings_bam_path will be */home/dnanexus/in/mappings_bam/my_mappings.bam*.  
<!--SECTION: Run samtools view -->

## Upload result

We use `dx-upload-all-outputs` to upload the data to the platform and associate
it as the output. `dx-upload-all-outputs` uses the folder pattern
*/home/dnanexus/out/[VARIABLE]*. It will upload the files in that folder
and then mark is as the output corresponding to *[VARIABLE]*. In this case,
the output is called `counts_txt`. [Above](#Create-output-directory), we have made the folder which
corresponds to this, and placed the output there.  
<!--SECTION: Upload result -->
