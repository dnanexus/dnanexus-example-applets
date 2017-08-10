This applet tutorial will perform a SAMtools count using parallel threads.

In order to take full advantage of the scalability that cloud computing offers, our scripts have to implement the correct methodologies. This applet tutorial will:
1. Install SAMtools
2. Download BAM file
3. Split workload
4. Count regions in parallel

This applet tutorial code is extremely similar to the [_Parallel Threads SAMtools count tutorial_](/python_parallel_tutorial.html#samtools_count_para_chr_subprocess_py), except `multiprocessing` is used instead of `multiprocessing.dummy`.

## How is the SAMtools dependency provided?
The SAMtools dependency is resolved by declaring an [Apt-Get](https://help.ubuntu.com/14.04/serverguide/apt-get.html) package in the `dxapp.json` `runSpec.execDepends` field.
```json
{
  "runSpec": {
    ...
    "execDepends": [
      {"name": "samtools"}
    ]
  }
```
For additional information, please refer to the [execDepends wiki page](https://wiki.dnanexus.com/Execution-Environment-Reference#Software-Packages).

## Download Inputs

This applet downloads all inputs at once using `dxpy.download_all_inputs`:
```python
inputs = dxpy.download_all_inputs()
# download_all_inputs returns a dictionary that contains mapping from inputs to file locations.
# Additionaly, helper keys, value pairs are added to the dicitonary, similar to bash helper functions
inputs
#     mappings_sorted_bam_path: [u'/home/dnanexus/in/mappings_sorted_bam/SRR504516.bam']
#     mappings_sorted_bam_name: u'SRR504516.bam'
#     mappings_sorted_bam_prefix: u'SRR504516'
#     mappings_sorted_bai_path: u'/home/dnanexus/in/mappings_sorted_bai/SRR504516.bam.bai'
#     mappings_sorted_bai_name: u'SRR504516.bam.bai'
#     mappings_sorted_bai_prefix: u'SRR504516'
```
## Split workload
We parallely process by using the python `multiprocessing` module. A rather simple pattern you can follow to compute in a parallel manner is shown below:
```python
print "Number of cpus: {0}".format(cpu_count())  # Get cpu count from multiprocessing
worker_pool = Pool(processes=cpu_count())  # Create a pool of workers, 1 for each core
results = worker_pool.map(run_cmd, view_cmds)  # map run_cmds to a collection
                                               # Pool.map will handle orchestrating the job
worker_pool.close()
worker_pool.join()  # Make sure to close and join workers when done
```
This convenient pattern allows you to quickly orchestrate jobs on a worker. For more detailed overview of the `multiprocessing` module, visit the [python2.7 docs](https://docs.python.org/2/library/multiprocessing.html).
<!-- INCLUDE: We create several helpers in our applet script to manage our workload. One helper you may have seen before is `run_cmd`; we use this function to manage or subprocess calls:-->
<!-- FUNCTION: run_cmd -->
<!-- INCLUDE: Before we can split our workload, we need to know what regions are present in our BAM input file. We handle this initial parsing in the `parse_sam_header_for_region` function:-->
<!-- FUNCTION: parse_sam_header_for_region -->

Once our workload is split and we've started processing, we wait and review the status of each `Pool` worker. Then, we merge and output our results.
```python
# Write results to file
resultfn = inputs['mappings_sorted_bam_name'][0]
resultfn = (
    resultfn[:-4] + '_count.txt'
    if resultfn.endswith(".bam")
    else resultfn + '_count.txt')
with open(resultfn, 'w') as f:
    sum_reads = 0
    for res, reg in izip(results, regions):
        read_count = int(res[0])
        sum_reads += read_count
        f.write("Region {0}: {1}\n".format(reg, read_count))
    f.write("Total reads: {0}".format(sum_reads))

count_file = dxpy.upload_local_file(resultfn)
output = {}
output["count_file"] = dxpy.dxlink(count_file)
return output
```
<!-- INCLUDE: {% include note.html content="The `run_cmd` function returns a tuple containing the stdout, stderr, and exit code of the subprocess call. We parse these outputs from our workers to determine whether the run failed or passed." %} -->
<!-- FUNCTION: verify_pool_status -->
