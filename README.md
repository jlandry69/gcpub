GCpipeline
==========

GCpipeline is a suite of tools taking as input illumina sequencing data (as bcl files) and producing fastq files as well as quality control metrics and plots. 
The different steps are:
* bcl to fastq convertion
* barcode demultiplexing (in-line or illumina multiplex scheme)
* alignment to the genome of interest
* quality control metrics and quality control plots

Additionally, we provide an archiving script which, first, cleans up the processed data, second writes a tarball of the run folder and finally archives QC related files (FastQC reports and a minimal set of files needed to use the Illumina Sequencing Analysis viewer).  

Installation
------------

You can build the GCpipeline from source. Dependencies are included as submodules so you need to do a recursive clone. 

`git clone --recursive git@github.com:jlandry69/gcpub.git`

`cd gcpub/`

`make all`

The config file (SampleSheetConverter.cfg) defines the different tools/volumes paths used by the pipeline.

Running GCpipeline
------------------

This is an example how to submit the pipeline (pipeline.sh) to a cluster using [LSF](http://www-03.ibm.com/systems/platformcomputing/products/lsf/processmanager.html title="LSF") manager. It takes one required argument: the run_folder that you want to process and the number of mismatch allowed during demultiplexing as one optional argument (default is 0).

Usage: $0 run_folder [N (allow N mismatches)]

Please find below an example of the bsub command to start the pipeline on a LSF cluster:
```
bsub -q ${QUEUE_NAME} -n ${NUMBER_CPU} -M ${MEMORY_REQUIREMENT} -J ${JOB_NAME}[1-8] -o pipeline_output.txt -e pipeline_error.txt path_to_pipeline.sh path_to_run_folder [N - number of mismatches allowed]"
```
Running Archiving
-----------------

This is an example how to submit the archiving script (archiving.sh) to a cluster using LSF manager. It takes one required argument: the run_folder that you want to archive, and the destination folder as one optional argument (default is current directory).

Usage: $0 run_folder [dest_folder]

PLease find below an example of the bsub command to start the archiving on a LSF cluster:
```
bsub -q ${QUEUE_NAME} -M ${MEMORY_REQUIREMENT} -J ${JOB_NAME} -o archive_output.txt -e archive_error.txt path_to_archive.sh path_to_run_folder [destination_folder]"
```
