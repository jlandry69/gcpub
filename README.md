GCpipeline
==========

GCpipeline is a suite of tools taking as input illumina sequencing data (as bcl files) and producing fastq files as well as quality control metrics and plots. 
The different steps are:
* bcl to fastq convertion
* barcode demultiplexing (in-line or illumina multiplex scheme)
* alignment to the genome of interest
* quality control metrics and quality control plots

Additionally, 

Installation
------------

You can build the GCpipeline from source. Dependencies are included as submodules so you need to do a recursive clone. 

`git clone --recursive https://github.com/tobiasrausch/gcpub.git`

`cd gcpub/`

`make all`

The config file (SampleSheetConverter.cfg) defines the different tool paths used by the pipeline.

Running GCpipeline
------------------

This is an example how to submit the pipeline (pipeline.sh) to a cluster using LSF manager. It takes one required argument: the run_folder that you want to process and the number of mismatch allowed during demultiplexing as one optional argument (default is 0).
Usage: $0 run_folder [N (allow N mismatches)]

```
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
MEM=memory_requirement (We use 48Gb)
QUEUE=queue_name
CPU=number_of_CPUs (We use 8)
MACHINE=name_of_the_submission_machine
d=$(date +'%a-%y%m%d-%H%M')

SL=`grep "^submit_log[^A-Za-z]" $1/SampleSheetConverter.cfg | sed 's/^.*=[ \t]*//'`

bsub_cmd="bsub -q ${QUEUE} -n ${CPU} -R \"span[hosts=1] select[mem>${MEM}]\" -M ${MEM} -J gcpip[1-8] -o ${SL}/pipeline_output_$d-PID$$ -e ${SL}/pipeline_error_$d-PID$$ -W 9000:00 ${BASEDIR}/pipeline.sh $1 $2"

echo $bsub_cmd | ssh${MACHINE} $(< /dev/fd/0)
```
Running Archiving
-----------------

This is an example how to submit the archiving script (archiving.sh) to a cluster using LSF manager. It takes one required argument: the run_folder that you want to archive, and the destination folder as  one optional argument (default is current directory).
Usage: $0 run_folder [dest_folder]

```
BASEDIR=$(dirname "$SCRIPT")
MEM=memory_requirement (We use 2Gb)
MACHINE=name_of_the_submission_machine
SL=`grep "^submit_log[^A-Za-z]" $1/SampleSheetConverter.cfg | sed 's/^.*=[ \t]*//'`
d=$(date +'%a-%y%m%d-%H%M')
bsub_cmd="bsub -M ${MEM} -J archive -o ${SL}/archive_output.$d-PID$$ -e ${SL}/archive_error.$d-PID$$ -W 9000:00 ${BASEDIR}/archive.sh $1 $2"
echo $bsub_cmd | ssh ${MACHINE} $(< /dev/fd/0)
```

