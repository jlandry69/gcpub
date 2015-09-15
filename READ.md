GCpipeline
==========

GCpipeline is a suite of tools taking as input illumina sequencing data (as bcl files). 
The different steps are:
* bcl to fastq convertion
* barcode demultiplexing (in-line or illumina multiplex scheme)
* alignment to the genome of interest
* quality control metrics and quality control plot

Installation
------------

You can build the GCpipeline from source. Dependencies are included as submodules so you need to do a recursive clone. 

`git clone --recursive https://github.com/tobiasrausch/gcpub.git`

`cd gcpub/`

`make all`

The config file (SampleSheetConverter.cfg) defines the different tool paths used by the pipeline.

Running GCpipeline
------------------

This is an example how to submit the pipeline to cluster using LSF manager:

#! /usr/bin/env bash

if [ $# -eq 0 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 <run_folder> [N (allow N mismatches)]"
    exit 1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
MEM=memory_requirement (We use 48Gb)
QUEUE=queue_name
CPU=number_of_CPUs (We use 8)
d=$(date +'%a-%y%m%d-%H%M')

SL=`grep "^submit_log[^A-Za-z]" $1/SampleSheet_lane${LANE}.cfg | sed 's/^.*=[ \t]*//'`

bsub_cmd="bsub -q ${QUEUE} -n ${CPU} -R \"span[hosts=1] select[mem>${MEM}]\" -M ${MEM} -J gcpip[1-8] -o ${SL}/pipeline_output_$d-PID$$ -e ${SL}/pipeline_error_$d-PID$$ -W 9000:00 ${BASEDIR}/pipeline.sh $1 $2"

echo $bsub_cmd | ssh solexa@submaster $(< /dev/fd/0)

