#!/bin/bash

if [ $# -ne 2 ]
then
    echo "Usage: $0 <runfolder> <lane>"
    exit -1
fi

FOLDER=$1
LANE=$2

# Pipeline base directory
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Check if private config file exists
if [ -f ${BASEDIR}/../gcpriv/SampleSheetConverter.cfg ]
then
    CFILE=${BASEDIR}/../gcpriv/SampleSheetConverter.cfg
else
    CFILE=${BASEDIR}/SampleSheetConverter.cfg
fi
# Tools
SAMTOOLS=${BASEDIR}/src/samtools/samtools
BWA=${BASEDIR}/src/bwa/bwa
DELLY=${BASEDIR}/src/delly
FASTQC_BIN=`grep "^fastqc[^A-Za-z]" ${CFILE} | sed 's/^.*=[ \t]*//'`
PY=python2.7
INSTRUMENT=`grep "^instrument[^A-Za-z]" ${FOLDER}/SampleSheet_lane${LANE}.cfg | sed 's/^.*=[ \t]*//'`

# Log folder and file prefix
LOGF="log_pipeline"
LOGPRE="gc_pipeline"

# Statistics folder
STATSFOLDER=${FOLDER}/Aligned_lane${LANE}/summary_stats

# Get all samples of this lane and move 
if [ -f ${FOLDER}/${LOGF}/lane${LANE}_barcodes.txt ]
then
    if [ -f ${FOLDER}/SampleSheet_lane${LANE}.cfg ]
    then
		SAMPLEARR=( `cut -f 1 ${FOLDER}/${LOGF}/lane${LANE}_barcodes.txt | tr '\n' ' '` )
		PAIRED=`grep "^paired[^A-Za-z]" ${FOLDER}/SampleSheet_lane${LANE}.cfg | sed 's/^.*=[ \t]*//'`

		# Create stats folder
		mkdir -p ${STATSFOLDER}
		if [ $? -ne 0 ]
		then
	    	echo "ERROR: Cannot create stats folder."
	    	exit 1;
		fi

		# Fastqc & Stats
		cd ${FOLDER}/Aligned_lane${LANE}
		for ((  i = 0 ;  i < ${#SAMPLEARR[@]};  i++  ))
		do
	    	if [ ${PAIRED} == "SE" ]
	    	then
				if [ -f ${SAMPLEARR[$i]}_sequence.txt.gz ]
				then
		    		${FASTQC_BIN} -o ${STATSFOLDER} --noextract ${SAMPLEARR[$i]}_sequence.txt.gz
				fi
	    	elif [ ${PAIRED} == "PE" ]
	    	then
				if [ -f ${SAMPLEARR[$i]}_1_sequence.txt.gz ]
				then
		    		${FASTQC_BIN} -o ${STATSFOLDER} --noextract ${SAMPLEARR[$i]}_1_sequence.txt.gz
				fi
				if [ -f ${SAMPLEARR[$i]}_2_sequence.txt.gz ]
				then
		    		${FASTQC_BIN} -o ${STATSFOLDER} --noextract ${SAMPLEARR[$i]}_2_sequence.txt.gz
				fi
	    	fi
	    	if [ -f ${SAMPLEARR[$i]}_sequence.bam ]
	    	then
				OUTPREFIX=${SAMPLEARR[$i]}
				${DELLY}/src/stats -o ${STATSFOLDER}/${OUTPREFIX}.summary.txt -i ${STATSFOLDER}/${OUTPREFIX}.insert.txt ${SAMPLEARR[$i]}_sequence.bam
				${DELLY}/src/stats -k -o ${STATSFOLDER}/${OUTPREFIX}.summary.cfg ${SAMPLEARR[$i]}_sequence.bam
				Rscript ${DELLY}/R/isize.R ${STATSFOLDER}/${OUTPREFIX}.insert.txt
	    	fi
		done
    fi
fi

# Create barcode summary
if [ -f ${FOLDER}/Unaligned_lane${LANE}/Basecall_Stats_*/Demultiplex_Stats.htm ]
then
    cd ${FOLDER}/Aligned_lane${LANE}
    ${PY} ${BASEDIR}/get_illumina_barcode_dist.py ${FOLDER}/Unaligned_lane${LANE}/Basecall_Stats_*/Demultiplex_Stats.htm > ${STATSFOLDER}/lane${LANE}.illumina.barcode_dist.txt
elif [ -f ${FOLDER}/Unaligned_lane${LANE}/Stats/DemultiplexingStats.xml ]
then
    if [[ "$INSTRUMENT" =~ ^ST ]]
    then
	cd ${FOLDER}/Aligned_lane${LANE}
	${PY} ${BASEDIR}/get_illumina_barcode_dist_ST.py ${FOLDER}/Unaligned_lane${LANE}/Stats/ConversionStats.xml > ${STATSFOLDER}/lane${LANE}.illumina.barcode_dist.txt
    else
	cd ${FOLDER}/Aligned_lane${LANE}
	${PY} ${BASEDIR}/get_illumina_barcode_dist_NS.py ${FOLDER}/Unaligned_lane${LANE}/Stats/DemultiplexingStats.xml > ${STATSFOLDER}/lane${LANE}.illumina.barcode_dist.txt
    fi
fi
if [ -f ${FOLDER}/${LOGF}/lane${LANE}.custom.barcode_dist.txt ]
then
    cp ${FOLDER}/${LOGF}/lane${LANE}.custom.barcode_dist.txt ${STATSFOLDER}/
	rm -f ${STATSFOLDER}/lane${LANE}.illumina.barcode_dist.txt
fi

# Plot barcode distribution
if [ -f ${FOLDER}/Aligned_lane${LANE}/summary_stats/*barcode_dist.txt ]
then
	cd ${FOLDER}/Aligned_lane${LANE}/summary_stats
	for f in *.barcode_dist.txt; do
    		outfn=`basename $f .txt`
    		echo "infile='"$f"'; plotTitle='"lane${LANE}"'; outfile='"$outfn"'" | cat - ${BASEDIR}/R/barcode_dist.R | R --vanilla --slave
    		ps2pdf *.barcode_dist.ps
	done
fi
