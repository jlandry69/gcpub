#! /usr/bin/env bash

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "*************************************************************************"
    echo "GeneCore Pipeline"
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "Version: 0.0.1"
    echo "Contact:"
    echo "Jonathan Landry (landry@embl.de)"
    echo "Markus Hsi-Yang Fritz (fritz@embl.de)"
    echo "Tobias Rausch (rausch@embl.de)"
    echo "*************************************************************************"
    echo ""
    echo "Usage: $0 <run_folder> [N (allow N mismatches)]"
    echo ""
    exit 1
fi

# Pipeline base directory
SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Check if private config file exists
if [ -f ${BASEDIR}/../gcpriv/SampleSheetConverter.cfg ]
then
    CFILE=${BASEDIR}/../gcpriv/SampleSheetConverter.cfg
else
    CFILE=${BASEDIR}/SampleSheetConverter.cfg

# Tools
SAMTOOLS=${BASEDIR}/src/samtools/samtools
BWA=${BASEDIR}/src/bwa/bwa
BCSPLIT=${BASEDIR}/src/bcsplit
PLANE=${BASEDIR}/lane.sh
CONVERT=`grep "^convert[^A-Za-z]" ${CFILE} | sed 's/^.*=[ \t]*//'`
CONVERTST=`grep "^convertst[^A-Za-z]" ${CFILE} | sed 's/^.*=[ \t]*//'`
MD5SUM=md5sum

# Log folder and file prefix
LOGF="log_pipeline"
LOGPRE="gc_pipeline"

# Python path
PY=`grep "^python[^A-Za-z]" ${CFILE} | sed 's/^.*=[ \t]*//'`

# LSF does not provide this. might be able to parse from LSB_MCPU_HOSTS
# if hardcoded like below, should match bsub -n (see above)
NCPUS=8

# User parameters: the run folder and lane
FOLDER=${1}
MISMATCH=0
if [ $# -eq 2 ]
then
    MISMATCH=${2}
    if [ ${2} != '1' ]; then
    	echo "Argument 2, if specified, has to be '1'"
    	exit 1
    fi
fi

# Check if LSB_JOBINDEX is set
if [ `echo ${LSB_JOBINDEX} | grep "." | wc -l | cut -f 1` -eq 1 ]
then
    LANE=${LSB_JOBINDEX}
else
    LANE=1
fi

# Check for the required sample sheet
if [ -f ${FOLDER}/SampleSheetOriginal.csv ]
then
    cd ${FOLDER}

    # Print out LSB_HOSTS in log file
    mkdir -p ${FOLDER}/${LOGF}
    echo "Job running on: ${LSB_HOSTS}" > ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt

    # Is lane present?
    if [ `cat ${FOLDER}/SampleSheetOriginal.csv | sed 's/,/\t/g' | awk '$2=="'${LANE}'"' 2> /dev/null | wc -l | cut -f 1` -eq 0 ]
    then
	echo "Lane: ${LANE}, Finished:" `date` >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
	exit;
    fi

    # Task1: SampleSheetConverter
    echo "SampleSheetConverter started." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
    ${PY} ${BASEDIR}/SampleSheetConverter.py --run ${FOLDER} --lane ${LANE} --config ${CFILE} >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt 2>&1
    if [ $? -ne 0 ]
    then
        echo "SampleSheetConverter failed." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
	exit 1
    fi
    echo "SampleSheetConverter done." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
    
    #Get configurations
    PAIRED=`grep "^paired[^A-Za-z]" ${FOLDER}/SampleSheet_lane${LANE}.cfg | sed 's/^.*=[ \t]*//'`
    USE_BASES=`grep "^use_bases[^A-Za-z]" ${FOLDER}/SampleSheet_lane${LANE}.cfg | sed 's/^.*=[ \t]*//'`
    INSTRUMENT=`grep "^instrument[^A-Za-z]" ${FOLDER}/SampleSheet_lane${LANE}.cfg | sed 's/^.*=[ \t]*//'`
    FLOWCELL=`grep "^flowcell[^A-Za-z]" ${FOLDER}/SampleSheet_lane${LANE}.cfg | sed 's/^.*=[ \t]*//'`
    BARCODE=`grep "^barcode_mode[^A-Za-z]" ${FOLDER}/SampleSheet_lane${LANE}.cfg | sed 's/^.*=[ \t]*//'`
    echo "Instrument: ${INSTRUMENT}, Flowcell: ${FLOWCELL}, Configuration: ${PAIRED}, ${USE_BASES}, ${BARCODE}" >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt

    # Task2: Do basecall conversion
    echo "Basecall conversion started." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
    if [ `echo ${INSTRUMENT} | egrep "^NS|^ST" | wc -l | cut -f 1` -ne 0 ]
    then
	if [ `echo ${INSTRUMENT} | grep "^NS" | wc -l | cut -f 1` -ne 0 ]
	then
            ${CONVERTST} --use-bases-mask ${USE_BASES} --barcode-mismatches ${MISMATCH} --runfolder-dir ${FOLDER} --output-dir ${FOLDER}/Unaligned_lane${LANE} > ${FOLDER}/Unaligned_lane${LANE}/nohup.out 2>&1
	else
	    ${CONVERTST} --tiles s_${LANE} --use-bases-mask ${USE_BASES} --barcode-mismatches ${MISMATCH} --runfolder-dir ${FOLDER} --output-dir ${FOLDER}/Unaligned_lane${LANE} --sample-sheet ${FOLDER}/SampleSheet_lane${LANE}.csv  > ${FOLDER}/Unaligned_lane${LANE}/nohup.out 2>&1
	fi

	if [ $? -ne 0 ]
	then
    	    echo "Basecall conversion failed." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
            exit 1
    	fi

	# Test if undetermined files exist
	ls -1 Unaligned_lane${LANE}/Undetermined_*.gz > /dev/null 2>&1
	if [ $? -eq 0 ]
	then
	    mv ${FOLDER}/Unaligned_lane${LANE}/Undetermined_*.gz ${FOLDER}/Unaligned_lane${LANE}/Undetermined_indices/
	fi
    else
	${CONVERT} --tiles s_${LANE} --force --use-bases-mask ${USE_BASES} --mismatches ${MISMATCH} --input-dir ${FOLDER}/Data/Intensities/BaseCalls/ --output-dir ${FOLDER}/Unaligned_lane${LANE} --sample-sheet ${FOLDER}/SampleSheet_lane${LANE}.csv --no-eamss &> ${FOLDER}/${LOGF}/${LOGPRE}.conversion_init.lane${LANE}.txt
    	cd ${FOLDER}/Unaligned_lane${LANE}    	
    	make -j ${NCPUS} > nohup.out 2>&1
    fi

    # Move files if necessary
    cat ${FOLDER}/SampleSheet_lane${LANE}.csv | sed -n -e '/SampleID/,$p' | grep "lane" | sed 's/^lane/,,lane/' | sed 's/,/\t/g' | cut -f 3,5 | sort | uniq > ${FOLDER}/${LOGF}/lane${LANE}_barcodes.txt
    for SID in `cut -f 1 ${FOLDER}/${LOGF}/lane${LANE}_barcodes.txt`
    do
  	if [ ! -d ${FOLDER}/Unaligned_lane${LANE}/Project_${FLOWCELL}/Sample_${SID} ]
	then
	    mkdir ${FOLDER}/Unaligned_lane${LANE}/Project_${FLOWCELL}/Sample_${SID}
	fi
	ls -1 ${FOLDER}/Unaligned_lane${LANE}/${SID}_* > /dev/null 2>&1
	if [ $? -eq 0 ]
	then 
	    mv ${FOLDER}/Unaligned_lane${LANE}/${SID}_* ${FOLDER}/Unaligned_lane${LANE}/Project_${FLOWCELL}/Sample_${SID}/
	fi
    done
    echo "Basecalling is done." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt

    # Task3: Merge FASTQ files
    if [ -d ${FOLDER}/Unaligned_lane${LANE}/Project_${FLOWCELL} ]
    then
	cd ${FOLDER}/Unaligned_lane${LANE}
	# Loop over Samples
	echo "Merging files is started." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
	for SAMPLEID in `cut -f 1 ${FOLDER}/${LOGF}/lane${LANE}_barcodes.txt`
	do  	
	    if [ ${PAIRED} == "PE" ]
	    then
	  	zcat ${FOLDER}/Unaligned_lane${LANE}/Project_${FLOWCELL}/Sample_${SAMPLEID}/*_R1_[0-9][0-9][0-9].fastq.gz | grep -A 3 -P '^@.* [^:]*:N:[^:]*:' | sed 's/^--$//' | sed '/^$/d' | gzip -c > ${SAMPLEID}_1_sequence.txt.gz
	  	zcat ${FOLDER}/Unaligned_lane${LANE}/Project_${FLOWCELL}/Sample_${SAMPLEID}/*_R2_[0-9][0-9][0-9].fastq.gz | grep -A 3 -P '^@.* [^:]*:N:[^:]*:' | sed 's/^--$//' | sed '/^$/d' | gzip -c > ${SAMPLEID}_2_sequence.txt.gz
	    else
		zcat ${FOLDER}/Unaligned_lane${LANE}/Project_${FLOWCELL}/Sample_${SAMPLEID}/*_[0-9][0-9][0-9].fastq.gz | grep -A 3 -P '^@.* [^:]*:N:[^:]*:' | sed 's/^--$//' | sed '/^$/d' | gzip -c > ${SAMPLEID}_sequence.txt.gz
    	    fi
	done
    	echo "Merging files is done." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
    fi

    # Task4: Split custom barcodes
    if [ ${BARCODE} == 'custom' ]
    then
	echo "Custom barcodes demultiplexing started." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
	cat ${FOLDER}/SampleSheetOriginal.csv | grep ",custom," | sed 's/,/\t/g'  | awk '$2=="'${LANE}'"' | cut -f 3,5 | sort | uniq > ${FOLDER}/${LOGF}/lane${LANE}_barcodes.txt
	if [ ${PAIRED} == "PE" ]
	then
	    ${BCSPLIT} pe -m ${MISMATCH} -o lane${LANE} ${FOLDER}/${LOGF}/lane${LANE}_barcodes.txt lane${LANE}_1_sequence.txt.gz lane${LANE}_2_sequence.txt.gz > ${FOLDER}/${LOGF}/lane${LANE}.custom.barcode_dist.txt 2>> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
	else
	    ${BCSPLIT} se -m ${MISMATCH} -o lane${LANE} ${FOLDER}/${LOGF}/lane${LANE}_barcodes.txt lane${LANE}_sequence.txt.gz > ${FOLDER}/${LOGF}/lane${LANE}.custom.barcode_dist.txt 2>> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
	fi
	if [ $? -ne 0 ]
	then
            echo "Custom barcodes demultiplexing failed." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
            exit 1
        else
	    echo "Custom barcodes demultiplexing done." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
	fi	      
    fi

    # Task5: Move and rename fastq from Unaligned_lane folder to Aligned_lane
    echo "Move from Unaligned_lane${LANE} folder to Aligned_lane${LANE} folder started." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
    SAMPLEARR=( `cut -f 1 ${FOLDER}/${LOGF}/lane${LANE}_barcodes.txt | tr '\n' ' '` )
    for ((  i = 0 ;  i < ${#SAMPLEARR[@]};  i++  ))
    do
	if [ ${PAIRED} == "SE" ]
	then
	    if [ -f ${FOLDER}/Unaligned_lane${LANE}/*_${SAMPLEARR[$i]}.fq.gz ]
   	    then
		mv ${FOLDER}/Unaligned_lane${LANE}/*_${SAMPLEARR[$i]}.fq.gz ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_sequence.txt.gz
		${MD5SUM} ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_sequence.txt.gz > ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_sequence.txt.gz.md5
	    elif [ -f ${FOLDER}/Unaligned_lane${LANE}/${SAMPLEARR[$i]}_sequence.txt.gz ]
	    then
		mv ${FOLDER}/Unaligned_lane${LANE}/${SAMPLEARR[$i]}_sequence.txt.gz ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_sequence.txt.gz
		${MD5SUM} ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_sequence.txt.gz > ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_sequence.txt.gz.md5
	    fi
	    # For custom barcoded lanes we also need to move the original unsplitted fq file
	    if [ -f ${FOLDER}/Unaligned_lane${LANE}/lane${LANE}_sequence.txt.gz ]
	    then
		mv ${FOLDER}/Unaligned_lane${LANE}/lane${LANE}_sequence.txt.gz ${FOLDER}/Aligned_lane${LANE}/lane${LANE}_sequence.txt.gz
		${MD5SUM} ${FOLDER}/Aligned_lane${LANE}/lane${LANE}_sequence.txt.gz > ${FOLDER}/Aligned_lane${LANE}/lane${LANE}_sequence.txt.gz.md5
	    fi
	elif [ ${PAIRED} == "PE" ]
	then
	    if [ -f ${FOLDER}/Unaligned_lane${LANE}/*_1_${SAMPLEARR[$i]}.fq.gz ]
	    then
		mv ${FOLDER}/Unaligned_lane${LANE}/*_1_${SAMPLEARR[$i]}.fq.gz ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_1_sequence.txt.gz
		${MD5SUM} ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_1_sequence.txt.gz > ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_1_sequence.txt.gz.md5
	    elif [ -f ${FOLDER}/Unaligned_lane${LANE}/${SAMPLEARR[$i]}_1_sequence.txt.gz ]
	    then
		mv ${FOLDER}/Unaligned_lane${LANE}/${SAMPLEARR[$i]}_1_sequence.txt.gz ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_1_sequence.txt.gz
		${MD5SUM} ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_1_sequence.txt.gz > ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_1_sequence.txt.gz.md5
	    fi
	    if [ -f ${FOLDER}/Unaligned_lane${LANE}/*_2_${SAMPLEARR[$i]}.fq.gz ]
	    then
		mv ${FOLDER}/Unaligned_lane${LANE}/*_2_${SAMPLEARR[$i]}.fq.gz ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_2_sequence.txt.gz
		${MD5SUM} ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_2_sequence.txt.gz > ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_2_sequence.txt.gz.md5
	    elif [ -f ${FOLDER}/Unaligned_lane${LANE}/${SAMPLEARR[$i]}_2_sequence.txt.gz ]
	    then
		mv ${FOLDER}/Unaligned_lane${LANE}/${SAMPLEARR[$i]}_2_sequence.txt.gz ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_2_sequence.txt.gz
		${MD5SUM} ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_2_sequence.txt.gz > ${FOLDER}/Aligned_lane${LANE}/${SAMPLEARR[$i]}_2_sequence.txt.gz.md5
	    fi
	    # For custom barcoded lanes we also need to move the original unsplitted fq file
	    if [ -f ${FOLDER}/Unaligned_lane${LANE}/lane${LANE}_1_sequence.txt.gz ]
	    then
		mv ${FOLDER}/Unaligned_lane${LANE}/lane${LANE}_1_sequence.txt.gz ${FOLDER}/Aligned_lane${LANE}/lane${LANE}_1_sequence.txt.gz
		${MD5SUM} ${FOLDER}/Aligned_lane${LANE}/lane${LANE}_1_sequence.txt.gz > ${FOLDER}/Aligned_lane${LANE}/lane${LANE}_1_sequence.txt.gz.md5
	    fi
	    if [ -f ${FOLDER}/Unaligned_lane${LANE}/lane${LANE}_2_sequence.txt.gz ]
	    then
               	mv ${FOLDER}/Unaligned_lane${LANE}/lane${LANE}_2_sequence.txt.gz ${FOLDER}/Aligned_lane${LANE}/lane${LANE}_2_sequence.txt.gz
		${MD5SUM} ${FOLDER}/Aligned_lane${LANE}/lane${LANE}_2_sequence.txt.gz > ${FOLDER}/Aligned_lane${LANE}/lane${LANE}_2_sequence.txt.gz.md5
            fi
	fi
    done
    echo "Move from Unaligned_lane${LANE} folder to Aligned_lane${LANE} folder done." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt

    # Task6: BWA Alignment
    GENOMES=`grep "^genomes[^A-Za-z]" ${CFILE} | sed 's/^.*=[ \t]*//'`
    cat ${FOLDER}/SampleSheetOriginal.csv | sed 's/,/\t/g' | awk '$2=="'${LANE}'"' | cut -f 3,4,5 > ${FOLDER}/${LOGF}/lane${LANE}_barcodes.txt
    while read SAMPLEID REF BARCODE
    do 
	if [ `echo ${REF} | egrep -w "na|NA|unknown|Unknown" | wc -l | cut -f 1` -eq 0 ]
    	then
	    echo "Alignment started." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
	    REFGEN=${GENOMES}/${REF}/${REF}.fa
	
	    if [ ${PAIRED} == "SE" ]
    	    then
		${BWA} mem -t ${NCPUS} ${REFGEN} ${FOLDER}/Aligned_lane${LANE}/${SAMPLEID}_sequence.txt.gz 2>> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.bwa.err | ${SAMTOOLS} view -bT ${REFGEN} - > ${FOLDER}/Aligned_lane${LANE}/${SAMPLEID}_sequence.bam 2>> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.samtools.err
	    else
		${BWA} mem -t ${NCPUS} ${REFGEN} ${FOLDER}/Aligned_lane${LANE}/${SAMPLEID}_1_sequence.txt.gz ${FOLDER}/Aligned_lane${LANE}/${SAMPLEID}_2_sequence.txt.gz 2>> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.bwa.err | ${SAMTOOLS} view -bT ${REFGEN} - > ${FOLDER}/Aligned_lane${LANE}/${SAMPLEID}_sequence.bam 2>> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.samtools.err
	    fi
	    ${SAMTOOLS} sort ${FOLDER}/Aligned_lane${LANE}/${SAMPLEID}_sequence.bam ${FOLDER}/Aligned_lane${LANE}/${SAMPLEID}_sequence.tmp 2>> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.samtools.err
	    mv ${FOLDER}/Aligned_lane${LANE}/${SAMPLEID}_sequence.tmp.bam ${FOLDER}/Aligned_lane${LANE}/${SAMPLEID}_sequence.bam
	    ${SAMTOOLS} index ${FOLDER}/Aligned_lane${LANE}/${SAMPLEID}_sequence.bam 2>> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.samtools.err
	    ${MD5SUM} ${FOLDER}/Aligned_lane${LANE}/${SAMPLEID}_sequence.bam > ${FOLDER}/Aligned_lane${LANE}/${SAMPLEID}_sequence.bam.md5
	fi
	echo "Alignment done." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
    done < ${FOLDER}/${LOGF}/lane${LANE}_barcodes.txt

    # Task6: Run post lane command
    echo "Post lane analysis started." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
    ${PLANE} ${FOLDER} ${LANE} > ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.postlane.txt 2>&1
    if [ $? -ne 0 ]
    then
    	echo "Post lane analysis failed." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
        exit 1
    fi
    echo "Post lane analysis finished." >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt

    # Finished analysis
    echo "Lane: ${LANE}, Finished:" `date` >> ${FOLDER}/${LOGF}/${LOGPRE}.lane${LANE}.txt
	echo "Lane: ${LANE}, Finished:" `date` >> ${FOLDER}/${LOGF}/${LOGPRE}_finished.txt

fi
