#! /usr/bin/env bash

if [ $# -eq 0 ] || [ $# -gt 2 ]; then
    echo "Usage: $0 <run_folder> [dest_folder]"
    exit 1
fi

# User parameter, the run folder
FOLDER=$1

# Force mount
ls ${FOLDER} >/dev/null

# Check if private config file exists
if [ -f ${BASEDIR}/../gcpriv/SampleSheetConverter.cfg ]
then
    CFILE=${BASEDIR}/../gcpriv/SampleSheetConverter.cfg
else
    CFILE=${BASEDIR}/SampleSheetConverter.cfg

# Paths
ARCVOL=`grep "^archiving_path[^A-Za-z]" ${CFILE} | sed 's/^.*=[ \t]*//'`
QCVOL=`grep "^fqc_path[^A-Za-z]" ${CFILE} | sed 's/^.*=[ \t]*//'`
DBCMD=`grep "^db_cmd[^A-Za-z]" ${CFILE} | sed 's/^.*=[ \t]*//'`
sendEmail=`grep "^send_email[^A-Za-z]" ${CFILE} | sed 's/^.*=[ \t]*//'`
# optionally pass a directory as 2nd arg to force writing the tarball to that directory
if [ $# -eq 2 ]; then
    DESTDIR=$2
fi

if [ -f ${FOLDER}/RunInfo.xml ]
    then
    RUN=`echo ${FOLDER} | sed 's/\/$//' |  sed 's/^.*\///'`
    DIRRUN=`echo ${FOLDER} | sed 's/\/Runs\/.*$//'`
    SUCCESS=0

    if [ -d ${DIRRUN} ]
	then
	ARCHIVESEARCH=( ${DESTDIR} ${ARCVOL} )
	for ((  i = 0 ;  i < ${#ARCHIVESEARCH[@]};  i++  ))
	  do
	  LOC=${ARCHIVESEARCH[$i]}
	  if [ -d ${LOC} ]
	      then
	      ls ${LOC} > /dev/null
	      SIZE=`df -B T | sed 's/[ \t][ \t]*/\t/g' | sed 's/^[ \t]*//' | cut -f 3,5 | grep "${LOC}" | cut -f 1 | sed 's/T$//'`
	      if [ ${SIZE} -gt 5 ]
		  then
		  # Keep the statistics
		  if [ -f ${FOLDER}/SampleSheetOriginal.csv ]
		      then
		      cd ${FOLDER}/../
		      tar -czf ${QC_VOL}/${RUN}.tar.gz ${RUN}/Aligned_lane*/summary_stats ${RUN}/Aligned_lane*/Project_*/Summary_Stats_*/Sample_Summary.htm ${RUN}/SampleSheetOriginal.* ${RUN}/InterOp/ ${RUN}/RunInfo.xml ${RUN}/runParameters.xml ${RUN}/Unaligned*/Basecall_Stats_*/Demultiplex_Stats.htm ${RUN}/Aligned_lane*/Flowcell_Summary_*.htm ${RUN}/Unaligned*/Basecall_Stats_*/IVC.htm ${RUN}/Unaligned*/Basecall_Stats_*/Plots/s_[1-8]_[A-Za-z]*.png
		      chmod go-wx ${QC_VOL}/${RUN}.tar.gz
		  fi

		  # Clean-up processed data
		  cd ${FOLDER}
		  #rm -rf Aligned* Unaligned* tobias_* config_*.txt
		  rm -rf Unaligned* tobias_* config_*.txt log_pipeline
		  rm -rf Data/Intensities/L001/C*/*.cif
		  rm -rf Data/Intensities/L002/C*/*.cif
		  rm -rf Data/Intensities/L003/C*/*.cif
		  rm -rf Data/Intensities/L004/C*/*.cif
		  rm -rf Data/Intensities/L005/C*/*.cif
		  rm -rf Data/Intensities/L006/C*/*.cif
		  rm -rf Data/Intensities/L007/C*/*.cif
		  rm -rf Data/Intensities/L008/C*/*.cif
		  cd ..
		  
		  # Upload Run metrics info to db
                  ${DBCMD}
		  
		  # Tar.gz run folder
		  echo ${LOC}/${RUN}.tar.gz > ${RUN}.email.txt
		  tar -czf ${LOC}/${RUN}.tar.gz ${RUN}/ 2>> ${RUN}.email.txt

		  ERROR=$?	      
		  if [ ${ERROR} -ne 0 ]
		      then
		      echo "Error! Exit code ${ERROR}."  >> ${RUN}.email.txt 
		  else
		      SUCCESS=1
		      echo ${RUN} `date +'%Y%m%d'` `du -hs ${LOC}/${RUN}.tar.gz | cut -f 1` >> ${DIRRUN}/archive_inventory.txt
		      echo " " >> ${RUN}.email.txt
		  fi
		  break
	      fi
	  fi
	done
    fi
    if [ ${SUCCESS} -ne 1 ]
	then
	echo "Archiving was NOT successful."  >> ${RUN}.email.txt 
	echo "Flowcell:" ${RUN} >> ${RUN}.email.txt
	echo " " >> ${RUN}.email.txt
        SUBJECT="WARNING: archiving not successful" 
	${sendEmail} "${SUBJECT}" ${RUN}.email.txt
    else
        SUBJECT="Please archive"
	${sendEmail} "${SUBJECT}" ${RUN}.email.txt 1
    fi
    rm ${RUN}.email.txt
fi

