#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd

source ../config.txt
moduledir="${mydir}/Module_2"

####################################
#
# Setting up log file
#
###################################
# Get date
date=$(date "+%Y-%m-%d")

# Define log file and redirect stdout and stderr to this file
if [ ! -d "${moduledir}/Logs/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/Logs/"
  mkdir ${moduledir}/Logs/
fi
log_file="${moduledir}/Logs/log_symlinks_$date"
exec &> >(tee -a "$log_file")

####################################
#
# Setting up samples array
#
####################################
# Declare samples array
declare -a samples

# Load experiment file, parse for samples and save them into the array
let p=0
while read -r f1 f2; do
	samples[p]="${f1}"
	((++p))	
done < $experiment_file

####################################
#
# Make folders
#
####################################

if [ ! -d "${moduledir}/readprocessing/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/"
  mkdir ${moduledir}/readprocessing/
fi

# Iterate over FASTQ files

for i in "${samples[@]}"; do		
			echo "Sample ${i}"
			if [ -e "${moduledir}/readprocessing/${i}_1.fq.gz" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
    		echo "$mytime FASTQ file link already created"
   	else
		ln -s ${readdir}/*${i}_R1.fastq.gz ${moduledir}/readprocessing/${i}_1.fq.gz
	fi			
			if [ -e "${moduledir}/readprocessing/${i}_2.fq.gz" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
    		echo "$mytime FASTQ file link already created"
   	else
		ln -s ${readdir}/*${i}_R2.fastq.gz ${moduledir}/readprocessing/${i}_2.fq.gz
	fi


done	
