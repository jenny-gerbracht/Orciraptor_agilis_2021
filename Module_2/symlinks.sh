#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd

# Define paths to read and working directory locations
readdir="/path/to/reads"
mydir="/path/to/wd"
moduledir="${mydir}/Module_2"

# Define paths to scripts
Rcorrector_dir="/scratch2/software/anaconda/envs/rcorrector/bin/"
Discard_dir="/scratch4/shess/Jenny/Scripts/"
	
# Log file
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

# Create necessary folders
if [ ! -d "${moduledir}/readprocessing/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/"
  mkdir ${moduledir}/readprocessing/
fi

if [ ! -d "${moduledir}/readprocessing/corrected/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/corrected/"
  mkdir ${moduledir}/readprocessing/corrected/
fi

if [ ! -d "${moduledir}/readprocessing/fastqc/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/fastqc/"
  mkdir ${moduledir}/readprocessing/fastqc/
fi

if [ ! -d "${moduledir}/readprocessing/fastqc/raw/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/fastqc/raw/"
  mkdir ${moduledir}/readprocessing/fastqc/raw/
fi

if [ ! -d "${moduledir}/readprocessing/fastqc/processed/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/fastqc/processed"
  mkdir ${moduledir}/readprocessing/fastqc/processed/
fi

# Declare samples array
declare -a samples

# Load experiment file, parse for samples and save them into the array
let p=0
while read -r f1 f2; do
	samples[p]="${f1}"
	((++p))	
done < $experiment_file

# Iterate over FASTQ files

for i in "${samples[@]}"; do		
			echo "Concatenate ${i}"
			if [ -e "${moduledir}/readprocessing/${i}_1.fq.gz" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
    		echo "$mytime FASTQ file link already created"
   	else
		ln -s ${readdir}/*${i}_R1.fastq.gz ${moduledir}/readprocessing/${i}_1.fq.gz
	fi			
			if [ -e "${moduledir}/readprocessing/${i}_1.fq.gz" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
    		echo "$mytime FASTQ file link already created"
   	else
		ln -s ${readdir}/*${i}_R2.fastq.gz ${moduledir}/readprocessing/${i}_2.fq.gz
	fi


done	
