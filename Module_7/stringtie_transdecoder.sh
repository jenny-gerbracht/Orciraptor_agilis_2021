#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd

# Define paths to read and working directory locations
mydir="/path/to/wd"
moduledir="${mydir}/Module_7"

# Log file
# Get date
date=$(date "+%Y-%m-%d")

# Define log file and redirect stdout and stderr to this file
if [ ! -d "${moduledir}/Logs/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/Logs/"
  mkdir ${moduledir}/Logs/
fi
log_file="${moduledir}/Logs/log_transdecoder_$date"
exec &> >(tee -a "$log_file")

# Create necessary folders
if [ ! -d "${moduledir}/transdecoder/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/transdecoder/"
  mkdir ${moduledir}/transdecoder/
fi

cd transdecoder

TransDecoder.LongOrfs -t ../SuperDuper_stringtie.fasta
TransDecoder.Predict -t ../SuperDuper_stringtie.fasta
