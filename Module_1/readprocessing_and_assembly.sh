#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe -threaded 10

# Define paths to read and working directory locations
readdir="/path/to/reads"
mydir="/path/to/wd"
moduledir="${mydir}/Module_1"

# Define paths to scripts
Rcorrector_dir="/scratch2/software/anaconda/envs/rcorrector/bin/"
Discard_dir="/scratch4/shess/Jenny/Scripts/"

# Log file
# Get date
date=$(date "+%Y-%m-%d")

# Define log file and redirect stdout and stderr to this file
if [ ! -d "${mydir}/Logs/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${mydir}/Logs/"
  mkdir ${mydir}/Logs/
fi
log_file="${mydir}/Logs/log_$date"
exec &> >(tee -a "$log_file")

# Create necessary folders
if [ ! -d "${moduledir}/readprocessing/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/"
  mkdir${moduledir}/readprocessing/
fi

if [ ! -d "${moduledir}/readprocessing/corrected/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/corrected/"
  mkdir ${moduledir}/readprocessing/corrected/
fi

# Make symlinks
ln -s $readdir/HI.3499.002.D704---D504.NA_R1.fastq ${moduledir}/readprocessing/NA_R1.fastq
ln -s $readdir/HI.3499.002.D704---D504.NA_R2.fastq ${moduledir}/readprocessing/NA_R2.fastq

# kmer-based error correction with Rcorrector
perl $Rcorrector_dir/run_rcorrector.pl \
  -p ${moduledir}/readprocessing/NA_R1.fastq ${moduledir}/readprocessing/NA_R2.fastq \
  -od ${moduledir}/corrected/ \
  -t 10
 
# Remove unfixable read pairs, compress output
python FilterUncorrectabledPEfastq.py \
-1 ${moduledir}/corrected/NA_R1.cor.fastq \
-2 ${moduledir}/corrected/NA_R2.cor.fastq \
-s NA
pigz unfixrm_NA_R1.cor.fq
pigz unfixrm_NA_R2.cor.fq

