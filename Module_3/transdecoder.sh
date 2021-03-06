#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd

###################################################################
#Script Name	:transdecoder.sh		                                                                                              
#Description	:Predicts ORFs for Orciraptor agilis de novo transcriptome                                                                               
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

source ../config.txt
moduledir="${mydir}/Module_3"

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
log_file="${moduledir}/Logs/log_transdecoder_$date"
exec &> >(tee -a "$log_file")

####################################
#
# Make folders
#
####################################

if [ ! -d "${moduledir}/transdecoder/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/transdecoder/"
  mkdir ${moduledir}/transdecoder/
fi

cd transdecoder

source activate transdecoder

TransDecoder.LongOrfs -t ${moduledir}/orciraptor_200_filtered.fasta
TransDecoder.Predict -t ${moduledir}/orciraptor_200_filtered.fasta

source deactivate
