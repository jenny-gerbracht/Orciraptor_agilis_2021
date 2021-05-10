#!/bin/bash

source ../config.txt
moduledir="${mydir}/Module_4"

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
  echo "${mytime} Make directory ${moduledir}/Logs/"
  mkdir ${moduledir}/Logs/
fi
log_file="${moduledir}/Logs/log_merops_$date"
exec &> >(tee -a "$log_file")

####################################
#
# Make folders
#
####################################

if [ ! -d "${moduledir}/merops/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/merops/"
  mkdir ${moduledir}/merops/
fi

# Annotation of peptidases with MEROPS database
echo ""
echo "###################"
echo "## MEROPS"
echo "###################"
echo ""

diamond \
blastp \
--sensitive \
--query ${mydir}/Module_3/transdecoder/orciraptor_transdecoder.pep_renamed.fasta \
--threads 10 \
--db /path/to/merops.dmnd \
--out ${moduledir}/merops/merops.out \
--outfmt 6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore \
--max-target-seqs 1
