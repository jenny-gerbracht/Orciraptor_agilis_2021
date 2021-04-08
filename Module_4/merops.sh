#!/bin/bash

# Define paths to working directory locations
mydir="/path/to/wd"
moduledir="${mydir}/Module_4"

# Log file
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

# Create necessary folders
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
--db /srv/Jenny/MEROPS/merops.dmnd \
--out ${moduledir}/merops/merops.out \
--outfmt 6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore \
--max-target-seqs 1
