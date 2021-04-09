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
log_file="${moduledir}/Logs/log_diamond_$date"
exec &> >(tee -a "$log_file")

# Create necessary folders
if [ ! -d "${moduledir}/diamond/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/diamond/"
  mkdir ${moduledir}/diamond/
fi

echo ""
echo "###################"
echo "## Diamond nr search"
echo "###################"
echo ""

diamond \
blastp \
--query ${mydir}/Module_3/transdecoder/orciraptor_transdecoder.pep_renamed.fasta \
--threads 10 \
--db /srv/Jenny/entap/bin/nr \
--out ${moduledir}/diamond/diamond_nr.out \
--outfmt 6 qseqid sseqid pident length evalue stitle \
--max-target-seqs 1
