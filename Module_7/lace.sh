#!/bin/bash

# Define paths to working directory locations
mydir="/path/to/wd"
moduledir="${mydir}/Module_7"

# Define paths to scripts
lace_dir="/home/jenny/Lace-1.14.1/Lace/"

# Define log file and redirect stdout and stderr to this file
if [ ! -d "${moduledir}/Logs/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/Logs/"
  mkdir ${moduledir}/Logs/
fi
log_file="${moduledir}/Logs/log_lace_$date"
exec &> >(tee -a "$log_file")

# Generate SuperTranscripts with Lace
echo ""
echo "###################"
echo "## Lace"
echo "###################"
echo ""

${lace_dir}/Lace_run.py \
--cores 10 \
--alternate \
${mydir}/Module_3/orciraptor_200_filtered.fasta \
${mydir}/Module_5/gene_trans_map
