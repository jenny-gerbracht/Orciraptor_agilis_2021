#!/bin/bash

# Define paths to working directory locations
mydir="/path/to/wd"
moduledir="${mydir}/Module_4"

# Define paths to scripts
interproscan_dir="/home/jenny/interproscan/interproscan-5.50-84.0"

# Log file
# Get date
date=$(date "+%Y-%m-%d")

# Define log file and redirect stdout and stderr to this file
if [ ! -d "${moduledir}/Logs/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/Logs/"
  mkdir ${moduledir}/Logs/
fi
log_file="${moduledir}/Logs/log_interproscan_$date"
exec &> >(tee -a "$log_file")

# Create necessary folders
if [ ! -d "${moduledir}/interproscan/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/interproscan/"
  mkdir ${moduledir}/interproscan/
fi

echo ""
echo "###################"
echo "## InterProScan"
echo "###################"
echo ""

# Remove asterisks from protein fasta
sed 's/*//g' ${mydir}/Module_3/transdecoder/orciraptor_transdecoder.pep_renamed.fasta \
> ${mydir}/Module_3/transdecoder/orciraptor_transdecoder.pep_renamed_nostop.fasta

# Run analysis
${interproscan_dir}/interproscan.sh \
--input ${mydir}/Module_3/transdecoder/orciraptor_transdecoder.pep_renamed.fasta \
--iprlookup -pa \
--output-dir ${moduledir}/interproscan \
--goterms \
--cpu 8
