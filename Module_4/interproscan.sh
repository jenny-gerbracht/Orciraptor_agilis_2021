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
log_file="${moduledir}/Logs/log_interproscan_$date"
exec &> >(tee -a "$log_file")

####################################
#
# Make folders
#
####################################

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
${Interproscan_dir}/interproscan.sh \
--input ${mydir}/Module_3/transdecoder/orciraptor_transdecoder.pep_renamed.fasta \
--iprlookup -pa \
--output-dir ${moduledir}/interproscan \
--goterms \
--cpu 8
