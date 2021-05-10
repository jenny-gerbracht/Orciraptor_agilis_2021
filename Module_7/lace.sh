#!/bin/bash

source ../config.txt
moduledir="${mydir}/Module_7"

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
log_file="${moduledir}/Logs/log_lace_$date"
exec &> >(tee -a "$log_file")

# Generate SuperTranscripts with Lace
echo ""
echo "###################"
echo "## Lace"
echo "###################"
echo ""

${Lace_dir}/Lace_run.py \
--cores 10 \
--alternate \
${mydir}/Module_3/orciraptor_200_filtered.fasta \
${mydir}/Module_5/gene_trans_map
