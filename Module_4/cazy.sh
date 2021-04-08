#!/bin/bash

# Define paths to read and working directory locations
mydir="/path/to/wd"
moduledir="${mydir}/Module_4"

# Define paths to scripts
cazy_dir="/srv/Jenny/CAZy_v9/run_dbcan/"

# Log file
# Get date
date=$(date "+%Y-%m-%d")

# Define log file and redirect stdout and stderr to this file
if [ ! -d "${moduledir}/Logs/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/Logs/"
  mkdir ${moduledir}/Logs/
fi
log_file="${moduledir}/Logs/log_cazy_$date"
exec &> >(tee -a "$log_file")

# Annotate carbohydrate-active enzymes with dbcan2 in HMM mode (http://bcb.unl.edu/dbCAN2/download/Tools/)
echo ""
echo "###################"
echo "## dbcan2"
echo "###################"
echo ""

python3 ${cazy_dir}/run_dbcan.py \
../Module_3/transdecoder/orciraptor_transdecoder.pep_renamed.fasta \
protein \
--db_dir ${cazy_dir}/db \
--dbCANFile dbCAN-HMMdb-V9.txt \
--hmm_eval 10 \
--hmm_cpu 10 \
--out_dir cazy \
--tools hmmer 
