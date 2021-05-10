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
log_file="${moduledir}/Logs/log_eggnog_$date"
exec &> >(tee -a "$log_file")

####################################
#
# Make folders
#
####################################

if [ ! -d "${moduledir}/eggnog/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/eggnog/"
  mkdir ${moduledir}/eggnog/
fi

# Annotation with eggnog-mapper
echo ""
echo "###################"
echo "## eggnog-mapper"
echo "###################"
echo ""

python ${Eggnog_dir}/emapper.py \
-i ${mydir}/Module_3/transdecoder/orciraptor_transdecoder.pep_renamed.fasta \
--output orciraptor_diamond \
--output_dir ${moduledir}/eggnog/ \
-m diamond \
--report_orthologs \
--pfam_transfer best_og \
--pfam_realign none \
--cpu 10 \
--data_dir /srv/Jenny/eggnog-db/

python ${Eggnog_dir}/emapper.py \
-i ${mydir}/Module_3/transdecoder/orciraptor_transdecoder.pep_renamed.fasta \
--output orciraptor_hmm \
--output_dir ${moduledir}/eggnog/ \
-m hmmer \
-d Eukaryota \
--report_orthologs \
--pfam_transfer best_og \
--pfam_realign none \
--cpu 10 \
--data_dir /srv/Jenny/eggnog-db/
