#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 10

###################################################################
#Script Name	:readprocessing_and_assembly.sh		                                                                                              
#Description	:Performs read processing (Rcorrector and Trim-Galore!), FastQC analysis of raw vs. processed reads,
#		 and de novo assembly of Mougeotia sp.                                                                               
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

source ../config.txt
moduledir="${mydir}/Module_1"

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
log_file="${moduledir}/Logs/log_assembly_$date"
exec &> >(tee -a "$log_file")

####################################
#
# Make folders
#
####################################

if [ ! -d "${moduledir}/readprocessing/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/readprocessing/"
  mkdir ${moduledir}/readprocessing/
fi

if [ ! -d "${moduledir}/readprocessing/corrected/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/readprocessing/corrected/"
  mkdir ${moduledir}/readprocessing/corrected/
fi

if [ ! -d "${moduledir}/readprocessing/fastqc/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/readprocessing/fastqc/"
  mkdir ${moduledir}/readprocessing/fastqc/
fi

if [ ! -d "${moduledir}/readprocessing/fastqc/raw/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/readprocessing/fastqc/raw/"
  mkdir ${moduledir}/readprocessing/fastqc/raw/
fi

if [ ! -d "${moduledir}/readprocessing/fastqc/processed/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/readprocessing/fastqc/processed"
  mkdir ${moduledir}/readprocessing/fastqc/processed/
fi

# Make symlinks
ln -s ${readdir_prey}/HI.3499.002.D704---D504.NA_R1.fastq ${moduledir}/readprocessing/NA_R1.fastq
ln -s ${readdir_prey}/HI.3499.002.D704---D504.NA_R2.fastq ${moduledir}/readprocessing/NA_R1.fastq

# kmer-based error correction with Rcorrector
echo ""
echo "###################"
echo "## Rcorrector"
echo "###################"
echo ""
perl ${Rcorrector_dir}/run_rcorrector.pl \
  -p ${moduledir}/readprocessing/NA_R1.fastq ${moduledir}/readprocessing/NA_R2.fastq \
  -od ${moduledir}/readprocessing/corrected/ \
  -t 10
 
# Remove unfixable read pairs, compress output
echo ""
echo "###################"
echo "## Discard"
echo "###################"
echo ""
cd ${moduledir}/readprocessing/corrected/
python ${Discard_dir}/FilterUncorrectabledPEfastq.py \
-1 ${moduledir}/readprocessing/corrected/NA_R1.cor.fq \
-2 ${moduledir}/readprocessing/corrected/NA_R2.cor.fq \
-s NA
pigz unfixrm_NA_R1.cor.fq
pigz unfixrm_NA_R2.cor.fq
rm NA_R1.cor.fq
rm NA_R2.cor.fq
cd ${moduledir}

# Perform quality and adapter trimming with Trim-Galore!
source activate trim-galore
echo ""
echo "###################"
echo "## Trim-Galore!"
echo "###################"
echo ""
echo -n "Trim-Galore! version: "
trim_galore -v
echo ""
trim_galore \
--fastqc \
-o ${moduledir}/readprocessing/trim_galore \
-j 4 \
--paired \
--length 49 \
${moduledir}/readprocessing/corrected/unfixrm_NA_R1.cor.fq.gz \
${moduledir}/readprocessing/corrected/unfixrm_NA_R2.cor.fq.gz
source deactivate

# Generate folder with fastqc analysis of raw vs. processed reads
source activate fastqc
echo ""
echo "###################"
echo "## FastQC"
echo "###################"
echo ""
echo -n "FastQC version: "
fastqc -v
echo ""
fastqc \
-o ${moduledir}/readprocessing/fastqc/raw/ \
${moduledir}/readprocessing/NA_R1.fastq \
${moduledir}/readprocessing/NA_R2.fastq
source deactivate
mv ${moduledir}/readprocessing/trim_galore/unfixrm_NA_R1.cor_val_1_fastqc* ${moduledir}/readprocessing/fastqc/processed/
mv ${moduledir}/readprocessing/trim_galore/unfixrm_NA_R2.cor_val_2_fastqc* ${moduledir}/readprocessing/fastqc/processed/

# Assembly with rnaSPAdes and generate assembly summary statistics
source activate spades
echo ""
echo "###################"
echo "## rnaSPAdes"
echo "###################"
echo ""
echo -n "spades version: "
spades.py -v
echo ""
spades.py \
--rna \
--ss rf \
--pe-1 1 ${moduledir}/readprocessing/trim_galore/unfixrm_NA_R1.cor_val_1.fq.gz \
--pe-2 1 ${moduledir}/readprocessing/trim_galore/unfixrm_NA_R2.cor_val_2.fq.gz \
-o NA_rnaspades
source deactivate
source activate trinity
TrinityStats.pl ${moduledir}/NA_rnaspades/transcripts.fasta > NA_assembly.txt
source deactivate
