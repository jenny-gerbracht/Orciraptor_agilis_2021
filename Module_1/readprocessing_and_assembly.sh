#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe -threaded 10

# Define paths to read and working directory locations
readdir="/path/to/reads"
mydir="/path/to/wd"
moduledir="${mydir}/Module_1"

# Define paths to scripts
Rcorrector_dir="/scratch2/software/anaconda/envs/rcorrector/bin/"
Discard_dir="/scratch4/shess/Jenny/Scripts/"

# Log file
# Get date
date=$(date "+%Y-%m-%d")

# Define log file and redirect stdout and stderr to this file
if [ ! -d "${moduledir}/Logs/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/Logs/"
  mkdir ${moduledir}/Logs/
fi
log_file="${moduledir}/Logs/log_$date"
exec &> >(tee -a "$log_file")

# Create necessary folders
if [ ! -d "${moduledir}/readprocessing/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/"
  mkdir${moduledir}/readprocessing/
fi

if [ ! -d "${moduledir}/readprocessing/corrected/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/corrected/"
  mkdir ${moduledir}/readprocessing/corrected/
fi

if [ ! -d "${moduledir}/readprocessing/fastqc/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/fastqc/"
  mkdir ${moduledir}/readprocessing/fastqc/
fi

if [ ! -d "${moduledir}/readprocessing/fastqc/raw/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/fastqc/raw/"
  mkdir ${moduledir}/readprocessing/fastqc/raw/
fi

if [ ! -d "${moduledir}/readprocessing/fastqc/processed/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/fastqc/processed"
  mkdir ${moduledir}/readprocessing/fastqc/processed/
fi

# Make symlinks
ln -s $readdir/HI.3499.002.D704---D504.NA_R1.fastq ${moduledir}/readprocessing/NA_R1.fastq
ln -s $readdir/HI.3499.002.D704---D504.NA_R2.fastq ${moduledir}/readprocessing/NA_R2.fastq

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
python ${Discard_dir}/FilterUncorrectabledPEfastq.py \
-1 ${moduledir}/readprocessing/corrected/NA_R1.cor.fastq \
-2 ${moduledir}/readprocessing/corrected/NA_R2.cor.fastq \
-s NA
pigz unfixrm_NA_R1.cor.fq
pigz unfixrm_NA_R2.cor.fq

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
--pe-1 1 ${moduledir}/readprocessing/trim_galore \
--pe-2 1 ${moduledir}/readprocessing/trim_galore \
-o NA_rnaspades
source deactivate
source activate trinity
perl TrinityStats.pl ${moduledir}/NA_rnaspades/transcripts.fasta NA_assembly.txt
source deactivate
