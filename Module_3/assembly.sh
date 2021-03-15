#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 10

# Define paths to working directory locations
mydir="/path/to/wd"
moduledir="${mydir}/Module_3"

# Log file
# Get date
date=$(date "+%Y-%m-%d")

# Define log file and redirect stdout and stderr to this file
if [ ! -d "${moduledir}/Logs/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/Logs/"
  mkdir ${moduledir}/Logs/
fi
log_file="${moduledir}/Logs/log_$date"
exec &> >(tee -a "$log_file")

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
--pe-1 1 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S1.1.fq.gz \
--pe-1 2 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S2.1.fq.gz \
--pe-1 3 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S3.1.fq.gz \
--pe-1 4 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S4.1.fq.gz \
--pe-1 5 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S5.1.fq.gz \
--pe-1 6 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S6.1.fq.gz \
--pe-1 7 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S7.1.fq.gz \
--pe-1 8 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S8.1.fq.gz \
--pe-1 9 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S9.1.fq.gz \
--pe-2 1 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S1.2.fq.gz \
--pe-2 2 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S2.2.fq.gz \
--pe-2 3 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S3.2.fq.gz \
--pe-2 4 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S4.2.fq.gz \
--pe-2 5 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S5.2.fq.gz \
--pe-2 6 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S6.2.fq.gz \
--pe-2 7 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S7.2.fq.gz \
--pe-2 8 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S8.2.fq.gz \
--pe-2 9 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_V1S9.2.fq.gz \
-o orciraptor_rnaspades
source deactivate

source activate trinity
TrinityStats.pl ${moduledir}/orciraptor_rnaspades/transcripts.fasta > orciraptor_assembly.txt
source deactivate
