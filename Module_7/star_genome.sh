#!/bin/bash

# Define paths to working directory locations
mydir="/path/to/wd"
moduledir="${mydir}/Module_7"

# Define log file and redirect stdout and stderr to this file
if [ ! -d "${moduledir}/Logs/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "${mytime} Make directory ${moduledir}/Logs/"
  mkdir ${moduledir}/Logs/
fi
log_file="${moduledir}/Logs/log_genome_$date"
exec &> >(tee -a "$log_file")

# Generate SuperTranscripts with Lace
echo ""
echo "###################"
echo "## Generate STAR genome"
echo "###################"
echo ""

STAR \
--runThreadN 10 \
--genomeSAindexNbases 11 \
--runMode genomeGenerate \
--genomeDir ${module_dir}/genome \
--genomeFastaFiles  ${module_dir}/SuperDuper.fasta \
--sjdbGTFfile ${module_dir}/SuperDuperTrans.gff \
--sjdbGTFfeatureExon exon \
--sjdbGTFtagExonParentTranscript trans_id \
--sjdbGTFtagExonParentGene gene_id
