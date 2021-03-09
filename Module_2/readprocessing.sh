#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 10

# Performs the following steps: a) k-mer based error correction using Rcorrector, b) discarding unfixable read pairs, c) quality and adapter trimming with Trim Galore!,
# d) read mapping and sorting to rRNA sequences of Orciraptor agilis and Mougeotia sp., e) read mapping and sorting to transcriptome of Mougeotia sp.

####################################
#
# Paths to reads and experiment file
#
####################################

# Paths to the experiment file and to the project folder
mydir="/path/to/wd"
moduledir="${mydir}/Module_2"
experiment_file="${moduledir}/experiment.txt"

####################################
#
# Paths to scripts
#
####################################

Rcorrector_dir="/scratch2/software/anaconda/envs/rcorrector/bin/"
Discard_dir="/scratch4/shess/Jenny/Scripts/"

####################################
#
# Setting up samples array
#
####################################

# Declare samples array
declare -a samples

# Load experiment file, parse for samples and save them into the array
let p=0
while read -r f1 f2; do
  samples[p]="${f1}"
  ((++p))
done < $experiment_file

# Declare condition arrays
declare -a cond

# Load experiment file, parse for conditions and save unique conditions into array
let i=1
cond[0]="control"
while read -r f1 f2; do
  if [[ " ${cond[*]} " == *"$f2"* ]];then
  continue
else
  cond[i]="${f2}"
  ((++i))
fi
done < $experiment_file

# Declare individual condition arrays
arr_length="$((${#cond[@]}-1))"
for i in $( eval echo {0..${arr_length}} );do
  declare -a cond${i}
done

# Load experiment file again, parse for conditions and save filenames into condition-specific arrays.
while read -r f1 f2; do
  for i in $( eval echo {0..${arr_length}} );do
  if [[ "$f2" == "${cond[i]}" ]];then
    eval cond${i}[cond${i}count]="${f1}"
    ((++cond${i}count))
  fi
done
done < $experiment_file

# State the conditions and samples for this analysis
echo "############################"
echo "## Conditions and samples ##"
echo "############################"
echo ""

arr_length="$((${#cond[@]}-1))"
for i in $( eval echo {0..${arr_length}} );do
  echo -e "cond${i} \t ${cond[i]} \t $(eval echo \${cond$i[*]})"
done

####################################
#
# Make folders
#
####################################

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

if [ ! -d "${moduledir}/rRNA/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/rRNA/"
  mkdir ${moduledir}/rRNA/
fi

if [ ! -d "${moduledir}/readprocessing/blacklist_rRNA/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/blacklist_rRNA/"
  mkdir ${moduledir}/readprocessing/blacklist_rRNA/
fi

if [ ! -d "${moduledir}/NA_ref/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/NA_ref/"
  mkdir ${moduledir}/NA_ref/
fi

if [ ! -d "${moduledir}/readprocessing/blacklist_NA/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/readprocessing/blacklist_NA/"
  mkdir ${moduledir}/readprocessing/blacklist_NA/
fi

# kmer-based error correction with Rcorrector
echo ""
echo "###################"
echo "## Rcorrector"
echo "###################"
echo ""

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Start rcorrector of sample ${i}"
  echo ""
  
  perl ${Rcorrector_dir}/run_rcorrector.pl \
  -p ${moduledir}/readprocessing/${i}_1.fq.gz ${moduledir}/readprocessing/${i}_2.fq.gz \
  -od ${moduledir}/readprocessing/corrected \
  -t 10
  
done

# Remove unfixable read pairs, compress output
echo ""
echo "###################"
echo "## Discard"
echo "###################"
echo ""

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Start discarding unfixable read pairs of sample ${i}"
  echo ""
  
  python ${Discard_dir}/FilterUncorrectabledPEfastq.py \
  -1 ${moduledir}/readprocessing/corrected/${i}_1.cor.fq.gz \
  -2 ${moduledir}/readprocessing/corrected/${i}_2.cor.fq.gz \
  -s ${i}
  pigz ${moduledir}/readprocessing/corrected/unfixrm_${i}_1.cor.fq
  pigz ${moduledir}/readprocessing/corrected/unfixrm_${i}_2.cor.fq
  
done

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

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Start trimming of sample ${i}"
  echo ""
  
  trim_galore \
  --fastqc \
  -o ${moduledir}/readprocessing/trim_galore \
  -j 4 \
  --paired \
  --length 49 \
  ${moduledir}/readprocessing/corrected/unfixrm_${i}_1.cor.fq.gz \
  ${moduledir}/readprocessing/corrected/unfixrm_${i}_2.cor.fq.gz
  
done
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

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Start fastqc analysis of sample ${i}"
  echo ""
  
  fastqc \
  -o ${moduledir}/readprocessing/fastqc/raw/ \
  ${moduledir}/readprocessing/${i}_1.fq.gz \
  ${moduledir}/readprocessing/${i}_2.fq.gz
  mv ${moduledir}/readprocessing/trim_galore/unfixrm_${i}_R1.cor_val_1_fastqc* ${moduledir}/readprocessing/fastqc/processed/
  mv ${moduledir}/readprocessing/trim_galore/unfixrm_${i}_R2.cor_val_2_fastqc* ${moduledir}/readprocessing/fastqc/processed/
  
done
source deactivate

# Map reads to Mougeotia and Orciraptor known rRNA sequences, sort in aligned and unaligned
# Create bowtie2 index for fasta file containing Mougeotia and Orciraptor rRNA sequences
source activate bowtie2
bowtie2-build ${moduledir}/ssu_T.fasta ${moduledir}/rRNA/ssu_T.fasta

echo ""
echo "###################"
echo "## rRNA mapping bowtie2"
echo "###################"
echo ""
echo -n "bowtie2 version: "
fastqc -v
echo ""

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Start rRNA mapping of sample ${i}"
  echo ""
  
  bowtie2 \
  -p 10 \
  --very-sensitive \
  --score-min C,0,0 \
  --phred33 \
  --al-conc-gz ${moduledir}/readprocessing/blacklist_rRNA/blacklist_rRNA_paired_aligned_${i}.fq.gz \
  --un-conc-gz ${moduledir}/readprocessing/blacklist_rRNA/blacklist_rRNA_paired_unaligned_${i}.fq.gz \
  -x ${moduledir}/rRNA/ssu_T.fasta \
  -1 ${moduledir}/readprocessing/trim_galore/unfixrm_${i}_1.cor_val_1.fq.gz \
  -2 ${moduledir}/readprocessing/trim_galore/unfixrm_${i}_2.cor_val_2.fq.gz
  
done

# Map reads to Mougeotia transcriptome created in Module_1
# Create bowtie2 index for Mougeotia transcriptome

bowtie2-build ${mydir}/Module_1/NA_rnaspades/transcripts.fasta ${moduledir}/NA_ref/NA_ref.fasta


  
