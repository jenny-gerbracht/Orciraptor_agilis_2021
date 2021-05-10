#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 10

# Performs the following steps: d) read mapping and sorting to rRNA sequences of Orciraptor agilis and Mougeotia sp., 
# e) read mapping and sorting to transcriptome of Mougeotia sp.

source ../config.txt
moduledir="${mydir}/Module_2"

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
let i=0
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

# Map reads to Mougeotia and Orciraptor known rRNA sequences, sort in aligned and unaligned
# Create bowtie2 index for fasta file containing Mougeotia and Orciraptor rRNA sequences from SILVA SSU r138.1 database: sequences from groups "Orciraptor" and "Mougeotia",
# U replaced with T

source activate bowtie2
bowtie2-build ${moduledir}/rRNA/ssu_T.fasta ${moduledir}/rRNA/ssu_T.fasta

echo ""
echo "###################"
echo "## rRNA mapping bowtie2"
echo "###################"
echo ""
echo -n "bowtie2 version: "
bowtie2 --version
echo ""

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Start rRNA mapping of sample ${i}"
  echo ""
  
  (bowtie2 \
  -p 10 \
  --very-sensitive \
  --score-min C,0,0 \
  --phred33 \
  --al-conc-gz ${moduledir}/readprocessing/blacklist_rRNA/blacklist_rRNA_paired_aligned_${i}.fq.gz \
  --un-conc-gz ${moduledir}/readprocessing/blacklist_rRNA/blacklist_rRNA_paired_unaligned_${i}.fq.gz \
  -x ${moduledir}/rRNA/ssu_T.fasta \
  -1 ${moduledir}/readprocessing/trim_galore/unfixrm_${i}_1.cor_val_1.fq.gz \
  -2 ${moduledir}/readprocessing/trim_galore/unfixrm_${i}_2.cor_val_2.fq.gz \
  --met-file bowtie2_alignment-metrics_rRNA.txt | samtools view -bS - > rRNA.BAM ) 3>&1 1>&2 2>&3 | tee rRNA.stderr.log 
  rm rRNA.BAM
  
done

# Map reads to Mougeotia transcriptome created in Module_1
# Create bowtie2 index for Mougeotia transcriptome

bowtie2-build ${mydir}/Module_1/NA_rnaspades/transcripts.fasta ${moduledir}/NA_ref/NA_ref.fasta

echo ""
echo "###################"
echo "## Prey organism mapping bowtie2"
echo "###################"
echo ""
echo -n "bowtie2 version: "
bowtie2 --version
echo ""

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Start mapping of sample ${i}"
  echo ""
  
  (bowtie2 \
  -p 10 \
  --phred33 \
  --al-conc-gz ${moduledir}/readprocessing/blacklist_NA/blacklist_NA_paired_aligned_${i}.fq.gz \
  --un-conc-gz ${moduledir}/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_${i}.fq.gz \
  -x ${moduledir}/NA_ref/NA_ref.fasta \
  -1 ${moduledir}/readprocessing/blacklist_rRNA/blacklist_rRNA_paired_unaligned_${i}.fq.1.gz \
  -2 ${moduledir}/readprocessing/blacklist_rRNA/blacklist_rRNA_paired_unaligned_${i}.fq.2.gz \
  --met-file bowtie2_alignment-metrics_mougeotia.txt | samtools view -bS - > mougeotia.BAM ) 3>&1 1>&2 2>&3 | tee mougeotia.stderr.log 
  rm mougeotia.BAM

done

# Rename output reads that are going to be used by the assembler

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Rename sample ${i}"
  echo ""
  
  mv ${moduledir}/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_${i}.fq.1.gz ${moduledir}/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_${i}.1.fq.gz
  mv ${moduledir}/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_${i}.fq.2.gz ${moduledir}/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_${i}.2.fq.gz
 
done
