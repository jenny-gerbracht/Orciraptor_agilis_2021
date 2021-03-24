#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd
#$ -pe threaded 10

####################################
#
# Paths to reads and experiment file
#
####################################

# Paths to the experiment file and to the project folder
mydir="/path/to/wd"
moduledir="${mydir}/Module_5"
experiment_file="${moduledir}/experiment.txt"

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

if [ ! -d "${moduledir}/BAM/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/BAM/"
  mkdir ${moduledir}/BAM/
fi

# Map processed reads to filtered Orciraptor transcriptome
# Create bowtie2 index for Orciraptor transcriptome

source activate bowtie2
bowtie2-build ${mydir}/Module_3/orciraptor_200_filtered.fasta ${moduledir}/BAM/orciraptor_200_filtered.fasta

echo ""
echo "###################"
echo "## mapping bowtie2"
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
  -x ${moduledir}/BAM/orciraptor_200_filtered.fasta \
  -1 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_${i}.1.fq.gz \
  -2 ${mydir}/Module_2/readprocessing/blacklist_NA/blacklist_NA_paired_unaligned_${i}.2.fq.gz \
  --met-file bowtie2_alignment-metrics.txt | samtools view -bS - > ${i}.BAM ) 3>&1 1>&2 2>&3 | tee stderr.log 
  
done
