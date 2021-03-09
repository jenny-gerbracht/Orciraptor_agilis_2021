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
srvdir="path/to/reads"
mydir="wd"
experiment_file="${mydir}/experiment.txt"

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

# Load experiment file, parse for conditions and save unique conditions into array. Note: Condition1 is always "control"
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

# Start rcorrector

for i in "${samples[@]}"; do
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo ""
  echo "$mytime Start rcorrector of sample ${i}"
  echo ""
  
  perl /scratch2/software/anaconda/envs/rcorrector/bin/run_rcorrector.pl \
  -p $mydir/rcorrector/${i}_1.fq.gz $mydir/rcorrector/${i}_2.fq.gz \
  -od $mydir/rcorrector/corrected \
  -t 10
  
done
