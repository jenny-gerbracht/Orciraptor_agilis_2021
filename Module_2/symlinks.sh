#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd

# Give absolute paths to the experiment file and to the project folder
srvdir="/scratch2/shess/TRANSCRIPTOMES/ORCIRAPTOR_RNASeq_FEB2017/READS/"
mydir="/scratch4/shess/Jenny/201109_ORC_readprocessing"
experiment_file="${mydir}/experiment.txt"

if [ ! -d "${mydir}/rcorrector/" ]; then		
	    	mkdir ${mydir}/rcorrector/
fi

	

# Declare samples array
declare -a samples

# Load experiment file, parse for samples and save them into the array
let p=0
while read -r f1 f2; do
	samples[p]="${f1}"
	((++p))	
done < $experiment_file

# Iterate over FASTQ files

for i in "${samples[@]}"; do		
			echo "Concatenate ${i}"
			if [ -e "$mydir/rcorrector/${i}_1.fq.gz" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
    		echo "$mytime FASTQ file link already created"
   	else
		ln -s $srvdir/*${i}_R1.fastq.gz $mydir/rcorrector/${i}_1.fq.gz
	fi			
			if [ -e "$mydir/rcorrector/${i}_2.fq.gz" ]; then
		mytime=$(date "+%Y-%m-%d %H:%M:%S")
    		echo "$mytime FASTQ file link already created"
   	else
		ln -s $srvdir/*${i}_R2.fastq.gz $mydir/rcorrector/${i}_2.fq.gz
	fi


done	
