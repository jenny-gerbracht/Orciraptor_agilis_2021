#!/bin/bash

###################################################################
#Script Name	:transdecoder_count.sh		                                                                                              
#Description	:Parses transdecoder output and summarises completeness and number of ORFs                                                                               
#Args: pep output from TransDecoder
#Usage: transdecoder_count.sh <TrandDecoder.pep> > output.txt
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

while read line;
do
	echo "Total_ORFs" $((`grep ">" -c`));
done < $1

while read line;
do
	echo "Complete" $((`grep "complete" -c`));
done < $1

while read line;
do
	echo "Internal" $((`grep "internal" -c`));
done < $1

while read line;
do
	echo "5prime_partial" $((`grep "5prime_partial" -c`));
done < $1

while read line;
do
	echo "3prime_partial" $((`grep "3prime_partial" -c`));
done < $1

while read line;
do
	echo "Forward" $((`grep "(+)" -c`));
done < $1

while read line;
do
	echo "Reverse" $((`grep "(-)" -c`));
done < $1
