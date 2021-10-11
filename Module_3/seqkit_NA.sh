#!/bin/bash

###################################################################
#Script Name	:seqkit_NA.sh		                                                                                              
#Description	:Removes ORFs that were identified in blastp search                                                                              
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

seqkit \
grep \
-n -f contaminants_NA.txt \
./transdecoder/orciraptor_200_filtered.fasta.transdecoder.pep \
-o ./transdecoder/orciraptor_200_filtered2.fasta.transdecoder.pep \
-v
