#!/bin/bash

###################################################################
#Script Name	:seqkit.sh		                                                                                              
#Description	:Removes contigs from fasta that were identified by the blastn search as contamination                                                                              
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

seqkit \
grep \
-n -f contaminants_orciraptor.txt \
orciraptor_200.fasta \
-o orciraptor_200_filtered.fasta \
-v
