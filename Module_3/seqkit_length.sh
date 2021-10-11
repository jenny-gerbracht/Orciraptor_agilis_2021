#!/bin/bash

###################################################################
#Script Name	:seqkit_length.sh		                                                                                              
#Description	:Removes contigs from fasta that are smaller than 200 nt                                                                      
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

seqkit \
seq \
-m 200 \
./orciraptor_rnaspades/transcripts.fasta > orciraptor_200.fasta
