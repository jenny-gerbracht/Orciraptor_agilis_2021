#!/bin/bash

###################################################################
#Script Name	:blastn.sh		                                                                                              
#Description	:Runs blastp search of Orciraptor vs. Mougeotia proteome,
#              returns one hit per ORF
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

source ../config.txt
moduledir="${mydir}/Module_3"

makeblastdb \
-in ${mydir}/Module_1/transdecoder/NA_transcripts.fasta.transdecoder.pep \
-dbtype prot

blastp \
-db ${mydir}/Module_1/transdecoder/NA_transcripts.fasta.transdecoder.pep \
-query ${moduledir}/transdecoder/orciraptor_200_filtered.fasta.transdecoder.pep \
-outfmt '6 qseqid sseqid pident length qcovhsp evalue bitscore' \
-max_target_seqs 1 \
-max_hsps 1 \
-evalue 1e-25 \
-num_threads 4 \
-out blastp_NA.hits
