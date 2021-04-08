#!/bin/bash

seqkit \
grep \
-n -f contaminants_NA.txt \
./transdecoder/orciraptor_200_filtered.fasta.transdecoder.pep \
-o /transdecoder/orciraptor_200_filtered2.fasta.transdecoder.pep \
-v
