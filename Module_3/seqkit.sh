#!/bin/bash

seqkit \
grep \
-n -f contaminants_orciraptor.txt \
orciraptor_200.fasta \
-o orciraptor_200_filtered.fasta \
-v
