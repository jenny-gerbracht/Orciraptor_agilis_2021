#!/bin/bash

seqkit \
seq \
-m 200 \
./orciraptor_rnaspades/transcripts.fasta > orciraptor_200.fasta
