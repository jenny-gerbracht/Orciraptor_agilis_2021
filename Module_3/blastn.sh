#!/bin/bash

blastn -db path/to/nt \
       -query orciraptor_200.fasta \
       -outfmt '6 qseqid sseqid pident length evalue bitscore sgi sacc staxids sskingdoms sscinames scomnames stitle' \
       -max_target_seqs 1 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 4 \
       -out Module_3_tax_1hit.hits
