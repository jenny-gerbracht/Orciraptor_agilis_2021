# Orciraptor_agilis_2021
Read processing and filtering, de novo transcriptome assembly, differential gene expression analysis and functional annotation of Orciraptor agilis RNA-seq data

## Module 1: Read processing and de novo transcriptome assembly of prey organism Mougeotia sp.

1) [MISSING, UPLOAD SRA!] Download of Mougeotia sp. RNA-seq data 
2) Run readprocessing_and_assembly.sh, output is de novo transcriptome assembly of Mougeotia sp.

## Module 2: Read processing of Orciraptor agilis

1) [MISSING, UPLOAD SRA!] Download Orciraptor agilis RNA-seq data
2) Run symlinks.sh
3) Run readprocessing.sh. Output are quality filtered and adapter trimmed reads that do not map to sequences from rRNA and/or Mougeotia sp.

## Module 3: De novo transcriptome assembly, decontamination, summary statistics, ORF prediction

1) Run assembly.sh 
2) Run blastn search with transcriptome
    -> Checked contigs with > 90% identity over a length of 100 nt, removed all bacterial, viral, ribosomal and algal contigs
4) ORF prediction with transdecoder.sh (1 strand)

## Module 4: Functional annotation
? entap ?
? Trinotate?
CAZY -> how to parse

## Module 4: Differential gene expression analysis
1) Mapping and counting
2) R script with DESeq2 and Figures

## Module 5: Assembly summary statistics
1) # Genes, isoforms, ORFs, 
