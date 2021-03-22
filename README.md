# Orciraptor_agilis_2021
Read processing and filtering, *de novo* transcriptome assembly, differential gene expression analysis and functional annotation of *Orciraptor agilis* RNA-seq data

## Module 1: Read processing and de novo transcriptome assembly of prey organism *Mougeotia* sp.

1. [MISSING, UPLOAD SRA!] Download of *Mougeotia* sp. RNA-seq data 
2. Run readprocessing_and_assembly.sh, output is *de novo* transcriptome assembly of *Mougeotia* sp.

## Module 2: Read processing of *Orciraptor agilis*

1. [MISSING, UPLOAD SRA!] Download *Orciraptor agilis* RNA-seq data
2. Run symlinks.sh
3. Run readprocessing.sh. Output are quality filtered and adapter trimmed reads.
4. Run mapping.sh. Output are reads that do not map to sequences from rRNA and/or *Mougeotia* sp.

## Module 3: *De novo* transcriptome assembly, decontamination, ORF prediction

1. Run assembly.sh to assemble the transcriptome from processed reads of all libraries. Output is *de novo* transcriptome assembly of *Orciraptor agilis* as a fasta.
2. Filter transcriptome for contigs larger than 200 nt with removesmalls.pl. Usage in folder Module_3:
```
perl removesmalls.pl 200 ${moduledir}/orciraptor_rnaspades/transcripts.fasta > orciraptor_200.fasta
```
3. Run blastn search with this transcriptome (nt database v5 updated on 2021-03-10)
  * Checked contigs with > 95% identity over a length of minimum 100 nt, saved contig identifiers of all bacterial, viral, ribosomal and algal contigs in contaminants.txt
  * Remove these sequences from transcriptome with seqkit.sh
4. ORF prediction with transdecoder.sh (1 strand)

## Module 4: Functional annotation
? entap ?
? Trinotate?
CAZY -> how to parse

## Module 4: Differential gene expression analysis
1) Mapping and counting
2) R script with DESeq2 and Figures

## Module 5: Assembly summary statistics
1) number of Genes, isoforms, ORFs, 
