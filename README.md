# Orciraptor_agilis_2021
Read processing and filtering, *de novo* transcriptome assembly, differential gene expression analysis and functional annotation of *Orciraptor agilis* RNA-seq data

## Module 1: Read processing and *de novo* transcriptome assembly of prey organism *Mougeotia* sp.

1. [MISSING, UPLOAD SRA!] Download of *Mougeotia* sp. RNA-seq data 
2. Run readprocessing_and_assembly.sh, output is *de novo* transcriptome assembly of *Mougeotia* sp.
3. Predict ORFs to use later for decontamination: transdecoder.sh

## Module 2: Read processing of *Orciraptor agilis*

1. [MISSING, UPLOAD SRA!] Download *Orciraptor agilis* RNA-seq data
2. Run symlinks.sh
3. Run readprocessing.sh. Output are quality filtered and adapter trimmed reads.
4. Run mapping.sh. Output are reads that do not map to sequences from rRNA and/or *Mougeotia* sp.

## Module 3: *De novo* transcriptome assembly, decontamination, ORF prediction

1. Run assembly.sh to assemble the transcriptome from processed reads of all libraries. Output is *de novo* transcriptome assembly of *Orciraptor agilis* as a fasta.
2. Filter transcriptome for contigs larger than 200 nt with removesmalls.pl. Usage in folder Module_3: [REPLACE WITH SOMETHING OWN]
```
perl removesmalls.pl 200 ${moduledir}/orciraptor_rnaspades/transcripts.fasta > orciraptor_200.fasta
```
3. Run blastn search with this transcriptome (nt database v5 updated on 2021-03-10): blastn.sh
  * Checked contigs with > 95% identity over a length of minimum 100 nt, saved contig identifiers of all bacterial, viral, ribosomal and algal contigs in contaminants.txt
  * Remove these sequences from transcriptome with seqkit.sh
4. ORF prediction with transdecoder.sh
5. Run blastp search of *Orciraptor* ORFs against *Mougeotia* predicted ORFs: blastp.sh  
  * Remove ORFs with a > 95% identity over a length of 150 aa from *Orciraptor* predicted proteome: seqkit_NA.sh
6. Rename ORFs to pattern "gx_iy.pz" (gene, isoform, peptide) with rename_transdecoder.py. Usage in folder Module_3/transdecoder: Output is "orciraptor_transdecoder.pep_renamed.fasta".
```
python rename_transdecoder.py orciraptor_200_filtered2.fasta.transdecoder.pep
```

## Module 4: Functional annotation
1. Run eggnog-mapper in diamond and hmm mode: eggnog.sh. Parse so that hmm annotation is used if there is no diamond annotation: eggnog_parse.R. 
2. Run InterProScan using interproscan.sh
3. Run a Diamond blastp search vs nr database (nr database downloaded in fasta format from https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/ on 2021-03-21 and formated as diamond db): diamond.sh
4. Annotation of carbohydrate-active enzymes (CAZymes) with with dbcan2 in HMM mode (database dbCAN-HMMdb-V9): cazy.sh
5. Annotation of peptidases with MEROPS database (database merops_scan.lib updated on 2019-05-19 from https://ftp.ebi.ac.uk/pub/databases/merops/current_release/): merops.sh

## Module 5: Differential gene expression analysis
1) Mapping the processed reads back to the newly generated transcriptome with bowtie2 and counting with salmon in alignment-mode (bowtie2.sh).
2) Perform 
3) Output table with parsed functional annotation and expression info
4) Generate Figures 

## Module 6: Assembly summary statistics
1) number of Genes, isoforms, ORFs (status), Ex90 

## Module 7: Supertranscripts
1) Run Lace to generate supertranscript fasta: lace.sh
2) Generate genome index of supertranscriptome for STAR mapping, perform STAR mapping of BAM files
3) Run stringtie and merge the gtfs
4) Convert gtf to fasta and predict ORFs
