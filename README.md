# Orciraptor_agilis_2021
Read processing and filtering, *de novo* transcriptome assembly (rnaSPAdes) of *Orciraptor agilis* RNA-seq data

## Module 1: Read processing and *de novo* transcriptome assembly of prey organism *Mougeotia* sp.

1. Run readprocessing_and_assembly.sh, output is *de novo* transcriptome assembly of *Mougeotia* sp.
2. Predict ORFs to use later for decontamination: transdecoder.sh
3. Filter transcriptome for contigs larger than 200 nt with seqkit_length.sh (only for upload)

## Module 2: Read processing of *Orciraptor agilis*

1. Run symlinks.sh
2. Run readprocessing.sh. Output are quality filtered and adapter trimmed reads.
3. Run mapping.sh. Output are reads that do not map to sequences from rRNA and/or *Mougeotia* sp.

## Module 3: *De novo* transcriptome assembly, decontamination, ORF prediction

1. Run assembly.sh to assemble the transcriptome from processed reads of all libraries. Output is *de novo* transcriptome assembly of *Orciraptor agilis* as a fasta.
2. Filter transcriptome for contigs larger than 200 nt with seqkit_length.sh
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
## Module 5: Generate gene_trans_map
1) Generate gene_trans_map file for Lace: gene_trans_map.R

## Module 6: Assembly summary statistics
1) number of Genes, isoforms, ORFs (status), Ex90 

## Module 7: Supertranscripts
1) Run Lace to generate supertranscript fasta: lace.sh
2) Generate genome index of supertranscriptome for STAR mapping (star_genome.sh), perform STAR mapping of processed reads (star_mapping.sh), index the bam files (index.sh)
3) Run stringtie and merge the gtfs (stringtie.sh)
4) Convert gtf to fasta and predict ORFs (stringtie_fasta.sh and stringtie_transdecoder.sh)
