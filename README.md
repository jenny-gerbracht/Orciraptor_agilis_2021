# Orciraptor_agilis_2021
Read processing and filtering, *de novo* transcriptome assembly (rnaSPAdes) of *Orciraptor agilis* RNA-seq data

## Module 1 + 2: See repository "Orciraptor_agilis_2021_Trinity"

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
2) Mapping the processed reads back to the newly generated transcriptome with bowtie2 and counting with salmon in alignment-mode (bowtie2.sh).

## Module 6: Assembly summary statistics
1) Number of number + length statistics of contigs was calculated with TrinityStats.pl script from Trinity toolkit
2) Number, completeness and orientation of ORFs is summarised with transdecoder_count.sh
3) ExN50 statistic is calculated with ExN50.sh 

## Module 7: Supertranscripts
1) Run Lace to generate supertranscript fasta: lace.sh
2) Generate genome index of supertranscriptome for STAR mapping (star_genome.sh), perform STAR mapping of processed reads (star_mapping.sh), index the bam files (index.sh)
3) Run stringtie and merge the gtfs (stringtie.sh)
4) Convert gtf to fasta and predict ORFs (stringtie_fasta.sh and stringtie_transdecoder.sh)
