#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -pe threaded 1

source ../../config.txt
moduledir="${mydir}/Module_5"

${Trinity_dir}/util/misc/fasta_seq_length.pl \
${mydir}/Module_3/orciraptor_200_filtered.fasta \
> transcripts.fasta.seq_lens

${Trinity_dir}/util/misc/TPM_weighted_gene_length.py \
--gene_trans_map ${moduledir}/gene_trans_map_trinity \
--trans_lengths transcripts.fasta.seq_lens \
--TPM_matrix ${moduledir}/salmon_isoform_TPM_matrix \
> gene_lengths.txt

${Trinity_dir}/Analysis/DifferentialExpression/run_GOseq.pl \
--factor_labeling /factor_labeling.txt \
--GO_assignments go_annotations.txt \
--lengths gene_lengths.txt \
--background background.txt
