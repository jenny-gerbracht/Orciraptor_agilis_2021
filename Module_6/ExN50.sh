#!/bin/bash
#$ -S /bin/bash
. /etc/profile
#$ -cwd

source ../config.txt

${Trinity_dir}/util/abundance_estimates_to_matrix.pl \
--est_method salmon \
${mydir}/Module_5/salmon/V1S1.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S2.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S3.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S4.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S5.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S6.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S7.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S8.salmon_quant/quant.sf \
${mydir}/Module_5/salmon/V1S9.salmon_quant/quant.sf \
