  #!/bin/bash
  
###################################################################
#Script Name	:stringtie_fasta.sh		                                                                                              
#Description	:Obtain sequences as fasta from stringtie gtf file                                                       
#Args:                                                                                           
#Author	:Jennifer Gerbracht                                               
#Email		:jennifer.gerbracht@gmx.de                                           
###################################################################

source ../config.txt
moduledir="${mydir}/Module_7"

gffread \
${moduledir}/stringtie_merged.gtf \
-g ${moduledir}/SuperDuper.fasta \
-w ${moduledir}/SuperDuper_stringtie.fasta \
-F
