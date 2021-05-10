  #!/bin/bash

source ../config.txt
moduledir="${mydir}/Module_7"

gffread \
${moduledir}/stringtie_merged.gtf \
-g ${moduledir}/SuperDuper.fasta \
-w ${moduledir}/SuperDuper_stringtie.fasta \
-F
