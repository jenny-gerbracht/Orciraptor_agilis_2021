  
#!/bin/bash

# Define paths to working directory locations
mydir="/path/to/wd"
moduledir="${mydir}/Module_7"
experiment_file="${moduledir}/experiment.txt"

gffread \
${moduledir}/stringtie_merged.gtf \
-g ${moduledir}/SuperDuper.fasta \
-w ${moduledir}/SuperDuper_stringtie.fasta \
-F
