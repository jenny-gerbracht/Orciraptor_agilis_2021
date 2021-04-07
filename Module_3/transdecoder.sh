#!/bin/bash

moduledir="${mydir}/Module_3"

if [ ! -d "${moduledir}/transdecoder/" ]; then
  mytime=$(date "+%Y-%m-%d %H:%M:%S")
  echo "$mytime Make directory ${moduledir}/transdecoder/"
  mkdir ${moduledir}/transdecoder/
fi

cd transdecoder

TransDecoder.LongOrfs -t ../orciraptor_200_filtered.fasta
TransDecoder.Predict -t ../orciraptor_200_filtered.fasta
