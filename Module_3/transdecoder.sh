#!/bin/bash

mkdir transdecoder
cd transdecoder

TransDecoder.LongOrfs -t ../orciraptor_200_filtered.fasta
TransDecoder.Predict -t ../orciraptor_200_filtered.fasta
