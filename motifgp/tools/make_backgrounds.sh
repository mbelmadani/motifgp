#!/bin/bash

FASTA=($(ls -d $1)) # | tr " " "\n")
ls $1 | parallel -j800% -I @ python motifgp.py -p @ -n 0

#for f in $FASTA; do
#    python motifgp.py -p "$fasta" -n 0
#done
