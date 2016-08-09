#!/bin/bash
"""
Take dump files from runtime_tmp and convert them back to fasta
"""


for sequences_ in $(ls runtime_tmp/ | grep -v ".fasta$"); do
    sequences="runtime_tmp/"$sequences_
    c=0
    echo "" > $sequences".fasta"
    while read seq; do
	echo ">"$c"_shuffled" >> $sequences".fasta"
	echo $seq >> $sequences".fasta"
	let c=c+1
    done < $sequences
done 
