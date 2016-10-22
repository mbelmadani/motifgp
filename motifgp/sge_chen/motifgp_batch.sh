#!/bin/bash
# Executes a single timed run for MotifGP

#$ -S /bin/bash
#$ -cwd
#$ -o OUT/output
#$ -e OUT/error
#$ -q abaqus.q 
#$ -l qname=abaqus.q

#DATABASES="/home/hpc3019/meme/db/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme,/home/hpc3019/meme/db/motif_databases/CIS-BP/Homo_sapiens.meme,/home/hpc3019/meme/db/motif_databases/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme"

DATABASES="/home/hpc3019/meme/db/motif_databases/JASPAR/JASPAR_CORE_2014_vertebrates.meme,/home/hpc3019/meme/db/motif_databases/MOUSE/uniprobe_mouse.meme"

echo "INPUTNAME=$1 MOO=$2 OUTPUT=$3 FASTA=$4 NGEN=$5 SEED=$6 FITNESS=$7 ARGS=$8"
(time python motifgp.py -p $4 -t $1 --moo $2 -o $3  -n $5 --random-seed $6 -f $7 $8 ) &> $3$1$2"."$6".time" && ./scripts/merge_nef_and_tomtom.sh $3/$1/$2/$6.nef "/home/hpc3019/meme/bin/tomtom" $DATABASES