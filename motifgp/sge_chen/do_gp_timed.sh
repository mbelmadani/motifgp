#!/bin/bash

TYPE="MOEADTEST"
#TYPE="AUTO"
VERSION="0.1.2"
POPSIZE="1000"
MATCHER="python"
FITNESS="FD"

indir="./CHEN/"
outdir="./CHEN_OUT_"$TYPE"_"$VERSION"_"$FITNESS"_"$POPSIZE"_"$MATCHER"/"

mkdir -p $outdir

STARTGEN=0
ENDGEN=100000
MAXSEED=109

TIMEFILE="sge_chen/chen_timelimits"

FASTAS=$(ls $indir/*.fasta)

#for MOO in "NSGAR" "SPEA2" "MOEAD ; do
for MOO in "NSGAR" "SPEA2" ; do
    for FASTA in $FASTAS; do
	i=100 # SEED
	target=$(basename $FASTA | sed s/.pos.fasta//g )
	timelimit=$(grep $target $TIMEFILE | cut -f2 -d" ")
	ARGS=" --timelimit=$timelimit --matcher=$MATCHER --popsize=$POPSIZE "
	while [ $i -le $MAXSEED ]
	do
	    FASTABASENAME=$(basename $FASTA)
	    mkdir -p $outdir	 
	    echo "Queuing motifgp" $FASTABASENAME $MOO $outdir $FASTA $ENDGEN $i $FITNESS "$ARGS"
	    qsub sge_chen/motifgp_batch.sh $FASTABASENAME $MOO $outdir $FASTA $ENDGEN $i $FITNESS "$ARGS"
	    i=$(( $i + 1 ))	    
	    #sleep 1
	done
    done
    break
done