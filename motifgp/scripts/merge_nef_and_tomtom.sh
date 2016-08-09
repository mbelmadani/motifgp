#!/bin/bash

# Takes an input .nef files and aligns the motifs with TOMTOM.
NEF_FILE=$1
ENGINE=$2
DATABASES=$3
echo "Processing $NEF_FILE with $ENGINE on $DATABASES"

# Create MEME files at nef subdirectory
python nef_to_neft.py $NEF_FILE $ENGINE $DATABASES

