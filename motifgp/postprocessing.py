#!/usr/bin/python
"""
Any kind of functions that can be considered post-processing on existing data
"""

PADLENGTH=8
ALPHABET=["A","C","G","T"]

def backpad(_sequences):
    """
    Pads a set of (background) sequences with length 8 strings of continuous nucleotides.
    This serves the purpose of filtering continues nucleotide sequences from occuring in predictions.
    """

    PADS = [x*PADLENGTH for x in ALPHABET]
    sequences = []
    i=0
    for seq in _sequences:
        sequences.append(PADS[i%len(PADS)]+seq)
        i+=1
    return sequences
