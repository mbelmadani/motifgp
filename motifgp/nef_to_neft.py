"""
Post-processing script.

Takes in a Network Expression Pareto Front (.nef) file
Creates a .MEME format file for each motifs.
"""

import fileinput
import os
import glob
import paretotopssm
import sys

if len(sys.argv) < 2:
    print "USAGE:", sys.argv, "filename", "motifevaluation (OPTIONAL), databases(comma-delimited)( OPTIONAL ) "
    print "Filename is a path. motifevaluation is an optional string for motif evaluation. Use 'tomtom', for example."
else:
    print "Converting NEF to MEME:", sys.argv[1]

r = paretotopssm.NEF_to_MEME(sys.argv[1])

# Process TOMTOM if 2nd argument is provided
if len(sys.argv) > 2 and sys.argv[2]: # So far
    TOMTOM_BIN = sys.argv[2]
    oldfiles = glob.glob(os.path.join(sys.argv[1], "*tomtom.*"))
    print "Deleting", len(oldfiles), "stale TOMTOM files"
    for oldfile in oldfiles:
        os.remove(oldfile)

    oldfiles = glob.glob(os.path.join(sys.argv[1], "*.neft"))
    print "Deleting", len(oldfiles), "stale NEFT files"
    for oldfile in oldfiles:
        os.remove(oldfile)

    count = 0
    if len(sys.argv) > 3:
        databases = " ".join(sys.argv[3].split(","))
    else:
        databases = " ~/meme/db/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme "

    if len(sys.argv) > 4:
        evalue_threshold = sys.argv[4]
    else:
        evalue_threshold = "0.5"

    for meme in r:
        tomtom_out = str(sys.argv[1])[:-4] #Trim .nef extension
        tomtom_out = os.path.join(tomtom_out, str(count))
        if not os.path.isdir(tomtom_out):
            os.mkdir(tomtom_out)
        print "Executing tomtom on", meme
        # SYSTEM CALL TO TOMTOM
        os.system(
            TOMTOM_BIN + " -no-ssc -oc " +tomtom_out +
            " -verbosity 1 -min-overlap 5 -mi 1 -dist pearson -evalue -thresh "+evalue_threshold+" " + meme + 
            " "+databases+" "
            ) # ~/meme/db/motif_databases/MOUSE/uniprobe_mouse.meme")
        print "Writing TOMTOM output to:", tomtom_out
        count+=1
    
    # Gather .nef files, match with TOMTOM output.
    with open(sys.argv[1], 'r') as f:
        nef_list = []
        for line in f:
            # Nef files
            nef_list.append(line)

        # Iterate nef files
        count = 0
        for meme in r:
            tomtom_out = str(sys.argv[1])[:-4] #Trim .nef extension
            tomtom_out = os.path.join(tomtom_out, str(count))
            tomtom_out = os.path.join(tomtom_out, "tomtom.txt")

            with open(tomtom_out, 'r') as tomtom_file:
                for line in tomtom_file:
                    tomtom_header = line
                    break
            
            alignment_line = None
            if os.path.isfile(tomtom_out):
                with open(tomtom_out, 'r') as tomtom_file:
                    tomtom_file.next() #Skip header
                    for line in tomtom_file:
                        alignment_line = line
                        break
            if alignment_line:
                nef_list[count + 1] = str(nef_list[count + 1]).strip() + ("\t" + alignment_line) 
            else:
                print "Skip;", tomtom_out, "is empty"
            count+=1
    nef_list[0] = str(nef_list[0]).strip() + ("\t" + tomtom_header) 
                
    # Creat a new version of the nef file with Tomtom
    NEFT_CONTENT = "".join(nef_list)
    NEFT_OUTPUT = sys.argv[1] + "t" #Adding 'neft' as the extension
    with open(NEFT_OUTPUT, 'w') as f:
        print "Writing .neft file at", NEFT_OUTPUT
        f.write(NEFT_CONTENT)
        

        
        
