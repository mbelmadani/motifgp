#!/usr/bin/python

import csv
import sys
import os.path

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-p", "--path", dest="path",
                  help="Path where the tomtom output subdirectories are.", default=None)

(options, args) = parser.parse_args()

if options.path == None:
    sys.err.write("Error, invalid path " + str(options.path) )
    sys.exit(-1)

counter = -1
digest = [] 
header = ""
while True:
    counter+=1
    current_path = options.path + str(counter) + "/tomtom.txt"
    if not os.path.isfile(current_path):
        print "Abort", current_path, "is not a file"
        break #Exit
    header = ["FILE"]
    with open(current_path, 'rb') as tomtom_file:
        reader = csv.reader( tomtom_file, delimiter='\t' )
        for row in reader:
            if row[0][0] == '#' and len(digest) > 0:
                header = header + row
            elif row[0][0] != '#':
                idx = [str(counter)]
                row =  idx + row # Add file number
                digest.append(row)

digest = sorted(digest, key=lambda tup: float(tup[4]))
digest = [header] + digest
for d in digest:
    print "\t".join(d)
