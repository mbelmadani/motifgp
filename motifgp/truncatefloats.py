"""
Read a file line by line. Whatever is a float, trunc the value to the 3rd decimal
"""

import sys

a_file = sys.argv[1]
with open(a_file) as f:
    for line in f:
        splitted = line.split(" ")
        def do_trunc(text):
            try:
                float(text)
                return str(round(float(text), 3))
            except ValueError:
                return text

        splitted = [do_trunc(x) for x in splitted]
        print " ".join(splitted)

