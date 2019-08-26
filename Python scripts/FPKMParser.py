#!/usr/bin/env python3

"""
Author: Erik IJland
Studentnr.: 960504980110

A parser to obtain only relevant genes from a cuffdiff output file.
"""
#####   DO NOT USE  #####
exit()
from sys import argv


InFileName = argv[1]
OutFileName = argv[2]
lines = open(InFileName, "r")
outputfile = open(OutFileName, "w")

for line in lines:
    if line.startswith("tracking_id"):
        outputfile.write(line + "\n")
    else:
        if "\t0\t0\t0" in line:
            continue
        else:
            outputfile.write(line + "\n")
lines.close()
outputfile.close()
