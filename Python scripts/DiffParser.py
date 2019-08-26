#!/usr/bin/env python3

"""
Author: Erik IJland
Studentnr.: 960504980110

A parser to obtain only relevant genes from a cuffdiff output file.
"""

from sys import argv


InFileName = argv[1]
OutFileName = "DiffParserOut.txt" #HARDCODED RIGHT NOW!!!
lines = open(InFileName, "r")
outputfile = open(OutFileName, "w")
dataset = ""
for line in lines:
    if line.startswith("test_id"):
        outputfile.write(line + "\n")
    elif "yes" not in line:
        continue
    else:
        outputfile.write(line + "\n")
lines.close()
outputfile.close()
