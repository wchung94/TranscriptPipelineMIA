#!/usr/bin/env python3
"""
Author: Nijveen group
Sript for running tool HISAT2
"""

from sys import argv
import subprocess
import os

def run_hisat2(index_file, input_single_strand, out_hisat2):
    """ Run HISAT2 program on FASTQ file
    
    index_file: string, filename of index file of the organism of the input files
    input_single_strand: string, filename of input FASTA/FASTQ file of single stranded sequence
    out_hisat2: string, filename of output .sam file of hisat2 
    """

    cmd = "hisat2 -p 3 --dta-cufflinks -x {} -U {} -S {}"\
    .format(index_file, input_single_strand, out_hisat2)
    if not os.path.exists (out_hisat2):
        process_hisat2 = subprocess.check_call(cmd, shell=True)
    else:
        print("File already exists")
    
if __name__ == "__main__":

    # run tool HISAT2
    index = argv[1]
    input_file = argv[2]
    out_file = argv[3]
    run_hisat2(index, input_file, out_file)
    

