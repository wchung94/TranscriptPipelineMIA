#!/usr/bin/env python3
"""
Author: Nijveen group
Sript for running tool samtools
"""

from sys import argv
import subprocess
import os

def run_samview(in_sam, unsorted_bam):
    """ Run samtools view program on .sam file
    
    in_sam: string, filename of input .sam file
    unsorted_bam: string, filename of output unsorted .bam file
    """
    cmd = " samtools view -bS {} > {}".format(in_sam, unsorted_bam)
    if not os.path.exists (unsorted_bam):
        process_sam = subprocess.check_call(cmd, shell=True)
    else:
        print("File already exists")


def run_samsort(unsorted_bam, sorted_bam):
    """ Run samtools sort program on unsorted .bam file
    
    unsorted_bam: string, filename of output unsorted .bam file
    sorted_bam: string, filename of output sorted .bam file
    """
    cmd = " samtools sort {} {}".format(unsorted_bam, sorted_bam)
    if not os.path.exists (sorted_bam):
        process_sam = subprocess.check_call(cmd, shell=True)
    else:
        print("File already exists")


if __name__ == "__main__":    
        
    # convert to .bam
    in_sam = argv[1]
    unsort_bam = argv[2]
    sort_bam = argv[3]
    run_samview(in_sam, unsort_bam)
    run_samsort(unsort_bam, sort_bam)
