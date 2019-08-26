#! /usr/bin/env python3
'''
Run cuffdiff 2 using gtf and sam files
'''

# GLOBAL FASTQ DATA LIST

from sys import argv
import os
import subprocess

def select_bam_files():
    sorted_bam_list = []
    for file in os.listdir('bam_files'):
        if file.endswith('.bam') and not file.endswith('_unsorted.bam'):
          sorted_bam_list.append('bam_files/' + file)
    sorted_bam_list = sorted(sorted_bam_list)
    return sorted_bam_list


def run_cuffdiff(gtf_file, sam_file_list):
    """Executes cuffdiff with an input_file ad save output in output_file
    sam_file, gtf_file: string, location and name of input file
    return: None
    """

    # check if the output_file already exists
    if not os.path.exists(gtf_file):
        print('gtf file does not exist')
        return
    for file in sam_file_list:
        if not os.path.exists(file):
            print('sam file not already exist')
            return

    cmd_cuffdiff = "cuffdiff -o {} -p 4 --library-type fr-firststrand -u -no-update-check {} {},{} {},{},{} {},{},{}"\
    .format('cuffdiff_output_directory_CR', gtf_file, sam_file_list[0], sam_file_list[7], sam_file_list[1], sam_file_list[2], sam_file_list[5], sam_file_list[3], sam_file_list[4], sam_file_list[6])
    e = subprocess.check_output(cmd_cuffdiff, shell=True) #e stores the output
    print('process start')
    subprocess.check_call(cmd_cuffdiff, shell=True)
    return

if __name__ == '__main__':
    #Select sorted bam files from directory
    sorted_bam_list = select_bam_files()
    print(sorted_bam_list)
    #run cuffdiff on selected bam files
    gtf_file = 'stringtie_gtfs/stringtie_merged.gtf'
    run_cuffdiff(gtf_file, sorted_bam_list)
