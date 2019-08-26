#! /usr/bin/env python3
'''
Run cuffdiff 2 using a reference gtf or gff file and multiple bam files
Erik, Charlie, Lukas
'''
from sys import argv
import os
import subprocess

def run_cuffdiff(gtf_reference, label_dictionary, bamdir, cuffdiff_output_dir):
    """Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do
    eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim
    ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut
    aliquip ex ea commodo consequat. Duis aute irure dolor in
    reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla
    pariatur. Excepteur sint occaecat cupidatat non proident, sunt in
    culpa qui officia deserunt mollit anim id est laborum.
    """
    #category_list = [s[0] for s in label_dictionary.values()]
    category_dict = {}
    bam_string = ""
    temp_list = []
    temp_string = ""
    
    
    for label, category in label_dictionary.items():
        category = category[0]
        category_dict[category] += label
    
    for category in category_dict.keys():
        for bam_label in category_dict[category]:
            bam_label = "../{}/" + s + ".bam".format(bamdir)
            temp_list += bam_label
        temp_string = ",".join(temp_list)
        bam_string += tempstring + " "
        temp_list = []
        temp_string = ""
    
    cmd_cuffdiff = "cuffdiff -o {} -p 16 -u -no-update-check {} {}"\
    .format(cuffdiff_output_dir, gtf_reference, bam_string)
    
    e = subprocess.check_output(cmd_cuffdiff, shell=True) #e stores the output
    print(datetime.now() + '\n cuffdiff process start')
    subprocess.check_call(cmd_cuffdiff, shell=True)
    print(datetime.now() + '\n cuffdiff process finished')
    return


if __name__ == '__main__':
    
    #cuffdiff setup
    os.chdir('..')
    try:
        os.mkdir('cuffdiff')
    except:
        #check if 
        pass #not propper python way. Should be written to a log file.
    os.chdir('cuffdiff')
    
    #cuffdiff
    bamdir = ""
    cuffdiff_out_dir = "cuffdiff_out"
    run_cuffdiff(gtf_reference, label_dictionary, bamdir, cuffdiff_out_dir)
