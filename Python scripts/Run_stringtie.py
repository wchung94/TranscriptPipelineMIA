#!/user/bin/env python3
"""
Author: Lukas van den Heuvel
Student number: 951218334040 Usage:
n.a.
Description: Subroutine template for command-line tools
"""

# Import statements
from datetime import datetime
from sys import argv
import subprocess
import os.path

def infile_check(filename):
    """Check the presence of your infile, return -1 if the file does not exist.

    Keywords:
    filename -- filename to check
    return: 0 if the script succesfull, other integers in case of an error """

    # check if the input files exist, if not return an error exit()
    if os.path.exists(filename):
        return(0)
    else:
        print("there is no such file: {}".format(filename))
        if auto_mode == "manual":
            return -1
        elif auto_mode == "overwrite" or auto_mode == "use_existing":
            return -1 #Could implement a system that fixes this user error.

def outfile_check(filename):
    """Get information on the filename availability. Overwrite it, use it or give an error.

    Keywords:
    filename -- filename to check
    return: 0 if the script succesfull, other integers in case of an error """

    # check if the output file already exists, if not run the program
    if os.path.exists(filename):
        print("{0} already exists\n".format(filename))
        question_loop = True
        if auto_mode == "manual":
            while question_loop == True:
                overwrite_question = input("Would you like to overwrite it? [y/n]")
                if overwrite_question == "y":
                    subprocess.check_call("rm {}".format(filename))
                    question_loop = False
                    return(0)
                elif overwrite_question == "n":
                    return(-1)
                else:
                    print("Please provide an answer as 'y' (yes) or 'n' (no)")
        elif auto_mode == "overwrite":
            subprocess.check_call("mv {} {}.old".format(filename,filename))
            return 0
        elif auto_mode == "use_existing":
            return 0
    else:
        return(0)

def call_error(codeblock,message):
    """Provide the user with a short error report, then exit the script.

    Keywords:
    codeblock -- block of code where the error occured
    message -- additional error message
    return: N/A """

    print(datetime.now(),"The pipeline was interrupted due to an error in the {} block".format(codeblock))
    print("ERROR: {}".format(message))
    exit()

def run_stringtie(label,gff_file, cores=8):
    """Run stringtie on a single .bam file, return a .gtf file.

    Keywords:
    label -- read-ID of the pipeline
    gff_file -- reference genome annotation file (or path)
    cores -- number of threads to give stringtie (default is 8)
    return: N/A, the function writes to a file """

    #set variables
    input_filename = label + ".bam"
    output_filename = label + ".gtf"

    if infile_check(input_filename) != 0:
        call_error("assembly","input file {} not found".format(input_filename))
    if outfile_check(output_filename) != 0:
        call_error("assembly","file {} already exists. [Overwrite = No]".format(output_filename))

    cmd_command = "stringtie {} -G {} -o {} -p {} -l {}".\
    format(input_filename,gff_file,output_filename,cores,label)

    if subprocess.check_call(cmd_command, shell=True) != 0:
        call_error("assembly","unexpected error while running 'stringtie'")
        
def run_stringtie(label_list):
    """Run stringtie --merge on a list of files

    Keywords:
    label -- read-ID of the pipeline
    gff_file -- reference genome annotation file (or path)
    cores -- number of threads to give stringtie (default is 8)
    return: N/A, the function writes to a file """

    #set variables
    
    
    gtf_list = []
    for label in label_list:
        gtflist.append("{}.gtf".format(label))
    
    ###FINISH THIS CODE YOU DUMB IDIOT!!!###

# main
if __name__ == '__main__':
    auto_mode = 'overwrite' #EITHER: "manual", "use_existing" or "overwrite"
    annotation_path = "/local/data/course/project/genomes/Catharanthus_roseus/genome_version2/cro_v2.gene_models.gff3"

    label_list = argv[1:]
    
    for label in label_list:
        print(datetime.now(), "Starting assembly of {}.".format(label))
        run_stringtie(label,annotation_path)
        print(datetime.now(), "Assembly of {} successful.".format(label))

    run_stringtieMerge(label_list)
