#!/user/bin/env python3
"""
Author: Vera, Nathalie, Erik, Lukas van den Heuvel, WY
Student number: Usage:
n.a.
Description: Work flow for command-line tools
"""

# Import statements
from datetime import datetime
from sys import argv
import os
import subprocess

#file checker

def infile_check(filename):
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
    print("The pipeline was interrupted due to an error in the {} block".\
    format(codeblock))
    print("ERROR: {}".format(message))
    exit()


#Fastq splitter

def select_fastq_files(directory='fastq_sequences'):
    sorted_fastq_list = []
    for fastq_file in os.listdir(directory):
        if fastq_file.endswith('.fastq') and not fastq_file.endswith('_1.fastq') and not fastq_file.endswith('_2.fastq'):
          sorted_fastq_list.append('fastq_sequences/' + fastq_file)
    sorted_fastq_list = sorted(sorted_fastq_list)
    return sorted_fastq_list

def get_fastq_rec(fastq_file):
    """Generator that yields FASTQ records

    fastq_file -- an opened FASTQ file
    record -- a list of 4 strings (lines) for a single FASTQ record
    """
    with open(fastq_file) as handle_file:
        record = []
        for counter, line in enumerate(handle_file):
            record.append(line.strip())
            if counter % 4 == 3:
                yield record
                record = []

def create_header(old_header,forward_or_reverse, length):
    header_list = old_header.strip().split(' ')
    part_one = header_list[0] + '.{}'.format(forward_or_reverse)
    full_header = part_one + ' length={}'.format(length)
    return full_header

def reverse_compliment(sequence):
    reverse = sequence[::-1]
    translate_table =  str.maketrans("ACGTN", "TGCAN")
    return reverse.translate(translate_table)

def fastq_split(input_file):
    file_name = input_file.split('.')[0]
    file_for = open(file_name + '_1.fastq', 'w')
    file_rev = open(file_name + '_2.fastq', 'w')
    for record in get_fastq_rec(input_file):
        file_for.write(create_header(record[0], 'for', '90') + '\n')    #header
        file_for.write(record[1][:90] + '\n')    #seq
        file_for.write(create_header(record[2], 'for', '90') + '\n')    #+header
        file_for.write(record[3][:90] + '\n')    #q-scores

        file_rev.write(create_header(record[0], 'rev', '90') + '\n') #header
        file_rev.write(reverse_compliment(record[1][90:]) + '\n') #reverse_complement (script)
        file_rev.write(create_header(record[2], 'rev', '90') + '\n') #+header
        file_rev.write(record[3][90:][::-1] + '\n') #reverse [::-1]
    return


#Hisat runner to create sam files

def select_split_fastq(directory='fastq_sequences', split):
    sorted_fastq_list = []
    for fastq_file in os.listdir(directory):
        if fastq_file.endswith(split + '.fastq'):
          sorted_fastq_list.append('fastq_sequences/' + fastq_file)
    sorted_fastq_list = sorted(sorted_fastq_list)
    return sorted_fastq_list


def hisat_index():
    index_genome = 'index_genome'
    hisat2_cmd = 'hisat2-build chrX.exon chrX_data/genome/chrX.fa chrX_tran'
    
    if not os.path.exists(index_genome):
        subprocess.check_call(hisat2_cmd, shell=True)
    else:
        return index_genome
    
def run_hisat2(index_file, input_for, input_rev, out_hisat2):
    """ Run HISAT2 program on FASTQ file
    
    index_file: string, filename of index file of the organism of the input files
    input_for: string, filename of input FASTA/FASTQ file of forward sequence 
    input_rev: string, filename of input FASTA/FASTQ file of reverse sequence 
    out_hisat2: string, filename of output .sam file of hisat2 
    """

    cmd = "hisat2 --dta -x {} -1 {} -2 {} -S {}"\
    .format(index_file, input_for, input_rev, out_hisat2)
    if not os.path.exists (out_hisat2):
        process_hisat2 = subprocess.check_call(cmd, shell=True)
    else:
        print("File already exists")
    

#stringtie runner

def run_stringtie(label,gff_file, cores=8):
    """Run the program ...
    output_filename -- the name the output should be saved to
    return: N/A, the function writes to a file
    """
    #set variables
    input_filename = label + ".bam"
    output_filename = label + ".gtf"

    if infile_check(input_filename) != 0:
        call_error("assembly","input file not found")
    if outfile_check(output_filename) != 0:
        call_error("assembly","file already exists. [Overwrite = No]")

    cmd_command = "stringtie {} -G {} -o {} -p {} -l {}".\
    format(input_filename,gff_file,output_filename,cores,label)

    if subprocess.check_call(cmd_command, shell=True) != 0:
        call_error("assembly","unexpected error while running 'stringtie'")


#cuffdiffrunner

def select_bam_files(directory='bam_files'):
    sorted_bam_list = []
    for bam_file in os.listdir(directory):
        if bam_file.endswith('.bam') and not bam_file.endswith('_unsorted.bam'):
          sorted_bam_list.append('bam_files/' + bam_file)
    sorted_bam_list = sorted(sorted_bam_list)
    return sorted_bam_list

def run_cuffdiff(gtf_file, bam_file_list):
    """Executes cuffdiff with an input_file ad save output in output_file
    bam_file, gtf_file: string, location and name of input file
    return: None
    """

    # check if the input files exists
    if infile_check(input_filename) != 0:
        call_error("assembly","input file not found")
    for bam_file in bam_file_list:
        if infile_check(bam_file) != 0:
            call_error("assembly","input file not found")
        

    cmd_cuffdiff = "cuffdiff -o {} -p 4 --library-type fr-firststrand -u -no-update-check {} {},{} {},{},{} {},{},{}"\
    .format('cuffdiff_output_directory_CR', gtf_file, bam_file_list[0], bam_file_list[7], bam_file_list[1], bam_file_list[2], bam_file_list[5], bam_file_list[3], bam_file_list[4], bam_file_list[6])
    e = subprocess.check_output(cmd_cuffdiff, shell=True) #e stores the output
    print('process start')
    subprocess.check_call(cmd_cuffdiff, shell=True)
    return


# main
if __name__ == '__main__':
 
    #initial defining of map with fastq sequences
    fastq_dir = 'fastq_sequences' # or define with ´fastq_dir = argv[1]´
    
    #fastQ parser and splitter
    fastq_files_list = select_fastq_files(fastq_dir)
    for input_file in fastq_files_list:
        fastq_split(input_file)
     
    #run tool HISAT2 to create sam files from splitted fastq files
    index = hisat_index() #index_genome/index_genome
    in_for_list = select_split_fastq(fastq_dir,'_1')
    in_rev_list = select_split_fastq(fastq_dir,'_2')
    out_file = argv[4]
    for input_file in range(len(in_for_list)):
        run_hisat2(index, in_for_list[input_file], in_rev_list[input_file], out_file)
    
    #run samtools to create bam files from sam files
    
    
    
    
    #run Stringtie  
    auto_mode = 'overwrite' #EITHER: "manual", "use_existing" or "overwrite"
    annotation_path = "/local/data/course/project/genomes/Catharanthus_roseus/genome_version2/cro_v2.gene_models.gff3"

    for label in argv[1:]:
        print(datetime.now(),"Starting assembly of {}.".format(label))
        run_stringtie(label,annotation_path)
        print(datetime.now(),"Assembly of {} successful.".format(label))


   #Select sorted bam files from directory
    sorted_bam_list = select_bam_files()
    print(sorted_bam_list)
    #run cuffdiff on selected bam files
    gtf_file = 'stringtie_gtfs/stringtie_merged.gtf'
    run_cuffdiff(gtf_file, sorted_bam_list)


