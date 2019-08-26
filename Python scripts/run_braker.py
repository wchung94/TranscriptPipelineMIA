#! /usr/bin/env python3
'''
Run Braker2 using pipeline B (Fa file and Bam files) 
'''

# GLOBAL FASTQ DATA LIST

from sys import argv
import os
import subprocess

def select_bam_files(bam_directory):
    """ Select all sorted BAM files from designated directory
   
    Keywords:
    bam_directory: string, directory name of the stored BAM files
    returns:
    sorted_bam_list: list of strings, file names of all sorted BAM files 
    
    WY Chung
    """    
    
    sorted_bam_list = []
    for file in os.listdir(bam_directory):
        if file.endswith('.bam') and not file.endswith('_unsorted.bam'):
          sorted_bam_list.append(bam_directory + '/' + file)
    sorted_bam_list = sorted(sorted_bam_list)
    
    return sorted_bam_list


def run_braker(fa_file, bam_file_list):
    """ Run Braker2 program on FASTA file and BAM files with the tools:
    GeneMark and Augustus.
   
    Keywords:
    fa_file: string, filename of FASTA of the reference genome
    bam_file_list: list of strings, filename of input BAM files from Hisat
    returns:
    out_gff: string, filename of output .gff file of Augustus from Braker2 
    
    WY Chung, Nathalie de Vries
    """
    
    # check if the output_file already exists
    if not os.path.exists(fa_file):
        print('FASTA file does not exist')
        return
    for file in bam_file_list:
        if not os.path.exists(file):
            print('bam file(s) do(es) not exist')
            return
                
    species = 'Catharanthusroseus'
    bam_files = ','.join(bam_file_list)
    #maybe make the working directory before this working dir? but braker makes a direcoctory in working directory by itself
    working_dir = '/local/data/BIF30806_2018_2/project/groups/nijveen/braker_results'
    
    #Changes to the copied local/prog/braker
    os.chdir('/local/data/BIF30806_2018_2/project/groups/nijveen/python_script/braker')
    
    cmd_braker = "./braker.pl -cores=8\
     --species={} --genome={} --bam={}\
     --AUGUSTUS_CONFIG_PATH=/local/data/course/project/groups/nijveen/python_script/augustus/config\
     --GENEMARK_PATH=/local/prog/gm_et_linux_64/gmes_petap\
     --BAMTOOLS_PATH=/usr/bin\
     --workingdir={} --useexisting"\
    .format(species, fa_file, bam_files, working_dir)
    subprocess.check_output(cmd_braker, shell=True) #e stores the output
    subprocess.check_call(cmd_braker, shell=True)

    out_gff = working_dir + '/braker/' + species + '/augustus.gff'
    os.chdir('/local/data/BIF30806_2018_2/project/groups/nijveen/python_script')
    
    #check if the augustus output gff file exists
    if not os.path.exists(out_gff):
        print('Braker failed.')
        return 
    
    return out_gff

if __name__ == '__main__':
    #Select sorted bam files from directory
    bam_directory = '/local/data/BIF30806_2018_2/project/groups/nijveen/Pipeline_new/bam_files'
    sorted_bam_list = select_bam_files(bam_directory)
    #print(sorted_bam_list)
    #run braker2 on selected bam files
    fa_file = '/local/data/BIF30806_2018_2/project/genomes/Catharanthus_roseus/genome_version2/cro_v2_asm.fasta'
    out_gff = run_braker(fa_file, sorted_bam_list)
    


