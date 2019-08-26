#! /usr/bin/env python3
'''
Run blast on the DEG  
'''

# GLOBAL FASTQ DATA LIST

from sys import argv
import os
import subprocess

def run_blast(seq_db_fasta, deg_fasta):
    """Executes blast by creating own database
    return: None
    """

    # check if the output_file already exists
    os.chdir('/local/data/BIF30806_2018_2/project/groups/nijveen/Pipeline_new/blast')
    
    if not os.path.exists(seq_db_fasta):
        print('fa file does not exist')
        return
    if not os.path.exists(deg_fa):
        print('fa file does not exist')
        return
        
    cmd_create_db = "makeblastdb -in {} -dbtype {} -out catha"/
    .format(seq_db_fasta, 'nuc1')
    cmd_blast = "blastn -db catha -evalue 1E-05 -num_alignments 1 -num_descriptions 1 -query {} -out results.out "/
    .format(deg_fasta)

    e = subprocess.check_output(cmd_blast, shell=True) #e stores the output
    subprocess.check_call(cmd_create_db, shell=True)
    subprocess.check_call(cmd_blast, shell=True)
    return

if __name__ == '__main__':
    #Select fasta files from directory
    sequence_db_fasta = '/local/data/BIF30806_2018_2/project/groups/nijveen/Pipeline_new/blast/sequence.fasta'
    deg_fasta = '/local/data/BIF30806_2018_2/project/groups/nijveen/Pipeline_new/blast/DEG_file.fasta'
        
    #run blast on selected
    run_blast(sequence_db_fasta, deg_fasta)



