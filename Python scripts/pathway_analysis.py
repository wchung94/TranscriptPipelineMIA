#! /usr/bin/env python3
'''
Pathway enrichment analysis on the significant expressed MIA genes.
'''

#Import statements
from sys import argv 
import ast
import os
import numpy
import scipy
import scipy.cluster.hierarchy as sch
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns


def parse_file(txt_file):
        
    data = pd.read_csv(txt_file, sep='\t')
    data = data.iloc[:,[0,4,5,9,11,13]]
    data.columns = ['gene','group_1','group_2','log_fold','p-val','significant']
    return data   
    
def dict_to_df(dict_file):
    with open(dict_file,'r') as data:
        for line in data:
            gene_dict = line
    
    gene_dict = ast.literal_eval(gene_dict)
    dict_df = pd.DataFrame.from_dict(gene_dict, orient='index')
    dict_df.transpose()
    dict_df.columns = ['gene','description','protein','significant']
    return dict_df



if __name__ == "__main__":
    

    
    

    

    txt_file = '/local/data/BIF30806_2018_2/project/groups/nijveen/Pipeline_new/blast/mia_lines2.txt'
    genes_df = parse_file(txt_file)
    print(genes_df)
    
    dict_txt = '/local/data/BIF30806_2018_2/project/groups/nijveen/Pipeline_new/blast/final_dictionary.txt'
    dict_df = dict_to_df(dict_txt) 
    print(dict_df)


