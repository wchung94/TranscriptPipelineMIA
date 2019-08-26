from sys import argv

def create_table(gene_exp_file_mia, found_genes, table_file):
    """Script to write a final tab delimited file as table for gene description of MIA genes

    :param gene_exp_file_mia: Self-created gene expression text file of only mia genes. (save_file in find_MIA_genes.py
    mia_lines2.txt)
    :param found_genes: tab delimited file with genename and description (mia_dictionary_table.txt)
    :param table_file: output file that is created
    """
    table_output = open(table_file, 'w')
    gene_exp_mia = open(gene_exp_file_mia,'r')
    found_genes_file = open(found_genes, 'r')
    table_output.write('Gene\tDescription\tUp/Downregulated\n')
    for line in found_genes_file:
        line_list1 = line.strip().split('\t')
        gene = line_list1[0]
        gene_exp_mia.seek(0)
        for line in gene_exp_mia:
            if line.startswith(gene):
                line_list = line.strip().split('\t')
                if line_list[-1] == 'yes':
                    sample1 = line_list[4].replace('q1','Control').replace('q2', 'MeJA')
                    sample2 = line_list[5].replace('q2','MeJA').replace('q3','Ethylene')
                    value1 = float(line_list[7])
                    value2 = float(line_list[8])
                    if value2 > value1:
                        table_output.write('{}\t{}\tUpregulated in {} compared to {}\n'.format(gene, line_list1[1], sample2,sample1))
                    else:
                        table_output.write('{}\t{}\tDownregulated in {} compared to {}\n'.format(gene, line_list1[1], sample2,sample1))

found_genes = argv[1]
gene_exp = argv[2]
table_out = argv[3]
create_table(gene_exp, found_genes, table_out)