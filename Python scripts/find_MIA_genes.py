from sys import argv

def create_gene_list(input_file):
    """Function to parse through gene expression file to find all DEG genes

    :param input_file: Gene_expression file
    :return: gene_list: Gene list of only significant genes (excluding T)
    """
    handle_file = open(input_file, 'r')
    gene_list = []
    for line in handle_file:
        line = line.strip().split('\t')
        if line[-1]=='yes':
            if ',' in line[2]:
                genes = line[2].split(',')
                gene_list += genes
            else:
                gene_list += [line[2]]
    handle_file.close()
    return gene_list

def create_proper_gene_list(gene_list):
    """Function to create the gene list with names containing T

    :param gene_list: List of genes that do not contain the letter T
    :return: gene_list2: List of genes that do contain the letter T
    """
    gene_list2 = []
    for gene in gene_list:
        position_dash = gene.index('_')
        new_gene_name = gene[:position_dash+1] + 'T' + gene[position_dash+1:]
        gene_list2 += [new_gene_name]
    return gene_list2

def parse_CRO_names(CRO_input_file):
    """Function to transform text file into a dictionary

    :param CRO_input_file: Self-uploaded file from MIA excel file. Contains gene name (with T) and description
    :return: CRO_dict: Dictionary of the gene name with description as value
    """
    CRO = open(CRO_input_file, 'r')
    CRO.readline()
    CRO_dict = {}
    for line in CRO:
        CRO_name = line.strip().split('\t')[0]
        gene_name = line.strip().split('\t')[1]
        CRO_dict[CRO_name] = gene_name
    CRO.close()
    return CRO_dict

def find_MIA_names(dic, gene_list):
    """Function to create dictionary of only MIA genes (without letter T) with description

    :param dic: CRO_dictionary earlier created
    :param gene_list: List of genes that do contain the letter T
    :return: final_MIA: Dictionary containing gene name without letter T and description
    """
    found_MIA = {}
    final_MIA = {}
    for gene in gene_list:
        if gene in dic:
            found_MIA[gene] = dic[gene]
    for key,items in found_MIA.items():
        position_T = key.index('T')
        new_key = key[:position_T]+key[position_T+1:]
        final_MIA[new_key] = items
    return final_MIA

def save_file(final_MIA_dictionary, MIA_file):
    """Function to create a file that contains gene name and description (later used for final table)

    :param final_MIA_dictionary: Dictionary of the gene name (without T) with description as value
    :param MIA_file: Name of the output file
    """
    handle_file = open(MIA_file,'w')
    for key, item in final_MIA_dictionary.items():
        handle_file.write('{}\t{}\n'.format(key,item))
    handle_file.close()

def find_exp(final_MIA, gene_exp):
    """Function to create list of cuffdiff gene_exp lines of MIA genes

    :param final_MIA: Dictionary of the gene name (without T) with description as value
    :param gene_exp: Cuffdiff output gene_expression file
    :return table: List containing lines of cuffdiff gene_exp of MIA genes
    """
    handle_file = open(gene_exp, 'r')
    table = []
    for key in final_MIA:
        handle_file.seek(0)
        for line in handle_file:
            if line.startswith(key):
                table += [line.strip()]
    handle_file.close()
    return table

def find_exp2(final_MIA, gene_exp, mia_output):
    """Function to write file in gene_exp cuffdiff format only of MIA genes (later used for final table)

    :param final_MIA: Dictionary of the gene name (without T) with description as value
    :param gene_exp: Cuffdiff output gene_expression file
    :param mia_output: output file with full gene_exp lines of MIA genes
    """
    handle_file = open(gene_exp, 'r')
    mia_output = open(mia_output,'w')
    mia_output.write(handle_file.readline())
    for key in final_MIA:
        handle_file.seek(0)
        for line in handle_file:
            if line.startswith(key):
                mia_output.write(line.strip() + '\n')
    handle_file.close()
    mia_output.close()

def create_new_table(table, mia_output):
    """Function to create simple table of gene name and the three significance values of CK/MJ, CK/ET and MJ/ET

    :param table: List containing lines of cuffdiff gene_exp of MIA genes
    :param mia_output: Name of output file that will contain a header and gene name and three significance values of
    CK/MJ, CK/ET and MJ/ET. Can be used for Venn Diagram
    """
    mia_file = open(mia_output, 'w')
    mia_dic = {}
    sig_list = []
    old_gene = None
    for line in table:
        line_list = line.split('\t')
        new_gene = line_list[0]
        if old_gene == new_gene or old_gene == None:
            sig_list += [line_list[-1].strip()]
            old_gene = new_gene
        else:
            mia_dic[old_gene] = sig_list
            sig_list = [line_list[-1].strip()]
            old_gene = line_list[0]
    if sig_list:
        mia_dic[old_gene] = sig_list
    mia_file.write('{}\t{}\t{}\t{}\n'.format('Gene', 'SigCK_MJ', 'SigCK_ET', 'SigMJ_ET'))
    for key, items in mia_dic.items():
        mia_file.write('{}\t{}\t{}\t{}\n'.format(key, items[0], items[1], items[2]))
    mia_file.close()


gene_exp_file = argv[1]
CRO_input = argv[2]
#Filename of file which will contain whole gene_exp data from MIA genes. USED IN FINAL TABLE
MIA_output1 = argv[3]
#Filename of file whichw ill contain only gene + significance
MIA_output2 = argv[4]
#File that will contain gene name and description. USED IN FINAL TABLE
MIA_dictionary_file = argv[5]

#Create list of genes in two formats
gene_list = create_gene_list(gene_exp_file)
final_gene_list = create_proper_gene_list(gene_list)

#Parse through CRO_file
CRO_dic = parse_CRO_names(CRO_input)
MIA_dic  = find_MIA_names(CRO_dic, final_gene_list)

#Create new files of mia genes
mia_table = find_exp(MIA_dic, gene_exp_file)
find_exp2(MIA_dic,gene_exp_file,MIA_output1)
create_new_table(mia_table, MIA_output2)
save_file(MIA_dic, MIA_dictionary_file)