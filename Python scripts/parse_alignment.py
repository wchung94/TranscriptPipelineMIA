from sys import argv

def create_gene_dic(input_file):
    gene_dic = {}
    alignment_file = open(input_file, 'r')
    for line in alignment_file:
        if line.startswith('Query='):
            label = line.strip().partition(' ')[2]
            gene_dic[label] = []
        if line.startswith('>'):
            if 'TSA' in line:
                continue
            else:
                gene_name = line.strip().partition(' ')[2]
                #print(gene_name)
                gene_dic[label] += [gene_name]
    return gene_dic

def final_dic(old_dic):
    final_dic = {}
    for key,item in old_dic.items():
        if item == []:
            continue
        else:
            final_dic[key] = item
    print(final_dic)

input_file = argv[1]
first_dic = create_gene_dic(input_file)
final_dic(first_dic)
