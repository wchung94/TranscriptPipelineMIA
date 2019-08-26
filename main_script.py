#!/user/bin/env python3
"""
Nijveen group
Authors: Vera Claasen, Wing Yu Chung, Lukas van den Heuvel, Erik ijland,
Nathalie de Vries
Student numbers: 970825157120, 940201157090, 951218334040, 960504980110,
920107916130
Usage: n.a.
Description: Work flow for command-line tools of ..
Version: 1.0 December 2018
"""

# Import statements
from datetime import datetime
import os
import subprocess
import argparse
import time
from argparse import RawTextHelpFormatter
from datetime import datetime
from multiprocessing import Process
from find_MIA_genes import create_gene_list, create_proper_gene_list, parse_CRO_names, find_MIA_names, save_file, find_exp, find_exp2, create_new_table

########################################################################
####################        Error handling        ######################
########################################################################

def infile_check(filename):
    """ Checks if the input file exists, if not returns an error message

    Kewords:
    filename: string, filename of input file

    Lukas van den Heuvel
    """
    if os.path.exists(filename):
        return(0)
    else:
        print(datetime.now(), "there is no such file: {}".format(filename))
        if auto_mode == "manual":
            return -1
        elif auto_mode == "overwrite" or auto_mode == "use_existing":
            return -1 #Could implement a system that fixes this user error.

def outfile_check(filename):
    """ Checks if the output file exists, if not runs the program

    Keywords:
    filename: string, filename of input file

    Lukas van den Heuvel
    """
    if os.path.exists(filename):
        print(datetime.now(), "{0} already exists".format(filename))
        question_loop = True
        if auto_mode == "manual":
            while question_loop == True:
                overwrite_question = input(str(datetime.now()) + " Would you like to overwrite it? [y/n]")
                if overwrite_question == "y":
                    subprocess.check_call("rm {}".format(filename), shell=True) #ADD SOME SORT OF PRINT/LOG FOR THIS!!!
                    question_loop = False
                    return 0
                elif overwrite_question == "n":
                    return -1
                else:
                    print(datetime.now(), "Please provide an answer as 'y' (yes) or 'n' (no)")
        elif auto_mode == "overwrite":
            subprocess.check_call("mv {} {}.old".format(filename,filename), shell=True)
            return 0
        elif auto_mode == "use_existing":
            return -1
    else:
        return(0)

def call_error(codeblock,message):
    """ Returns location of an error and error message

    Keywords:
    codeblock: string, block of code were error occurs
    message: string, text message of an error

    Lukas van de Heuvel
    """
    print(datetime.now(), "The pipeline was interrupted due to an error in the {} block".format(codeblock))
    print(datetime.now(), "ERROR: {}".format(message))
    exit()

########################################################################
####################        User Interface        ######################
########################################################################

def select_input_files(directory):
### labelfile as parameter!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    """ Let's the user select and label datasets out of given directory
    (manualy or automatted with labelfile)

    Keywords:
    directory: string, name of input directory
    labelfile: string, name of text file with dataset name and nicknames
    Retuns:
    label_dict: dictionary of {label: dataset name, nickname}

    Lukas van den Heuvel
    """
    filelist = os.listdir(directory)
    selected_list = []

    if "dataset_labels.txt" in filelist:
        label_dict = {}
        print("found 'dataset_labels.txt' in your directory. Using this instead of manual input.")
        with open(directory+'/dataset_labels.txt') as dataset_labels:
            for line in dataset_labels:
                line = line.strip().split('\t')
                label_dict[line[0]] = [line[1],line[2]]
    else:
        print("There are {} files in this directory:".format(len(filelist)))
        print('\t'.join(filelist))

        for filename in filelist:
            if yes_no("include {} in the analysis? [y/n] ".format(filename)):
                selected_list.append(filename)

        print("You selected {} files for the analysis:".format(len(selected_list)))

        label_dict = {}
        dataset_count = {}
        for selected_name in selected_list:
            dataset_name = input("Which dataset does input {} belong to? ".format(selected_name))

            if dataset_name in dataset_count.keys():
                dataset_count[dataset_name] += 1
            else:
                dataset_count[dataset_name] = 1

            label = selected_name.replace('.fastq','').strip()
            label_dict[label] = [dataset_name, dataset_name + str(dataset_count[dataset_name])]

            print("added {} to dataset {} with nickname {}".format(selected_name, dataset_name, dataset_name + str(dataset_count[dataset_name])))
    return label_dict

def yes_no(inputstring):
    """ Checks if user input is y/n, if not returns error message

    Keywords:
    inputstring: string, user input (y/n)

    Lukas van den Heuvel
    """
    question_loop = True
    while question_loop == True:
        answer = input(inputstring)
        if answer.strip() == "y":
            question_loop = False
            return(True)
        elif answer.strip() == "n":
            return(False)
        else:
            print("\nERROR: Please provide an answer as 'y' (yes) or 'n' (no)")

def print_dictionary(mydict):
    """ Prints each key and value(s) of dictionary on tab delimeted line

    Keywords:
    mydict: dictionary of {key: value(s)}

    Lukas van den Heuvel
    """
    for key, value in mydict.items():
        print(key,"\t","\t".join(value))
    print("")

def set_arguments():
    """ Set arguments for command line of main script

    Lukas van den Heuvel
    """
    parser.add_argument("inputdirectory",\
    help="The directory containing your RNAseq transcripts")

    parser.add_argument("--labelfile",\
    help="list containing files-of-interest and their categories/groups",\
    metavar=".txt") #NOT YET IMPLEMENTED

    parser.add_argument("--genome",\
    help="Path to your reference genome from your present working directory.\n\tDefaults to '/local/data/course/project/genomes/Catharanthus_roseus/genome_version2/cro_v2_asm.fasta'",\
    default= '/local/data/course/project/genomes/Catharanthus_roseus/genome_version2/cro_v2_asm.fasta',\
    metavar=".fasta")

    parser.add_argument("--annotation",\
    help="Path to your reference annotation from your present working directory.\n\tDefaults to '/local/data/course/project/groups/nijveen/FINAL/augustus.gff'",\
    default= '/local/data/course/project/groups/nijveen/FINAL/augustus.gff',\
     metavar=".gff")

    parser.add_argument("--mode",\
    help="When encountering an existing output file:\n\t'manual'\tasks user per incident.[default]\n\t'overwrite'\toverwrites existing files.\n\t'use_existing'\tuses existing files.",\
    default="manual",\
    metavar="")

    parser.add_argument("--paired",\
    help="use paired reads instead of forward reads",\
    action="store_true")

    parser.add_argument("--hisatindex",\
    help="Path to hisat2 index files.\n\tDefaults to '/local/data/course/project/groups/nijveen/index_genome/'",\
    default="/local/data/course/project/groups/nijveen/index_genome/index_genome",\
    metavar=".ht2")

    #parser.add_argument("--version", help="print the version number of this script",action="store_true")

########################################################################
####################        FASTQ Handling        ######################
########################################################################

def fastq_split(label):
    """ Splits off forward strand (90 bp) of paired FASTQ sequence (180 bp)

    Keywords:
    label: string, filename of paired FASTQ sequence
    Returns:
    fastq_out: string, filename of forward FASTQ sequence

    Lukas van den Heuvel
    """
    fastq_in = open(label + '.fastq', 'r')
    fastq_out = open(label + '_1.fastq', 'w')

    for linenumber, line in enumerate(fastq_in):
        mod = linenumber%2
        if mod == 0:
            fastq_out.write("{}.for length=90\n".format(line.strip().split(" ")[0]))
        if mod == 1:
            fastq_out.write("{}\n".format(line[:90]))
    print(datetime.now(), "FASTQ-splitting complete for", label)

def check_presence(label_dict, extension):
    """ Checks if the input file exists, if not returns an error exit()

    Keywords:
    label_dict: dictionary of {label: dataset name, nickname}
    extension: string, part added to label which completes filename

    Lukas van de Heuvel
    """
    for label in label_dict:
        if label + extension in os.listdir("."):
            print(datetime.now(), "Self-check succesfull. Continuing.")
            return 0
        else:
            print(datetime.now(), "ERROR: Files not found.")
            return -1

########################################################################
####################     Commandline functions    ######################
########################################################################

def run_hisat2(index_file, label_dict, fastqdir):
#### change filename to hisat_out (better readability
    """ Runs HISAT2 program on FASTQ sequence file

    Keywords:
    index_file: string, filename of index file of the used organism
    label_dict: dictionary of {label: dataset name, nickname}
    fastqdir: directory of forward FASTQ sequences
    Returns:
    filename:

    Nathalie de Vries, Vera Claasen, Lukas van den Heuvel
    """
    temp_infile_list = [s for s in label_dict.keys()]
    #"../{}/"+ s + "_1.fastq".format(fastqdir)

    for label in label_dict:
        filename = label + ".sam"
        if outfile_check(filename) != 0:
            temp_infile_list.pop(temp_infile_list.index(label))

    for templabel in temp_infile_list:
        cmd = "hisat2 -p 8 --dta-cufflinks -x {} -U {} -S {}".format(index_file, "{}/{}/".format(starting_directory,fastqdir) + templabel + "_1.fastq", templabel + ".sam")
        hisat_return = subprocess.check_call(cmd, shell=True)
        if hisat_return != 0:
            call_error("FASTQ Allignment","an unexpected error occured while running HISAT2 on {}".format(templabel))

def run_samview(label):
### label dictionary as parameter!!!!!!!!!!!! part out of main in function
    """ Run samtools view program on .sam file

    Keywords:
    label_dict: dictionary of {label: dataset name, nickname}
    Returns:
    unsorted_bam: string, filename of output unsorted .bam file

    Nathalie de Vries, Vera Claasen, Lukas van den Heuvel
    """
    in_sam = label + ".sam"
    unsorted_bam = label + "_unsorted.bam"
    cmd = " samtools view -bS {} > {}".format(in_sam, unsorted_bam)     #PRINT: RUNNING SAMVIEW FOR [LABEL] ON PROCESS [ProcessID]
    samview_return = subprocess.check_call(cmd, shell=True)
    if samview_return != 0:
            call_error("SAM-BAM Conversion","an unexpected error occured while converting .sam to .bam for {}".format(label))


def run_samsort(label):
### label dictionary as parameter!!!!!!!!!!!! part out of main in function
    """ Run samtools sort program on unsorted .bam file

    Keywords:
    unsorted_bam: string, filename of output unsorted .bam file
    sorted_bam: string, filename of output sorted .bam file

    """
    unsorted_bam = label + "_unsorted.bam"
    sorted_bam = label
    cmd = " samtools sort {} {}".format(unsorted_bam, sorted_bam)
    samsort_return = subprocess.check_call(cmd, shell=True)
    if samsort_return != 0:
            call_error("BAM Sorting","an unexpected error occured while sorting allignment {}".format(label))

def run_cuffdiff(gtf_reference, label_dict, bampath):
    """Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do
    eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim
    ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut
    aliquip ex ea commodo consequat. Duis aute irure dolor in
    reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla
    pariatur. Excepteur sint occaecat cupidatat non proident, sunt in
    culpa qui officia deserunt mollit anim id est laborum.
    """
    category_dict = {}
    bam_string = ""
    temp_list = []
    temp_string = ""
    cuffdif_outfiles = ["bias_params.info", "cds_exp.diff",\
    "genes.fpkm_tracking", "isoforms.count_tracking", "promoters.diff",\
    "splicing.diff", "tss_groups.fpkm_tracking", "cds.count_tracking",\
    "cds.fpkm_tracking", "gene_exp.diff", "genes.read_group_tracking",\
    "isoforms.fpkm_tracking", "read_groups.info", "tss_group_exp.diff",\
    "tss_groups.read_group_tracking", "cds.diff",\
    "cds.read_group_tracking", "genes.count_tracking",\
    "isoform_exp.diff", "isoforms.read_group_tracking", "run.info",\
     "tss_groups.count_tracking", "var_model.info"]

    for label, category in label_dict.items():
        category = category[0]
        if category in category_dict.keys():
            category_dict[category] += [label]
        else:
            category_dict[category] = [label]

    for category in category_dict.keys():
        for bam_label in category_dict[category]:
            bam_label = bampath + bam_label + ".bam"
            temp_list += [bam_label]
        temp_string = ",".join(temp_list)
        bam_string += temp_string + " "
        temp_list = []
        temp_string = ""

    for filename in cuffdif_outfiles:
        if os.path.exists(filename):
            print(datetime.now(),"CUFFDIFF OUTPUT",filename,"ALREADY EXISTS")
            print(datetime.now(),"Continuing the pipeline regardless.")


    cmd_cuffdiff = "cuffdiff -p 16 -u -no-update-check {} {}"\
    .format(gtf_reference, bam_string)

    subprocess.check_call(cmd_cuffdiff, shell=True)

    dir_list = os.listdir()
    for filename in cuffdif_outfiles:
        if filename not in dir_list:
            print(datetime.now(), "Cuffdiff output", filename, "is missing")

def diffparser(InFileName):
    OutFileName = "DiffParserOut.txt" #HARDCODED RIGHT NOW!!!
    lines = open(InFileName, "r")
    outputfile = open(OutFileName, "w")
    dataset = ""

    for line in lines:
        if line.startswith("test_id"):
            outputfile.write(line + "\n")
        elif "yes" not in line:
            continue
        else:
            outputfile.write(line + "\n")
    lines.close()
    outputfile.close()

def run_MIA_analysis(gene_exp_file,CRO_input):
    #Filename of file which will contain whole gene_exp data from MIA genes. USED IN FINAL TABLE
    MIA_output1 = 'MIA_expression_data.txt'
    #Filename of file whichw ill contain only gene + significance
    MIA_output2 = 'MIA_significance_description.txt'
    #File that will contain gene name and description. USED IN FINAL TABLE
    MIA_dictionary_file = 'MIA_dictionary.txt'

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
    table_output.close()



########################################################################
####################              MAIN            ######################
########################################################################

if __name__ == '__main__':

    starting_directory = str(os.getcwd())
    start_time = time.time()

    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    set_arguments()
    args = parser.parse_args()

    auto_mode = args.mode
    if auto_mode not in ["manual","overwrite","use_existing"]:
        print("please set mode to 'manual', 'overwrite' or 'use_existing'")
        exit()

    label_dict = select_input_files(args.inputdirectory)
    print_dictionary(label_dict)

    #fastq-split prep
    print(datetime.now(), "\t>Start FASTQ-splitting<")
    os.chdir(args.inputdirectory)

    temp_label_dict = label_dict.copy()
    for label in label_dict:
        if outfile_check(label + "_1.fastq") != 0:
            temp_label_dict.pop(label)

    # FASTQ-SPLIT
    for label in temp_label_dict:
        p = Process(target=fastq_split, args=(label,))
        p.start()
    try:
        p.join()
    except:
        pass

    #check succes of fastq-split
    if check_presence(label_dict,"_1.fastq") != 0:
        call_error("FASTQ-handling","one or more split FASTQ files are missing")
    print(datetime.now(), "\t>Finished FASTQ-splitting<\n")


    #hisat2 prep
    print(datetime.now(), "\t>Start FASTQ Allignment<")
    os.chdir(starting_directory)
    try:
        os.mkdir('hisat2_sam')
    except:
        pass #not propper python way. Should be written to a log file.
    os.chdir('hisat2_sam')

    # HISAT
    run_hisat2(args.hisatindex,label_dict,args.inputdirectory)

    #check succes of hisat
    if check_presence(label_dict,".sam") != 0:
        call_error("FASTQ Allignment","one or more '.sam' files are missing")
    print(datetime.now(), "\t>Finished FASTQ Allignment<\n")

    #samtools-view prep
    print(datetime.now(), "\t>Start SAM-BAM Conversion<")
    temp_label_dict = label_dict.copy()
    for label in label_dict:
        if outfile_check(label + "_unsorted.bam") != 0:
            temp_label_dict.pop(label)

    # SAMTOOLS-view
    for label in temp_label_dict:
        p = Process(target=run_samview, args=(label,))
        p.start()
    try: #IF: temp_list empty
        p.join()
    except:
        pass

    #samtools-sort prep
    temp_label_dict = label_dict.copy()
    for label in label_dict:
        if outfile_check(label + ".bam") != 0:
            temp_label_dict.pop(label)

    # SAMTOOLS-sort
    for label in temp_label_dict:
        p = Process(target=run_samsort, args=(label,))
        p.start()
    try: #IF: temp_list empty
        p.join()
    except:
        pass

    #check succes of samtools
    if check_presence(label_dict,".bam") != 0:
        call_error("SAM-BAM Conversion","one or more '.bam' files are missing")
    print(datetime.now(), "\t>Finished SAM-BAM Conversion<\n")

    #cuffdiff prep
    print(datetime.now(), "\t>Start Cuffdiff<")
    os.chdir(starting_directory)
    try:
        os.mkdir('cuffdiff')
    except:
        pass #not propper python way. Should be written to a log file.
    os.chdir('cuffdiff')

    # CUFFDIFF
    if os.path.exists("gene_exp.diff") == False:
        run_cuffdiff(args.annotation, label_dict, starting_directory + "/hisat2_sam/")
    else:
        print(datetime.now(), "cuffdiff/gene_exp.diff already exists")
    print(datetime.now(), "\t>Finished Cuffdiff<\n")

    os.chdir(starting_directory)

    #output-parse prep
    print(datetime.now(), "\t>Start Output Analysis<")
    MIA_analysis_dir = "output_" + str(datetime.now())[:19].replace(' ','_')
    try:
        os.mkdir(MIA_analysis_dir)
    except:
        pass #not propper python way. Should be written to a log file.
    os.chdir(MIA_analysis_dir)

    # CUFFDIFF-PARSE
    diffparser(starting_directory + "/cuffdiff/gene_exp.diff")
    print(datetime.now(), "Parsing of cuffdiff output completed.")
    print(datetime.now(), "Outputs in DiffParserOut.txt")
    print(datetime.now(), "\t>Finished Output Analysis<\n")


    # PARSE OUTPUTS [MIA SPECIFIC]
    #run_MIA_analysis(starting_directory + "/cuffdiff/gene_exp.diff", starting_directory + '/CRO_names.txt') #COULD turn this into an argument variable
    # FINAL OUTPUT TABLE [MIA SPECIFIC]
    #create_table('MIA_dictionary.txt', 'MIA_expression_data.txt', 'final_table.txt')

    #FINISHING UP
    elapsed_time = time.time() - start_time
    elapsed_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    
    # INFORM USER:
    print(datetime.now(), "\tSTOPPED LOGGING\n")
    print("Congratulations! You have completed your run of the Nijveen pipeline.")
    print("Elapsed time: {}".format(elapsed_time))
    
    print("If you wish to continue your analysis run the provided R scripts with:")
    print("Rscript --default-packages=ggplot2 convert_output.R gene_exp.diff -graph_title.png\n")
    
    print("Thank you! For any feedback please contact the offices of Nijveen group at:")
    print("e-mail: lukas.vandenheuvel@wur.nl")
    print("e-mail: vera.claassen@wur.nl")
    print("e-mail: wing.chung@wur.nl")
    print("e-mail: nathalie.devries@wur.nl")
    print("e-mail: erik.ijland@wur.nl")
    print("")
    



