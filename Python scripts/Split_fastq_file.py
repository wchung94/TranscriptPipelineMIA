#! /usr/bin/env python3

from sys import argv

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

if __name__ == '__main__':
    input_file = argv[1]
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