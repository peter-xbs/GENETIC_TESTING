# _*_ coding:utf-8 _*_

import os
import re
import argparse

parser = argparse.ArgumentParser(description='INFO: Screen the qualified gene from GTR websites')
parser.add_argument('-f', '--file', dest='gtr_file', help='Provide the gtr spider file', required=True)
parser.add_argument('-k', '--key_words', dest='key_words', help='Provide the key words, wrapped with double quotes', required=True)

args = parser.parse_args()
gtr_file = args.gtr_file
key_words = args.key_words

key_words = key_words.strip().split(' ')

cwd = os.path.dirname(os.path.realpath(__file__))
cwd_p = os.path.dirname(cwd)
config_dir = os.path.join(cwd_p, 'Config_files')
omim_gene_file = os.path.join(config_dir, 'OMIM_GENE_TEXT.txt')

with open(omim_gene_file, errors='ignore') as f:
    omim_gene_dic = {}
    for line in f:
        if line.strip():
            line_list = line.strip().split('\t')
            if len(line_list) == 2: 
                key = line_list[0]
                omim_gene_dic[key] = line_list[1]
            else:
                print(line)
with open(gtr_file) as f:
    final_gene_list = []
    f_str = f.read()
    f_list = re.split('\*{5,100}', f_str)
    gene_info = f_list[1].strip()
    gene_list = gene_info.split('\n')
    gene_list = [item.split('\t')[0] for item in gene_list]
    for gene in gene_list:
        if gene in omim_gene_dic:
            gene_info_str = omim_gene_dic[gene]
            for key in key_words:
                if key in gene_info_str or key.lower() in gene_info_str:
                    final_gene_list.append(gene)

print('\n'.join(final_gene_list))
