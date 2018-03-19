# _*_ coding:utf-8 _*_

import os
import re
import json
import argparse

parser = argparse.ArgumentParser(description='INFO: Get the gene-related mim and filter them by key words!', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-g', dest='gene_config', help='Provide the dis-related all gene lists file', required=True)
parser.add_argument('-k', dest='key_words', help='Provide the key words list, all in double quote and seperated by space!', required=True)
parser.add_argument('-o', dest='output', help='Provide the output gene_mim configure file!', required=True)

args = parser.parse_args()
gene_file = args.gene_config
key = args.key_words
output = args.output

key_list = re.split('\s|\t|;|,', key)
cwd = os.path.dirname(os.path.realpath(__file__))
cwd_p = os.path.dirname(cwd)
gene2mim = os.path.join(cwd_p, 'Config_files/databases/gene2mim.json')
omimpheno = os.path.join(cwd_p, 'Config_files/OMIM_PHENO_TEXT.txt')
mim2chpo = os.path.join(cwd_p, 'Config_files/OMIM2CHPO_Integrate_YiHao.txt')

with open(mim2chpo) as f:
    mim2chpo_dic = {}
    mim2disname_dic = {}
    for line in f:
        if line.strip():
            line_list = line.strip().split('\t')
            chpo_terms = re.sub('CHPO:', '', line_list[-1]).strip()
            if not chpo_terms:
                yihao = line_list[-2].strip()
                if yihao != 'Null':
                    chpo_terms = yihao
            if len(chpo_terms) >= 80:
                chpo_terms_trunc = chpo_terms[:75]+'...'
            else:
                chpo_terms_trunc = chpo_terms
            omim = line_list[0]
            mim2chpo_dic[omim] = (chpo_terms, chpo_terms_trunc)
            mim2disname_dic[omim] = line_list[1]
with open(gene2mim) as f:
    gene2mim_dic = json.load(f)

with open(omimpheno, errors='ignore') as f:
    omimpheno_dic = {}
    for line in f:
        line_list = line.strip().split('\t')
        omimpheno_dic[line_list[0]] = line_list[1]

with open(output, 'w') as fo:
    with open(gene_file) as f:
        for line in f:
            line = line.strip()
            if line in gene2mim_dic:
                mim_list = gene2mim_dic[line]
                for mim in mim_list:
                    if mim in omimpheno_dic:
                        pheno_text = omimpheno_dic[mim].lower()
                        for item in key_list:
                            if item.lower() in pheno_text:
                                if mim in mim2chpo_dic:
                                    chpo = mim2chpo_dic[mim]
                                    dis_name = mim2disname_dic[mim]
                                    if not chpo[0].startswith('N'):
                                        new_line = '\t'.join([line, mim, dis_name, chpo[0], chpo[1]])+'\n'
                                        fo.write(new_line)
