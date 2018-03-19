# _*_ coding:utf-8 _*_

import os
import argparse
from collections import defaultdict

cwd = os.path.dirname(os.path.realpath(__file__))
cwd_p = os.path.dirname(cwd)
config_dir = os.path.join(cwd_p, 'Config_files')
exon_config_default = os.path.join(config_dir, 'OMIM_genes_hg19_exons.bed')
db_config_default = os.path.join(config_dir, 'IntegratedPathoSitesDatabases_modify.txt')

parser = argparse.ArgumentParser(description='Description: Calculate the bam file Coverage info for all the disease_related_genes and sites',
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-o', dest='output', help='Provide the output file name',
                    required=True)
parser.add_argument('-dp', dest='depth_file', help='Provide the depth file name',
                    required=True)
parser.add_argument('-e', dest='MIM_gene_exon', help='Provide OMIM genes merged hg19 exons bed file, default: ../Config_files/OMIM_genes_hg19_exons.bed',
                    default=exon_config_default)
parser.add_argument('-db', dest='db_file', help='Provide db_file, default: ../Config_files/IntegratedPathoSitesDatabases_modify.txt',
                    default=db_config_default)
args = parser.parse_args()

db_file = args.db_file
hg19_exon_gene = args.MIM_gene_exon
depth_file = args.depth_file
output = args.output

with open(db_file) as f:
    sites_id_dic = defaultdict(list)
    sites_db_dic = defaultdict(list)
    sites_mim_dic = defaultdict(list)
    db_gene_set = set()
    for line in f:
        if line.startswith('#'):
            continue
        line_list = line.strip().split('\t')
        gene = line_list[7]
        chrid = line_list[0]
        sites_id = line_list[5]
        pos = line_list[-1]
        item = ':'.join([chrid, pos])
        db = line_list[8]
        mim = line_list[6]
        sites_db_dic[item].append(db)
        sites_id_dic[item].append(sites_id)
        sites_mim_dic[item].append(mim)
        db_gene_set.add(gene)

with open(depth_file) as f:
    depth_dic = {}
    for line in f:
        if line.strip():
            line_list = line.strip().split('\t')
            chrid = line_list[0]
            pos = line_list[1]
            dp = line_list[2]
            if chrid not in depth_dic:
                depth_dic[chrid] = {}
            depth_dic[chrid][pos] = dp

with open(hg19_exon_gene) as f:
    exon_gene_dic = defaultdict(list)
    for line in f:
        line_list = line.strip().split('\t')
        gene = line_list[-2]
        start = line_list[1]
        end = line_list[2]
        chrid = line_list[0]
        region = (chrid, start, end)
        if gene in db_gene_set:
            exon_gene_dic[gene].append(region)

with open(output, 'w') as fo:
    header_line = '\t'.join(['#Chr', 'Pos', 'RelativePos', 'Depth', 'ID', 'ExonID','DB_Source', 'MIM', 'Gene'])+'\n'
    fo.write(header_line)
    for gene in exon_gene_dic:
        exon_id = 0
        relative_coordinate = 0
        for region in exon_gene_dic[gene]:
            exon_id += 1
            start = int(region[1])+1
            end = int(region[2])+1
            chrid = region[0]
            for pos in range(start, end):
                relative_coordinate += 1
                key = ':'.join([chrid, str(pos)])
                pos = str(pos)
                if pos in depth_dic[chrid]:
                    dp = depth_dic[chrid][pos]
                else:
                    dp = '0'
                if key in sites_db_dic:
                    db_flag = ';'.join(sites_db_dic[key])
                else:
                    db_flag = '---'
                if key in sites_id_dic:
                    sites_id = ';'.join(sites_id_dic[key])
                else:
                    sites_id = '---'
                if key in sites_mim_dic:
                    mim_id = ';'.join(sites_mim_dic[key])
                else:
                    mim_id = '---'
                new_line = '\t'.join([chrid, pos, str(relative_coordinate), dp, sites_id, str(exon_id), db_flag, mim_id, gene]) + '\n'
                fo.write(new_line)
            relative_coordinate += 100
# print('------------------Sites Not in the Exon Regions!!!----------------------')
# for key in sites_id_dic:
#     chrid = key.split(':')[0]
#     if chrid == 'chrM':
#         continue
#     pos = key.split(':')[1]
#     chr_dic = depth_dic[chrid]
#     if pos in chr_dic:
#         continue
#     else:
#         print(key, sites_id_dic[key])
# print('---------------------------Finished!!!-------------------------------------')

