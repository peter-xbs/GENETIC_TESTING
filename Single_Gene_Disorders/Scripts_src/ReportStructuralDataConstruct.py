# _*_ coding:utf-8 _*_

import os
import json
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='INFO: Make the structural data for the gene test report, including disease info and all sites info!', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-s', dest='site_stat_file', help='Provide the all sites statistic file produced by the DB_sites_stat.py!', required=True)
parser.add_argument('-g', dest='geneomim_config', help='Provide the gene omim configure file produced by the Gene2MimFilter.py!', required=True)
parser.add_argument('-o', dest='output_dir', help='Provide the output_dir!', required=True)
parser.add_argument('-p', dest='prefix', help='Provide the prefix name, usually the sample name!', required=True)
args = parser.parse_args()
site_stat_file = args.site_stat_file
genemim_config = args.geneomim_config
output_dir = args.output_dir
prefix = args.prefix

with open(genemim_config) as f:
    mim2gene_dic = defaultdict(set)
    mim2chpo_dic = {}
    for line in f:
        line_list = line.strip().split('\t')
        if line_list:
            gene, mim, dis_name, chpo, chpo_trunc = line_list[0], line_list[1], line_list[2], line_list[3], line_list[4]
            mim2chpo_dic[mim] = (dis_name, chpo, chpo_trunc)
            mim2gene_dic[mim].add(gene)

with open(site_stat_file) as f:
    gene_dic = defaultdict(list)
    mim2site_num_dic = defaultdict(int)
    for line in f:
        line_dic = {}
        line_list = line.strip().split('\t')
        gene = line_list[3]
        mim = line_list[4]
        site_id = line_list[2]
        mut_type = line_list[11]
        patho_status = line_list[9]
        hetero_status = line_list[-1]
        db = line_list[5]
        mut_expression = line_list[7]+'('+line_list[8]+')'

        if mim in mim2gene_dic:
            if gene in mim2gene_dic[mim]:
                mim2site_num_dic[mim] += 1

        if mut_expression.startswith('-'):
            mut_expression = '--'

        if mut_type.startswith('-'):
            mut_type = '未变异'
        elif 'nonsynonymous' in mut_type:
            mut_type = '非同义突变'
        elif mut_type == 'synonymous SNV':
            mut_type = '同义突变'
        elif mut_type == 'nonframeshift insertion':
            mut_type = '非移码插入'
        elif mut_type == 'frameshift insertion':
            mut_type ='移码插入'
        elif mut_type == 'stopgain':
            mut_type = '无义突变'
        elif mut_type == 'stoploss':
            mut_type = 'stoploss'
        elif mut_type == 'nonframeshift deletion':
            mut_type = '非移码缺失'
        elif mut_type == 'frameshift deletion':
            mut_type = '移码缺失'
        else:
            mut_type = '其它'

        label = ''
        if patho_status.startswith('-'):
            patho_status = '--'
            label = ''
        elif 'ncertain' in patho_status:
            patho_status = '未知'
            label = '空心圆'
        elif 'Likely benign' in patho_status:
            patho_status = '可能良性'
            label = ''
        elif 'Likely pathogenic' in patho_status:
            patho_status = '可能致病'
            label = '实心圆'
        elif 'Pathogenic' in patho_status:
            patho_status = '致病'
            label = '实心圆'
        elif 'Benign' in patho_status:
            patho_status = '良性'
            label = ''

        if hetero_status.startswith('-'):
            hetero_status = '--'
        elif hetero_status == '0.5':
            hetero_status = '杂合'
        else:
            hetero_status = '纯合'

        if db.startswith('-'):
            DB = '--'
        else:
            DB = db+': '+site_id

        line_dic['基因'] = gene
        line_dic['变异位点'] = mut_expression
        line_dic['变异类型'] = mut_type
        line_dic['致病性'] = patho_status
        line_dic['标记'] = label
        line_dic['杂合性'] = hetero_status
        line_dic['数据库收录'] = DB

        gene_dic[gene].append(line_dic)
all_sites_list = []
for gene in gene_dic:
    for item in gene_dic[gene]:
        all_sites_list.append(item)
all_sites_json = os.path.join(output_dir, prefix+'_all_sites.json')
with open(all_sites_json, 'w') as fo:
    json.dump(all_sites_list, fo, ensure_ascii=False)

dis_info_list = []
#first_empty_dic = {}
#first_empty_dic['OMIM'] = ' '
#first_empty_dic['数据库收录位点'] = ' '
#first_empty_dic['疾病名称'] = ' '
#first_empty_dic['疾病临床特征'] = ' '
#first_empty_dic['全长内容'] = ' '
#first_empty_dic['致病基因'] = ' '
#dis_info_list.append(first_empty_dic)

for mim in mim2gene_dic:
    dis_info_dic = {}
    dis_info_dic['OMIM'] = mim
    dis_info_dic['数据库收录位点'] = mim2site_num_dic[mim]
    dis_info_dic['疾病名称'] = mim2chpo_dic[mim][0]
    dis_info_dic['疾病临床特征'] = mim2chpo_dic[mim][2]
    dis_info_dic['全长内容'] = mim2chpo_dic[mim][1]
    dis_info_dic['致病基因'] = '、'.join(list(mim2gene_dic[mim]))
    dis_info_list.append(dis_info_dic)
all_dis_info_json = os.path.join(output_dir, prefix+'_all_dis_info.json')
with open(all_dis_info_json, 'w') as fo:
    json.dump(dis_info_list, fo, ensure_ascii=False)
