# _*_ coding:utf-8 _*_

import re
import os
import json
import argparse
from copy import deepcopy

parser = argparse.ArgumentParser(description='INFO: calculate the db sites and patho/vus sites in the vcf/bam file', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--bam', dest='bam_file', help='Provide the recal.bam file', required=True)
parser.add_argument('--conf', dest='gene_mim_file', help='Provide the gene omim configure file', required=True)
parser.add_argument('--anno', dest='intervar_file', help='Provide the intervar annotated file', required=True)
args = parser.parse_args()
recal_bam_file = args.bam_file
gene_omim_file = args.gene_mim_file
vcf_txt_file = args.intervar_file

# aachange2varid configure file
cwd = os.path.dirname(os.path.realpath(__file__))
cwd_p = os.path.dirname(cwd)
config_dir = cwd_p + '/Config_files'
aachange2varid = os.path.join(config_dir, 'disease_aa_change2var_id.json')
with open(aachange2varid) as f:
    aachange2varid_dic = json.load(f)

# db_file_modify file path
db_file_modify = os.path.join(config_dir, 'IntegratedPathoSitesDatabases_modify.txt')
# recal.bam file path
output_path = os.path.dirname(recal_bam_file)
# site_bed file for samtools
sites_bed = os.path.join(output_path, os.path.basename(recal_bam_file)+'.sites.bed')
sites_bed_fh = open(sites_bed, 'w')
# initial temporary config file
tmp1_config = os.path.join(output_path, os.path.basename(recal_bam_file)+'.tmp1.config')
tmp1_config_fh = open(tmp1_config, 'w')
tmp2_config = os.path.join(output_path, os.path.basename(recal_bam_file)+'.tmp2.config')
tmp2_config_fh = open(tmp2_config, 'w')
tmp3_config = os.path.join(output_path,os.path.basename(recal_bam_file)+'.sites.stat')
tmp3_config_fh = open(tmp3_config, 'w')
tmp1_vcf_txt = os.path.join(output_path, os.path.basename(recal_bam_file)+'.tmp1.vcf_txt')
tmp1_vcf_txt_fh = open(tmp1_vcf_txt, 'w')
# sites depth output file
depth_output = os.path.join(output_path, os.path.basename(recal_bam_file)+'.sites.depth')

# produce intermediate configure files
with open(gene_omim_file, 'r') as f:
    gene_mim_set = set()
    gene_set = set()
    for line in f:
        line_list = re.split(';|:|\t|\s', line.strip())
        gene = line_list[0]
        omim = line_list[1]
        term = gene+':'+omim
        gene_mim_set.add(term)
        gene_set.add(gene)

with open(db_file_modify, 'r') as f:
    header_dic = {}
    header = f.readline()
    header_list = header.strip().split('\t')
    tmp_set = set()
    for item in header_list:
        header_dic[item] = header_list.index(item)
    for line in f:
        line_list = line.strip().split('\t')
        chrid, ID, dis_mim_id, gene, database, site_start, site_end = line_list[header_dic['#chr']], \
                                                                      line_list[header_dic['ID']], \
                                                                      line_list[header_dic['dis_mim_id']], \
                                                                      line_list[header_dic['gene']], \
                                                                      line_list[header_dic['database']], \
                                                                      line_list[header_dic['site_start']], \
                                                                      line_list[header_dic['site_end']]
        gene_mim = gene+':'+dis_mim_id
        if gene_mim in gene_mim_set:
            sites_line = '\t'.join([chrid, site_start, site_end])+'\n'
            if sites_line not in tmp_set:
                tmp_set.add(sites_line)
                sites_bed_fh.write(sites_line)
            new_line = '\t'.join([chrid, site_end, ID, gene, dis_mim_id, database])+'\n'
            if new_line not in tmp_set:
                tmp_set.add(new_line)
                tmp1_config_fh.write(new_line)
sites_bed_fh.close()
tmp1_config_fh.close()
# depth calculation process

depth_calc_cmd = 'samtools depth -b {} {} > {}'.format(sites_bed, recal_bam_file, depth_output)
print(depth_calc_cmd)
#os.system(depth_calc_cmd)

# 加入为数据库位点加入depth数据
with open(depth_output, 'r') as f:
    depth_dic = {}
    for line in f:
        line_list = line.strip().split('\t')
        key = line_list[0]+'_'+line_list[1]
        depth_dic[key] = line_list[2]
with open(tmp1_config) as f:
    for line in f:
        line_list = line.strip().split('\t')
        key = line_list[0]+'_'+line_list[1]
        if key in depth_dic:
            new_line = line.strip()+'\t'+depth_dic[key]+'\n'
            tmp2_config_fh.write(new_line)
        else:
            new_line = line.strip()+'\t'+'0'+'\n'
            tmp2_config_fh.write(new_line)

tmp2_config_fh.close()
# 为数据库位点加入vcf中的信息，同时提取出vcf中VUS级别以上的位点另存为一个文件之中
with open(vcf_txt_file, encoding='utf8', errors='ignore') as f:
    ID_dic = {}
    header_line = f.readline().strip()
    header_dic = {}
    header_list = header_line.split('\t')
    for item in header_list:
        header_dic[item] = header_list.index(item)

    for line in f:
        ID = '.'
        line_list = line.strip().split('\t')
        rsid = line_list[header_dic['avsnp150']]
        gene = line_list[header_dic['Ref.Gene']]
        func = line_list[header_dic['Func.refGene']]
        exonic_func = line_list[header_dic['ExonicFunc.refGene']]
        genotype = line_list[header_dic['Otherinfo']]

        chrid = line_list[0]
        site_start = line_list[1]
        intervar = line_list[13]
        status = re.split(':|PVS1', intervar)[1].strip()
        if gene in gene_set:
            aachange_info = line_list[header_dic['AAChange.refGene']]
            aachange_pos_max = 0
            aachange = '----'
            cchange = '----'
            if aachange_info != '.':
                print(aachange_info)
                aachange_list = re.split(',', aachange_info)
                for elem in aachange_list:
                    elem_list = elem.split(':')
                    AAchange = elem_list[-1]
                    Cchange = elem_list[-2]
                    pos = re.findall('p\.\w+?(\d{1,5})\w+', AAchange)
                    print(pos)
                    if pos:
                        pos = pos[0]
                        if int(pos) > aachange_pos_max:
                            aachange_pos_max = int(pos)
                            cchange = Cchange
                            aachange = AAchange
                            gene_aachange = gene + '.' + aachange.replace('p.', '')
                            if gene_aachange in aachange2varid_dic:
                                ID = aachange2varid_dic[gene_aachange]
            if rsid != '.':
                ID = rsid

            new_line = '\t'.join([chrid, site_start, ID, cchange, aachange, status, gene, func, exonic_func, genotype])+'\n'
            tmp1_vcf_txt_fh.write(new_line)
tmp1_vcf_txt_fh.close()

# vcf_txt和数据库位点取交集
with open(tmp2_config) as f1:
    with open(tmp1_vcf_txt) as f2:
        f2_list = f2.readlines()
        f2_cp_list = deepcopy(f2_list)
        for line1 in f1:
            for line2 in f2_list:
                line1_list = line1.strip().split('\t')
                ID1 = line1_list[2]
                line2_list = line2.strip().split('\t')
                ID2 = line2_list[2]
                if ID1 == ID2:
                    f2_cp_list.remove(line2)
                    new_line = '\t'.join([line1.strip(), line2_list[3], line2_list[4], line2_list[5], line2_list[7], line2_list[8], line2_list[9]])+'\n'
                    break
            else:
                new_line = '\t'.join([line1.strip(), '----', '----', '----', '----', '----', '----'])+'\n'
            tmp3_config_fh.write(new_line)
        for line3 in f2_cp_list:
            line3_list = line3.strip().split('\t')
            status = line3_list[5]
            func = line3_list[7]
            if (func == 'splicing' or func == 'exonic') and (not 'enign' in status):
                new_line = '\t'.join(['chr'+line3_list[0], line3_list[1], line3_list[2], line3_list[6], '----', '----', '----',
                                  line3_list[3], line3_list[4], line3_list[5], line3_list[7], line3_list[8], line3_list[9]])+'\n'
                tmp3_config_fh.write(new_line)
tmp3_config_fh.close()

# 统计并提供报告所用数据

