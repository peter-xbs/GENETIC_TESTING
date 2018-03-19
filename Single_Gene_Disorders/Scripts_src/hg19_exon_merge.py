# _*_coding:utf-8 _*_

import re
import os
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='INFO: Merge genes exon region, excludes UTR region!', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--ref', dest='hg19_refgene', help='Provide the hg19_refGene.txt', required=True)
parser.add_argument('--output', dest='output', help='Provide the finally merged file')
args = parser.parse_args()

hg19_refgene = args.hg19_refgene
output = args.output

output_dir = os.path.dirname(output)
hg19_exon_merge = os.path.join(output_dir, 'hg19_refGene_intermediate_merged.txt')
hg19_exon_merge_tmp1 = os.path.join(output_dir, 'hg19_refGene_1.tmp')
hg19_exon_merge_tmp2 = os.path.join(output_dir, 'hg19_refGene_2.tmp')
hg19_exon_merge_tmp3 = os.path.join(output_dir, 'hg19_refGene_3.tmp')

with open(hg19_refgene) as f:
    gene_dic = defaultdict(list)
    gene_chro_strand_dic = {}
    for line in f:
        if re.search('chr(\d){1,2}_', line) or re.search('chrU', line):
            continue
        else:
            line_list = line.strip().split('\t')
            gene = line_list[12]
            transcript = line_list[1]
            chrom = line_list[2]
            strand = line_list[3]
            cds_start = line_list[6]
            cds_end = line_list[7]
            if transcript.startswith('NR_'):
                continue
            gene_chro_strand_dic[gene] = (chrom, strand)
            exon_start = [x for x in line_list[9].strip().split(',') if x.strip()]
            exon_end = [x for x in line_list[10].strip().split(',') if x.strip()]
            exon_regions = list(zip(exon_start, exon_end))
            exon_regions = [item for item in exon_regions if item]
            exon_modify_regions = []
            if len(exon_regions) > 1:
                for item in exon_regions:
                    if cds_start >= item[1]:
                        continue
                    elif item[0] <= cds_start <= item[1]:
                        region = (cds_start, item[1])
                        exon_modify_regions.append(region)
                    elif item[0] <= cds_end <= item[1]:
                        region = (item[0], cds_end)
                        exon_modify_regions.append(region)
                        break
                    else:
                        region = (item[0], item[1])
                        exon_modify_regions.append(region)
            else:
                region = [cds_start, cds_end]
                exon_modify_regions.append(region)
            if gene not in gene_dic:
                for region in exon_modify_regions:
                    region = list(region)
                    gene_dic[gene].append(region)
            else:
                for region in exon_modify_regions:
                    for item in gene_dic[gene]:
                        if region[0] >= item[0] and region[1] <= item[1]:
                            break
                        elif (region[0] < item[0]) and (region[1] >= item[0]):
                            gene_dic[gene][gene_dic[gene].index(item)][0] = region[0]
                            if region[1] > item[1]:
                                gene_dic[gene][gene_dic[gene].index(item)][1] = region[1]
                            break
                        elif (item[0] <= region[0] <= item[1]) and (region[1] > item[1]):
                            gene_dic[gene][gene_dic[gene].index(item)][1] = region[1]
                            break
                    else:
                        gene_dic[gene].append(list(region))

with open(hg19_exon_merge, 'w') as fo:
    for gene in gene_dic:
        chrid, strand = gene_chro_strand_dic[gene]
        exon_list = gene_dic[gene]
        for exon in exon_list:
            new_line = '\t'.join([chrid, ] + exon + [gene, strand]) + '\n'
            fo.write(new_line)
cmd1 = "grep -v chrX {}|grep -v chrY|sort -k 1.4n -k 4 > {}".format(hg19_exon_merge, hg19_exon_merge_tmp1)                       
cmd2 = 'grep chrX {}|sort -k1.4n -k 4 > {}'.format(hg19_exon_merge, hg19_exon_merge_tmp2)
cmd3 = 'grep chrY {}|sort -k1.4n -k 4 > {}'.format(hg19_exon_merge, hg19_exon_merge_tmp3)
cmd4 = 'cat {} {} {} > {}'.format(hg19_exon_merge_tmp1, hg19_exon_merge_tmp2, hg19_exon_merge_tmp3, output)
os.system(cmd1)
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)
