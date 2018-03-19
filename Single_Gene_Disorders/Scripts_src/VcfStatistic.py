from collections import defaultdict
import sys

annotated_txt = sys.argv[1]
output = sys.argv[2]
fo = open(output, 'w')
sys.stdout = fo

indel_dic = {'frameshift insertion': 'Frameshift_insertion', 'frameshift deletion': 'Frameshift_deletion', 
             'frameshift substitution': 'Frameshift_substitution', 'nonframeshift deletion': 'Nonframeshift_deletion', 
             'nonframeshift substitution': 'Nonframeshift_substitution', 'unknown': 'unknown', 'nonframeshift insertion': 'Nonframeshift_insertion'}
snp_dic = {'stopgain': 'stopgain_SNV', 'stoploss': 'stoploss_SNV', 'nonsynonymous SNV': 'missense_SNV', 
           'synonymous SNV': 'synonymous_SNV', 'unknown': 'unknown'}
region_dic = {'exonic': 'Exonic', 'exonic;splicing': 'Exonic and splicing', 'splicing': 'splicing', 
              'ncRNA_exonic;splicing': 'ncRNA_exonic and splicing',  'ncRNA_splicing': 'ncRNA_splicing', 
              'ncRNA_exonic': 'ncRNA_exonic', 'ncRNA_UTR5;ncRNA_UTR3': 'ncRNA_UTR5 and ncRNA_UTR3', 
              'ncRNA_UTR5': 'ncRNA_UTR5', 'ncRNA_UTR3': 'ncRNA_UTR3', 'ncRNA_intronic': 'ncRNA_intronic', 
              'UTR5;UTR3': 'UTR5 and UTR3', 'UTR5': 'UTR5', 'UTR3': 'UTR3', 'intronic': 'intronic', 
              'upstream;downstream': 'upstream and downstream', 'upstream': 'upstream', 'downstream': 'downstream', 'intergenic': 'intergenic'}
snp_order = ['stopgain', 'stoploss', 'nonsynonymous SNV', 'synonymous SNV', 'unknown']
indel_order = ['frameshift insertion', 'frameshift deletion', 'frameshift substitution', 'nonframeshift insertion', 'nonframeshift deletion', 'nonframeshift substitution', 'unknown']
region_order = ['ncRNA_UTR5;ncRNA_UTR3', 'ncRNA_UTR5', 'ncRNA_UTR3', 'ncRNA_intronic', 'UTR5;UTR3', 'UTR5', 'UTR3', 'intronic', 'upstream;downstream', 'upstream', 'downstream', 'intergenic', 'exonic;splicing', 'exonic', 'ncRNA_exonic;splicing', 'ncRNA_exonic', 'splicing', 'ncRNA_splicing']
def vcf_stat(file_path):
    with open(file_path) as f:
        header_line = f.readline()
        header_list = header_line.strip().split('\t')
        header_dic = {}
        snp_region_dic = defaultdict(int)
        snp_type_dic = defaultdict(int)
        indel_region_dic = defaultdict(int)
        indel_type_dic = defaultdict(int)
        for header in header_list:
            header_dic[header] = header_list.index(header)
        for line in f:
            line_list = line.strip().split('\t')
            func = line_list[header_dic['Func.refGene']]
            exonic_func = line_list[header_dic['ExonicFunc.refGene']]
            ref = line_list[header_dic['Ref']]
            alt = line_list[header_dic['Alt']]
            if ref in ['A', 'T', 'C', 'G'] and alt in ['A', 'T', 'C', 'G']:
                snp_region_dic[func] += 1
                snp_type_dic[exonic_func] += 1

            else:
                indel_region_dic[func] += 1
                indel_type_dic[exonic_func] += 1

        return snp_region_dic, snp_type_dic, indel_region_dic, indel_type_dic
            
def calc_total(input_dic1, input_dic2, order_list):
    total_num = 0
    for item in order_list:
        print('\t'.join([input_dic2[item],str(input_dic1[item])]))
        total_num += input_dic1[item]
    print('\t'.join(['Total', str(total_num)]))
    
    
results = vcf_stat(annotated_txt)
snp_region_dic, snp_type_dic, indel_region_dic, indel_type_dic = results

print('*********************************************************')
calc_total(snp_region_dic, region_dic, region_order)
print('*********************************************************')
calc_total(snp_type_dic, snp_dic, snp_order)
print('*********************************************************')
calc_total(indel_region_dic, region_dic, region_order)
print('*********************************************************')
calc_total(indel_type_dic, indel_dic, indel_order)
print('*********************************************************')
fo.close()
