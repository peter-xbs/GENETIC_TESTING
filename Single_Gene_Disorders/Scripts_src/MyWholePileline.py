# _*_ coding:utf-8 _*_

import json
import re
import copy
import pdfkit
from Haplotype_construct import *
import subprocess


class ExtractSampleId:
    """
    从输入vcf文件中提取出sample_id
    """
    def __init__(self, vcf_path):
        """
        :param vcf_path: 输入的单个病人的vcf文件
        """
        self.vcf_file = vcf_path
    def extract_info(self):
        with open(self.vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    line_list = line.strip().split('\t')
                    sample_id = line_list[-1].strip()
                    break
            else:
                sample_id = 'Null'
                print('Error: No sample_id info in the vcf file, please check it!!')
            return sample_id


class FilterVcf:
    """
    VCF文件在进入流程之前对其按照0.3 <alt/total<0.8, 20<depth<400进行过滤，保证位点的准确性
    """
    def __init__(self, raw_vcf, filtered_vcf):
        """
        :param raw_vcf: 单个人的vcf文件，待过滤
        :param filtered_vcf: 过滤后的输出文件
        """
        self.raw_vcf = raw_vcf
        self.filtered_vcf = filtered_vcf

    def filter_vcf(self):
        with open(self.raw_vcf, 'r') as f_in:
            with open(self.filtered_vcf, 'w') as fo:
                for line in f_in:
                    line = line.strip()
                    if line.startswith('#'):
                        new_line = line + '\n'
                        fo.write(new_line)
                    else:
                        line_list = line.split('\t')
                        info = line_list[9]
                        if info.startswith('0/0') or info.startswith('./.'):
                            pass
                        elif info.startswith('1/1'):
                            DP = info.split(':')[2]
                            if DP == '.':
                                pass
                            elif int(DP) >20:
                                new_line = line + '\n'
                                fo.write(new_line)
                        else:
                            AD = info.split(':')[1]
                            DP = info.split(':')[2]
                            ref = int(AD.split(',')[0])
                            alt = int(AD.split(',')[1])
                            if not (ref + alt) == 0 and not DP == '.':
                                fraction = alt/(ref + alt)
                                if 0.3 < fraction < 0.8 and 40 < int(DP) < 300:
                                    new_line = line + '\n'
                                    fo.write(new_line)


class FilterVcfByDisorderGene:
    """
    将按位点深度过滤后的vcf文件再次进行过滤，按照单基因遗传病对应的基因进行过滤，若某个位点存在于基因区域内，予以保留，否则去除，输入到一个新文件中，供InterVar调用
    """

    def __init__(self, filter_vcf_input_path, data_file_path, filter_vcf_out_path):
        """
        :param filter_vcf_input_path: 按位点深度过滤后的病人vcf文件
        :param data_file_path: 使用的疾病位点信息的配置文件
        :param filter_vcf_out_path: 按致病基因过滤后的病人vcf输出文件
        """
        self.filter_vcf_input_path = filter_vcf_input_path
        self.data_file_path = data_file_path
        self.filter_vcf_out_path = filter_vcf_out_path

    def filter_vcf(self):
        with open(self.filter_vcf_input_path) as f_in:
            with open(self.data_file_path, encoding='utf-8') as f_gene:
                with open(self.filter_vcf_out_path, 'w') as fo:
                    gene_set = set()
                    gene_pat = re.compile('Gene.refGene=(.*?);')
                    for line in f_gene:
                        line_list = line.strip().split('\t')
                        dis_gene = line_list[2].strip()
                        gene_set.add(dis_gene)
                    for line in f_in:
                        if line.startswith('#'):
                            fo.write(line)
                        else:
                            line_gene = re.findall(gene_pat, line)
                            if line_gene:
                                if line_gene[0] != '.':
                                    for gene in line_gene:
                                        if gene in gene_set:
                                            fo.write(line)
                                            break


class VcfToInterVar:
    """
    主要尝试取代PatVcfPreCope这个class的功能，返回一个dict，该dict内部的位点全是经InterVar推断为pathogenic/likely pathogenic状态的点
    """

    def __init__(self, filter_vcf_input, anno_vcf_out, intervar='../InterVar/Intervar.py', output_dict=None):
        """
        :param filter_vcf_input: 输入vcf文件，已经过位点深度以及gene过滤，位点数目相对较少，有利于加快速度
        :param anno_vcf_out: 输出vcf文件, intervar注释后产生文件，内部包含pathogenic信息
        :param intervar: intervar工具的路径(后续如果报错，就调整config.ini文件)
        """
        self.filter_vcf_input = filter_vcf_input
        self.anno_vcf_out = anno_vcf_out
        self.intervar = intervar
        if output_dict is None:
            self.output_dict = {}
        else:
            self.output_dict = output_dict

    def intervar_annotation(self):
        """
        执行intervar，生成判别结果(pathogenic or benign)
        :return: None
        """
        os.environ['intervar'] = self.intervar
        os.environ['input'] = self.filter_vcf_input
        os.environ['output'] = self.anno_vcf_out
        cmd = 'python2.6 $intervar -b hg19 -i $input --input_type=VCF -o $output'
        subprocess.call(cmd, shell=True)

    def contruct_vcf_dict(self):
        """
        主要目标是构建pathogenic/likely pathogenic位点的vcf_dict，目标取代PatVcfPreCope class的输出结果
        :return: vcf_dic for VcfToDisease class使用
        """
        intervar_output_file = self.anno_vcf_out+'.hg19_multianno.txt.intervar'
        intervar_dict = {}
        if os.path.isfile(intervar_output_file):
            with open(intervar_output_file, 'r') as f_intervar:
                for line in f_intervar:
                    line_list = line.split('\t')
                    chr_id = line_list[0]
                    site_pos = line_list[1]
                    intervar_raw_res = line_list[13]
                    intervar_res = re.split(':|PVS1=', intervar_raw_res)[1].strip()
                    coordinate = chr_id+'_'+site_pos
                    if intervar_res == 'Likely pathogenic' or intervar_res == 'Pathogenic':
                        intervar_dict[coordinate] = intervar_res
        else:
            print('There must be something wrong with the intervar steps! Please check it')

        with open(self.filter_vcf_input, 'r') as f:
            aachange_pat = re.compile(r'AAChange\.refGene=(.*?);')
            for line in f:
                if line.startswith('#'):
                    pass
                else:
                    line_list = line.strip().split('\t')
                    line_list = line.split('\t')
                    ref_N = line_list[3]
                    alt_N = line_list[4]
                    site_id = line_list[2]
                    site_genotype = line_list[-1].split(':')[0]
                    site_chr_id = re.sub('chr', '', line_list[0])
                    site_position = line_list[1]
                    coordinate = site_chr_id+'_'+site_position
                    if coordinate in intervar_dict:
                        if site_id.startswith('rs') and site_genotype != '0/0':
                            site_dict = {}
                            site_dict['base_change'] = ref_N + '>' + alt_N
                            site_dict['genotype'] = site_genotype
                            site_dict['site_chr_id'] = site_chr_id
                            site_dict['site_position'] = site_position
                            site_dict['acmg_conclusion'] = intervar_dict[coordinate]
                            site_dict['all_info'] = line
                            self.output_dict[site_id] = site_dict
                        elif site_id == '.' and site_genotype != '0/0':
                            info_str = line_list[7]
                            AAChange_str = re.findall(aachange_pat, info_str)
                            if AAChange_str:
                                AAChange_list = AAChange_str[0].split(',')
                                for elem in AAChange_list:
                                    elem_list = elem.split(':')
                                    gene = elem_list[0]
                                    aachange = elem_list[-1]
                                    site_id = gene + '.' + re.sub('p\.', '', aachange)
                                    if not site_id.startswith('.'):
                                        site_dict = {}
                                        site_dict['base_change'] = ref_N + '>' + alt_N
                                        site_dict['genotype'] = site_genotype
                                        site_dict['site_chr_id'] = site_chr_id
                                        site_dict['site_position'] = site_position
                                        site_dict['acmg_conclusion'] = intervar_dict[coordinate]
                                        site_dict['all_info'] = line
                                        self.output_dict[site_id] = site_dict


class SingleGeneDisorderAPI:
    """
    将单基因遗传病数据结构化后，生成一个嵌套dict用于使用
    """

    def __init__(self, data_file_path, output_dict=None):

        """
        :param data_file_path: 使用的疾病位点信息的配置文件
        :param output_dict: 输出的供使用的嵌套词典
        """
        self.data = data_file_path
        if output_dict is None:
            self.output_dict = {}
        else:
            self.output_dict = output_dict

    def construct_dict(self):
        with open(self.data, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                line_list = [x.strip() for x in line.split('\t')]
                dis_ch_name = line_list[0]
                dis_omim_id = line_list[1]
                dis_gene = line_list[2]
                site_id = line_list[3]
                dis_inheritance = line_list[4]
                site_aa_change = line_list[5]

                site_dict = {}

                site_dict['dis_ch_name'] = dis_ch_name
                site_dict['dis_omim_id'] = dis_omim_id
                site_dict['dis_inheritance'] = dis_inheritance
                site_dict['dis_gene'] = dis_gene
                site_dict['aa_change'] = site_aa_change
                if site_id not in self.output_dict:
                    self.output_dict[site_id] = []
                    self.output_dict[site_id].append(site_dict)
                else:
                    self.output_dict[site_id].append(site_dict)

            return self.output_dict


class PatVcfPreCope:
    """
    目的是将病人vcf文件进行预处理，提取出各潜在致病位点，及其基因型信息，碱基变化信息
    """

    def __init__(self, vcf_file_path, output_dict=None):
        """
        :param vcf_file_path: 病人vcf文件，一个vcf文件中只包含一个病人
        :param output_dict: 将vcf文件中有用信息整合成结构化较好的一个dict
        """
        if output_dict is None:
            self.output_dict = {}
        else:
            self.output_dict = output_dict

        self.data = vcf_file_path

    def construct_dict(self):
        with open(self.data, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    continue
                else:
                    line_list = line.split('\t')
                    ref_N = line_list[3]
                    alt_N = line_list[4]
                    site_id = line_list[2]
                    site_genotype = line_list[-1].split(':')[0]
                    site_chr_id = re.sub('chr', '', line_list[0])
                    site_position = line_list[1]

                    if site_id.startswith('rs') and site_genotype != '0/0':
                        site_dict = {}
                        site_dict['base_change'] = ref_N + '>' + alt_N
                        site_dict['genotype'] = site_genotype
                        site_dict['site_chr_id'] = site_chr_id
                        site_dict['site_position'] = site_position
                        site_dict['all_info'] = line
                        self.output_dict[site_id] = site_dict
                    elif site_id == '.' and site_genotype != '0/0':
                        re_pat = re.compile(r'AAChange\.refGene\=(.*?);')
                        info_str = line_list[7]
                        AAChange_str = re.findall(re_pat, info_str)
                        if AAChange_str:
                            AAChange_list = AAChange_str.split(',')
                            for elem in AAChange_list:
                                elem_list = elem.split(':')
                                gene = elem_list[0]
                                aachange = elem_list[-1]
                                site_id = gene + '.' + re.sub('p\.', '', aachange)
                                if not site_id.startswith('.'):
                                    site_dict = {}
                                    site_dict['base_change'] = ref_N + '>' + alt_N
                                    site_dict['genotype'] = site_genotype
                                    site_dict['site_chr_id'] = re.sub('chr', '', site_chr_id)
                                    site_dict['site_position'] = site_position
                                    site_dict['all_info'] = line
                                    self.output_dict[site_id] = site_dict

            return self.output_dict


class VcfToDisease:
    """
    主要是将vcf中提取出来的位点信息映射到疾病，并输出每个疾病对应的位点详细信息
    """

    def __init__(self, vcf_dic, dis_dic, haplotype_API, sample_id, dis_site_dic=None, disease_dict=None, report_dict=None):
        """
        :param vcf_dic: 从Patient vcf文件中构造的结构化字典
        :param dis_dic: disease和位点关系的配置文件字典
        :param haplotype_API: 传入的haplotype_dict字典，可按照dic['chr_id'][sample_id][position]去查看某个位点的基因型(phase后的)
        :param sample_id: 输入文件的sample ID信息
        :param dis_site_dic: vcf_dic和dis_dic交集字典
        :param disease_dict: 最终转换的以disease为key的字典，表明该病人的哪些位点导致了哪些疾病
        :report_dict: 最终提供给检测报告的数据格式，详情见接口设计文档中的第二部分
        """
        self.vcf_dic = vcf_dic
        self.dis_dic = dis_dic
        self.haplotype_dic = haplotype_API
        self.sample_id = sample_id
        if dis_site_dic is None:
            self.dis_site_dic = {}
        else:
            self.dis_site_dic = dis_site_dic
        if disease_dict is None:
            self.disease_dict = {}
        else:
            self.disease_dict = disease_dict
        if report_dict is None:
            self.report_dict = {}
        else:
            self.report_dict = report_dict

    def integrate_dics(self):
        for key in self.vcf_dic:
            if key in self.dis_dic:
                self.dis_site_dic[key] = []
                for item in self.dis_dic[key]:
                    value_dic = {}
                    value_dic['base_change'] = self.vcf_dic[key]['base_change']
                    value_dic['genotype'] = self.vcf_dic[key]['genotype']
                    value_dic['site_chr_id'] = self.vcf_dic[key]['site_chr_id']
                    value_dic['site_position'] = self.vcf_dic[key]['site_position']
                    value_dic['dis_ch_name'] = item['dis_ch_name']
                    value_dic['dis_omim_id'] = item['dis_omim_id']
                    value_dic['dis_inheritance'] = item['dis_inheritance']
                    value_dic['dis_gene'] = item['dis_gene']
                    value_dic['aa_change'] = item['aa_change']
                    self.dis_site_dic[key].append(value_dic)

        return self.dis_site_dic

    def construct_dis2site_dict(self):
        dis_site_dic = self.integrate_dics()
        for key in dis_site_dic:
            for item in dis_site_dic[key]:
                dis_omim_id = item['dis_omim_id']
                if dis_omim_id not in self.disease_dict:
                    self.disease_dict[dis_omim_id] = set()
                    self.disease_dict[dis_omim_id].add(key)
                else:
                    self.disease_dict[dis_omim_id].add(key)

        return self.disease_dict

    def construct_result_for_report(self):
        dis_result = self.construct_dis2site_dict()
        dis_site_info = self.integrate_dics()
        for key in dis_result:
            key_dict = {}
            key_dict['site_info'] = {}
            if len(dis_result[key]) == 1:
                site_id = list(dis_result[key])[0]
                tmp_list = dis_site_info[site_id]
                # 临时用的一个列表变量，还没想好如何取名
                for elem in tmp_list:
                    if elem['dis_omim_id'] == key:
                        heterozygosity_status = '杂合' if elem ['genotype'] == '0/1' else '纯合'
                        key_dict['dis_omim_id'] = key
                        key_dict['dis_ch_name'] = elem['dis_ch_name']
                        key_dict['dis_inheritance'] = elem['dis_inheritance']
                        key_dict['dis_gene'] = elem['dis_gene']
                        if elem['genotype'] == '0/1':
                            site_genotype = '杂合'
                        else:
                            site_genotype = '纯合'
                        key_dict['site_info'][site_id] = {'genotype': site_genotype,
                                                          'inheritance': elem['dis_inheritance'],
                                                          'base_change': elem['base_change'],
                                                          'aa_change': elem['aa_change'],
                                                          'dis_inheritance': elem['dis_inheritance'],
                                                          'site_position': elem['site_position'],
                                                          'site_chr_id': elem['site_chr_id']}
                        if key_dict['dis_inheritance'] == 'AD':
                            key_dict['haplotype_inheritance_pattern'] = 'AD'
                            key_dict['heterozygosity_status'] = heterozygosity_status
                            key_dict['dis_conclusion'] = '患病'
                        elif key_dict['dis_inheritance'] == 'AR':
                            key_dict['haplotype_inheritance_pattern'] = 'AR'
                            key_dict['heterozygosity_status'] = heterozygosity_status
                            if heterozygosity_status == '杂合':
                                key_dict['dis_conclusion'] = '携带'
                            else:
                                key_dict['dis_conclusion'] = '患病'
                        elif key_dict['dis_inheritance'] == 'XLR':
                            key_dict['haplotype_inheritance_pattern'] = 'XLR'
                            key_dict['heterozygosity_status'] = heterozygosity_status
                            if heterozygosity_status == '杂合':
                                key_dict['dis_conclusion'] = '携带'
                            else:
                                key_dict['dis_conclusion'] = '患病'
                        elif key_dict['dis_inheritance'] == 'XLD':
                            key_dict['haplotype_inheritance_pattern'] = 'XLD'
                            key_dict['heterozygosity_status'] = heterozygosity_status
                            key_dict['dis_conclusion'] = '患病'
                        else:
                            key_dict['haplotype_inheritance_pattern'] = 'YL'
                            key_dict['heterozygosity_status'] = heterozygosity_status
                            key_dict['dis_conclusion'] = '患病'
                        self.report_dict[key] = key_dict
                        break
                else:
                    print('There must be something wrong in this Programe, Please check it at line before 273!!')
                    # 正常情况下，改语句不会输出，上述for循环中的If语句一定可以执行成功，如果不行，一定是系统问题
            else:
                # 该部分脚本待整合的更精炼一些，目前的逻辑有些重复繁琐
                site_id_list = list(dis_result[key])
                dis_sites_genotype_left = []
                dis_sites_genotype_right = []

                for site_id in site_id_list:
                    site_id_chr = dis_site_info[site_id][0]['site_chr_id']
                    site_id_position = dis_site_info[site_id][0]['site_position']
                    site_haplotype = self.haplotype_dic[self.sample_id][site_id_chr][site_id_position].split('|')
                    dis_sites_genotype_left.append(site_haplotype[0])
                    dis_sites_genotype_right.append(site_haplotype[1])

                    tmp_list = dis_site_info[site_id]
                    for elem in tmp_list:
                        if elem['dis_omim_id'] == key:
                            key_dict['dis_omim_id'] = key
                            key_dict['dis_ch_name'] = elem['dis_ch_name']
                            key_dict['dis_inheritance'] = elem['dis_inheritance']
                            key_dict['dis_gene'] = elem['dis_gene']
                            if elem['genotype'] == '0/1':
                                site_genotype = '杂合'
                            else:
                                site_genotype = '纯合'
                            key_dict['site_info'][site_id] = {'genotype': site_genotype,
                                                              'inheritance': elem['dis_inheritance'],
                                                              'base_change': elem['base_change'],
                                                              'aa_change': elem['aa_change'],
                                                              'dis_inheritance': elem['dis_inheritance'],
                                                              'site_position': elem['site_position'],
                                                              'site_chr_id': elem['site_chr_id']}
                            break
                    else:
                        print('There must be something wrong in this Programe, Please check it at line before 341!!')
                if '1' in dis_sites_genotype_left and '1' in dis_sites_genotype_right:
                    heterozygosity_status = '纯合'
                else:
                    heterozygosity_status = '杂合'
                if key_dict['dis_inheritance'] == 'AD':
                    key_dict['haplotype_inheritance_pattern'] = 'AD'
                    key_dict['heterozygosity_status'] = heterozygosity_status
                    key_dict['dis_conclusion'] = '患病'
                elif key_dict['dis_inheritance'] == 'AR':
                    key_dict['haplotype_inheritance_pattern'] = 'AR'
                    key_dict['heterozygosity_status'] = heterozygosity_status
                    if key_dict['heterozygosity_status'] == '杂合':
                        key_dict['dis_conclusion'] = '携带'
                    else:
                        key_dict['dis_conclusion'] = '患病'
                elif key_dict['dis_inheritance'] == 'XLR':
                    key_dict['haplotype_inheritance_pattern'] = 'XLR'
                    key_dict['heterozygosity_status'] = heterozygosity_status
                    if key_dict['heterozygosity_status'] == '杂合':
                        key_dict['dis_conclusion'] = '携带'
                    else:
                        key_dict['dis_conclusion'] = '患病'
                elif key_dict['dis_inheritance'] == 'XLD':
                    key_dict['haplotype_inheritance_pattern'] = 'XLD'
                    key_dict['heterozygosity_status'] = heterozygosity_status
                    key_dict['dis_conclusion'] = '患病'
                else:
                    key_dict['haplotype_inheritance_pattern'] = 'YL'
                    key_dict['heterozygosity_status'] = heterozygosity_status
                    key_dict['dis_conclusion'] = '患病'
                self.report_dict[key] = key_dict
                print(key_dict)

        return self.report_dict


class GetInfoFromSencondAPI:
    """
    该class主要用于承接VCFTODISEASE模块过来的数据，进一步从第二个API中提取所需要的数据，供基因检测报告用
    """

    def __init__(self, vcf_to_dis_dict, second_API, vcf_to_report_dict=None):
        """
        :param vcf_to_dis_dict: VcfToDisease Class中第三个函数construct_result_for_report返回的结果，即模块接口设计v2中第二步的内容
        :param vcf_to_report_dict: 此步输出的report_dict最终直接用于html文档的构建
        """
        self.vcf_to_dis_dict = vcf_to_dis_dict
        self.second_API = second_API
        if vcf_to_report_dict is None:
            self.vcf_to_report_dict = {}
        else:
            self.vcf_to_report_dict = vcf_to_report_dict

    def add_new_info(self):
        with open(self.second_API, 'r') as f:
            second_API_dict = json.load(f)
            reference_list = []
            self.vcf_to_report_dict['diseases'] = {}
            for key in self.vcf_to_dis_dict:
                self.vcf_to_report_dict['diseases'][key] = copy.deepcopy(self.vcf_to_dis_dict[key])
                self.vcf_to_report_dict['diseases'][key]['dis_description'] = second_API_dict['OMIM_API'][key]['description']
                self.vcf_to_report_dict['diseases'][key]['dis_en_name'] = second_API_dict['OMIM_API'][key]['en_name']
                self.vcf_to_report_dict['diseases'][key]['dis_type'] = second_API_dict['OMIM_API'][key]['dis_type']
                self.vcf_to_report_dict['diseases'][key]['treat_recomm'] = second_API_dict['OMIM_API'][key]['treat_recomm']
                rs_IDs = list(self.vcf_to_dis_dict[key]['site_info'].keys())
                for rs_ID in rs_IDs:
                    if rs_ID in second_API_dict['rsID_ref_API']:
                        reference_list.extend(second_API_dict['rsID_ref_API'][rs_ID])
            self.vcf_to_report_dict['reference_list'] = reference_list

            return self.vcf_to_report_dict


class GetBasicInfoOfPatient:
    """
    用于获取受检者基本信息，可填充至检测报告首部
    """

    def __init__(self, basic_info_file_path, output_info_dict=None):
        self.basic_info_data = basic_info_file_path
        if output_info_dict is None:
            self.output_info_dict = {}
        else:
            self.output_info_dict = output_info_dict

    def get_basic_info(self):
        with open(self.basic_info_data, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                line_list = line.split('\t')
                if line.startswith('#'):
                    pass
                else:
                    name, sex, birthday, tel, initial_diagnose, father_dis, mother_dis, familiy_dis_inheritance, sample_id, sample_type, \
                    sample_get_date, sample_receive_date, sample_report_date, sample_department, sample_doctor, sample_seq_com = \
                    line_list[0], line_list[1], line_list[2], line_list[3], line_list[4], line_list[5], line_list[6], line_list[7], line_list[8],\
                    line_list[9], line_list[10], line_list[11], line_list[12], line_list[13], line_list[14], line_list[15]

                    self.output_info_dict[sample_id] = {}

                    self.output_info_dict[sample_id]['sex'] = sex
                    self.output_info_dict[sample_id]['name'] = name
                    self.output_info_dict[sample_id]['birthday'] = birthday
                    self.output_info_dict[sample_id]['tel'] = tel
                    self.output_info_dict[sample_id]['initial_diagnose'] = initial_diagnose
                    self.output_info_dict[sample_id]['father_dis'] = father_dis
                    self.output_info_dict[sample_id]['mother_dis'] = mother_dis
                    self.output_info_dict[sample_id]['familiy_dis_inheritance'] = familiy_dis_inheritance
                    self.output_info_dict[sample_id]['sample_type'] = sample_type
                    self.output_info_dict[sample_id]['sample_get_date'] = sample_get_date
                    self.output_info_dict[sample_id]['sample_receive_date'] = sample_receive_date
                    self.output_info_dict[sample_id]['sample_report_date'] = sample_report_date
                    self.output_info_dict[sample_id]['sample_department'] = sample_department
                    self.output_info_dict[sample_id]['sample_doctor'] = sample_doctor
                    self.output_info_dict[sample_id]['sample_seq_com'] = sample_seq_com

            return self.output_info_dict
    def get_pat_parent_dis_pheno(self, sample_id):
        basic_info_dict = self.get_basic_info()
        return basic_info_dict[sample_id]['father_dis'], basic_info_dict[sample_id]['mother_dis']


class BasicInfoHtml:
    """
    基因检测报告头部信息，即基本信息构建，接受一个基本信息的字典
    """
    def __init__(self, sample_id, basic_info_dict):
        self.sample_id = sample_id
        self.basic_info_dict = basic_info_dict

    def construct_html(self):
        header_line = '<tr style="background-color:#6699FF; border:2px"><td>受检者信息</td><td>报告信息</td><td>送检信息</td></tr>'
        second_line_style = '<tr style="font-size:12px;background-color:#CCC;"><td>姓名：{}</td><td>样本编号：{}</td><td>送检单位：{}</td></tr>'
        second_line = second_line_style.format(self.basic_info_dict['name'], sample_id, self.basic_info_dict['sample_department'])
        third_line_style = '<tr style="font-size:12px;background-color:#FFF;"><td>性别：{}</td><td>样本类型：{}</td><td>送检医生：{}</td></tr>'
        third_line = third_line_style.format(self.basic_info_dict['sex'], self.basic_info_dict['sample_type'], self.basic_info_dict['sample_doctor'])
        fourth_line_style = '<tr style="font-size:12px;background-color:#CCC;"><td>出生日期：{}</td><td>样本采集日期：{}</td><td>检测单位：{}</td></tr>'
        fourth_line = fourth_line_style.format(self.basic_info_dict['birthday'], self.basic_info_dict['sample_get_date'], self.basic_info_dict['sample_department'])
        fifth_line_style = '<tr style="font-size:12px;background-color:#FFF;"><td>联系电话：{}</td><td>样本接受日期：{}</td><td></td></tr>'
        fifth_line = fifth_line_style.format(self.basic_info_dict['tel'], self.basic_info_dict['sample_receive_date'])
        sixth_line_style =  '<tr style="font-size:12px;background-color:#CCC;"><td>性别：{}</td><td></td></tr>'
        six_line = sixth_line_style.format(self.basic_info_dict['sex'])
        seventh_line_style = '<tr style="font-size:12px;background-color:#FFF;"><td>父母疾病表型：{}</td></tr>'
        seventh_line = seventh_line_style.format(self.basic_info_dict['father_dis']+';'+self.basic_info_dict['mother_dis'])
        content = header_line + second_line + third_line + fourth_line + fifth_line + six_line + seventh_line
        table_frame = '<table style="border:2px solid #6699FF;rules:rows; cellspacing:0.5px">{}</table>'.format(content)
        return table_frame


class DisBackGroundHtml:
    """
    用来构建第二部分，遗传背景介绍以及基因检测意义
    """
    def __init__(self, background_text_path):
        """
        :param background_text_path: 为背景资料准备的文本文件路径
        """
        self.background_text_path = background_text_path

    def construct_html(self):
        file = open(self.background_text_path, 'r', encoding='utf-8')
        f = [x.strip() for x in file.readlines()]
        title = '遗传背景知识和基因检测意义'
        line_style = '<p text-indent:2em; line-height:15px>{}</p>'
        line1 = line_style.format(f[0])
        line2 = line_style.format(f[1])
        line3 = line_style.format(f[2])
        text = ''.join([line1, line2, line3])
        file.close()
        return text


class TestResultSummaryHtml:
    """
    用来构建第三部分，展示个人单基因遗传病检测结果的汇总情况
    """
    def __init__(self, text_file_path, report_dict, pat_parent_pheno_status):
        """
        :param text_file_path: 第三部分，表格前面的总结性信息，提供描述模板
        :param report_dict: 受检者所有疾病信息汇总的字典
        :param pat_parent_pheno_status: 受检者父母疾病表型，来源于基本信息配置文件，可由class GetBasicInfoOfPatient中的get_pat_parent_dis_pheno函数获得
        """

        self.text_file_path = text_file_path
        self.report_dict = report_dict
        self.pat_parent_pheno_status = pat_parent_pheno_status

    def construct_html(self):
        dis_nums = len(self.report_dict['diseases'])
        with open(self.text_file_path, 'r', encoding='utf-8') as f1:
            pre_text = f1.read()
            text = pre_text.format(dis_nums)
        text_line =  '<p text-indent:2em; line-height:15px><strong>{}</strong></p>'.format(text)
        line_style = '<tr style="background-color:#3399FF; border:1px; height:20px"><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td width=55>{}</td><td>{}</td><td>{}</td></tr>'
        dis_line_style = '<tr style="font-size:13px;text-height:2px;"><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td style="text-align:center;">{}</td><td>' \
                         '<b style="color:red">{}</b></td></tr>'
        carrier_line_style = '<tr style="font-size:13px;text-height:2px"><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td style="text-align:center;">{}</td><td>' \
                             '<b style="color:blue">{}</b></td></tr>'
        dis_line_list = []
        carrier_line_list = []
        i = 0
        j = 0
        header_line = line_style.format('遗传病名称', 'OMIM ID', '致病基因名称', '遗传方式', '类型', '受检者父母疾病表型', '检测结论')
        for key in self.report_dict['diseases']:
            if self.pat_parent_pheno_status[0] == self.report_dict['diseases'][key]['dis_ch_name']:
                father_pheno = '(+)'
            else:
                father_pheno = '(-)'
            if self.pat_parent_pheno_status[1] == self.report_dict['diseases'][key]['dis_ch_name']:
                mother_pheno = '(+)'
            else:
                mother_pheno = '(-)'
            # 上述代码有修改余地，暂未考虑父母患有多种疾病的情况
            if father_pheno == '(+)' and mother_pheno == '(+)':
                parent_pheno = '<font style="font-size:8px;color:red;">父亲：'+father_pheno + \
                               '；母亲：'+mother_pheno+'</font>'
            elif father_pheno == '(-)' and mother_pheno == '(+)':
                parent_pheno = '<font style="font-size:8px;color:blue;">父亲：'+father_pheno + \
                               '；母亲：'+mother_pheno+'</font>'
            elif father_pheno == '(+)' and mother_pheno == '(-)':
                parent_pheno = '<font style="font-size:8px; color:red">父亲：'+father_pheno + \
                               '；母亲：'+mother_pheno+'</font>'
            else:
                parent_pheno = '<font style="font-size:8px;color:blue;">父亲：'+father_pheno + \
                               '；母亲：'+mother_pheno+'</font>'

            if self.report_dict['diseases'][key]['dis_conclusion'] == '患病':
                i += 1
                line = dis_line_style.format(self.report_dict['diseases'][key]['dis_ch_name'],
                                             self.report_dict['diseases'][key]['dis_omim_id'],
                                             self.report_dict['diseases'][key]['dis_gene'],
                                             self.report_dict['diseases'][key]['haplotype_inheritance_pattern'],
                                             self.report_dict['diseases'][key]['dis_type'],
                                             parent_pheno,
                                             self.report_dict['diseases'][key]['dis_conclusion'])
                dis_line_list.append(line)
            else:
                j += 1
                line = carrier_line_style.format(self.report_dict['diseases'][key]['dis_ch_name'],
                                                 self.report_dict['diseases'][key]['dis_omim_id'],
                                                 self.report_dict['diseases'][key]['dis_gene'],
                                                 self.report_dict['diseases'][key]['haplotype_inheritance_pattern'],
                                                 self.report_dict['diseases'][key]['dis_type'],
                                                 parent_pheno,
                                                 self.report_dict['diseases'][key]['dis_conclusion'])
                carrier_line_list.append(line)

        dis_remained_line = ''.join(dis_line_list)
        dis_table_frame = '<table style="border:2px solid #6699FF;rules:rows; cellspacing:0.5px">' + header_line + dis_remained_line + '</table>'
        carrier_remained_line = ''.join(carrier_line_list)
        carrier_table_frame = '<table style="border:2px solid #6699FF;rules:rows; cellspacing:0.5px">' + header_line + carrier_remained_line + '</table>'
        if i == 0 and j == 0:
            return '<p><strong>恭喜您，基因检测未显示您患有任何单基因遗传病</strong></p>'
        elif i == 0 and j != 0:
            para = '<p><strong>恭喜您，基因检测未显示您患有任何单基因遗传病，但您是单基因遗传病致病基因携带者</strong></p>'
            return para + carrier_table_frame
        elif i != 0 and j == 0:
            return text_line + dis_table_frame
        else:
            para = '<p><strong>检测您携带如下疾病的相关致病基因，结果如下表：</strong></p>'
            return text_line + dis_table_frame + para + carrier_table_frame


class TestResultDetailHtml:
    """
    此模块主要用于展示每种疾病及其详细位点信息
    """
    def __init__(self, report_dict):
        """
        :param report_dict:接收详细report_dict字典中的信息
        """
        self.report_dict = report_dict

    def construct_site_map_to_figure(self):
        site_map_to_figure_dict = {}
        dis_list = []
        carrier_list = []
        for key in self.report_dict['diseases']:
            site_info = self.report_dict['diseases'][key]['site_info']
            if self.report_dict['diseases'][key]['dis_conclusion'] == '患病':
                for item in site_info:
                    dis_list.append(item)
            else:
                for item in site_info:
                    carrier_list.append(item)
        total_list = dis_list + carrier_list
        i = 0
        for elem in total_list:
            i += 1
            site_map_to_figure_dict[elem] = '附图'+str(i)
        return site_map_to_figure_dict

    def construct_html(self):
        site_map_to_figure_dict = self.construct_site_map_to_figure()
        dis_nums = len(self.report_dict['diseases'])
        background_para = '<p text-indent:2em><strong>通过检测变异位点遗传信息，获得了多个致病突变变异位点，发现受检者{}， {}详细变异位点信息如下表：</strong></p>'
        sub_header_line_style = '<td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td>'
        sub_line_style = '<td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td style="color:#0099FF;font-size:9px;">{}</td>'
        table_line_style = '<tr height:10px><td rowspan="{}">{}</td><td rowspan="{}">{}</td>{}</tr>'
        header_line_style = '<tr style="background-color:#3399FF; border:1px; height:20px"><td>{}</td><td>{}</td>{}</tr>'
        sub_header_line = sub_header_line_style.format('核酸改变', '氨基酸改变', 'rs ID', '杂合/纯合', '遗传方式', '变异位点IGV图')
        header_line = header_line_style.format('遗传病名称', '致病基因名称', sub_header_line)
        dis_line_list = [header_line]
        carrier_line_list = [header_line]
        dis = 0
        carrier = 0
        for key in self.report_dict['diseases']:
            site_info = self.report_dict['diseases'][key]['site_info']
            site_num = len(site_info)
            if self.report_dict['diseases'][key]['dis_conclusion'] == '患病':
                dis += 1
                sub_line_list = []
                site_td_num = 0
                # 该变量用于控制当一种疾病的多个位点均发生变异的情况
                for item in site_info:
                    site_td_num += 1
                    sub_line = sub_line_style.format(
                        self.report_dict['diseases'][key]['site_info'][item]['base_change'],
                        self.report_dict['diseases'][key]['site_info'][item]['aa_change'],
                        item,
                        self.report_dict['diseases'][key]['site_info'][item]['genotype'],
                        self.report_dict['diseases'][key]['site_info'][item]['dis_inheritance'],
                        site_map_to_figure_dict[item])
                    if site_td_num == 1:
                        sub_line_list.append(sub_line)
                    else:
                        sub_line_list.append('<tr>'+sub_line+'</tr>')

                line = table_line_style.format(site_num,
                                               self.report_dict['diseases'][key]['dis_ch_name'],
                                               site_num,
                                               self.report_dict['diseases'][key]['dis_gene'],
                                               ''.join(sub_line_list))
                dis_line_list.append(line)
            else:
                carrier += 1
                sub_line_list = []
                site_td_num = 0
                for item in site_info:
                    site_td_num += 1
                    sub_line = sub_line_style.format(
                        self.report_dict['diseases'][key]['site_info'][item]['base_change'],
                        self.report_dict['diseases'][key]['site_info'][item]['aa_change'],
                        item,
                        self.report_dict['diseases'][key]['site_info'][item]['genotype'],
                        self.report_dict['diseases'][key]['site_info'][item]['dis_inheritance'],
                        site_map_to_figure_dict[item])
                    if site_td_num == 1:
                        sub_line_list.append(sub_line)
                    else:
                        sub_line_list.append('<tr>' + sub_line + '</tr>')
                line = table_line_style.format(site_num,
                                               self.report_dict['diseases'][key]['dis_ch_name'],
                                               site_num,
                                               self.report_dict['diseases'][key]['dis_gene'],
                                               ''.join(sub_line_list))
                carrier_line_list.append(line)
        table1 = background_para.format('患有'+str(dis)+'种单基因遗传病', '所患疾病') + '<table style="border:2px solid #6699FF;rules:rows; cellspacing:0.5px">'+''.join(dis_line_list)+'</table>'
        table2 = background_para.format('携带'+str(carrier)+'种单基因遗传病的相关致病位点', '相关疾病') + '<table style="border:2px solid #6699FF;rules:rows; cellspacing:0.5px">'+''.join(carrier_line_list)+'</table>'
        if dis == 0 and carrier == 0:
            return '<p style="color:green">恭喜您，经基因组测序检测，未发现致病相关变异位点!</p>'
        elif dis != 0 and carrier == 0:
            return table1
        elif dis == 0 and carrier != 0:
            return table2
        else:
            return table1 + table2


class DiseaseDescriptionHtml:
    """
    该模块主要用来介绍受检者相关的各种疾病的介绍(必须是患病的介绍)
    """
    def __init__(self, report_dict):
        """
        :param report_dict: 从中抽取各疾病相关的介绍
        """
        self.report_dict = report_dict
    def construct_html(self):
        line_style = '<p style="text-indent:0em;color:#3399FF;"><strong>{}.{}</strong></p><p><font face="Times New Roman">{}</font></p>'
        line_list = []
        i = 0
        for key in self.report_dict['diseases']:
            if self.report_dict['diseases'][key]['dis_conclusion'] == '患病':
                i += 1
                line = line_style.format(str(i), self.report_dict['diseases'][key]['dis_ch_name'], self.report_dict['diseases'][key]['dis_description'])
                line_list.append(line)
        return ''.join(line_list)


class TreatRecommendationHtml:
    """
    此模块提供受检者的治疗建议，若是患病，提供详细的治疗方案；若是携带者，提醒优生优育
    """
    def __init__(self, report_dict):
        self.report_dict = report_dict

    def construct_html(self):
        line_style = '<p style="text-indent:0em;color:#3399FF;"><strong>{}.  {}<b style="color:{}">{}</b>建议：</strong></p><p>{}</p>'
        dis_line_list = []
        carrier_line_list = []
        i = 0
        j = 0
        for key in self.report_dict['diseases']:
            if self.report_dict['diseases'][key]['dis_conclusion'] == '患病':
                i += 1
                line = line_style.format(str(i), self.report_dict['diseases'][key]['dis_ch_name'], 'red', '患病治疗', self.report_dict['diseases'][key]['treat_recomm'])
                dis_line_list.append(line)
            else:
                j += 1
                line = line_style.format(str(j), self.report_dict['diseases'][key]['dis_ch_name'], 'blue', '致病基因携带', self.report_dict['diseases'][key]['treat_recomm'])
                carrier_line_list.append(line)
        total_list = dis_line_list + carrier_line_list

        return ''.join(total_list)


class ReferenceHtml:
    """
    该模块用于获取选取位点的参考文献信息，可以精准匹配
    """
    def __init__(self, report_dict):
        self.report_dict = report_dict

    def construct_html(self):
        line_style = '<p height=5px, style="font-size:10px"><font face="Times New Roman">{}. {}<strong>{}</strong>{}<font style="color:blue">{}</font></font></p>'
        references = self.report_dict['reference_list']
        i = 0
        line_list = []
        for elem in references:
            i += 1
            pmid = elem[3] if len(elem) == 4 else ''
            pmid = re.sub('related citations|images|,', '', pmid).strip()
            reference_line = line_style.format(str(i), elem[0], elem[1], elem[2], pmid)
            line_list.append(reference_line)

        return ''.join(line_list)


class WholeDocHtml:
    """
    该模块可将以上所有模块进行汇总，生成检测报告全文文档
    """
    def __init__(self, css_config, basic_info, background_info, summary_info, detail_info, specification_info, recommendation_info, reference_info):
        """
        :css_config: 用于提供部分样式格式
        :param basic_info: 检测报告基本信息部分
        :param background_info: 检测报告遗传背景知识和基因检测意义
        :param summary_info: 检测报告单基因遗传病检测结果概览
        :param detail_info: 检测报告致病变异位点检测结果
        :param specification_info: 检测报告单基因遗传病检测详细说明
        :param recommendation_info: 检测报告给予受测者健康建议
        :param reference_info: 检测报告参考文献部分
        """
        self.basic_info = basic_info
        self.background_info = background_info
        self.summary_info = summary_info
        self.detail_info = detail_info
        self.specification_info = specification_info
        self.recommendation_info = recommendation_info
        self.reference_info = reference_info
        self.css_config= css_config

    def construct_html(self, input_file, output_file):
        with open(output_file, 'w', encoding='utf-8') as fo:
            with open(input_file, 'r', encoding='utf-8') as f:
                doc_string_blank = f.read()
                doc_string_fill = doc_string_blank.format(self.css_config[0], self.css_config[1],
                                                          self.css_config[2], self.css_config[3],
                                                          self.css_config[4], self.css_config[5],
                                                          self.basic_info, self.background_info,
                                                          self.summary_info, self.detail_info,
                                                          self.specification_info, self.recommendation_info,
                                                          self.reference_info)
                fo.write(doc_string_fill)

# ###################具体流程###########################

# 1. 提取出单个病人的sample_id，由VCF文件中提取
sample_id = ExtractSampleId('./0116D01Q200088A01_2101.vcf').extract_info()

# 2. 构造或加载haplotype API
## 生成Haplotype API字典
# input_vcf = '/opt/mnt/biogfs2/wanghengtao/Fuwai/Cardiovascular/WGS/8-NORMAL_from_lung_kdiney_liver/GVCF_list.txt.combine.vcf'
# output_path = '../'
# config_file = './tmp'
# haplotype = ConstructHaplotypeRef(input_vcf, output_path, config_file)
# haplotype_dict = test.construct_position_to_genotypes_dict()
with open('./Config_files/haplotype.json', 'r') as f:
    haplotype_dict = json.load(f)
print(haplotype_dict['X'].keys())

dic1_instance = SingleGeneDisorderAPI('基因检测报告v0.1_5.txt')
Dic1 = dic1_instance.construct_dict()
filter = FilterVcf('./0116D01Q200088A01_2101.vcf', './基因检测报告v0.1_single_pat_filter.vcf').filter_vcf()
dic2_instance = PatVcfPreCope('././基因检测报告v0.1_single_pat_filter.vcf')
Dic2 = dic2_instance.construct_dict()
test = VcfToDisease(Dic2, Dic1, haplotype_dict, sample_id)
dis_site_info = test.integrate_dics()
dis_result = test.construct_dis2site_dict()
Dic3 = test.construct_result_for_report()
test2 = GetInfoFromSencondAPI(Dic3, './基因检测报告v0.1_Second_API.json')
test3 = GetBasicInfoOfPatient('./基因检测报告v0.1_基本信息配置文件')
report_dict = test2.add_new_info()
for key in report_dict['diseases']:
    print(key, report_dict['diseases'][key])
# for key in dis_site_info:
#     print(key, dis_site_info[key])


sample_basic_dict = test3.get_basic_info()
sample_id = '0116D01D600015C_2101'
single_pat_basic_info_dict = sample_basic_dict[sample_id]
single_pat_parent_pheno = test3.get_pat_parent_dis_pheno(sample_id)


basic_html = BasicInfoHtml(sample_id, single_pat_basic_info_dict).construct_html()
background_html = DisBackGroundHtml('./基因检测报告v0.1_遗传背景介绍文本').construct_html()
summary_html = TestResultSummaryHtml('./基因检测报告v0.1_疾病汇总部分简介文本', report_dict, single_pat_parent_pheno).construct_html()
description_html = DiseaseDescriptionHtml(report_dict).construct_html()
recomm_html = TreatRecommendationHtml(report_dict).construct_html()
reference_html = ReferenceHtml(report_dict).construct_html()
test_result_detail_html = TestResultDetailHtml(report_dict).construct_html()
css_config = ['{color:#0099FF; font-size:20px}', '{height:20px; text-align:left; border:blue; line-type:gt;text-height:2px; font-size:12px}',
              '{background: #CCC; font-size:12px}', '{background: #FFF; font-size:12px}', '{line-height:15px; font-size:13px; text-indent:2em}',
              '{border: 3px solid orange;border-top-left-radius: 50px;border-top-right-radius: 50px;border-bottom-left-radius: 50px;border-bottom-right-radius: 50px;padding:10px;}']
doc = WholeDocHtml(css_config, basic_html, background_html, summary_html, test_result_detail_html, description_html, recomm_html, reference_html)
doc.construct_html('./whole_doc_template.html', './Single_Gene_Disorder_Test_Report.html')

config = pdfkit.configuration(wkhtmltopdf='C:/Program Files/wkhtmltopdf/bin/wkhtmltopdf.exe')
with open('./Single_Gene_Disorder_Test_Report.html', 'r', encoding='utf-8') as f:
    pdfkit.from_file(f, './Single_Gene_Disorder_Test_Report.pdf', configuration=config)











