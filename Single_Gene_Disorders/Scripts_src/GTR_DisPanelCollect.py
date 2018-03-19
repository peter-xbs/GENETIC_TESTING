# 收集疾病相关基因，通过GTR网站
# 分析所有参数（包括过滤参数均存在于url中，可统一配置到base_url中）
# 所有翻页操作也均在url中，通过解析首页获得全部页面后，构造Url，从而可以避免动态网页爬取
import requests
import re
import os
from bs4 import BeautifulSoup as BS
from lxml import etree
import time
import sys
import json
from collections import defaultdict
from operator import itemgetter

cwd = os.path.dirname(os.path.realpath(__file__))
cwd_p = os.path.dirname(cwd)
config_db_dir = os.path.join(cwd_p, 'Config_files/databases')
gtr2mim = os.path.join(config_db_dir, "GTRConditions2MIM.json")
with open(gtr2mim) as f:
    condition2mim_dic = json.load(f)

base_url = 'https://www.ncbi.nlm.nih.gov/gtr/all/tests/?term={}&filter=testtype:clinical;testpurpose:diagnosis;method:2_7;certification:clia'
if len(sys.argv)>1:
    term = sys.argv[1]
else:
    term = ''
if term == '-h' or term == '--help' or term == '-H':
    print('USAGE: GTR_DisPanelCollect.py "disease name"')
    sys.exit()
elif not term:
    print('USAGE: GTR_DisPanelCollect.py "disease name"')
    sys.exit()
initial_url = base_url.format(term)
headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/62.0.3202.94 Safari/537.36'}
html = requests.get(initial_url, headers=headers)

html_tree = etree.HTML(html.content)

pages_node = html_tree.xpath('//*[@id="gtr_maincontent"]/div/div/div[3]/div/div[1]/h3')
page_text = pages_node[0].text
total_pages = re.findall('of.*?(\d){1,2}', page_text)[0].strip()

# 获得首页表格中的信息
output_name = re.sub('\s+','_',term)+'_panel.txt'
fo = open(output_name, 'w')
dis_genes_dic = defaultdict(int)
dis_omim_dic = defaultdict(int)
omim2term_dic = {}
omimterm_dic = defaultdict(int)

header_line = '\t'.join(['Title', 'Company', 'Country', 'Conditions', 'Gene Num', 'Panel URL','Conditions', 'Genes','Condition2MIM'])+'\n'
fo.write(header_line)

def get_detail_info(url):
    condition_mim_list = []
    html = requests.get(url, headers=headers).content
    tree = etree.HTML(html)
    conditions_whole = tree.xpath('//*[@id="tabs-content"]/div[3]/div[1]/ul/li/a/text()')
    genes = tree.xpath('//*[@id="tabs-content"]/div[3]/div[2]/ul/li/a/text()')
    if conditions_whole:
        condition_str = ';'.join(conditions_whole)
        gene_str = ';'.join(genes)
    else:
        conditions_whole = tree.xpath('//*[@id="tabs-content"]/div[2]/div[1]/ul/li/a/text()')
        genes = tree.xpath('//*[@id="tabs-content"]/div[2]/div[2]/ul/li/a/text()')
        condition_str = ';'.join(conditions_whole)
        gene_str = ';'.join(genes)
    for condition in conditions_whole:
        if condition in condition2mim_dic:
            omim = condition2mim_dic[condition]
            condition_mim_list.append(omim+':'+'"'+condition+'"')
    condition_mim_str = ';'.join(condition_mim_list)
    return [condition_str, gene_str, condition_mim_str]

def save2file(trs, tr_coeff):
    for tr in trs:
        panel_name = tr.select('td p')[0].text
        panel_url = tr.select('td p a')[0]
        panel_url = re.findall('href\=\"(.*?)"', str(panel_url))[0]
        panel_url = 'https://www.ncbi.nlm.nih.gov'+panel_url
        panel_country = tr.select('td div span')[0].text
        panel_comp = tr.select('td div')[0].text
        panel_comp = re.sub(str(panel_country),'',panel_comp)
        panel_conditions = tr.select('td a')[1].text
        panel_genes_num = tr.select('td a')[2].text
        score = 1*tr_coeff
        if int(panel_conditions) <= 50:
            if int(panel_conditions) > 5:
                score = 0.75*tr_coeff
            if int(panel_conditions) > 10:
                score = 0.5*tr_coeff
            elif int(panel_conditions) >30:
                score = 0.25*tr_coeff
            detail_info = get_detail_info(panel_url)
            new_line = '\t'.join([panel_name, panel_comp, panel_country, panel_conditions, panel_genes_num, panel_url, detail_info[0], detail_info[1],detail_info[2]])+'\n'
            gene_list = detail_info[1].split(';')
            for gene in gene_list:
                dis_genes_dic[gene] += score
            conditions_list = detail_info[0].split(';')            
            for condition in conditions_list:
                omimterm_dic[condition] += score
            omim_list = detail_info[2].split(';')
            omim_list = [item.split(':') for item in omim_list]
            for omim in omim_list:
                dis_omim_dic[omim[0]] += score
                if len(omim) > 1:
                    omim2term_dic[omim[0]] = omim[1]
                else:
                    omim2term_dic[omim[0]] = 'NA'
            time.sleep(0.2)
            fo.write(new_line)

if int(total_pages) > 1:
    tr_coeff = 10
    max_page = 8 if int(total_pages)>=7 else int(total_pages)+1
    for i in range(1,int(max_page)):
        new_url = initial_url+'&page={}'.format(i)
        print("Now starting to crawl page: {}".format(i))
        html = requests.get(new_url, headers=headers).content
        soup = BS(html, "html5lib")
        table = soup.find('table')
        trs = table.select('tbody tr')
        save2file(trs, tr_coeff)
        if tr_coeff > 2:
            tr_coeff -= 2
        else:
            tr_coeff = tr_coeff * 0.5
        
    fo.write('\n\n******************************************************\n\n')
    gene_info_list = list(dis_genes_dic.items())
    gene_info_list = sorted(gene_info_list, key=itemgetter(1), reverse=True)
    for gene in gene_info_list:
        new_line = '\t'.join([gene[0], str(gene[1])])+'\n'
        fo.write(new_line)
    fo.write('\n\n******************************************************\n\n')
    omim_info_list = list(dis_omim_dic.items())
    omim_info_list = sorted(omim_info_list, key=itemgetter(1), reverse=True)
    for omim in omim_info_list:
        if omim[0]:
            new_line = '\t'.join([omim[0], omim2term_dic[omim[0]],str(omim[1])])+'\n'
            fo.write(new_line)       
    fo.write('\n\n******************************************************\n\n')
    omim_term_list = list(omimterm_dic.items())
    omim_term_list = sorted(omim_term_list, key=itemgetter(1), reverse=True)
    for omim in omim_term_list:
        new_line = '\t'.join([omim[0], str(omim[1])])+'\n'
        fo.write(new_line)

fo.close()
