# _*_ coding:utf-8 _*_

import sys

with open('./rs_id_4.avinput', 'r') as f:
    rs_id_dic = {}
    for line in f:
        line_list = line.strip().split('\t')
        key = line_list[0] + '_' + line_list[1]
        rs_id_dic[key] = line_list[-1]

with open(sys.argv[1], 'r') as f:
    for line in f:
        line_list = line.strip().split('\t')
        key = line_list[0] + '_' + line_list[1]
        if key in rs_id_dic:
            print('\t'.join([rs_id_dic[key], line]))
