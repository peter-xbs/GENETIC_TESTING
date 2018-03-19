# _*_ coding:utf-8 _*_

import os
import argparse


parser = argparse.ArgumentParser(description="INFO: calculate all sites in database file's capture efficiency")
parser.add_argument('-b', '--bam', dest='bam_file', help='Provide the recal.bam file', required=True)
parser.add_argument('-t', '--threshold', dest='depth_threshold', help='Provide the depth threshold', required=True)
parser.add_argument('-d', '--depth_file', dest='depth_file', help='Provide the output depth file name', required=True)
args = parser.parse_args()
depth_threshold = args.depth_threshold
bam_file = args.bam_file
depth_file = args.depth_file

cwd = os.path.dirname(os.path.realpath(__file__))
cwd_p = os.path.dirname(cwd)
config_dir = os.path.join(cwd_p, 'Config_files')
sites_bed_file = os.path.join(config_dir, 'all_sites.bed')

cmd = "samtools -b {} {} > {}".format(sites_bed_file, bam_file, depth_file)
# os.system(cmd)

with open(sites_bed_file) as f:
    all_sites_set = set()
    for line in f:
        if line.startswith('#'):
            continue
        else:
            line_list = line.strip().split('\t')
            key = ':'.join([line_list[0], line_list[2]])
            all_sites_set.add(key)

with open(depth_file) as f:
    qualified_set = set()
    for line in f:
        line_list = line.strip().split('\t')
        key = ':'.join([line_list[0],line_list[1]])
        depth = int(line_list[2])
        if depth >= int(depth_threshold):
            if key in all_sites_set:
                qualified_set.add(key)
            else:
                print(line)

percent = round(len(qualified_set)/len(all_sites_set), 3)*100
percent = str(percent)+'%'
print('Capture efficency of {} is {}'.format(depth_file, percent))
