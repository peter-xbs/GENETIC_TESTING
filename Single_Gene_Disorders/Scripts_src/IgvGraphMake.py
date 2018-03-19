# _*_ coding:utf-8 _*_

import os
import sys
import getopt

def usage():
    print('USAGE INFO: python3 IgvGraphMake.py -B <bam_file> -C <site config file> -H <figure height,also could not provide, default 500)-O <output dir> -h <help>')

opts, args = getopt.getopt(sys.argv[1:], "B:C:H:O:h")
if len(opts) < 2:
    usage()
    sys.exit()

for op, value in opts:
    if op == '-B':
        bam_file = value
    elif op == '-C':
        site_info_file = value
    elif op == '-O':
        output_dir = value
    elif op == '-H':
        height = value
    elif op == '-h':
        usage()
        sys.exit()
    else:
        print('ERROR: Please provide the right arguments!')
        usage()
        sys.exit()


def check_files(bam_file, site_info_file, output_dir):
    if not os.path.isfile(site_info_file):
        print('ERROR: You did not provide the disease site info file! Please check it!')
        sys.exit()
    if not os.path.isdir(output_dir):
        try:
            os.makedirs(output_dir)
            print('Create output directions {}.....'.format(output_dir))
        except Exception as e:
            print('ERROR: The {} does not exist! And you seemed not have the permission to create it! Please check that'.format(output_dir))
            sys.exit()

make_igv_snapshots_tool_path = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[0:-1])+'/IGV'
make_igv_snapshots_tool = os.path.join(make_igv_snapshots_tool_path, 'make_IGV_snapshots.py')

if not os.path.isfile(make_igv_snapshots_tool):
    print("ERROR: Can't find the make_IGV_snapshots.py tool under the {} direction!".format(make_igv_snapshots_tool_path))
    sys.exit()


class IgvGraphMaking:
    """
    该模块输入bam文件和位点信息配置文件，调用主文件夹下的IGV内工具，自动生成截图以及相关pdf附件文档
    """
    def __init__(self, igv_tool, bam_files, site_info_file, output_dir):
        """
        :param igv_tool: 执行命令的工具，来源于github，已完全阅读源码，鉴定可用
        :param bam_file: 单个个体的bam文件
        :param site_info_file: 当个个体的致病位点信息，包括4列，分别是chr, pos_start, pos_end, id_name
        """
        self.igv_tool = igv_tool
        self.bam_files = bam_files.split('+')
        self.site_info_file = site_info_file
        self.output_dir = output_dir

    def make_graph(self):
        if height:
            cmd = 'python3 {} -r {} -ht {} -o {} {}'.format(self.igv_tool, self.site_info_file, height, self.output_dir, ' '.join(self.bam_files))
        else:
            cmd = 'python3 {} -r {} -o {} {}'.format(self.igv_tool, self.site_info_file, self.output_dir, ' '.join(self.bam_files))
        os.system(cmd)

def main():
    check_files(bam_file, site_info_file, output_dir)
    make_graph_instance = IgvGraphMaking(make_igv_snapshots_tool, bam_file, site_info_file, output_dir)
    make_graph_instance.make_graph()
if __name__ == '__main__':
    main()

