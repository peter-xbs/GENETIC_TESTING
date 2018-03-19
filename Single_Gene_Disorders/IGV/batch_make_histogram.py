import os
import sys

file_name = sys.argv[2]
with open(file_name, 'w') as fo:
    line = 'track type="bedGraph" name="{}" color=0,128,0 visiblity=full'.format(file_name) + '\n'
    fo.write(line)
cmd_str = 'python3 ./make_histogram.py {} {} {} {} >> {}'
with open('./histogram_region', 'r') as f:
    for line in f:
        line_list = line.strip().split('\t')
        cmd = cmd_str.format(sys.argv[1], line_list[0].strip(), line_list[1].strip(), line_list[2].strip(), file_name)
        os.system(cmd)
