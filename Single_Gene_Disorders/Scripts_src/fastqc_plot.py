# _*_ coding:utf-8 _*_

import os
import re
import sys
import argparse


def get_info(files, pos=0, pos2=0):
    with open(files) as f:
        fq_1 = f.read().split('>>')
        for item in fq_1:
            n_pos_content_list = []
            if "base N content" in item:
                elems = item.split('#')[1].split('\n')[1:-1]
                for s in elems:
                    n_pos_content_list.append(s.split('\t')[1])
                break
        else:
            print('Error: N-content info are missing in the {}'.format(fq_1_summary))
            sys.exit()
        for item in fq_1:
            if "sequence content" in item:
                elems = item.split('#')[1].split('\n')
                header_dic = {}
                for elem in elems:
                    if elem.startswith('B'):
                        elem_list = elem.split('\t')
                        for s in elem_list:
                            header_dic[s] = elem_list.index(s)
                    else:
                        if elem.strip():
                            elem_list = elem.split('\t')
                            new_line = '\t'.join(['A', 'NULL', elem_list[header_dic['A']], 'T', 'NULL',
                                                  elem_list[header_dic['T']], 'G', 'NULL', elem_list[header_dic['G']],
                                                  'C', 'NULL', elem_list[header_dic['C']], 'N', 'NULL', n_pos_content_list[pos]])
                            pos += 1
                            if pos < 9 or pos>36:
                                fo.write(new_line+'\n')
                            else:
                                for k in range(5):
                                    fo.write(new_line+'\n')
        for item in fq_1:
            if "base sequence quality" in item:
                elems = item.split('\n')[2:-1]
                for elem in elems:
                    elem_list = elem.split('\t')
                    new_line = elem_list[1]+'\n'
                    pos2 += 1
                    if pos2 < 9 or pos2 >36:
                        fo2.write(new_line)
                    else:
                        for k in range(5):
                            fo2.write(new_line)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Extract QC information from fastqc results', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--input1', dest='fq_zip_1', help='Provide the fq1.zip', required=True)
    parser.add_argument('--input2', dest='fq_zip_2', help='Provide the fq2.zip', required=True)
    parser.add_argument('--output', dest='output_dir', help='Provide output file dir', required=True)
    parser.add_argument('--prefix', dest='prefix', help='Provide the plot prefix name', default='output')
    args = parser.parse_args()

    fq_zip_1 = args.fq_zip_1
    fq_zip_2 = args.fq_zip_2
    output_dir = args.output_dir
    prefix = args.prefix

    os.environ['fq_zip_1'] = fq_zip_1
    os.environ['fq_zip_2'] = fq_zip_2

    fq_1_dir = re.sub('\.zip', '', fq_zip_1)
    fq_2_dir = re.sub('\.zip', '', fq_zip_2)
    fq_p_dir = os.path.dirname(fq_1_dir)
    os.environ['fq_p_dir'] = fq_p_dir


    if not os.path.isdir(fq_1_dir):
        os.system("unzip ${fq_zip_1} -d ${fq_p_dir}")
    if not os.path.isdir(fq_2_dir):
        os.system("unzip ${fq_zip_2} -d ${fq_p_dir}")

    fq_1_summary = os.path.join(fq_1_dir, "fastqc_data.txt")
    fq_2_summary = os.path.join(fq_2_dir, "fastqc_data.txt")

    gc_output_tmp = os.path.join(output_dir, 'atcg_content.tmp')
    quality_output_tmp = os.path.join(output_dir, 'quality.tmp')
    gc_output = os.path.join(output_dir, prefix+'_gc_content')
    quality_output = os.path.join(output_dir, prefix+'_quality')
    with open(quality_output_tmp, 'w') as fo2:
        with open(gc_output_tmp, 'w') as fo:
            get_info(fq_1_summary)
            get_info(fq_2_summary)

    with open(gc_output, 'w') as fo3:
        with open(gc_output_tmp) as f:
            i = 0
            for line in f:
                new_line = str(i) + '\t' + line
                i += 1
                fo3.write(new_line)
    with open(quality_output, 'w') as fo4:
        with open(quality_output_tmp) as f:
            i = 0
            for line in f:
                new_line = str(i) + '\t' + line
                i += 1
                fo4.write(new_line)

    # GC_plot
    read_length = 150
    total_length = read_length*2
    half_length = read_length/2
    output1_png = gc_output+'.'+'png'
    vertical_bar = "abline(v={}, col='grey', lty=2)".format(read_length)

    gc_plot_cmd = """
gc<- read.table({})
site<- gc[, 1]
base_a <- gc[, 4]
base_t <- gc[, 7]
base_g <- gc[, 10]
base_c <- gc[, 13]
base_n <- gc[, 16]
total_sites <- {}
half_sites<- {}
png({}, width=800, height=600, units="px")
sp_a <- spline(site, base_a, n=1000)
sp_t <- spline(site, base_t, n=1000)
sp_g <- spline(site, base_g, n=1000)
sp_c <- spline(site, base_c, n=1000)
sp_n <- spline(site, base_n, n=1000)
plot(sp_a, xlim=c(0, total_sites), ylim=c(0, 50), axes=FALSE, col="red", type="l", xlab="Position along reads",
 ylab="percent", main="Base percentage composition along reads", lty=1, lwd=3.3)
lines(sp_t, col="magenta", type="l", lty=2, lwd=3)
lines(sp_g, col="darkblue", type="l", lty=3, lwd=3)
lines(sp_c, col="green", type="l", lty=4, lwd=3)
lines(sp_n, col="cyan3", type="l", lty=5, lwd=3)
legend("topright", legend=c("A", "T", "G", "C", "N"), col=c("red", "magenta", "darkblue", "green", "cyan3"),lty=c(1, 2, 3, 4, 5))
{}
axis(side=1, at=seq(from=0, to = total_sites, by = half_sites))
axis(side=2, at=seq(from=0, to = 50, by = 10))
dev.off()
    """.format('"'+gc_output+'"', total_length, half_length, '"'+output1_png+'"', vertical_bar)
    print(gc_plot_cmd)
    # mean quality plot
    output2_png = quality_output+'.png'
    mean_quality_plot_cmd = """
table<- read.table({})
site<- table[, 1]
quality<- table[, 2]
#error<- table[, 3]
total_sites<- {}
png({}, width=800, height=600, units="px")
plot(site, quality, xlim=c(0, total_sites), ylim=c(0, 40), axes=FALSE, col="red", type="p", pch=20, cex=0.35,
 xlab="Position along reads", ylab="Quality", main="Distribution of qualities")
axis(side=1, at=seq(from=0, to = total_sites, by = 20))
axis(side=2, at=seq(from=0, to = 40, by = 10))
abline(h=20, col="darkblue", lty=2)
abline(v=seq(0, total_sites, by=10), col="darkblue", lty=3)
dev.off()""".format('"'+quality_output+'"', total_length, "'"+output2_png+"'")

    os.environ['gc_cmd'] = gc_plot_cmd
    os.environ['quality_cmd'] = mean_quality_plot_cmd
    Rscripts = os.path.join(output_dir, 'fq_QC_plot.tmp.R')
    with open(Rscripts, 'w') as fo:
        fo.write(mean_quality_plot_cmd+'\n')
        fo.write(gc_plot_cmd+'\n')
    os.environ['Rscripts'] = Rscripts
    os.system("Rscript ${Rscripts}")
    os.system("rm ${Rscripts}")
    os.environ['gc_output_tmp'] = gc_output_tmp
    os.environ['quality_output_tmp'] = quality_output_tmp
    os.system("rm ${quality_output_tmp} ${gc_output_tmp}")
