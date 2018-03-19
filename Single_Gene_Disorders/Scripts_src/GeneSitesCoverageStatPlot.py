# _*_ coding:utf-8 _*_

import re
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='INFO: Calculate the coverage of disease related gene and sites and plot them!', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-c', dest='coverage_file', help='Provide the coverage_file', required=True)
parser.add_argument('-g', dest='gene_omim_file', help='Provide the gene omim configure file!', required=True)
parser.add_argument('-o', dest='output_dir', help='Provide the output_dir!', required=True)
parser.add_argument('-p', dest='prefix', help='Provide the prefix name for the intermediate config files or dirs!', required=True)
args = parser.parse_args()
coverage_file = args.coverage_file
gene_omim_file = args.gene_omim_file
output_dir = args.output_dir
prefix = args.prefix

plot_dir = os.path.join(output_dir, prefix)
if not os.path.isdir(plot_dir):
    try:
        os.makedirs(plot_dir)
    except OSError as e:
        print('ERROR: You seemed not have the permission to create the dir: {}, Check it!'.format(plot_dir))
        sys.exit()

with open(gene_omim_file) as f:
    gene_mim_dic = {}
    for line in f:
        line_list = re.split(':|;|\t|,', line.strip())
        gene = line_list[0]
        mim = line_list[1]
        gene_mim = '_'.join([gene, mim])
        gene_mim_dic[gene] = mim

# 构建位点和基因统计文件
output_config = os.path.join(plot_dir, os.path.basename(coverage_file)+'.plot.config')
with open(output_config, 'w') as fo1:
    with open(coverage_file) as f:
        output_config_dic = {}
        gene_exon_num_dic = {}
        col_dic = {}
        for line in f:
            line_list = line.strip().split('\t')
            if line.startswith('#'):
                for item in line_list:
                    col_dic[item] = line_list.index(item)
            else:
                site_id, exon_id = line_list[col_dic['ID']], line_list[col_dic['ExonID']]
                dp, db = line_list[col_dic['Depth']], line_list[col_dic['DB_Source']]
                chrid, pos, relative_pos = line_list[0], line_list[col_dic['Pos']], line_list[col_dic['RelativePos']]
                gene, mim = line_list[col_dic['Gene']], line_list[col_dic['MIM']]
                gene_output = os.path.join(output_dir, prefix+'_'+gene+'.gene.stat')
                sites_output = os.path.join(output_dir, prefix+'_'+gene+'.db.sites.stat')
                if gene in gene_mim_dic:
                    gene_mim = '_'.join([gene, gene_mim_dic[gene]])
                    config_line1 = '\t'.join([gene_output, sites_output, gene_mim])

                    output_config_dic[gene] = config_line1
                    gene_exon_num_dic[gene] = exon_id

                    fh1 = open(gene_output, 'a')
                    if int(exon_id) % 2 == 1:
                        colour = 'green4'
                    else:
                        colour = 'blue'
                    new_line1 = '\t'.join([gene, exon_id, chrid, pos, dp, relative_pos, colour]) + '\n'
                    fh1.write(new_line1)

                    if (mim in gene_mim_dic[gene]) and (not site_id.startswith('-')):
                        fh2 = open(sites_output, 'a')
                        new_line2 = '\t'.join([gene, exon_id, chrid, pos, dp, relative_pos, 'red']) + '\n'
                        fh2.write(new_line2)

    # 画图配置文件
    for gene in output_config_dic:
        part1 = output_config_dic[gene]
        part2 = gene_exon_num_dic[gene]
        new_line3 = '\t'.join([part1, part2])+'\n'
        fo1.write(new_line3)

plot_cmd1 = 'cov <- as.matrix(read.table("{}", header=F, sep="\t"))'.format(output_config)
plot_cmd2 = """
for (i in 1: dim(cov)[1]) {
"""
plot_cmd3 = """
    pdf(paste("{}/", cov[i, 3], ".pdf", sep=""), height=3, width=8)
""".format(plot_dir)
plot_cmd4 = """
    par(mfrow=c(1, 1), mar=c(3.5, 4, 2, 1), oma=c(0, 0, 0, 0), xaxs="i", yaxs="i", bty="l")

    data1 <- read.table(cov[i, 1], header=T, sep="\t")
    data2 <- read.table(cov[i, 2], header=T, sep="\t")
    km <- which(data1[, 5] >= 10)
    rate <- length(km) / dim(data1)[1] * 100
    rate <- round(rate, 2)
    da1 <- as.matrix(data1)
    da2 <- as.matrix(data2)
    plot(data1[, 6], data1[, 5], main=cov[i, 3], col=da1[, 7], col.axis="black", col.lab="black", col.main="darkred", xlab="", ylab="Depth", xaxt="n", type='h', cex.lab=0.8, ylim=c(-1, max(data1[, 5])), lwd=0.8)
    # mtext(paste("(DP>=10: ",rate,"%)",sep=""),side=3,line=-1.5,cex=0.85,las=1,col="black")
    text(x=max(data1[, 6])*0.9, y=max(data1[, 5])*1.1, adj=1, labels=paste("(DP>=10: ", rate, "%)", sep=""), srt=0, xpd=TRUE, cex=0.82)
    points(data2[, 6], data2[, 5], col = da2[, 7], type = 'p', pch = 19, cex = 1)
    xm1 <- c(min(data1[, 6]), max(data1[, 6]))
    xm2 <- c(min(data1[, 4]), max(data1[, 4]))
    # n.x<-length(xm1)
    # par(xaxt="s")
    # for (j in 1:n.x) {
    # axis(side=1,at=xm1[j],tck=-0.02,labels=xm2[j],col.axis="black",cex.axis=0.8,las=1,hadj=0.5,padj=-1,line=0.05)
    # }
    tmm <- data.frame(table(data1[, 2]))
    pp1 <- paste("Exon", tmm[, 1])
    pp2 <- c()
    for (i in 1: dim(tmm)[1]){
        pk <- data1[which(data1[, 2] == tmm[i, 1]), 6]
        lines(pk, rep(10, length(pk)), type='l', lwd=1.5, col="grey", lty=2)
        x0 <- (min(pk) + max(pk)) / 2
        pp2 <- c(pp2, x0)
        }
    text(x=pp2, y=-5, adj=1, labels=pp1, srt=45, xpd=TRUE, cex=0.82)
}
dev.off()
"""

RplotScripts = os.path.join(plot_dir, prefix+'.Stat.plot.R')
with open(RplotScripts, 'w') as fo2:
    fo2.write(plot_cmd1+'\n')
    fo2.write(plot_cmd2+'\n')
    fo2.write(plot_cmd3+'\n')
    fo2.write(plot_cmd4+'\n')
os.environ['Rscripts'] = RplotScripts
os.system("Rscript ${Rscripts}")
