#! /usr/bin/env python3
# _*_ coding:utf-8 _*_

import sys
import os
import re
import gzip
import glob
import time
import argparse
import configparser

# -----------------构造默认配置文件目录路径-------------------------
cwd = os.path.dirname(os.path.realpath(__file__))
cwd_p = os.path.dirname(cwd)
config_dir = os.path.join(cwd_p, 'Config_files')
tmp_dir = os.path.join(cwd_p, 'tmp')
metrics_statistics_scripts = os.path.join(cwd, 'Metrics_statistic.sh')
vcf_muttype_stat_scripts = os.path.join(cwd, 'VcfStatistic.py')
fastqc_plot_scripts = os.path.join(cwd, 'fastqc_plot.py')
omim_genes_bed = os.path.join(config_dir, 'OMIM_genes_hg19_exons.bed')

# -----------------参数选项设置------------------------------------
parser = argparse.ArgumentParser(description='APPLICATION: Germline Calling Process From Raw Fq files and Basic data statistics!', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-I', dest='fastq_dir', help='Provide the fastq files dir!', required=True)
parser.add_argument('-P', dest='patient_id', help='Provide the patient_id !', required=True)
parser.add_argument('-L', dest='lane_id', help='Provide the lane_id, or data batch_id!', required=True)
parser.add_argument('-C', dest='config_files_dir', help='Provide the config files dir, default ../Config_files/!', default=config_dir)
parser.add_argument('-O', dest='out_path', help='Provide the output direction!')
args = parser.parse_args()
fastq_dir = args.fastq_dir
patient_id = args.patient_id
lane_id = args.lane_id
config_files_dir = args.config_files_dir
out_path = args.out_path.rstrip('/')

# ---------------------建立文件目录系统----------------------
dirs_list = ['0-LOG', '1-FASTQC', '2-MAPPING', '3-PICARD', '4-GERMLINE', '5-STATISTICS']
for path_ in dirs_list:
    abs_path = os.path.join(out_path, path_)
    try:
        os.makedirs(abs_path)
    except OSError:
        pass

# -------------------组学工具路径, 输入变量及log日志文件配置--------------------------
# genomic tools 路径配置
genomic_tools_config = os.path.join(config_files_dir, 'genomic_tools_env.ini')
config = configparser.ConfigParser()
config.read(genomic_tools_config)
my_config = config['conf']
for item in my_config:
    os.environ[item] = my_config[item]
# 外部输入变量转换为系统环境变量
env_conf = dict()
env_conf['SM'] = patient_id
env_conf['lane'] = lane_id
env_conf['out_path'] = out_path.rstrip('/')
env_conf['fastq_dir'] = fastq_dir.rstrip('/')
env_conf['config_files_dir'] = config_files_dir
# paired-end fastq文件绝对路径确定及转换为系统环境变量
fq_1_regex = fastq_dir.rstrip('/')+'/'+patient_id+'/'+lane_id+'/'+patient_id+'*1*fq*'
fq_2_regex = fastq_dir.rstrip('/')+'/'+patient_id+'/'+lane_id+'/'+patient_id+'*2*fq*'
fq_1_glob = glob.glob(fq_1_regex)
fq_2_glob = glob.glob(fq_2_regex)
if fq_1_glob:
    fq_1 = fq_1_glob[0]
else:
    print('ERROR: NOT FIND THE Forward fq file, Please check that!')
    sys.exit()
if fq_2_glob:
    fq_2 = fq_2_glob[0]
else:
    print('ERROR: NOT FIND THE Reverse fq file, Please check that!')
    sys.exit()
env_conf['fq_1'] = fq_1
env_conf['fq_2'] = fq_2
# 获取fastq文件的第一行信息，并提取出LB, ID, PU信息
if fq_1.endswith('gz'):
    first_line = gzip.open(fq_1).readline()
else:
    first_line = open(fq_1).readline()
LB = first_line.split(':')[2]
ID = LB + '.' + first_line.split(':')[3]
PU = ID + '.' + patient_id
PL = 'Illumina'
env_conf['LB'] = LB
env_conf['ID'] = ID
env_conf['PU'] = PU
env_conf['PL'] = PL

# 各步骤log文件环境变量配置
env_conf['bwa_log'] = os.path.join(out_path, '0-LOG/'+patient_id+'.bwa.log')
env_conf['picard_log'] = os.path.join(out_path, '0-LOG/'+patient_id+'.picard.log')
env_conf['haplotypecaller_log'] = os.path.join(out_path, '0-LOG/'+patient_id+'.haplotypecaller.log')
env_conf['fastqc_log'] = os.path.join(out_path, '0-LOG/'+patient_id+'.fastqc.log')
env_conf['trimmomatic_log'] = os.path.join(out_path, '0-LOG/'+patient_id+'.trimmomatic.log')
env_conf['anno_log'] = os.path.join(out_path, '0-LOG/'+patient_id+'.anno.log')
env_conf['stat_log'] = os.path.join(out_path, '0-LOG/'+patient_id+'.stat.log')
env_conf['error_log'] = os.path.join(out_path, '0-LOG/'+patient_id+'.ERROR.log')

whole_log_file = open(os.path.join(out_path, '0-LOG/'+patient_id+'.pipeline.log'), 'w')
error_log_file = open(env_conf['error_log'], 'w')

# In house统计脚本环境变量配置
env_conf['metric_stat_script'] = metrics_statistics_scripts
env_conf['vcf_stat_script'] = vcf_muttype_stat_scripts
env_conf['fastqc_plot_script'] = fastqc_plot_scripts

# omim-genes bed regions file
env_conf['omim_genes_bed'] = omim_genes_bed

for elem in env_conf:
    os.environ[elem] = env_conf[elem]

# --------------构建log判断函数，判断某一步骤是正常完成还是处于异常状态------------------


def filter_log(error_list, tail_list, log_file):
    status = 1
    if os.path.isfile(log_file):
        with open(log_file) as f:
            log_str = f.read()
            for error_info in error_list:
                if error_info in log_str:
                    status = 0
                    break
            if status:
                for tail_info in tail_list:
                    if tail_info not in log_str:
                        status = 0
                        break
    else:
        status = 0
    return status

# ------------------定义vcf位点质量过滤标准函数---------------------------
class FilterVcf:
    """
    VCF文件在进入流程之前对其按照0.2 <alt/total<0.8, depth>10进行过滤，保证位点的准确性
    """
    def __init__(self, raw_vcf, filtered_vcf):
        """
        :param raw_vcf: 单个人的vcf文件，待过滤
        :param output_dir: 过滤后的输出文件所在目录，输出文件为后缀为site_quality_filtered.vcf
        """
        self.raw_vcf = raw_vcf
        self.filtered_vcf = filtered_vcf

    def filter_vcf(self):
        with open(self.raw_vcf, 'r') as f_in:
            with open(self.filtered_vcf, 'w') as fo:
                for line in f_in:
                    if line.startswith('#'):
                        fo.write(line)
                    else:
                        line_list = line.strip().split('\t')
                        info_field = line_list[9]
                        if info_field.startswith('0/0') or info_field.startswith('./.'):
                            pass
                        elif info_field.startswith('1/1'):
                            DP = info_field.split(':')[2]
                            if DP == '.':
                                pass
                            elif int(DP) > 10:
                                fo.write(line)
                        else:
                            AD = info_field.split(':')[1]
                            DP = info_field.split(':')[2]
                            ref = int(AD.split(',')[0])
                            alt = int(AD.split(',')[1])
                            if not (ref + alt) == 0 and not DP == '.':
                                fraction = alt/(ref + alt)
                                if 0.2 < fraction < 0.8 and int(DP) > 10:
                                    fo.write(line)
# ------------------QC STEP，包括fastqc > trimmomatic > fastqc----------------------

print('INFO: Raw fq files FASTQC started at: {}'.format(time.ctime()), file=whole_log_file, flush=True)
os.system("fastqc -t 2 -o ${out_path}/1-FASTQC/ ${fq_1} ${fq_2} 2> ${fastqc_log} 1>&2")
print('INFO: Raw fq files FASTQC finished at: {}'.format(time.ctime()), file=whole_log_file, flush=True)
trimmomatic_cmd_str_file = os.path.join(config_files_dir, 'trimmomatic_cmd_str')
trimmomatic_output_path = os.path.join(out_path, '1-FASTQC')
try:
    os.makedirs(trimmomatic_output_path)
except OSError:
    pass

with open(trimmomatic_cmd_str_file) as f:
    trimmomatic_cmd_str = f.read()

fq_1_trimmed = os.path.join(trimmomatic_output_path, patient_id + '.'+lane_id+'.R1.clean.fq')
fq_2_trimmed = os.path.join(trimmomatic_output_path, patient_id + '.'+lane_id+'.R2.clean.fq')
fq_1_unpaired = os.path.join(trimmomatic_output_path, patient_id + '.'+lane_id+'.R1.unpaired.fq')
fq_2_unpaired = os.path.join(trimmomatic_output_path, patient_id + '.'+lane_id+'.R2.unpaired.fq')

trimmomatic_cmd = trimmomatic_cmd_str.format(fq_1, fq_2, fq_1_trimmed, fq_1_unpaired, fq_2_trimmed, fq_2_unpaired, env_conf['trimmomatic_log'])
os.system(trimmomatic_cmd)
print('INFO: Raw fq files TRIMMOMATIC step finished at: {}'.format(time.ctime()), file=whole_log_file, flush=True)

os.environ['fq_1_trimmed'] = fq_1_trimmed
os.environ['fq_2_trimmed'] = fq_2_trimmed

os.system("fastqc -t 2 -o ${out_path}/1-FASTQC/ ${fq_1_trimmed} ${fq_2_trimmed} 2>> ${fastqc_log} 1>&2")
print('INFO: Clean fq files FASTQC finished at: {}'.format(time.ctime()), file=whole_log_file, flush=True)

# -----------------------bwa > samtools > picard ------------------------------------
# mapping and sort
bwa_cmd = "bwa mem -t 2 -R '@RG\\tID:{}\\tSM:{}\\tPU:{}\\tLB:{}\\tPL:{}' -M {}/hg19 {} {} 2>{}|samtools sort -T {} -O BAM  -o {}/2-MAPPING/{}.{}.sorted.bam 2>{} 1>&2".format(patient_id, patient_id, PU, LB, PL, my_config['hgRef'], fq_1_trimmed, fq_2_trimmed, env_conf['bwa_log'], tmp_dir+'/'+env_conf['SM']+'.'+env_conf['lane'], out_path, patient_id, lane_id, env_conf['bwa_log'])
os.system(bwa_cmd)
bwa_error_list = ['ERROR', 'error', 'Error', 'truncated file']
filter_bwa = filter_log(bwa_error_list, ['[main] Real time:'], env_conf['bwa_log'])
if filter_bwa:
    print('INFO: bwa and sort procedure finished at {}'.format(time.ctime()), file=whole_log_file, flush=True)
else:
    print('ERROR: Sample {} is collapsed at the bwa & sort step on {}'.format(env_conf['SM'], time.ctime()), file=error_log_file, flush=True)
    sys.exit()

# markduplicates
sorted_bam = os.path.join(out_path, '2-MAPPING/'+patient_id+'.'+lane_id+'.sorted.bam')
marked_bam = os.path.join(out_path, '3-PICARD/'+patient_id+'.'+lane_id+'.marked.bam')
markdup_metrics = os.path.join(out_path, '3-PICARD/'+patient_id+'.'+lane_id+'.Mkdup.metrics')
markdup_cmd = "java -Xmx5g -jar {} MarkDuplicates VALIDATION_STRINGENCY=SILENT INPUT={} OUTPUT={} METRICS_FILE={} 2> {} 1>&2".format(my_config['picardjar'], sorted_bam, marked_bam, markdup_metrics, env_conf['picard_log'])
os.system(markdup_cmd)
filter_picard = filter_log(['ERROR', 'problem', 'invalid'], ['picard.sam.markduplicates.MarkDuplicates done'], env_conf['picard_log'])
if filter_picard:
    print('INFO: picard marduplicates finished at {}'.format(time.ctime()), file=whole_log_file, flush=True)
else:
    print('ERROR: Sample {} is collapsed at the picard marduplicates step on {}'.format(env_conf['SM'], time.ctime()), file=error_log_file, flush=True)
    sys.exit()

os.system("samtools index ${out_path}/3-PICARD/${SM}.${lane}.marked.bam 2>> ${picard_log} 1>&2")
print('INFO: samtools index step finished at {}'.format(time.ctime()), file=whole_log_file, flush=True)

os.system("samtools view -b ${out_path}/3-PICARD/${SM}.${lane}.marked.bam chrY > ${out_path}/3-PICARD/${SM}.${lane}.chrY.bam")
print('INFO: samtools pick chrY mapping results at {}'.format(time.ctime()), file=whole_log_file, flush=True)

# germline calling
# os.system("java -Xmx5g -jar ${gatkjar} -T HaplotypeCaller -R ${hgref}/hg19.fa -L ${intervallist} --dbsnp ${dbsnp}  -I ${out_path}/3-PICARD/${SM}.${lane}.marked.bam -o ${out_path}/4-GERMLINE/${SM}.${lane}.germline.vcf 2> ${haplotypecaller_log} 1>&2")
filter_gatk = filter_log(['ERROR', 'problem', 'invalid'], ['ProgressMeter - Total runtime'], env_conf['haplotypecaller_log'])
if filter_gatk:
    print('INFO: GATK HaplotypeCaller finished at {}'.format(time.ctime()), file=whole_log_file, flush=True)
else:
    print('ERROR: Sample {} is collapsed at the GATK HaplotypeCaller step on {}'.format(env_conf['SM'], time.ctime()), file=error_log_file, flush=True)
    sys.exit()
# vcf sites quality filtering
raw_vcf = os.path.join(out_path, '4-GERMLINE/'+env_conf['SM']+'.'+env_conf['lane']+'.germline.vcf')
filter_vcf = os.path.join(out_path, '4-GERMLINE/'+env_conf['SM']+'.'+env_conf['lane']+'.germline_filter.vcf')
filter_func = FilterVcf(raw_vcf, filter_vcf)
filter_func.filter_vcf()

# annotation
os.system("${annovar}/table_annovar.pl ${out_path}/4-GERMLINE/${SM}.${lane}.germline_filter.vcf ${humandb} --buildver hg19 --outfile ${out_path}/4-GERMLINE/${SM}.${lane}.germline_filter --remove --protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp150,dbnsfp33a,clinvar_20170905,exac03,dbscsnv11,dbnsfp31a_interpro,rmsk,ensGene,knownGene  --operation  g,f,f,f,f,f,f,f,f,r,g,g --nastring . -vcfinput 2>${anno_log} 1>&2")
filter_annovar = filter_log([], ['VCF output is written to', 'Multianno output file is written to'], env_conf['anno_log'])
if filter_annovar:
    print('INFO: Annovar finished at {}'.format(time.ctime()), file=whole_log_file, flush=True)
else:
    print('ERROR: Sample {} is collapsed at the annovar step on {}'.format(env_conf['SM'], time.ctime()), file=error_log_file, flush=True)
    sys.exit()

# intervar evaluation
os.system("python2 ${intervar} -b hg19 -i ${out_path}/4-GERMLINE/${SM}.${lane}.germline_filter.vcf --input_type=VCF -o ${out_path}/4-GERMLINE/${SM}.${lane}.germline_filter --skip_annovar 2>>${anno_log} 1>&2")
filter_intervar = filter_log([], ['The InterVar is finished'], env_conf['anno_log'])
if filter_intervar:
    print('INFO: InterVar finished at {}'.format(time.ctime()), file=whole_log_file, flush=True)
else:
    print('ERROR: Sample {} is collapsed at the intervar step on {}'.format(env_conf['SM'], time.ctime()), file=error_log_file, flush=True)
    sys.exit()

# data statistics
os.system("java -Xmx5g -jar ${picardjar} CalculateHsMetrics I=${out_path}/3-PICARD/${SM}.${lane}.marked.bam O=${out_path}/5-STATISTICS/${SM}.${lane}.marked.bam.metrics BAIT_INTERVALS=${intervallist} TARGET_INTERVALS=${intervallist} 2>> ${picard_log} 1>&2")
filter_metrics = filter_log(['ERROR:', 'USAGE: CalculateHsMetrics'], [], env_conf['picard_log'])
if filter_metrics:
    print('INFO: Picard CalculateHsMetrics finished at {}'.format(time.ctime()), file=whole_log_file, flush=True)
else:
    print('ERROR: Sample {} is collapsed at the Picard CalculateHsMetrics step on {}'.format(env_conf['SM'], time.ctime()), file=error_log_file, flush=True)
    sys.exit()
os.system("${metric_stat_script} -M ${out_path}/5-STATISTICS/${SM}.${lane}.marked.bam.metrics -O ${out_path}/5-STATISTICS/ -P ${SM}.${lane} 2>${stat_log} 1>&2")
os.system("python ${vcf_stat_script} ${out_path}/4-GERMLINE/${SM}.${lane}.germline_filter.hg19_multianno.txt ${out_path}/5-STATISTICS/${SM}.${lane}.stat.txt 2>${stat_log} 1>&2")
fq_1_zip = os.path.join(os.path.dirname(fq_1_trimmed), os.path.basename(fq_1))
fq_2_zip = os.path.join(os.path.dirname(fq_2_trimmed), os.path.basename(fq_2))
raw_fq_1_fastqc_zip = re.sub('\.fq|\.fastq|\.fq\.gz|\.fastq\.gz', '_fastqc.zip', fq_1_zip)
raw_fq_2_fastqc_zip = re.sub('\.fq|\.fastq|\.fq\.gz|\.fastq\.gz', '_fastqc.zip', fq_2_zip)
clean_fq_1_fastqc_zip = re.sub('\.fq|\.fastq|\.fq\.gz|\.fastq\.gz', '.clean_fastqc.zip', fq_1_zip)
clean_fq_2_fastqc_zip = re.sub('\.fq|\.fastq|\.fq\.gz|\.fastq\.gz', '.clean_fastqc.zip', fq_1_zip)
os.environ['raw_fq_1_fastqc_zip'] = raw_fq_1_fastqc_zip
os.environ['raw_fq_2_fastqc_zip'] = raw_fq_2_fastqc_zip
os.environ['clean_fq_1_fastqc_zip'] = clean_fq_1_fastqc_zip
os.environ['clean_fq_2_fastqc_zip'] = clean_fq_2_fastqc_zip
os.system('python3 ${fastqc_plot_script} --input1 ${raw_fq_1_fastqc_zip} --input2 ${raw_fq_2_fastqc_zip} --output ${out_path}/5-STATISTICS --prefix ${SM}.${lane}.raw 2>>${stat_log} 1>&2')
os.system('python3 ${fastqc_plot_script} --input1 ${clean_fq_1_fastqc_zip} --input2 ${clean_fq_2_fastqc_zip} --output ${out_path}/5-STATISTICS --prefix ${SM}.${lane}.clean 2>>${stat_log} 1>&2')
os.system('samtools depth -b ${omim_genes_bed} ${out_path}/3-PICARD/${SM}.${lane}.marked.bam > ${out_path}/5-STATISTICS/${SM}.${lane}.omim_genes.depth')
print('INFO: Statistic step finished at {}'.format(time.ctime()), file=whole_log_file, flush=True)

# 关闭文件句柄
whole_log_file.close()
error_log_file.close()

