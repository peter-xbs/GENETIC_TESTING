#! /usr/bin/python
# -*- coding: utf-8 -*-
import sys
import os
import gzip
import time
import re
import commands
#import send_mail
import ConfigParser
import getopt

def usage():
    print 'USAGE INFO: python2 BasicProcess.py -I <fastq files dir> -P <patient sample id> -L <batch or lane id> -C <configure files dir> -O <output_path> -h <help>'

opts, args = getopt.getopt(sys.argv[1:], "I:P:L:C:O:h")
if len(opts) < 2:
    usage()
    sys.exit()
for op, value in opts:
    if op == '-I':
        fastq_dir = value
    elif op == '-P':
        patient_id = value
    elif op == '-L':
        lane_id = value
    elif op == '-C':
        config_files_dir = value
    elif op == '-O':
        out_path = value
    elif op == '-h':
        usage()
        sys.exit()
    else:
        print('Please Provide the right arguments:')
        usage()
        sys.exit()

class myconf(ConfigParser.ConfigParser):  
    def __init__(self,defaults=None):  
        ConfigParser.ConfigParser.__init__(self,defaults=None)  
    def optionxform(self, optionstr):  
        return optionstr

def getDocSize(filenam):
    try:
        size = os.path.getsize(filenam)
        return size
    except Exception as err:
        return 0

def filter_log(error_info_list,tail_info,log_file):
    if (getDocSize(log_file)>0):
        for error_info in error_info_list:
            if filter(lambda x:error_info in x,file(log_file)):
	        return
        else:
            return 1
        if not filter(lambda x:tail_info in x,file(log_file)):
            return 
        else:
            return 1
    else:
        return 0
conf={}
conf['SM']= patient_id 
conf['lane'] = lane_id
genomic_config_file = config_files_dir.strip('/')+'/genomic_tools_env.ini'
cf=myconf()
cf.read(genomic_config_file)
v = cf.items("conf")
for i in range(0,len(v)):
    conf[v[i][0]]=v[i][1]

conf['out_path'] = out_path
conf['Fq_position'] = fastq_dir
conf['file1_1'] = conf['Fq_position']+"/"+conf['SM']+"/"+conf['lane']+"/"+conf['SM']+".R1.fq.gz"
conf['file1_2'] = conf['Fq_position']+"/"+conf['SM']+"/"+conf['lane']+"/"+conf['SM']+".R1.fq"
conf['file1_3'] = conf['Fq_position']+"/"+conf['SM']+"/"+conf['lane']+"/"+conf['SM']+"_1.fq.gz"
conf['file1_4'] = conf['Fq_position']+"/"+conf['SM']+"/"+conf['lane']+"/"+conf['SM']+"_1.fq"
conf['file2_1'] = conf['Fq_position']+"/"+conf['SM']+"/"+conf['lane']+"/"+conf['SM']+".R2.fq.gz"
conf['file2_2'] = conf['Fq_position']+"/"+conf['SM']+"/"+conf['lane']+"/"+conf['SM']+".R2.fq"
conf['file2_3'] = conf['Fq_position']+"/"+conf['SM']+"/"+conf['lane']+"/"+conf['SM']+"_2.fq.gz"
conf['file2_4'] = conf['Fq_position']+"/"+conf['SM']+"/"+conf['lane']+"/"+conf['SM']+"_2.fq"
conf['bwa_log'] = conf['out_path'] + '/0-LOG/' + conf['SM'] + '.bwa.log'
conf['samtools_sort_log']=conf['out_path']+'/0-LOG/'+conf['SM']+'.samtools_sort.log'
conf['picard_merge_log']=conf['out_path']+'/0-LOG/'+conf['SM']+'.picard_merge.log'
conf['markduplicates_log']=conf['out_path']+'/0-LOG/'+conf['SM']+'.MarkDuplicates.log'
conf['indexmarkduplicates_log']=conf['out_path']+'/0-LOG/'+conf['SM']+'.Index.markduplicates.log'
conf['indelrealigner_log']=conf['out_path']+'/0-LOG/'+conf['SM']+'.GATK.IndelRealigner.log'
conf['baserecalibrator_log']=conf['out_path']+'/0-LOG/'+conf['SM']+'.GATK.BaseRecalibrator.log'
conf['realignertargetcreator_log']=conf['out_path']+'/0-LOG/'+conf['SM']+'.GATK.RealignerTargetCreator.log'
conf['printreads_log']=conf['out_path']+'/0-LOG/'+conf['SM']+'.GATK.PrintReads.log'
conf['haplotypecaller_log']=conf['out_path']+'/0-LOG/'+conf['SM']+'.GATK.HaplotypeCaller.log'
conf['calculatehsmetrics_log']=conf['out_path']+'/0-LOG/'+conf['SM']+'.CalculateHsMetrics.log'

for key in conf:
    os.environ[key] = conf[key]

out_path_list = ['/0-LOG','/1-FASTQC','/2-MAPPING','/3-PICARD','/4-GATK','/4-GATK/' + conf['lane'],'/5-GERMLINE','/6-HsMetrics', '/7-QC_Stat']
for path_bar in out_path_list:
    path=conf['out_path']+path_bar
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass

logfile = open(conf['out_path']+'/0-LOG/'+conf['SM']+'.log','w')
print >> logfile,'Start at ' + time.ctime()
error_log_fq = open(conf['out_path']+'/0-LOG/error_fq.log','a')
error_log_bwa = open(conf['out_path']+'/0-LOG/error_bwa.log','a')
error_log_mergesamfiles = open(conf['out_path']+'/0-LOG/error_mergesamfiles.log','a')
error_log_markduplicate = open(conf['out_path']+'/0-LOG/error_MarkDuplicates.log','a')
error_log_indelrealigner = open(conf['out_path']+'/0-LOG/error_IndelRealigner.log','a')
error_log_baserecalibrator = open(conf['out_path']+'/0-LOG/error_BaseRecalibrator.log','a')
error_log_realignertargetcreator = open(conf['out_path']+'/0-LOG/error_RealignerTargetCreator.log','a')
error_log_printreads = open(conf['out_path']+'/0-LOG/error_PrintReads.log','a')
error_log_haplotypecaller = open(conf['out_path']+'/0-LOG/error_HaplotypeCaller.log','a')
error_log_calculatehsmetrics = open(conf['out_path']+'/0-LOG/error_CalculateHsMetrics.log','a')

if os.path.exists(conf['file1_1']):
    fq_1 = conf['file1_1']
    fq_2 = conf['file2_1']
    line = gzip.open(fq_1).readline()
elif os.path.exists(conf['file1_2']):
    fq_1 = conf['file1_2']
    fq_2 = conf['file2_2']
    line = open(fq_1).readline()   
elif os.path.exists(conf['file1_3']):
    fq_1 = conf['file1_3']
    fq_2 = conf['file2_3']
    line = gzip.open(fq_1).readline()
elif os.path.exists(conf['file1_4']):
    fq_1 = conf['file1_4']
    fq_2 = conf['file2_4']
    line = open(fq_1).readline()
else:
    print >> error_log_fq,conf['SM']
    print 'Please provide the right input file dirs/files!'
    #send_mail.send_log_mail('Wrong Fq_position : does not have fastq file','Wrong_Fq_'+conf['SM'])
    os._exit(1)

os.environ['fq_1'] = fq_1
os.environ['fq_2'] = fq_2
print(fq_1, fq_2)
os.system("fastqc -t 2 -o ${out_path}/1-FASTQC/ ${fq_1} ${fq_2} 2> ${out_path}/0-LOG/${SM}.${lane}.fastqc.log 1>&2")

trimmomatic_cmd_str_file = config_files_dir.strip('/') + '/trimmomatic_cmd_str'
with open(trimmomatic_cmd_str_file) as f:
    trimmomatic_cmd_str = f.read()
fq_1_trimmed = os.path.join(out_path.rstrip('/')+'/1-FASTQC/'+conf['SM']+'/'+conf['lane'] , conf['SM'] + '.'+conf['lane']+'.R1.clean.fq')
fq_2_trimmed = os.path.join(out_path.rstrip('/')+'/1-FASTQC/'+conf['SM']+'/'+conf['lane'] , conf['SM'] + '.'+conf['lane']+'.R2.clean.fq')
fq_1_unpaired = os.path.join(out_path.rstrip('/')+'/1-FASTQC/'+conf['SM']+'/'+conf['lane'] , conf['SM'] + '.'+conf['lane']+'.R1.unparied.fq')
fq_2_unpaired = os.path.join(out_path.rstrip('/')+'/1-FASTQC/'+conf['SM']+'/'+conf['lane'] , conf['SM'] + '.'+conf['lane']+'.R2.unparied.fq')
try:
    os.system('mkdir -p ${out_path}/1-FASTQC/${SM}/${lane}')
except Exceptions as e:
    pass
trimmomatic_cmd = trimmomatic_cmd_str %(fq_1, fq_2, fq_1_trimmed, fq_1_unpaired, fq_2_trimmed, fq_2_unpaired) 
print(trimmomatic_cmd)
os.system(trimmomatic_cmd)
os.environ['fq_1_trimmed'] = fq_1_trimmed
os.environ['fq_2_trimmed'] = fq_2_trimmed
print(fq_1_trimmed, fq_2_trimmed)
os.system("fastqc -t 2 -o ${out_path}/1-FASTQC/ ${fq_1_trimmed} ${fq_2_trimmed} 2> ${out_path}/0-LOG/${SM}.${lane}.fastqc.log 1>&2")

print >> logfile,'fastqc finished at ' + time.ctime()
error_log_fq.close()
LB = line.split(":")[2]
ID = LB + '.' + line.split(":")[3]
PU = ID + '.' + conf['SM']
os.environ['LB'] = LB
os.environ['ID'] = ID
os.environ['PU'] = PU
os.system("bwa mem -t 2 -R \'@RG\\tID:"+conf['SM']+"\\tSM:"+conf['SM']+"\\tPU:"+PU+"\\tLB:"+LB+"\\tPL:"+conf['PL']+"\' -M ${hgRef}/hg19 ${fq_1_trimmed} ${fq_2_trimmed} 2>${bwa_log} | samtools view -bS -o ${out_path}/2-MAPPING/${SM}.${lane}.bam - 2>${bwa_log} 1>&2")
bwa_error_list=['ERROR','error','Error','truncated file']
filter_bwa=filter_log(bwa_error_list,'[main] Real time:',conf['bwa_log'])
if filter_bwa:
    print >> logfile,'bwa finished at ' + time.ctime()
else:
    print >> error_log_bwa,conf['SM']
    #send_mail.send_log_mail('Wrong in bwa','Wrong_bwa_'+conf['SM'])
    os._exit(1)
os.system("samtools sort ${out_path}/2-MAPPING/${SM}.${lane}.bam -T /opt/tmp/${SM}.${lane} -o ${out_path}/2-MAPPING/${SM}.${lane}.sorted.bam 2> ${out_path}/0-LOG/${SM}.${lane}.bamsort.log 1>&2")
print >> logfile,'bwa_sort finished at ' + time.ctime()
error_log_bwa.close()

os.system("java -Xmx5g -jar ${PICARDJAR} MarkDuplicates VALIDATION_STRINGENCY=SILENT INPUT=${out_path}/2-MAPPING/${SM}.${lane}.sorted.bam OUTPUT=${out_path}/3-PICARD/${SM}.${lane}.marked.bam METRICS_FILE=${out_path}/3-PICARD/${SM}.${lane}.Mkdup.metrics 2> ${markduplicates_log} 1>&2")
filter_markduplicate=filter_log(['ERROR','problem','invalid'],'picard.sam.markduplicates.MarkDuplicates done',conf['markduplicates_log'])
if filter_markduplicate:
    print >> logfile,'MarkDuplicates finished at ' + time.ctime()
else:
    print >> error_log_markduplicate,conf['SM']
    #send_mail.send_log_mail('Wrong MarkDuplicates','Wrong_MarkDuplicates_'+conf['SM'])
    os._exit(1)
error_log_markduplicate.close()

os.system("samtools index ${out_path}/3-PICARD/${SM}.${lane}.marked.bam 2> ${indexmarkduplicates_log} 1>&2")
print >> logfile,'MarkDuplicates_samtools_index finished at ' + time.ctime()

os.system("java -Xmx5g -jar ${GATKJAR} -T RealignerTargetCreator -U ALLOW_N_CIGAR_READS -R ${hgRef}/hg19.fa -I ${out_path}/3-PICARD/${SM}.${lane}.marked.bam -o ${out_path}/4-GATK/${lane}/${SM}.${lane}.marked.intervals -known ${known1000G_indels} -known ${GoldStandard_indels} 2> ${realignertargetcreator_log} 1>&2")
filter_realignertargetcreator=filter_log(['ERROR','problem','invalid'],'ProgressMeter - Total runtime',conf['realignertargetcreator_log'])
if filter_realignertargetcreator:
    print >> logfile,'GATK_RealignerTargetCreator finished at ' + time.ctime()
else:
    print >> error_log_realignertargetcreator,conf['SM']
    #send_mail.send_log_mail('Wrong GATK RealignerTargetCreator','Wrong_GATK_RealignerTargetCreator_'+conf['SM'])
    os._exit(1)
error_log_realignertargetcreator.close()

os.system("java -Xmx5g -jar ${GATKJAR} -T IndelRealigner -R ${hgRef}/hg19.fa -I ${out_path}/3-PICARD/${SM}.${lane}.marked.bam -targetIntervals ${out_path}/4-GATK/${lane}/${SM}.${lane}.marked.intervals -o ${out_path}/4-GATK/${lane}/${SM}.${lane}.marked.realn.bam -known ${known1000G_indels} -known ${GoldStandard_indels} 2> ${indelrealigner_log} 1>&2")
filter_indelrealigner=filter_log(['ERROR','problem','invalid'],'ProgressMeter - Total runtime',conf['indelrealigner_log'])
if filter_indelrealigner:
    print >> logfile,'GATK_IndelRealigner finished at ' + time.ctime()
else:
    print >> error_log_indelrealigner,conf['SM']
    #send_mail.send_log_mail('Wrong GATK IndelRealigner','Wrong_GATK_IndelRealigner_'+conf['SM'])
    os._exit(1)
error_log_indelrealigner.close()

os.system("java -Xmx5g -jar ${GATKJAR} -nct 2 -T BaseRecalibrator -R ${hgRef}/hg19.fa -I ${out_path}/4-GATK/${lane}/${SM}.${lane}.marked.realn.bam -o ${out_path}/4-GATK/${lane}/${SM}.${lane}.marked.realn.recal -knownSites ${known1000G_indels} -knownSites ${GoldStandard_indels} -knownSites ${dbSNP} 2> ${baserecalibrator_log} 1>&2")
filter_baserecalibrator=filter_log(['ERROR','problem','invalid'],'ProgressMeter - Total runtime',conf['baserecalibrator_log'])
if filter_baserecalibrator:
    print >> logfile,'GATK_BaseRecalibrator finished at ' + time.ctime()
else:
    print >> error_log_baserecalibrator,conf['SM']
    #send_mail.send_log_mail('Wrong GATK BaseRecalibrator','Wrong_GATK_BaseRecalibrator_'+conf['SM'])
    os._exit(1)
error_log_baserecalibrator.close()

os.system("java -Xmx5g -jar ${GATKJAR} -nct 2 -T PrintReads -R ${hgRef}/hg19.fa -I ${out_path}/4-GATK/${lane}/${SM}.${lane}.marked.realn.bam -o ${out_path}/4-GATK/${lane}/${SM}.${lane}.marked.realn.recal.bam --BQSR ${out_path}/4-GATK/${lane}/${SM}.${lane}.marked.realn.recal  2> ${printreads_log} 1>&2")
filter_printreads=filter_log(['ERROR','problem','invalid'],'ProgressMeter - Total runtime',conf['printreads_log'])
if filter_printreads:
    print >> logfile,'GATK_PrintReads finished at ' + time.ctime()
else:
    print >> error_log_printreads,conf['SM']
    #send_mail.send_log_mail('Wrong GATK PrintReads','Wrong_GATK_PrintReads_'+conf['SM'])
    os._exit(1)
error_log_printreads.close()

try:
    os.mknod(conf['out_path'] + '/4-GATK/' + conf['lane'] + "/" + conf['SM'] + '_finished')
except OSError as exc:
    pass
os.system("samtools view -b ${out_path}/4-GATK/${lane}/${SM}.${lane}.marked.realn.recal.bam chrY > ${out_path}/4-GATK/${lane}/${SM}.${lane}.chrY.bam")

os.system("java -Xmx5g -jar ${GATKJAR} -T HaplotypeCaller -R ${hgRef}/hg19.fa -L ${IntervalList} --dbsnp ${dbSNP}  -I ${out_path}/4-GATK/${lane}/${SM}.${lane}.marked.realn.recal.bam -o ${out_path}/5-GERMLINE/${SM}.${lane}.germline.vcf 2> ${haplotypecaller_log}")
print "okay"
os.system("${ANNOVAR}/table_annovar.pl ${out_path}/5-GERMLINE/${SM}.${lane}.germline.vcf ${HUMANDB} --buildver hg19 --outfile ${out_path}/5-GERMLINE/${SM}.${lane}.germline --remove --protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp150,dbnsfp33a,clinvar_20170905,exac03,dbscsnv11,dbnsfp31a_interpro,rmsk,ensGene,knownGene  --operation  g,f,f,f,f,f,f,f,f,r,g,g --nastring . -vcfinput")

os.system("python2 ${INTERVAR} -b hg19 -i ${out_path}/5-GERMLINE/${SM}.${lane}.germline.vcf --input_type=VCF -o ${out_path}/5-GERMLINE/${SM}.${lane}.germline --skip_annovar")
filter_printreads=filter_log(['error''ERROR','problem','invalid'],'Total runtime',conf['haplotypecaller_log'])
if filter_printreads:
    print >> logfile,'GATK_HaplotypeCaller finished at ' + time.ctime()
else:
    print >> error_log_printreads,conf['SM']
    #send_mail.send_log_mail('Wrong GATK HaplotypeCaller','Wrong_GATK_PrintReads_'+conf['SM'])
    os._exit(1)
error_log_haplotypecaller.close()

os.system("java -Xmx5g -jar ${PICARDJAR} CalculateHsMetrics I=${out_path}/4-GATK/${lane}/${SM}.${lane}.marked.realn.recal.bam O=${out_path}/6-HsMetrics/${SM}.${lane}.marked.realn.recal.bam.metrics BAIT_INTERVALS=${IntervalList} TARGET_INTERVALS=${IntervalList} 2> ${calculatehsmetrics_log} 1>&2")
os.system("./Metrics_statistic.sh -D ${out_path}/6-HsMetrics/ -M ${SM}.${lane}.marked.realn.recal.bam.metrics -O ${out_path}/6-HsMetrics/ -P ${SM}")
os.system("python ./VcfStatistic.py ${out_path}/5-GERMLINE/${SM}.${lane}.germline.hg19_multianno.txt ${out_path}/7-QC_Stat/${SM}.${lane}.stat.txt")
metrics = os.path.join(out_path.rstrip('/')+'/6-HsMetrics/',conf['SM']+'.Formated')
os.environ["metrics"] = metrics
if os.path.isfile(metrics):
    os.system("cp ${metrics} ${out_path}/7-QC_Stat/") 
filter_printreads=filter_log(['error''ERROR','problem','invalid'],'picard.analysis.directed.CalculateHsMetrics done', conf['calculatehsmetrics_log'])
if filter_printreads:
    print >> logfile,'GATK_CalculateHsMetrics finished at ' + time.ctime()
else:
    print >> error_log_calculatehsmetrics,conf['SM']
    #send_mail.send_log_mail('Wrong GATK CalculateHsMetrics','Wrong_GATK_PrintReads_'+conf['SM'])
    os._exit(1)
error_log_haplotypecaller.close()

print >> logfile,'End at ' + time.ctime()
logfile.close()

log_finished = open(conf['out_path']+'/0-LOG/GATK_finished.log','a')
print >> log_finished,conf['SM']
log_finished.close()
