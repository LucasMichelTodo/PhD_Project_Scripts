import subprocess as sp
import os

gatk = '/home/lucas/Programs/gatk-4.1.9.0/gatk'
wd = '/mnt/Disc4T/Projects/PhD_Project/Variant_Calling/'
os.chdir(wd)

def AddOrReplaceReadGroups(bam):
    output = bam.replace('.bam', '_withRG.bam')
    name = bam.rsplit('.')[0]
    cmd = ("java -jar /home/lucas/Programs/picard.jar "
           "AddOrReplaceReadGroups "
           "-INPUT {} "
           "-OUTPUT {} "
           "-RGID group_{} "
           "-RGLB lib_{} "
           "-RGPL illumina "
           "-RGPU unit1 "
           "-RGSM {}_sample") .format(bam, output, name, name, name)
    sp.call(cmd, shell=True)

def mark_duplicates(bam):
    outfile = bam.replace('.bam', '_markedDuplicates.bam')
    mtrcsfile = bam.replace('.bam', '_metrics.txt')
    args = [gatk, bam, outfile, mtrcsfile]
    cmd = '{} MarkDuplicates -I {} -O {} -M {}' .format(*args)
    sp.call(cmd, shell=True)

    cmd = 'samtools sort {} -o {}' .format(outfile, outfile)
    sp.call(cmd, shell=True)

def base_recalibration(bam):

    outfile = bam.replace('.bam', '_baserecal_table.table')
    known = './empty.bed'
    ref = './ref.fasta'
    args = [gatk, bam, ref, known, outfile]

    cmd = ('{} BaseRecalibrator '
           '-I {} -R {} '
           '--known-sites {} '
           '-O {}') .format(*args)

    sp.call(cmd, shell=True)

def applyBQSR(bam):

    outfile = bam.replace('.bam', '_BQSR.bam')
    recal_table = bam.replace('.bam', '_baserecal_table.table')
    ref = './ref.fasta'
    args = [gatk, bam, ref, recal_table, outfile]

    cmd = ('{} ApplyBQSR '
           '-I {} -R {} '
           '--bqsr-recal-file {} '
           '-O {}') .format(*args)

    sp.call(cmd, shell=True)

def mergeBams(*bams, out):
    nbams = len(bams)
    inputs = '-I {} '*nbams
    cmd = 'java -jar /home/lucas/Programs/picard.jar ' \
        'MergeSamFiles ' + \
        inputs .format(*bams) + \
        '-O {}.bam' .format(out)
    sp.call(cmd, shell=True)

def call_variants(bam):
    outfile = bam.replace('.bam', '_variants.vcf')
    ref = './ref.fasta'
    args = [gatk, ref, bam, outfile]

    cmd = ('{} --java-options "-Xmx4g" HaplotypeCaller '
           '-R {} -I {} -O {} -ploidy 1') .format(*args)

    sp.call(cmd, shell=True)

def call_VEP(vcf, gff, fasta):

    out = vcf.replace('.vcf', '_VEPannotated.txt')
    args = [vcf, out, gff, fasta]

    cmd = ("/home/lucas/Programs/ensembl-vep/vep "
           "-i {} "
           "-o {} "
           "--gff {} "
           "--fasta {} "
           "--force_overwrite "
           "--vcf") .format(*args)

    sp.call(cmd, shell=True)

import os
import subprocess as sp

wd = '/mnt/Disc4T/Projects/PhD_Project/Variant_Calling/'
os.chdir(wd)

gatk = '/home/lucas/Programs/gatk-4.1.9.0/gatk'
indir = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/Bams/'

os.listdir(indir)

bams = ['1.2B_in_sort_q5.bam',
        '10G_in_sort_q5.bam',
        'A7K9_in_sort_q5.bam',
        'E5K9_in_sort_q5.bam',
        'B11_in_sort_q5.bam'
        #'NF54_in_renamed_q5_sort.bam',
        ]

for bam in bams:
    bam = indir+bam

    AddOrReplaceReadGroups(bam)
    bam = bam.replace('.bam', '_withRG.bam')

    mark_duplicates(bam)
    bam = bam.replace('.bam', '_markedDuplicates.bam')

    base_recalibration(bam)
    applyBQSR(bam)
    bam = bam.replace('.bam', '_BQSR.bam')


bamlist = [f for f in os.listdir(indir) if f.endswith('_withRG_markedDuplicates_BQSR.bam')]

os.chdir(indir)
mergeBams(*bamlist, out = 'merged_12B_10G_A7_E5_B11')

bam = 'merged_12B_10G_A7_E5_B11.bam'
sp.call('samtools index {}' .format(bam), shell=True)
call_variants(bam)

os.chdir(wd)
vcf = './merged_12B_10G_A7_E5_B11_variants.vcf'
gff = './PlDB-52_Pfalciparum3D7_vep_changetypes.gff.gz'
fasta = './ref.fasta'
call_VEP(vcf, gff, fasta)
