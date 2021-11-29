import subprocess as sp
import os

## Functions

def call_BBDUK(in1, in2, out1, out2, outm, ref, params):
    cmd = ("bbduk.sh in={} in2={} "
           "out={} out2={} outm={} "
           "ref={}") .format(in1, in2, out1, out2, outm, ref)

    cmd = cmd+" "+params
    sp.call(cmd, shell=True)

## Calls

params = "ktrim=r k=22 mink=6 overwrite=t"
ref = "/home/lucas/Programs/bbmap/resources/adapters.fa"

root_path = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Raw_Data/"
read1s = []
read2s = []
for path, subdirs, files in os.walk(root_path):
    for f in files:
        if all(x in f for x in ["read1", ".fastq"]):
            read1s.append(os.path.join(path, f))
        elif all(x in f for x in ["read2", ".fastq"]):
            read2s.append(os.path.join(path, f))
        else:
            print(f)


read1s = sorted(read1s)
read2s = sorted(read2s)
outpath = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Clean_Reads/"

for pair in zip(read1s, read2s):
    in1, in2 = pair[0], pair[1]
    out1 = outpath + pair[0].rsplit("/", 1)[1].replace(".fastq", "_clean.fastq")
    out2 = outpath + pair[1].rsplit("/", 1)[1].replace(".fastq", "_clean.fastq")
    outm = outpath + pair[0].rsplit("/", 1)[1].replace(".fastq", "_badreads.fastq")
    call_BBDUK(in1, in2, out1, out2, outm, ref, params)

## Funtions

def call_Bowtie2(in1, in2, out, params):
    cmd = "bowtie2 -1 {} -2 {} -S {}" .format(in1, in2, out)
    cmd = cmd+" "+params
    print(cmd)
    sp.call(cmd, shell=True)

## Calls
params = ("-p 4 --very-sensitive --local "
          "-5 4 -3 4 -I 50 -X 200 "
          "-x /home/lucas/Programs/bowtie2-2.3.0-legacy/Pf3D7")

inpath = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/Clean_Reads/"
# files = [f for f in os.listdir(inpath)if f.endswith("_clean.fastq.gz")]

# read1s = sorted([f for f in files if "read1" in f])
# read2s = sorted([f for f in files if "read2" in f])

read1s = ['i10-E5HAme_26899_TAGCTT_read1_clean.fastq.gz',
          'i2-E5HAin_26891_CGATGT_read1_clean.fastq.gz',
          'i6-E5HAac_26895_GCCAAT_read1_clean.fastq.gz']

read2s = ['i10-E5HAme_26899_TAGCTT_read2_clean.fastq.gz',
          'i2-E5HAin_26891_CGATGT_read2_clean.fastq.gz',
          'i6-E5HAac_26895_GCCAAT_read2_clean.fastq.gz']

outpath = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/Current/"

for pair in zip(read1s, read2s):
    name = pair[0].split("_")[0]

    in1, in2 = inpath+pair[0], inpath+pair[1]
    #out = inpath+name+".sam"
    out = outpath+name+".sam"
    call_Bowtie2(in1, in2, out, params)

import subprocess as sp
import sys
from tqdm import tqdm

## Functions

def from_sam_to_bam(samfile):
    name = samfile.rsplit(".")[0]
    cmd = "samtools view -bS {} > {}" .format(samfile, name+".bam")
    sp.call(cmd, shell=True)

    ### Erase SAM after creating BAM
    # cmd = "rm {}" .format(samfile)
    # sp.call(cmd, shell=True)

    cmd = "samtools sort {} > {}" .format(name+".bam", name+"_sort.bam")
    sp.call(cmd, shell=True)

    ### Erase bam after creating sortedBAM
    cmd = "rm {}" .format(name+".bam")
    sp.call(cmd, shell=True)

    cmd = "samtools index {} > {}" .format(name+"_sort.bam", name+"_sort.bam.bai")
    sp.call(cmd, shell=True)

    ## Filter only >=q5 reads
    cmd = "samtools view -b -q 5 {} > {}" .format(name+"_sort.bam", name+"_q5_sort.bam")
    sp.call(cmd, shell=True)

## Calls

#indir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/Bams/"
indir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/New_Bams/"
samfiles = [f for f in os.listdir(indir) if f.endswith(".sam")]

for f in tqdm(samfiles):
    sam = indir+f
    from_sam_to_bam(sam)

import os

def remove_duplicates(indir, outdir, bam):

    i = indir+bam
    o = outdir+bam.replace(".bam", "_noDup.bam")
    m = outdir+bam.replace(".bam", "_metrics.txt")

    cmd = ("java -jar /home/lucas/Programs/picard.jar MarkDuplicates "
           "REMOVE_DUPLICATES=true I={} O={} M={}") .format(i, o, m)

    sp.call(cmd, shell=True)

# indir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/Bams/"
# bams = sorted([b for b in os.listdir(indir) if b.endswith(".bam")])
# outdir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/New_Coverage/NoDupBams/"


indir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/Current/"
bams = sorted([b for b in os.listdir(indir) if b.endswith("q5_sort.bam")])
outdir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/Current/"

for bam in bams:
    remove_duplicates(indir, outdir, bam)
