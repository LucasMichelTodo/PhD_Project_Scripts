import subprocess as sp
import os

def normCoverage(bam, binsize, process, suffix, inpath="./", outpath="./"):

    inbam = inpath+bam
    out = outpath+bam.rsplit(".", 1)[0]+suffix

    cmd = ["bamCoverage",
           "--bam", inbam,
           "--outFileName", out,
           "--binSize", str(binsize),
           "--numberOfProcessors", str(process),
           "--extendReads",
           "--outFileFormat", "bigwig",
           "--normalizeUsing", "RPKM"]

    sp.run(cmd,
           stdout=open(outpath+'normCov_log.txt', 'w'))


def normbyInput(b1, b2, ps, binsize, p, suffix,  inpath="./", outpath="./"):

    out = outpath+b1.rsplit(".", 1)[0]+suffix
    b1 = inpath+b1
    b2 = inpath+b2

    cmd = ["bigwigCompare",
           "-b1", b1,
           "-b2", b2,
           "--pseudocount", str(ps),
           "-bs", str(binsize),
           "-p", str(p),
           "-o", out,
           "-of", "bedgraph"]

    sp.run(cmd,
           stdout=open(outpath+'normInp_log.txt', 'w'))

# inpath = "/media/lucas/Disc4T/Projects/Cristina_ChIP_All/New_Coverage/NoDupBams/"
# outpath = "/media/lucas/Disc4T/Projects/Cristina_ChIP_All/New_Coverage/"

# Set filenames
# bams = ["1.2B_me_sort_q5_noDup.bam", "1.2B_in_sort_q5_noDup.bam",
#         "10G_me_sort_q5_noDup.bam", "10G_in_sort_q5_noDup.bam",
#         "B11_me_sort_q5_noDup.bam", "B11_in_sort_q5_noDup.bam",
#         "A7K9_me_sort_q5_noDup.bam", "A7K9_in_sort_q5_noDup.bam",
#         "E5K9_me_sort_q5_noDup.bam", "E5K9_in_sort_q5_noDup.bam"]

# bwigs = [bam.replace(".bam", "_RPKM.bw") for bam in bams]

inpath = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/Current/"
outpath = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/Current/"

bams = [b for b in os.listdir(inpath) if b.endswith('_q5_sort_noDup.bam')]
bwigs = [('i10-E5HAme_q5_sort_noDup_RPKM.bw', 'i2-E5HAin_q5_sort_noDup_RPKM.bw'),
         ('i6-E5HAac_q5_sort_noDup_RPKM.bw', 'i2-E5HAin_q5_sort_noDup_RPKM.bw')]


# Calls
for bam in bams:
    sp.call(f"samtools index {inpath+bam}", shell=True)
    normCoverage(bam, 150, 8, "_RPKM.bw", inpath=inpath, outpath=outpath)

for x in bwigs:
    normbyInput(x[0], x[1], 1, 150, 8, "_normInp_ps1.bdg", inpath, outpath)
