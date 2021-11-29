#### Functions ####
import subprocess
import os

## Functions
def get_RPKMs(bam, bs, smooth, out):

    cmd = ('bamCoverage -b {} '
       '--outFileFormat bedgraph '
       '--normalizeUsing RPKM '
       '-p 8 '
       '-bs {} '
       '--smoothLength {} '
       '-o {}')

    print(cmd)
    subprocess.call(cmd .format(bam, bs, smooth, out), shell=True)

indir = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/Bams/'
outdir = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_1bp/' 

bams = [
    'E5K9_in_sort_q5.bam',
    'E5K9_me_sort_q5.bam',
    'A7K9_me_sort_q5.bam',
    'A7K9_in_sort_q5.bam',
    'NF54_me_renamed_q5_sort.bam',
    'NF54_in_renamed_q5_sort.bam',
]

for bam in bams:
    infile = indir + bam
    outfile = outdir + bam.replace('.bam', '_RPKMs_1bp_smooth150.bdg')
    get_RPKMs(infile, 1, 150, outfile)


    
