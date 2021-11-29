import pybedtools as pb
import os

indir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/New_Coverage/Binsize_150/"
fls = os.listdir(indir)
beds = [indir+bdg for bdg in fls if bdg.endswith(".bdg")]

for bed in beds:
    cov = pb.BedTool(bed)
    cov = cov.tabix()
