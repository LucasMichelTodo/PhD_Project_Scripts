import pybedtools as pb
import numpy as np
import pandas as pd
import os
from collections import defaultdict

indir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/New_Coverage/Binsize_150/"

#bin_bed = pb.BedTool("/mnt/Disc4T/Projects/PhD_Project/Data/PlasmoDB-45_Pfalciparum3D7_bin5_2prevGenes.bed")

bin_bed = pb.BedTool("/mnt/Disc4T/Projects/PhD_Project/Data/PlasmoDB-46_Pfalciparum3D7_withGDV1_ncRNAs_bin5_2prevGenes.bed")

cov_files = [f for f in os.listdir(indir) if f.endswith("_normInp_ps1.bdg.gz")]

#cov_files = ['NF54_me_q5_sort_noDup_RPKM_normInp_ps1.bdg.gz']

for cov in cov_files:

    flstr = indir+cov
    cov = pb.BedTool(flstr)
    print(flstr, "Converted to bed!")
    outfile = flstr.replace(".bdg.gz", "_5binned_cov_2prevpost.csv")
    #logfile = outfile.replace(".bed", ".log")
    genevals = defaultdict(list)

    for interval in bin_bed:

        gene = interval.name
        pos = interval.score

        # Not all of them have a match!
        try:
            match = cov.tabix_intervals(interval)
            val = np.mean([float(x.fields[3]) for x in match])
            genevals[gene].append((val, pos))

        except:
            pass

    # Rearrange values deppending on strandness (we have to "flip" genes on "-" strand)
    sorted_genevals = {}
    for gene, val in genevals.items():
        svals = sorted(val, key = lambda x:int(x[1]))
        vals = [x[0] for x in svals]
        sorted_genevals[gene] = vals

    # Write output
    df = pd.DataFrame.from_dict(sorted_genevals, orient='index')
    df.to_csv(outfile)
    print("Done with file: {}" .format(flstr))

import numpy as np
import pandas as pd
import os
import pybedtools as pb
from collections import defaultdict

wd = '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_normInput/'
os.chdir(wd)

indir = './'
bin_bed = pb.BedTool("/mnt/Disc4T/Projects/PhD_Project/Data/PlasmoDB-46_Pfalciparum3D7_withGDV1_ncRNAs_bin5_2prevGenes.bed")

cov_files = [f for f in os.listdir(indir) if f.endswith(".bdg.gz")]
#cov = cov_files[0]

for cov in cov_files:

    flstr = indir+cov
    cov = pb.BedTool(flstr)
    print(flstr, "Converted to bed!")
    outfile = flstr.replace(".bdg.gz", "_5binned_cov_2prevpost.csv")
    #logfile = outfile.replace(".bed", ".log")
    genevals = defaultdict(list)

    for interval in bin_bed:

        #interval = bin_bed[0]
        gene = interval.name
        pos = interval.score

        # Not all of them have a match!
        try:
            match = cov.tabix_intervals(interval)
            val = np.mean([float(x.fields[3]) for x in match])
            genevals[gene].append((val, pos))

        except:
            pass

    # Rearrange values deppending on strandness (we have to "flip" genes on "-" strand)
    sorted_genevals = {}
    for gene, val in genevals.items():
        svals = sorted(val, key = lambda x:int(x[1]))
        vals = [x[0] for x in svals]
        sorted_genevals[gene] = vals

    # Write output
    df = pd.DataFrame.from_dict(sorted_genevals, orient='index')
    df.to_csv(outfile)
    print("Done with file: {}" .format(flstr))
