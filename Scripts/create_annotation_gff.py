# Import packages
import pybedtools as py
from pybedtools.featurefuncs import TSS
from pybedtools.featurefuncs import three_prime
import subprocess as sp

# Parameters to set
upstream = 500
dwstream = 500

ref = '/mnt/Disc4T/Projects/PhD_Project/Data/PlasmoDB-46_Pfalciparum3D7_withGDV1_ncRNAs.gff'
out_name = '/mnt/Disc4T/Projects/PhD_Project/Data/PlDB-46_Pf3D7_GDV1_ncRNAs_500bp_bothends.gff'

# Create genome dict
genome={
    'Pf3D7_01_v3': 640851,
    'Pf3D7_02_v3': 947102,
    'Pf3D7_03_v3': 1067971,
    'Pf3D7_04_v3': 1200490,
    'Pf3D7_05_v3': 1343557,
    'Pf3D7_06_v3': 1418242,
    'Pf3D7_07_v3': 1445207,
    'Pf3D7_08_v3': 1472805,
    'Pf3D7_09_v3': 1541735,
    'Pf3D7_10_v3': 1687656,
    'Pf3D7_11_v3': 2038340,
    'Pf3D7_12_v3': 2271494,
    'Pf3D7_13_v3': 2925236,
    'Pf3D7_14_v3': 3291936,
    'Pf3D7_API_v3': 34250,
    'Pf_M76611': 5967
}


# Filter GFF to retain only "gene" entries (we don't want to anotate each peak many times).
gff = py.BedTool(ref)
sel = ["gene"]
gene_gff = gff.filter(lambda x: x.fields[2] in sel)

def elongate_feat(feat, up, down):
    if feat.start >= up:
        feat.start = feat.start-up
    else:
        feat.start = 0

    feat.stop = feat.stop+down
    if feat.stop > genome[feat.chrom]:
        feat.stop = genome[feat.chrom]
    return(feat)

out = gff.each(elongate_feat, up=upstream, down=dwstream)
out.saveas(out_name)
