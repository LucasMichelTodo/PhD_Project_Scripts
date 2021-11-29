import pybedtools as pb
import pandas as pd

gff = pb.BedTool("/mnt/Disc4T/Projects/PhD_Project/Data/PlasmoDB-46_Pfalciparum3D7_withGDV1_ncRNAs.gff")
gff = gff.filter(lambda x: x[2] == "gene")

df = []
for feat in gff:
    df.append([feat.attrs['ID'], feat.attrs['description']])

anot = pd.DataFrame(df)
anot.to_csv("/mnt/Disc4T/Projects/PhD_Project/Data/r_annot.csv", index=False, header=False)
