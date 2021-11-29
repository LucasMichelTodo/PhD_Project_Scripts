import os
import pandas as pd

wd = "/mnt/Disc4T/Projects/PhD_Project/Variant_Calling/"
os.chdir(wd)

var_df = pd.read_csv("parsed_variants.csv", index_col=0)
var_df.head()

# Drop empty columns
empty = var_df.columns[var_df.isna().all(axis=0)].tolist()
var_df = var_df.drop(columns=empty)
var_df.columns

# Set ratio difference and read depth thresholds
r_thld = 0.5
d_thld = 20

# Subset DF
NF54_p18 = abs(var_df["RefRatio_NF54"] - var_df["RefRatio_P18"]) >= r_thld
NF54_p63 = abs(var_df["RefRatio_NF54"] - var_df["RefRatio_P63"]) >= r_thld
p18_p63 = abs(var_df["RefRatio_P63"] - var_df["RefRatio_P18"]) >= r_thld

d_NF54 = var_df["depth_NF54"] >= d_thld
d_p18 = var_df["depth_P18"] >= d_thld
d_p63 = var_df["depth_P63"] >= d_thld

# df = var_df[(NF54_p18 | NF54_p63 | p18_p63) &
#                    (d_NF54 & d_p18 & d_p63)]

df = var_df[(NF54_p18 & d_NF54 & d_p18) |
            (NF54_p63 & d_NF54 & d_p63) |
            (p18_p63 & d_p18 & d_p63)]

vep_cols = ['Consequence',
            'IMPACT',
            'Feature_type',
            'Feature',
            'BIOTYPE',
            'EXON',
            'INTRON',
            'cDNA_position',
            'CDS_position',
            'Protein_position',
            'Amino_acids',
            'Codons',
            'DISTANCE',
            'STRAND']

df = df[['Ref', 'Alt', 'Chrom', 'Pos',
         'RefRatio_NF54', 'RefRatio_P18', 'RefRatio_P63',
         'depth_NF54', 'depth_P18', 'depth_P63',
         'Gene', 'Annot']+vep_cols]

df.to_csv("parsed_epireset_annotated_newfilter.csv")
