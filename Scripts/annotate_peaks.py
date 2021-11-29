import pybedtools as pb
import pandas as pd
import os

## Function

def annotate_bed(bedfile):

    bed = pb.BedTool(bedfile)
    anot = bed.intersect(ref, wao=True)

    parsed_anot = []
    for interval in anot:
        originalfields = interval.fields[0:5]
        if interval.fields[5] == ".":
            gid, anot = "intergenic", "NA"
        else:
            gid = interval.fields[8]
            anot = interval.fields[9]
        parsed_anot.append(originalfields + [gid, anot])

    df = pd.DataFrame(parsed_anot)
    outfile = bedfile.replace(".bed", "_annotated.csv")
    df.to_csv(outfile, sep="\t", header=False, index=False)

## Calls

gff = "/mnt/Disc4T/Projects/PhD_Project/Binned_Coverage/final_binned_1000tss_wholegene_parsed.bed"
peaksdir = "/mnt/Disc4T/Projects/Chip_Seq_Data_2021/BetaFit_DuplDel_Filtered/"
peaks = [f for f in os.listdir(peaksdir) if f.endswith("_filtered.bed")]

ref = pb.BedTool(gff)
ref = ref.sort()

for peak in peaks:
    bedfile = peaksdir+peak
    annotate_bed(bedfile)
