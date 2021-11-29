import os

indir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/New_Coverage/New_Peaks/"
macs2_files = [f for f in os.listdir(indir) if f.endswith(".narrowPeak")]
macs2_files = ['NF54_Macspeaks_peaks.narrowPeak', 'E5K9_Macspeaks_peaks.narrowPeak']

def narrowPeak_to_bed(macs2file):
    with open(macs2file, "r+") as infile:
        outfile = macs2file.replace(".narrowPeak", ".bed")
        with open(outfile, "w+") as output:
            for line in infile:
                linelist = line.strip().split("\t")
                bedlist = linelist[0:3]+[linelist[6]]
                output.write("\t".join(bedlist)+"\n")

for f in macs2_files:
    fl = indir+f
    narrowPeak_to_bed(fl)
