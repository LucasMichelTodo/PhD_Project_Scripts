import subprocess
import os
import itertools

## Functions

def getName(filename):
    outname = filename.rsplit("/", 1)[1] #Remove path
    out = outname.split("_", 1)[0] #Remove added names
    return(out)


def getDepth(peaksfile):
    dlist = []
    with open(peaksfile) as f:
        for line in f:
            if line.startswith("# total fragments"):
                d = line.split(":")[1].strip()
                dlist.append(int(d))
    return(str(min(dlist)))


def macs2callpeak(t, c, params):

    # Make sure we are using apropiate MACS2 version (2.1.2)
    print("You are using MACS version:")
    cmd = "macs2 --version"
    subprocess.call(cmd, shell=True)
    print("\n")

    cmd = ("macs2 callpeak -t {} -c {} ") .format(t, c) + params

    subprocess.call(cmd, shell=True)


def macs2DifPeaks(t1, c1, t2, c2):

    # Make sure we are using apropiate MACS2 version (2.1.2)
    print("You are using MACS version:")
    cmd = "macs2 --version"
    subprocess.call(cmd, shell=True)
    print("\n")

    t1pile = getName(t1) + "_treat_pileup.bdg"
    c1pile = getName(c1) + "_control_lambda.bdg"

    t2pile = getName(t2) + "_treat_pileup.bdg"
    c2pile = getName(c2) + "_control_lambda.bdg"

    peaks1 = getName(t1) + "_peaks.xls"
    peaks2 = getName(t2) + "_peaks.xls"

    d1 = getDepth(peaks1)
    d2 = getDepth(peaks2)

    prefix = getName(t1)+"_vs_"+getName(t2)+"_g140_l150_c50"

    cmd = ("macs2 bdgdiff "
           "--t1 {} --c1 {} "
           "--t2 {} --c2 {} "
           "--d1 {} --d2 {} "
           "--o-prefix {} "
           "-g 140 -l 150 --cutoff 50") .format(t1pile, c1pile,
                                                t2pile, c2pile,
                                                d1, d2, prefix)

    subprocess.call(cmd, shell=True)

def manormDifPeaks(peaks1, peaks2, reads1, reads2, params):

    #Check which MaNorm version we are using:
    print("\n")
    print("You are using MANorm version:")
    cmd = "manorm --version"
    subprocess.call(cmd, shell=True)
    print("\n")

    cmd = ("manorm --p1 {} --p2 {} --r1 {} --r2 {} ") .format(peaks1, peaks2,
                                                              reads1, reads2)
    cmd = cmd + params

    subprocess.call(cmd, shell=True)

indir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/New_Coverage/NoDupBams/"
outdir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/New_Coverage/New_Peaks/"

# Select chip and input tracks from a given folder:
def is_chip(x):
    isbam = x.endswith("_noDup.bam")
    ismet = x.split("_")[1] == "me"
    return(isbam and ismet)

def is_input(x):
    isbam = x.endswith("_noDup.bam")
    isin = x.split("_")[1] == "in"
    return(isbam and isin)

chips = [b for b in os.listdir(indir) if is_chip(b)]
inputs = [b for b in os.listdir(indir) if is_input(b)]

chips = [b for b in chips if "NF54" in b]
inputs = [b for b in inputs if "NF54" in b]


# Ensure same ordering (by name)
chips = sorted(chips)
inputs = sorted(inputs)

# Call PeakCalling
params_form = ("-f BAMPE -B "
               "-g 2.41e7 "
               "--keep-dup all "
               "--fe-cutoff 1.5 "
               "--nomodel "
               "--extsize 150 "
               "-n {} "
               f"--outdir {outdir}")

for pair in zip(chips, inputs):

    t = indir + pair[0]
    c = indir + pair[1]
    name = pair[0].split("_")[0]+"_Macspeaks"
    params = params_form .format(name)

    macs2callpeak(t, c, params)

    print("==============================")
    print("Finished {}!" .format(name))
    print("==============================\n\n\n")

import os
import itertools

## Original Call
peaksdir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/New_Coverage/New_Peaks/"
peakfiles = os.listdir(peaksdir)
peakfiles = sorted([p for p in peakfiles if p.endswith(".narrowPeak")])

readsdir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/New_Coverage/NoDupBams/"
readfiles = os.listdir(readsdir)
readfiles = [b for b in readfiles if b.endswith(".bam")]
readfiles = sorted([b for b in readfiles if b.split("_")[1] == "me"])

outdir = "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/New_Coverage/MaNormPeaks"
pairs = zip(peakfiles, readfiles)
difs = list(itertools.combinations(pairs, 2))


# difs_NF54 = []
# for d in difs:
#     if any("NF54" in fl for fl in [d[0][0],
#                                    d[0][1],
#                                    d[1][0],
#                                    d[1][1]]):
#         difs_NF54.append(d)


for dif in difs_NF54:

    p1 = peaksdir+dif[0][0]
    p2 = peaksdir+dif[1][0]

    r1 = readsdir+dif[0][1]
    r2 = readsdir+dif[1][1]

    params = ("--pf narrowpeak "
              "--rf bam --pe -o {}") .format(outdir)

    manormDifPeaks(p1, p2, r1, r2, params)
