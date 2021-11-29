  import pybedtools as pb
  import pandas as pd
  import os

  wd = "/mnt/Disc4T/Projects/PhD_Project/Subtelomers/"
  os.chdir(wd)

  # Load telomeres and peak files
  subtel = pd.read_csv(
      "subtelomers.csv",
      header=0,
      names=["Chrom", "Right", "Left"],
      sep="\t"
  )

  subtel_dict = dict(zip(subtel["Chrom"], zip(subtel["Right"], subtel["Left"])))

  # ## Search for peaks outside telomers

  peaks = pb.BedTool("/mnt/Disc4T/Projects/Cristina_ChIP_All/New_Coverage/New_Peaks/all_peaks.bed")

  def checkSubtel(peak):

      """ Check if a region (in bed format) is outside subtelomeric regions.
      Returns True if a peak is outside subtelomers (as defined in a dict
      called subtel_dict) and False otherwise."""

      left = peak.start > subtel_dict[peak.chrom][0]*1000
      right = peak.end < subtel_dict[peak.chrom][1]*1000
      return(all([left, right]))

  island_peaks = peaks.filter(lambda peak: checkSubtel(peak))

  island_peaks.saveas("/mnt/Disc4T/Projects/PhD_Project/Island_Peaks/island_peaks.bed")

  Create a list of subtelomeric genes

  def isSubtel(gene):

      """ Check if a region (in bed format) is in the subtelomeric region.
      Returns True if a region is in subtelomers (as defined in a dict
      called subtel_dict) and False otherwise."""

      if gene.chrom in subtel_dict.keys():
          left = gene.start < subtel_dict[gene.chrom][0]*1000
          right = gene.end > subtel_dict[gene.chrom][1]*1000
      else:
          left, right = False, False

      return(any([left, right]))

  genes = pb.BedTool("/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-45_Pfalciparum3D7.gff")
  genes = genes.filter(lambda x: x[2] == "gene")
  genes = genes.sort()

  subtel_genes = genes.filter(lambda gene: isSubtel(gene))
  subtel_genes = set([x[8].split(";")[0].replace("ID=", "") for x in subtel_genes])

  outfile = "/mnt/Disc4T/Projects/PhD_Project/Island_Peaks/subtel_genes.csv"

  with open(outfile, "w+") as output:
      for gene in subtel_genes:
          output.write(gene+"\n")

  gff = "/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-45_Pfalciparum3D7.gff"
  bed = "/mnt/Disc4T/Projects/PhD_Project/Island_Peaks/island_peaks.bed"

  ref = pb.BedTool(gff)
  ref = ref.filter(lambda x: x[2] == "gene")
  ref = ref.sort()

  bed = pb.BedTool(bed)
  bed.count()

  raw_anot = bed.intersect(ref, wao=True)
  anot.count()

  def parseRawAnot(bed):

      entries = []

      for entry in bed:

          chrom, start, stop = entry.chrom, entry.start, entry.stop
          info = entry.fields[11]

          if info == ".":
              gid = "Non-coding"
              anot = "NA"
          else:
              gid = info.split(";")[0].replace("ID=", "")
              anot = info.split(";")[1].replace("description=", "")

          entries.append([chrom, start, stop, gid, anot])

      df = pd.DataFrame(entries)
      df.columns = ["Chrom", "Start", "Stop", "GeneID", "Annot"]
      return(df)


  anot = parseRawAnot(raw_anot)
  anot.to_csv("/mnt/Disc4T/Projects/PhD_Project/Island_Peaks/island_peaks_annotated.csv", index = None)
