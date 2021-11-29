  import pandas as pd
  #pd.DataFrame(het_genes).to_csv("/media/lucas/Disc4T/Projects/PhD_Project/het_genes.csv", header = None, index = False)

  frascka = pd.read_csv("/media/lucas/Disc4T/Projects/PhD_Project/Fraschka/fraschka_all.csv", header = None)
  fra = set(frascka[0].tolist())

  both = het_genes & fra
  ours = het_genes - fra
  theirs = fra - het_genes

  len(het_genes)
  len(fra)

  len(both)
  len(ours)
  len(theirs)

  for x in ours:
      print(x)
