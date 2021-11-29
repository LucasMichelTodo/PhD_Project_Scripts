  setwd("/media/lucas/Disc4T/Projects/PhD_Project/Microarrays/Results/")

  anep1 <- read.csv("./Anastasia_Epireset_R1/aMAFC_geneLevel.csv")
  anep2 <- read.csv("./Anastasia_Epireset_R2/aMAFC_geneLevel.csv")

  anli1 <- read.csv("./Anastasia_Lipids_R1/aMAFC_geneLevel.csv")
  anli2 <- read.csv("./Anastasia_Lipids_R2/aMAFC_geneLevel.csv")
  anli3 <- read.csv("./Anastasia_Lipids_R3/aMAFC_geneLevel.csv")

  elhs <- read.csv("./eli_heatshock/aMAFC_geneLevel.csv")
  elhs <- head(elhs[,c(1,2,11,16,17,26,31)])

  orind1 <- read.csv("./oriol_inductions_R1/aMAFC_geneLevel.csv")
  orind2 <- read.csv("./oriol_inductions_R2/aMAFC_geneLevel.csv")


  all_aMAFC <- plyr::join_all(list(anep1,
                                   anep2,
                                   anli1,
                                   anli2,
                                   anli3,
                                   elhs,
                                   orind1,
                                   orind2),
                              by = "X", type = "full")

  head(all_aMAFC)
  colnames(all_aMAFC)
