library(ggplot2)
library(tsne)
library(reshape2)
library(tidyverse)
library(scales)
require(gridExtra)

#### Load Data ####

wd <- "/mnt/Disc4T/Projects/Cristina_ChIP_All/Data/New_Coverage/Binsize_150/"
setwd(wd)

cov12B <- read.csv("./1.2B_me_sort_q5_noDup_RPKM_normInp_ps1_5binned_cov_2prevpost.csv")
cov10G <- read.csv("./10G_me_sort_q5_noDup_RPKM_normInp_ps1_5binned_cov_2prevpost.csv")
covA7K9 <- read.csv("./A7K9_me_sort_q5_noDup_RPKM_normInp_ps1_5binned_cov_2prevpost.csv")
covE5K9 <- read.csv("./E5K9_me_sort_q5_noDup_RPKM_normInp_ps1_5binned_cov_2prevpost.csv")
covB11 <- read.csv("./B11_me_sort_q5_noDup_RPKM_normInp_ps1_5binned_cov_2prevpost.csv")

strains <- list(df12B=cov12B,
                df10G=cov10G,
                dfA7K9=covA7K9,
                dfE5K9=covE5K9,
                dfB11=covB11
                )

## Change dfs in a list
strains <- lapply(strains, function(df) {
    colnames(df)[1] <- "Gene_id"
    df
})

## Convert list into individual objects again
list2env(strains, envir=.GlobalEnv)

##Create numeric non-na mtxs
nona.mtxs <- lapply(strains, function(df) {
    mtx <- df[complete.cases(df),-1]
    rownames(mtx) <- df[complete.cases(df),1]
    mtx
})

names(nona.mtxs) <- c("mtx12B",
                      "mtx10G",
                      "mtxA7K9",
                      "mtxE5K9",
                      "mtxB11"
                      )

##Create side-by-side matrices
mtx_all <- cbind(nona.mtxs$mtx12B,
                 nona.mtxs$mtx10G,
                 nona.mtxs$mtxA7K9,
                 nona.mtxs$mtxE5K9,
                 nona.mtxs$mtxB11
                 )

colnames(mtx_all) <- 1:115

## Load Gene Characteristics
variant <- read.csv("/home/lucas/ISGlobal/Gen_Referencies/Gens_variants_extended.txt",
                    header = TRUE, sep = "\t")

het_genes = read.csv("/mnt/Disc4T/Projects/PhD_Project/het_genes.csv",
                     header = F, stringsAsFactors = F)[,1]

fraschka = read.csv("/mnt/Disc4T/Projects/PhD_Project/Fraschka/fraschka_all.csv",
                    header = F, stringsAsFactors = F)[,1]

both = intersect(het_genes, fraschka)
ours = setdiff(het_genes, fraschka)
theirs = setdiff(fraschka, het_genes)

subtelomeric = read.csv("/mnt/Disc4T/Projects/PhD_Project/Island_Peaks/subtel_genes.csv",
                        header = F, stringsAsFactors = F)[,1]

## Create Gene info DF

## anot <- read.csv("/mnt/Disc4T/Projects/PhD_Project/chip_seq_genes_annotated.csv",
##                  sep = "\t", header = FALSE)

anot <- read.csv("/mnt/Disc4T/Projects/PhD_Project/Data/r_annot.csv", header=FALSE)
colnames(anot) <- c("Gene_id", "Annot")

geneDF = as.data.frame(strains$df12B$Gene_id)
colnames(geneDF) <- "Gene_id"

geneDF <- left_join(geneDF, anot, by = "Gene_id")
geneDF["Annot"] <- gsub("Plasmodium", "Pl.", geneDF$Annot)
geneDF["Annot"] <- gsub("protein", "prot.", geneDF$Annot)
geneDF["Annot"] <- gsub("membrane", "memb.", geneDF$Annot)
geneDF["Annot"] <- gsub("conserved", "cvd.", geneDF$Annot)
geneDF["Annot"] <- gsub("function", "func.", geneDF$Annot)
geneDF["Annot"] <- gsub("unknown", "ukwn.", geneDF$Annot)
geneDF["Annot"] <- gsub("exported", "xptd.", geneDF$Annot)
geneDF["Annot"] <- gsub("pseudogene", "pseudo", geneDF$Annot)
geneDF["Annot"] <- gsub("putative", "put.", geneDF$Annot)
geneDF["Annot"] <- gsub("%2C", "", geneDF$Annot)


geneDF["Epi"] <- "non-variant"
geneDF[geneDF$Gene_id %in% both, "Epi"] <- "Both"
geneDF[geneDF$Gene_id %in% ours, "Epi"] <- "Ours"
geneDF[geneDF$Gene_id %in% theirs, "Epi"] <- "Theirs"

geneDF["Subtelomeric"] <- "Normal"
geneDF[geneDF$Gene_id %in% subtelomeric, "Subtelomeric"] <- "Subtelomeric"
island <- geneDF$Epi != "non-variant" & geneDF$Subtelomeric != "Subtelomeric"
geneDF[island, "Subtelomeric"] <- "Island"

geneDFnona <- geneDF[geneDF$Gene_id %in% rownames(nona.mtxs$mtx12B),]

## Add transcription data

#transDF <- read.csv("/mnt/Disc4T/Projects/PhD_Project/External_Data/trasncription_parsed.csv")
transDF <- read.csv("/mnt/Disc4T/Projects/PhD_Project/Data/variantome_parsed.csv")

#### Create HeatmapDF ####
difpeaksDF <- data.frame(matrix(NA, nrow = length(geneDF$Gene_id), ncol = 72))

## Select ON and OFF examples for each dif. gene
for (i in 1:length(geneDF$Gene_id)){
  gene <- as.character(geneDF$Gene_id[i])
  means <- c()

  for (df in strains){
    ## We take into acount to decide "on"/"of" only the last 3 bins of 5' and first 2 of ORF.
    vect <- as.numeric(df[df$Gene_id == gene, 8:12])
    means <- c(means, mean(vect, na.rm = T))
  }

  onidx <- which.min(means)
  ondf <- strains[[onidx]]
  onstrain <- gsub("df", "", names(strains)[onidx])
  onvect <- ondf[ondf$Gene_id == gene, -1]

  offidx <- which.max(means)
  offdf <- strains[[offidx]]
  offstrain <- gsub("df", "", names(strains)[offidx])
  offvect <- offdf[offdf$Gene_id == gene, -1]

  difpeaksDF[i,] <- c(gene,
                      onvect,
                      offvect,
                      offvect-onvect,
                      onstrain,
                      offstrain)
}

prefix <- "On_"
seq <- 1:23
nameON <- paste0(prefix, seq)
prefix <- "Off_"
nameOFF <- paste0(prefix, seq)
prefix <- "Dif_"
nameDIFF <- paste0(prefix, seq)
names <- c(nameON, nameOFF, nameDIFF)

colnames(difpeaksDF)[1] <- "Gene_id"
colnames(difpeaksDF)[2:70] <- names
colnames(difpeaksDF)[71:72] <- c("On_strain", "Off_strain")

difpeaksDF <- left_join(difpeaksDF, geneDF, by="Gene_id")

difpeak_mtx <- as.matrix(difpeaksDF[,2:47])
rownames(difpeak_mtx) <- difpeaksDF$Gene_id

labels <- apply(difpeaksDF, 1, function(x){
  paste0(x[1], ": ", x[73])
})

difpeaksDF["Labels"] <- labels
melting_vars = c("Gene_id",
                 "On_strain",
                 "Off_strain",
                 "Annot",
                 "Epi",
                 "Subtelomeric",
                 "Labels",
                 "Cluster")

## Add transcription Data ##

completeDF <- left_join(difpeaksDF, transDF, by="Gene_id")

length(transDF$Gene_id)
length(unique(transDF$Gene_id))

on_trans <- rep(NA, (dim(completeDF)[1]))
on_trans[completeDF$On_strain == "10G"] <- completeDF$X10G[completeDF$On_strain == "10G"]
on_trans[completeDF$On_strain == "12B"] <- completeDF$X1.2B[completeDF$On_strain == "12B"]
on_trans[!completeDF$On_strain %in% c("12B", "10G")] <- completeDF$X3D7.B[!completeDF$On_strain %in% c("12B", "10G")]

off_trans <- rep(NA, (dim(completeDF)[1]))
off_trans[completeDF$Off_strain == "10G"] <- completeDF$X10G[completeDF$Off_strain == "10G"]
off_trans[completeDF$Off_strain == "12B"] <- completeDF$X1.2B[completeDF$Off_strain == "12B"]
off_trans[!completeDF$Off_strain %in% c("12B", "10G")] <- completeDF$X3D7.B[!completeDF$Off_strain %in% c("12B", "10G")]

completeDF["On"] <- on_trans
completeDF["Off"] <- off_trans

customMax <- function(vect){
  if (all(is.na(vect))) {
    x <- NA
  } else {
    x <- vect[which.max(abs(vect))]
  }
  return(x)
}

completeDF['Dif_12B.10G'] <- completeDF["X1.2B"] - completeDF["X10G"]
completeDF['Dif_12B.3D7B'] <- completeDF["X1.2B"] - completeDF["X3D7.B"]
completeDF['Dif_10G.3D7B'] <- completeDF["X10G"] - completeDF["X3D7.B"]

completeDF["MaxTrans_Dif"] <- apply(select(completeDF,
                                           'Dif_12B.10G',
                                           'Dif_12B.3D7B',
                                           'Dif_10G.3D7B'),
                                    1, customMax)

completeDF["TransDif"] <- completeDF["On"] - completeDF["Off"]
#write.csv(completeDF, "/mnt/Disc4T/Projects/PhD_Project/Data/r_complete_df.csv", row.names = F)

#### Heatmap Function ####

myHeatmap <- function(mdf){
  p <- ggplot(mdf, aes(x = variable, y = Labels, fill = value)) +

    geom_tile(colour="snow3") +

    theme(
      legend.position='bottom',
      legend.title = element_blank(),

      strip.background = element_blank(),
      panel.background=element_blank(),
      panel.border=element_blank(),panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),plot.background=element_blank(),

      axis.title = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank())
  p
}

myHeatmap_het <- function(mdf){
  p <- myHeatmap(mdf)
  p <- p + scale_fill_gradient(low = "white",
                               high = "red",
                               na.value="grey90") +

        scale_y_discrete(position = "right") +
        theme(strip.text.y = element_text(angle = 0)) +

        facet_grid(Cluster ~ Strain,
                   scales="free",
                   space="free")
  p
}

myHeatmap_Dif <- function(mdf){
  p <- myHeatmap(mdf)
  p <- p + scale_fill_gradient2(low = "#536DFE",
                                mid = "white",
                                high = "#00796B",
                                na.value="grey90") +

    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          strip.text.y=element_blank()) +

    facet_grid(Cluster ~ Strain,
               scales="free",
               space="free")
  p
}

myHeatmap_trans <- function(mdf){
  p <- myHeatmap(mdf)
  p <- p + scale_fill_gradient2(low = "blue",
                                mid = "white",
                                high = "yellow",
                                limits = c(-3,3),
                                oob=squish,
                                #labels=c("<-1",-0.5,0,0.5,">1"),
                                na.value="grey90") +
    scale_y_discrete(position = "left") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.y = element_text(angle = 180)) +

    facet_grid(Cluster ~ variable, scales="free", space="free", switch = "y")
  p
}

#### Clusterings ####

## ## Kmeans clustering
## kc_df <- noVARS[complete.cases(noVARS),]
## kclust <- kmeans(kc_df[,hetcols], nclust, nstart = 100)
## kc_df["Cluster"] <- kclust$cluster
## kc_dfsort <- arrange(kc_df, Cluster)
## kc_df$Labels <- factor(kc_df$Labels, levels = kc_dfsort$Labels)

customHClust <- function(df, nclust = 12){
  hetcols <- df %>% select(starts_with("On_"),
                           starts_with("Off_"),
                           -"On_strain",
                           -"Off_strain")
  dmtx <- dist(hetcols, method = "maximum")
  cl <- hclust(dmtx, method = 'complete')
  clusters <- cutree(cl, nclust)
  df$Labels <- factor(df$Labels, levels = df$Labels[cl$order])
  df["Cluster"] <- clusters
  return(df)
}

#### Call Het Heatmap ####

plotHclust <- function(df){
  hdf <- df %>% select(starts_with("On_"),
                       starts_with("Off_"),
                       melting_vars)
  mdf <- melt(hdf, id.vars = melting_vars)
  mdf["Strain"] <-sapply(mdf$variable, function(x) {strsplit(as.character(x), "_")[[1]][1]})
  phet <- myHeatmap_het(mdf)
}

## ggsave(phet,
##        file="/mnt/Disc4T/Projects/PhD_Project/Plots/heat_all_hc.png",
##        width = 50,
##        height = 40,
##        units="cm" )

#### Call Dif Heatmap ####

plotDif <- function(df){
  difDF <- df %>% select(starts_with("Dif_"),
                         -"Dif_12B.10G",
                         -"Dif_12B.3D7B",
                         -"Dif_10G.3D7B",
                         melting_vars)

  mdf <- melt(difDF, id.vars = melting_vars)
  mdf["Strain"] <- "Diff"
  pdif <- myHeatmap_Dif(mdf)
  return(pdif)
}

#### Plot transcription ####

plotTrans <- function(df){
  trans_var_heat <- df %>% select("On", "Off", "Labels", "Cluster")
  mdf <- melt(trans_var_heat, id.vars = c("Labels", "Cluster"))
  ptrans <- myHeatmap_trans(mdf)
  return(ptrans)
}

plotTransOneCol <- function(df){
  trans_var_heat <- df %>% select("TransDif", "Labels", "Cluster")
  mdf <- melt(trans_var_heat, id.vars = c("Labels", "Cluster"))
  ptrans <- myHeatmap_trans(mdf)
  return(ptrans)
}


## df <- completeDF[mask_noVARS & mask_difTrans & mask_difPeaks,]
## cdf <- customHClust(df)

## plotTrans(cdf)
## plotTransOneCol(cdf)

wholeHeatmap <- function(df, nclust = 12){
  clust <- customHClust(df, nclust)
  hplot <- plotHclust(clust)
  dplot <- plotDif(clust)
  tplot <- plotTrans(clust)

  all_plot <- grid.arrange(tplot, dplot, hplot, nrow = 1, widths = c(1,2.8,8))
  return(all_plot)
}

wholeHeatmap2 <- function(df, nclust = 12){
  clust <- customHClust(df, nclust)
  hplot <- plotHclust(clust)
  dplot <- plotDif(clust)
  tplot <- plotTransOneCol(clust)

  all_plot <- grid.arrange(tplot, dplot, hplot, nrow = 1, widths = c(1,2,8))
  return(all_plot)
}

#### Filtering ####

## Load data-set with genes that have diffrenential het. peaks
dif_genes = read.csv("/mnt/Disc4T/Projects/PhD_Project/diff_genes.csv",
                     header=F, stringsAsFactors=F)[,1]

## Load dataset to remove non-expressed vars, rifins, stevors and pfmc-2tm
expresed_vars <- read.csv("/mnt/Disc4T/Projects/PhD_Project/Data/var_rif_stevor_mc2tm_expression_rosetta_annotated_NOVA.txt", header=FALSE, sep="\t")
colnames(expresed_vars) <- c("Old_id", "Gene_id", "Anot", "Family")

## Filters
mask_difPeaks <- completeDF$Gene_id %in% dif_genes
mask_noVARS <- sapply(completeDF$Labels, function(x) grepl("rifin", x) |
                                                     grepl("PfEMP1", x) |
                                                     grepl("stevor", x) |
                                                     grepl("Pfmc-2TM", x))
mask_noVARS <- !mask_noVARS | completeDF$Gene_id %in% expresed_vars$Gene_id

mask_difTrans <- !is.na(completeDF$TransDif) & completeDF$TransDif > 2

## Only het.difPeaks, noVARS
x <- completeDF[mask_noVARS & mask_difPeaks,]

ggsave(wholeHeatmap(x),
       file="/mnt/Disc4T/Projects/PhD_Project/Plots/het_noVAR.png",
       width = 50,
       height = 40,
       units="cm" )

## Same but filtered by expression difference
x <- completeDF[mask_noVARS & mask_difTrans & mask_difPeaks,]

ggsave(wholeHeatmap(x, 4),
       file="/mnt/Disc4T/Projects/PhD_Project/Plots/het_noVAR_trans.png",
       width = 50,
       height = 40,
       units="cm" )

## Filter only with Trans

x <- completeDF[mask_difTrans,]

ggsave(wholeHeatmap(x),
       file="/mnt/Disc4T/Projects/PhD_Project/Plots/trans.png",
       width = 50,
       height = 40,
       units="cm" )

## Filter only with Trans and noVAR
## Aquest!
x <- completeDF[mask_difTrans & mask_noVARS,]

ggsave(wholeHeatmap(x, 4),
       file="/mnt/Disc4T/Projects/PhD_Project/Plots/trans_noVAR_nou.png",
       width = 40,
       height = 20,
       units="cm" )

#### Single Gene Plots ####

plot_gene <- function(gene_id, df, nclust = 12){
  clust <- customHClust(df, nclust)
  clust <- clust[clust$Gene_id == gene_id,]
  hplot <- plotHclust(clust)
  dplot <- plotDif(clust)
  tplot <- plotTrans(clust)

  all_plot <- grid.arrange(tplot, dplot, hplot, nrow = 1, widths = c(1,2,8))
  return(all_plot)
}


plot_gene("PF3D7_0320400", noVARS)
completeDF %>%
  filter(Gene_id == "PF3D7_0320400") %>%
  select("Gene_id", "On_strain", "Off_strain")

plot_gene("PF3D7_1241000", noVARS)
completeDF %>%
  filter(Gene_id == "PF3D7_1241000") %>%
  select("Gene_id", "On_strain", "Off_strain")

plot_gene("PF3D7_0601600", noVARS)
completeDF %>%
  filter(Gene_id == "PF3D7_0601600") %>%
  select("Gene_id", "On_strain", "Off_strain")

plot_gene("PF3D7_1149300", noVARS)
completeDF %>%
  filter(Gene_id == "PF3D7_1149300") %>%
  select("Gene_id", "On_strain", "Off_strain")

plot_gene("Custom_PF3D7_0935400_intergenic_lnc2", noVARS)
completeDF %>%
  filter(Gene_id == "Custom_PF3D7_0935400_intergenic_lnc2") %>%
  select("Gene_id", "On_strain", "Off_strain")

#### Tendency Plots ####

cluster_plot <- function(df, clust, nclust = 12){
  df <- customHClust(df, nclust)
  df <- df %>%
    filter(Cluster == clust) %>%
    select(starts_with("On_"),
           starts_with("Off_"),
           melting_vars)
  m <- melt(df, id.vars = melting_vars)
  m["State"] <- sapply(m$variable, function(x) strsplit(as.character(x), "_")[[1]][1])
  bins <- lapply(c(1:23), function(x) rep(x, dim(df)[1]))
  bins <- unlist(bins)
  bins <- c(bins, bins)
  m["Bin"] <- bins
  p <- ggplot(m, aes(x = bins, y = value, group = State))
  p <- p + geom_smooth(aes(color = State, fill = State))
  ## ggsave(p,
  ##        file = sprintf("/mnt/Disc4T/Projects/PhD_Project/Plots/Cluster_plots/cluster_%s.png",
  ##                       cluster))
  return(p)
}


for (clust in 1:12){
  p <- cluster_plot(noVARS, clust)
  print(p)
}

## Call heatmap
pk <- myDifHeatmap(hc_df)
print(pk)

## ggsave(pk,
##        file="/mnt/Disc4T/Projects/PhD_Project/Plots/heat_all_kc.png",
##        width = 50,
##        height = 40,
##        units="cm" )

## Create PCA DF
all_pca <- prcomp(noVARS[,2:47])

all_pca_df <- as.data.frame(all_pca$x[,c(1,2)])
all_pca_df <- cbind(noVARS[,c(1,48:51)], all_pca_df)
all_pca_df["HCluster"] <- hc_df$Cluster
all_pca_df["KCluster"] <- kc_df$Cluster

p <- ggplot(all_pca_df, aes(x=PC1, y=PC2, color = factor(KCluster))) +
  geom_point() +
  scale_fill_discrete()

p

p <- ggplot(all_pca_df, aes(x=PC1, y=PC2, color = factor(HCluster)))
p <- p + geom_point()
p <- p + scale_fill_discrete()
print(p)
