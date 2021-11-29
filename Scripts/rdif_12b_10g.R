library(ggplot2)
library(tsne)
library(reshape2)
library(scales)

#### Load Data ####
wd <- "/mnt/Disc4T/Projects/Cristina_ChIP_All/New_Coverage/Binsize50_deepTools/"
setwd(wd)

## "Old" coverage
## cov12B <- read.csv("./1.2B_me_sort_q5_noDup_RPKM_normInp_ps1_binned_cov_new.bed",
##                    stringsAsFactors = F)
## cov10G <- read.csv("./10G_me_sort_q5_noDup_RPKM_normInp_ps1_binned_cov_new.bed",
##                    stringsAsFactors = F)

## New coverage
cov12B <- read.csv("./1.2B_me_sort_q5_noDup_RPKM_normInp_ps1_5binned_cov_2prevpost.csv",
                                       stringsAsFactors = F)
cov10G <- read.csv("./10G_me_sort_q5_noDup_RPKM_normInp_ps1_5binned_cov_2prevpost.csv",
                                       stringsAsFactors = F)



strains <- list(df12B=cov12B, df10G=cov10G)

trans.data <- read.csv("/mnt/Disc4T/Projects/PhD_Project/External_Data/trasncription_parsed.csv",
                       stringsAsFactors = F)


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

names(nona.mtxs) <- c("mtx12B", "mtx10G")

##Create side-by-side matrices
mtx_12b10g <- cbind(nona.mtxs$mtx12B, nona.mtxs$mtx10G)
colnames(mtx_12b10g) <- 1:46

##Create diferences matrices
dif_12b10g <- nona.mtxs$mtx12B-nona.mtxs$mtx10G

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

anot <- read.csv("/mnt/Disc4T/Projects/PhD_Project/chip_seq_genes_annotated.csv",
                 sep = "\t", header = FALSE)

geneDF = as.data.frame(strains$df12B$Gene_id)
colnames(geneDF) <- "Gene_id"

geneDF["Annot"] <- anot$V2
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

geneDF["Trans"] <- trans.data$Dif_12B.10G
geneDF["Subtelomeric"] <- "Normal"
geneDF[geneDF$Gene_id %in% subtelomeric, "Subtelomeric"] <- "Subtelomeric"
island <- geneDF$Epi != "non-variant" & geneDF$Subtelomeric != "Subtelomeric"
geneDF[island, "Subtelomeric"] <- "Island"

geneDFnona <- geneDF[geneDF$Gene_id %in% rownames(nona.mtxs$mtx12B),]

cvgDF <- geneDFnona[geneDFnona$Epi != "non-variant",]
cvg_mtx <- dif_12b10g[rownames(dif_12b10g) %in% cvgDF$Gene_id,]

dif_pca <- prcomp(dif_12b10g)
cvg_pca <- prcomp(cvg_mtx)

dif_pca_df <- as.data.frame(dif_pca$x[, c(1,2)])
dif_pca_df <- cbind(dif_pca_df, geneDFnona[,-1])

cvg_pca_df <- as.data.frame(cvg_pca$x[, c(1,2)])
cvg_pca_df <- cbind(cvg_pca_df, cvgDF[,-1])


alpha <- sapply(dif_pca_df$Epi == "non-variant", function(x) if (x) {0.1} else {1})

#### PCA plots ####
## All-genes
p <- ggplot(dif_pca_df, aes(x=PC1, y=PC2, color = Epi))
p <- p + geom_point(alpha=alpha)
p

p <- ggplot(dif_pca_df, aes(x=PC1, y=PC2, color = Subtelomeric))
p <- p + geom_point(alpha=alpha)
p

p <- ggplot(dif_pca_df, aes(x=PC1, y=PC2, color = Trans))
p <- p + geom_point(alpha=alpha)
p <- p + scale_color_gradient2(midpoint = 0, low="red", mid = "black", high="green")
p

## CVGs
p <- ggplot(cvg_pca_df, aes(x=PC1, y=PC2, color = Epi))
p <- p + geom_point()
p

p <- ggplot(cvg_pca_df, aes(x=PC1, y=PC2, color = Subtelomeric))
p <- p + geom_point()
p

p <- ggplot(cvg_pca_df, aes(x=PC1, y=PC2, color = Trans))
p <- p + geom_point()
p <- p + scale_color_gradient2(midpoint = 0, low="red", mid = "black", high="green")
p

difpeaks12B <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Narrow_fe15/XLS_contrasts/Overlapped_and_filetred/calcFE/bed/annotated2/1.2B_10G_peaks_overlappandfilter_calcFE_annotated.csv",
                        sep = "\t",
                        stringsAsFactors = FALSE)

difpeaks10G <- read.csv("/home/lucas/ISGlobal/Chip_Seq/DATA/Aligns/q5/Narrow_fe15/XLS_contrasts/Overlapped_and_filetred/calcFE/bed/annotated2/10G_1.2B_peaks_overlappandfilter_calcFE_annotated.csv",
                        sep = "\t",
                        stringsAsFactors = FALSE)

hetdifgenes <-  unique(c(difpeaks12B$Gene, difpeaks10G$Gene))

hetdifDF <- geneDFnona[geneDFnona$Gene_id %in% hetdifgenes,]

hetdifDF
hetdif_mtx <- dif_12b10g[rownames(dif_12b10g) %in% hetdifDF$Gene_id,]

hetdif_pca <- prcomp(hetdif_mtx)

hetdif_pca_df <- as.data.frame(hetdif_pca$x[, c(1,2)])
hetdif_pca_df <- cbind(hetdif_pca_df, hetdifDF[,-1])

#### PCA plots ####
## All-genes
p <- ggplot(hetdif_pca_df, aes(x=PC1, y=PC2, color = Trans))
p <- p + geom_point()
p <- p + scale_color_gradient2(midpoint = 0, low="red", mid = "black", high="green")
p


hist(hetdifDF$Trans, breaks = 20)

customHeatmap <- function(df, limits){
  df <- melt(df)
  p <- ggplot(df, aes(x = variable, y = Gene_id, fill = value)) +
    geom_tile(colour="snow3",size=0.10) +
    scale_fill_gradient2(midpoint = 0,
                         low = "green",
                         mid = "black",
                         high = "red",
                         limits = limits,
                         oob=squish) +
    theme(
      strip.background = element_blank(),
      axis.title = element_blank(),
      strip.text.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank())
  p
}

## Absolute value of difference
heatDF <- cbind(hetdifDF[,-3], abs(hetdif_mtx))

## Make hierarquical Clustering
dmtx <- dist(abs(hetdif_mtx), method = "euclidean")
cl <- hclust(dmtx, method = 'average')
heatDF$Gene_id <- factor(heatDF$Gene_id, levels = heatDF$Gene_id[cl$order])
customHeatmap(heatDF, c(0,2))

## Difference
heatDF <- cbind(hetdifDF[,-3], hetdif_mtx)

## Make hierarquical Clustering
dmtx <- dist(hetdif_mtx, method = "euclidean")
cl <- hclust(dmtx, method = 'average')
heatDF$Gene_id <- factor(heatDF$Gene_id, levels = heatDF$Gene_id[cl$order])
customHeatmap(heatDF, c(-2,2))

difpeak_mtx <- mtx_12b10g[rownames(mtx_12b10g) %in% hetdifgenes,]
heatDF <- cbind(hetdifDF[,-4], difpeak_mtx)


heatDF[heatDF$Gene_id == "PF3D7_1220700",]

head(heatDF)
labels <- apply(heatDF, 1, function(x){
  paste0(x[1], ": ", x[2])
})

heatDF["Labels"] <- labels

## Make hierarquical Clustering
dmtx <- dist(difpeak_mtx, method = "euclidean")
cl <- hclust(dmtx, method = 'complete')
heatDF$Labels <- factor(heatDF$Labels, levels = heatDF$Labels[cl$order])
#customHeatmap(heatDF, c(-4,4))

#head(heatDF)
mdf <- melt(heatDF)
mdf["Strain"] <- sapply(mdf$variable, function(x) if (as.numeric(x) < 24) {"12B"} else {"10G"})
#mdf["Strain"] <- strain

p <- ggplot(mdf, aes(x = variable, y = Labels, fill = value)) +

  geom_tile(colour="snow3") +
            #size=0.10,
            #height=.9) +

  scale_fill_gradient2(midpoint = 0,
                       low = "white",
                       high = "red") +

  scale_y_discrete(position = "right") +

  theme(
    strip.background = element_blank(),
    axis.title = element_blank(),
    #strip.text.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()) +

  facet_grid(~Strain,
             scales="free_x",
             space="free")
p

##ggsave("/mnt/Disc4T/Projects/PhD_Project/dif12b_10g_heatmap.png", p,
##       device = "png", width = 40, height = 20,  units = "cm")


head(heatDF)

## Labeling each row of the melted df
## n <- dim(heatDF)[1]
## strain <- factor(c(rep("12Bpre2", n),
##                    rep("12Bpre2", n),
##                    rep("12Bpre1", n),
##                    rep("12Bpre1", n),
##                    rep(rep("12Bprom", 5), n),
##                    rep(rep("12Bbody", 5), n),
##                    rep(rep("12Bterm", 5), n),
##                    rep("12Bpost1", n),
##                    rep("12Bpost1", n),
##                    rep("12Bpost2", n),
##                    rep("12Bpost2", n),

##                    rep("10Gpre2", n),
##                    rep("10Gpre2", n),
##                    rep("10Gpre1", n),
##                    rep("10Gpre1", n),
##                    rep(rep("10Gprom", 5), n),
##                    rep(rep("10Gbody", 5), n),
##                    rep(rep("10Gterm", 5), n),
##                    rep("10Gpost1", n),
##                    rep("10Gpost1", n),
##                    rep("10Gpost2", n),
##                    rep("10Gpost2", n)),

##                  levels = c("12Bpre2", "12Bpre1",
##                             "12Bprom", "12Bbody", "12Bterm",
##                             "12Bpost1", "12Bpost2",
##                             "10Gpre2", "10Gpre1",
##                             "10Gprom", "10Gbody", "10Gterm",
##                             "10Gpost1", "10Gpost2"))
