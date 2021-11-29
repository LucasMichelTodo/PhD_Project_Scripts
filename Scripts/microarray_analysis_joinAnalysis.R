wd <- '/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/Join_Analysis/'
setwd(wd)

library(tidyverse)
library(eulerr)
library(readxl)
library(viridis)
library(ggdendro)
library(gridExtra)

## Load Old Dif Genes
old_fld <- '../Old_Arrays/R_results_OldArrays_Variantome/'
suffix <- '_log2FC2_red15_maxtime2.csv'

v12B_10G <- read_csv(paste0(old_fld, '12B_vs_10G', suffix)) %>%
         select(Gene_id) %>%
         pull()

v12B_3D7B <- read_csv(paste0(old_fld, '12B_vs_3D7B', suffix)) %>%
  select(Gene_id) %>%
  pull()

v10G_3D7B <- read_csv(paste0(old_fld, '10G_vs_3D7B', suffix)) %>%
  select(Gene_id) %>%
  pull()

old_gids <- unique(c(v12B_10G, v12B_3D7B, v10G_3D7B))


## Load New Dif Genes
new_fld <- '../New_Arrays/R_results_NewArray/'
suffix <- '_log2FC2_red15_maxtime2.csv'

vA7_B11 <- read_csv(paste0(new_fld, 'A7_vs_B11', suffix)) %>%
  select(Gene_id) %>%
  pull()

vA7_E5 <- read_csv(paste0(new_fld, 'A7_vs_E5', suffix)) %>%
  select(Gene_id) %>%
  pull()

vB11_E5 <- read_csv(paste0(new_fld, 'B11_vs_E5', suffix)) %>%
  select(Gene_id) %>%
  pull()

new_gids <- unique(c(vA7_B11, vA7_E5, vB11_E5))


## Plot
A <- old_gids
B <- new_gids
AB <- intersect(A, B)

ab <- length(AB)
a <- length(A[!A %in% AB])
b <- length(B[!B %in% AB])


fit <- euler(c(A=a, B=b, "A&B"=ab))

scales::viridis_pal()(2)

d <- plot(fit, fills = list(fill = c('#440154FF', "#FDE725FF"), alpha = 0.5),
          edges = list(lwd = 0.1),
          quantities = list(quantities = T),
          labels = list(labels=c("Old Arrays", "New Arrays")))

ggsave(d, filename = './join_Difs_Venn.pdf', device = "pdf")

plot(d)
print(fit)

library(viridis)

## Load Areas
exp_old <- read_csv('../Old_Arrays/R_results_OldArrays_Variantome/geneLevel_exp.csv') %>%
  select(-Name, -Variant, -Annot)

exp_new <- read_csv('../New_Arrays/R_results_NewArray/geneLevel_exp.csv') %>%
  select(-name, -Variant, -Annot)

exp_df <- full_join(exp_old, exp_new)
exp_df[exp_df$Gene_id == 'PF3D7_0935400_as',]$Gene_id <- 'PF3D7_0935390'

## Load Red Signal
red_old <- read_csv('../Old_Arrays/R_results_OldArrays_Variantome/geneLevel_redSignalexp.csv') %>%
  select(-Name, -Variant, -Annot)

red_new <- read_csv('../New_Arrays/R_results_NewArray/geneLevel_redSignal_exp.csv') %>%
  select(-name, -Variant, -Annot)

red_df <- full_join(red_old, red_new)
red_df[red_df$Gene_id == 'PF3D7_0935400_as',]$Gene_id <- 'PF3D7_0935390'

## Load info_df
info_df <- read_csv('../../../Binned_Coverage/info_df.csv')

## Load dif_genes
dif_df <- read_csv('../../../Binned_Coverage/max_log2FC2_filters_passed.csv')

## PCA Dif Genes

pca_df <- exp_df %>%
  filter(Gene_id %in% dif_df$Gene_id) %>%
  filter(complete.cases(.)) %>%
  pivot_longer(-Gene_id) %>%
  pivot_wider(names_from = Gene_id, values_from = value)

pca <- prcomp(pca_df %>% select(-name))
cmp1 <- format(round(summary(pca)$importance[2,1]*100, 2), nsmall = 2)
cmp2 <- format(round(summary(pca)$importance[2,2]*100, 2), nsmall = 2)

df_pca <- as.data.frame(pca$x)
head(df_pca)

strains <- sapply(pca_df$name, function(x) str_split(x, fixed('_'))[[1]][1])
tps <- sapply(pca_df$name, function(x) str_split(x, fixed('_'))[[1]][2])

?str_split

df_pca$Strain <- strains
df_pca$TimePoint <- tps

p <- ggplot(df_pca, aes(x=PC1,y=PC2, col = Strain, group = Strain))
p <- p + geom_point(aes(size= TimePoint))
p <- p + geom_path()
p <- p + scale_x_continuous(name=paste0("PC1: ", cmp1, "%"))
p <- p + scale_y_continuous(name=paste0("PC2: ", cmp2, "%"))
#p <- p + scale_color_viridis(discrete=TRUE, option="viridis")
p

ggsave(p, filename = "PCA_dif_genes.pdf", device = "pdf")

## PCA colapse genes into 1 val

maxexp_df <- exp_df %>%
  rowwise() %>%
  mutate(Max12B = max(c_across(contains('12B')))) %>%
  mutate(Max10G = max(c_across(contains('10G')))) %>%
  mutate(Max3D7B = max(c_across(contains('3D7B')))) %>%
  mutate(MaxA7 = max(c_across(contains('A7')))) %>%
  mutate(MaxE5 = max(c_across(contains('E5')))) %>%
  mutate(MaxB11 = max(c_across(contains('B11')))) %>%
  ungroup() %>%
  select(Gene_id, contains('Max'))

pca_df <- maxexp_df %>%
  #filter(Gene_id %in% dif_df$Gene_id) %>%
  filter(complete.cases(.)) %>%
  pivot_longer(-Gene_id) %>%
  pivot_wider(names_from = Gene_id, values_from = value)

pca <- prcomp(pca_df %>% select(-name))
cmp1 <- format(round(summary(pca)$importance[2,1]*100, 2), nsmall = 2)
cmp2 <- format(round(summary(pca)$importance[2,2]*100, 2), nsmall = 2)

df_pca <- as.data.frame(pca$x)
head(df_pca)

strains <- sapply(pca_df$name, function(x) gsub('Max', '', x, fixed = T))
df_pca$Strain <- strains

p <- ggplot(df_pca, aes(x=PC1,y=PC2, col = Strain, group = Strain))
p <- p + geom_point()
p <- p + geom_path()
p <- p + scale_x_continuous(name=paste0("PC1: ", cmp1, "%"))
p <- p + scale_y_continuous(name=paste0("PC2: ", cmp2, "%"))
##p <- p + scale_color_viridis(discrete=TRUE, option="viridis")
p

ggsave(p, filename = "PCA_collapsed_genes.pdf", device = "pdf")

my_percentile <- function(vector){
  ecdf(vector)(vector)*100
}

red_perc_df <- red_df %>%
  mutate(across(.cols = -Gene_id, .fns = my_percentile, .names = "Perc_{.col}")) %>%
  select(Gene_id, contains('Perc'))


## PCA Red Signal Perc

pca_df <- red_perc_df %>%
  #filter(Gene_id %in% dif_df$Gene_id) %>%
  filter(complete.cases(.)) %>%
  pivot_longer(-Gene_id) %>%
  pivot_wider(names_from = Gene_id, values_from = value)

pca <- prcomp(pca_df %>% select(-name))
cmp1 <- format(round(summary(pca)$importance[2,1]*100, 2), nsmall = 2)
cmp2 <- format(round(summary(pca)$importance[2,2]*100, 2), nsmall = 2)

df_pca <- as.data.frame(pca$x)

strains <- sapply(pca_df$name, function(x) str_split(x, fixed('_'))[[1]][2])
tps <- sapply(pca_df$name, function(x) str_split(x, fixed('_'))[[1]][3])

df_pca$Strain <- strains
df_pca$TimePoint <- tps

p <- ggplot(df_pca, aes(x=PC1,y=PC2, col = Strain, group = Strain))
p <- p + geom_point(aes(size= TimePoint))
p <- p + geom_path()
p <- p + scale_x_continuous(name=paste0("PC1: ", cmp1, "%"))
p <- p + scale_y_continuous(name=paste0("PC2: ", cmp2, "%"))
#p <- p + scale_color_viridis(discrete=TRUE, option="viridis")
p

ggsave(p, filename = "PCA_redperc_genes.pdf", device = "pdf")

## Collapse into 1 val

maxperc_df <- red_perc_df %>%
  rowwise() %>%
  mutate(Max12B = max(c_across(contains('12B')))) %>%
  mutate(Max10G = max(c_across(contains('10G')))) %>%
  mutate(Max3D7B = max(c_across(contains('3D7B')))) %>%
  mutate(MaxA7 = max(c_across(contains('A7')))) %>%
  mutate(MaxE5 = max(c_across(contains('E5')))) %>%
  mutate(MaxB11 = max(c_across(contains('B11')))) %>%
  ungroup() %>%
  select(Gene_id, contains('Max'))

pca_df <- maxexp_df %>%
  #filter(Gene_id %in% dif_df$Gene_id) %>%
  filter(complete.cases(.)) %>%
  pivot_longer(-Gene_id) %>%
  pivot_wider(names_from = Gene_id, values_from = value)

pca <- prcomp(pca_df %>% select(-name))
cmp1 <- format(round(summary(pca)$importance[2,1]*100, 2), nsmall = 2)
cmp2 <- format(round(summary(pca)$importance[2,2]*100, 2), nsmall = 2)

df_pca <- as.data.frame(pca$x)
head(df_pca)

strains <- sapply(pca_df$name, function(x) gsub('Max', '', x, fixed = T))
df_pca$Strain <- strains

p <- ggplot(df_pca, aes(x=PC1,y=PC2, col = Strain, group = Strain))
p <- p + geom_point()
p <- p + geom_path()
p <- p + scale_x_continuous(name=paste0("PC1: ", cmp1, "%"))
p <- p + scale_y_continuous(name=paste0("PC2: ", cmp2, "%"))
##p <- p + scale_color_viridis(discrete=TRUE, option="viridis")
p

ggsave(p, filename = "PCA_redperc_collapsed_genes.pdf", device = "pdf")

old_path <-'../Old_Arrays/R_results_OldArrays_Variantome/'
new_path <-'../New_Arrays/R_results_NewArray/'

old_df <- read_csv(paste0(old_path, 'old_arrays_final_df.csv'))
new_df <- read_csv(paste0(new_path, 'new_arrays_final_df.csv'))

old_areas <- read_csv(paste0(old_path, 'area_geneLevel.csv'))
new_areas <- read_csv(paste0(new_path, 'area_geneLevel.csv'))
new_areas[new_areas$Gene_id == 'PF3D7_0935400_as',]$Gene_id <- 'PF3D7_0935390'

old_max <- read_csv(paste0(old_path, 'all_aMAFC.csv'))
new_max <- read_csv(paste0(new_path, 'all_aMAFC.csv'))
new_max[new_max$Gene_id == 'PF3D7_0935400_as',]$Gene_id <- 'PF3D7_0935390'

## Load info_df
info_df <- read_csv('../../../Binned_Coverage/info_df.csv')

## OLD TOP DIFERENCES
## Quin timepoint agafem? El de màxima diferència. Entre quines soques?
## Què passa si un mateix gen té diferències a TPs diferents en contrastos diferents?

## Agafem per a cada gen el timepoint de màxima diferència del contrast amb màxima diferència.

dif_12B_10G <- read_csv(paste0(old_path, '12B_vs_10G_log2FC2_red15_maxtime.csv'))
dif_12B_3D7B <- read_csv(paste0(old_path, '12B_vs_3D7B_log2FC2_red15_maxtime.csv'))
dif_10G_3D7B <- read_csv(paste0(old_path, '10G_vs_3D7B_log2FC2_red15_maxtime.csv'))

table(dif_12B_10G$Dupl_Del)
table(dif_12B_3D7B$Dupl_Del)
table(dif_10G_3D7B$Dupl_Del)

tp_12B_10G <- old_max %>%
  select(Gene_id, `12B-10G_MaxTime`) %>%
  filter(Gene_id %in% dif_12B_10G$Gene_id)

tp_12B_3D7B <- old_max %>%
  select(Gene_id, `12B-3D7B_MaxTime`) %>%
  filter(Gene_id %in% dif_12B_3D7B$Gene_id)

tp_10G_3D7B <- old_max %>%
  select(Gene_id, `10G-3D7B_MaxTime`) %>%
  filter(Gene_id %in% dif_10G_3D7B$Gene_id)

x <- dif_12B_10G %>%
  select(Gene_id, `12B-10G_MaxVal`, Variant, Gam_specific, Dupl_Del) %>%
  rename(MaxVal = `12B-10G_MaxVal`) %>%
  mutate(Contrast = '12B_10G') %>%
  left_join(tp_12B_10G) %>%
  rename(MaxTime = `12B-10G_MaxTime`)

y <- dif_12B_3D7B %>%
  select(Gene_id, `12B-3D7B_MaxVal`, Variant, Gam_specific, Dupl_Del) %>%
  rename(MaxVal = `12B-3D7B_MaxVal`) %>%
  mutate(Contrast = '12B_3D7B') %>%
  left_join(tp_12B_3D7B) %>%
  rename(MaxTime = `12B-3D7B_MaxTime`)

z <- dif_10G_3D7B %>%
  select(Gene_id, `10G-3D7B_MaxVal`, Variant, Gam_specific, Dupl_Del) %>%
  rename(MaxVal = `10G-3D7B_MaxVal`) %>%
  mutate(Contrast = '10G_3D7B') %>%
  left_join(tp_10G_3D7B) %>%
  rename(MaxTime = `10G-3D7B_MaxTime`)

old_maxdifs <- bind_rows(x, y, z) %>%
  arrange(-abs(MaxVal)) %>%
  distinct(Gene_id, .keep_all = TRUE)

old_heat_areas <- NULL
new_heat_areas <- NULL
for (gid in old_maxdifs$Gene_id) {
  tp <- old_maxdifs %>%
    filter(Gene_id == gid) %>%
    select(MaxTime) %>%
    pull()

  old_a <- old_areas %>%
    filter(Gene_id == gid) %>%
    select(contains(tp))

  new_a <- new_areas %>%
    filter(Gene_id == gid) %>%
    select(contains(tp))

  if (dim(new_a)[1] == 0){
    new_a = tibble(
      c1 = as.numeric(NA),
      c2 = as.numeric(NA),
      c3 = as.numeric(NA)
    )
    }
  names(old_a) <- c('1.2B', '10G', '3D7-B')
  names(new_a) <- c('A7', 'B11', 'E5')

  old_heat_areas <- bind_rows(old_heat_areas, old_a)
  new_heat_areas <- bind_rows(new_heat_areas, new_a)
}
old_heat_df <- bind_cols(old_maxdifs, old_heat_areas)
new_heat_df <- bind_cols(old_maxdifs %>% select(Gene_id, Variant, Gam_specific, Dupl_Del), new_heat_areas)

old_heat_df <- old_heat_df %>%
  rowwise() %>%
  mutate(rwmean_centered_12B = `1.2B` - mean(c(`1.2B`, `10G`, `3D7-B`))) %>%
  mutate(rwmean_centered_10G = `10G` - mean(c(`1.2B`, `10G`, `3D7-B`))) %>%
  mutate(rwmean_centered_3D7B = `3D7-B` - mean(c(`1.2B`, `10G`, `3D7-B`)))

new_heat_df <- new_heat_df %>%
  rowwise() %>%
  mutate(rwmean_centered_A7 = A7 - mean(c(A7, E5, B11))) %>%
  mutate(rwmean_centered_E5 = E5 - mean(c(A7, E5, B11))) %>%
  mutate(rwmean_centered_B11 = B11 - mean(c(A7, E5, B11)))

old_pre_dupl_filtering <- old_heat_df
new_pre_dupl_filtering <- new_heat_df

old_heat_df <- old_heat_df %>%
  filter(!Dupl_Del)

new_heat_df <- new_heat_df %>%
  filter(!Dupl_Del)

## Ordering
mtx <- old_heat_df %>%
  select(contains('rwmean'))

##Make hierarquical Clustering
dmtx <- dist(scale(mtx), method = "euclidean")
cl <- hclust(dmtx, method = 'average')
old_heat_df$Gene_id <- factor(old_heat_df$Gene_id, levels = old_heat_df$Gene_id[cl$order])
new_heat_df$Gene_id <- factor(new_heat_df$Gene_id, levels = new_heat_df$Gene_id[cl$order])

old_tree <- ggdendrogram(cl, rotate = T)

old_m_infodf <- old_heat_df %>%
  select(Gene_id, Variant, Gam_specific) %>%
  gather(variable, value, -Gene_id)

mold_df <- old_heat_df %>%
  select(Gene_id, contains('rwmean')) %>%
  gather(variable, value, -Gene_id)

mnew_df <- new_heat_df %>%
  select(Gene_id, contains('rwmean')) %>%
  gather(variable, value, -Gene_id)

## old_heat_df %>%
##   left_join(info_df, by = 'Gene_id') %>%
##   select(
##     Gene_id,
##     Variant,
##     Gam_specific,
##     contains('rwmean')
##     ) %>%
##   write_tsv('old_arrays_maxFC_heatmap_df.tsv')

## NEW TOP DIFERENCES
## Quin timepoint agafem? El de màxima diferència. Entre quines soques?
## Què passa si un mateix gen té diferències a TPs diferents en contrastos diferents?

## Agafem per a cada gen el timepoint de màxima diferència del contrast amb màxima diferència.

dif_A7_E5 <- read_csv(paste0(new_path, 'A7_vs_E5_log2FC2_red15_maxtime.csv'))
dif_A7_B11 <- read_csv(paste0(new_path, 'A7_vs_B11_log2FC2_red15_maxtime.csv'))
dif_B11_E5 <- read_csv(paste0(new_path, 'B11_vs_E5_log2FC2_red15_maxtime.csv'))

tp_A7_E5 <- new_max %>%
  select(Gene_id, `A7-E5_MaxTime`) %>%
  filter(Gene_id %in% dif_A7_E5$Gene_id)

tp_A7_B11 <- new_max %>%
  select(Gene_id, `A7-B11_MaxTime`) %>%
  filter(Gene_id %in% dif_A7_B11$Gene_id)

tp_B11_E5 <- new_max %>%
  select(Gene_id, `B11-E5_MaxTime`) %>%
  filter(Gene_id %in% dif_B11_E5$Gene_id)

x <- dif_A7_E5 %>%
  select(Gene_id, `A7-E5_MaxVal`, Variant, Gam_specific, Dupl_Del) %>%
  rename(MaxVal = `A7-E5_MaxVal`) %>%
  mutate(Contrast = 'A7_E5') %>%
  left_join(tp_A7_E5) %>%
  rename(MaxTime = `A7-E5_MaxTime`)

y <- dif_A7_B11 %>%
  select(Gene_id, `A7-B11_MaxVal`, Variant, Gam_specific, Dupl_Del) %>%
  rename(MaxVal = `A7-B11_MaxVal`) %>%
  mutate(Contrast = 'A7_B11') %>%
  left_join(tp_A7_B11) %>%
  rename(MaxTime = `A7-B11_MaxTime`)

z <- dif_B11_E5 %>%
  select(Gene_id, `B11-E5_MaxVal`, Variant, Gam_specific, Dupl_Del) %>%
  rename(MaxVal = `B11-E5_MaxVal`) %>%
  mutate(Contrast = 'B11_E5') %>%
  left_join(tp_B11_E5) %>%
  rename(MaxTime = `B11-E5_MaxTime`)

new_maxdifs <- bind_rows(x, y, z) %>%
  arrange(-abs(MaxVal)) %>%
  distinct(Gene_id, .keep_all = TRUE)

new_old_heat_areas <- NULL
new_new_heat_areas <- NULL
for (gid in new_maxdifs$Gene_id) {

  tp <- new_maxdifs %>%
    filter(Gene_id == gid) %>%
    select(MaxTime) %>%
    pull()

  old_a <- old_areas %>%
    filter(Gene_id == gid) %>%
    select(contains(tp))

  new_a <- new_areas %>%
    filter(Gene_id == gid) %>%
    select(contains(tp))

  if (dim(old_a)[1] == 0){
    old_a = tibble(
      c1 = as.numeric(NA),
      c2 = as.numeric(NA),
      c3 = as.numeric(NA)
    )
    }
  names(old_a) <- c('1.2B', '10G', '3D7-B')
  names(new_a) <- c('A7', 'B11', 'E5')

  new_old_heat_areas <- bind_rows(new_old_heat_areas, old_a)
  new_new_heat_areas <- bind_rows(new_new_heat_areas, new_a)
}
new_old_heat_df <- bind_cols(new_maxdifs %>% select(Gene_id, Variant, Gam_specific, Dupl_Del), new_old_heat_areas)
new_new_heat_df <- bind_cols(new_maxdifs, new_new_heat_areas)

new_old_heat_df <- new_old_heat_df %>%
  rowwise() %>%
  mutate(rwmean_centered_12B = `1.2B` - mean(c(`1.2B`, `10G`, `3D7-B`))) %>%
  mutate(rwmean_centered_10G = `10G` - mean(c(`1.2B`, `10G`, `3D7-B`))) %>%
  mutate(rwmean_centered_3D7B = `3D7-B` - mean(c(`1.2B`, `10G`, `3D7-B`)))

new_new_heat_df <- new_new_heat_df %>%
  rowwise() %>%
  mutate(rwmean_centered_A7 = A7 - mean(c(A7, E5, B11))) %>%
  mutate(rwmean_centered_E5 = E5 - mean(c(A7, E5, B11))) %>%
  mutate(rwmean_centered_B11 = B11 - mean(c(A7, E5, B11)))


new_old_pre_dupl_filtering <- new_old_heat_df
new_new_pre_dupl_filtering <- new_new_heat_df

new_old_heat_df <- new_old_heat_df %>%
  filter(!Dupl_Del)

new_new_heat_df <- new_new_heat_df %>%
  filter(!Dupl_Del)



## Ordering
mtx <- new_new_heat_df %>%
  select(contains('rwmean'))

##Make hierarquical Clustering
dmtx <- dist(scale(mtx), method = "euclidean")
cl <- hclust(dmtx, method = 'average')
new_old_heat_df$Gene_id <- factor(new_old_heat_df$Gene_id, levels = new_old_heat_df$Gene_id[cl$order])
new_new_heat_df$Gene_id <- factor(new_new_heat_df$Gene_id, levels = new_new_heat_df$Gene_id[cl$order])

new_tree <- ggdendrogram(cl, rotate = T)

new_m_infodf <- new_new_heat_df %>%
  select(Gene_id, Variant, Gam_specific) %>%
  gather(variable, value, -Gene_id)

mnew_old_df <- new_old_heat_df %>%
  select(Gene_id, contains('rwmean')) %>%
  gather(variable, value, -Gene_id)

mnew_new_df <- new_new_heat_df %>%
  select(Gene_id, contains('rwmean')) %>%
  gather(variable, value, -Gene_id)

## Heatmap function
my_heatmap <- function(m_df){
  p <- ggplot(m_df, aes(x = variable, y = Gene_id, fill = value))
  p <- p + geom_tile(colour="snow3")
  p <- p + theme(
             ##text=element_text(size=24, family="Roboto"),
             legend.position='bottom',
             legend.title = element_blank(),

             panel.background=element_blank(),
             panel.grid.minor=element_blank(),
             plot.background=element_blank(),

             axis.title = element_blank(),
             axis.line.x = element_blank(),
             axis.ticks.x = element_blank(),
             axis.title.x = element_blank()
           )
  return(p)
}

my_areas_heatmap <- function(m_df){
  p <- my_heatmap(m_df)
  p <- p + scale_fill_gradient2(
             low = "#3da4ab",
             high = "#f6cd61",
             mid = 'black',
             na.value="grey",
             )
  return(p)
}

my_info_heatmap <- function(m_df){
  p <- my_heatmap(m_df)
  p <- p + scale_fill_manual(
             values = c('white', '#fe8a71')
           )
  p <- p + theme(
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),

               panel.border=element_blank(),
               panel.grid.major=element_blank(),

               strip.background = element_blank(),
               strip.text.x = element_blank(),
               strip.text.y = element_blank(),
             )
  return(p)
}

library(Cairo)
install.packages('Cairo')

## Old maxFC
p_old <- my_areas_heatmap(mold_df)
p_new <- my_areas_heatmap(mnew_df)
i <- my_info_heatmap(old_m_infodf)
whole_plot <- arrangeGrob(p_old, p_new, i, old_tree, nrow = 1, widths = c(2, 2, 1, 1))
plot(whole_plot)
ggsave('./Heatmaps/old_topdif_whole_new.pdf', whole_plot, device = 'pdf')

## New maxFC
p_new <- my_areas_heatmap(mnew_new_df)
p_old <- my_areas_heatmap(mnew_old_df)
i <- my_info_heatmap(new_m_infodf)
whole_plot <- arrangeGrob(p_new, p_old, i, new_tree, nrow = 1, widths = c(2, 2, 1, 1))
plot(whole_plot)
ggsave('./Heatmaps/new_topdif_whole_new.pdf', whole_plot, device = 'pdf')

## Load Old Dif genes
old_fld <- '../Old_Arrays/R_results_OldArrays_Variantome/'
suffix <- '_log2FC2_red15_maxtime.csv'

v12B_10G <- read_csv(paste0(old_fld, '12B_vs_10G', suffix))
v12B_3D7B <- read_csv(paste0(old_fld, '12B_vs_3D7B', suffix))
v10G_3D7B <- read_csv(paste0(old_fld, '10G_vs_3D7B', suffix))

## Load New Dif Genes
new_fld <- '../New_Arrays/R_results_NewArray/'
suffix <- '_log2FC2_red15_maxtime.csv'

vA7_B11 <- read_csv(paste0(new_fld, 'A7_vs_B11', suffix))
vA7_E5 <- read_csv(paste0(new_fld, 'A7_vs_E5', suffix))
vB11_E5 <- read_csv(paste0(new_fld, 'B11_vs_E5', suffix))

## Mix
n <- unique(c(vA7_B11$Gene_id, vB11_E5$Gene_id))
o <- unique(c(v12B_3D7B$Gene_id, v10G_3D7B$Gene_id))
intersect(n,o)


x <- v12B_3D7B %>%
  select(Gene_id, contains('MaxVal')) %>%
  full_join(v10G_3D7B %>% select(Gene_id, contains('MaxVal'))) %>%
  full_join(vA7_B11 %>% select(Gene_id, contains('MaxVal'))) %>%
  full_join(vB11_E5 %>% select(Gene_id, contains('MaxVal')))

x %>%
  filter(Gene_id %in% intersect(n,o)) %>%
  print(n = 100)


info_df %>%
  filter(Gene_id %in% intersect(n,o)) %>%
  select(Gene_id, Name, Annot, Variant, Gam_specific, Family, SubFamily) %>%
  write_csv('similarities_3D7A_B11.csv')

library(tidyverse)

getwd()

## Annot
info_df <- read_csv('./info_df.csv')

## Old arrays
old <- read_csv('../Microarrays/New_Old_separate_approach/Old_Arrays/R_results_OldArrays_Variantome/all_aMAFC.csv')

## New arrays
new <- read_csv('../Microarrays/New_Old_separate_approach/New_Arrays/R_results_NewArray/all_aMAFC.csv')

## Set thresholds
fc_th <- 1 ## FC considered to be a DE gene
dif_th <- 1 ## difference we alow between strains

### Get 3D7B DE genes
old
s3d7b_difs <- old %>%
  filter(abs(`12B-3D7B_MaxVal`) > fc_th | abs(`10G-3D7B_MaxVal`) > fc_th) %>%
  select(Gene_id, `12B-3D7B_MaxVal`, `10G-3D7B_MaxVal`)

s12B10G_3d7b <- s3d7b_difs %>%
  filter(sign(`12B-3D7B_MaxVal`) == sign(`10G-3D7B_MaxVal`)) %>% ## Check for different sign
  rowwise() %>%
  filter(abs(`12B-3D7B_MaxVal` - `10G-3D7B_MaxVal`) < dif_th) %>%
  mutate(Meandif_12B10G_3D7B = mean(c(`12B-3D7B_MaxVal`, `10G-3D7B_MaxVal`))) %>%
  ungroup()

### Get B11 DE genes
b11_difs <- new %>%
  filter(abs(`A7-B11_MaxVal`) > fc_th | abs(`B11-E5_MaxVal`) > fc_th) %>%
  select(Gene_id, `A7-B11_MaxVal`, `B11-E5_MaxVal`)

b11_a7e5 <- b11_difs %>%
  filter(sign(`A7-B11_MaxVal`) != sign(`B11-E5_MaxVal`)) %>% ## Check for different sign
  rowwise() %>%
  filter(abs(sum(`A7-B11_MaxVal`, `B11-E5_MaxVal`)) < dif_th) %>%
  mutate(Meandif_A7E5_B11 = mean(c(`A7-B11_MaxVal`, -`B11-E5_MaxVal`))) %>%
  ungroup()

## Join differences
join_meandifs <- b11_a7e5 %>%
  select(Gene_id, contains('Mean')) %>%
  left_join(s12B10G_3d7b %>% select(Gene_id, contains('Mean')), by = 'Gene_id') %>%
  filter(complete.cases(.))

final_df <- join_meandifs %>%
  filter(sign(Meandif_A7E5_B11) != sign(Meandif_12B10G_3D7B)) %>% ## Check for different sign
  rowwise() %>%
  filter(abs(Meandif_A7E5_B11 + Meandif_12B10G_3D7B) < 1) %>%
  #mutate(Meandif_12B10G_3D7B = mean(c(`12B-3D7B_MaxVal`, `10G-3D7B_MaxVal`))) %>%
  ungroup() %>%
  left_join(info_df, by = 'Gene_id')

final_df

final_df %>%
  write_tsv('B11_12B10G_similarities.tsv')

library(readxl)
library(tidyverse)

gdv1_dd <- read_xlsx('./GDV1_del_effects/gdv1_DD_supptable.xlsx')

gdv1_paper <- gdv1_dd %>%
  rename(Gene_id = `Gene Id`, GDV1_paper = `sig. up or down`) %>%
  select(Gene_id, GDV1_paper) %>%
  filter(GDV1_paper %in% c('u', 'd'))

our_list <- read_tsv('./GDV1_del_effects/B11_12B10G_similarities.tsv')

cross_df <- our_list %>%
  left_join(gdv1_paper)

cross_df %>%
  count(GDV1_paper)

cross_df %>%
  filter(!is.na(GDV1_paper)) %>%
  select(Gene_id, contains('Mean'), GDV1_paper, Annot)
