#### Imports and Dirs ####

library(ggplot2)
library(tidyverse)
library(readxl)
library(reshape2)
require(gridExtra)
library(readxl)
library(ggh4x)
library(ggrepel)
library(tsne)
library(scales)
library(viridis)
library(cluster)
library(NbClust)
library(factoextra)
library(class)

wd <- '/mnt/Disc4T/Projects/PhD_Project/Binned_Coverage/'
setwd(wd)

difpeaks_dir <- '/mnt/Disc4T/Projects/Chip_Seq_Data_2021/BetaFit_DuplDel_Filtered/'

info_df <- read_tsv('/mnt/Disc4T/Projects/PhD_Project/Data/PlasmoDB-52_Pfalciparum3D7_parsed_annotation.tsv')

## Flag tRNAs
info_df <- info_df %>%
  mutate(Is_tRNA = grepl('tRNA', Annot, fixed = T) & Type == 'ncRNA_gene')

#### Load gene lists ####
gene_lists <- read_excel('../Microarrays/New_Old_separate_approach/All gene lists_160719.xlsx', sheet = 2)

#### Load gene families data ####

gene_fam <- read_excel('./cvgfamilylist/Supplementary_table_2_CVG_list_161120_ap.xlsx', sheet = 2)
gene_fam <- gene_fam %>%
  rename(Gene_id = `Gene ID`,
         Gene_name = `Gene Name or Symbol`,
         SubFamily = `Family Detail`) %>%
  mutate(SubFamily = case_when(SubFamily == 'var pseudo,truncated or -like.' ~ 'var-like',
                               TRUE ~ SubFamily)) %>%
  mutate(Gene_name = ifelse(Gene_name == 'N/A', NA, Gene_name)) %>%
  select(Gene_id, Gene_name, Family, SubFamily)

gene_fam %>%
  filter(Family == 'OTHER') %>%
  mutate(NewFam = case_when(is.na(SubFamily) ~ Gene_name,
                            !is.na(SubFamily) ~ SubFamily))

bigfams <- c(
  'VAR',
  'FIKK',
  'HYP',
  'PHIST',
  'RIFIN',
  'STEVOR',
  'PFMC-2TM'
)


info_df <- info_df %>%
  left_join(gene_fam) %>%
  mutate(Name = ifelse(is.na(Name), Gene_name, Name)) %>%
  mutate(Name = ifelse(is.na(Name) & Family != 'OTHER', Family, Name)) %>%
  mutate(Name = ifelse(Gene_id == 'PF3D7_0935390', 'GDV1as', Name)) %>%
  mutate(Label = ifelse(is.na(Name), Gene_id, paste(Gene_id, Name, sep = ': '))) %>%
  mutate(Family_Grouped = case_when(
           Family %in% bigfams ~ Family,
           !is.na(Family) & !Family %in% bigfams ~ 'Other CVGs',
           is.na(Family) ~ 'Not CVGs'))

table(info_df$Name == info_df$Gene_name)

info_df %>%
  filter(Name != Gene_name) %>%
  select(Gene_id, Name, Gene_name) %>%
  print(n = 40)

#### Load variant genes data ####
cvgs <- read.csv2('../Microarrays/New_Old_separate_approach/New_Arrays/Files/taula_CVG_final.csv', stringsAsFactors = F)
cvgs <- cvgs %>%
  select(Gene_id = Gene.ID, Variant = Final.Customized) %>%
  mutate(Variant = ifelse(Variant == 'YES', TRUE, FALSE))

varlist <- filter(cvgs, Variant) %>% select(Gene_id) %>% pull()

info_df <- info_df %>%
  mutate(Variant = Gene_id %in% varlist)

## Probably Variant
## Some genes are from variant families but we have them as non-variant
info_df %>%
  filter(!Variant & !is.na(Family)) %>%
  select(Gene_id, Annot, Family) %>%
  print(n = 100)

## Check no variant gene has NA Family
info_df %>%
  filter(Variant & is.na(Family))

gams <- read_csv('../Microarrays/New_Old_separate_approach/Oriol_Suplementary/gam_table.csv') %>%
  select(Gene_id, Gam_specific)

info_df <- info_df %>%
  left_join(gams)

#### Load Dif.Peaks data ####

get_difpeak_gids <- function(df){
  df <- df %>%
    mutate(Gene_id = na_if(Gene_id, 'intergenic')) %>%
    select(Gene_id) %>%
    drop_na() %>%
    pull() %>%
    unique()
}
add_difpeaks_col <- function(df, strains){
  cn <- c('Chrom', 'Start', 'Stop', 'Peak_id', 'Score', 'Gene_id', 'Annot')
  file_str1 <- paste0(strains[1], '_vs_', strains[2], '_g200_l150_c1_c1.0_cond1_beta_cdf_099_IndelDup_filtered_annotated.csv')
  file_str2 <- paste0(strains[1], '_vs_', strains[2], '_g200_l150_c1_c1.0_cond2_beta_cdf_099_IndelDup_filtered_annotated.csv')
  peaks1 <- read_tsv(paste0(difpeaks_dir, file_str1), cn)
  peaks2 <- read_tsv(paste0(difpeaks_dir,file_str2), cn)
  peaks1 <- get_difpeak_gids(peaks1)
  peaks2 <- get_difpeak_gids(peaks2)
  df %>%
    mutate('Difpeaks_{strains[3]}_{strains[4]}' := Gene_id %in% peaks1 | Gene_id %in% peaks2)
}

strains_list <- list(
  c('1.2B', '10G', '12B', '10G'),
  c('A7K9', 'E5K9', 'A7', 'E5'),
  c('A7K9', 'B11', 'A7', 'B11'),
  c('B11', 'E5K9', 'B11', 'E5')
)

for (strains in strains_list) {info_df <- add_difpeaks_col(info_df, strains)}

## Change B11-E5 to E5-B11
info_df <- info_df %>%
  rename(Difpeaks_E5_B11 = Difpeaks_B11_E5)

write_csv(info_df, 'info_df.csv')

info_df %>%
  select(Gene_id, contains('Difpeaks')) %>%
  count(Difpeaks_E5_B11)

info_df %>%
  count(Family)

novar_family <- info_df %>%
  filter((!Variant & !is.na(Family)) | (Variant & is.na(Family)))

write_csv(novar_family, 'noVariant_with_family.csv')

#### Load transcription data ####
trans_df_old <- read_csv('../Microarrays/New_Old_separate_approach/Old_Arrays/R_results_OldArrays_Variantome/final_summary_table.csv')
trans_df_new <- read_csv('../Microarrays/New_Old_separate_approach/New_Arrays/R_results_NewArray/final_summary_table.csv')
trans_df <- full_join(trans_df_old, trans_df_new) %>%
  select(-Name, -Annot, -GamGene)

## Change GDV1as ID from the arrays to new id
trans_df[trans_df$Gene_id == 'PF3D7_0935400_as',]$Gene_id <- 'PF3D7_0935390'

## Subset trans_df to genes that appear on info_df
trans_df <- trans_df %>%
  filter(Gene_id %in% info_df$Gene_id)

trans_df <- trans_df %>%
  mutate(`E5-B11_MaxVal` = -`B11-E5_MaxVal`) %>%
  mutate(`E5-B11_MaxTime` = `B11-E5_MaxTime`) %>%
  select(-`B11-E5_MaxVal`, -`B11-E5_MaxTime`)


## Load Areas Data
arees_old <- read_csv('../Microarrays/New_Old_separate_approach/Old_Arrays/R_results_OldArrays_Variantome/areaDiferences_geneLevel.csv') %>%
  rename(Gene_id = ...1)

arees_new <- read_csv('../Microarrays/New_Old_separate_approach/New_Arrays/R_results_NewArray/areaDiferences_geneLevel.csv') %>%
  rename(Gene_id = ...1)

arees_df <- full_join(arees_old, arees_new)
arees_df[arees_df$Gene_id == 'PF3D7_0935400_as',]$Gene_id <- 'PF3D7_0935390'

## Load Red Percentile Data
red_percent_old <- read_csv('../Microarrays/New_Old_separate_approach/Old_Arrays/R_results_OldArrays_Variantome/red_percentiles.csv')
red_percent_new <- read_csv('../Microarrays/New_Old_separate_approach/New_Arrays/R_results_NewArray/red_percentiles.csv')

red_df <- full_join(red_percent_old, red_percent_new, by = 'Gene_id')
red_df[red_df$Gene_id == 'PF3D7_0935400_as',]$Gene_id <- 'PF3D7_0935390'

red_df

## Load Max-Time data
old_maxtime <- read_csv('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/Old_Arrays/R_results_OldArrays_Variantome/old_arrays_maxtime.csv')

new_maxtime <- read_csv('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/New_Arrays/R_results_NewArray/new_arrays_maxtime.csv')


maxtime_df <- full_join(old_maxtime, new_maxtime, by = 'Gene_id', suffix = c('_Old', '_New'))
maxtime_df[maxtime_df$Gene_id == 'PF3D7_0935400_as',]$Gene_id <- 'PF3D7_0935390'


old_breaks <- read_csv('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/Old_Arrays/R_results_OldArrays_Variantome/old_area_breaks.csv')

new_breaks <- read_csv('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/New_Arrays/R_results_NewArray/new_area_breaks.csv')

breaks_df <- bind_cols(old_breaks, new_breaks)
colnames(breaks_df) <- c('Old_Area_Breaks', 'New_Area_Breaks')

## Load DuplDel Data
f_path <- '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/Duplication_Deletion_Regions_Mean_Separate_DuplDel/Crossed_with_genes/'
suffix <- '_in_sort_q5_RPKMs_bymean_dupl_del_fact_1.75up_0.1dw_minlen500_mergelen_200_filtered_genes.tsv'
prefix <- c('1.2B', '10G', 'A7K9', 'E5K9', 'B11')

get_dupl_del <- function(prefix){
  fname <- paste0(prefix, suffix)
  df <- read_tsv(paste0(f_path, fname), col_names = F) %>%
    rename(Gene_id = X1, Name = X2, Annot = X3)
}

dupl_del <- lapply(prefix, get_dupl_del)
names(dupl_del) <- prefix

#### Load heterochromatin data 1000bp+500CDS ####
cov_12b <- read_tsv('./Genewise_Coverages/cov_1000tss_500orf_1.2B.bed', col_names = F)
cov_10g <- read_tsv('./Genewise_Coverages/cov_1000tss_500orf_10G.bed', col_names = F)
#cov_3d7b <- read_tsv('./Genewise_Coverages/cov_1000tss_500orf_3D7.bed', col_names = F)
cov_a7 <- read_tsv('./Genewise_Coverages/cov_1000tss_500orf_A7K9.bed', col_names = F)
cov_e5 <- read_tsv('./Genewise_Coverages/cov_1000tss_500orf_E5K9.bed', col_names = F)
cov_b11 <- read_tsv('./Genewise_Coverages/cov_1000tss_500orf_B11.bed', col_names = F)

join_cols = c('X1', 'X2', 'X3', 'X4')

het_df <- cov_12b %>%
  full_join(cov_10g, by = join_cols) %>%
  #full_join(cov_3d7b, by = join_cols) %>%
  full_join(cov_a7, by = join_cols) %>%
  full_join(cov_e5, by = join_cols) %>%
  full_join(cov_b11, by = join_cols)

colnames(het_df) <- c('Chrom', 'Start', 'Stop', 'Gene_id',
                      'Het_12B', 'Het_10G', 'Het_A7', 'Het_E5', 'Het_B11')

het_df <- het_df %>% select(Gene_id, contains('Het'), everything())

hdf <- het_df %>% select(Gene_id, contains('Het'))
cnames <- str_replace(colnames(hdf), 'Het', 'Cov_1000fp_500orf')
hdf <- hdf %>% set_names(cnames)

#### Load heterochromatin data 5pORF3p ####
## Load 3prime, ORF, 5prime Coverage

load3ORF5 <- function(cov_5ORF3_file){
  cov_5ORF3 <- read_tsv(cov_5ORF3_file, col_names = F)
  cov_5ORF3 <- cov_5ORF3 %>%
    setNames(c('Chrom', 'Start', 'Stop', 'Gene_id', 'Intensity', 'Strand', 'Cov')) %>%
    select(-Intensity)
  prime5 <- cov_5ORF3 %>% filter(grepl('5prime', fixed = T, Gene_id))
  ORF <- cov_5ORF3 %>% filter(!grepl('5prime', fixed = T, Gene_id) & !grepl('3prime', fixed = T, Gene_id))
  prime3 <- cov_5ORF3 %>% filter(grepl('3prime', fixed = T, Gene_id))
  ORF['Cov_5prime'] <- prime5$Cov
  ORF['Cov_3prime'] <- prime3$Cov
  ORF <- ORF %>% rename(Cov_ORF = Cov)
  return(ORF)
}

e5_5ORF3 <- load3ORF5('./Genewise_Coverages/cov_10005p_orf_10003p_E5K9.bed')
a7_5ORF3 <- load3ORF5('./Genewise_Coverages/cov_10005p_orf_10003p_A7K9.bed')
b11_5ORF3 <- load3ORF5('./Genewise_Coverages/cov_10005p_orf_10003p_B11.bed')
x12b_5ORF3 <- load3ORF5('./Genewise_Coverages/cov_10005p_orf_10003p_1.2B.bed')
x10g_5ORF3 <- load3ORF5('./Genewise_Coverages/cov_10005p_orf_10003p_10G.bed')
#x3d7b_5ORF3 <- load3ORF5('./Genewise_Coverages/cov_10005p_orf_10003p_3D7.bed')

join_cols = c('Chrom', 'Start', 'Stop', 'Gene_id', 'Strand')

cov_5ORF3_df <- x12b_5ORF3 %>%
  full_join(x10g_5ORF3, by = join_cols, suffix = c('_12B', '_10G')) %>%
  #full_join(x3d7b_5ORF3, by = join_cols, suffix = c('', '_3D7B')) %>%
  full_join(a7_5ORF3, by = join_cols, suffix = c('', '_A7')) %>%
  full_join(e5_5ORF3, by = join_cols, suffix = c('', '_E5')) %>%
  full_join(b11_5ORF3, by = join_cols, suffix = c('', '_B11')) %>%
  rename(Cov_ORF_A7 = Cov_ORF, Cov_5prime_A7 = Cov_5prime, Cov_3prime_A7 = Cov_3prime) %>%
  select(all_of(join_cols), contains('5prime'), contains('ORF'), contains('3prime'))

cov_5orf3 <- cov_5ORF3_df %>% select(Gene_id, contains('Cov'))
cnames <- str_replace(colnames(cov_5orf3), '5prime', '1000fp') %>%
  str_replace('ORF', 'allorf') %>%
  str_replace('3prime', '1000tp')
cov_5orf3 <- cov_5orf3 %>% set_names(cnames)

#### Load ORF abs-bp Coverages ####

read_abs_orf_cov <- function(strain, len){
  col_names <- c('Chrom', 'Start', 'Stop', 'Gene_id', 'Cov')
  path <- paste0('./Genewise_Coverages/cov_', len, 'orf_allowoverlap_', strain, '.bed')
  df <- read_tsv(path, col_names = col_names) %>%
    select(Gene_id, Cov)
  return(df)
}

cov_500_12B <- read_abs_orf_cov('1.2B', '500')
cov_500_10G <- read_abs_orf_cov('10G', '500')
cov_500_A7 <- read_abs_orf_cov('A7K9', '500')
cov_500_E5 <- read_abs_orf_cov('E5K9', '500')
cov_500_B11 <- read_abs_orf_cov('B11', '500')

cov_1000_12B <- read_abs_orf_cov('1.2B', '1000')
cov_1000_10G <- read_abs_orf_cov('10G', '1000')
cov_1000_A7 <- read_abs_orf_cov('A7K9', '1000')
cov_1000_E5 <- read_abs_orf_cov('E5K9', '1000')
cov_1000_B11 <- read_abs_orf_cov('B11', '1000')

cov_1500_12B <- read_abs_orf_cov('1.2B', '1500')
cov_1500_10G <- read_abs_orf_cov('10G', '1500')
cov_1500_A7 <- read_abs_orf_cov('A7K9', '1500')
cov_1500_E5 <- read_abs_orf_cov('E5K9', '1500')
cov_1500_B11 <- read_abs_orf_cov('B11', '1500')

cov_orf_abs_df <- cov_500_12B %>%
  full_join(cov_500_10G, by = 'Gene_id') %>%
  full_join(cov_500_A7, by = 'Gene_id') %>%
  full_join(cov_500_E5, by = 'Gene_id') %>%
  full_join(cov_500_B11, by = 'Gene_id') %>%

  full_join(cov_1000_12B, by = 'Gene_id') %>%
  full_join(cov_1000_10G, by = 'Gene_id') %>%
  full_join(cov_1000_A7, by = 'Gene_id') %>%
  full_join(cov_1000_E5, by = 'Gene_id') %>%
  full_join(cov_1000_B11, by = 'Gene_id') %>%

  full_join(cov_1500_12B, by = 'Gene_id') %>%
  full_join(cov_1500_10G, by = 'Gene_id') %>%
  full_join(cov_1500_A7, by = 'Gene_id') %>%
  full_join(cov_1500_E5, by = 'Gene_id') %>%
  full_join(cov_1500_B11, by = 'Gene_id') %>%

  set_names('Gene_id',
            'Cov_500orf_12B', 'Cov_500orf_10G', 'Cov_500orf_A7', 'Cov_500orf_E5', 'Cov_500orf_B11',
            'Cov_1000orf_12B', 'Cov_1000orf_10G', 'Cov_1000orf_A7', 'Cov_1000orf_E5', 'Cov_1000orf_B11',
            'Cov_1500orf_12B', 'Cov_1500orf_10G', 'Cov_1500orf_A7', 'Cov_1500orf_E5', 'Cov_1500orf_B11')

cov_orf_abs_df %>%
  filter(!complete.cases(.))


cov_orf_abs_df %>%
  select(contains('12B'))

cov_orf_abs_df %>%
  select(contains('10G'))

load_pDB <- function(cov_pDB_file){

  col_names <- c('Chrom', 'Start', 'Stop', 'Gene_id', 'Type', 'Cov')
  cov_pDB <- read_tsv(cov_pDB_file, col_names = col_names) %>%
    mutate(Gene_id = str_replace(Gene_id, '.*(PF3D7_\\d{7}).*', '\\1'))

  prime5 <- cov_pDB %>%
    filter(Type == 'five_prime_UTR') %>%
    rename(p5utr_Cov = Cov) %>%
    select(Gene_id, p5utr_Cov)

  prime3 <- cov_pDB %>%
    filter(Type == 'three_prime_UTR') %>%
    rename(p3utr_Cov = Cov) %>%
    select(Gene_id, p3utr_Cov)

  outdf <- prime5 %>%
    full_join(prime3, by = 'Gene_id') %>%
    group_by(Gene_id) %>%
    summarize(p5utr_Cov = mean(p5utr_Cov), p3utr_Cov = mean(p3utr_Cov))

  return(outdf)
}

cov_plasmoDB_12B <- load_pDB('./Genewise_Coverages/cov_plasmoDB_UTRs_1.2B.bed')
cov_plasmoDB_10G <- load_pDB('./Genewise_Coverages/cov_plasmoDB_UTRs_10G.bed')
#cov_plasmoDB_3D7B <- load_pDB('./Genewise_Coverages/cov_plasmoDB_UTRs_3D7.bed')
cov_plasmoDB_A7 <- load_pDB('./Genewise_Coverages/cov_plasmoDB_UTRs_A7K9.bed')
cov_plasmoDB_E5 <- load_pDB('./Genewise_Coverages/cov_plasmoDB_UTRs_E5K9.bed')
cov_plasmoDB_B11 <- load_pDB('./Genewise_Coverages/cov_plasmoDB_UTRs_B11.bed')

cov_pDB_df <- cov_plasmoDB_12B %>%
  full_join(cov_plasmoDB_10G, by = 'Gene_id', suffix = c('_12B', '_10G')) %>%
  #full_join(cov_plasmoDB_3D7B, by = 'Gene_id', suffix = c('', '_3D7B')) %>%
  full_join(cov_plasmoDB_A7, by = 'Gene_id', suffix = c('', '_A7')) %>%
  full_join(cov_plasmoDB_E5, by = 'Gene_id', suffix = c('', '_E5')) %>%
  full_join(cov_plasmoDB_B11, by = 'Gene_id', suffix = c('', '_B11')) %>%
  rename(p5utr_Cov_A7 = p5utr_Cov, p3utr_Cov_A7 = p3utr_Cov)


load_pDB_TSS <- function(cov_pDB_file){
  col_names <- c('Chrom', 'Start', 'Stop', 'Gene_id', 'TSS_Cov')
  cov_pDB <- read_tsv(cov_pDB_file, col_names = col_names) %>%
    mutate(Gene_id = str_replace(Gene_id, '.*(PF3D7_\\d{7}).*', '\\1')) %>%
    select(Gene_id, TSS_Cov) %>%
    group_by(Gene_id) %>%
    summarize(TSS_Cov = mean(TSS_Cov))

  return(cov_pDB)
}

cov_plasmoDB_TSS_12B <- load_pDB_TSS('./Genewise_Coverages/cov_plasmoDB_TSSs_1.2B.bed')
cov_plasmoDB_TSS_10G <- load_pDB_TSS('./Genewise_Coverages/cov_plasmoDB_TSSs_10G.bed')
#cov_plasmoDB_TSS_3D7B <- load_pDB_TSS('./Genewise_Coverages/cov_plasmoDB_TSSs_3D7.bed')
cov_plasmoDB_TSS_A7 <- load_pDB_TSS('./Genewise_Coverages/cov_plasmoDB_TSSs_A7K9.bed')
cov_plasmoDB_TSS_E5 <- load_pDB_TSS('./Genewise_Coverages/cov_plasmoDB_TSSs_E5K9.bed')
cov_plasmoDB_TSS_B11 <- load_pDB_TSS('./Genewise_Coverages/cov_plasmoDB_TSSs_B11.bed')

cov_pDB_TSS_df <- cov_plasmoDB_TSS_12B %>%
  full_join(cov_plasmoDB_TSS_10G, by = 'Gene_id', suffix = c('_12B', '_10G')) %>%
  #full_join(cov_plasmoDB_TSS_3D7B, by = 'Gene_id', suffix = c('', '_3D7B')) %>%
  full_join(cov_plasmoDB_TSS_A7, by = 'Gene_id', suffix = c('', '_A7')) %>%
  full_join(cov_plasmoDB_TSS_E5, by = 'Gene_id', suffix = c('', '_E5')) %>%
  full_join(cov_plasmoDB_TSS_B11, by = 'Gene_id', suffix = c('', '_B11')) %>%
  rename(TSS_Cov_A7 = TSS_Cov)

cov_pDB_df <- cov_pDB_df %>%
  full_join(cov_pDB_TSS_df)

#### Load Many-Bin coverages ####

get_bins <- function(dfin, bin_txt, outcol){
  dfout <- dfin %>% filter(grepl(bin_txt, fixed = T, Gene_id)) %>%
    mutate(Gene_id = str_replace(Gene_id, bin_txt, '')) %>%
    rename(!!outcol := Cov) %>%
    select(Gene_id, Strand, !!outcol)
}

get_manybins <- function(strain){

  fpath <- paste0('./Genewise_Coverages/cov_3ORF5_manybins_allowoverlap_', strain, '.bed')
  cov_5orf3_manybins <- read_tsv(fpath, col_names = F)%>%
    setNames(c('Chrom', 'Start', 'Stop', 'Gene_id', 'Intensity', 'Strand', 'Cov')) %>%
    select(-Intensity)

  cov_2000fp_manybins <- get_bins(cov_5orf3_manybins, '_5prime_2000', 'Cov_2000fp_manybins')
  cov_1500fp_manybins <- get_bins(cov_5orf3_manybins, '_5prime_1500', 'Cov_1500fp_manybins')
  cov_1000fp_manybins <- get_bins(cov_5orf3_manybins, '_5prime_1000', 'Cov_1000fp_manybins')
  cov_500fp_manybins <- get_bins(cov_5orf3_manybins, '_5prime_500', 'Cov_500fp_manybins')
  cov_1qorf_manybins <- get_bins(cov_5orf3_manybins, '_q1', 'Cov_1qorf_manybins')
  cov_2qorf_manybins <- get_bins(cov_5orf3_manybins, '_q2', 'Cov_2qorf_manybins')
  cov_3qorf_manybins <- get_bins(cov_5orf3_manybins, '_q3', 'Cov_3qorf_manybins')
  cov_4qorf_manybins <- get_bins(cov_5orf3_manybins, '_q4', 'Cov_4qorf_manybins')
  cov_1500tp_manybins <- get_bins(cov_5orf3_manybins, '_3prime_1500', 'Cov_1500tp_manybins')
  cov_1000tp_manybins <- get_bins(cov_5orf3_manybins, '_3prime_1000', 'Cov_1000tp_manybins')
  cov_500tp_manybins <- get_bins(cov_5orf3_manybins, '_3prime_500', 'Cov_500tp_manybins')

  join_cols <- c('Strand', 'Gene_id')
  df <- cov_2000fp_manybins %>%
    full_join(cov_1500fp_manybins, by = join_cols) %>%
    full_join(cov_1000fp_manybins, by = join_cols) %>%
    full_join(cov_500fp_manybins, by = join_cols) %>%
    full_join(cov_1qorf_manybins, by = join_cols) %>%
    full_join(cov_2qorf_manybins, by = join_cols) %>%
    full_join(cov_3qorf_manybins, by = join_cols) %>%
    full_join(cov_4qorf_manybins, by = join_cols) %>%
    full_join(cov_1500tp_manybins, by = join_cols) %>%
    full_join(cov_1000tp_manybins, by = join_cols) %>%
    full_join(cov_500tp_manybins, by = join_cols)

  return(df)
}

cov_manybins_12B <- get_manybins('1.2B')
cov_manybins_10G <- get_manybins('10G')
#cov_manybins_3D7B <- get_manybins('3D7')
cov_manybins_A7 <- get_manybins('A7K9')
cov_manybins_E5 <- get_manybins('E5K9')
cov_manybins_B11 <- get_manybins('B11')


## cov_manybins_B11 %>%
##   filter(Gene_id == 'PF3D7_1130800') %>%
##   print(width = 800)


join_cols <- c('Gene_id', 'Strand')
cov_manybins <- cov_manybins_12B %>%
  full_join(cov_manybins_10G, by = join_cols, suffix = c('_12B', '_10G')) %>%
  #full_join(cov_manybins_3D7B, by = join_cols, suffix = c('', '_3D7B')) %>%
  full_join(cov_manybins_A7, by = join_cols, suffix = c('', '_A7')) %>%
  full_join(cov_manybins_E5, by = join_cols, suffix = c('', '_E5')) %>%
  full_join(cov_manybins_B11, by = join_cols, suffix = c('', '_B11')) %>%
  rename(Cov_2000fp_manybins_A7 = Cov_2000fp_manybins,
         Cov_1500fp_manybins_A7 = Cov_1500fp_manybins,
         Cov_1000fp_manybins_A7 = Cov_1000fp_manybins,
         Cov_500fp_manybins_A7 = Cov_500fp_manybins,
         Cov_1qorf_manybins_A7 = Cov_1qorf_manybins,
         Cov_2qorf_manybins_A7 = Cov_2qorf_manybins,
         Cov_3qorf_manybins_A7 = Cov_3qorf_manybins,
         Cov_4qorf_manybins_A7 = Cov_4qorf_manybins,
         Cov_1500tp_manybins_A7 = Cov_1500tp_manybins,
         Cov_1000tp_manybins_A7 = Cov_1000tp_manybins,
         Cov_500tp_manybins_A7 = Cov_500tp_manybins,
         ) %>%
  select(-Strand)

## Convert to numeric
cov_manybins <- cov_manybins %>%
  mutate(across(contains('Cov'), as.numeric))

cov_manybins %>%
  filter(!complete.cases(.))

## get_correlations <- function(strains, trans_filter, difpeaks_filter){
##     df1 <- paste0('cov_manybins_', strains[1])
##     df2 <- paste0('cov_manybins_', strains[2])

##     dif_2000fp_manybins <- get(df1)$Cov_2000fp_manybins - get(df2)$Cov_2000fp_manybins
##     dif_1500fp_manybins <- get(df1)$Cov_1500fp_manybins - get(df2)$Cov_1500fp_manybins
##     dif_1000fp_manybins <- get(df1)$Cov_1000fp_manybins - get(df2)$Cov_1000fp_manybins
##     dif_500fp_manybins <- get(df1)$Cov_500fp_manybins - get(df2)$Cov_500fp_manybins

##     dif_1q_manybins <- get(df1)$Cov_1qorf_manybins - get(df2)$Cov_1qorf_manybins
##     dif_2q_manybins <- get(df1)$Cov_2qorf_manybins - get(df2)$Cov_2qorf_manybins
##     dif_3q_manybins <- get(df1)$Cov_3qorf_manybins - get(df2)$Cov_3qorf_manybins
##     dif_4q_manybins <- get(df1)$Cov_4qorf_manybins - get(df2)$Cov_4qorf_manybins

##     dif_1500tp_manybins <- get(df1)$Cov_1500tp_manybins - get(df2)$Cov_1500tp_manybins
##     dif_1000tp_manybins <- get(df1)$Cov_1000tp_manybins - get(df2)$Cov_1000tp_manybins
##     dif_500tp_manybins <- get(df1)$Cov_500tp_manybins - get(df2)$Cov_500tp_manybins

##     cor_cov_df <- bind_cols(list(
##       Gene_id = get(df1)$Gene_id,
##       dif_2000fp_manybins = dif_2000fp_manybins,
##       dif_1500fp_manybins = dif_1500fp_manybins,
##       dif_1000fp_manybins = dif_1000fp_manybins,
##       dif_500fp_manybins = dif_500fp_manybins,
##       dif_1q_manybins = dif_1q_manybins,
##       dif_2q_manybins = dif_2q_manybins,
##       dif_3q_manybins = dif_3q_manybins,
##       dif_4q_manybins = dif_4q_manybins,
##       dif_1500tp_manybins = dif_1500tp_manybins,
##       dif_1000tp_manybins = dif_1000tp_manybins,
##       dif_500tp_manybins = dif_500tp_manybins
##       ))


##     full_df <- inner_join(cor_cov_df, trans_df, by = 'Gene_id')  %>%
##       left_join(info_df)

##     trans_col <- paste0(strains[1], '-', strains[2], '_MaxVal')
##     difpeaks_col <- paste0('Difpeaks_', strains[1], '_', strains[2])

##     cor_df <- full_df %>%
##       filter(abs(get(trans_col)) > trans_filter) %>%
##       select(Gene_id, matches(trans_col), contains('dif_'))

##     if (difpeaks_filter){df_dif <- df_dif %>% filter(get(difpeaks_col))}

##     cormtx <- cor(cor_df %>% select(-Gene_id) %>% drop_na(), method = 'pearson')
##     print(paste0('Strains: ', strains[1], ' vs ', strains[2]))
##     print(paste0('Number of genes: ', dim(cor_df)[1]))
##     print(cormtx)
##     print('+++++++++++++++++++++++++++++++++++++++++++++++')
##     outdf <- as.data.frame(cormtx)
##     outdf['Corr_with'] <- rownames(outdf)
##     outdf <- outdf %>% select(Corr_with, everything())
##     outname = paste0('corr_', strains[1], '_', strains[2], '_transth_', trans_filter, '_difpeaks_', difpeaks_filter, '.csv')
##     write_csv(outdf, file = outname)
## }

## tf <- 1.5
## df <- FALSE

## contrasts <- list(
##   c('12B', '10G'),
##   c('A7', 'E5'),
##   c('A7', 'B11'),
##   c('B11', 'E5')
## )

## for (c in contrasts) {get_correlations(c, tf, df)}


## Don't, computer can't handle it!

load_abs_ORF <- function(cov_abs_orf_file){
  cov_abs_ORF <- read_tsv(cov_abs_orf_file, col_names = F)
  cov_abs_ORF <- cov_abs_ORF %>%
    setNames(c('Chrom', 'Start', 'Stop', 'Gene_id', 'Intensity', 'Strand', 'Cov')) %>%
    select(-Intensity)
  ORF <- cov_abs_ORF %>%
    filter(!grepl('5prime', fixed = T, Gene_id) & !grepl('3prime', fixed = T, Gene_id))
  ORF <- ORF %>%
    rename(Cov_ORF = Cov) %>%
    mutate(Cov_ORF = as.double(Cov_ORF)) %>%
    select(Gene_id, Cov_ORF)
  return(ORF)
}

abs_ORF_500_12B <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_0_500_1.2B.bed')
abs_ORF_500_10G <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_0_500_10G.bed')
abs_ORF_500_A7 <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_0_500_A7K9.bed')
abs_ORF_500_E5 <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_0_500_E5K9.bed')
abs_ORF_500_B11 <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_0_500_B11.bed')

abs_ORF_1000_12B <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_500_1000_1.2B.bed')
abs_ORF_1000_10G <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_500_1000_10G.bed')
abs_ORF_1000_A7 <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_500_1000_A7K9.bed')
abs_ORF_1000_E5 <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_500_1000_E5K9.bed')
abs_ORF_1000_B11 <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_500_1000_B11.bed')

abs_ORF_1500_12B <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_1000_1500_1.2B.bed')
abs_ORF_1500_10G <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_1000_1500_10G.bed')
abs_ORF_1500_A7 <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_1000_1500_A7K9.bed')
abs_ORF_1500_E5 <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_1000_1500_E5K9.bed')
abs_ORF_1500_B11 <- load_abs_ORF('./Genewise_Coverages/cov_abs_ORF_1000_1500_B11.bed')

join_cols = c('Gene_id')

cov_abs_ORF <- abs_ORF_500_12B %>%
  full_join(abs_ORF_500_10G, by = join_cols, suffix = c('_0_500_12B', '_0_500_10G')) %>%
  full_join(abs_ORF_500_A7, by = join_cols, suffix = c('', '_0_500_A7')) %>%
  full_join(abs_ORF_500_E5, by = join_cols, suffix = c('', '_0_500_E5')) %>%
  full_join(abs_ORF_500_B11, by = join_cols, suffix = c('', '_0_500_B11')) %>%

  full_join(abs_ORF_1000_12B, by = join_cols, suffix = c('', '_500_1000_12B')) %>%
  full_join(abs_ORF_1000_10G, by = join_cols, suffix = c('', '_500_1000_10G')) %>%
  full_join(abs_ORF_1000_A7, by = join_cols, suffix = c('', '_500_1000_A7')) %>%
  full_join(abs_ORF_1000_E5, by = join_cols, suffix = c('', '_500_1000_E5')) %>%
  full_join(abs_ORF_1000_B11, by = join_cols, suffix = c('', '_500_1000_B11')) %>%

  full_join(abs_ORF_1500_12B, by = join_cols, suffix = c('', '_1000_1500_12B')) %>%
  full_join(abs_ORF_1500_10G, by = join_cols, suffix = c('', '_1000_1500_10G')) %>%
  full_join(abs_ORF_1500_A7, by = join_cols, suffix = c('', '_1000_1500_A7')) %>%
  full_join(abs_ORF_1500_E5, by = join_cols, suffix = c('', '_1000_1500_E5')) %>%
  full_join(abs_ORF_1500_B11, by = join_cols, suffix = c('', '_1000_1500_B11')) %>%

  rename(Cov_ORF_0_500_A7 = Cov_ORF) %>%
  select(Gene_id, contains('Cov_'))

## Create 'fixed' color palette
col_factor <- factor(unique(info_df$Family_Grouped),
                     levels=c('Not CVGs',
                              'Other CVGs',
                              'FIKK',
                              'HYP',
                              'PHIST',
                              'STEVOR',
                              'RIFIN',
                              'VAR',
                              'PFMC-2TM'))

my_colors <- c('gray', scales::viridis_pal(begin = 0, end = 1)(8))
names(my_colors) <- levels(col_factor)
my_scale <- scale_fill_manual(name = "Family_Grouped", values = my_colors)

cor_cov_df <- cov_orf_abs_df %>%
  full_join(hdf, by = 'Gene_id') %>%
  full_join(cov_5orf3, by = 'Gene_id') %>%
  full_join(cov_pDB_df, by = 'Gene_id') %>%
  full_join(cov_manybins, by = 'Gene_id') %>%
  full_join(cov_abs_ORF, by = 'Gene_id')

positive_cov_df <- cor_cov_df %>%
  mutate(across(-Gene_id, ~ ifelse(.x > 0, .x, 0)))

full_df <- inner_join(cor_cov_df, trans_df)  %>%
  left_join(info_df)

#### Correlations ####

get_cor_mtx <- function(df, contrast, trans_filter, difpeaks_filter){

  ## contrast <- c('12B', '10G')
  ## df <- full_df

  trans_col <- paste0(contrast[1], '-', contrast[2], '_MaxVal')
  difpeaks_col <- paste0('Difpeaks_', contrast[1], '_', contrast[2])


  ## df %>%
  ##   select(contains('Cov_ORF_'))

  df_dif <- df %>%
  ## Absolute bp ORF
    mutate(Cov_500orf_Dif = get(paste0('Cov_500orf_', contrast[1])) - get(paste0('Cov_500orf_', contrast[2]))) %>%
    mutate(Cov_1000orf_Dif = get(paste0('Cov_1000orf_', contrast[1])) - get(paste0('Cov_1000orf_', contrast[2]))) %>%
    mutate(Cov_1500orf_Dif = get(paste0('Cov_1500orf_', contrast[1])) - get(paste0('Cov_1500orf_', contrast[2]))) %>%

  ## Absolute bp ORF, by interval
    mutate(Cov_0_500orf_Dif = get(paste0('Cov_ORF_0_500_', contrast[1])) - get(paste0('Cov_ORF_0_500_', contrast[2]))) %>%
    mutate(Cov_500_1000orf_Dif = get(paste0('Cov_ORF_500_1000_', contrast[1])) - get(paste0('Cov_ORF_500_1000_', contrast[2]))) %>%
    mutate(Cov_1000_1500orf_Dif = get(paste0('Cov_ORF_1000_1500_', contrast[1])) - get(paste0('Cov_ORF_1000_1500_', contrast[2]))) %>%

  ## 100bp 5prime + 500bp ORF
    mutate(Cov_1000fp_500orf_Dif = get(paste0('Cov_1000fp_500orf_', contrast[1])) - get(paste0('Cov_1000fp_500orf_', contrast[2]))) %>%

  ## 100bp 5prime/ ORF / 1000bp 3prime
    mutate(Cov_1000fp_Dif = get(paste0('Cov_1000fp_', contrast[1])) - get(paste0('Cov_1000fp_', contrast[2]))) %>%
    mutate(Cov_allorf_Dif = get(paste0('Cov_allorf_', contrast[1])) - get(paste0('Cov_allorf_', contrast[2]))) %>%
    mutate(Cov_1000tp_Dif = get(paste0('Cov_1000tp_', contrast[1])) - get(paste0('Cov_1000tp_', contrast[2]))) %>%

  ## PlsmoDB 5prime and TSS
    mutate(Cov_plasmoDB_5p_Dif = get(paste0('p5utr_Cov_', contrast[1])) - get(paste0('p5utr_Cov_', contrast[2]))) %>%
    mutate(Cov_plasmoDB_TSS_Dif = get(paste0('TSS_Cov_', contrast[1])) - get(paste0('TSS_Cov_', contrast[2]))) %>%
    mutate(Cov_plasmoDB_3p_Dif = get(paste0('p3utr_Cov_', contrast[1])) - get(paste0('p3utr_Cov_', contrast[2]))) %>%

  ## Manybins (2000, 1500, 1000, 500, 1q, 2q, 3q, 4q, 500, 1000, 1500)
    mutate(Cov_2000fp_manybins_Dif = get(paste0('Cov_2000fp_manybins_', contrast[1])) - get(paste0('Cov_2000fp_manybins_', contrast[2]))) %>%
    mutate(Cov_1500fp_manybins_Dif = get(paste0('Cov_1500fp_manybins_', contrast[1]))  - get(paste0('Cov_1500fp_manybins_', contrast[2]))) %>%
    mutate(Cov_1000fp_manybins_Dif = get(paste0('Cov_1000fp_manybins_', contrast[1])) - get(paste0('Cov_1000fp_manybins_', contrast[2]))) %>%
    mutate(Cov_500fp_manybins_Dif = get(paste0('Cov_500fp_manybins_', contrast[1]))- get(paste0('Cov_500fp_manybins_', contrast[2]))) %>%

    mutate(Cov_1qorf_manybins_Dif = get(paste0('Cov_1qorf_manybins_', contrast[1])) - get(paste0('Cov_1qorf_manybins_', contrast[2]))) %>%
    mutate(Cov_2qorf_manybins_Dif = get(paste0('Cov_2qorf_manybins_', contrast[1]))  - get(paste0('Cov_2qorf_manybins_', contrast[2]))) %>%
    mutate(Cov_3qorf_manybins_Dif = get(paste0('Cov_3qorf_manybins_', contrast[1])) - get(paste0('Cov_3qorf_manybins_', contrast[2]))) %>%
    mutate(Cov_4qorf_manybins_Dif = get(paste0('Cov_4qorf_manybins_', contrast[1]))- get(paste0('Cov_4qorf_manybins_', contrast[2]))) %>%

    mutate(Cov_1500tp_manybins_Dif = get(paste0('Cov_1500tp_manybins_', contrast[1]))  - get(paste0('Cov_1500tp_manybins_', contrast[2]))) %>%
    mutate(Cov_1000tp_manybins_Dif = get(paste0('Cov_1000tp_manybins_', contrast[1])) - get(paste0('Cov_1000tp_manybins_', contrast[2]))) %>%
    mutate(Cov_500tp_manybins_Dif = get(paste0('Cov_500tp_manybins_', contrast[1]))- get(paste0('Cov_500tp_manybins_', contrast[2]))) %>%

    filter(abs(get(trans_col)) > trans_filter) %>%
    select(Gene_id, matches(trans_col), contains('_Dif'))

  ## df_dif %>%
  ##   select(Gene_id, Cov_0_500orf_Dif, Cov_500_1000orf_Dif, Cov_1000_1500orf_Dif) %>%
  ##   summary()

  ## table(df_dif$Cov_500_1000orf_Dif)

  if (difpeaks_filter){df_dif <- df_dif %>% filter(get(difpeaks_col))}

  cormtx <- cor(df_dif %>% select(-Gene_id) %>% drop_na(), method = 'pearson')
  print(paste0('Contrast ', contrast[1], ' vs ', contrast[2]))
  print(paste0('Number of genes: ', dim(df_dif %>% select(-Gene_id) %>% drop_na())[1]))
  print(as.data.frame(cormtx[,1]))
  print('+++++++++++++++++++++++++++++')

  outdf <- as.data.frame(cormtx)
  outdf['Corr_with'] <- rownames(outdf)
  outdf <- outdf %>% select(Corr_with, everything())
  outname = paste0('corr_ALL_', contrast[1], '_', contrast[2],
                   '_transth_', trans_filter, '_difpeaks_', difpeaks_filter, '.csv')
  write_csv(outdf, file = outname)
}

contrasts <- list(
  c('12B', '10G'),
  c('A7', 'E5'),
  c('A7', 'B11'),
  c('E5', 'B11')
)

trans_filter <- 0
difpeaks_filter <- F

for (contrast in contrasts) {get_cor_mtx(full_df, contrast, trans_filter, difpeaks_filter)}

cor_12B_10G <- get_cor_mtx(full_df, c('12B', '10G'), trans_filter, difpeaks_filter)
cor_A7_E5 <- get_cor_mtx(full_df, c('A7', 'E5'), trans_filter, difpeaks_filter)
cor_A7_B11 <- get_cor_mtx(full_df, c('A7', 'B11'), trans_filter, difpeaks_filter)
cor_B11_E5 <- get_cor_mtx(full_df, c('E5', 'B11'), trans_filter, difpeaks_filter)

all_Corr <- bind_cols(Corr_With = cor_12B_10G[-1,1],
                      Cor_12B_10G = cor_12B_10G[-1,2],
                      Cor_A7_E5 = cor_A7_E5[-1,2],
                      Cor_A7_B11 = cor_A7_B11[-1,2],
                      Cor_B11_E5 = cor_B11_E5[-1,2])

write_csv(all_Corr, file = paste0('./Correlations_transth_', trans_filter, '_allowoverlaps.csv'))

#### Correlation Plots ####

myCorPlot <- function(df, region, s1, s2, trans_th){

  df <- plot_5ORF3_df

  trans <- paste0(s1,'-',s2,'_MaxVal')
  difpeak <- paste0('Difpeaks_', s1, '_', s2)
  title_str <- paste0(s1, ' vs ', s2)
  outfile <- paste0('./Plots/corplot_', s1, '_', s2, '_', region, '.png')

  plot_df <- df %>%
    mutate(Dif_het = get(paste0('Cov_', region, '_', s1)) - get(paste0('Cov_', region, '_', s2))) %>%
    select(Gene_id, matches(trans), Dif_het, matches(difpeak)) %>%
    setNames(c('Gene_id', 'Trans', 'Het', 'Difpeak')) %>%
    filter(abs(Trans) > trans_th)

  p <- ggplot(plot_df, aes(x=Het,
                           y=Trans)) +
    scale_alpha_discrete(range = c(0.2, 1)) +
    ggtitle(title_str) +
    ylab('aMAFC') +
    xlab(paste0('Heterochromatin difference at ', region)) +
    geom_point()
  ggsave(outfile, p, device = 'png')
  print(p)
}

plot_5ORF3_df <- inner_join(cov_5ORF3_df, trans_df)  %>%
  left_join(info_df)

regions <- c('5prime', 'ORF', '3prime')

contrasts <- list(
  c('12B', '10G'),
  c('A7', 'E5'),
  c('A7', 'B11'),
  c('E5', 'B11')
)

## ## Call corr plots for 5p/ORF/3p coverage
## for (region in regions){
##   for (contrast in contrasts){
##     s1 <- contrast[1]
##     s2 <- contrast[2]
##     trans_th <- 0
##     myCorPlot(plot_5ORF3_df, region, s1, s2, trans_th)
##   }
## }

region <- 'ORF'
s1 <- 'A7'
s2 <- 'B11'
trans_th <- 1
myCorPlot(plot_5ORF3_df, region, s1, s2, trans_th)

#### Correlation Plots ####

myCorPlot_all <- function(df, region, s1, s2, trans_th){

  trans <- paste0(s1,'-',s2,'_MaxVal')
  difpeak <- paste0('Difpeaks_', s1, '_', s2)
  title_str <- paste0(s1, ' vs ', s2)
  outfile <- paste0('./Plots/Strain_Corr/corplot_',
                    s1, '_', s2, '_', region,
                    'transth_', trans_th,
                    '.png')

  plot_df <- df %>%
    mutate(Dif_het = get(paste0(region, s1)) - get(paste0(region, s2))) %>%
    select(Gene_id, matches(trans), Dif_het, matches(difpeak)) %>%
    setNames(c('Gene_id', 'Trans', 'Het', 'Difpeak')) %>%
    filter(abs(Trans) > trans_th)

  p <- ggplot(plot_df, aes(x=Het,
                      y=Trans,
                      color = Difpeak,
                      alpha = Difpeak)) +
    scale_alpha_discrete(range = c(0.2, 1)) +
    ggtitle(title_str) +
    ylab('aMAFC') +
    xlab(paste0('Heterochromatin difference at ', region)) +
    geom_point()
  ggsave(outfile, p, device = 'png')
  print(p)
}

## Make plots for all strains and coverages
cnms <- str_replace(colnames(full_df %>% select(contains('Cov') & contains('12B'))), '12B', '')

contrasts <- list(
  c('12B', '10G'),
  c('A7', 'E5'),
  c('A7', 'B11'),
  c('E5', 'B11')
)

## ## Plot correlations for all Coverages
## for (c in cnms){
##   for (contrast in contrasts){
##     s1 <- contrast[1]
##     s2 <- contrast[2]
##     trans_th <- 1
##     myCorPlot_all(full_df, c, s1, s2, trans_th)
##   }
## }

region <- 'Cov_allorf_'
s1 <- 'A7'
s2 <- 'E5'
trans_th <- 1
myCorPlot_all(full_df, region, s1, s2, trans_th)

hist(cor_cov_df$Cov_500orf_12B - cor_cov_df$Cov_500orf_10G, breaks = 500)
hist(cor_cov_df$Cov_500orf_12B - cor_cov_df$Cov_500orf_A7, breaks = 500)
hist(cor_cov_df$Cov_500orf_A7 - cor_cov_df$Cov_500orf_B11, breaks = 500)

x <- positive_cov_df %>%
  filter(across(contains('Cov')) != 0)

hist(x$Cov_500orf_12B - x$Cov_500orf_10G, breaks = 5000)
hist(x$Cov_500orf_12B - x$Cov_500orf_A7, breaks = 5000)
hist(x$Cov_500orf_A7 - x$Cov_500orf_B11, breaks = 5000)

## General Approach

chi_sqared_aproach <- function(contrast, hetdif_th, transdif_th, cov_col){

  ## cat('\n##########################\n')

  ## header <- paste0(
  ##   'Calculating Chi-Sqared for: ',
  ##   cov_col, ' in ',
  ##   contrast[1], ' vs ', contrast[2]
  ## )
  ## print(header)

  trans_col <- paste0(contrast[1], '-', contrast[2], '_MaxVal')

  hetdif_df <- positive_cov_df %>%
    mutate(het_dif = get(paste0(cov_col, contrast[1])) - get(paste0(cov_col, contrast[2]))) %>%
    select(Gene_id, het_dif)

  transdif_df <- trans_df %>%
    rename(trans_dif = matches(trans_col)) %>%
    select(Gene_id, trans_dif)

  final_df <- hetdif_df %>%
    inner_join(transdif_df) %>%
    filter(complete.cases(.)) %>%
    #filter(abs(het_dif) > hetdif_th) %>%
    filter(abs(trans_dif) > transdif_th) %>%
    mutate(logical_het = ifelse(het_dif > 0, 'Positive', 'Negative')) %>%
    mutate(logical_het = ifelse(abs(het_dif) > hetdif_th, logical_het, 'Neutral')) %>%
    mutate(logical_trans = ifelse(trans_dif > 0, 'Positive', 'Negative'))
    #mutate(logical_trans = ifelse(abs(trans_dif) > transdif_th, logical_trans, 'Neutral'))

  ngenes <- dim(final_df)[1]
  ## print(paste0('Number of genes: ', as.character(ngenes)))

  cont_table <- table(
    final_df$logical_het,
    final_df$logical_trans,
    dnn = c('Het_dif', 'Trans_dif')
  )

  ## print(cont_table)

  chisq <- chisq.test(cont_table)
  ## print(chisq)

  frac_stat <- (cont_table[1,2] + cont_table[3,1]) / sum(cont_table)
  ## print(paste0('Custom fraction stat: ', as.character(frac_stat)))
  ## cat('\n##########################\n\n\n\n')

  return(list(
    pval = chisq$p.value,
    frac_stat = frac_stat,
    ngenes = ngenes,
    ##table = cont_table
    het_neg_trans_neg = cont_table[1,1],
    het_neut_trans_neg = cont_table[2,1],
    het_pos_trans_neg = cont_table[3,1],
    het_neg_trans_pos = cont_table[1,2],
    het_neut_trans_pos = cont_table[2,2],
    het_pos_trans_pos = cont_table[3,2]
    ))
}

## Manual 1 contrast 1 coverage column

contrast <- c('12B', '10G')
hetdif_th <- 0.5
transdif_th <- 1
cov_col <- 'Cov_500orf_'

chi_sqared_aproach(contrast, hetdif_th, transdif_th, cov_col)

cov_cols <- unique(sapply(colnames(positive_cov_df), function(x) gsub('12B|10G|A7|E5|B11', '', x)))
cov_cols <- cov_cols[2:length(cov_cols)]

## Create tables by strain

hetdif_th <- 0
transdif_th <- 1

contrasts <- list(
  c('12B', '10G'),
  c('A7', 'E5'),
  c('A7', 'B11'),
  c('E5', 'B11')
)

for (contrast in contrasts){

  xsq_l <- lapply(cov_cols, function(x) chi_sqared_aproach(contrast, hetdif_th, transdif_th, x))
  xsq_df <- tibble(xsq_l) %>%
    unnest_wider(xsq_l) %>%
    mutate(Cov_col = cov_cols) %>%
    select(Cov_col, everything())

  str_contrast <- paste(contrast[1], contrast[2], sep = '_')
  outname <- paste0(
    './Contingency_Table_Corr/xsquared_',
    str_contrast,
    '_transth_', as.character(transdif_th),
    '_hetth_', as.character(hetdif_th),
    '.tsv')
  write_tsv(xsq_df, outname)
}

xsq_df <- as.data.frame(t(xsq_df))
xsq_df['Cov_Col'] <- rownames(xsq_df)
as_tibble(xsq_df)

## Create Table for Pval and Frac_Stat for all strains together
hetdif_th <- 0.5
transdif_th <- 1

pvals_12B_10G <- sapply(cov_cols, function(x) chi_sqared_aproach(c('12B', '10G'), hetdif_th, transdif_th, x)$pval)
fracstats_12B_10G <- sapply(cov_cols, function(x) chi_sqared_aproach(c('12B', '10G'), hetdif_th, transdif_th, x)$frac_stat)

pvals_A7_E5 <- sapply(cov_cols, function(x) chi_sqared_aproach(c('A7', 'E5'), hetdif_th, transdif_th, x)$pval)
fracstats_A7_E5 <- sapply(cov_cols, function(x) chi_sqared_aproach(c('A7', 'E5'), hetdif_th, transdif_th, x)$frac_stat)

pvals_A7_B11 <- sapply(cov_cols, function(x) chi_sqared_aproach(c('A7', 'B11'), hetdif_th, transdif_th, x)$pval)
fracstats_A7_B11 <- sapply(cov_cols, function(x) chi_sqared_aproach(c('A7', 'B11'), hetdif_th, transdif_th, x)$frac_stat)

pvals_E5_B11 <- sapply(cov_cols, function(x) chi_sqared_aproach(c('E5', 'B11'), hetdif_th, transdif_th, x)$pval)
fracstats_E5_B11 <- sapply(cov_cols, function(x) chi_sqared_aproach(c('E5', 'B11'), hetdif_th, transdif_th, x)$frac_stat)

contingency_results <- tibble(
  Cov_Col = names(pvals_12B_10G),
  Pvals_12B_10G = pvals_12B_10G,
  Pvals_A7_E5 = pvals_A7_E5,
  Pvals_A7_B11 = pvals_A7_B11,
  Pvals_E5_B11 = pvals_E5_B11,

  Frac_12B_10G = fracstats_12B_10G,
  Frac_A7_E5 = fracstats_A7_E5,
  Frac_A7_B11 = fracstats_A7_B11,
  Frac_E5_B11 = fracstats_E5_B11) %>%
  print(n = 2000)

write_tsv(contingency_results, './Contingency_Table_Corr/contingency_aproach_results.tsv')

#### Heatmap Funtions ####

myHeatmap <- function(df, family_facet){
  p <- ggplot(df, aes(x = variable, y = Label, fill = value)) +
    geom_tile(colour="snow3") +
    theme(
      text=element_text(size=24),

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
  if (family_facet){p <- p+facet_grid(Family ~., scales = "free_y", space = "free")}
  return(p)
}
hetHeatmap <- function(df, family_facet){
  p <- myHeatmap(df, family_facet)
  p <- p + scale_fill_gradient(low = "white",
                               high = "orange",
                               na.value="white",
                               limits = c(0,NA)) +
    scale_y_discrete(position = "right") +
    theme(
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
hetDifHeatmap <- function(df, family_facet){
  p <- myHeatmap(df, family_facet)
  p <- p + scale_fill_gradient2(low = "chartreuse3",
                                mid = "white",
                                high = "darkred",
                                na.value="grey90") +
    scale_y_discrete(position = "left") +
    theme(
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
transHeatmap <- function(df, family_facet){
  p <- myHeatmap(df, family_facet)
  p <- p + scale_fill_gradient2(low = "#536DFE",
                                mid = "white",
                                high = "yellow2",
                                na.value="grey90") +
    scale_y_discrete(position = "left") +
    theme(
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

family_heatmap <- function(mdf){
  p <- myHeatmap(mdf, family_facet)
  ##p <- p + geom_text(aes(label=Label))
  p <- p + my_scale
  p <- p + scale_y_discrete(limits=(rev(levels(mdf$Label))))
  p <- p + theme(
             axis.text.x = element_blank(),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank(),

             panel.border=element_blank(),
             panel.grid.major=element_blank(),

             strip.background = element_blank(),
             strip.text.x = element_blank(),
             strip.text.y = element_blank()
           )
  return(p)
}

#### A first plot ####

trans_het_df <- trans_df %>%
  left_join(het_df, by = 'Gene_id') %>%
  left_join(info_df, by = 'Gene_id')

e5_b11 <- trans_het_df %>%
  mutate(Het_dif = Het_E5 - Het_B11) %>%
  filter(abs(`E5-B11_MaxVal`) > 1 & abs(Het_dif) > 0) %>%
  arrange(desc(`E5-B11_MaxVal`)) %>%
  mutate(Label = factor(Label, levels = Label))

df_e5_b11_trans <- e5_b11 %>%
  select(Gene_id, `E5-B11_MaxVal`, Annot, Family, SubFamily, Label)

df_e5_b11_het <- e5_b11 %>%
  select(Gene_id, Het_dif, Het_B11, Het_E5, Annot, Family, SubFamily, Label)


mdf_e5_b11_het <- melt(df_e5_b11_het, id.vars = c('Gene_id', 'Annot', 'Family', 'SubFamily', 'Label'))
mdf_e5_b11_trans <- melt(df_e5_b11_trans, id.vars = c('Gene_id', 'Annot', 'Family', 'SubFamily', 'Label'))

family_facet <- FALSE

het_plot <- hetHeatmap(mdf_e5_b11_het, family_facet)
trans_plot <- transHeatmap(mdf_e5_b11_trans, family_facet)
all_plot <- grid.arrange(trans_plot, het_plot, nrow = 1, widths = c(1,3))

#### General Plots ####

all_transcov_df <- trans_df %>%
  left_join(positive_cov_df, by = 'Gene_id') %>%
  left_join(info_df, by = 'Gene_id')

finalHeatmap <- function(contrast, het_col, abs_trans_filter, abs_het_filter, family, family_facet) {

  trans_col <- paste0(contrast[1], '-', contrast[2], '_MaxVal')
  het_col_1 <- paste0(het_col, '_', contrast[1])
  het_col_2 <- paste0(het_col, '_', contrast[2])

  subset_all <- all_transcov_df %>%
    mutate(Het_dif = get(het_col_1) - get(het_col_2)) %>%
    filter(abs(get(trans_col)) > abs_trans_filter &
           (abs(Het_dif) > abs_het_filter | is.na(Het_dif))) %>%
    select(Gene_id,
           matches(trans_col),
           matches(het_col_1),
           matches(het_col_2),
           Het_dif,
           Label, Name, Annot, Family, SubFamily
           )

  ## Subset by family is needed
  if (!is.na(family)){
    subset_all <- subset_all %>%
      filter(Family == family)
  }

  ## Ordering
  if (sorting == 'clust') {
    mtx <- subset_all %>%
      select(where(is.numeric))

    ## Make hierarquical Clustering
    dmtx <- dist(scale(mtx), method = "euclidean")
    cl <- hclust(dmtx, method = 'average')
    subset_all$Label <- factor(subset_all$Label, levels = subset_all$Label[cl$order])

  } else if (sorting == 'trans') {
    subset_all <- subset_all %>%
      arrange(desc(get(trans_col))) %>%
      mutate(Label = factor(Label, levels = Label))

  } else if (sorting == 'het') {
    subset_all <- subset_all %>%
      arrange(desc(Het_dif)) %>%
      mutate(Label = factor(Label, levels = Label))
  }

  subset_trans <- subset_all %>%
    select(Gene_id, matches(trans_col), Label, Name, Annot, Family, SubFamily)

  subset_het <- subset_all %>%
    select(Gene_id, matches(het_col_1), matches(het_col_2), Label, Name, Annot, Family, SubFamily)

  subset_het %>%
    print(width = 200)

  subset_hetDif <- subset_all %>%
    select(Gene_id, Het_dif, Label, Name, Annot, Family, SubFamily)

  melt_vars <- c('Gene_id', 'Label', 'Name', 'Annot', 'Family', 'SubFamily')

  mdf_het <- melt(subset_het, id.vars = melt_vars)
  mdf_trans <- melt(subset_trans, id.vars = melt_vars)
  mdf_hetDif <- melt(subset_hetDif, id.vars = melt_vars)

  het_plot <- hetHeatmap(mdf_het, family_facet)
  trans_plot <- transHeatmap(mdf_trans, family_facet)
  hetDif_plot <- hetDifHeatmap(mdf_hetDif, family_facet)
  all_plot <- grid.arrange(trans_plot, hetDif_plot, het_plot, nrow = 1, widths = c(1,1,3))

  outname <- paste0('./Plots/Trans_Het_by_Strain/',
                    contrast[1], '_',
                    contrast[2], '_',
                    het_col, '_',
                    'transth_', as.character(abs_trans_filter),
                    '.pdf')
  print(outname)
  ggsave(outname, all_plot, device = 'pdf')
}

colnames(cor_cov_df)

contrast <- c('12B', '10G')
het_col <- 'Cov_500fp_manybins'

abs_trans_filter <- 1
abs_het_filter <- 0

family_facet <- F
family <- NA
## trans, het or clust
sorting <- 'trans'

plt <- finalHeatmap(contrast, het_col, abs_trans_filter, abs_het_filter, family, family_facet)

analysis_df <- trans_df %>%
  left_join(info_df, by = 'Gene_id') %>%
  left_join(cor_cov_df %>% select(Gene_id, contains('Cov_500fp_manybins')), by = 'Gene_id')

## analysis_df %>%
##   filter(abs(`A7-E5_MaxVal`) > 2) %>%
##   select(Gene_id, `A7-E5_MaxVal`)

## analysis_df %>%
##   filter(Gene_id == 'PF3D7_1372600') %>%
##   select(Gene_id, `A7-E5_MaxVal`)

##### Plotting part #####
contrast <- c('E5', 'B11')
het_col <- 'Cov_500fp_manybins'

trans_col <- paste0(contrast[1], '-', contrast[2], '_MaxVal')
het_col_1 <- paste0(het_col, '_', contrast[1])
het_col_2 <- paste0(het_col, '_', contrast[2])

infocols <- colnames(info_df)

x <- analysis_df %>%
  select(Gene_id,
         Family_Grouped,
         matches(trans_col),
         matches(het_col_1),
         matches(het_col_2),
         one_of(infocols)
         ) %>%
  filter(abs(get(trans_col)) > 1)

table(x$Gam_specific)
table(x$Family_Grouped)
table(x$Variant)

## Test for Gam Genes enrichment
count_mtx <- matrix(
  c(sum(x$Gam_specific),
    sum(info_df$Gam_specific),
    sum(!x$Gam_specific),
    sum(!info_df$Gam_specific)),
  nrow=2,ncol=2)

fisher.test(count_mtx, alternative="greater")

## Donut Plot

contrasts_donuts <- function(contrast, het_col, df){
  trans_col <- paste0(contrast[1], '-', contrast[2], '_MaxVal')
  het_col_1 <- paste0(het_col, '_', contrast[1])
  het_col_2 <- paste0(het_col, '_', contrast[2])

  infocols <- colnames(info_df)

  x <- df %>%
    select(Gene_id,
           Family_Grouped,
           matches(trans_col),
           matches(het_col_1),
           matches(het_col_2),
           one_of(infocols)
           ) %>%
    filter(abs(get(trans_col)) > 1)
  data <- as.data.frame(table(x$Family_Grouped))
  colnames(data) <- c('category', 'count')

  ## Compute percentages
  data$fraction = data$count / sum(data$count)

  ## Compute the cumulative percentages (top of each rectangle)
  data$ymax = cumsum(data$fraction)

  ## Compute the bottom of each rectangle
  data$ymin = c(0, head(data$ymax, n=-1))

  ## Compute label position
  data$labelPosition <- (data$ymax + data$ymin) / 2

  ## Compute a good label
  data$label <- paste0(data$category, "\n value: ", data$count)

  ## Make the plot
  p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect(color="white") +
                                        #geom_label( x=4.2, aes(y=labelPosition, label=label), size=3) +
    coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
    xlim(c(2, 4)) +# Try to remove that to see how to make a pie chart
    theme_void() + my_scale#scale_fill_viridis(discrete=TRUE)

  print(p)
  outname <- paste0('./Plots/Donuts/',
                    contrast[1], '_',
                    contrast[2], '_',
                    'families',
                    '.svg')

  ggsave(outname, p, device = 'svg')
}

contrasts <- list(
  c('12B', '10G'),
  c('A7', 'E5'),
  c('A7', 'B11'),
  c('E5', 'B11')
  )
het_col <- 'Cov_500fp_manybins'
contrast <-   c('12B', '10G')
df <- analysis_df

for(contrast in contrasts){
  print(contrast)
  contrasts_donuts(contrast, het_col, analysis_df)
}

info_df %>%
  filter(Family_Grouped == 'Not CVGs')

table(info_df$Variant)

info_df %>%
  filter(!Variant) %>%
  filter(Family_Grouped != 'Not CVGs')

#### 5/ORF/3 Plots ####

finalHeatmap <- function(trans_col, het_col_1, het_col_2, abs_trans_filter, abs_het_filter, family, family_facet) {

  subset_all <- trans_het_df %>%
    filter(abs(get(trans_col)) > abs_trans_filter) %>%
    arrange(desc(get(trans_col))) %>%
    mutate(Gene_id = factor(Gene_id, levels = Gene_id))

  if (!is.na(family)){
    subset_all <- subset_all %>%
      filter(Family == family)
  }

  subset_trans <- subset_all %>%
    select(Gene_id, Label, matches(trans_col), Annot, Family, SubFamily)

  subset_het <- cov_5ORF3_df %>%
    left_join(info_df %>% select(Gene_id), by = 'Gene_id') %>%
    select(Gene_id, contains(het_col_1), contains(het_col_2))

  hetcol_names <- c(paste(rep(het_col_1, 3), c('5p', 'ORF', '3p'), sep = '_'),
                    paste(rep(het_col_2, 3), c('5p', 'ORF', '3p'), sep = '_'))

  subset_het <- subset_het %>% setNames(c('Gene_id', hetcol_names))

  subset_het <- subset_trans %>%
    left_join(subset_het, by = 'Gene_id') %>%
    select(-matches(trans_col)) %>%
    mutate(Label = factor(Label, levels = Label))

  table(is.na(subset_het$Label))
  subset_het %>%
    filter(is.na(Label))

  info_df %>%
    filter(Gene_id == 'PF3D7_0112400')


  cols_5p <- colnames(subset_het %>% select(contains('5p')))
  cols_ORF <- colnames(subset_het %>% select(contains('ORF')))
  cols_3p <- colnames(subset_het %>% select(contains('3p')))

  subset_hetDif <- subset_het %>%
    mutate(Dif_5p = get(cols_5p[1]) - get(cols_5p[2])) %>%
    mutate(Dif_ORF = get(cols_ORF[1]) - get(cols_ORF[2])) %>%
    mutate(Dif_3p = get(cols_3p[1]) - get(cols_3p[2])) %>%
    select(Gene_id, Label, Annot, Family, SubFamily, contains('Dif'))

  melt_vars <- c('Gene_id', 'Label', 'Annot', 'Family', 'SubFamily')

  mdf_het <- melt(subset_het, id.vars = melt_vars)
  mdf_trans <- melt(subset_trans, id.vars = melt_vars)
  mdf_hetDif <- melt(subset_hetDif, id.vars = melt_vars)

  het_plot <- hetHeatmap(mdf_het, family_facet)
  trans_plot <- transHeatmap(mdf_trans, family_facet)
  hetDif_plot <- hetDifHeatmap(mdf_hetDif, family_facet)
  all_plot <- grid.arrange(trans_plot, hetDif_plot, het_plot, nrow = 1, widths = c(0.7,1,3))
}

trans_col <- 'A7-E5_MaxVal'
het_col_1 <- 'A7'
het_col_2 <- 'E5'

abs_trans_filter <- 1.5
abs_het_filter <- 0

family_facet <- T

family <- NA


finalHeatmap(trans_col, het_col_1, het_col_2, abs_trans_filter, abs_het_filter, family, family_facet)

old_arrays <- '../Microarrays/New_Old_separate_approach/Old_Arrays/R_results_OldArrays_Variantome/'

filtered_12B_10G <- read_tsv(paste0(old_arrays, '12B_10G_final_df.tsv')) %>%
  filter(PassRed & PassDuplDel)
filtered_12B_3D7B <- read_tsv(paste0(old_arrays, '12B_3D7B_final_df.tsv')) %>%
  filter(PassRed & PassDuplDel)
filtered_10G_3D7B <- read_tsv(paste0(old_arrays, '10G_3D7B_final_df.tsv')) %>%
  filter(PassRed & PassDuplDel)

new_arrays <- '../Microarrays/New_Old_separate_approach/New_Arrays/R_results_NewArray/'

filtered_A7_B11 <- read_tsv(paste0(new_arrays, 'A7_B11_final_df.tsv')) %>%
  filter(PassRed & PassDuplDel)
filtered_A7_E5 <- read_tsv(paste0(new_arrays, 'A7_E5_final_df.tsv')) %>%
  filter(PassRed & PassDuplDel)
filtered_B11_E5 <- read_tsv(paste0(new_arrays, 'B11_E5_final_df.tsv')) %>%
  filter(PassRed & PassDuplDel)

## Swich B11vsE5 for E5vsB11
filtered_E5_B11 <- filtered_B11_E5 %>%
  rename(`E5-B11_MaxVal` = `B11-E5_MaxVal`, `E5-B11_MaxTime` = `B11-E5_MaxTime`) %>%
  mutate(`E5-B11_MaxVal` = -`E5-B11_MaxVal`)

trans_df %>%
  filter(!Gene_id %in% filtered_12B_10G$Gene_id) %>%
  filter(!Gene_id %in% filtered_12B_3D7B$Gene_id) %>%
  filter(!Gene_id %in% filtered_10G_3D7B$Gene_id) %>%
  filter(!Gene_id %in% filtered_A7_B11$Gene_id) %>%
  filter(!Gene_id %in% filtered_A7_E5$Gene_id) %>%
  filter(!Gene_id %in% filtered_E5_B11$Gene_id) %>%
  write_tsv('filtered_out_genes_by_Red_DuplDel.tsv')

## Select MaxDif for each gene (once filtered by redfilter and dupl_del)
maxDif_df <- tibble()
for (gid in trans_df$Gene_id){

  ## Create a list with each contrast difference df
  difs <- list(
    filtered_12B_10G,
    filtered_12B_3D7B,
    filtered_10G_3D7B,
    filtered_A7_B11,
    filtered_A7_E5,
    filtered_E5_B11
  )

  ## Filter by Gene_id
  dif_dfs <- lapply(difs, function(x) x %>% filter(Gene_id == gid))

  ## Function to get MaxVal of each df and join them in a vector
  get_MaxVal <- function(x){
    maxVal <- x %>%
      select(contains('MaxVal')) %>%
      pull()
    if (identical(maxVal, numeric(0))) {maxVal <- NA}
    return(maxVal)
  }

  maxVect <- sapply(dif_dfs, get_MaxVal)
  max_idx <- which.max(abs(maxVect))

  ## Handle genes that don't pass filters (set to NA)
  if (identical(maxVect, rep(NA, 6))) {
    out_row <- tibble(
      Gene_id = gid,
      Max_aMAFC = NA,
      Max_Time = NA,
      On_trans = NA,
      Off_trans = NA,
      PassMaxtime = FALSE
    )
    ## Handle rest of genes
  } else {
    ## Get on/off strain names
    maxDif <- dif_dfs[[max_idx]] %>% select(contains('_MaxVal'))
    maxVal <- maxDif %>% pull()

    maxStrains <- colnames(maxDif)
    strains <- strsplit(maxStrains, split = '_', fixed = T)[[1]][1]
    strains <- strsplit(strains, split = '-', fixed = T)[[1]]
    if (maxVal >= 0){
      On <- strains[1]
      Off <- strains[2]
    } else {
      On <- strains[2]
      Off <- strains[1]
    }

    ## Collect MaxDif row from the appropiate df
    out_row <- dif_dfs[[max_idx]] %>%
      select(Gene_id, contains('_MaxVal'), contains('_MaxTime'), PassMaxtime) %>%
      setNames(c('Gene_id', 'Max_aMAFC', 'Max_Time', 'PassMaxtime')) %>%
      mutate(
        On_trans = On,
        Off_trans = Off
      )

  }
  maxDif_df <- bind_rows(maxDif_df, out_row)
}

maxDif_df <- maxDif_df %>%
  arrange(-abs(Max_aMAFC)) %>%
  left_join(trans_df, by = 'Gene_id')

maxDif_df %>%
  print(width = 400)

pre3D7_substitution <- maxDif_df
new_array_areas <- read_csv('/mnt/Disc4T/Projects/PhD_Project/Microarrays/New_Old_separate_approach/New_Arrays/R_results_NewArray/area_geneLevel.csv') %>%
  select(-Name, -Annot, -Variant)

max_trans_newarray <- maxDif_df %>%
  left_join(new_array_areas)

get_new_onoff <- function(gid) {

  #gid <- 'PF3D7_0421500'
  on <- max_trans_newarray %>%
    filter(Gene_id == gid) %>%
    select(On_trans) %>%
    pull()

  off <- max_trans_newarray %>%
    filter(Gene_id == gid) %>%
    select(Off_trans) %>%
    pull()

  time <- max_trans_newarray %>%
    filter(Gene_id == gid) %>%
    select(Max_Time) %>%
    pull()

  onoff <- TRUE
  if (is.na(on) | is.na(off)){
    onoff <- FALSE
  } else {
    if (on == '3D7B') {func <- which.max}
    if (off == '3D7B') {func <- which.min}
    if (on != '3D7B' & off != '3D7B') {onoff <- FALSE}
  }

  if (is.na(time)){
    out <- NA
  } else if (!onoff){
    out <- NA
  } else {
    strains <- c('A7', 'B11', 'E5')

    vect <- max_trans_newarray %>%
      filter(Gene_id == gid) %>%
      select(contains(time))

    ## Filter by red percent if 3D7B is 'On' strain
    red_pcnts <- red_df %>%
      filter(Gene_id == gid) %>%
      select(A7, B11, E5)

    red_mask <- red_pcnts > 15
    if (on != '3D7B'){red_mask <- c(TRUE, TRUE, TRUE)}

    ## Filter by dupl/del
    dupl_del_mask <- c(
      !gid %in% dupl_del$A7K9$Gene_id,
      !gid %in% dupl_del$B11$Gene_id,
      !gid %in% dupl_del$E5K9$Gene_id
    )

    ## Apply both filters
    whole_mask <- red_mask & dupl_del_mask

    vect <- vect[whole_mask]
    strains <- strains[whole_mask]

    ## Final output
    ifelse(all(is.na(vect)) | !any(whole_mask), out <- NA, out <- strains[func(vect)])
  }
  return(out)
}

newonoffs <- sapply(max_trans_newarray$Gene_id, get_new_onoff)
maxDif_df['New_OnOffs'] <- newonoffs
max_trans_newarray['New_OnOffs'] <- newonoffs

max_trans_newarray %>%
  select(Gene_id, On_trans, Off_trans, New_OnOffs) %>%
  print(n = 50)


maxDif_df <- maxDif_df %>%
  mutate(On_trans = ifelse(On_trans == '3D7B',
                           New_OnOffs,
                           On_trans)) %>%
  mutate(Off_trans = ifelse(Off_trans == '3D7B',
                           New_OnOffs,
                           Off_trans)) %>%
  mutate(Is_3D7B = !is.na(New_OnOffs)) %>%
  select(-New_OnOffs)

maxDif_df %>%
  select(-contains('_MaxVal'), -contains('_MaxTime'))

#old_maxTrans_df <- maxTrans_df
maxTrans_df <- maxDif_df

## Remove tRNAs
transFC_df <- maxTrans_df %>%
  left_join(info_df, by = 'Gene_id') %>%
  filter(!Is_tRNA)

## Add red percentile
get_red_percent_onstrain <- function(gid){
  on_strain <- transFC_df %>%
    filter(Gene_id == gid) %>%
    select(On_trans) %>%
    pull()

  if (is.na(on_strain) | !gid %in% red_df$Gene_id){
    pcnt <- NA
  } else {
    pcnt <- red_df %>%
      filter(Gene_id == gid) %>%
      select(matches(on_strain)) %>%
      pull()
  }
  return(pcnt)
}

percs <- sapply(transFC_df$Gene_id, get_red_percent_onstrain)

transFC_df <- transFC_df %>%
  mutate(Red_Pcnt_Onstrain = percs)

transFC_df %>%
  filter(abs(Max_aMAFC) > 2) %>%
  filter(Red_Pcnt_Onstrain < 15) %>%
  select(Gene_id, Label, Max_aMAFC, On_trans, Off_trans, Red_Pcnt_Onstrain, Variant) %>%
  write_csv('max_red_percentile_15_filter_failed.csv')

transFC_df %>%
  filter(abs(Max_aMAFC) > 2) %>%
  filter(Red_Pcnt_Onstrain < 15) %>%
  select(Gene_id, Label, Max_aMAFC, On_trans, Off_trans, Red_Pcnt_Onstrain, Variant)


transFC_df %>%
  mutate(Red_filter = case_when(
           abs(Max_aMAFC) >= 1 & Red_Pcnt_Onstrain > 20 ~ 'Trans_Red',
           abs(Max_aMAFC) >= 1 & Red_Pcnt_Onstrain <= 20 ~ 'Trans_NoRed',
           abs(Max_aMAFC) < 1 & Red_Pcnt_Onstrain > 20 ~ 'NoTrans_Red',
           abs(Max_aMAFC) < 1 & Red_Pcnt_Onstrain <= 20 ~ 'NoTrans_NoRed',
           )) %>%
  count(Red_filter)

## New approach
## Check which areas does maxtimepoint overlapp -> check if aAFC > th at this areas

point_overlap <- function(point, interval){
  point >= interval[1] & point <= interval[2]
}

transFC_df_noNAs <- transFC_df %>%
  filter(!is.na(On_trans) & !is.na(Off_trans))

aFC_in_maxtime <- c()
gids <- c()
th <- 2
for (gid in transFC_df_noNAs$Gene_id){

  #gid <- 'PF3D7_0900200'
  ## Get max contrast
  strains <- max_trans_newarray %>% ## this is the one where 3D7B appears as such
    filter(Gene_id == gid) %>%
    select(On_trans, Off_trans) %>%
    as.character()

  if (!any(is.na(strains))) {
    ## Check wehter its old/new df
    batch <- ifelse(any(c('A7', 'E5', 'B11') %in% strains), 'New', 'Old')

    ## Create time-regions
    breaks <- breaks_df %>%
      select(contains(batch)) %>%
      pull()

    left <- c(breaks[1], breaks[3])
    right <- c(breaks[3], breaks[5])
    mid <- c(breaks[2], breaks[4])
    sides_l <- c(breaks[1], breaks[2])
    sides_r <- c(breaks[4], breaks[5])

    ## Get maxtime
    maxtime <- maxtime_df %>%
      filter(Gene_id == gid) %>%
      select(contains(batch)) %>%
      pull()

    ## Ensure maxtime is in the areas intervals
    if (maxtime < breaks[1]) {maxtime <- breaks[1]}
    if (maxtime > breaks[5]) {maxtime <- breaks[5]}

    ## Get overlappped areas
    areas <- list('left' = left, 'right' = right, 'mid' = mid,
                  'sides' = sides_l, 'sides' = sides_r)
    overlaps <- sapply(areas, function(x) point_overlap(maxtime, x))

    ## Get aAFC in overlapping areas
    aFCs <- arees_df %>%
      filter(Gene_id == gid) %>%
      select(contains(names(areas[overlaps])) &
             contains(strains[1]) &
             contains(strains[2])) %>%
      as.numeric()

    aFC_in_maxtime <- c(aFC_in_maxtime, any(abs(aFCs) > th))
    gids <- c(gids, gid)
  } else {
    aFC_in_maxtime <- c(aFC_in_maxtime, NA)
  }
}

maxtime_aAFC_df <- tibble(Gene_id = gids, aAFC_at_maxtime = aFC_in_maxtime)

transFC_df_noNAs <- transFC_df_noNAs %>%
  left_join(maxtime_aAFC_df)

transFC_df_noNAs %>%
  filter(abs(Max_aMAFC) > 2) %>%
  filter(!aAFC_at_maxtime) %>%
  select(Label, Max_aMAFC, On_trans, Off_trans, Annot) %>%
  write_csv('areaFC_in_maxtp_filter_failed.csv')

transFC_df_noNAs %>%
  filter(abs(Max_aMAFC) > 2) %>%
  filter(!aAFC_at_maxtime) %>%
  select(Label, Max_aMAFC, On_trans, Off_trans, Annot)

transFC_df_noNAs

check_dupldel <- function(gid){

  #gid <- 'PF3D7_0831700'
  strains <- transFC_df_noNAs %>%
    filter(Gene_id == gid) %>%
    select(On_trans, Off_trans) %>%
    as.character()
  #print(strains)

                                        #f_path <- '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/Duplication_Deletion_Regions/Crossed_with_genes/'
  ## file_list <- list.files(path=f_path)
  ## genes_tsv <- file_list[sapply(file_list, function(x) grepl(strains[1], x) & grepl(strains[2], x))]
  ## #print(genes_tsv)
  ## dd_genes <- read_tsv(paste0(f_path, genes_tsv), col_names = F, col_types = cols())
  ## gid %in% dd_genes$X1

  gid %in% dupl_del[[strains[1]]]$Gene_id | gid %in% dupl_del[[strains[2]]]$Gene_id
}


dupl_del <- sapply(transFC_df_noNAs$Gene_id, check_dupldel)
transFC_df_noNAs['Dupl_Del_in_Contrast'] <- dupl_del


transFC_df_noNAs %>%
  filter(Dupl_Del_in_Contrast) %>%
  write_csv('genes_with_dupl_del.csv')

transFC_df_noNAs %>%
  filter(Dupl_Del_in_Contrast) %>%
  select(Gene_id, Max_aMAFC, On_trans, Off_trans, Dupl_Del_in_Contrast) %>%
  print(n = 100)

filtered_10G_3D7B %>%
  filter(Gene_id == 'PF3D7_0935400')

#### Check CGH data matches out input derived dupl/del data ####

## Load old_array for translation
old <- read_csv('../Microarrays/New_Old_separate_approach/Old_Arrays/gene_level_raw_table.csv') %>%
  select(Old_id, Gene_id)

cgh <- read_xls('../Microarrays/New_Old_separate_approach/Old_Arrays/Variantome_Original/3D7_Variantome_AllData_withGam.xls', sheet = 4)

cgh <- cgh %>%
  select(
    `...1`,
    X1.2b,
    X10g,
    X3d7b
  ) %>%
  rename(
    Old_id = `...1`,
    CGH_12B = X1.2b,
    CGH_10G = X10g,
    CGH_3D7B = X3d7b
  ) %>%
  mutate(across(contains('CGH'), as.numeric)) %>%
  left_join(old) %>%
  relocate(Gene_id, Old_id) %>%
  filter(!is.na(Gene_id))

cgh_th <- 0.7
cgh <- cgh %>%
  mutate(CGH_12B_10G = abs(CGH_12B - CGH_10G) > cgh_th) %>%
  mutate(CGH_12B_3D7B = abs(CGH_12B - CGH_3D7B) > cgh_th) %>%
  mutate(CGH_10G_3D7B = abs(CGH_10G - CGH_3D7B) > cgh_th)

f_path <- '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/Duplication_Deletion_Regions/Crossed_with_genes/'
file_list <- c(
  "12B_minus_10G_100bp_500smth_RPKM_cov_norm_pdf_0999999_minlen500_genes.tsv",
  "12B_minus_3D7_100bp_500smth_RPKM_cov_norm_pdf_0999999_minlen500_genes.tsv",
  "10G_minus_3D7_100bp_500smth_RPKM_cov_norm_pdf_0999999_minlen500_genes.tsv"
)

dupl_dl_12B_10G <- read_tsv(paste0(f_path, file_list[1]), col_names = F) %>%
  select(X1) %>% pull()

dupl_dl_12B_3D7B <- read_tsv(paste0(f_path, file_list[2]), col_names = F) %>%
  select(X1) %>% pull()

dupl_dl_10G_3D7B <- read_tsv(paste0(f_path, file_list[3]), col_names = F) %>%
  select(X1) %>% pull()

## 3D7B no és el que toca!!!!

cgh %>%
  filter(CGH_12B_10G) %>%
  select(Gene_id) %>%
  pull() %in% dupl_dl_12B_10G


cgh <- cgh %>%
  mutate(In_Dupl_Del_12B_10G = Gene_id %in% dupl_dl_12B_10G) %>%
  mutate(In_Dupl_Del_12B_3D7B = Gene_id %in% dupl_dl_12B_3D7B) %>%
  mutate(In_Dupl_Del_10G_3D7B = Gene_id %in% dupl_dl_10G_3D7B)

cgh %>%
  select(
    Gene_id, Old_id,
    CGH_12B, CGH_10G,
    CGH_12B_10G, In_Dupl_Del_12B_10G
  ) %>%
  write_csv('cgh_inputs_12B_10G.csv')

## cgh %>%
##   select(
##     Gene_id, Old_id,
##     CGH_12B, CGH_10G, CGH_3D7B,
##     CGH_12B_3D7B, In_Dupl_Del_12B_3D7B
##   ) %>%
##   filter(CGH_12B_3D7B) %>%
##   write_csv('cgh_inputs_12B_3D7B.csv')

## cgh %>%
##   select(
##     Gene_id, Old_id,
##     CGH_12B, CGH_10G, CGH_3D7B,
##     CGH_10G_3D7B, In_Dupl_Del_10G_3D7B
##   ) %>%
##   filter(CGH_10G_3D7B) %>%
##   write_csv('cgh_inputs_10G_3D7B.csv')

finalFC <- transFC_df_noNAs %>%
  filter(abs(Max_aMAFC) > 2) %>%
  filter(Red_Pcnt_Onstrain > 15) %>%
  filter(aAFC_at_maxtime) %>%
  filter(!Dupl_Del_in_Contrast)

write_csv(finalFC, 'max_log2FC2_filters_passed.csv')

## transFC_df_noNAs %>%
##   filter(abs(Max_aMAFC) > 2) %>%
##   filter(Red_Pcnt_Onstrain > 15)

## transFC_df_noNAs %>%
##   select(Red_Pcnt_Onstrain) %>%
##   summary()


## colnames(transFC_df_noNAs)

## transFC_df_noNAs %>%
##   select(aAFC_at_maxtime) %>%
##   summary()

## transFC_df_noNAs %>%
##   filter(abs(Max_aMAFC) > 1) %>%
##   filter(Red_Pcnt_Onstrain > 15) %>%
##   filter(!aAFC_at_maxtime) %>%
##   write_csv('maxFC_lg2FC1_fail_aAFCmaxtime.csv')

## transFC_df_noNAs %>%
##   filter(abs(Max_aMAFC) > 1)

finalFC %>%
  filter(On_trans == 'B11' |
         Off_trans == 'B11')

finalFC %>%
  count(Variant)

## From Peak-Calling
transFC_df_noNAs %>%
  filter(!Gene_id %in% finalFC$Gene_id) %>%
  filter(Difpeaks_12B_10G | Difpeaks_A7_E5 | Difpeaks_A7_B11 | Difpeaks_E5_B11) %>%
  select(
    Gene_id, Max_aMAFC,
    On_trans, Off_trans,
    #contains('Difpeaks'),
    Red_Pcnt_Onstrain,
    aAFC_at_maxtime,
    Dupl_Del_in_Contrast
  ) %>%
  left_join(info_df, by = 'Gene_id') %>%
  write_tsv('genes_with_difpeaks_not_in_FC_tables.tsv')

on_off_donut <- function(col){
  data <- finalFC %>%
    count(get(col)) %>%
    set_names('category', 'count')

  ## Compute percentages
  data$fraction = data$count / sum(data$count)

  ## Compute the cumulative percentages (top of each rectangle)
  data$ymax = cumsum(data$fraction)

  ## Compute the bottom of each rectangle
  data$ymin = c(0, head(data$ymax, n=-1))

  ## Compute label position
  data$labelPosition <- (data$ymax + data$ymin) / 2

  ## Compute a good label
  data$label <- paste0(data$category, "\n value: ", data$count)

  ## Make the plot
  p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category))
  p <- p +  geom_rect(color="white")
  ##p <- p + geom_label( x=4.2, aes(y=labelPosition, label=label), size=3)
  p <- p + coord_polar(theta="y") # Try to remove that to understand how the chart is built initially
  p <- p + xlim(c(2, 4)) # Try to remove that to see how to make a pie chart
  p <- p + theme_void()
  p <- p + scale_fill_viridis(discrete=TRUE)
  p
}

ggsave('../Plots/on_difgenes_donut.pdf', on_off_donut('On_trans'), device = 'pdf')
ggsave('../Plots/off_difgenes_donut.pdf', on_off_donut('Off_trans'), device = 'pdf')

## Select transcription cols

tdf <- transFC_df_noNAs %>%
  select(Gene_id, Max_aMAFC, On_trans, Off_trans, Is_3D7B)

#print(maxTrans_df[63,], width = 400)

## Create empty DF
hetcols <- colnames(cor_cov_df %>% select(contains('12B')))
hetcols <- str_replace(hetcols, '_12B', '')

col_names <- c('Gene_id',
               paste0(hetcols, '_On'),
               paste0(hetcols, '_Off'),
               paste0(hetcols, '_Dif'))

cn <- setNames(rep('', length(col_names)), col_names)
df <- bind_rows(cn)[0,]

## Traverse max trans. df and get coverage cols

## which(tdf$Gene_id == 'PF3D7_0100300')
## i <- 63
## tdf[63,]

for (i in 1:dim(tdf)[1]){

  gid <- as.character(tdf$Gene_id[i])
  on <- tdf$On_trans[i]
  off <- tdf$Off_trans[i]

  ## Select contrast in which difference sign must be switched
  contrast <- paste0(on, '-', off)
  positive_contrasts <- c(
    '12B-10G', '12B-A7', '12B-E5', '12B-B11',
    '10G-A7', '10G-E5', '10G-B11',
    'A7-E5', 'A7-B11',
    'E5-B11'
  )

  ifelse(
    contrast %in% positive_contrasts,
    neg_dif <- FALSE,
    neg_dif <- TRUE
    )

  ## Main Loop
  if (gid %in% cor_cov_df$Gene_id & !is.na(on) & !is.na(off)){

    onvect <- cor_cov_df %>%
      filter(Gene_id == gid) %>%
      select(contains(on))

    offvect <- cor_cov_df %>%
      filter(Gene_id == gid) %>%
      select(contains(off))

    diffvect <- onvect-offvect
    if (neg_dif) {diffvect <- -diffvect}

    row <- c(gid,
             unlist(onvect, use.names = F),
             unlist(offvect, use.names = F),
             unlist(diffvect, use.names = F))
    row <- setNames(row, col_names)
    df <- df %>% add_row(bind_rows(row))
  }
}

signed_maxtrans_cov <- df %>%
  mutate(across(-Gene_id,  as.numeric)) %>%
  full_join(tdf, by = 'Gene_id') %>%
  #mutate(Max_aMAFC = abs(Max_aMAFC)) %>%
  left_join(info_df)

## genes_nocov <- unsigned_maxtrans_cov %>%
##      select(Gene_id, contains('Cov_500fp_manybins')) %>%
##      filter(!complete.cases(.)) %>%
##      select(Gene_id) %>% pull()

## maxTrans_df %>%
##   filter(is.na(On_trans) | is.na(Off_trans))

## cor_cov_df %>%
##   select(Gene_id, contains('Cov_500fp_manybins')) %>%
##   filter(!complete.cases(.)) %>%
##   print(width = 400)

## Select transcription cols

tdf <- transFC_df_noNAs %>%
  select(Gene_id, Max_aMAFC, On_trans, Off_trans, Is_3D7B)

#print(maxTrans_df[63,], width = 400)

## Create empty DF
hetcols <- colnames(cor_cov_df %>% select(contains('12B')))
hetcols <- str_replace(hetcols, '_12B', '')

col_names <- c('Gene_id',
               paste0(hetcols, '_On'),
               paste0(hetcols, '_Off'),
               paste0(hetcols, '_Dif'))

cn <- setNames(rep('', length(col_names)), col_names)
df <- bind_rows(cn)[0,]

## Traverse max trans. df and get coverage cols

## which(tdf$Gene_id == 'PF3D7_0100300')
## i <- 63
## tdf[63,]

for (i in 1:dim(tdf)[1]){

  gid <- as.character(tdf$Gene_id[i])
  on <- tdf$On_trans[i]
  off <- tdf$Off_trans[i]

  ## Select contrast in which difference sign must be switched
  ## contrast <- paste0(on, '-', off)
  ## positive_contrasts <- c(
  ##   '12B-10G', '12B-A7', '12B-E5', '12B-B11',
  ##   '10G-A7', '10G-E5', '10G-B11',
  ##   'A7-E5', 'A7-B11',
  ##   'B11-E5'
  ## )

  ## ifelse(
  ##   contrast %in% positive_contrasts,
  ##   neg_dif <- FALSE,
  ##   neg_dif <- TRUE
  ##   )

  ## Main Loop
  if (gid %in% positive_cov_df$Gene_id & !is.na(on) & !is.na(off)){

    onvect <- positive_cov_df %>%
      filter(Gene_id == gid) %>%
      select(contains(on))

    offvect <- positive_cov_df %>%
      filter(Gene_id == gid) %>%
      select(contains(off))

    diffvect <- onvect-offvect
    #if (neg_dif) {diffvect <- -diffvect}

    row <- c(gid,
             unlist(onvect, use.names = F),
             unlist(offvect, use.names = F),
             unlist(diffvect, use.names = F))
    row <- setNames(row, col_names)
    df <- df %>% add_row(bind_rows(row))
  }
}

unsigned_maxtrans_cov <- df %>%
  mutate(across(-Gene_id,  as.numeric)) %>%
  full_join(tdf, by = 'Gene_id') %>%
  mutate(Max_aMAFC = abs(Max_aMAFC)) %>%
  left_join(info_df)

genes_nocov <- unsigned_maxtrans_cov %>%
     select(Gene_id, contains('Cov_500fp_manybins')) %>%
     filter(!complete.cases(.)) %>%
     select(Gene_id) %>% pull()

maxTrans_df %>%
  filter(is.na(On_trans) | is.na(Off_trans))

unsigned_maxtrans_cov %>%
  select(Gene_id, contains('Cov_500fp_manybins')) %>%
  filter(!complete.cases(.)) %>%
  print(width = 400)

## Non-Variant With transFC analysis
non_variant_transDif <- unsigned_maxtrans_cov %>%
  filter(!Variant) %>%
  filter(!Is_tRNA) %>%
  filter(abs(Max_aMAFC) > 1) %>%
  select(Label, Max_aMAFC, On_trans, Off_trans, Family, Gam_specific, Annot, Gene_id)

write_csv(non_variant_transDif, 'non_variant_transth_1.csv')

non_variant_transDif %>%
  print(n = 50)

non_variant_transDif %>%
  count(Family)

non_variant_transDif %>%
  count(Gam_specific)

non_variant_transDif %>%
  filter(!is.na(Family))

non_variant_transDif %>%
  filter(On_trans == 'B11' | Off_trans == 'B11')


non_variant_transDif %>%
  filter(On_trans != '12B' & Off_trans != '12B' &
         On_trans != '10G' & Off_trans != '10G') %>%
  select(Gene_id, Max_aMAFC) %>%
  write_csv('non_variant_log2FC1_newarrays.csv')


non_variant_transDif

unsigned_maxtrans_cov %>%
  filter(abs(Max_aMAFC) > 1) %>%
  filter(Variant) %>%
  filter(On_trans != '12B' & Off_trans != '12B' &
         On_trans != '10G' & Off_trans != '10G') %>%
  write_csv('variant_log2fc1_newarrays.csv')

unsigned_maxtrans_cov %>%
  filter(abs(Max_aMAFC) > 1) %>%
  count(Family_Grouped)

unsigned_maxtrans_cov %>%
  filter(abs(Max_aMAFC) > 1.5)

hist(abs(unsigned_maxtrans_cov$Max_aMAFC), breaks = 1000)


unsigned_maxtrans_cov %>%
  filter(abs(Max_aMAFC) > 3) %>%
  count(Variant) %>%
  filter(Variant) %>%
  pull()

unsigned_maxtrans_cov %>%
  filter(abs(Max_aMAFC) > 2)


## percs <- c()
## for (x in seq(0, 7, by = 0.1)){

##   total <- dim(unsigned_maxtrans_cov %>%
##     filter(abs(Max_aMAFC) > x))[1]

##   var <- unsigned_maxtrans_cov %>%
##     filter(abs(Max_aMAFC) > x) %>%
##     count(Variant) %>%
##     filter(Variant) %>%
##     pull()

##   perc <- (var/total)*100
##   percs <- c(percs, perc)
## }
## percs

## df <- tibble('Thresolds' = seq(0, 7, by = 0.1), 'Perc.Variant' = percs)


## p <- ggplot(df, aes(x = Thresholds, y = Perc.Variant))
## p <- p + geom_point()
## p <- p + scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
## p <- p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
## p

## Get Correlations

cor_signed_maxtrans <- signed_maxtrans_cov %>%
  filter(!Is_3D7B)

## Subset column by column
clnms <- colnames(cor_signed_maxtrans %>% select(contains('_Dif')))
for (c in clnms) {
  print(c)
  cor_genes <- cor_signed_maxtrans %>%
    select(Max_aMAFC, matches(c)) %>%
    filter(abs(Max_aMAFC) > 1.5) %>%
    drop_na()

  cormtx <- cor(cor_genes, method = 'pearson')
  print(paste0('Number of genes: ', dim(cor_genes)[1]))
  print(as.data.frame(cormtx[,1]))
  print('------------------')
}

## All toghether

trans_th <- 1

cor_genes <- cor_signed_maxtrans %>%
  select(Max_aMAFC, contains('_Dif')) %>%
  filter(abs(Max_aMAFC) > trans_th) %>%
  drop_na()

cor_genes

cormtx <- cor(cor_genes, method = 'pearson')
print(paste0('Number of genes: ', dim(cor_genes)[1]))
print(as.data.frame(cormtx[,1]))

outdf <- as.data.frame(cormtx)
outdf['Corr_with'] <- rownames(outdf)
outdf <- outdf %>% select(Corr_with, everything())
outname = paste0('corr_MaxFC_', trans_th, '_allowoverlaps.csv')
write_csv(outdf, outname)

## By subsets

get_corr <- function(df){
  cor_genes <- df %>%
    select(-Gene_id) %>%
    filter(abs(Max_aMAFC) > trans_th) %>%
    drop_na()

  cormtx <- cor(cor_genes, method = 'pearson')
  print(paste0('Number of genes: ', dim(cor_genes)[1]))

  outdf <- as.data.frame(cormtx)
  outdf['Corr_with'] <- rownames(outdf)
  outdf <- outdf %>% select(Corr_with, everything())
  outdf <- as_tibble(outdf[,c(1,2)])
  print(outdf)
}

trans_th <- 1

## Abs bins
abs_bins <- cor_signed_maxtrans %>%
  select(Gene_id, Max_aMAFC,
         Cov_1500fp_manybins_Dif,
         Cov_1000fp_manybins_Dif,
         Cov_500fp_manybins_Dif,
         contains('ORF', ignore.case = F) & contains('Dif'),
         Cov_500tp_manybins_Dif,
         Cov_1000tp_manybins_Dif
         )

abs_bins_cor <- get_corr(abs_bins)

## Rel bins

rel_bins <- cor_signed_maxtrans %>%
  select(Gene_id, Max_aMAFC,
         p5utr_Cov_Dif,
         contains('qorf') & contains('Dif'),
         p3utr_Cov_Dif
         )

rel_bins_cor <- get_corr(rel_bins)

## TSS

tss_bins <- cor_signed_maxtrans %>%
  select(Gene_id, Max_aMAFC,
         TSS_Cov_Dif,
         Cov_500fp_manybins_Dif,
         p5utr_Cov_Dif
         )

tss_bins_cor <- get_corr(tss_bins)


## Plots

df <- cor_signed_maxtrans %>% filter(abs(Max_aMAFC) > 1)

for (c in clnms){
  p <- ggplot(df,
              aes(x = Max_aMAFC, y = get(c), color = Gam_specific)) +
    geom_point() +
    geom_label_repel(data=subset(df, get(c) > 1),
                     aes(label = Gene_id),
                     box.padding   = 0.35,
                     point.padding = 0.5,
                     segment.color = 'grey50')
  ggsave(paste0('./Plots/MaxFC_Cor_Plots/cor_signed_maxFC_', c, '_allowoverlaps.png'), p, device = 'png')
}

corr_data <- as_tibble(outdf[,1:2])

corr_data %>%
  arrange(Max_aMAFC) %>%
  print(n = 40)

bins <- c(
  'Cov_2000fp_manybins_Dif',
  'Cov_1500fp_manybins_Dif',
  'Cov_1000fp_manybins_Dif',
  'Cov_500fp_manybins_Dif',
  'Cov_500orf_Dif',
  'Cov_1000orf_Dif',
  'Cov_1500orf_Dif',
  'Cov_500tp_manybins_Dif',
  'Cov_1000tp_manybins_Dif',
  'Cov_1500tp_manybins_Dif'
  )


## Geom_step approach

pos <- c(-2000, -1500, -1000, -500, 0, 500, 1000, 2500, 3000, 3500)

cor_bins <- corr_data %>%
  filter(Corr_with %in% bins) %>%
  arrange(factor(Corr_with, levels = bins)) %>%
  mutate(Pos = pos)

cor_bins

cor_bins2 <- cor_bins %>%
  add_row(Corr_with = 'start', Max_aMAFC = 0, Pos = -2500) %>%
  add_row(Corr_with = 'end_of_gene', Max_aMAFC = 0, Pos = 1500) %>%
  add_row(Corr_with = 'end', Max_aMAFC = 0, Pos = 4000) %>%
  add_row(Corr_with = 'endtail', Max_aMAFC = 0, Pos = 4500) %>%
  arrange(Pos)

xbreaks <- c(
  '',
  'ATG-2000',
  'ATG-1500',
  'ATG-1000',
  'ATG-500',
  'ATG',
  'ATG+500',
  'ATG+1000',
  'ATG+1500',
  'End',
  'End+500',
  'End+1000',
  'End+1500',
  ''
  )

p <- ggplot(cor_bins2, aes(x = Pos, y = -Max_aMAFC, group = 1)) +
  geom_step() +
  ylim(0, 1) +
  ylab('- Pearson Corr') + xlab('Region') +
  theme_minimal() +
  scale_x_continuous(breaks = cor_bins2$Pos, labels = xbreaks) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'gray') +
  geom_vline(xintercept = 2500, linetype = 'dashed', color = 'gray')

p
ggsave('corr_genemodel_steps_allowoverlaps.png', p, device = 'png')

## Geom Bar approach

cor_bins_plot <- function(df, labs, regions){
  corplot_df <- df %>%
    filter(Corr_with != 'Max_aMAFC') %>%
    mutate(Corr_with = factor(Corr_with, levels = Corr_with)) %>%
    arrange(Corr_with) %>%
    mutate(Labs = factor(labs, levels = labs)) %>%
    mutate(Region = regions)

  p <- ggplot(corplot_df, aes(x = Labs, y = -Max_aMAFC, color = Region)) +
    geom_boxplot(key_glyph = "point") +
    ylim(0, 1) +
    ylab('- Pearson Corr\n') + xlab('\nRegion') +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      text = element_text(size=30),
            legend.position = "none"
          )

  ggsave(paste0('./Plots/Correlations/', deparse(substitute(df)), '_allowoverlaps.pdf'), p, device = 'pdf')
  p
}


## Abs

labs = c(
  '-1500 to -1000',
  '-1000 to -500',
  '-500 to ATG',
  'ATG to 500',
  '500 to 1000',
  '1000 to 1500',
  'End to +500',
  '+500 to +1000'
)

regions = c(rep('5prime', 3), rep('ORF', 3), rep('3prime', 2))

cor_bins_plot(abs_bins_cor, labs, regions)

## Rel

labs = c(
  '5\'UTR',
  'ORF 1/4',
  'ORF 2/4',
  'ORF 3/4',
  'ORF 4/4',
  '3\'UTR'
)

regions = c(rep('5prime', 1), rep('ORF', 4), rep('3prime', 1))

cor_bins_plot(rel_bins_cor, labs, regions)

## TSS

labs = c(
  'TSS',
  '-500 to ATG',
  '5\'UTR'
)

regions = c(rep('5prime', 3), rep('ORF', 0), rep('3prime', 0))

cor_bins_plot(tss_bins_cor, labs, regions) + geom_boxplot(color = 'green')



## cor_bins <- corr_data %>%
##   filter(Corr_with %in% bins) %>%
##   mutate(Corr_with = factor(Corr_with, levels = bins)) %>%
##   arrange(Corr_with)

## labs = c(
##   'ATG-2000',
##   'ATG-1500',
##   'ATG-1000',
##   'ATG-500',
##   'ATG+500',
##   'ATG+1000',
##   'ATG+1500',
##   'End+500',
##   'End+1000',
##   'End+1500'
## )

## cor_bins2 <- cor_bins %>%
##   mutate(Region = c(rep('5prime', 4), rep('ORF', 3), rep('3prime', 3))) %>%
##   mutate(Labs = factor(labs, levels = labs))

## p <- ggplot(cor_bins2, aes(x = Labs, y = -Max_aMAFC, color = Region)) +
##   geom_boxplot(key_glyph = "point") +
##   ylim(0, 1) +
##   ylab('- Pearson Corr') + xlab('Region') +
##   theme_minimal() +
##   theme(axis.text.x = element_text(angle = 45, hjust = 1))

## p

## ggsave('corr_genemodel_lines_abs.png', p, device = 'png')


## bins <- c(
##   'TSS_Cov_Dif',
##   'p5utr_Cov_Dif',
##   'Cov_1qorf_manybins_Dif',
##   'Cov_2qorf_manybins_Dif',
##   'Cov_3qorf_manybins_Dif',
##   'Cov_4qorf_manybins_Dif',
##   'p3utr_Cov_Dif'
## )

## cor_bins2 <- corr_data %>%
##   filter(Corr_with %in% bins) %>%
##   arrange(factor(Corr_with, levels = bins)) %>%
##   mutate(Corr_with = factor(Corr_with, levels = bins)) %>%
##   mutate(Region = c(rep('5prime', 2), rep('ORF', 4), rep('3prime', 1)))

## p <- ggplot(cor_bins2, aes(x = Corr_with, y = -Max_aMAFC, color = Region)) +
##   geom_boxplot(key_glyph = "point") +
## #  geom_bar(stat='identity', aes(fill = Region), width = 1, color = 'darkslategray') +
##   ylim(0, 1) +
##   ylab('- Pearson Corr') + xlab('Region') +
##   theme_minimal() +
##   theme(axis.text.x = element_text(angle = 45, hjust = 1))
## p

## ggsave('corr_genemodel_lines_rel.png', p, device = 'png')

## corr_data$Corr_with

trans_th <- 1

cor_genes <- cor_signed_maxtrans %>%
  select(Max_aMAFC, contains('_Dif')) %>%
  filter(abs(Max_aMAFC) > trans_th) %>%
  drop_na()

cor_genes

cormtx <- cor(cor_genes, method = 'pearson')
print(paste0('Number of genes: ', dim(cor_genes)[1]))
print(as.data.frame(cormtx[,1]))

outdf <- as.data.frame(cormtx)
outdf['Corr_with'] <- rownames(outdf)
outdf <- outdf %>% select(Corr_with, everything()) %>% as_tibble()

plotdf <- outdf %>% select(Corr_with, Max_aMAFC) %>% arrange(Max_aMAFC)
mplotdf <- melt(plotdf)
mplotdf$Corr_with <- factor(mplotdf$Corr_with, levels = mplotdf$Corr_with)

ggplot(mplotdf, aes(x = variable, y = Corr_with)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label=value)) +
  scale_fill_gradient2(low = "chartreuse3",
                       mid = "white",
                       high = "darkred",
                       na.value="grey90")

plot_df <- unsigned_maxtrans_cov %>%
  arrange(desc(Max_aMAFC)) %>%
  mutate(Gene_id = factor(Gene_id, levels = Gene_id))

heatMap_allstrains_trans <- function(aMAFC_th, cov_col, family, family_facet){

  on_col <- paste0(cov_col, '_On')
  off_col <- paste0(cov_col, '_Off')
  dif_col <- paste0(cov_col, '_Dif')

  subset_all <- plot_df %>%
    filter(abs(Max_aMAFC) > aMAFC_th) %>%
    select(Gene_id,
           Label,
           Max_aMAFC,
           matches(on_col),
           matches(off_col),
           matches(dif_col),
           Annot,
           Family,
           SubFamily
           ) %>%
    filter(complete.cases(.))

  if (!is.na(family)){
    subset_all <- subset_all %>%
      filter(Family == family)
  }

  ## Ordering
  if (sorting == 'clust') {
    mtx <- subset_all %>%
      select(matches(on_col),
             matches(off_col),
             matches(dif_col))
      #select(where(is.numeric))

    ##Make hierarquical Clustering
    dmtx <- dist(scale(mtx), method = "euclidean")
    cl <- hclust(dmtx, method = 'average')
    subset_all['Cluster'] <- cutree(cl, nclust)
    subset_all$Label <- factor(subset_all$Label, levels = subset_all$Label[cl$order])

    ## ## Make k-means Clustering
    ## cl <- kmeans(scale(mtx), nclust)
    ## subset_all$Label <- factor(subset_all$Label, levels = subset_all$Label[cl$cluster])

  } else if (sorting == 'trans') {
    subset_all <- subset_all %>%
      arrange(desc(Max_aMAFC)) %>%
      mutate(Label = factor(Label, levels = Label))

  } else if (sorting == 'het') {
    subset_all <- subset_all %>%
      arrange(desc(dif_col)) %>%
      mutate(Label = factor(Label, levels = Label))
  }


  subset_trans <- subset_all %>%
    select(Gene_id, Max_aMAFC, Label, Annot, Family, SubFamily)

  subset_het <- subset_all %>%
    select(Gene_id, matches(on_col), matches(off_col), Label, Annot, Family, SubFamily)

  subset_hetDif <- subset_all %>%
    select(Gene_id, matches(dif_col), Label, Annot, Family, SubFamily)

  melt_vars <- c('Gene_id', 'Label', 'Annot', 'Family', 'SubFamily')
  mdf_het <- melt(subset_het, id.vars = melt_vars)
  mdf_trans <- melt(subset_trans, id.vars = melt_vars)
  mdf_hetDif <- melt(subset_hetDif, id.vars = melt_vars)

  het_plot <- hetHeatmap(mdf_het, family_facet)
  trans_plot <- transHeatmap(mdf_trans, family_facet)
  hetDif_plot <- hetDifHeatmap(mdf_hetDif, family_facet)
  all_plot <- grid.arrange(trans_plot, hetDif_plot, het_plot, nrow = 1, widths = c(1,1,3))
}

aMAFC_th <- 1.5
family_facet <- F
family <- NA
cov_col <- 'Cov_500fp_manybins'
sorting <- 'trans'
nclust <- 5

##colnames(plot_df)


## Inspect NAs
plot_df %>%
  filter(is.na(Cov_500orf_Dif)) %>%
  select(Gene_id, Name, contains('Cov_500orf')) %>%
  print(n = 300)

heatMap_allstrains_trans(aMAFC_th, cov_col, family, family_facet)

plot_df %>%
  filter(Name == 'CLAG3.1') %>%
  select(Gene_id)

plot_df <- unsigned_maxtrans_cov %>%
  arrange(desc(Max_aMAFC)) %>%
  mutate(Gene_id = factor(Gene_id, levels = Gene_id))

length(unsigned_maxtrans_cov$Label)
unique(length(unsigned_maxtrans_cov$Label))

plot_df %>%
  filter(abs(Max_aMAFC) > 1)

plot_df %>%
  filter(abs(Max_aMAFC) > 1) %>%
  select(Family_Grouped) %>%
  table()

plot_df %>%
  filter(abs(Max_aMAFC) > 1) %>%
  filter(Family_Grouped == 'Not CVGs') %>%
  select(Gene_id, Name, Max_aMAFC, Gam_specific) %>%
  arrange(-abs(Max_aMAFC)) %>%
  print(n=40)

plot_df$Chrom


plot_df %>%
  filter(abs(Max_aMAFC) > 1) %>%
  filter(Family_Grouped == 'Not CVGs') %>%
  filter(Annot == 'conserved Plasmodium protein, unknown function')

plot_df %>%
  filter(abs(Max_aMAFC) > 1) %>%
  filter(Family_Grouped == 'Not CVGs') %>%
  filter(Annot == 'Plasmodium exported protein, unknown function')

plot_df %>%
  filter(abs(Max_aMAFC) > 1) %>%
  filter(Family_Grouped == 'Not CVGs') %>%
  select(Gam_specific) %>%
  table()


table(info_df$Difpeaks_E5_B11)

heatMap_allstrains_trans <- function(aMAFC_th, cov_col, family, family_facet){

  on_col <- paste0(cov_col, '_On')
  off_col <- paste0(cov_col, '_Off')
  dif_col <- paste0(cov_col, '_Dif')

  subset_all <- plot_df %>%
    filter(abs(Max_aMAFC) > aMAFC_th) %>%
    select(Gene_id,
           Label,
           Max_aMAFC,
           matches(on_col),
           matches(off_col),
           matches(dif_col),
           Annot,
           Family,
           Family_Grouped,
           SubFamily
           ) %>%
    filter(complete.cases(.))

  if (!is.na(family)){
    subset_all <- subset_all %>%
      filter(Family == family)
  }


  ## Hierarchical Clustering
  mtx <- subset_all %>%
    select(matches(on_col),
           matches(off_col))

  ##Make hierarquical Clustering
  dmtx <- dist(scale(mtx), method = "euclidean")
  cl <- hclust(dmtx, method = 'average')
  subset_all['Cluster'] <- cutree(cl, nclust)

  ## Ordering
  if (sorting == 'clust') {
    subset_all$Label <- factor(subset_all$Label, levels = subset_all$Label[cl$order])

    ## ## Make k-means Clustering
    ## cl <- kmeans(scale(mtx), nclust)
    ## subset_all$Label <- factor(subset_all$Label, levels = subset_all$Label[cl$cluster])

  } else if (sorting == 'trans') {
    subset_all <- subset_all %>%
      arrange(desc(Max_aMAFC)) %>%
      mutate(Label = factor(Label, levels = Label))

  } else if (sorting == 'het') {
    subset_all <- subset_all %>%
      arrange(desc(dif_col)) %>%
      mutate(Label = factor(Label, levels = Label))
  }


  subset_trans <- subset_all %>%
    select(Gene_id, Max_aMAFC, Label, Annot, Family, Family_Grouped, SubFamily, Cluster)

  subset_het <- subset_all %>%
    select(Gene_id, matches(on_col), matches(off_col), Label, Annot, Family, Family_Grouped, SubFamily, Cluster)

  subset_hetDif <- subset_all %>%
    select(Gene_id, matches(dif_col), Label, Annot, Family, Family_Grouped, SubFamily, Cluster)

  melt_vars <- c('Gene_id', 'Label', 'Annot', 'Family', 'Family_Grouped', 'SubFamily', 'Cluster')
  mdf_het <- melt(subset_het, id.vars = melt_vars)
  mdf_trans <- melt(subset_trans, id.vars = melt_vars)
  mdf_hetDif <- melt(subset_hetDif, id.vars = melt_vars)

  het_plot <- hetHeatmap(mdf_het, family_facet)

  trans_plot <- transHeatmap(mdf_trans, family_facet)

  hetDif_plot <- hetDifHeatmap(mdf_hetDif, family_facet)
  hetDif_plot <- hetDif_plot + scale_fill_gradient(low = "chartreuse3",
                                                   high = "white",
                                                   na.value="white")

  #all_plot <- grid.arrange(trans_plot, hetDif_plot, het_plot, nrow = 1, widths = c(1,1,3))
  all_plot <- grid.arrange(hetDif_plot, het_plot, nrow = 1, widths = c(1,3))
  #ggsave('./Plots/MaxAMFC_Heatmaps/max_mafc_cov.svg', all_plot, device = 'svg')

}

aMAFC_th <- 3
family_facet <- T
family <- NA
cov_col <- 'Cov_500fp_manybins'
sorting <- 'clust'
nclust <- 5

##colnames(plot_df)
heatMap_allstrains_trans(aMAFC_th, cov_col, family, family_facet)

## ## Inspect NAs
## plot_df %>%
##   filter(is.na(Cov_500orf_Dif)) %>%
##   select(Gene_id, Name, contains('Cov_500orf')) %>%
##   print(n = 300)

#### Load Data ####
cov_dir <- '/home/lucas/ISGlobal/Projects/Phd_Project/ChIP_Seq/RPKMs_normInput/'

## Coverage
gm_12B <- read_csv(paste0(cov_dir, '1.2B_me_sort_q5_RPKMs_normInput_5binned_cov_2prevpost.csv'))
gm_10G <- read_csv(paste0(cov_dir, '10G_me_sort_q5_RPKMs_normInput_5binned_cov_2prevpost.csv'))
gm_A7 <- read_csv(paste0(cov_dir, 'A7K9_me_sort_q5_RPKMs_normInput_5binned_cov_2prevpost.csv'))
gm_E5 <- read_csv(paste0(cov_dir, 'E5K9_me_sort_q5_RPKMs_normInput_5binned_cov_2prevpost.csv'))
gm_B11 <- read_csv(paste0(cov_dir, 'B11_me_sort_q5_RPKMs_normInput_5binned_cov_2prevpost.csv'))


gm_strains <- list(df12B=gm_12B,
                df10G=gm_10G,
                dfA7=gm_A7,
                dfE5=gm_E5,
                dfB11=gm_B11
                )

## Change dfs in a list
gm_strains <- lapply(gm_strains, function(df) {
    colnames(df)[1] <- "Gene_id"
    df
})

## Convert list into individual objects again
list2env(gm_strains, envir=.GlobalEnv)

##Create numeric non-na mtxs
nona.mtxs <- lapply(gm_strains, function(df) {
    mtx <- as.matrix(df[complete.cases(df),-1])
    rownames(mtx) <- df[complete.cases(df),1] %>% pull()
    mtx
})



nona_gm_strains <- lapply(gm_strains, function(tibble) {
  tibble %>%
    filter(complete.cases(tibble))
  tibble
})


names(nona.mtxs) <- c("mtx12B", "mtx10G", 'mtxA7', 'mtxE5', 'mtxB11')
names(nona_gm_strains) <- c('nona_12B', 'nona_10G', 'nona_A7', 'nona_E5', 'nona_B11')

## Create side-by-side matrices
#mtx_12b10g <- cbind(nona.mtxs$mtx12B, nona.mtxs$mtx10G)

mtx_all <- cbind(nona.mtxs$mtx12B,
                 nona.mtxs$mtx10G,
                 nona.mtxs$mtxA7,
                 nona.mtxs$mtxE5,
                 nona.mtxs$mtxB11
                 )

colnames(mtx_all) <- 1:dim(mtx_all)[2]

bins <- paste0('bin', c(1:23))
all_cols <- c(paste0(bins, '_12B'),
              paste0(bins, '_10G'),
              paste0(bins, '_A7'),
              paste0(bins, '_E5'),
              paste0(bins, '_B11')
              )

gmodel_all <- gm_strains$df12B %>%
  full_join(gm_strains$df10G, by = 'Gene_id') %>%
  full_join(gm_strains$dfA7, by = 'Gene_id') %>%
  full_join(gm_strains$dfE5, by = 'Gene_id') %>%
  full_join(gm_strains$dfB11, by = 'Gene_id') %>%
  setNames(c('Gene_id', all_cols))

nona_gm_all <- gmodel_all %>%
  filter(complete.cases(gmodel_all))

## Clonally Variant genes

#cvgDF <- geneDFnona[geneDFnona$Epi != "non-variant",]
cvgDF <- nona_gm_all %>%
  left_join(info_df)



#cvg_mtx <- dif_12b10g[rownames(dif_12b10g) %in% cvgDF$Gene_id,]
cvg_mtx <- cvgDF %>%
  filter(Variant) %>%
  select(contains('bin'))

cvg_pca <- prcomp(cvg_mtx)
cvg_pca_df <- as.data.frame(cvg_pca$x[, c(1,2)])
cvg_pca_df <- cbind(cvg_pca_df, cvgDF %>% filter(Variant))
#alpha <- sapply(dif_pca_df$Epi == "non-variant", function(x) if (x) {0.1} else {1})


## All genes, filtered by transcription

trans_th <- 0

gm_trans <- maxTrans_df %>%
  left_join(gmodel_all) %>%
  left_join(info_df) %>%
  filter(abs(Max_aMAFC) > trans_th) %>%
  filter(across(contains('bin'), complete.cases))

gm_trans_mtx <- gm_trans %>%
  select(contains('bin'))

gm_trans_pca <- prcomp(gm_trans_mtx)
gm_trans_df <- as_tibble(gm_trans_pca$x[, c(1,2)])
gm_trans <- bind_cols(gm_trans_df,  gm_trans)

## Create 'fixed' color palette
col_factor <- factor(unique(analysis_df$Family_Grouped),
                     levels=c('Not CVGs',
                              'Other CVGs',
                              'FIKK',
                              'HYP',
                              'PHIST',
                              'STEVOR',
                              'RIFIN',
                              'VAR'))

my_colors <- c('gray', scales::viridis_pal()(7))
names(my_colors) <- levels(col_factor)
my_scale2 <- scale_color_manual(name = "Family_Grouped", values = my_colors)

#### PCA plots ####
## All-genes, transcription filter
p <- ggplot(gm_trans, aes(x=PC1, y=PC2, color = Family_Grouped))
p <- p + geom_point() + my_scale2
p

ggsave('./Plots/pca_families.svg', p, device = 'svg')

p <- ggplot(gm_trans, aes(x=PC1, y=PC2, color = Gam_specific))
p <- p + geom_point()
p

p <- ggplot(gm_trans, aes(x=PC1, y=PC2, color = Max_aMAFC))
p <- p + geom_point()
p <- p + scale_color_gradient2(midpoint = 0, low="red", mid = "black", high="green")
p

## CVGs
p <- ggplot(cvg_pca_df, aes(x=PC1, y=PC2, color = Family_Grouped))
p <- p + geom_point() + my_scale2
p

p <- ggplot(cvg_pca_df, aes(x=PC1, y=PC2, color = Gam_specific))
p <- p + geom_point()
p

p <- ggplot(cvg_pca_df, aes(x=PC1, y=PC2, color = Variant))
p <- p + geom_point()
p <- p
p

customHeatmap <- function(df, family_facet, limits){
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
  if (family_facet){p <- p+facet_grid(vars(Family), scales = "free", space = "free")}
  p
}

## Absolute value of difference
#heatDF <- cbind(hetdifDF[,-3], abs(hetdif_mtx))

gm_dif_heatmap <- function(contrast, trans_th, difpeaks, family_facet){
  difpeaks_col <- paste0('Difpeaks_', contrast[1], '_', contrast[2])
  trans_col <- paste0(contrast[1], '-', contrast[2], '_MaxVal')

  gm_dif <- nona_gm_all %>%
    left_join(info_df) %>%
    left_join(trans_df) %>%
    filter(abs(get(trans_col)) > trans_th)

  if (difpeaks) {gm_dif <- gm_dif %>% filter(get(difpeaks_col))}

  x <- gm_dif %>%
    select(contains(contrast[1]) & contains('bin'))

  y <- gm_dif %>%
    select(contains(contrast[2]) & contains('bin'))

  mtx <- abs(x-y)

  out_df <- gm_dif %>%
    select(-contains('bin'), -contains('MaxVal'), -contains('MaxTime')) %>%
    bind_cols(mtx)

  vals <- out_df %>%
    select(contains('bin')) %>%
    pull()

  ## Make hierarquical Clustering
  dmtx <- dist(mtx, method = "euclidean")
  cl <- hclust(dmtx, method = 'average')
  out_df$Gene_id <- factor(out_df$Gene_id, levels = out_df$Gene_id[cl$order])
  customHeatmap(out_df, family_facet, c(min(vals),max(vals)))
}

contrast <- c('12B', '10G')
trans_th <- 1
difpeaks <- F
family_facet <- F

gm_dif_heatmap(contrast, trans_th, difpeaks, family_facet)

#difpeak_mtx <- mtx_12b10g[rownames(mtx_12b10g) %in% hetdifgenes,]
#heatDF <- cbind(hetdifDF[,-4], difpeak_mtx)

gm_sidebyside_heatmap <- function(contrast, family_facet, trans_th){
  trans_col <- paste0(contrast[1], '-', contrast[2], '_MaxVal')

  gm_df <- nona_gm_all %>%
    select(Gene_id, contains(contrast)) %>%
    left_join(info_df, by = 'Gene_id') %>%
    left_join(trans_df, by = 'Gene_id') %>%
    filter(abs(get(trans_col)) > trans_th)

  gm_mtx <- gm_df %>%
    select(contains('bin'))


  ## Make hierarquical Clustering
  dmtx <- dist(gm_mtx, method = "euclidean")
  cl <- hclust(dmtx, method = 'complete')
  gm_df$Gene_id <- factor(gm_df$Gene_id, levels = gm_df$Gene_id[cl$order])

  #customHeatmap(gm_df, c(min(gm_mtx), max(gm_mtx)))

  #head(gm_df)
  mdf <- melt(gm_df %>% select(-c(contains('MaxVal'), contains('MaxTime'))))
  mdf["Strain"] <- sapply(mdf$variable, function(x) strsplit(as.character(x), split = "_", fixed = TRUE)[[1]][2])

  p <- ggplot(mdf, aes(x = variable, y = Gene_id, fill = value)) +

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

  if (family_facet){p <- p+facet_grid(vars(Family), vars(Strain), scales = "free", space = "free")}

  p

}

contrast <- c('E5', 'B11')
family_facet <- F
trans_th <- 1.5

gm_sidebyside_heatmap(contrast, family_facet, trans_th)


##ggsave("/mnt/Disc4T/Projects/PhD_Project/dif12b_10g_heatmap.png", p,
##       device = "png", width = 40, height = 20,  units = "cm")


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

## tdf <- maxTrans_df %>%
##   select(Gene_id, Max_aMAFC, On_trans, Off_trans, Is_3D7B)

tdf <- transFC_df_noNAs %>%
  select(Gene_id, Max_aMAFC, On_trans, Off_trans, Is_3D7B)

gmodel_all_pos <- gmodel_all %>%
  mutate(across(-Gene_id, ~ ifelse(.x > 0, .x, 0)))

## Create empty DF
hetcols <- gmodel_all %>%
  select(contains('12B') & contains('bin')) %>%
  colnames(.) %>%
  gsub('12B', '', ., fixed=TRUE)

col_names <- c('Gene_id',
               paste0(hetcols, 'On'),
               paste0(hetcols, 'Off'),
               paste0(hetcols, 'Dif'))

cn <- setNames(rep('', length(col_names)), col_names)
df <- bind_rows(cn)[0,]

## Traverse max trans. df and get coverage cols

for (i in 1:dim(tdf)[1]){

  gid <- as.character(tdf$Gene_id[i])
  on <- tdf$On_trans[i]
  off <- tdf$Off_trans[i]

  ## Select contrast in which difference sign must be switched
  ## contrast <- paste0(on, '-', off)
  ## positive_contrasts <- c(
  ##   '12B-10G', '12B-A7', '12B-E5', '12B-B11',
  ##   '10G-A7', '10G-E5', '10G-B11',
  ##   'A7-E5', 'A7-B11',
  ##   'B11-E5'
  ## )

  ## ifelse(
  ##   contrast %in% positive_contrasts,
  ##   neg_dif <- FALSE,
  ##   neg_dif <- TRUE
  ## )

  ### Main Loop
  if (gid %in% gmodel_all_pos$Gene_id & !is.na(on) & !is.na(off)){

    onvect <- gmodel_all_pos %>%
      filter(Gene_id == gid) %>%
      select(contains(on))

    offvect <- gmodel_all_pos %>%
      filter(Gene_id == gid) %>%
      select(contains(off))

    diffvect <- onvect-offvect
    #if (neg_dif) {diffvect <- -diffvect}

    row <- c(gid,
             unlist(onvect, use.names = F),
             unlist(offvect, use.names = F),
             unlist(onvect-offvect, use.names = F))
    row <- setNames(row, col_names)
    df <- df %>% add_row(bind_rows(row))
  }
}

df <- df %>%
  mutate(across(-Gene_id, as.numeric))

df %>% select(contains('Off'))
tdf %>% filter(Gene_id == 'PF3D7_0832200')
gmodel_all %>%
  filter(Gene_id == 'PF3D7_0832200') %>%
  select(contains('B11'))

gm_pos_maxtrans <- df %>%
  full_join(tdf, by = 'Gene_id') %>%
  left_join(info_df) %>%
  mutate(Max_aMAFC = abs(Max_aMAFC))

heatMap_allstrains_trans <- function(df, aMAFC_th, family, family_facet){

  ## Returns a list object with list = (df, plot)

  subset_all <- df %>%
    filter(abs(Max_aMAFC) > aMAFC_th)

  if (!is.na(family)){
    subset_all <- subset_all %>%
      filter(Family == family)
  }

  subset_all <- subset_all %>%
    filter(across(
      contains('bin'),
      complete.cases
    ))

  ## mtx <- subset_all %>%
  ##   select(contains('bin'))

  mtx <- subset_all %>%
    select(contains('bin'),
           -contains('bin1_'),
           -contains('bin2_'),
           -contains('bin3_'),
           -contains('bin4_'),
           -contains('bin20_'),
           -contains('bin21_'),
           -contains('bin22_'),
           -contains('bin23_'),
           )

  ## Make hierarquical Clustering
  dmtx <- dist(mtx, method = "euclidean")
  cl <- hclust(dmtx, method = 'complete')
  subset_all$Label <- factor(subset_all$Label, levels = subset_all$Label[cl$order])
  subset_all['Cluster'] <- cutree(cl, nclust)
  subset_all$Label <- factor(subset_all$Label, levels = subset_all$Label[cl$order])

  subset_trans <- subset_all %>%
    select(Gene_id, Max_aMAFC, Label, Name, Annot, Family, SubFamily, Cluster)

  subset_het_on <- subset_all %>%
    select(Gene_id, Label,
           contains('On') & contains('bin'),
           Name, Annot, Family, SubFamily, Cluster)

  subset_het_off <- subset_all %>%
    select(Gene_id, Label,
           contains('Off') & contains('bin'),
           Name, Annot, Family, SubFamily, Cluster)

  ##write_csv(subset_het %>% select(-contains('bin')) %>% arrange(Cluster), './het_trans_genemodel.csv')

  subset_hetDif <- subset_all %>%
    select(Gene_id, contains('_Dif'), Label, Name, Annot, Family, SubFamily, Cluster)

  melt_vars <- c('Gene_id', 'Label', 'Name', 'Annot', 'Family', 'SubFamily', 'Cluster')
  mdf_het_on <- melt(subset_het_on, id.vars = melt_vars)
  mdf_het_off <- melt(subset_het_off, id.vars = melt_vars)
  mdf_trans <- melt(subset_trans, id.vars = melt_vars)
  mdf_hetDif <- melt(subset_hetDif, id.vars = melt_vars)

  het_plot_on <- hetHeatmap(mdf_het_on, family_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),

      panel.border=element_blank(),
      panel.grid.major=element_blank(),

      strip.background = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      )


  het_plot_off <- hetHeatmap(mdf_het_off, family_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")

  trans_plot <- transHeatmap(mdf_trans, family_facet) +
    scale_fill_gradient(low = "white",
                        high = "blue",
                        na.value="gray",
                        limits = c(0,NA))  +
    facet_grid(Cluster ~., scales = "free_y", space = "free")

  hetDif_plot <- hetDifHeatmap(mdf_hetDif, family_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")

  all_plot <- grid.arrange(hetDif_plot, het_plot_on, het_plot_off, nrow = 1, widths = c(1,1,2))

  result <- list(df = subset_all %>% arrange(Label), plot = all_plot)
  return(result)
}


##info_df %>% select(Family) %>% pull() %>% unique()

aMAFC_th <- 2
family_facet <- F
family <- NA
nclust <- 5
df <- gm_pos_maxtrans ##%>% filter(is.na(Family))

x <- heatMap_allstrains_trans(df, aMAFC_th, family, family_facet)
ggsave('./Plots/current_genemodel_plot_with_neighbors.pdf', x$plot, device = 'pdf', height = 80, width = 60, units = 'cm')

x$df %>%
  select(-contains('bin')) %>%
  select(Label, Annot, Family)

all_gmodel_heat <- x$df

trans_df %>%
  filter(Gene_id == 'PF3D7_0713100') %>%
  print(width = 400)

## Plot by Family

for (family in info_df %>% select(Family) %>% pull() %>% unique()){

  aMAFC_th <- 1
  family_facet <- F
  nclust <- 3
  df <- gm_pos_maxtrans
  n <- dim(gm_pos_maxtrans %>% filter(Family == family & abs(Max_aMAFC) > aMAFC_th))[1]

  if (n >= nclust & !is.na(family)){
    plotname <- paste0('./Plots/Gene_Model/genemodel_positive_complete_', family, '.svg')
    csvname <- paste0('./Gene_Model_Tables/', family, '_log2trans_', aMAFC_th, '.csv')

    x <- heatMap_allstrains_trans(gm_pos_maxtrans, aMAFC_th, family, family_facet)

    ggsave(plotname, x$plot, device = 'svg')
    write_csv(x$df %>% select(-contains('bin')), csvname)
  }
}

## Plot the "NA" Family

aMAFC_th <- 2
family_facet <- F
family <- NA
nclust <- 3
df <- gm_pos_maxtrans %>%
  filter(is.na(Family) & !Is_tRNA)
x <- heatMap_allstrains_trans(df, aMAFC_th, family, family_facet)
plotname <- paste0('./Plots/Gene_Model/genemodel_positive_complete_', family, '.svg')
csvname <- paste0('./Gene_Model_Tables/', family, '_log2trans_', aMAFC_th, '.csv')
ggsave(plotname, x$plot, device = 'svg')
write_csv(x$df %>% select(-contains('bin')), csvname)

test_kmeans <- function(mtx){

  ## Set max k for k analysis
  maxk <- mtx %>% distinct() %>% nrow()
  maxk <- min(maxk, 20)

  ## Elbow method, keep reducing maxk if it fails
  while (maxk > 0) {
    test <- try(
      ## elbow plot
      elbow <- fviz_nbclust(mtx, kmeans, method = "wss", k.max = maxk) +
        labs(subtitle = "Elbow method")
    , silent = T)
    ## If kmax throws an error, reduce it by 1
    if (class(test)[1] != "gg") {
      maxk <- maxk -1
    } else {
      break
    }
  }

  ## Silhouette method (if kmax works for elbow it will work for silhouette)
  silhouette <- fviz_nbclust(mtx, kmeans, method = "silhouette", k.max = maxk) +
    labs(subtitle = "Silhouette method")

  result = list(elbow = elbow, silhouette = silhouette)
}

heatMap_allstrains_kmeans <- function(df, aFC_th, fam, fam_facet, nclu, tbar, fbar){

  ## Returns a list object with list = (df, plot)

  subset_all <- df %>%
    filter(abs(Max_aMAFC) > aFC_th)

  if (!is.na(fam)){
    subset_all <- subset_all %>%
      filter(Family == fam)
  }

  subset_all <- subset_all %>%
    filter(across(
      contains('bin'),
      complete.cases
    ))

  mtx <- subset_all %>%
    select(contains('bin'),
           -contains('bin1_'),
           -contains('bin2_'),
           -contains('bin3_'),
           -contains('bin4_'),
           -contains('bin20_'),
           -contains('bin21_'),
           -contains('bin22_'),
           -contains('bin23_'),
           )
    #select(contains('On'), contains('Off'))
    #select(contains('Dif'))

  ## Make K-means Clustering

  km_fit <- kmeans(mtx, nclu, iter.max = 1000)
  cl <- km_fit$cluster
  subset_all['Cluster'] <- cl

  ## Test optimal k's
  k_tests <- test_kmeans(mtx)
  elbow <- k_tests$elbow
  silhouette <- k_tests$silhouette

  ## Subset DF for different heatmaps

  subset_trans <- subset_all %>%
    select(Gene_id, Max_aMAFC, Label, Name, Annot, Family, SubFamily, Cluster)

  subset_het_on <- subset_all %>%
    select(Gene_id, Label,
           contains('On') & contains('bin'),
           Name, Annot, Family, SubFamily, Cluster)

  subset_het_off <- subset_all %>%
    select(Gene_id, Label,
           contains('Off') & contains('bin'),
           Name, Annot, Family, SubFamily, Cluster)

  subset_hetDif <- subset_all %>%
    select(Gene_id, contains('_Dif'), Label, Name, Annot, Family, SubFamily, Cluster)

  subset_fambar <- subset_all %>%
    select(Gene_id, Label,
           Family_Grouped,
           Name, Annot, Family, SubFamily, Cluster) %>%
    #mutate(Label = factor(Label, levels = Label)) %>%
    replace_na(list(Family_Grouped = 'Non CVGs'))

  ## Melt Data-Frames

  melt_vars <- c('Gene_id', 'Label', 'Name', 'Annot', 'Family', 'SubFamily', 'Cluster')
  mdf_het_on <- melt(subset_het_on, id.vars = melt_vars)
  mdf_het_off <- melt(subset_het_off, id.vars = melt_vars)
  mdf_trans <- melt(subset_trans, id.vars = melt_vars)
  mdf_hetDif <- melt(subset_hetDif, id.vars = melt_vars)
  mdf_fambar <- melt(subset_fambar, id.vars = melt_vars)

  head(mdf_fambar)

  ## Create individual Heatmaps

  het_plot_on <- hetHeatmap(mdf_het_on, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free") +
    theme(axis.text.x = element_blank())

  het_plot_off <- hetHeatmap(mdf_het_off, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")+
    theme(axis.text.x = element_blank())

  trans_plot <- transHeatmap(mdf_trans, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free") +
    theme(axis.text.x = element_blank()) +
    scale_fill_gradient(low = "white",
                        high = "blue",
                        na.value="gray",
                        limits = c(0,NA))

  hetDif_plot <- hetDifHeatmap(mdf_hetDif, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free") +
    theme(axis.text.x = element_blank())

  fambar_plot <- family_heatmap(mdf_fambar) +
    facet_grid(Cluster ~., scales = "free_y", space = "free")

  ## Create fake Heatmap for Labels

  labels_df <- mdf_trans %>%
    select(Label, Family, Cluster) %>%
    mutate(value = 1, variable = 'fake_val')

  labels_plot <- hetHeatmap(labels_df, fam_facet) +
    facet_grid(Cluster ~., scales = "free_y", space = "free") +
    theme(
      axis.text.y = element_text(),
      axis.ticks.y = element_line(),
      axis.text.x = element_blank()
    )

  ## Arrange plots (with or without transcription heatmap and fambar)

  pl_wd <- tibble(P_names = c('trans', 'het_dif', 'het_on', 'het_off', 'fambar', 'labels'),
                  Plots = list(
                       trans_plot,
                       hetDif_plot,
                       het_plot_on,
                       het_plot_off,
                       fambar_plot,
                       labels_plot),
                  widths = c(0.1,1,1,1,0.1,1))

  if (!tbar){pl_wd <- pl_wd %>% filter(P_names != 'trans')}
  if (!fbar){pl_wd <- pl_wd %>% filter(P_names != 'fambar')}

  all_plot <- grid.arrange(grobs = pl_wd$Plots, nrow = 1, widths = pl_wd$widths)

  ## Create output
  result <- list(df = subset_all, plot = all_plot, elbow = elbow, silhouette = silhouette)
  return(result)
}



## set.seed(123)

## x <- heatMap_allstrains_kmeans(
##   df = gm_pos_maxtrans,
##   aFC_th = 2,
##   fam = NA,
##   fam_facet = F,
##   nclu = 7,
##   tbar = F,
##   fbar = T
## )

## Current Heatmap

set.seed(123)

x <- heatMap_allstrains_kmeans(
  df = gm_pos_maxtrans,
  aFC_th = 2,
  fam = NA,
  fam_facet = F,
  nclu = 7,
  tbar = F,
  fbar = T
)

## Save resulting Heatmap and DF

ggsave('./Plots/current_genemodel_plot_with_neighbors_Kmeans_021121.pdf', x$plot, device = 'pdf', height = 80, width = 60, units = 'cm')

all_gmodel_kmeans <- x$df

all_gmodel_kmeans %>%
  select(-contains('bin')) %>%
  write_tsv('current_kmeans_clusters.tsv')

## Check output
x$elbow
x$silhouette
plot(x$plot)

## Some tests
cl1 <- all_gmodel_kmeans %>%
  filter(Cluster == 1)

cl1$Label[order(cl1$Label)]

cl1$Label
cl1 %>%
  arrange(desc(Label)) %>%
  select(Label)


all_gmodel_kmeans %>%
  select(-contains('bin')) %>%
  group_by(Cluster) %>%
  summarize(Mean_aMAFC = mean(Max_aMAFC))

## famdf <- all_gmodel_kmeans %>%
##   arrange(Cluster, desc(Label))

## famdf <- famdf %>%
##   mutate(Label = factor(Label, levels = Label)) %>%
##   select(Label, Family_Grouped) %>%
##   replace_na(list(Family = 'Non CVGs'))

## mfamdf <- melt(famdf, id='Label')

## fam_bar <- ggplot(mfamdf, aes(x = variable, y = Label)) +
##   geom_tile(aes(fill = value)) +
##   #geom_text(aes(label=value)) +
##   my_scale +
##   scale_y_discrete(limits=(rev(levels(famdf$Label)))) +
##   theme(
##       panel.border=element_blank(),
##       panel.grid.major=element_blank(),
##       #strip.background = element_blank(),
##       strip.text.x = element_blank(),
##     )

## fam_bar

## ggsave('./Plots/Kmeans/families_bar.pdf', fam_bar, device = 'pdf')

## Plot by Family

set.seed(123)

info_df %>%
  count(Family)

fam_plot <- heatMap_allstrains_kmeans(
  df = gm_pos_maxtrans %>%
    mutate(Label = paste(Gene_id, SubFamily, sep = ': ')),
  aFC_th = 0,
  fam = 'CLAG',
  fam_facet = F,
  nclu = 3,
  tbar = T,
  fbar = F
)


fam_plot$elbow
fam_plot$silhouette
plot(fam_plot$plot)

## Current famillies and k's
fams_ks <- list(
  c('ACS', 3),
  c('CLAG', 2),
  c('HYP', 3),
  c('OTHER', 7),
  c('PFMC-2TM', 3),
  c('PHIST', 6),
  c('RIFIN', 5),
  c('STEVOR', 3),
  c('VAR', 5)
)

for (fam_k in fams_ks){

  print(fam_k)

  fam_plot <- heatMap_allstrains_kmeans(
    df = gm_pos_maxtrans %>%
      mutate(Label = paste(Gene_id, SubFamily, sep = ': ')),
    aFC_th = 0,
    fam = fam_k[1],
    fam_facet = F,
    nclu = as.numeric(fam_k[2]),
    tbar = T,
    fbar = F
  )

  ggsave(paste0('./Plots/Kmeans/Families/', fam_k[1], '_aMAFC', aMAFC_th, '_k', fam_k[2], '.pdf'),
         fam_plot$plot, width = 60, units = 'cm', device = 'pdf')
  ggsave(paste0('./Plots/Kmeans/Families/', fam_k[1], '_elbow', '.pdf'),
         fam_plot$elbow, device = 'pdf')
  ggsave(paste0('./Plots/Kmeans/Families/', fam_k[1], '_silhouette', '.pdf'),
         fam_plot$silhouette, device = 'pdf')
}



for (family in info_df %>% select(Family) %>% pull() %>% unique()){

  #family <- 'RIFIN'
  aMAFC_th <- 0
  family_facet <- F
  df <- gm_pos_maxtrans
  n <- dim(gm_pos_maxtrans %>% filter(Family == family & abs(Max_aMAFC) > aMAFC_th))[1]
  nclust <- n%/%10

  if (n >= nclust & !is.na(family)){
    plotname <- paste0('./Plots/Gene_Model/genemodel_positive_complete_', family, '.svg')
    csvname <- paste0('./Gene_Model_Tables/', family, '_log2trans_', aMAFC_th, '.csv')

    x <- heatMap_allstrains_kmeans(gm_pos_maxtrans, aMAFC_th, family, family_facet, nclust)
    #ggsave(plotname, x$plot, device = 'svg')
    #write_csv(x$df %>% select(-contains('bin')), csvname)
  }
}

## Plot the "NA" Family

aMAFC_th <- 2
family_facet <- F
family <- NA
nclust <- 3
df <- gm_pos_maxtrans %>%
  filter(is.na(Family) & !Is_tRNA)
x <- heatMap_allstrains_trans(df, aMAFC_th, family, family_facet)
plotname <- paste0('./Plots/Gene_Model/genemodel_positive_complete_', family, '.svg')
csvname <- paste0('./Gene_Model_Tables/', family, '_log2trans_', aMAFC_th, '.csv')
ggsave(plotname, x$plot, device = 'svg')
write_csv(x$df %>% select(-contains('bin')), csvname)

tendency_plot_loess <- function(df, cluster){

  max_y <- max(df %>% select(contains('_On') | contains('_Off')))
  min_y <- min(df %>% select(contains('_On') | contains('_Off')))

  plt_df <-  df %>%
    select(Gene_id, contains('On'), contains('Off'), Cluster) %>%
    filter(Cluster == cluster) %>%
    select(-Cluster)

  plot_mdf <- melt(plt_df)
  head(plot_mdf)

  plot_mdf['State'] <- sapply(plot_mdf$variable, function(x) str_split(x, '_')[[1]][2])
  plot_mdf <- plot_mdf %>%
    mutate(variable = as.numeric(sub("bin(\\d+).*", "\\1", variable)))

  head(plot_mdf)

  p <- ggplot(plot_mdf, aes(x = variable, y = value, group = State))
  p <- p + geom_smooth(aes(color = State), method = 'loess', level = 0.95)
  p <- p + scale_color_manual(values = c('red', 'green'))
  p <- p + ggtitle(paste0('Cluster ', cluster))
  p <- p + geom_vline(xintercept = 4.5, linetype="solid")
  p <- p + geom_vline(xintercept = 9.5, linetype="dotted")
  p <- p + geom_vline(xintercept = 14.5, linetype="dotted")
  p <- p + geom_vline(xintercept = 19.5, linetype="solid")
  p <- p + theme_classic()
  p <- p + theme(axis.title.x = element_blank())
  p <- p + ylab('H3K9me3 enrichment')
  p <- p + scale_x_continuous(
             breaks=c(2.5,4.5,7,12,17,19.5, 21.5),
             labels=c('Prev\nGenes', "ATG", "5'UTR", "CDS", "3'UTR", "END", 'Post\nGenes'),
             expand = c(0, 0)
           )
  p <- p + scale_y_continuous(limits = c(min_y-0.3, max_y+0.3), expand = c(0, 0))
  p
}

tendency_plot_mean_sd <- function(df, cluster){

  onoff_df <- df %>%
    select(-contains('_Dif'))

  mean_df <- onoff_df %>%
    group_by(Cluster) %>%
    summarize(across(contains('bin'), list(mean))) %>%
    pivot_longer(!Cluster, values_to = 'mean')

  sd_df <- onoff_df %>%
    group_by(Cluster) %>%
    summarize(across(contains('bin'), list(sd))) %>%
    pivot_longer(!Cluster, values_to = 'sd')

  plot_df <- mean_df %>%
    full_join(sd_df, by = c('Cluster', 'name'))

  max_y <- max(plot_df$mean + plot_df$sd)
  min_y <- min(plot_df$mean - plot_df$sd)

  plot_df <- plot_df %>%
    filter(Cluster == cluster) %>%
    select(-Cluster)

  plot_df['State'] <- sapply(plot_df$name, function(x) str_split(x, '_')[[1]][2])
  plot_mdf <- plot_df %>%
    mutate(variable = as.numeric(sub("bin(\\d+).*", "\\1", name)))

  #print(min_y)

  p <- ggplot(plot_mdf, aes(x = variable, group = State))
  p <- p + geom_line(aes(y = mean, color = State), size = 1)
  p <- p + geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = State),
                       alpha = .2)
  p <- p + scale_color_manual(values = c('red', 'green'))
  p <- p + ggtitle(paste0('Cluster ', cluster))
  p <- p + geom_vline(xintercept = 4.5, linetype="solid")
  p <- p + geom_vline(xintercept = 9.5, linetype="dotted")
  p <- p + geom_vline(xintercept = 14.5, linetype="dotted")
  p <- p + geom_vline(xintercept = 19.5, linetype="solid")
  p <- p + theme_classic()
  p <- p + theme(axis.title.x = element_blank())
  p <- p + ylab('H3K9me3 enrichment')
  p <- p + scale_x_continuous(
             breaks=c(2.5,4.5,7,12,17,19.5, 21.5),
             labels=c('Prev\nGenes', "ATG", "5'UTR", "CDS", "3'UTR", "END", 'Post\nGenes'),
             expand = c(0, 0)
           )
  p <- p + scale_y_continuous(limits = c(min_y, max_y), expand = c(0, 0))
  p
}


## 1 cluster plots
tendency_plot_loess(all_gmodel_heat, 6)
tendency_plot_loess(all_gmodel_kmeans, 6)
tendency_plot_mean_sd(all_gmodel_kmeans, 6)

## All cluster plots, overlapped

## Hierarchical clustering loess
nclust <- length(unique(all_gmodel_heat$Cluster))

plots <- lapply(1:nclust, function(x) tendency_plot_loess(all_gmodel_heat, x))
clusters_plot <- grid.arrange(grobs = plots, nrow = length(plots), ncol = 1)

ggsave(
  paste0('./Plots/clusters_plot_overlapped_with_neighbors.pdf'),
  clusters_plot,
  height = 80, width = 15, units = 'cm',
  device = 'pdf'
)


## K-means loess
plots <- lapply(1:nclust, function(x) tendency_plot_loess(all_gmodel_kmeans, x))
clusters_plot <- grid.arrange(grobs = plots, nrow = length(plots), ncol = 1)

ggsave(
  paste0('./Plots/clusters_plot_overlapped_with_neighbors_kmeans.pdf'),
  clusters_plot,
  height = 80, width = 15, units = 'cm',
  device = 'pdf'
)

## K-means means +- sd
plots <- lapply(1:nclust, function(x) tendency_plot_mean_sd(all_gmodel_kmeans, x))
clusters_plot <- grid.arrange(grobs = plots, nrow = length(plots), ncol = 1)

ggsave(
  paste0('./Plots/clusters_plot_overlapped_with_neighbors_kmeans_means_sd.pdf'),
  clusters_plot,
  height = 80, width = 15, units = 'cm',
  device = 'pdf'
)

box_df <- all_gmodel_kmeans %>%
  select(Gene_id, Max_aMAFC, Cluster)

p <- ggplot(all_gmodel_kmeans, aes(x = as.character(Cluster), y= Max_aMAFC))
p <- p + stat_boxplot(geom = "errorbar", width = 0.5) #add perpendicular whiskers
p <- p + geom_boxplot(outlier.colour="red")
p <- p + theme_classic()
p <- p + xlab('Cluster')
p <- p + ylab('Level of expression (aMAFC)')
p

ggsave(paste0('./Plots/cluster_boxplots.pdf'), p, device = 'pdf')

## Donut Plot

## ## Create 'fixed' color palette
## all_fams <- c('Not CVGs',
##               'Other CVGs',
##               '6-CYS',
##               'ACBP',
##               'ACS',
##               'CLAG',
##               'DNAJ',
##               'GBP',
##               'PFMC-2TM',
##               'SURFIN',
##               'FIKK',
##               'HYP',
##               'PHIST',
##               'STEVOR',
##               'RIFIN',
##               'VAR')

## col_factor <- factor(all_fams, levels=all_fams)

## my_colors <- c('gray', scales::viridis_pal()(15))
## names(my_colors) <- levels(col_factor)
## my_scale_dount <- scale_fill_manual(name = "Family_Grouped", values = my_colors)


## donut_df <- all_gmodel_heat %>%
##   mutate(Family_Dount = case_when(
##            Family %in% col_factor ~ Family,
##            Family == 'OTHER' ~ 'Other CVGs',
##            is.na(Family) ~ 'Not CVGs'))


donut_plot <- function(df, cluster) {
    x <- df %>%
      filter(Cluster == cluster)

    data <- as.data.frame(table(x$Family_Grouped))
    colnames(data) <- c('category', 'count')

    count(info_df, Family)

    ## Compute percentages
    data$fraction = data$count / sum(data$count)

    ## Compute the cumulative percentages (top of each rectangle)
    data$ymax = cumsum(data$fraction)

    ## Compute the bottom of each rectangle
    data$ymin = c(0, head(data$ymax, n=-1))

    ## Compute label position
    data$labelPosition <- (data$ymax + data$ymin) / 2

    ## Compute a good label
    data$label <- paste0(data$category, "\n value: ", data$count)

    ## Make the plot
    p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
      geom_rect(color="white") +
      #geom_label( x=4.2, aes(y=labelPosition, label=label), size=3) +
      coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
      xlim(c(2, 4)) +# Try to remove that to see how to make a pie chart
      theme_void() + my_scale

    outname <- paste0('./Plots/Donuts/gmodel_hetmap_cluster_',
                      as.character(cluster),
                      '.svg')
    #ggsave(outname, p, device = 'svg')
    p
}

count(all_gmodel_kmeans, Cluster)

p1 <- donut_plot(all_gmodel_heat, 1)
p2 <- donut_plot(all_gmodel_heat, 2)
p3 <- donut_plot(all_gmodel_heat, 3)
p4 <- donut_plot(all_gmodel_heat, 4)
p5 <- donut_plot(all_gmodel_heat, 5)
p6 <- donut_plot(all_gmodel_heat, 6)
p7 <- donut_plot(all_gmodel_heat, 7)

all_donuts <- grid.arrange(p1, p2, p3, p4, p5, p6, p7,  nrow = 7, ncol = 1)
ggsave('./Plots/Donuts_with_neighbors/gmodel_heatmap_alldonuts.svg', all_donuts, device = 'svg')


## Plot all genes to get legend
## all_genes_donut <- info_df %>%
##   mutate(Family_Dount = case_when(
##            Family %in% col_factor ~ Family,
##            Family == 'OTHER' ~ 'Other CVGs',
##            is.na(Family) ~ 'Not CVGs'))


data <- as.data.frame(table(info_df$Family_Grouped))
colnames(data) <- c('category', 'count')

count(info_df, Family)

## Compute percentages
data$fraction = data$count / sum(data$count)

## Compute the cumulative percentages (top of each rectangle)
data$ymax = cumsum(data$fraction)

## Compute the bottom of each rectangle
data$ymin = c(0, head(data$ymax, n=-1))

## Compute label position
data$labelPosition <- (data$ymax + data$ymin) / 2

## Compute a good label
data$label <- paste0(data$category, "\n value: ", data$count)

## Make the plot
p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect(color="white") +
                                        #geom_label( x=4.2, aes(y=labelPosition, label=label), size=3) +
  coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
  xlim(c(2, 4)) +# Try to remove that to see how to make a pie chart
  theme_void() + my_scale

ggsave('./Plots/Donuts/all_genes_donut.svg', p, device = 'svg')
p


## Create 'fixed' color palette
col_factor <- factor(unique(analysis_df$Family_Grouped),
                     levels=c('Not CVGs',
                              'Other CVGs',
                              'FIKK',
                              'HYP',
                              'PHIST',
                              'STEVOR',
                              'RIFIN',
                              'VAR'))

my_colors <- c('gray', scales::viridis_pal(begin = 0, end = 1)(7))
names(my_colors) <- levels(col_factor)
my_scale <- scale_fill_manual(name = "Family_Grouped", values = my_colors)

all_gmodel_heat %>%
  filter(Cluster == 1) %>%
  select(Gene_id) %>%
  left_join(transFC_df_noNAs) %>%
  left_join(red_df) %>%
  select(Gene_id, Max_aMAFC, On_trans, Off_trans, Is_3D7B, A7, E5, B11) %>%
  arrange(Gene_id) %>%
  print(n = 50)



red_df



transFC_df_noNAs

## Load gene state tables

state_old <- read_csv('/mnt/Disc4T/Projects/Active_gene_list/Results_Tables/gene_state_final.csv') %>%
  select(Gene_id, category_12B, category_10G)

state_new <- read_tsv('/mnt/Disc4T/Projects/Active_gene_list/New_Arrays_Results/state_df_rna25_red25_redresc40_reddw15_areaFC1_redpcntdif30_allstrains.csv')%>%
  select(Gene_id, category_A7, category_E5, category_B11)

state_df <- state_old %>%
  full_join(state_new, by = 'Gene_id') %>%
  rename_with(~ gsub('category_', 'Gene_State_', .x))


v <- c('Var_Repressed', 'Var_Repressed')
my_difstate_filter(v)

silenced <- any(grepl('Var_Repressed', state_vect)) | any(grepl('Inactive', state_vect))



## Check for which genes we have a strain in an active state and a strain in a silenced state
my_difstate_filter <- function(state_vect){
  active <- any(grepl('Active', state_vect))
  silenced <- any(grepl('Var_Repressed', state_vect)) | any(grepl('Inactive', state_vect))
  return(active & silenced)
}



state_df <- state_df %>%
  mutate(DifState = apply(select(., contains('State')), 1, my_difstate_filter))

stinfo_df <- state_df %>%
  left_join(info_df, by='Gene_id')

stinfo_df %>%
  filter(DifState) %>%
  select(Gene_id, contains('State'), Variant, Annot) %>%
  print(width = 200)


stinfo_df %>%
  filter(DifState) %>%
  count(Variant)

vars <- stinfo_df %>%
  filter(DifState) %>%
  filter(Variant) %>%
  select(Gene_id, contains('State'), Annot) %>%
  print(width = 200)

novars <- stinfo_df %>%
  filter(DifState) %>%
  filter(!Variant) %>%
  select(Gene_id, contains('State'), Annot) %>%
  print(width = 200)

vh <- gm_pos_maxtrans %>%
  filter(Gene_id %in% vars$Gene_id)

nvh <- gm_pos_maxtrans %>%
  filter(Gene_id %in% novars$Gene_id)


x <- heatMap_allstrains_kmeans(
  df = nvh,
  aFC_th = 0,
  fam = NA,
  fam_facet = F,
  nclu = 1,
  tbar = F,
  fbar = T
)

x$silhouette

#### Load Data ####

## ## Remove neighboring genes bins
## gmodel_noneighbors <- gmodel_all %>%
##   select(
##     -contains('bin1_'),
##     -contains('bin2_'),
##     -contains('bin3_'),
##     -contains('bin4_'),
##     -contains('bin20_'),
##     -contains('bin21_'),
##     -contains('bin22_'),
##     -contains('bin23_'),
##     )

## bins <- paste0('bin', c(1:15))
## all_cols <- c(paste0(bins, '_12B'),
##               paste0(bins, '_10G'),
##               paste0(bins, '_A7'),
##               paste0(bins, '_E5'),
##               paste0(bins, '_B11')
##               )
## gmodel_noneighbors <- gmodel_noneighbors %>%
##   set_names('Gene_id', all_cols)

## nona_gm_noneighbors <- gmodel_noneighbors %>%
##   filter(complete.cases(gmodel_noneighbors))

## gm_dif_heatmap <- function(contrast, trans_th, difpeaks, family_facet){
##   difpeaks_col <- paste0('Difpeaks_', contrast[1], '_', contrast[2])
##   trans_col <- paste0(contrast[1], '-', contrast[2], '_MaxVal')

##   gm_dif <- nona_gm_noneighbors %>%
##     left_join(info_df) %>%
##     left_join(trans_df) %>%
##     filter(abs(get(trans_col)) > trans_th)

##   if (difpeaks) {gm_dif <- gm_dif %>% filter(get(difpeaks_col))}

##   x <- gm_dif %>%
##     select(contains(contrast[1]) & contains('bin'))

##   y <- gm_dif %>%
##     select(contains(contrast[2]) & contains('bin'))

##   mtx <- abs(x-y)

##   out_df <- gm_dif %>%
##     select(-contains('bin'), -contains('MaxVal'), -contains('MaxTime')) %>%
##     bind_cols(mtx)

##   vals <- out_df %>%
##     select(contains('bin')) %>%
##     pull()

##   ## Make hierarquical Clustering
##   dmtx <- dist(mtx, method = "euclidean")
##   cl <- hclust(dmtx, method = 'average')
##   out_df$Gene_id <- factor(out_df$Gene_id, levels = out_df$Gene_id[cl$order])
##   customHeatmap(out_df, family_facet, c(min(vals),max(vals)))
## }

## contrast <- c('12B', '10G')
## trans_th <- 1
## difpeaks <- F
## family_facet <- F

## gm_dif_heatmap(contrast, trans_th, difpeaks, family_facet)

#difpeak_mtx <- mtx_12b10g[rownames(mtx_12b10g) %in% hetdifgenes,]
#heatDF <- cbind(hetdifDF[,-4], difpeak_mtx)

## gm_sidebyside_heatmap <- function(contrast, family_facet, trans_th){
##   trans_col <- paste0(contrast[1], '-', contrast[2], '_MaxVal')

##   gm_df <- nona_gm_noneighbors %>%
##     select(Gene_id, contains(contrast)) %>%
##     left_join(info_df, by = 'Gene_id') %>%
##     left_join(trans_df, by = 'Gene_id') %>%
##     filter(abs(get(trans_col)) > trans_th)

##   gm_mtx <- gm_df %>%
##     select(contains('bin'))


##   ## Make hierarquical Clustering
##   dmtx <- dist(gm_mtx, method = "euclidean")
##   cl <- hclust(dmtx, method = 'complete')
##   gm_df$Gene_id <- factor(gm_df$Gene_id, levels = gm_df$Gene_id[cl$order])

##   #customHeatmap(gm_df, c(min(gm_mtx), max(gm_mtx)))

##   #head(gm_df)
##   mdf <- melt(gm_df %>% select(-c(contains('MaxVal'), contains('MaxTime'))))
##   mdf["Strain"] <- sapply(mdf$variable, function(x) strsplit(as.character(x), split = "_", fixed = TRUE)[[1]][2])

##   p <- ggplot(mdf, aes(x = variable, y = Gene_id, fill = value)) +

##     geom_tile(colour="snow3") +
##               #size=0.10,
##               #height=.9) +

##     scale_fill_gradient2(midpoint = 0,
##                          low = "white",
##                          high = "red") +

##     scale_y_discrete(position = "right") +

##     theme(
##       strip.background = element_blank(),
##       axis.title = element_blank(),
##       #strip.text.x = element_blank(),
##       axis.line.x = element_blank(),
##       axis.ticks.x = element_blank(),
##       axis.text.x = element_blank()) +

##     facet_grid(~Strain,
##                scales="free_x",
##                space="free")

##   if (family_facet){p <- p+facet_grid(vars(Family), vars(Strain), scales = "free", space = "free")}

##   p

## }

## contrast <- c('E5', 'B11')
## family_facet <- F
## trans_th <- 1.5

## gm_sidebyside_heatmap(contrast, family_facet, trans_th)

## tdf <- transFC_df_noNAs %>%
##   select(Gene_id, Max_aMAFC, On_trans, Off_trans, Is_3D7B)

## gmodel_noneighbors_pos <- gmodel_noneighbors %>%
##   mutate(across(-Gene_id, ~ ifelse(.x > 0, .x, 0)))

## ## Create empty DF
## hetcols <- gmodel_noneighbors %>%
##   select(contains('12B') & contains('bin')) %>%
##   colnames(.) %>%
##   gsub('12B', '', ., fixed=TRUE)

## col_names <- c('Gene_id',
##                paste0(hetcols, 'On'),
##                paste0(hetcols, 'Off'),
##                paste0(hetcols, 'Dif'))

## cn <- setNames(rep('', length(col_names)), col_names)
## df <- bind_rows(cn)[0,]

## ## Traverse max trans. df and get coverage cols

## for (i in 1:dim(tdf)[1]){

##   gid <- as.character(tdf$Gene_id[i])
##   on <- tdf$On_trans[i]
##   off <- tdf$Off_trans[i]

##   ## Select contrast in which difference sign must be switched
##   ## contrast <- paste0(on, '-', off)
##   ## positive_contrasts <- c(
##   ##   '12B-10G', '12B-A7', '12B-E5', '12B-B11',
##   ##   '10G-A7', '10G-E5', '10G-B11',
##   ##   'A7-E5', 'A7-B11',
##   ##   'B11-E5'
##   ## )

##   ## ifelse(
##   ##   contrast %in% positive_contrasts,
##   ##   neg_dif <- FALSE,
##   ##   neg_dif <- TRUE
##   ## )

##   ### Main Loop
##   if (gid %in% gmodel_noneighbors_pos$Gene_id & !is.na(on) & !is.na(off)){

##     onvect <- gmodel_noneighbors_pos %>%
##       filter(Gene_id == gid) %>%
##       select(contains(on))

##     offvect <- gmodel_noneighbors_pos %>%
##       filter(Gene_id == gid) %>%
##       select(contains(off))

##     diffvect <- onvect-offvect
##     #if (neg_dif) {diffvect <- -diffvect}

##     row <- c(gid,
##              unlist(onvect, use.names = F),
##              unlist(offvect, use.names = F),
##              unlist(onvect-offvect, use.names = F))
##     row <- setNames(row, col_names)
##     df <- df %>% add_row(bind_rows(row))
##   }
## }

## df <- df %>%
##   mutate(across(-Gene_id, as.numeric))

## gm_noneighbors_pos_maxtrans <- df %>%
##   full_join(tdf, by = 'Gene_id') %>%
##   left_join(info_df) %>%
##   mutate(Max_aMAFC = abs(Max_aMAFC))

## heatMap_allstrains_trans <- function(df, aMAFC_th, family, family_facet){

##   ## Returns a list object with list = (df, plot)

##   df
##   subset_all <- df %>%
##     filter(abs(Max_aMAFC) > aMAFC_th)

##   if (!is.na(family)){
##     subset_all <- subset_all %>%
##       filter(Family == family)
##   }

##   subset_all <- subset_all %>%
##     filter(across(
##       contains('bin'),
##       complete.cases
##     ))

##   mtx <- subset_all %>%
##     select(contains('_On') | contains('_Off'))

##   ## Make hierarquical Clustering
##   dmtx <- dist(mtx, method = "euclidean")
##   cl <- hclust(dmtx, method = 'complete')
##   subset_all$Label <- factor(subset_all$Label, levels = subset_all$Label[cl$order])
##   subset_all['Cluster'] <- cutree(cl, nclust)
##   subset_all$Label <- factor(subset_all$Label, levels = subset_all$Label[cl$order])

##   subset_trans <- subset_all %>%
##     select(Gene_id, Max_aMAFC, Label, Name, Annot, Family, SubFamily, Cluster)

##   subset_het_on <- subset_all %>%
##     select(Gene_id, Label,
##            contains('On') & contains('bin'),
##            Name, Annot, Family, SubFamily, Cluster)

##   subset_het_off <- subset_all %>%
##     select(Gene_id, Label,
##            contains('Off') & contains('bin'),
##            Name, Annot, Family, SubFamily, Cluster)

##   ##write_csv(subset_het %>% select(-contains('bin')) %>% arrange(Cluster), './het_trans_genemodel.csv')

##   subset_hetDif <- subset_all %>%
##     select(Gene_id, contains('_Dif'), Label, Name, Annot, Family, SubFamily, Cluster)

##   melt_vars <- c('Gene_id', 'Label', 'Name', 'Annot', 'Family', 'SubFamily', 'Cluster')
##   mdf_het_on <- melt(subset_het_on, id.vars = melt_vars)
##   mdf_het_off <- melt(subset_het_off, id.vars = melt_vars)
##   mdf_trans <- melt(subset_trans, id.vars = melt_vars)
##   mdf_hetDif <- melt(subset_hetDif, id.vars = melt_vars)

##   het_plot_on <- hetHeatmap(mdf_het_on, family_facet) +
##     facet_grid(Cluster ~., scales = "free_y", space = "free") +
##     theme(
##       axis.text.y = element_blank(),
##       axis.ticks.y = element_blank(),

##       panel.border=element_blank(),
##       panel.grid.major=element_blank(),

##       strip.background = element_blank(),
##       strip.text.x = element_blank(),
##       strip.text.y = element_blank(),
##       )


##   het_plot_off <- hetHeatmap(mdf_het_off, family_facet) +
##     facet_grid(Cluster ~., scales = "free_y", space = "free")

##   trans_plot <- transHeatmap(mdf_trans, family_facet) +
##     scale_fill_gradient(low = "white",
##                         high = "blue",
##                         na.value="gray",
##                         limits = c(0,NA))  +
##     facet_grid(Cluster ~., scales = "free_y", space = "free")

##   hetDif_plot <- hetDifHeatmap(mdf_hetDif, family_facet) +
##     facet_grid(Cluster ~., scales = "free_y", space = "free")

##   all_plot <- grid.arrange(hetDif_plot, het_plot_on, het_plot_off, nrow = 1, widths = c(1,1,2))

##   result <- list(df = subset_all %>% arrange(Label), plot = all_plot)
##   return(result)
## }


## ##info_df %>% select(Family) %>% pull() %>% unique()

## aMAFC_th <- 2
## family_facet <- F
## family <- NA
## nclust <- 8
## df <- gm_noneighbors_pos_maxtrans ##%>% filter(is.na(Family))

## x <- heatMap_allstrains_trans(df, aMAFC_th, family, family_facet)

## x$df %>%
##   select(Gene_id, Annot, Max_aMAFC, On_trans, Off_trans, Variant, Cluster) %>%
##   filter(Cluster == 6)

## #ggsave('./Plots/current_genemodel_plot.pdf', x$plot, device = 'pdf', height = 40, width = 60, units = 'cm')

## x$df %>%
##   select(-contains('bin')) %>%
##   select(Label, Annot, Family)

## all_noneighbors_gmodel_heat <- x$df

## trans_df %>%
##   filter(Gene_id == 'PF3D7_0713100') %>%
##   print(width = 400)

## ## Plot by Family

## for (family in info_df %>% select(Family) %>% pull() %>% unique()){

##   aMAFC_th <- 1
##   family_facet <- F
##   nclust <- 3
##   df <- gm_pos_maxtrans
##   n <- dim(gm_pos_maxtrans %>% filter(Family == family & abs(Max_aMAFC) > aMAFC_th))[1]

##   if (n >= nclust & !is.na(family)){
##     plotname <- paste0('./Plots/Gene_Model/genemodel_positive_complete_noneighbors_', family, '.svg')
##     csvname <- paste0('./Gene_Model_Tables/', family, '_log2trans_', aMAFC_th, '.csv')

##     x <- heatMap_allstrains_trans(gm_pos_maxtrans, aMAFC_th, family, family_facet)

##     ggsave(plotname, x$plot, device = 'svg')
##     write_csv(x$df %>% select(-contains('bin')), csvname)
##   }
## }

## ## Plot the "NA" Family

## aMAFC_th <- 2
## family_facet <- F
## family <- NA
## nclust <- 3
## df <- gm_pos_maxtrans %>%
##   filter(is.na(Family) & !Is_tRNA)
## x <- heatMap_allstrains_trans(df, aMAFC_th, family, family_facet)
## plotname <- paste0('./Plots/Gene_Model/genemodel_positive_complete_noneighbors_', family, '.svg')
## csvname <- paste0('./Gene_Model_Tables/', family, '_log2trans_', aMAFC_th, '.csv')
## ggsave(plotname, x$plot, device = 'svg')
## write_csv(x$df %>% select(-contains('bin')), csvname)

## my_tendency_plot_sides <- function(df, cluster){

##   max_y <- max(df %>% select(contains('bin')))
##   min_y <- min(df %>% select(contains('bin')))

##   on <-  df %>%
##     select(Gene_id, contains('On'), Cluster) %>%
##     filter(Cluster == cluster) %>%
##     select(-Cluster)

##   plot_on_df <- melt(on) %>%
##     mutate(variable = as.numeric(gsub('_On', '', gsub('bin', '', variable))))

##   p_on <- ggplot(plot_on_df, aes(x = variable, y = value))
##   p_on <- p_on + geom_smooth(method = 'loess', level = 0.95, color = 'green')
##   p_on <- p_on + ylim(min_y, max_y)


##   off <-  df %>%
##     select(Gene_id, contains('Off'), Cluster) %>%
##     filter(Cluster == cluster) %>%
##     select(-Cluster)

##   plot_off_df <- melt(off) %>%
##     mutate(variable = as.numeric(gsub('_Off', '', gsub('bin', '', variable))))

##   p_off <- ggplot(plot_off_df, aes(x = variable, y = value))
##   p_off <- p_off + geom_smooth(method = 'loess', level = 0.95, color = 'red')
##   p_off <- p_off + ylim(min_y, max_y)

##   title = paste0('Cluster ', as.character(cluster))
##   all_plot <- grid.arrange(p_on, p_off, nrow = 1, widths = c(1,1), top = title)
## }

## my_tendency_plot_overlapped <- function(df, cluster){

##   max_y <- max(df %>% select(contains('bin')))
##   min_y <- min(df %>% select(contains('bin')))

##   plt_df <-  df %>%
##     select(Gene_id, contains('On'), contains('Off'), Cluster) %>%
##     filter(Cluster == cluster) %>%
##     select(-Cluster)

##   plot_mdf <- melt(plt_df)
##   head(plot_mdf)

##   plot_mdf['State'] <- sapply(plot_mdf$variable, function(x) str_split(x, '_')[[1]][2])
##   plot_mdf <- plot_mdf %>%
##     mutate(variable = as.numeric(sub("bin(\\d+).*", "\\1", variable)))

##   head(plot_mdf)

##   p <- ggplot(plot_mdf, aes(x = variable, y = value, group = State))
##   p <- p + geom_smooth(aes(color = State), method = 'loess', level = 0.95)
##   p <- p + ylim(min_y, max_y)
##   p <- p + scale_color_manual(values = c('red', 'green'))
##   p <- p + ggtitle(paste0('Cluster ', cluster))
##   p
## }

## ## 1 cluster plots
## my_tendency_plot_sides(all_noneighbors_gmodel_heat, 5)
## my_tendency_plot_overlapped(all_noneighbors_gmodel_heat, 5)

## ## All cluster plots, side-by-side
## nclust <- length(unique(all_noneighbors_gmodel_heat$Cluster))

## plots <- lapply(1:nclust, function(x) my_tendency_plot_sides(all_noneighbors_gmodel_heat, x))
## clusters_plot <- grid.arrange(grobs = plots, nrow = length(plots), ncol = 1)

## ggsave(
##   paste0('./Plots/clusters_plot_sides.pdf'),
##   clusters_plot,
##   height = 80, width = 15, units = 'cm',
##   device = 'pdf'
## )


## plots <- lapply(1:nclust, function(x) my_tendency_plot_overlapped(all_noneighbors_gmodel_heat, x))
## clusters_plot <- grid.arrange(grobs = plots, nrow = length(plots), ncol = 1)

## ggsave(
##   paste0('./Plots/clusters_plot_overlapped.pdf'),
##   clusters_plot,
##   height = 80, width = 15, units = 'cm',
##   device = 'pdf'
## )

## Donut Plot

## ## Create 'fixed' color palette
## all_fams <- c('Not CVGs',
##               'Other CVGs',
##               '6-CYS',
##               'ACBP',
##               'ACS',
##               'CLAG',
##               'DNAJ',
##               'GBP',
##               'PFMC-2TM',
##               'SURFIN',
##               'FIKK',
##               'HYP',
##               'PHIST',
##               'STEVOR',
##               'RIFIN',
##               'VAR')

## col_factor <- factor(all_fams, levels=all_fams)

## my_colors <- c('gray', scales::viridis_pal()(15))
## names(my_colors) <- levels(col_factor)
## my_scale_dount <- scale_fill_manual(name = "Family_Grouped", values = my_colors)


## donut_df <- all_gmodel_heat %>%
##   mutate(Family_Dount = case_when(
##            Family %in% col_factor ~ Family,
##            Family == 'OTHER' ~ 'Other CVGs',
##            is.na(Family) ~ 'Not CVGs'))


## donut_plot <- function(df, cluster) {
##     x <- df %>%
##       filter(Cluster == cluster)

##     data <- as.data.frame(table(x$Family_Grouped))
##     colnames(data) <- c('category', 'count')

##     count(info_df, Family)

##     ## Compute percentages
##     data$fraction = data$count / sum(data$count)

##     ## Compute the cumulative percentages (top of each rectangle)
##     data$ymax = cumsum(data$fraction)

##     ## Compute the bottom of each rectangle
##     data$ymin = c(0, head(data$ymax, n=-1))

##     ## Compute label position
##     data$labelPosition <- (data$ymax + data$ymin) / 2

##     ## Compute a good label
##     data$label <- paste0(data$category, "\n value: ", data$count)

##     ## Make the plot
##     p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
##       geom_rect(color="white") +
##       #geom_label( x=4.2, aes(y=labelPosition, label=label), size=3) +
##       coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
##       xlim(c(2, 4)) +# Try to remove that to see how to make a pie chart
##       theme_void() + my_scale

##     outname <- paste0('./Plots/Donuts/gmodel_hetmap_cluster_noneighbors',
##                       as.character(cluster),
##                       '.svg')
##     #ggsave(outname, p, device = 'svg')
##     p
## }

## count(all_noneighbors_gmodel_heat, Cluster)

## p1 <- donut_plot(all_noneighbors_gmodel_heat, 1)
## p2 <- donut_plot(all_noneighbors_gmodel_heat, 2)
## p3 <- donut_plot(all_noneighbors_gmodel_heat, 3)
## p4 <- donut_plot(all_noneighbors_gmodel_heat, 4)
## p5 <- donut_plot(all_noneighbors_gmodel_heat, 5)
## p6 <- donut_plot(all_noneighbors_gmodel_heat, 6)
## p7 <- donut_plot(all_noneighbors_gmodel_heat, 7)
## p8 <- donut_plot(all_noneighbors_gmodel_heat, 8)

## all_donuts <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8,  nrow = 8, ncol = 1)
## ggsave('./Plots/Donuts/gmodel_heatmap_alldonuts.svg', all_donuts, device = 'svg')


## ## Plot all genes to get legend
## ## all_genes_donut <- info_df %>%
## ##   mutate(Family_Dount = case_when(
## ##            Family %in% col_factor ~ Family,
## ##            Family == 'OTHER' ~ 'Other CVGs',
## ##            is.na(Family) ~ 'Not CVGs'))


## data <- as.data.frame(table(info_df$Family_Grouped))
## colnames(data) <- c('category', 'count')

## count(info_df, Family)

## ## Compute percentages
## data$fraction = data$count / sum(data$count)

## ## Compute the cumulative percentages (top of each rectangle)
## data$ymax = cumsum(data$fraction)

## ## Compute the bottom of each rectangle
## data$ymin = c(0, head(data$ymax, n=-1))

## ## Compute label position
## data$labelPosition <- (data$ymax + data$ymin) / 2

## ## Compute a good label
## data$label <- paste0(data$category, "\n value: ", data$count)

## ## Make the plot
## p <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
##   geom_rect(color="white") +
##                                         #geom_label( x=4.2, aes(y=labelPosition, label=label), size=3) +
##   coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
##   xlim(c(2, 4)) +# Try to remove that to see how to make a pie chart
##   theme_void() + my_scale

## ggsave('./Plots/Donuts/all_genes_donut.svg', p, device = 'svg')
## p


## ## Create 'fixed' color palette
## col_factor <- factor(unique(analysis_df$Family_Grouped),
##                      levels=c('Not CVGs',
##                               'Other CVGs',
##                               'FIKK',
##                               'HYP',
##                               'PHIST',
##                               'STEVOR',
##                               'RIFIN',
##                               'VAR'))

## my_colors <- c('gray', scales::viridis_pal(begin = 0, end = 1)(7))
## names(my_colors) <- levels(col_factor)
## my_scale <- scale_fill_manual(name = "Family_Grouped", values = my_colors)

#### Save/load environtment ####
#save.image('binned_coverage_021121.RData')
setwd('/mnt/Disc4T/Projects/PhD_Project/Binned_Coverage/')
#load('binned_coverage_061021.RData')
