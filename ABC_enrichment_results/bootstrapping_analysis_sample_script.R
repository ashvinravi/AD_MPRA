## Bootstrapping Analysis 

# This is a sample bootstrapping analysis designed to test variant enrichment in cell-type specific enhancer data. 
#It takes in two inputs - summary statistcs (fine-mapping results, GWAS, common variants) and a bed file (ATAC-seq, ChIP-seq, ABC scores), containing annotation information. 


# Load bed files for cell-type specific data (ABC Model Results/ATAC-seq data)
library(valr)
library(readr)
library(dplyr)

set.seed(252)

setwd('/Users/ashvinravi/Desktop/finemapping_MPRA_project/AD_MPRA/')
abc_enhancer_coords <- read_tsv('ABC_enrichment_results/ABC_enhancers/abc_enhancer_coords_hg19.tsv')

# Load list of common variants from 1000G Phase 3 Files for each chromosome - concatenate into one file 
common_variants <- data.frame()
setwd('/Users/ashvinravi/Desktop/finemapping_MPRA_project/AD_MPRA/ABC_enrichment_results/1000G_panel_hg19/bim_files/')
for (file in list.files('/Users/ashvinravi/Desktop/finemapping_MPRA_project/AD_MPRA/ABC_enrichment_results/1000G_panel_hg19/bim_files/', pattern = '*.bim')) {
  f <- read.table(file)
  print(file)
  common_variants <- rbind(common_variants, f)
}
colnames(common_variants) = c('chrom', 'RSID', 'Z', 'start', 'A1', 'A2')
common_variants$chrom <- paste('chr', common_variants$chrom, sep = '')
common_variants$end <- common_variants$start + 1 
microglia_enhancer_coords <- abc_enhancer_coords[abc_enhancer_coords$CellType == 'microglia',]
microglia_abc_coords <- microglia_enhancer_coords


# Load fine-mapped results here 

fine_mapping_results <- read.table('/Users/ashvinravi/Desktop/finemapping_MPRA_project/AD_MPRA/full_results/susie_ld_gwas_with_rsids_hg19.tsv.gz', sep = '\t', header = TRUE)
names(microglia_abc_coords)[names(microglia_abc_coords) == 'chr'] = 'chrom'
names(fine_mapping_results)[names(fine_mapping_results) == 'CHR'] = 'chrom'
names(fine_mapping_results)[names(fine_mapping_results) == 'BP'] = 'start'
fine_mapping_results$end <- fine_mapping_results$start + 1

fine_mapping_results <- fine_mapping_results %>% filter(PIP >= 0.1)

fine_mapping_results$chrom <- paste('chr', fine_mapping_results$chrom, sep = '')
microglia_intersect <- bed_intersect(microglia_abc_coords, fine_mapping_results)
microglia_intersect_results <- microglia_intersect %>% distinct(chrom, start.y, end.y, LOCUS.y, PIP.y, CREDIBLE_SET.y)

monocyte_enhancers <- abc_enhancer_coords[abc_enhancer_coords$Category == 'MNP',]
names(monocyte_enhancers)[names(monocyte_enhancers) == 'chr'] = 'chrom'

monocyte_intersect <- bed_intersect(fine_mapping_results, monocyte_enhancers)
monocyte_intersect_results <- monocyte_intersect %>% distinct(chrom, start.x, end.x, LOCUS.x, PIP.x, CellType.y)

# Calculates number of fine-mapped SNPs overlapping in a given bed file 
calculate_finemap_overlap <- function(enhancer, c) {
  finemap_overlap <- bed_intersect(fine_mapping_results, enhancer)
  counts <- data.frame(nrow(finemap_overlap %>% distinct(chrom, start.x, end.x)))
  counts$celltype <- c
  colnames(counts) = c('number_of_snps', 'CellType')
  return(counts)
}

# 1000 simulations of caluclating overlap between sampling random number of common variants (same amount as the number of fine-mapped SNPs) 
# and calculating overlap between selected enhancer 
calculate_commonvar_overlap <- function(enhancer) {
  counts = data.frame()
  for (i in 1:1000) {
    d = sample_n(common_variants, nrow(fine_mapping_results))
    overlap <- bed_intersect(d, enhancer)
    number_overlap = data.frame(nrow(overlap %>% distinct(chrom, start.x, end.x)))
    counts = rbind(counts, number_overlap)
    print(i)
  }
  return(counts)
}

# Bootstrap MNP 
MNP_enhancers <- abc_enhancer_coords[abc_enhancer_coords$Category == 'MNP',]
MNP_enhancer_categories <- data.frame(table(MNP_enhancers$CellType))

MNP_gwas = data.frame()
for (category in MNP_enhancer_categories$Var1) {
  overlap_MNP <- calculate_finemap_overlap(MNP_enhancers[MNP_enhancers$CellType == category,], category)
  MNP_gwas = rbind(MNP_gwas, overlap_MNP)
}

MNP_common_var <- data.frame()
for (category in MNP_enhancer_categories$Var1) {
  temp <- MNP_enhancers[MNP_enhancers$CellType == category,]
  overlap_MNP_common <- calculate_commonvar_overlap(temp)
  overlap_MNP_common$CellType <- category 
  colnames(overlap_MNP_common) <- c('number_of_snps', 'CellType')
  MNP_common_var = rbind(MNP_common_var, overlap_MNP_common)
}
summary_MNP <- MNP_common_var %>% group_by(CellType) %>% summarise_at(vars(number_of_snps), list(name = mean))
summary_MNP <- merge(summary_MNP, MNP_gwas, by = "CellType") 
summary_MNP$OR <- summary_MNP$number_of_snps / summary_MNP$name

write_tsv(MNP_common_var, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/MNP_common_variant_overlap.tsv')
write_tsv(summary_MNP, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/summary_MNP_results.tsv')

# Bootstrap Fibroblasts 
fibroblast_enhancers <- abc_enhancer_coords[abc_enhancer_coords$Category == 'Fibroblasts',]
fibroblast_enhancer_categories <- data.frame(table(fibroblast_enhancers$CellType))

fibroblast_gwas = data.frame()
for (category in fibroblast_enhancer_categories$Var1) {
  overlap_fibroblast <- calculate_finemap_overlap(fibroblast_enhancers[fibroblast_enhancers$CellType == category,], category)
  fibroblast_gwas = rbind(fibroblast_gwas, overlap_fibroblast)
}

fibroblast_common_var <- data.frame()
for (category in fibroblast_enhancer_categories$Var1) {
  temp <- fibroblast_enhancers[fibroblast_enhancers$CellType == category,]
  overlap_fibroblast_common <- calculate_commonvar_overlap(temp)
  overlap_fibroblast_common$CellType <- category 
  colnames(overlap_fibroblast_common) <- c('number_of_snps', 'CellType')
  fibroblast_common_var = rbind(fibroblast_common_var, overlap_fibroblast_common)
}
summary_fibroblasts <- fibroblast_common_var %>% group_by(CellType) %>% summarise_at(vars(number_of_snps), list(name = mean))
summary_fibroblasts <- merge(summary_fibroblasts, fibroblast_gwas, by = "CellType") 
summary_fibroblasts$OR <- summary_fibroblasts$number_of_snps / summary_fibroblasts$name
write_tsv(fibroblast_common_var, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/fibroblast_common_variant_overlap.tsv')
write_tsv(summary_fibroblasts, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/summary_fibroblast_results.tsv')

# Bootstrap B-cells
b_cell_enhancers <- abc_enhancer_coords[abc_enhancer_coords$Category == 'B-cell',]
b_cell_enhancer_categories <- data.frame(table(b_cell_enhancers$CellType))

b_cell_gwas = data.frame()
for (category in b_cell_enhancer_categories$Var1) {
  overlap_b_cell <- calculate_finemap_overlap(b_cell_enhancers[b_cell_enhancers$CellType == category,], category)
  b_cell_gwas = rbind(b_cell_gwas, overlap_b_cell)
}

b_cell_common_var <- data.frame()
for (category in b_cell_enhancer_categories$Var1) {
  temp <- b_cell_enhancers[b_cell_enhancers$CellType == category,]
  overlap_b_cell_common <- calculate_commonvar_overlap(temp)
  overlap_b_cell_common$CellType <- category 
  colnames(overlap_b_cell_common) <- c('number_of_snps', 'CellType')
  b_cell_common_var = rbind(b_cell_common_var, overlap_b_cell_common)
}
summary_b_cell <- b_cell_common_var %>% group_by(CellType) %>% summarise_at(vars(number_of_snps), list(name = mean))
summary_b_cell <- merge(summary_b_cell, b_cell_gwas, by = "CellType") 
summary_b_cell$OR <- summary_b_cell$number_of_snps / summary_b_cell$name

write_tsv(b_cell_common_var, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/b_cell_common_variant_overlap.tsv')
write_tsv(summary_b_cell, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/summary_b_cell_results.tsv')

# Bootstrap T-cells
t_cell_enhancers <- abc_enhancer_coords[abc_enhancer_coords$Category == 'T-cell',]
t_cell_enhancer_categories <- data.frame(table(t_cell_enhancers$CellType))

t_cell_gwas = data.frame()
for (category in t_cell_enhancer_categories$Var1) {
  overlap_t_cell <- calculate_finemap_overlap(t_cell_enhancers[t_cell_enhancers$CellType == category,], category)
  t_cell_gwas = rbind(t_cell_gwas, overlap_t_cell)
}

t_cell_common_var <- data.frame()
for (category in t_cell_enhancer_categories$Var1) {
  temp <- t_cell_enhancers[t_cell_enhancers$CellType == category,]
  overlap_t_cell_common <- calculate_commonvar_overlap(temp)
  overlap_t_cell_common$CellType <- category 
  colnames(overlap_t_cell_common) <- c('number_of_snps', 'CellType')
  t_cell_common_var = rbind(t_cell_common_var, overlap_t_cell_common)
}
summary_t_cell <- t_cell_common_var %>% group_by(CellType) %>% summarise_at(vars(number_of_snps), list(name = mean))
summary_t_cell <- merge(summary_t_cell, t_cell_gwas, by = "CellType") 
summary_t_cell$OR <- summary_t_cell$number_of_snps / summary_t_cell$name

write_tsv(t_cell_common_var, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/t_cell_common_variant_overlap.tsv')
write_tsv(summary_t_cell, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/summary_t_cell_results.tsv')

# Bootstrap Epithelial 
epithelial_enhancers <- abc_enhancer_coords[abc_enhancer_coords$Category == 'Epithelial',]
epithelial_enhancer_categories <- data.frame(table(epithelial_enhancers$CellType))

epithelial_gwas = data.frame()
for (category in epithelial_enhancer_categories$Var1) {
  overlap_epithelial <- calculate_finemap_overlap(epithelial_enhancers[epithelial_enhancers$CellType == category,], category)
  epithelial_gwas = rbind(epithelial_gwas, overlap_epithelial)
}

epithelial_common_var <- data.frame()
for (category in epithelial_enhancer_categories$Var1) {
  temp <- epithelial_enhancers[epithelial_enhancers$CellType == category,]
  overlap_epithelial_common <- calculate_commonvar_overlap(temp)
  overlap_epithelial_common$CellType <- category 
  colnames(overlap_epithelial_common) <- c('number_of_snps', 'CellType')
  epithelial_common_var = rbind(epithelial_common_var, overlap_epithelial_common)
}
summary_epithelial <- epithelial_common_var %>% group_by(CellType) %>% summarise_at(vars(number_of_snps), list(name = mean))
summary_epithelial <- merge(summary_epithelial, epithelial_gwas, by = "CellType") 
summary_epithelial$OR <- summary_epithelial$number_of_snps / summary_epithelial$name

write_tsv(epithelial_common_var, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/epithelial_common_variant_overlap.tsv')
write_tsv(summary_epithelial, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/summary_epithelial_results.tsv')

# Bootstrap Hepatocytes
hepatocyte_enhancers <- abc_enhancer_coords[abc_enhancer_coords$Category == 'Other_blood',]
hepatocyte_enhancer_categories <- data.frame(table(hepatocyte_enhancers$CellType))

hepatocyte_gwas = data.frame()
for (category in hepatocyte_enhancer_categories$Var1) {
  overlap_hepatocyte <- calculate_finemap_overlap(hepatocyte_enhancers[hepatocyte_enhancers$CellType == category,], category)
  hepatocyte_gwas = rbind(hepatocyte_gwas, overlap_hepatocyte)
}

hepatocyte_common_var <- data.frame()
for (category in hepatocyte_enhancer_categories$Var1) {
  temp <- hepatocyte_enhancers[hepatocyte_enhancers$CellType == category,]
  overlap_hepatocyte_common <- calculate_commonvar_overlap(temp)
  overlap_hepatocyte_common$CellType <- category 
  colnames(overlap_hepatocyte_common) <- c('number_of_snps', 'CellType')
  hepatocyte_common_var = rbind(hepatocyte_common_var, overlap_hepatocyte_common)
}
summary_hepatocyte <- hepatocyte_common_var %>% group_by(CellType) %>% summarise_at(vars(number_of_snps), list(name = mean))
summary_hepatocyte <- merge(summary_hepatocyte, hepatocyte_gwas, by = "CellType") 
summary_hepatocyte$OR <- summary_hepatocyte$number_of_snps / summary_hepatocyte$name

write_tsv(hepatocyte_common_var, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/hepatocyte_common_variant_overlap.tsv')
write_tsv(summary_hepatocyte, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/summary_hepatocyte_results.tsv')

# Bootstrap Miscellaneous
other_enhancers <- abc_enhancer_coords[abc_enhancer_coords$Category == 'Other',]
other_enhancer_categories <- data.frame(table(other_enhancers$CellType))

other_gwas = data.frame()
for (category in other_enhancer_categories$Var1) {
  overlap_other <- calculate_finemap_overlap(other_enhancers[other_enhancers$CellType == category,], category)
  other_gwas = rbind(other_gwas, overlap_other)
}

other_common_var <- data.frame()
for (category in other_enhancer_categories$Var1) {
  temp <- other_enhancers[other_enhancers$CellType == category,]
  overlap_other_common <- calculate_commonvar_overlap(temp)
  overlap_other_common$CellType <- category 
  colnames(overlap_other_common) <- c('number_of_snps', 'CellType')
  other_common_var = rbind(other_common_var, overlap_other_common)
}
summary_other <- other_common_var %>% group_by(CellType) %>% summarise_at(vars(number_of_snps), list(name = mean))
summary_other <- merge(summary_other, other_gwas, by = "CellType") 
summary_other$OR <- summary_other$number_of_snps / summary_other$name

write_tsv(other_common_var, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/other_common_variant_overlap.tsv')
write_tsv(summary_other, '/sc/arion/projects/ad-omics/ashvin/cell_type_enrichment_results/gwas/summary_other_results.tsv')



