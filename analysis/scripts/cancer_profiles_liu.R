# Library
library(DESeq2); library(tidyverse); library(msigdbr); library(fgsea); library(pheatmap)
library(RColorBrewer); library(readxl); library(writexl); library(ggpubr);library(ggsci)
library(biomaRt);library(openxlsx)
source("src/custom_plots.R")
source("src/analysis_functions.R")

# Title: Cancer-specific Differential Expression Analysis Pipeline (Liu data)
# Description: Analysis of RNA-seq data with tumor purity design formula
# Author: Donagh
# Date: 2025-01-22
# Version: 1.1

# Loud counts - also available from Liu et al 2019 Github
counts <- read_delim("data/processed/rnaseq_rawcounts_liu.txt")
counts <- counts %>% column_to_rownames("...1")

# Load meta data and deconvolution
meta_data <- readRDS("data/meta/meta_data_liu.Rds")
TME_liu <- data.frame(t(readRDS("data/processed/consensus_tme_liu.Rds")))

# Subset counts according to meta data
counts <- counts[,colnames(counts) %in% meta_data$PATIENT_ID]

# Define APC and NON_APC  
APC_CELLS <- c("Dendritic_cells" , "Monocytes","B_cells",
               "Macrophages",  "Macrophages_M1", "Macrophages_M2")

NON_APC_CELLS <- c("Cytotoxic_cells", "Eosinophils", "Mast_cells", "NK_cells","Neutrophils",
                   "T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta", "T_regulatory_cells",
                   "Plasma_cells")

# Average APC and NON_APC 
TME_liu$APC_CELLS <- rowMeans(TME_liu[, APC_CELLS, drop=FALSE]) 
TME_liu$NON_APC_CELLS <- rowMeans(TME_liu[, NON_APC_CELLS, drop = FALSE]) 

# Merge meta and tme
meta_tme <- merge(TME_liu, meta_data, by.x = 0, by.y="SAMPLE_ID") %>% column_to_rownames("PATIENT_ID")
meta_tme$group <- as.factor(meta_tme$group) # group as factor 

# Scale model covariates 
meta_tme_scale <- scale(meta_tme[,c(2:22,32), drop=F])
meta_tme_scale <- merge(meta_tme_scale, meta_tme[,"group", drop =F], by =0) %>% column_to_rownames("Row.names")

## unify order of patients 
patient_order <- rownames(meta_tme)
counts_order <- counts[, patient_order]

# Run deseq analysis
result <- run_deseq_analysis(counts_order, meta_tme_scale, min_count = 5,
                             min_samples = 15, comparison = "groupNofailure.PURITY",
                             purity_var = "PURITY")

# Create msigdb objects - Reactome
msigdb_reactome <- msigdbr(species = "human", category = "C2")
msigdb_reactome <- msigdb_reactome[grepl(c("REACTOME"), msigdb_reactome$gs_name),] # reactome

# Create msigdb objects - Hallmarks
msigdb_hallmarks <- msigdbr(species = "human", category = "H") # hallmarks

# Create pathways
pathways_R <- split(x = msigdb_reactome$gene_symbol, f = msigdb_reactome$gs_name)
pathways_H <- split(x = msigdb_hallmarks$gene_symbol, f = msigdb_hallmarks$gs_name)

# Pull and order genes by stat
genes <- result$results %>% pull(stat, gene) %>% sort(decreasing = T)

# Pathway enrichment
fgsea_R_nf_ar <- fgsea::fgseaSimple(pathways = pathways_R, stats = genes, nperm = 100000)
fgsea_H_nf_ar <- fgsea::fgseaSimple(pathways = pathways_H, stats = genes, nperm = 100000)
