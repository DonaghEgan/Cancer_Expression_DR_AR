# Title: Cancer-specific Differential Expression Analysis Pipeline (HMF data)
# Description: Analysis of RNA-seq data with tumor purity design formula
# Author: Donagh
# Date: 2025-01-22
# Version: 1.1

# Library
library(DESeq2); library(tidyverse); library(msigdbr); library(fgsea); library(pheatmap)
library(RColorBrewer); library(readxl); library(writexl); library(ggpubr);library(ggsci)
library(biomaRt);library(openxlsx)
source("src/custom_plots.R")
source("src/analysis_functions.R")

# Load counts + meta
counts <- readRDS("insert_path") ## available from the HMF. 
meta_data <- readRDS("data/meta/meta_data.Rds")

# Overlap counts and meta and remove non-assigned samples
meta_data <- meta_data[rownames(meta_data) %in% colnames(counts), ] 
meta_data <- meta_data[!meta_data$group == "NA",] 

# ConsensusTME output
consensus_tme <- data.frame(t(readRDS("data/processed/consensus_tme.Rds")))

# Unify meta and and patients with rna seq 
counts <- counts[, colnames(counts) %in% rownames(meta_data)]

# Define APC and NON_APC  
APC_CELLS <- c("Dendritic_cells" , "Monocytes","B_cells",
               "Macrophages",  "Macrophages_M1", "Macrophages_M2")

NON_APC_CELLS <- c("Cytotoxic_cells", "Eosinophils", "Mast_cells", "NK_cells","Neutrophils",
                   "T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta", "T_regulatory_cells",
                   "Plasma_cells")

# Average APC and NON_APC 
consensus_tme$APC_CELLS <- rowMeans(consensus_tme[, APC_CELLS, drop=FALSE]) 
consensus_tme$NON_APC_CELLS <- rowMeans(consensus_tme[, NON_APC_CELLS, drop = FALSE]) 

# Add tme to meta data 
meta_tme <- merge(consensus_tme, meta_data, by = 0) %>% column_to_rownames("Row.names")

# Scale variables in model
meta_tme_scale <- scale(meta_tme[,c(1:21,48), drop=F])
meta_tme_scale <- merge(meta_tme_scale, meta_tme[,"group", drop =F], by = 0) %>% column_to_rownames("Row.names")

# Run deseq analysis
result <- run_deseq_analysis(counts, meta_tme_scale, min_count = 10,
                             min_samples = 40, comparison = "groupNofailure.tumorPurity",
                             purity_var = "tumorPurity")

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
