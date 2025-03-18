# Title: Spatial Signatures 
# Description: Create signatures for pathways based 
# on leading-edge genes and for AR,DR, and PR based on DEGs
# Author: Donagh
# Date: 2025-01-22
# Version: 1.1

source("src/analysis_functions.R")

# Directory containing your fgsea
my_dir <- "data/processed/"

# List all files that contain "fgsea" in their names
files_fgsea <- list.files(
  path        = my_dir, 
  pattern     = "fgsea", 
  full.names  = TRUE
)

# Get names
file_names <- basename(files_fgsea)

# Read them into a list
my_data_list <- lapply(files_fgsea, function(x) {
  readRDS(x)
})

names(my_data_list) <- file_names

# Extract and filter genes for each pathway
IFNG <- find_common_genes(list(
  my_data_list[["fgsea_hmf_H_cancer.Rds"]]$leadingEdge[[27]],
  my_data_list[["fgsea_liu_H_cancer.Rds"]]$leadingEdge[[27]]
))

IFNA <- find_common_genes(list(
  my_data_list[["fgsea_hmf_H_cancer.Rds"]]$leadingEdge[[26]],
  my_data_list[["fgsea_liu_H_cancer.Rds"]]$leadingEdge[[26]]
))

MYC <- find_common_genes(list(
  my_data_list[["fgsea_hmf_H_cancer.Rds"]]$leadingEdge[[32]],
  my_data_list[["fgsea_liu_H_cancer.Rds"]]$leadingEdge[[32]]
))

MHCII <- find_common_genes(list(
  my_data_list[["fgsea_hmf_R_cancer.Rds"]]$leadingEdge[[781]],
  my_data_list[["fgsea_liu_R_cancer.Rds"]]$leadingEdge[[783]]
))

# Print the results
print("Genes in IFNG appearing at least twice:")
print(IFNG)

print("Genes in IFNA appearing at least twice:")
print(IFNA)

print("Genes in MYC appearing at least twice:")
print(MYC)

print("Genes in MHCII appearing at least twice:")
print(MHCII)

# Create a named list of pathways and genes
pathway_list <- list(
  IFNG = IFNG,
  IFNA = IFNA,
  MYC = MYC,
  MHCII = MHCII
)

# Convert to a dataframe
gene_pathway_df <- do.call(rbind, lapply(names(pathway_list), function(pathway) {
  data.frame(
    Pathway = pathway,
    Gene = pathway_list[[pathway]]
  )
}))

write_csv(gene_pathway_df, "data/processed/pathway_genes_2datasets.csv")

#-------------------------------------------------------------------------------
# AR _ DR _ Signatures
#-------------------------------------------------------------------------------

DR_AR_hmf <- readRDS("analysis/output/results_data/hmf_deseq_purity_nf_aq_deg.Rds")
AR_PR_hmf <- readRDS("analysis/output/results_data/hmf_deseq_purity_ar_pr_deg.Rds")
DR_PR_hmf <- readRDS("analysis/output/results_data/hmf_deseq_purity_nf_pr_deg.Rds")

DR_AR_liu <- readRDS("analysis/output/results_data/liu_deseq_purity_nf_aq_deg.Rds")
AR_PR_liu <- readRDS("analysis/output/results_data/liu_deseq_purity_ar_pr_deg.Rds")
DR_PR_liu <- readRDS("analysis/output/results_data/liu_deseq_purity_nf_pr_deg.Rds")

# HMF dataset comparisons
DR_AR_hmf_sig <- DR_AR_hmf[DR_AR_hmf$padj < 0.1 & DR_AR_hmf$log2FoldChange > 0, ]$gene
AR_DR_hmf_sig <- DR_AR_hmf[DR_AR_hmf$padj < 0.1 & DR_AR_hmf$log2FoldChange < 0, ]$gene

AR_PR_hmf_sig <- AR_PR_hmf[AR_PR_hmf$padj < 0.1 & AR_PR_hmf$log2FoldChange > 0, ]$gene
PR_AR_hmf_sig <- AR_PR_hmf[AR_PR_hmf$padj < 0.1 & AR_PR_hmf$log2FoldChange < 0, ]$gene

DR_PR_hmf_sig <- DR_PR_hmf[DR_PR_hmf$padj < 0.1 & DR_PR_hmf$log2FoldChange > 0, ]$gene
PR_DR_hmf_sig <- DR_PR_hmf[DR_PR_hmf$padj < 0.1 & DR_PR_hmf$log2FoldChange < 0, ]$gene

# Liu dataset comparisons
DR_AR_liu_sig <- DR_AR_liu[DR_AR_liu$padj < 0.1 & DR_AR_liu$log2FoldChange > 0, ]$gene
AR_DR_liu_sig <- DR_AR_liu[DR_AR_liu$padj < 0.1 & DR_AR_liu$log2FoldChange < 0, ]$gene

AR_PR_liu_sig <- AR_PR_liu[AR_PR_liu$padj < 0.1 & AR_PR_liu$log2FoldChange > 0, ]$gene
PR_AR_liu_sig <- AR_PR_liu[AR_PR_liu$padj < 0.1 & AR_PR_liu$log2FoldChange < 0, ]$gene

DR_PR_liu_sig <- DR_PR_liu[DR_PR_liu$padj < 0.1 & DR_PR_liu$log2FoldChange > 0, ]$gene
PR_DR_liu_sig <- DR_PR_liu[DR_PR_liu$padj < 0.1 & DR_PR_liu$log2FoldChange < 0, ]$gene

# Gene unions for forward comparisons
DR_AR_union <- union(DR_AR_hmf_sig, DR_AR_liu_sig)
AR_DR_union <- union(AR_DR_hmf_sig, AR_DR_liu_sig)

AR_PR_union <- union(AR_PR_hmf_sig, AR_PR_liu_sig)
PR_AR_union <- union(PR_AR_hmf_sig, PR_AR_liu_sig)

DR_PR_union <- union(DR_PR_hmf_sig, DR_PR_liu_sig)
PR_DR_union <- union(PR_DR_hmf_sig, PR_DR_liu_sig)

# Create lists for each comparison
pathway_lists <- list(
  DR_AR = DR_AR_union,
  AR_DR = AR_DR_union,
  AR_PR = AR_PR_union,
  PR_AR = PR_AR_union,
  DR_PR = DR_PR_union,
  PR_DR = PR_DR_union
)

ar_dr_pathway_df <- do.call(rbind, lapply(names(pathway_lists), function(pathway) {
  data.frame(
    Pathway = pathway,
    Gene = pathway_lists[[pathway]]
  )
}))

# save as CSV
write_csv(ar_dr_pathway_df, "data/processed/group_comparisons.csv")
