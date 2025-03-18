# Title: CIN signature 
# Description: Score CIN genes from Bakjoum et al 
# on leading-edge genes and for AR,DR, and PR based on DEGs
# Author: Donagh
# Date: 2025-01-22
# Version: 1.1

library(readr)
library(tidyverse)

# load signature - available from supplementary figure of source paper
CIN_vector <- c("PELI2", "BMP2", "SHH", "TNS4", "RAB3B", "ROBO1", 
                 "ARHGAP28", "CHN2", "CST1", "F13A1", "CPVL", "SEMA6D", 
                 "NHSL2", "GTF2IP7", "DPYSL3", "PCDH7", "KHDRBS3", "TRAC", 
                 "TMEM156", "CST4", "CD24", "FGF5", "NTN4")

# as a list
cin_list <- list(CIN = CIN_vector)

# load ar vs pr genes  
AR_PR_hmf <- readRDS("analysis/output/results_data/hmf_deseq_purity_ar_pr_deg.Rds")
AR_PR_liu <- readRDS("analysis/output/results_data/liu_deseq_purity_ar_pr_deg.Rds")

# Pull and order genes by stat
genes_liu <- AR_PR_liu %>% pull(stat, gene) %>% sort(decreasing = T)
genes_hmf <- AR_PR_hmf %>% pull(stat, gene) %>% sort(decreasing = T)

# Pathway enrichment
fgsea_ar_pr_liu <- fgsea::fgseaSimple(pathways = cin_list, stats = genes_liu, nperm = 1000)
fgsea_ar_pr_hmf <- fgsea::fgseaSimple(pathways = cin_list, stats = genes_hmf, nperm = 1000)

fgsea_ar_pr_liu$pathway <- paste(fgsea_ar_pr_liu$pathway, "LIU", sep = "_")
fgsea_ar_pr_hmf$pathway <- paste(fgsea_ar_pr_hmf$pathway, "HMF", sep = "_")

combined_cin <- rbind(fgsea_ar_pr_liu, fgsea_ar_pr_hmf)
saveRDS(combined_cin, "analysis/output/results_datacin_fgsea_pr_ar.Rds")

pdf("analysis/output/figures/cin_sig_hmf_liu.pdf", width = 3, height = 2)
ggplot(combined_cin, aes(x = NES, y = reorder(pathway, NES), color = padj)) +
  geom_point(size = 5) +
  scale_color_viridis_c() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(x = "Normalized Enrichment Score (NES)",
       y = "Pathway",
       color = "padj") +
  theme_bw(base_size = 9) +
  xlim(-2,1)
dev.off()