# Title: Baegev Validation
# Description: This script analyzes gene expression data from purified populations to validate
# differential expression results from univariate and cancer-cell specific models.
# Author: Donagh
# Date: 2025-05-22
# Version: 1.1

# Load required libraries
library(tidyverse)   
library(ggpubr)     
library(readxl)      
library(reshape2)    

# Define input and output directories (relative paths for portability)
input_dir <- "data/processed/"
output_dir <- "output/figures/"

# Load purified population data
baegev_sig <- read_xlsx(file.path(input_dir, "bagaev_cell_exp.xlsx"), skip = 1)

# Remove rows with missing values and set 'Gene' as row names
baegev_sig <- na.omit(baegev_sig) %>% column_to_rownames("Gene")

# Load differential expression results
purity_model_deseq <- readRDS(file.path(output_dir, "hmf_deseq_purity_res_nonres_deg.Rds"))
standard_model_deseq <- readRDS(file.path(output_dir, "hmf_deseq_standard_res_nonres_deg.Rds"))

# Select significant genes (adjusted p-value < 0.1)
pur_sig <- purity_model_deseq[purity_model_deseq$padj < 0.1, ]
std_sig <- standard_model_deseq[standard_model_deseq$padj < 0.1, ]

# Subset purified population data for significant genes
baegev_std <- baegev_sig[rownames(baegev_sig) %in% rownames(std_sig), ]
baegev_pur <- baegev_sig[rownames(baegev_sig) %in% rownames(pur_sig), ]

# Scale and center the data (Z-score normalization)
# This standardizes gene expression values for comparison across cell types
baegev_std <- data.frame(t(scale(t(baegev_std))))
baegev_pur <- data.frame(t(scale(t(baegev_pur))))

# Melt the data for plotting using tidyr::pivot_longer (modern alternative to reshape2::melt)
baegev_std_melt <- baegev_std %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "CellType", values_to = "Expression") %>%
  filter(!is.na(Expression))  # Remove any NA or NaN values

baegev_pur_melt <- baegev_pur %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "CellType", values_to = "Expression") %>%
  filter(!is.na(Expression))  # Remove any NA or NaN values

# Create boxplot for Univariate Model
pdf(file.path(output_dir, "baegev_non_model_val.pdf"), height = 4, width = 4)
ggplot(baegev_std_melt, aes(x = reorder(CellType, Expression, FUN = median), y = Expression)) +
  geom_boxplot(fill = "#919BC6") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.09)) +
  ylab("Gene Expression (Z score)") +
  xlab("Cell Type") +
  stat_compare_means(label = "p.format", label.x.npc = 0.05) +  # Default: Wilcoxon test
  ggtitle("Univariate Model")
dev.off()

# Create boxplot for Cancer-cell Specific Model
pdf(file.path(output_dir, "baegev_model_val.pdf"), height = 4, width = 4)
ggplot(baegev_pur_melt, aes(x = reorder(CellType, Expression, FUN = median), y = Expression)) +
  geom_boxplot(fill = "#B6E76A") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.09)) +
  ylab("Gene Expression (Z score)") +
  xlab("Cell Type") +
  stat_compare_means(label = "p.format", label.x.npc = 0.05) +  # Default: Wilcoxon test
  ggtitle("Cancer-cell Specific Model")
dev.off()