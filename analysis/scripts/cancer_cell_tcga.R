# Title: Cancer cell-specific analysis in TCGA
# Description: Exploring the AR and DR signatures
# in cancer cell expression profiles from TCGA,
# inferred using CIBSERSORTx
# Author: Donagh
# Date: 2025-01-22
# Version: 1.1

# Library
library(GSVA)
library(data.table)
library(ggpubr)
library(tidyverse)
library(survminer)
library(survMisc)
library(ggsci)
library(survival)
library(matrixStats)
library(car)
source("src/custom_plots.R")

# Global variables.
TOP_N_GENES <- 500  # Number of top genes for each signature slice
TMB_CUTOFF  <- 10   # TMB threshold for High vs. Low

# Output file paths
PATH_MALIG_EXP <- "data/processed/CIBERSORTxHiRes_Job44_Malignant_Window32._deseq_tcga.txt"
PATH_DEG_HMF <- "data/processed/hmf_deseq_purity_nf_aq_deg.Rds"
PATH_SAMPLE_META <- "data/processed/data_clinical_sample.txt"
PATH_PATIENT_META <- "data/processed/data_clinical_patient.txt"
PATH_CLONE_COUNTS <- "data/processed/wolf_et_al_clones.txt"

# Read data
malig_exp <- fread(PATH_MALIG_EXP) %>%
  as.data.frame() %>%
  column_to_rownames("GeneSymbol")

# Convert '.' to '-' in column names
colnames(malig_exp) <- gsub("\\.", "-", colnames(malig_exp))

# Log2 transform
malig_exp <- log2(malig_exp + 1)

# Read in DEGs (HMF) from RDS
deg_hmf <- readRDS(PATH_DEG_HMF)

# Read sample metadata - TMB
sample_meta <- read.delim(PATH_SAMPLE_META, skip = 4)

# Read patient metadata - survival info
patient_meta <- read.delim(PATH_PATIENT_META, header = TRUE, comment.char = "#")

# Read in clone count estimates
clone_counts <- read.delim(PATH_CLONE_COUNTS)

# Filter and prepare genes
hmf_up <- deg_hmf %>%
  arrange(desc(stat)) %>%
  slice(seq_len(TOP_N_GENES))

hmf_down <- deg_hmf %>%
  arrange(stat) %>%  # ascending order
  slice(seq_len(TOP_N_GENES))

# Subset malignant expression by selected gene lists
malig_hmf_up <- malig_exp[rownames(malig_exp) %in% hmf_up$gene, ]
malig_hmf_down <- malig_exp[rownames(malig_exp) %in% hmf_down$gene, ]

# Combined Signatures
combined_sig <- data.frame(
  hmf_up   = colMedians(t(scale(t(malig_hmf_up))),   na.rm = TRUE),
  hmf_down = colMedians(t(scale(t(malig_hmf_down))), na.rm = TRUE),
  SAMPLE_ID = colnames(malig_hmf_down)
)

# Merge with sample meta data - contains TMB estimates
sample_meta_sig <- merge(combined_sig, sample_meta, by = "SAMPLE_ID", all.x = TRUE)

# ITH analysis
# format for merging
#clone_counts <- clone_counts %>%
#  mutate(Patient = gsub("-", ".", Patient))
clone_counts_sig <- merge(
  clone_counts,
  sample_meta_sig,
  by.x = "Patient",
  by.y = "PATIENT_ID",
  all.x = F  # Keep all intra_hetero rows if you want
)

# Remove Samples with no TMB
sample_meta_sig <- sample_meta_sig[!is.na(sample_meta_sig$TMB_NONSYNONYMOUS),]

# Group by TMB
sample_meta_sig$tmb_group <- ifelse(sample_meta_sig$TMB_NONSYNONYMOUS < 10, "Low\n(<10 mut/Mb)", "High\n(>10 mut/Mb)")

# make clones numeric
clone_counts_sig$clones <- as.numeric(clone_counts_sig$clones)

# calculate cor
hmf_cor <- cor.test(
  x = clone_counts_sig$hmf_up,
  y = clone_counts_sig$clones,
  method = "spearman"
)

# Extract correlation coefficient and p-value
rho     <- hmf_cor$estimate
p_value <- hmf_cor$p.value

# PLOT
p1 <- ggscatter(combined_sig, x = "hmf_up", y = "hmf_down",
          add = "reg.line",
          conf.int = TRUE,
          cor.coef = TRUE,
          cor.method = "spearman",
          color = pal_npg()(4)[4],
          add.params = list(color = "red",
                            fill = "lightgray"),
          size = 0.9,
          alpha = 0.6,
          title = "TCGA",
          xlab = "AR Expression",
          ylab = "DR Expression")

# Get comparisons 
comp <- combn(levels(as.factor(sample_meta_sig$tmb_group)), 2, simplify = FALSE)

# order 
p2 <- plot_boxplot(sample_meta_sig, x_var = "tmb_group", y_var = "hmf_up", alpha = 0.8,
  y_label = "Exp of DR signature", x_label = "TMB", comparisons = comp, color_pal = rev(pal_frontiers()(10)[c(3,1)]))

clone_counts_sig$clones <- as.factor(clone_counts_sig$clones)
p3 <- ggplot(clone_counts_sig, aes(x = clones, y = hmf_up, color = clones)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, na.rm = T) +
  theme_classic() +
  ylab("Mean Exp NF signature") +
  theme(legend.position = "none") +
  stat_summary(
    geom = "crossbar",
    fun.y = "median",
    col = "black",
    size = 0.3, width = 0.35, na.rm = T
  ) +
  scale_color_npg() +
  annotate("text", x = Inf, y = Inf, label = paste("Spearman's rho =", round(rho, 2),
  "\np-value =", format.pval(p_value, digits = 2)), 
  hjust = 1.7, vjust = 1.1, size = 4, color = "black")

pdf("analysis/output/figures/tcga_cancer_analysis.pdf", onefile = TRUE, height = 3, width = 3)
p1
p2
p3
dev.off()
