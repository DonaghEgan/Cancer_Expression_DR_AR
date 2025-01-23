# Title: Cell line differentiation analysis. https://doi.org/10.1016/j.ccell.2018.03.017
# Author: Donagh Egan
# Date: 23-01-25
# Description: Analysis of RNA-Seq data from  Tsoi et al

# Load required libraries with error handling
required_packages <- c("GEOquery", "DESeq2", "tidyverse", "data.table",
                       "msigdbr", "ggpubr", "GSVA", "ggsci", "readxl",
                       "patchwork")

source("src/data_download.R")
source("src/custom_plots.R")

# Function to install and load packages
load_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing package:", pkg))
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

# Load all required packages
tryCatch({
  load_packages(required_packages)
}, error = function(e) {
  stop("Failed to load required packages. Please check your R environment.")
})

## load degs
deg_hmf <- readRDS("Z:/working/GEL_MELANOMA/External_Datasets/HMF_GCP/RNA_seq_processed/OUTPUT/hmf_deseq_purity_nf_aq_deg.Rds")
hmf_up <- deg_hmf[deg_hmf$padj < 0.1 & deg_hmf$log2FoldChange > 0,] # AR
hmf_down <- deg_hmf[deg_hmf$padj < 0.1 & deg_hmf$log2FoldChange < 0,] #DR

# Set seed for reproducibility
set.seed(123)

# Read the gzipped text file
data <- read.table(gzfile("data/processed/GSE80829_Melanoma_cell_lines_RNASeq_FPKM.txt.gz"),
                   header = TRUE, sep = "\t") %>% column_to_rownames("Gene")
colnames(data) <- gsub("_DMSO", "", colnames(data))

# Meta
meta <- read_xlsx("data/meta/dediff_meta.xlsx", skip = 1) %>%
  rename(cell_line = 1)

# Data preprocessing
data <- data[, colnames(data) %in% meta$cell_line]
meta <- meta[meta$cell_line %in% colnames(data), ]
data <- round(data)

# create deseq object
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = meta,
                              design = ~ Subtype)
# filtering 
keep <- rowSums(counts(dds) >= 1) >= 10 
table(keep)
dds <- dds[keep,]

# run comparison  
dds <- DESeq(dds)
res <- results(dds)

# normalize
vsd <- vst(dds, blind=FALSE)
normalised_data <- data.frame(assay(vsd))

# merge
exp_meta <- merge(t(normalised_data), meta, by.x = 0, by.y = "cell_line")
exp_meta$Subtype <- factor(exp_meta$Subtype, levels = c("Undifferentiated", "Neural crest like",
                                                        "Transitory", "Melanocytic"))

## score 
pathways <- list("dr" = hmf_up$gene, "ar" = hmf_down$gene)
path_score <- GSVA::gsva(as.matrix(normalised_data), gset.idx.list = pathways, method = "zscore")

## merge
score_meta <- merge(t(path_score), meta, by.x = 0, by.y = "cell_line")
score_meta$Subtype <- factor(score_meta$Subtype, levels = c("Undifferentiated", "Neural crest like",
                                                            "Transitory", "Melanocytic"))

# Define comparisons for each plot
comparisons_mitf <- list(c("Transitory", "Undifferentiated"),
                         c("Transitory", "Neural crest like"),
                         c("Transitory", "Melanocytic"), 
                         c("Neural crest like", "Melanocytic"))

comparisons_ar <- list(c("Neural crest like", "Undifferentiated"),
                       c("Neural crest like", "Melanocytic"))

# x axis labels
labels <- c("U", "N", "T", "M")

# Generate the individual plots using the function
plot1 <- plot_boxplot(exp_meta, "MITF", "MITF Expression",
                      comparisons_mitf, x_labels = labels)
plot2 <- plot_boxplot(exp_meta, "MYC", "MYC Expression", comparisons_mitf,
                      x_labels = labels)
plot3 <- plot_boxplot(score_meta, "dr", "DR Signature", comparisons_mitf,
                      x_labels = labels)
plot4 <- plot_boxplot(score_meta, "ar", "AR Signature", comparisons_ar,
                      x_labels = labels)

# Combine the four plots into a 2x2 grid with a shared legend
combined_plot <- (plot1 | plot2) / (plot3 | plot4) +
  plot_layout(guides = "collect")

# Save the combined plot to a PDF file
pdf("analysis/output/figures/myc_dediff_cell_lines.pdf", height = 4.5, width = 5)
combined_plot
dev.off()

