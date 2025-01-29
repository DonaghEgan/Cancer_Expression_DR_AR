# Title: Tumor Mutation Burden (TMB) Analysis and Visualization
# Description: This script reads TMB data and compares log10(TMB) across different outcome groups.
# Author: Donagh
# Date: 2025-01-28
# Version: 1.1

# Load Required Libraries
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)
library(readr)
source("src/custom_plots.R")

# Define relative paths for reproducibility
tmb_file <- "data/processed/totMuts_sampleIDs.txt"
meta_data_file <- "data/meta/meta_data.Rds"  
counts_file <- "insert_path" ## available from the HMF. 

# Read Tumor Mutation Burden (TMB) data
tmb <- read_delim(tmb_file, delim = "\t", col_types = cols())

# Read Metadata
meta_data <- readRDS(meta_data_file)
counts <- readRDS(counts_file)
meta_data <- meta_data[rownames(meta_data) %in% colnames(counts), ]## load data 

# Define Color Palette   
pal_colors <- rev(ggsci::pal_npg("nrc")(3)[c(2, 1, 3)])  

# Data Processing
tmb_groups <- meta_data %>%
  rownames_to_column(var = "SampleID") %>%  # Add rownames as a column
  filter(group != "NA") %>%  # Filter out rows with 'NA' in 'group'
  dplyr::inner_join(tmb, by = "SampleID") %>%  # Join on 'SampleID'
  # Calculate log10(TMB) and set factor levels for 'group'
  mutate(
    log_TMB = log10(TMB),
    group = factor(group, levels = c("Primary_Resistance", "Acq_Resistance", "Nofailure"))
  )

# Define pairwise comparisons for statistical annotations
comparisons <- list(
  c("Acq_Resistance", "Nofailure"),
  c("Acq_Resistance", "Primary_Resistance"),
  c("Nofailure", "Primary_Resistance")
)

#Generate Plot
tmb_plot <- plot_boxplot(data = tmb_groups, x_var = "group", y_var = "log_TMB", 
                         y_label = "Log(10) TMB", comparisons = comparisons,
                         x_labels =  c("Primary Res", "Acquired Res", "Durable Rsp"),
                         color_pal = pal_colors)

pdf("analysis/output/figures/tmb_hmf_groups.pdf", height = 3, width = 3)
tmb_plot
dev.off()
