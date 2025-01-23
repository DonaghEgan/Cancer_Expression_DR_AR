# Title: IFNG stimulation in Melanoma cell lines
# Author: Donagh Egan
# Date: 23-01-25
# Description: Analysis of RNA-Seq data from GEO dataset GSE152755

# Load required libraries with error handling
required_packages <- c("GEOquery", "DESeq2", "tidyverse", "data.table",
                       "msigdbr", "ggpubr", "GSVA")

source("src/data_download.R")

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

# Set seed for reproducibility
set.seed(123)

# Define global variables
GSE_ID <- "GSE152755"
TAR_FILE_PATH <- file.path("data", "processed", "GSE152755_RAW.tar")

# Get files from GEO
gse <- download_geo_direct(gse_id = GSE_ID)
gse_data <- gse[[1]]

# Extract tar files
raw_files <- extract_tar_files(TAR_FILE_PATH)

# Process sample files
data_list <- process_sample_files(raw_files$files)

# Merge count data
merged_counts <- Reduce(function(...) merge(..., all = TRUE), data_list) %>%
  column_to_rownames("GeneId")

# Prepare experimental design
exp_design <- data.frame(
  sample = rownames(gse_data@phenoData@data),
  condition = gse_data@phenoData@data$title,
  cell_line = gse_data@phenoData@data$`cell line:ch1`
) %>%
  mutate(condition = sub("^[^_]*_", "", condition))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = round(merged_counts),
  colData = exp_design,
  design = ~ condition
)

# Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

#Normalize data
vsd <- vst(dds)
normalized_data <- data.frame(assay(vsd))

# Final data preparation
exp_data <- merge(exp_design, t(normalized_data), by.x = "sample", by.y = 0) %>%
  filter(condition != "TNF_RNA")

# Create pathways
msigdb <- msigdbr(species = "human", category = "H") # hallmarks
pathways <- split(x = msigdb$gene_symbol, f = msigdb$gs_name)

# Perform GSVA analysis
gsva_results <- GSVA::gsva(
  as.matrix(normalized_data), 
  gset.idx.list = pathways, 
  method = "zscore"
)

# Transpose and merge with experimental design
score_df <- merge(
  exp_design, 
  t(gsva_results),  
  by.x = "sample", 
  by.y = 0
) %>%
  filter(condition != "TNF_RNA")

# plot
pdf("analysis/output/figures/myc_ifng_trx.pdf", height = 3, width = 3)
ggplot(score_df, aes(x=condition, y = HALLMARK_MYC_TARGETS_V1, fill=condition)) + geom_boxplot(alpha = 0.9) +
  theme(axis.text.x = element_text(angle = 90)) +
  stat_compare_means(paired = T) + geom_point(alpha = 0.5) + theme_classic() +
  scale_fill_manual(values = c("#D94F4F", "#6BAF92")) +
  theme(legend.position = "none") + xlab("") +
  scale_x_discrete(labels = c("D0", "IFNG")) +
  ylab("MYC Activity")
dev.off()
