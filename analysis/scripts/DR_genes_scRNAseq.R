# Title: Bulk-to-single cell validation
# Description: Analyze differentially expressed genes between DRvsAR patients in scRNAseq
# Author: Donagh Egan
# Date: 2025-01-22
# Version: 1.1

# Library
library(Seurat);library(limma);library(fgsea);library(tidyverse);library(ggsci)
library(ggpubr);library(msigdbr);library(ggpubr);library(readxl);library(sjPlot)
library(sjmisc);library(interactions)
source("src/data_download.R")
source("src/analysis_functions.R")

# Load differentially expressed genes
deg_hmf <- readRDS("data/processed/hmf_deseq_purity_nf_aq_deg.Rds")
hmf_up <- deg_hmf[deg_hmf$padj < 0.1 & deg_hmf$log2FoldChange > 0,] # DR
hmf_down <- deg_hmf[deg_hmf$padj < 0.1 & deg_hmf$log2FoldChange < 0,] # AR

# Immune system genes -> reactome
immune_genes <- qusage::read.gmt("data/processed/REACTOME_IMMUNE_SYSTEM.v2024.1.Hs.gmt")

# Download data from Arnon et al
data_path <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978%5Fcounts.csv.gz"

# Clean - removes after download 
sc_data <- download_geo_data(ftp_url = data_path, rownames_col = "X", cleanup = TRUE)

# Load cell and patient info 
cell_annot <- read_csv("data/meta/GSE115978_cell.annotations.csv")
patient_info <- readxl::read_xlsx("data/meta/GSE115978_patient_info.csv", skip = 2)

# Join patient and cell data 
processed_annot <- merge(patient_info, cell_annot, by.x = "Sample", by.y= "samples")

# Filter data
processed_annot <- processed_annot %>%
  filter(
    `Treatment group` == "Untreated",  # Using backticks instead of [[ ]]
    cell.types != "?"
  ) %>% column_to_rownames("cells")

# Unify meta and counts
sc_filter <- sc_data[,colnames(sc_data) %in% rownames(processed_annot)]

# Process seurat object
seurat_obj <- seurat_wrapper(counts = sc_filter,
                             metadata = processed_annot,
                             project_name = "hmf", min_cells = 10, 
                             min_features = 100,scale_factor = 10000,
                             n_var_features = 2000, scale_all_features = TRUE)

# Values to subset
sample_id <- "Sample"  
celltype1 <- "Mal"  
celltype2 <- c("T.cell", "T.CD8", "T.CD4")  

# Get common samples 
common_samples <- intersect(
  unique(seurat_obj@meta.data[seurat_obj@meta.data$cell.types == celltype1, sample_id]),
  unique(seurat_obj@meta.data[seurat_obj@meta.data$cell.types %in% celltype2, sample_id])
)

# Single subset operation with pre-filtered samples
subsets <- list(
  celltype1 = subset(seurat_obj, subset = cell.types == celltype1 & Sample %in% common_samples),
  celltype2 = subset(seurat_obj, subset = cell.types %in% celltype2 & Sample %in% common_samples)
)

# Calculate average expression for both cell types at once
exprs_list <- lapply(subsets, function(x) {
  AverageExpression(x, group.by = sample_id, return.seurat = TRUE, assays = "RNA")
})

# Extract expression matrices with aligned samples
exprs_mal <- exprs_list$celltype1[, common_samples]
exprs_tcell <- exprs_list$celltype2[, common_samples]

# Score HMF genes in samples of cancer cells only
cancer_bulk_score <- AddModuleScore(
  object = exprs_mal,
  features = list(score = c(hmf_up$gene)),
  ctrl = 100,  
  name = "ModuleScore"  
)

# Extract t cell expression data
tcell_data <- data.frame(t(GetAssayData(exprs_tcell, slot = "data")))
tcell_data <- tcell_data[,colnames(tcell_data) %in% immune_genes$REACTOME_IMMUNE_SYSTEM, drop = F]
tcell_data <- tcell_data[,colSums(tcell_data) > 10 & apply(tcell_data, 2, function(x) sum(x == 0)) < nrow(tcell_data)]

# Extract module scores from seurat_object1
hmf_scores <- cancer_bulk_score@meta.data[ , grep("ModuleScore", colnames(cancer_bulk_score@meta.data)), drop = F]

# Merge t cell expression data & HMF scores
hmf_tcell_data <- merge(hmf_scores, tcell_data, by = 0) %>% column_to_rownames("Row.names")

# module score to be correlated
target_column <- "ModuleScore1"

# Initialize a data frame to store the results
results <- data.frame(column = character(),
                      estimate = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)

# Loop through each column in the data frame
for (col in colnames(hmf_tcell_data)) {
  
  if (col != target_column) {
    
    # Perform correlation test
    test_result <- cor.test(hmf_tcell_data[[target_column]], hmf_tcell_data[[col]], 
                            use = "complete.obs")
    
    # Append the results to the data frame
    results <- rbind(results, data.frame(column = col,
                                         estimate = test_result$estimate,
                                         p_value = test_result$p.value))
  }
}

results$FDR <- p.adjust(results$p_value, method = "BH")
plot_results <- results %>% arrange(estimate) %>% dplyr::slice(c(1:3, (n()-2):n()))
plot_results$label <- ifelse(plot_results$FDR < 0.05, "0.05", "non_sig")

#-------------------------------------------------------------------------------
# Analyze HMF genes in cancer cells and the relationship with MHCII:IFNG
#-------------------------------------------------------------------------------

# load mitf program
mitf_program <- read_xlsx("data/processed/tirosh_et_al_mitf.xlsx", col_names = F)

# create pathways - reactome and hallmarks
msigdbr_h <- msigdbr(species = "human", category = "H")
msigdbr_r <- msigdbr(species = "human", category = "C2")

pathways_R = split(x = msigdbr_r$gene_symbol, f = msigdbr_r$gs_name)
pathways_H = split(x = msigdbr_h$gene_symbol, f = msigdbr_h$gs_name)

# Calculate module scores
cancer_cell_score <- AddModuleScore(
  object = subsets$celltype1, # cancer cells 
  features = list(mhc = c(pathways_R[["REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION"]]),
                  ifng = c(pathways_H[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]),
                  dr_sig = c(hmf_up$gene),
                  mitf = c(mitf_program$...1),
                  axl = c(mitf_program$...2)),
  ctrl = 300,  
  name = "path_Score"  
)

# Calculate ratios using residuals
cancer_cell_score@meta.data$ratio_mhc_ifng <- residuals(lm(path_Score1 ~ path_Score2,
                                                           cancer_cell_score@meta.data))

cancer_cell_score@meta.data$ratio_mitf_axl <- residuals(lm(path_Score4 ~ path_Score5,
                                                           cancer_cell_score@meta.data))

# fit linear model - predict ratio mhc:ifng with mitf as predictor and dr sig as
# interaction term
lm_cancer <- lm(ratio_mhc_ifng ~ ratio_mitf_axl*path_Score3, data = cancer_cell_score@meta.data)
summary(lm_cancer)


