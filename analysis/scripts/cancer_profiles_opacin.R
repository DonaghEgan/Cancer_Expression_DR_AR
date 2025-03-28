# Author: Donagh Egan
# Date: 23-01-25
# Description: Infer cancer cell-specific expression profiles
# from the Opacin-Neo dataset

# Load required libraries with error handling
required_packages <- c("GEOquery", "DESeq2", "tidyverse", "data.table",
                       "msigdbr", "ggpubr", "GSVA")

source("src/data_download.R")
source("src/analysis_functions.R")

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

# Load counts and meta
counts <- readRDS("counts.Rds")
meta_data <- read_csv("meta_data")
colnames(meta_data)

# available from source publication
meta_data <- read.csv("insert_path")
meta_data$TumorPerCent <- as.numeric(meta_data$TumorPerCent)
counts <- counts[,colnames(counts) %in% meta_data$sample] 

# Define APC and NON_APC  
APC_CELLS <- c("Dendritic_cells" , "Monocytes","B_cells",
               "Macrophages",  "Macrophages_M1", "Macrophages_M2")

NON_APC_CELLS <- c("Cytotoxic_cells", "Eosinophils", "Mast_cells", "NK_cells","Neutrophils",
                   "T_cells_CD4", "T_cells_CD8", "T_cells_gamma_delta", "T_regulatory_cells",
                   "Plasma_cells")

# Load TME
consensus_tme <- data.frame(t(readRDS("data/processed/consensus_tme_opacin.Rds")))

# Average APC and NON_APC 
consensus_tme$APC_CELLS <- rowMeans(consensus_tme[, APC_CELLS, drop=FALSE]) 
consensus_tme$NON_APC_CELLS <- rowMeans(consensus_tme[, NON_APC_CELLS, drop = FALSE]) 

# Merge meta and TME
meta_tme <- merge(consensus_tme, meta_data, by.x = 0, by.y = "sample") %>% column_to_rownames("Row.names")

# Scale
meta_tme_scale <- scale(meta_tme[,c(1:21,44), drop = F])
meta_tme_scale <- merge(meta_tme_scale, meta_tme[,"group_option3", drop =F], by =0) %>% column_to_rownames("Row.names")

# Run DESeq2
result <- run_deseq_analysis(counts, meta_tme_scale, min_count = 4,
                             min_samples = 49, purity_var = "TumorPerCent",
                             comparison = "group_option3PathResponse.TumorPerCent",
                             group_var = "group_option3")

# Create msigdb objects - Hallmarks
msigdb_H <- msigdbr(species = "human", category = "H") 
pathways_H <- split(x = msigdb_H$gene_symbol, f = msigdb_H$gs_name)

# Create msigdb objects - Reactome
msigdb_reactome <- msigdbr(species = "human", category = "C2")
msigdb_reactome <- msigdb_reactome[grepl(c("REACTOME"), msigdb_reactome$gs_name),] # reactome
pathways_R <- split(x = msigdb_reactome$gene_symbol, f = msigdb_reactome$gs_name)

# Pull Genes
genes <- result$results %>% pull(stat, gene) %>% sort(decreasing = T)

# FGSEA
fgsea_H <- fgsea::fgseaSimple(pathways = pathways_H,
                              stats = genes, nperm = 100000)

fgsea_R <- fgsea::fgseaSimple(pathways = pathways_R,
                              stats = genes, nperm = 100000)

# Plot
# ------------------------------------------------------------------------------

# Combine for plotting
combined_pathways <- rbind(fgsea_R, fgsea_H)

# Pathways of interest 
path_interest <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                   "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                   "HALLMARK_MYC_TARGETS_V1",
                   "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION")

combined_pathways <- combined_pathways[combined_pathways$pathway %in% 
                                         path_interest,]

combined_pathways$pathway <- gsub("REACTOME_", "", combined_pathways$pathway)
combined_pathways$pathway <- gsub("HALLMARK_", "", combined_pathways$pathway)

pdf("analysis/output/figures/enriched_pathways_opacin.pdf", width = 5.5, height = 3)
ggplot(combined_pathways, aes(x = NES, y = reorder(pathway, NES), fill = padj)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = "Normalized Enrichment Score (NES)", y = "Pathway") +
  ggtitle("Enriched Pathways: Response vs Non Response") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  viridis::scale_fill_viridis(discrete = FALSE, direction = -1, option = "D")
dev.off()
