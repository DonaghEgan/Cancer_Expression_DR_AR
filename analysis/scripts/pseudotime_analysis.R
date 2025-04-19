library(Seurat)
library(patchwork)
library(readr)
library(readxl)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)
library(interactions)
library(viridis)
library(monocle3)
library(AUCell)
library(SingleCellExperiment)
library(SeuratWrappers)
library(RColorBrewer)
library(harmony)
library(monocle3)
library(AUCell)
library(msigdbr)

pal_custom <- c(pal_nejm()(8), pal_npg()(10))

# load seurat object
seurat_obj <- readRDS("Malignant_cells.rds") # https://doi.org/10.48804/GSAXBN

# select pre treatment only
seurat_obj <- subset(seurat_obj, subset = Timepoint == "BT")

# file paths
base_path <- "setwd"
signatures_file <- file.path(base_path, "pathway_genes_2datasets.csv")
group_comparisons_file <- file.path(base_path, "group_comparisons.csv") # DR AR sig
benci_signature_file <- file.path(base_path, "ifng_signature_benci_et_al.xlsx")
mitf_file <- file.path(base_path, "tirosh_et_al_mitf.xlsx")

# relevant signatures
signatures <- read_csv(signatures_file)
signatures_ar_dr <- read_csv(group_comparisons_file)
benci_sig <- read_xlsx(benci_signature_file, col_names = FALSE)
mitf <- read_xlsx(mitf_file, col_names = FALSE)

# Combine signatures data
combined_signatures <- rbind(signatures, signatures_ar_dr)

# Create pathways list
pathways <- split(combined_signatures$Gene, combined_signatures$Pathway)

# Add IFNG-related signatures
pathways[["IFNG_GS"]] <- unlist(benci_sig[, c(2, 4, 5, 6, 7)], use.names = FALSE)
pathways[["IFNG_RS"]] <- benci_sig[[1]]

# Setting up pathways
msigdbr_df <- msigdbr(species = "human", category = "H")
pathways_H = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

msigdbr_df <- msigdbr(species = "human", category = "C2")
pathways_R = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

# Convert Seurat object to a Monocle 3 CDS object
cds <- as.cell_data_set(seurat_obj)

# let's use the clustering information have
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 
list_cluster <- seurat_obj@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assign UMAP coordinate - cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- seurat_obj@reductions$umap@cell.embeddings

# Learn trajectory graph 
cds <- learn_graph(cds, use_partition = FALSE)

p1 <- plot_cells(cds,
           color_cells_by = 'cluster', 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 4)

# Order the cells in pseudotime 
root_cells <- colnames(cds)[cds$Malignant_clusters == "Mesenchymal_like"] # set mesenchymal as root node

# Order cells in the trajectory
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = root_cells)

p2 <- plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)


# cells ordered by monocle3 pseudotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
seurat_obj$Malignant_clusters
  
# score signatures - use AUcell because certain genes in DR sig have low counts
# Extract the normalized count matrix from Seurat object
expr_matrix <- GetAssayData(seurat_obj, layer = "scale.data", assay = "SCT")

# Define pathways as gene sets for AUCell
gene_sets <- list(
  AR = pathways[["AR_DR"]],
  DR = pathways[["DR_AR"]],
  PR = c(pathways[["PR_DR"]]),  
  IFNG_GS = pathways[["IFNG_GS"]],
  IFNG_RS = pathways[["IFNG_RS"]],
  MITF = mitf$...1,
  AXL = mitf$...2,
  IFNG_H = pathways_H[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
  MHCII = pathways_R[["REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION"]]
)

# Rank genes for each cell
cells_rankings <- AUCell_buildRankings(expr_matrix, nCores = 3, plotStats = FALSE)

# Compute AUC scores for each gene set
cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings, aucMaxRank = nrow(expr_matrix) * 0.75, normAUC = F)
seurat_obj <- AddMetaData(seurat_obj, as.data.frame(t(getAUC(cells_AUC))))

# plot pseudotime scores / cell cluster
ggplot(data.pseudo, aes(x = reorder(Malignant_clusters, monocle3_pseudotime, median), y = monocle3_pseudotime, fill = Malignant_clusters)) +
  geom_boxplot()+ theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") + ylab("Pseudotime") + xlab("") + 
  scale_fill_manual(values = pal_custom)

# Add signatures 
data.pseudo$DR <- seurat_obj$DR
data.pseudo$AR <- seurat_obj$AR
data.pseudo$AXL <- seurat_obj$AXL
data.pseudo$IFNG_RS <- seurat_obj$IFNG_RS
data.pseudo$IFNG_GS <- seurat_obj$IFNG_GS

# Define a theme with a compact legend
compact_legend_theme <- theme(
  legend.position = "right",         # Place legend on the right
  legend.key.size = unit(0.4, "cm"), # Reduce legend key size
  legend.text = element_text(size = 6),  # Reduce legend text size
  legend.title = element_text(size = 6, face = "bold"), # Smaller bold title
  legend.spacing.y = unit(5, "mm"),
  legend.spacing.x = unit(5, "mm")
)

# Plot
p3 <- ggplot(data.pseudo, aes(x = monocle3_pseudotime, y = AR)) +
  geom_point(aes(color = ident), alpha = 0.5, size = 0.3) +   # reduce point opacity to minimize overplotting
  geom_smooth(method = "loess", span = 1, se = FALSE, size = 1, color = "#ca2323") +  # increase smooth line width
  scale_color_manual(values = pal_custom) +
  labs(
    x = "Pseudotime",
    y = "AR signature",
    color = "Cell Type",
  ) +
  theme_classic(base_size = 8) + compact_legend_theme 

p4 <- ggplot(data.pseudo, aes(x = monocle3_pseudotime, y = DR)) +
  geom_point(aes(color = ident), alpha = 0.5, size = 0.3) +   # reduce point opacity to minimize overplotting
  geom_smooth(method = "loess", span = 0.76, se = TRUE, size = 1, color = "#ca2323") +  # increase smooth line width
  scale_color_manual(values = pal_custom) +
  labs(
    x = "Pseudotime",
    y = "DR signature",
    color = "Cell Type",
  ) +
  theme_classic(base_size = 8) + compact_legend_theme 

p5 <- ggplot(data.pseudo, aes(x = monocle3_pseudotime, y = AXL)) +
  geom_point(aes(color = ident), alpha = 0.5, size = 0.3) +   # reduce point opacity to minimize overplotting
  geom_smooth(method = "loess", span = 1, se = FALSE, size = 1, color = "#ca2323") +  # increase smooth line width
  scale_color_manual(values = pal_custom) +
  labs(
    x = "Pseudotime",
    y = "AXL signature",
    color = "Cell Type",
  ) +
  theme_classic(base_size = 8) + compact_legend_theme

# combine plot
combined_plot <- (p3 + p5 + p4) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")  

pdf("signatures_pseudo_plots.pdf", height = 2.8, width = 4.3)
combined_plot
dev.off()

################################################################################
## IFNG_GS vs IFNG_RS
################################################################################

# Define a function to plot with the Spectral color scale
custom_feature_plot <- function(seurat_obj, feature_name) {
  # Extract feature values
  feature_values <- seurat_obj@meta.data[[feature_name]]
  
  # Calculate min, median, and max values
  min_val <- min(feature_values, na.rm = TRUE)
  mean_val <- median(feature_values, na.rm = TRUE)
  max_val <- max(feature_values, na.rm = TRUE)
  
  # Generate FeaturePlot with dynamic color scale
  FeaturePlot(seurat_obj, features = feature_name) +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red", 
                           midpoint = mean_val, limits = c(min_val, max_val)) +
    ggtitle(paste("Feature Plot:", feature_name))
}

# Identify key malignant clusters
key_clusters <- c("Interferon_Alpha_Beta_Response", "Antigen_Presentation")

# Add a column to distinguish key clusters
seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(is_key_cluster = ifelse(Malignant_clusters %in% key_clusters, TRUE, FALSE))

# Create the plot
pdf("/mnt/8TB/Projects/POIAZ/Donagh/jc_marine/pseuodtime_infg_ifngrs.pdf", height = 3, width = 5.2)
ggplot(seurat_obj@meta.data, aes(x = IFNG_RS, y = IFNG_GS, color = Malignant_clusters)) +
  
  # Plot background points first (non-key clusters)
  geom_point(data = seurat_obj@meta.data %>% filter(!is_key_cluster), alpha = 0.3, size = 0.5) +
  
  # Highlight key clusters on top
  geom_point(data = seurat_obj@meta.data %>% filter(is_key_cluster), alpha = 0.7, size = 0.8) + 
  
  scale_color_manual(values = c(pal_nejm()(8), (pal_npg()(8)))) +
  theme_classic() +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8,
              color = "#252222", alpha = 0.5, linetype = "dashed") +
  labs(
    x = "IFNG_RS",  # Italicize labels if needed
    y = "IFNG_GS",
    color = "Malignant Cluster"
  )
dev.off()








