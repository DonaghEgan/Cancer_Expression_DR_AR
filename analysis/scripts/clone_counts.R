library(readr)
library(data.table)
library(dplyr)
library(entropy)
library(ggplot2)
library(ggpubr)
library(ggsci)

# load colors
pal_colors <- rev(pal_npg()(3)[c(2,1,3)])

# Define the directory containing the files and the pattern to match
file_directory <- "Z:/working/GEL_MELANOMA/External_Datasets/HMF_GCP/analysis/PyClone/pyclone_output/"  
file_pattern <- "CPC.*\\.tsv$"  # Replace with your file pattern, e.g., "PyClone_output_.*\\.csv$"

# List all files matching the pattern in the specified directory
file_list <- list.files(path = file_directory, pattern = file_pattern, full.names = TRUE)

# Read each file into a list of dataframes
data_list <- lapply(file_list, read_tsv)

# Combine all dataframes into one dataframe
combined_data <- bind_rows(data_list)
dim(combined_data)

# Join cluster ID and Sample
combined_data$samp_clust_id <- paste(combined_data$sample_id, 
                                     combined_data$cluster_id, sep = "_")

count_table <- data.frame(samp_clust_id_count = table(combined_data$samp_clust_id))

combined_data_final <- merge(combined_data, count_table,
                             by.x = "samp_clust_id", 
                             by.y = "samp_clust_id_count.Var1")

pts <- data.frame(unique(combined_data_final$sample_id))

clone_counts_l<-list()
for(a in 1:nrow(pts)){
  
  these_muts <-combined_data_final[combined_data_final$sample_id == pts[a,1],4]
  these_clone_counts<-table(these_muts)
  
  print(table(these_muts))
  
  clone_counts_l[[a]]<-as.data.frame(these_clone_counts)
  clone_counts_l[[a]]$sample<-pts[a,1]
  
}

meta_data <- read.csv("Z:/working/GEL_MELANOMA/External_Datasets/HMF_GCP/RNA_seq_processed/OUTPUT/meta_filtered.csv",row.names = 1)
counts <- readRDS("Z:/working/GEL_MELANOMA/External_Datasets/HMF_GCP/RNA_seq_processed/OUTPUT/counts.Rds")

# intersect meta data with counts
meta_data <- meta_data[rownames(meta_data) %in% colnames(counts), ]## load data
meta_data$group <- factor(meta_data$group, levels = c("Primary_Resistance", "Acq_Resistance", "Nofailure"))

meta_pts <- merge(meta_data, pts, by.x = 0, by.y = "sample")

# max number of clones per sample
clone_numbers <- clone_counts %>%
  dplyr::group_by(sample) %>%
  summarize(clone_count = max(these_muts))

comp <- list(c("Acq_Resistance", "Nofailure"), c("Acq_Resistance", "Primary_Resistance"), c("Nofailure", "Primary_Resistance"))

clone_meta <- merge(meta_data, clone_numbers, by.x = 0, by.y = "sample")
pdf("Z:/working/GEL_MELANOMA/External_Datasets/HMF_GCP/RNA_seq_processed/OUTPUT/figures/clone_counts.pdf", height = 3, width = 3)
ggplot(clone_meta, aes(x = group, y = clone_count, color = group)) + 
  geom_point(position=position_jitter(width=0.05), alpha = 0.6) +  # Adding jittered points for visibility
  theme_classic() + 
  ylab("Clone Counts") + 
  theme(legend.position = "none") +
  stat_compare_means(comparisons = comp, method =  "wilcox.test", method.args = list(alternative = "greater")) +  stat_summary(
    geom = "crossbar",
    fun.y = "median",
    col = "black",
    size = 0.3, width = 0.3) + scale_color_manual(values = pal_colors) + ggtitle("HMF") + 
  scale_x_discrete(labels=c("Primary Res","Acquired Res","Durable Rsp")) + xlab("")
dev.off()













# Load libraries
library(readr)
library(data.table)
library(dplyr)
library(entropy)
library(ggplot2)
library(ggpubr)
library(ggsci)

# Load colors (using ggpubr/ggsci)
pal_colors <- rev(pal_npg()(3)[c(2,1,3)])

# Define directory and file pattern
file_directory <- "inset_path"  
file_pattern <- "CPC.*\\.tsv$"  

# List all files matching the pattern
file_list <- list.files(path = file_directory, pattern = file_pattern, full.names = TRUE)

# Load meta data and counts
meta_data <- read.csv("insert_path", row.names = 1)
counts <- readRDS("insert_path")

# Intersect meta_data with counts (keep only samples present in both)
meta_data <- meta_data[rownames(meta_data) %in% colnames(counts), ]
meta_data$group <- factor(meta_data$group, levels = c("Primary_Resistance", "Acq_Resistance", "Nofailure"))

# only keep samples with matched RNA
keep_names <- rownames(meta_data)

# Get the base filename without extension for each file
file_names <- tools::file_path_sans_ext(basename(file_list))

# Filter file_list where the extracted file name is in keep_names
filtered_file_list <- file_list[file_names %in% keep_names]

# Read each file into a list of dataframes
data_list <- lapply(filtered_file_list, read_tsv)

# Combine all dataframes into one dataframe
combined_data <- bind_rows(data_list)
dim(combined_data)

# Create a new identifier by pasting sample and cluster IDs together
combined_data$samp_clust_id <- paste(combined_data$sample_id, 
                                     combined_data$cluster_id, sep = "_")

# Create a count table for each sample-cluster combination
count_table <- data.frame(samp_clust_id_count = table(combined_data$samp_clust_id))

# Merge counts back into the combined data
combined_data_final <- merge(combined_data, count_table,
                             by.x = "samp_clust_id", 
                             by.y = "samp_clust_id_count.Var1")

# filter to ensure clones supported by > 0 mutations
combined_data_final<-combined_data_final[combined_data_final$samp_clust_id_count > 0, ]

# Calculate number of clones and number of corresponding mutations / sample
pts <- data.frame(unique(combined_data_final$sample_id))
clone_counts_l<-list()

for(a in 1:nrow(pts)){
  
  these_muts <-combined_data_final[combined_data_final$sample_id == pts[a,1],4]
  these_clone_counts<-table(these_muts)
  
  print(table(these_muts))
  
  clone_counts_l[[a]]<-as.data.frame(these_clone_counts)
  clone_counts_l[[a]]$sample<-pts[a,1]  
}

# renmame for merging
colnames(pts)[1] <- "sample"

# join as dataframe -> as make numeric
clone_counts <- do.call(rbind.data.frame,clone_counts_l)
clone_counts$these_muts <- as.numeric(clone_counts$these_muts)

# determine maximum amount of clones per sample
clone_numbers <- clone_counts %>%
  dplyr::group_by(sample) %>%
  dplyr::summarize(clone_count = max(these_muts))

dim(clone_numbers)
rownames(meta_data)

clone_meta <- merge(meta_data, clone_numbers, by.x = 0, by.y = "sample")
dim(clone_meta)

# comparisons 
comp <- list(c("Acq_Resistance", "Nofailure"),
c("Acq_Resistance", "Primary_Resistance"),
c("Primary_Resistance", "Nofailure"))

pdf("analysis/output/figures/clone_counts.pdf", height = 3, width = 3)
ggplot(clone_meta, aes(x = group, y = clone_count, color = group)) + 
  geom_point(position=position_jitter(width=0.05), alpha = 0.6) +  # Adding jittered points for visibility
  theme_classic() + 
  ylab("Clone Counts") + 
  theme(legend.position = "none") +
  stat_compare_means(comparisons = comp, method =  "wilcox.test") +  stat_summary(
    geom = "crossbar",
    fun.y = "median",
    col = "black",
    size = 0.3, width = 0.3) + scale_color_manual(values = pal_colors) + ggtitle("HMF") + 
  scale_x_discrete(labels=c("Primary Res","Acquired Res","Durable Rsp")) + xlab("")
dev.off()







# Load meta data and counts
meta_data <- read.csv("Z:/working/GEL_MELANOMA/External_Datasets/HMF_GCP/RNA_seq_processed/OUTPUT/meta_filtered.csv", row.names = 1)
counts <- readRDS("Z:/working/GEL_MELANOMA/External_Datasets/HMF_GCP/RNA_seq_processed/OUTPUT/counts.Rds")

# Intersect meta_data with counts (keep only samples present in both)
meta_data <- meta_data[rownames(meta_data) %in% colnames(counts), ]
meta_data$group <- factor(meta_data$group, levels = c("Primary_Resistance", "Acq_Resistance", "Nofailure"))

# Merge clone numbers with meta_data.
# Here we assume the rownames of meta_data correspond to sample IDs.
meta_data <- meta_data %>% 
  mutate(sample = rownames(meta_data))
clone_meta <- left_join(meta_data, clone_numbers, by = c("sample" = "sample_id"))

# Define comparisons for stat_compare_means
comp <- list(c("Acq_Resistance", "Nofailure"),
             c("Acq_Resistance", "Primary_Resistance"),
             c("Nofailure", "Primary_Resistance"))

# Create and save the plot
pdf("Z:/working/GEL_MELANOMA/External_Datasets/HMF_GCP/RNA_seq_processed/OUTPUT/figures/clone_counts.pdf", height = 3, width = 3)
ggplot(clone_meta, aes(x = group, y = clone_count, color = group)) + 
  geom_point(position = position_jitter(width = 0.05), alpha = 0.6) +  # Adding jittered points
  theme_classic() + 
  ylab("Clone Counts") + 
  theme(legend.position = "none") +
  stat_compare_means(comparisons = comp, 
                     method = "wilcox.test", 
                     method.args = list(alternative = "greater")) +  
  stat_summary(geom = "crossbar",
               fun.y = "median",
               col = "black",
               size = 0.3, width = 0.3) + 
  scale_color_manual(values = pal_colors) + 
  ggtitle("HMF") + 
  scale_x_discrete(labels = c("Primary Res", "Acquired Res", "Durable Rsp")) + 
  xlab("")
dev.off()
