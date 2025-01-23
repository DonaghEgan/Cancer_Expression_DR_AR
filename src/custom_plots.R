# volcano Plot 
custom_volcano_plot <- function(volcano_data, subset_data, fc_threshold, padj_threshold, title = "NF vs AR") {
  ggplot(volcano_data, aes(x = logFC, y = negLogPadj)) +
    geom_point(aes(color = ifelse(abs(logFC) > fc_threshold & negLogPadj > -log10(padj_threshold), "red", "grey"))) +
    scale_color_identity() +
    xlim(c(-max(abs(volcano_data$logFC)), max(abs(volcano_data$logFC)))) +
    ylim(c(0, max(-log10(volcano_data$padj), na.rm = TRUE))) +
    geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "red") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "red") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_text(data = subset_data, aes(label = rownames(subset_data)),
              hjust = 1, vjust = -0.4, size = 2.5, color = "black", check_overlap = TRUE) +
    labs(x = "Log2 Fold Change", y = "-log10(Adjusted p-value)", title = title) +
    theme_bw()
}

# Define a reusable function for plotting
plot_boxplot <- function(data, y_var, y_label, comparisons, x_labels) {
  ggplot(data, aes(x = Subtype, y = !!sym(y_var), fill = Subtype)) + 
    geom_boxplot(alpha = 0.8) + 
    labs(x = "", y = y_label) + 
    theme_classic() + 
    scale_x_discrete(labels = x_labels) +
    scale_fill_npg() +
    stat_compare_means(label = "p.signif", comparisons = comparisons, size = 2.5) + 
    theme(legend.key.size = unit(0.3, "cm"), 
          axis.title.y = element_text(size = 9),
          axis.text.x = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8))
}
