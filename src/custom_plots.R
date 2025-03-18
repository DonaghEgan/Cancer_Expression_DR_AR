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
plot_boxplot <- function(data, x_var, y_var, y_label = NULL, x_label = NULL, alpha = NULL, 
                         comparisons = NULL, color_pal, test_method = "wilcox.test") {
  # Ensure required packages are loaded
  require(ggplot2)
  require(ggpubr)
  require(ggsci) # For scale_fill_npg()
  
  ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]])) + 
      geom_boxplot(outlier.shape = NA, alpha = alpha) +
      stat_compare_means(
      comparisons = comparisons,
      method = test_method, # More appropriate for non-normal distributions
      label = "p.format", # Show exact p-values instead of symbols
      bracket.size = 0.3,
      tip.length = 0.02,
      step.increase = 0.12,
      size = 3
    ) +
    scale_fill_manual(values = color_pal) +
    labs(
      x = x_label,
      y = y_label) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none")
}
