# =============================================================================
# VOLCANO PLOTS MODULE
# =============================================================================
# This module contains functions for creating volcano plots

# Function to create enhanced volcano plots
create_enhanced_volcano <- function(results, title, stat_column, stat_threshold, logFC_threshold, top_n = 10) {
  
  # Add gene names for labeling
  results$gene <- rownames(results)
  
  # Calculate y-values
  if (stat_column == "FDR") {
    results$y_values <- -log10(pmax(results$FDR, 1e-300))  # Avoid -Inf
    y_label <- "-log₁₀(FDR)"
    threshold_line <- -log10(stat_threshold)
  } else {
    results$y_values <- -log10(pmax(results$PValue, 1e-300))  # Avoid -Inf
    y_label <- "-log₁₀(P-value)"
    threshold_line <- -log10(stat_threshold)
  }
  
  # Get top genes for labeling
  significant_genes <- results[results$regulation != "Not significant", ]
  
  if (nrow(significant_genes) > 0) {
    significant_genes <- significant_genes[order(significant_genes[[stat_column]]), ]
    if (nrow(significant_genes) > top_n) {
      top_genes <- significant_genes[1:top_n, ]
    } else {
      top_genes <- significant_genes
    }
  } else {
    top_genes <- data.frame()
  }
  
  # Create plot
  p <- ggplot(results, aes(x = logFC, y = y_values)) +
    geom_point(aes(color = regulation), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up-regulated" = "#E31A1C", 
                                  "Down-regulated" = "#1F78B4", 
                                  "Not significant" = "gray70")) +
    geom_hline(yintercept = threshold_line, linetype = "dashed", color = "black", alpha = 0.7) +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black", alpha = 0.7) +
    labs(title = title,
         subtitle = paste("Total:", nrow(results), "| Up:", sum(results$regulation == "Up-regulated"),
                          "| Down:", sum(results$regulation == "Down-regulated")),
         x = "log₂ Fold Change", y = y_label, color = "Regulation") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "bottom")
  
  # Add gene labels if available
  if (nrow(top_genes) > 0) {
    p <- p + geom_text_repel(data = top_genes, 
                             aes(x = logFC, y = y_values, label = gene), 
                             size = 2.5, max.overlaps = 20, 
                             box.padding = 0.3, point.padding = 0.3,
                             show.legend = FALSE)
  }
  
  return(p)
}

# Function to create and save all volcano plots
create_all_volcano_plots <- function(de_results, params) {  # Fixed: accepts params object
  
  cat("Creating volcano plots...\n")
  
  # Extract parameters from params object
  stat_column <- params$stat_column
  stat_threshold <- params$stat_threshold
  logFC_threshold <- params$LOGFC_THRESHOLD
  top_n <- params$TOP_N_VOLCANO_LABELS
  
  # Create individual volcano plots
  volcano_plots <- list()
  for (comparison in names(de_results)) {
    clean_name <- gsub("Group_", "", comparison)
    clean_name <- gsub("_vs_", " vs ", clean_name)
    
    volcano_plots[[comparison]] <- create_enhanced_volcano(
      de_results[[comparison]], 
      clean_name, 
      stat_column, 
      stat_threshold, 
      logFC_threshold,
      top_n  # Added missing parameter
    )
    
    # Save individual plot
    filename <- paste0("volcano_", gsub("[^A-Za-z0-9]", "_", comparison), "_", stat_column, "_", stat_threshold, ".png")
    ggsave(filename, volcano_plots[[comparison]], width = 8, height = 6, dpi = 300)
    cat(paste("Saved:", filename, "\n"))
  }
  
  # Create combined volcano plot
  if (length(volcano_plots) > 0) {
    combined_plot <- do.call(grid.arrange, c(volcano_plots, ncol = ceiling(sqrt(length(volcano_plots)))))
    combined_filename <- paste0("combined_volcano_plots_", stat_column, "_", stat_threshold, ".png")
    ggsave(combined_filename, combined_plot, width = 12, height = 8, dpi = 300)
    cat(paste("Saved combined plot:", combined_filename, "\n"))
  }
  
  return(volcano_plots)
}