# =============================================================================
# DIAGNOSTIC MODULE - modules/diagnostic_functions.R
# Enhanced miRNA Differential Expression Analysis - Diagnostic Functions
# =============================================================================

# Function to create Create BCV plots for quality control
#' Create BCV plots for quality control
#' 
#' @param dge_filtered DGEList object (filtered and normalized)
#' @param filename Character. Output filename for the plot
#' @return NULL (saves plot to file)
#' @export
create_bcv_plots <- function(dge_filtered, filename) {
  
  png(filename, width = 12, height = 8, units = "in", res = 300)
  par(mfrow = c(1, 2))
  
  # BCV plot
  plotBCV(dge_filtered, main = "Biological Coefficient of Variation")
  
  # MDS plot
  plotMDS(dge_filtered, main = "Multi-Dimensional Scaling Plot",
          col = as.numeric(factor(dge_filtered$samples$group)))
  
  dev.off()
  cat(paste("BCV plots saved to", filename, "\n"))
  
  return(invisible(NULL))
}

# Function to create safe color palettes
create_safe_color_palette <- function(categories, palette_name = "Set3") {
  n_categories <- length(unique(categories))
  
  if (n_categories == 1) {
    colors <- c("#8DD3C7")  # Single color
  } else if (n_categories == 2) {
    colors <- c("#8DD3C7", "#FFFFB3")  # Two colors
  } else if (n_categories <= 12) {
    colors <- RColorBrewer::brewer.pal(max(3, n_categories), palette_name)[1:n_categories]
  } else {
    # For more than 12 categories, use rainbow
    colors <- rainbow(n_categories)
  }
  
  # Ensure colors are named
  names(colors) <- unique(categories)
  return(colors)
}

# Function to run data structure diagnostic
run_data_diagnostic <- function(dge_obj) {
  cat("\n=== DATA STRUCTURE DIAGNOSTIC ===\n")
  cat("Sample groups:\n")
  print(table(dge_obj$samples$group))
  cat("\nTissue types:\n")
  print(table(dge_obj$samples$Tissue))
  cat("\nStage types:\n")
  print(table(dge_obj$samples$Stage))
  
  # Quick color palette test
  cat("\n=== COLOR PALETTE TEST ===\n")
  test_group_colors <- create_safe_color_palette(dge_obj$samples$group, "Set3")
  test_tissue_colors <- create_safe_color_palette(dge_obj$samples$Tissue, "Set1") 
  test_stage_colors <- create_safe_color_palette(dge_obj$samples$Stage, "Set2")
  
  cat("Group colors:\n")
  print(test_group_colors)
  cat("Tissue colors:\n")
  print(test_tissue_colors)
  cat("Stage colors:\n")
  print(test_stage_colors)
  
  # Test if all colors are properly named
  cat("All group colors named?", all(names(test_group_colors) %in% unique(dge_obj$samples$group)), "\n")
  cat("All tissue colors named?", all(names(test_tissue_colors) %in% unique(dge_obj$samples$Tissue)), "\n")
  cat("All stage colors named?", all(names(test_stage_colors) %in% unique(dge_obj$samples$Stage)), "\n")
  
  return(list(
    group_colors = test_group_colors,
    tissue_colors = test_tissue_colors,
    stage_colors = test_stage_colors
  ))
}