# =============================================================================
# MA PLOTS MODULE
# =============================================================================
# This module contains functions for creating MA plots

# Function to create MA plots for all contrasts
create_ma_plots <- function(fit, contrast_matrix, params) {  # Fixed: accepts params object
  
  cat("Creating MA plots...\n")
  
  # Extract parameters from params object
  stat_column <- params$stat_column
  stat_threshold <- params$stat_threshold
  logFC_threshold <- params$LOGFC_THRESHOLD
  
  filename <- paste0("MA_plots_", stat_column, "_", stat_threshold, ".png")
  png(filename, width = 3600, height = 2400, res = 300)
  
  # Calculate layout
  n_contrasts <- ncol(contrast_matrix)
  n_cols <- ceiling(sqrt(n_contrasts))
  n_rows <- ceiling(n_contrasts / n_cols)
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 2))
  
  for (i in 1:n_contrasts) {
    qlf <- glmQLFTest(fit, contrast = contrast_matrix[, i])
    clean_title <- gsub("Group_", "", colnames(contrast_matrix)[i])
    clean_title <- gsub("_vs_", " vs ", clean_title)
    plotMD(qlf, main = clean_title)
    abline(h = c(-logFC_threshold, logFC_threshold), col = "red", lty = 2)
  }
  
  dev.off()
  cat(paste("MA plots saved as:", filename, "\n"))
}


# Function to create BCV plots
create_bcv_plots <- function(dge_obj, filename = "BCV_plots_comprehensive.png") {
  
  cat("Creating BCV plots...\n")
  
  png(filename, width = 3600, height = 2400, res = 300)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2), oma = c(1, 1, 2, 1))
  
  # 1. Standard BCV plot
  plotBCV(dge_obj, main = "Biological Coefficient of Variation")
  
  # 2. Mean-variance plot
  plotMeanVar(dge_obj, show.raw.vars = TRUE, show.tagwise.vars = TRUE,
              show.ave.raw.vars = TRUE, main = "Mean-Variance Relationship")
  
  # 3. Dispersion estimates
  plot(dge_obj$AveLogCPM, sqrt(dge_obj$tagwise.dispersion), 
       pch = 16, cex = 0.5, col = "blue", 
       xlab = "Average log CPM", ylab = "Sqrt(Tagwise Dispersion)",
       main = "Tagwise Dispersion Estimates")
  abline(h = sqrt(dge_obj$common.dispersion), col = "red", lwd = 2)
  if (!is.null(dge_obj$trended.dispersion)) {
    lines(dge_obj$AveLogCPM[order(dge_obj$AveLogCPM)], 
          sqrt(dge_obj$trended.dispersion[order(dge_obj$AveLogCPM)]), 
          col = "green", lwd = 2)
  }
  legend("topright", legend = c("Tagwise", "Common", "Trended"), 
         col = c("blue", "red", "green"), lwd = c(1, 2, 2), pch = c(16, NA, NA))
  
  # 4. QL dispersion plot - Fixed: create design matrix here
  design <- model.matrix(~ 0 + group, data = dge_obj$samples)
  colnames(design) <- make.names(levels(dge_obj$samples$group))
  fit_temp <- glmQLFit(dge_obj, design)
  plotQLDisp(fit_temp, main = "QL Dispersion Estimates")
  
  mtext("Dispersion Analysis", outer = TRUE, cex = 1.5, font = 2)
  dev.off()
  
  cat(paste("BCV plots saved as:", filename, "\n"))
}