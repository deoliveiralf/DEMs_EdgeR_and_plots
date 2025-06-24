# =============================================================================
# NORMALIZATION AND FILTERING MODULE
# =============================================================================
# This module handles normalization, dispersion estimation, and BCV plots

#' Normalize Data and Estimate Dispersions
#' 
#' @param dge_filtered DGEList object (filtered)
#' @param design Design matrix
#' @return List containing normalized DGE object, fitted model, AND CPM matrix
#' @export
normalize_and_estimate_dispersions <- function(dge_filtered, design) {
  
  # Estimate dispersions
  cat("Estimating dispersions...\n")
  dge_filtered <- estimateDisp(dge_filtered, design)
  
  # Calculate CPM
  norm_counts <- cpm(dge_filtered, normalized.lib.sizes = TRUE, log = FALSE)
  
  # Fit the model
  fit <- glmQLFit(dge_filtered, design)
  
  return(list(
    dge_filtered = dge_filtered,
    design = design,
    fit = fit,
    norm_counts = norm_counts  # Add normalized counts to output
  ))
}

#' Create Comprehensive BCV Plots
#' 
#' @param dge_obj DGEList object with estimated dispersions
#' @param filename Character. Output filename
#' @param design Design matrix (optional, will try to create from group if not provided)
#' @return NULL (saves plot to file)
#' @export
create_bcv_plots <- function(dge_obj, filename = "BCV_plots_comprehensive.png", design = NULL) {
  
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
  
  # 4. QL dispersion plot - Fixed section
  tryCatch({
    if (!is.null(design)) {
      # Use provided design matrix
      fit_temp <- glmQLFit(dge_obj, design)
    } else if ("group" %in% colnames(dge_obj$samples)) {
      # Create design matrix from group factor if not provided
      design_temp <- model.matrix(~0 + group, data = dge_obj$samples)
      colnames(design_temp) <- levels(dge_obj$samples$group)
      fit_temp <- glmQLFit(dge_obj, design_temp)
    } else {
      # Fallback: create simple intercept model
      design_temp <- model.matrix(~1, data = data.frame(row.names = colnames(dge_obj)))
      fit_temp <- glmQLFit(dge_obj, design_temp)
    }
    plotQLDisp(fit_temp, main = "QL Dispersion Estimates")
  }, error = function(e) {
    # If QL plot fails, create a placeholder
    plot(1, 1, type = "n", xlab = "", ylab = "", main = "QL Dispersion Plot - Error")
    text(1, 1, paste("Could not create QL plot:\n", e$message), cex = 0.8)
  })
  
  mtext("Dispersion Analysis", outer = TRUE, cex = 1.5, font = 2)
  dev.off()
  
  cat(paste("BCV plots saved as:", filename, "\n"))
  return(invisible(NULL))
}