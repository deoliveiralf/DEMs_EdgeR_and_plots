# =============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS MODULE
# =============================================================================
# This module handles the differential expression analysis and result processing

#' Perform Comprehensive Differential Expression Analysis
#' 
#' @param fit Fitted model object from glmQLFit
#' @param contrast_matrix Contrast matrix
#' @param params List. Analysis parameters containing thresholds
#' @return List containing DE results and summary statistics
#' @export
perform_comprehensive_DE_analysis <- function(fit, contrast_matrix, params) {
  
  all_results <- list()
  summary_stats <- data.frame()
  
  for (i in 1:ncol(contrast_matrix)) {
    contrast_name <- colnames(contrast_matrix)[i]
    cat(paste("Processing contrast:", contrast_name, "\n"))
    
    # Perform QL F-test
    qlf <- glmQLFTest(fit, contrast = contrast_matrix[, i])
    results <- topTags(qlf, n = Inf)$table
    
    # Add regulation categories
    results$regulation <- "Not significant"
    results$regulation[results$logFC > params$LOGFC_THRESHOLD & 
                       results[[params$stat_column]] < params$stat_threshold] <- "Up-regulated"
    results$regulation[results$logFC < -params$LOGFC_THRESHOLD & 
                       results[[params$stat_column]] < params$stat_threshold] <- "Down-regulated"
    
    # Store results
    all_results[[contrast_name]] <- results
    
    # Calculate summary statistics
    n_up <- sum(results$regulation == "Up-regulated")
    n_down <- sum(results$regulation == "Down-regulated")
    n_total <- nrow(results)
    n_significant <- n_up + n_down
    
    # Add to summary table
    summary_stats <- rbind(summary_stats, 
                           data.frame(
                             Comparison = contrast_name,
                             Total_miRNAs = n_total,
                             Significant = n_significant,
                             Up_regulated = n_up,
                             Down_regulated = n_down,
                             Percent_DE = round((n_significant/n_total) * 100, 2),
                             Max_logFC_up = ifelse(n_up > 0, round(max(results$logFC[results$regulation == "Up-regulated"]), 2), 0),
                             Max_logFC_down = ifelse(n_down > 0, round(min(results$logFC[results$regulation == "Down-regulated"]), 2), 0),
                             stringsAsFactors = FALSE
                           ))
    
    # Save individual results
    output_file <- paste0("DE_results_", gsub("[^A-Za-z0-9]", "_", contrast_name), 
                          "_", params$stat_column, "_", params$stat_threshold, ".txt")
    write.table(results, output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  }
  
  # Save summary statistics
  summary_file <- paste0("DE_summary_", params$stat_column, "_", params$stat_threshold, ".csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  
  return(list(results = all_results, summary = summary_stats))
}

#' Extract Significantly Differentially Expressed miRNAs
#' 
#' @param de_results List. Results from differential expression analysis
#' @param params List. Analysis parameters
#' @return List of significantly DE miRNAs for each comparison
#' @export
extract_significant_mirnas <- function(de_results, params) {
  
  significant_lists <- list()
  
  for (comparison in names(de_results)) {
    results <- de_results[[comparison]]
    
    # Extract significantly DE miRNAs
    significant_mirnas <- rownames(results[abs(results$logFC) > params$LOGFC_THRESHOLD & 
                                            results[[params$stat_column]] < params$stat_threshold, ])
    
    significant_lists[[comparison]] <- significant_mirnas
  }
  
  return(significant_lists)
}

#' Create DE Results Summary Table
#' 
#' @param de_results List. Results from differential expression analysis
#' @param params List. Analysis parameters
#' @return Data frame with comprehensive summary
#' @export
create_de_summary_table <- function(de_results, params) {
  
  summary_data <- data.frame()
  
  for (comparison in names(de_results)) {
    results <- de_results[[comparison]]
    
    # Calculate statistics
    total_genes <- nrow(results)
    significant_genes <- sum(abs(results$logFC) > params$LOGFC_THRESHOLD & 
                             results[[params$stat_column]] < params$stat_threshold)
    up_regulated <- sum(results$logFC > params$LOGFC_THRESHOLD & 
                        results[[params$stat_column]] < params$stat_threshold)
    down_regulated <- sum(results$logFC < -params$LOGFC_THRESHOLD & 
                          results[[params$stat_column]] < params$stat_threshold)
    
    # Statistical summaries
    median_pvalue <- median(results$PValue, na.rm = TRUE)
    median_fdr <- median(results$FDR, na.rm = TRUE)
    
    # Add to summary
    summary_data <- rbind(summary_data, data.frame(
      Comparison = comparison,
      Total_miRNAs = total_genes,
      Significant_DEMs = significant_genes,
      Up_regulated = up_regulated,
      Down_regulated = down_regulated,
      Percent_DE = round((significant_genes/total_genes) * 100, 2),
      Median_PValue = round(median_pvalue, 4),
      Median_FDR = round(median_fdr, 4),
      Max_logFC = round(max(abs(results$logFC), na.rm = TRUE), 2),
      stringsAsFactors = FALSE
    ))
  }
  
  return(summary_data)
}

#' Get Top Differentially Expressed miRNAs
#' 
#' @param de_results List. Results from differential expression analysis
#' @param comparison Character. Name of the comparison
#' @param n_top Integer. Number of top miRNAs to return
#' @param params List. Analysis parameters
#' @return Data frame with top DE miRNAs
#' @export
get_top_de_mirnas <- function(de_results, comparison, n_top = 20, params) {
  
  if (!comparison %in% names(de_results)) {
    stop("Comparison not found in results")
  }
  
  results <- de_results[[comparison]]
  
  # Filter significant results
  significant_results <- results[abs(results$logFC) > params$LOGFC_THRESHOLD & 
                                 results[[params$stat_column]] < params$stat_threshold, ]
  
  if (nrow(significant_results) == 0) {
    cat(paste("No significant DEMs found for", comparison, "\n"))
    return(data.frame())
  }
  
  # Sort by statistical significance and fold change
  significant_results <- significant_results[order(significant_results[[params$stat_column]], 
                                                   -abs(significant_results$logFC)), ]
  
  # Return top n
  top_results <- head(significant_results, n_top)
  top_results$miRNA <- rownames(top_results)
  
  return(top_results)
}