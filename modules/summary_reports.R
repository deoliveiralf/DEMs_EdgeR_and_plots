# =============================================================================
# SUMMARY REPORT MODULE - modules/summary_report.R
# Enhanced miRNA Differential Expression Analysis - Summary Report Functions
# =============================================================================

# Function to generate comprehensive summary (main entry point)
generate_comprehensive_summary <- function(raw_data = NULL, filtered_data, de_analysis = NULL, 
                                           de_results = NULL, params = NULL,
                                           stat_column = "FDR", stat_threshold = 0.05,
                                           LOGFC_THRESHOLD = 1, MIN_CPM = 1, MIN_SAMPLES = 3,
                                           dem_heatmaps = NULL, venn_results = NULL, ...) {
  
  cat("Generating comprehensive summary report...\n")
  
  # Input validation and debugging
  cat("Debugging input data...\n")
  cat("- raw_data:", ifelse(is.null(raw_data), "NULL", paste("exists, dim:", paste(dim(raw_data), collapse="x"))), "\n")
  cat("- filtered_data:", ifelse(is.null(filtered_data), "NULL", paste("exists, dim:", paste(dim(filtered_data), collapse="x"))), "\n")
  cat("- de_analysis:", ifelse(is.null(de_analysis), "NULL", "exists"), "\n")
  cat("- de_results:", ifelse(is.null(de_results), "NULL", "exists"), "\n")
  cat("- params:", ifelse(is.null(params), "NULL", "exists"), "\n")
  
  # Enhanced debugging for filtered_data
  if (!is.null(filtered_data)) {
    cat("- filtered_data class:", paste(class(filtered_data), collapse = ", "), "\n")
    cat("- filtered_data structure:\n")
    str(filtered_data, max.level = 1)
  }
  
  # Validate essential inputs
  if (is.null(filtered_data)) {
    stop("Error: filtered_data is NULL. Cannot proceed with summary report.")
  }
  
  # Improved validation logic
  is_valid_data <- FALSE
  data_info <- list(type = "unknown", nrow = 0, ncol = 0)
  
  # Check for DGEList object
  if (inherits(filtered_data, "DGEList")) {
    is_valid_data <- TRUE
    data_info$type <- "DGEList"
    if (!is.null(filtered_data$counts)) {
      data_info$nrow <- nrow(filtered_data$counts)
      data_info$ncol <- ncol(filtered_data$counts)
    }
    cat("- Detected DGEList object\n")
  }
  # Check for data.frame
  else if (is.data.frame(filtered_data)) {
    is_valid_data <- TRUE
    data_info$type <- "data.frame"
    data_info$nrow <- nrow(filtered_data)
    data_info$ncol <- ncol(filtered_data)
    cat("- Detected data.frame object\n")
  }
  # Check for matrix
  else if (is.matrix(filtered_data)) {
    is_valid_data <- TRUE
    data_info$type <- "matrix"
    data_info$nrow <- nrow(filtered_data)
    data_info$ncol <- ncol(filtered_data)
    cat("- Detected matrix object\n")
  }
  # Check for list with known elements
  else if (is.list(filtered_data)) {
    cat("- List elements:", paste(names(filtered_data), collapse = ", "), "\n")
    
    # Check for common expression data elements
    if (!is.null(filtered_data$counts)) {
      is_valid_data <- TRUE
      data_info$type <- "list_with_counts"
      data_info$nrow <- nrow(filtered_data$counts)
      data_info$ncol <- ncol(filtered_data$counts)
      cat("- Detected list object with counts element\n")
    } else if (!is.null(filtered_data$E)) {
      is_valid_data <- TRUE
      data_info$type <- "list_with_E"
      data_info$nrow <- nrow(filtered_data$E)
      data_info$ncol <- ncol(filtered_data$E)
      cat("- Detected list object with E element\n")
    } else if (!is.null(filtered_data$logCPM)) {
      is_valid_data <- TRUE
      data_info$type <- "list_with_logCPM"
      data_info$nrow <- nrow(filtered_data$logCPM)
      data_info$ncol <- ncol(filtered_data$logCPM)
      cat("- Detected list object with logCPM element\n")
    } else if (!is.null(filtered_data$normalized)) {
      is_valid_data <- TRUE
      data_info$type <- "list_with_normalized"
      data_info$nrow <- nrow(filtered_data$normalized)
      data_info$ncol <- ncol(filtered_data$normalized)
      cat("- Detected list object with normalized element\n")
    } else if (!is.null(filtered_data$data)) {
      is_valid_data <- TRUE
      data_info$type <- "list_with_data"
      data_info$nrow <- nrow(filtered_data$data)
      data_info$ncol <- ncol(filtered_data$data)
      cat("- Detected list object with data element\n")
    } else {
      # Look for any matrix-like elements
      matrix_elements <- sapply(filtered_data, function(x) is.matrix(x) || is.data.frame(x))
      if (any(matrix_elements)) {
        first_matrix_name <- names(filtered_data)[which(matrix_elements)[1]]
        first_matrix <- filtered_data[[first_matrix_name]]
        is_valid_data <- TRUE
        data_info$type <- paste0("list_with_", first_matrix_name)
        data_info$nrow <- nrow(first_matrix)
        data_info$ncol <- ncol(first_matrix)
        cat("- Detected list object with matrix element:", first_matrix_name, "\n")
      }
    }
  }
  # Check if it's a numeric structure with dimensions
  else if (is.numeric(filtered_data) && !is.null(dim(filtered_data))) {
    is_valid_data <- TRUE
    data_info$type <- "numeric_matrix"
    data_info$nrow <- nrow(filtered_data)
    data_info$ncol <- ncol(filtered_data)
    cat("- Detected numeric object with dimensions\n")
  }
  
  if (!is_valid_data) {
    cat("Warning: Could not recognize data format. Data info:\n")
    cat("- Class:", paste(class(filtered_data), collapse = ", "), "\n")
    cat("- Length:", length(filtered_data), "\n")
    if (is.list(filtered_data)) {
      cat("- List names:", paste(names(filtered_data), collapse = ", "), "\n")
    }
    cat("Attempting to proceed with best guess...\n")
  }
  
  # Handle different input formats - prioritize de_results if provided
  if (!is.null(de_results)) {
    de_analysis <- de_results
    cat("Using de_results as de_analysis\n")
  }
  
  # Extract parameters from params object if provided
  if (!is.null(params)) {
    cat("Extracting parameters from params object...\n")
    if (!is.null(params$stat_column)) stat_column <- params$stat_column
    if (!is.null(params$stat_threshold)) stat_threshold <- params$stat_threshold
    if (!is.null(params$LOGFC_THRESHOLD)) LOGFC_THRESHOLD <- params$LOGFC_THRESHOLD
    if (!is.null(params$MIN_CPM)) MIN_CPM <- params$MIN_CPM
    if (!is.null(params$MIN_SAMPLES)) MIN_SAMPLES <- params$MIN_SAMPLES
  }
  
  # Get original count from raw data with safe handling
  original_dge_count <- 0
  tryCatch({
    if (!is.null(raw_data)) {
      if (is.list(raw_data) && !is.null(raw_data$counts)) {
        original_dge_count <- nrow(raw_data$counts)
      } else if (!is.null(dim(raw_data))) {
        original_dge_count <- nrow(raw_data)
      }
    } else if (data_info$nrow > 0) {
      original_dge_count <- data_info$nrow
    }
    cat("Original DGE count:", original_dge_count, "\n")
  }, error = function(e) {
    cat("Warning: Could not determine original count, using 0:", e$message, "\n")
    original_dge_count <- 0
  })
  
  # Create the main summary report with error handling
  tryCatch({
    create_summary_report(
      dge_filtered = filtered_data,
      de_analysis = de_analysis,
      stat_column = stat_column,
      stat_threshold = stat_threshold,
      LOGFC_THRESHOLD = LOGFC_THRESHOLD,
      MIN_CPM = MIN_CPM,
      MIN_SAMPLES = MIN_SAMPLES,
      dem_heatmaps = dem_heatmaps,
      venn_results = venn_results,
      original_dge_count = original_dge_count,
      data_info = data_info
    )
  }, error = function(e) {
    cat("Error in create_summary_report:", e$message, "\n")
    cat("Attempting simplified summary...\n")
    create_simple_summary(filtered_data, de_analysis, original_dge_count, data_info)
  })
  
  # Save detailed results summary to file with error handling
  if (!is.null(de_analysis)) {
    tryCatch({
      save_results_summary(
        de_analysis = de_analysis,
        stat_column = stat_column,
        stat_threshold = stat_threshold,
        output_file = "comprehensive_results_summary.txt"
      )
    }, error = function(e) {
      cat("Error saving results summary:", e$message, "\n")
    })
  } else {
    cat("No DE analysis results to save\n")
  }
  
  cat("Summary report generation completed!\n")
  
  # Additional debug info
  cat("DEBUG: Exiting generate_comprehensive_summary function normally\n")
}

# Function to create comprehensive summary report
create_summary_report <- function(dge_filtered, de_analysis, stat_column, stat_threshold, 
                                  LOGFC_THRESHOLD, MIN_CPM, MIN_SAMPLES, dem_heatmaps = NULL, 
                                  venn_results = NULL, original_dge_count, data_info = NULL) {
  
  cat("Creating summary report...\n")
  
  # Validate inputs
  if (is.null(dge_filtered)) {
    cat("Error: dge_filtered is NULL\n")
    return(invisible(NULL))
  }
  
  # Create separator line safely using strrep
  separator_line <- strrep("=", 80)
  
  # Print comprehensive summary
  cat("\n")
  cat(separator_line, "\n")
  cat("COMPREHENSIVE miRNA DIFFERENTIAL EXPRESSION ANALYSIS SUMMARY\n")
  cat(separator_line, "\n")
  
  cat("\nANALYSIS PARAMETERS:\n")
  cat(paste("- Statistical threshold:", stat_column, "<", stat_threshold, "\n"))
  cat(paste("- Log fold change threshold: |logFC| >", LOGFC_THRESHOLD, "\n"))
  cat(paste("- Minimum CPM for filtering:", MIN_CPM, "\n"))
  cat(paste("- Minimum samples with CPM >", MIN_CPM, ":", MIN_SAMPLES, "\n"))
  
  cat("\nDATA SUMMARY:\n")
  cat(paste("- Total miRNAs loaded:", ifelse(is.null(original_dge_count), "Unknown", original_dge_count), "\n"))
  
  # Use data_info if available, otherwise try to extract dimensions
  if (!is.null(data_info) && data_info$nrow > 0) {
    cat(paste("- miRNAs after filtering:", data_info$nrow, "\n"))
    cat(paste("- Total samples:", data_info$ncol, "\n"))
    cat(paste("- Data type:", data_info$type, "\n"))
  } else {
    # Fallback to original dimension extraction logic
    tryCatch({
      if (inherits(dge_filtered, "DGEList")) {
        if (!is.null(dge_filtered$counts)) {
          filtered_nrow <- nrow(dge_filtered$counts)
          filtered_ncol <- ncol(dge_filtered$counts)
        } else {
          filtered_nrow <- nrow(dge_filtered)
          filtered_ncol <- ncol(dge_filtered)
        }
      } else if (is.list(dge_filtered)) {
        # Try different common elements for list objects
        if (!is.null(dge_filtered$counts)) {
          filtered_nrow <- nrow(dge_filtered$counts)
          filtered_ncol <- ncol(dge_filtered$counts)
        } else if (!is.null(dge_filtered$E)) {
          filtered_nrow <- nrow(dge_filtered$E)
          filtered_ncol <- ncol(dge_filtered$E)
        } else if (!is.null(dge_filtered$logCPM)) {
          filtered_nrow <- nrow(dge_filtered$logCPM)
          filtered_ncol <- ncol(dge_filtered$logCPM)
        } else if (!is.null(dge_filtered$normalized)) {
          filtered_nrow <- nrow(dge_filtered$normalized)
          filtered_ncol <- ncol(dge_filtered$normalized)
        } else if (!is.null(dge_filtered$data)) {
          filtered_nrow <- nrow(dge_filtered$data)
          filtered_ncol <- ncol(dge_filtered$data)
        } else {
          # If it's a list but no recognizable elements, try to find any matrix-like elements
          matrix_elements <- sapply(dge_filtered, function(x) is.matrix(x) || is.data.frame(x))
          if (any(matrix_elements)) {
            first_matrix <- dge_filtered[[which(matrix_elements)[1]]]
            filtered_nrow <- nrow(first_matrix)
            filtered_ncol <- ncol(first_matrix)
          } else {
            filtered_nrow <- NULL
            filtered_ncol <- NULL
          }
        }
      } else {
        filtered_nrow <- nrow(dge_filtered)
        filtered_ncol <- ncol(dge_filtered)
      }
      
      cat(paste("- miRNAs after filtering:", ifelse(is.null(filtered_nrow), "Unknown", filtered_nrow), "\n"))
      cat(paste("- Total samples:", ifelse(is.null(filtered_ncol), "Unknown", filtered_ncol), "\n"))
    }, error = function(e) {
      cat("- miRNAs after filtering: Error determining count -", e$message, "\n")
      cat("- Total samples: Error determining count\n")
    })
  }
  
  # Safe handling of sample groups
  tryCatch({
    if (inherits(dge_filtered, "DGEList") && !is.null(dge_filtered$samples) && !is.null(dge_filtered$samples$group)) {
      groups <- levels(as.factor(dge_filtered$samples$group))
      cat(paste("- Sample groups:", paste(groups, collapse = ", "), "\n"))
    } else if (is.list(dge_filtered) && !is.null(dge_filtered$samples) && !is.null(dge_filtered$samples$group)) {
      groups <- levels(as.factor(dge_filtered$samples$group))
      cat(paste("- Sample groups:", paste(groups, collapse = ", "), "\n"))
    } else {
      cat("- Sample groups: Not available\n")
    }
  }, error = function(e) {
    cat("- Sample groups: Error determining groups -", e$message, "\n")
  })
  
  cat("\nDIFFERENTIAL EXPRESSION RESULTS:\n")
  if (!is.null(de_analysis)) {
    tryCatch({
      if (!is.null(de_analysis$summary)) {
        print(de_analysis$summary)
      } else if (!is.null(de_analysis$results)) {
        cat("DE results available for", length(de_analysis$results), "comparisons\n")
        for (comp_name in names(de_analysis$results)) {
          comp_results <- de_analysis$results[[comp_name]]
          if (!is.null(comp_results) && nrow(comp_results) > 0) {
            n_sig <- sum(comp_results$regulation != "Not significant", na.rm = TRUE)
            cat(paste("- ", comp_name, ": ", n_sig, " significant DEMs\n"))
          }
        }
      } else {
        cat("DE analysis structure not recognized\n")
      }
    }, error = function(e) {
      cat("Error displaying DE results:", e$message, "\n")
    })
  } else {
    cat("No DE analysis results available\n")
  }
  
  cat("\nHEATMAP SUMMARY:\n")
  if (!is.null(dem_heatmaps) && length(dem_heatmaps) > 0) {
    tryCatch({
      for (comp in names(dem_heatmaps)) {
        if (!is.null(dem_heatmaps[[comp]]$n_dems)) {
          cat(paste("-", comp, ":", dem_heatmaps[[comp]]$n_dems, "DEMs\n"))
        } else {
          cat(paste("-", comp, ": Heatmap created\n"))
        }
      }
    }, error = function(e) {
      cat("Error displaying heatmap summary:", e$message, "\n")
    })
  } else {
    cat("No DEM heatmaps created\n")
  }
  
  cat("\nOUTPUT FILES CREATED:\n")
  cat("- BCV_plots_comprehensive.png\n")
  cat("- global_expression_heatmap_fixed.png\n")
  cat("- global_expression_heatmap_fixed_complex.png\n")
  cat("- volcano plots for each comparison\n")
  cat(paste("- MA_plots_", stat_column, "_", stat_threshold, ".png\n"))
  cat("- DEM heatmaps for each comparison (if DEMs found)\n")
  cat("- DE results tables for each comparison\n")
  cat("- Summary CSV files\n")
  
  if (!is.null(venn_results)) {
    cat("- venn_diagram_DEMs.png\n")
    cat("- venn_diagram_summary.csv\n")
  }
  
  cat("\n")
  cat(separator_line, "\n")
  cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
  cat(separator_line, "\n")
  
  # Save session info
  tryCatch({
    sessionInfo_output <- capture.output(sessionInfo())
    writeLines(sessionInfo_output, "session_info.txt")
    cat("\nSession information saved to session_info.txt\n")
  }, error = function(e) {
    cat("\nWarning: Could not save session info:", e$message, "\n")
  })
}

# Function to create simple summary when main function fails
create_simple_summary <- function(filtered_data, de_analysis, original_count, data_info = NULL) {
  # Create separator line safely using strrep
  simple_separator <- strrep("=", 25)
  
  cat("\n", simple_separator, "\n")
  cat("SIMPLIFIED SUMMARY\n")
  cat(simple_separator, "\n")
  cat("Analysis completed with some limitations\n")
  
  if (!is.null(data_info) && data_info$nrow > 0) {
    cat("Filtered data dimensions:", data_info$nrow, "x", data_info$ncol, "\n")
    cat("Data type:", data_info$type, "\n")
  } else if (!is.null(filtered_data)) {
    tryCatch({
      # Try to get dimensions safely
      if (inherits(filtered_data, "DGEList") && !is.null(filtered_data$counts)) {
        dims <- dim(filtered_data$counts)
      } else if (is.list(filtered_data)) {
        # Try different common elements for list objects
        if (!is.null(filtered_data$counts)) {
          dims <- dim(filtered_data$counts)
        } else if (!is.null(filtered_data$E)) {
          dims <- dim(filtered_data$E)
        } else if (!is.null(filtered_data$logCPM)) {
          dims <- dim(filtered_data$logCPM)
        } else if (!is.null(filtered_data$normalized)) {
          dims <- dim(filtered_data$normalized)
        } else if (!is.null(filtered_data$data)) {
          dims <- dim(filtered_data$data)
        } else {
          dims <- NULL
        }
      } else {
        dims <- dim(filtered_data)
      }
      
      if (!is.null(dims)) {
        cat("Filtered data dimensions:", paste(dims, collapse = " x "), "\n")
      } else {
        cat("Could not determine filtered data dimensions\n")
      }
      
      # Print list structure for debugging
      if (is.list(filtered_data)) {
        cat("List elements:", paste(names(filtered_data), collapse = ", "), "\n")
      }
      
    }, error = function(e) {
      cat("Could not determine filtered data dimensions:", e$message, "\n")
    })
  }
  
  if (!is.null(de_analysis)) {
    cat("DE analysis results are available\n")
  } else {
    cat("No DE analysis results\n")
  }
  
  cat(simple_separator, "\n")
}

# Function to save comprehensive results summary to file
save_results_summary <- function(de_analysis, stat_column, stat_threshold, output_file = "comprehensive_results_summary.txt") {
  
  tryCatch({
    # Create separator line safely using strrep
    file_separator <- strrep("=", 77)
    sub_separator <- strrep("-", 32)
    
    # Create summary text
    summary_text <- c(
      file_separator,
      "COMPREHENSIVE miRNA DIFFERENTIAL EXPRESSION ANALYSIS RESULTS SUMMARY",
      file_separator,
      "",
      paste("Analysis Date:", Sys.Date()),
      paste("Statistical Threshold:", stat_column, "<", stat_threshold),
      "",
      "DIFFERENTIAL EXPRESSION SUMMARY:",
      sub_separator
    )
    
    # Add detailed results for each comparison
    if (!is.null(de_analysis$results)) {
      for (comparison in names(de_analysis$results)) {
        results <- de_analysis$results[[comparison]]
        
        # Count significant genes
        n_up <- sum(results$regulation == "Up-regulated", na.rm = TRUE)
        n_down <- sum(results$regulation == "Down-regulated", na.rm = TRUE)
        n_total <- nrow(results)
        
        summary_text <- c(summary_text,
                          "",
                          paste("Comparison:", comparison),
                          paste("  Total miRNAs:", n_total),
                          paste("  Up-regulated:", n_up),
                          paste("  Down-regulated:", n_down),
                          paste("  Total significant:", n_up + n_down),
                          paste("  Percentage DE:", round(((n_up + n_down)/n_total) * 100, 2), "%"))
        
        # Add top DE genes if any
        if (n_up > 0 || n_down > 0) {
          significant_genes <- results[results$regulation != "Not significant", ]
          significant_genes <- significant_genes[order(significant_genes[[stat_column]]), ]
          
          summary_text <- c(summary_text,
                            paste("  Top 5 most significant:"))
          
          top_genes <- head(significant_genes, 5)
          for (i in 1:nrow(top_genes)) {
            gene_name <- rownames(top_genes)[i]
            logfc <- round(top_genes$logFC[i], 3)
            stat_val <- formatC(top_genes[[stat_column]][i], format = "e", digits = 2)
            regulation <- top_genes$regulation[i]
            
            summary_text <- c(summary_text,
                              paste("    ", gene_name, "| logFC:", logfc, "|", stat_column, ":", stat_val, "|", regulation))
          }
        }
      }
    }
    
    # Write to file
    writeLines(summary_text, output_file)
    cat(paste("Comprehensive results summary saved to:", output_file, "\n"))
    
  }, error = function(e) {
    cat("Error saving results summary:", e$message, "\n")
  })
}

# Function to print analysis parameters summary
print_analysis_parameters <- function(stat_column, stat_threshold, LOGFC_THRESHOLD, MIN_CPM, MIN_SAMPLES) {
  # Create separator line safely using strrep
  param_separator <- strrep("-", 24)
  
  cat("\nANALYSIS PARAMETERS USED:\n")
  cat(param_separator, "\n")
  cat(paste("Statistical threshold:", stat_column, "<", stat_threshold, "\n"))
  cat(paste("Log fold change threshold: |logFC| >", LOGFC_THRESHOLD, "\n"))
  cat(paste("Minimum CPM for filtering:", MIN_CPM, "\n"))
  cat(paste("Minimum samples with CPM >", MIN_CPM, ":", MIN_SAMPLES, "\n"))
}