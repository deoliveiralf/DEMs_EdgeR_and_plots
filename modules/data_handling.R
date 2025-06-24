# =============================================================================
# DATA HANDLING MODULE - FIXED VERSION WITH ROBUST SAMPLE MATCHING
# =============================================================================
# This module handles data loading, validation, and initial processing

#' Safely Read Files with Error Handling
#' 
#' @param filename Character. Path to the file
#' @param type Character. File type ("csv" or "tsv")
#' @return Data frame with file contents
#' @export
safe_read_file <- function(filename, type = "csv") {
  if (!file.exists(filename)) {
    stop(paste("File not found:", filename))
  }
  
  tryCatch({
    if (type == "csv") {
      return(read.csv(filename, stringsAsFactors = FALSE))
    } else if (type == "tsv") {
      return(read.delim(filename, row.names = 1, check.names = FALSE))
    }
  }, error = function(e) {
    stop(paste("Error reading file", filename, ":", e$message))
  })
}

#' Load and Validate Data Files with Robust Sample Matching
#' 
#' @param data_files List. Named list with sample_info and count_matrix file paths
#' @return List containing validated sample_info and count_matrix
#' @export
load_and_validate_data <- function(data_files) {
  
  # Read the data files
  cat("Loading data files...\n")
  sample_info <- safe_read_file(data_files$sample_info, "csv")
  count_matrix <- safe_read_file(data_files$count_matrix, "tsv")
  
  # Validate data
  if (nrow(sample_info) == 0 || ncol(count_matrix) == 0) {
    stop("Error: Empty data files detected")
  }
  
  cat(paste("Loaded", nrow(count_matrix), "miRNAs and", ncol(count_matrix), "samples from count matrix\n"))
  cat(paste("Loaded", nrow(sample_info), "samples from sample info\n"))
  
  # Check for sample matching
  samples_in_info <- sample_info$Sample
  samples_in_matrix <- colnames(count_matrix)
  
  # Find common samples
  common_samples <- intersect(samples_in_info, samples_in_matrix)
  missing_in_info <- setdiff(samples_in_matrix, samples_in_info)
  missing_in_matrix <- setdiff(samples_in_info, samples_in_matrix)
  
  # Report sample matching status
  cat("\n=== SAMPLE MATCHING REPORT ===\n")
  cat(paste("Samples in both files:", length(common_samples), "\n"))
  cat(paste("Samples in count matrix only:", length(missing_in_info), "\n"))
  cat(paste("Samples in sample info only:", length(missing_in_matrix), "\n"))
  
  if (length(missing_in_info) > 0) {
    cat("Samples in count matrix but not in sample info:\n")
    cat(paste("  ", missing_in_info, collapse = "\n"), "\n")
  }
  
  if (length(missing_in_matrix) > 0) {
    cat("Samples in sample info but not in count matrix:\n")
    cat(paste("  ", missing_in_matrix, collapse = "\n"), "\n")
  }
  
  # Check if we have enough common samples
  if (length(common_samples) < 6) {
    stop(paste("ERROR: Only", length(common_samples), "matching samples found. Need at least 6 for analysis."))
  }
  
  # Filter to common samples only
  cat(paste("\nFiltering to", length(common_samples), "common samples...\n"))
  
  # Filter and reorder sample_info to match count_matrix
  sample_info_filtered <- sample_info[sample_info$Sample %in% common_samples, ]
  count_matrix_filtered <- count_matrix[, common_samples]
  
  # Ensure sample order matches between files
  sample_info_filtered <- sample_info_filtered[match(colnames(count_matrix_filtered), sample_info_filtered$Sample), ]
  
  # Final validation
  if (!all(colnames(count_matrix_filtered) == sample_info_filtered$Sample)) {
    stop("ERROR: Sample order mismatch after filtering")
  }
  
  # Validate that required columns exist
  required_sample_cols <- c("Sample", "Group")
  missing_cols <- setdiff(required_sample_cols, colnames(sample_info_filtered))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in sample_info:", paste(missing_cols, collapse = ", ")))
  }
  
  # Check for numeric count data
  if (!all(sapply(count_matrix_filtered, is.numeric))) {
    warning("Non-numeric values detected in count matrix")
  }
  
  cat(paste("Final dataset:", nrow(count_matrix_filtered), "miRNAs x", ncol(count_matrix_filtered), "samples\n"))
  
  return(list(
    sample_info = sample_info_filtered,
    count_matrix = count_matrix_filtered
  ))
}

#' Create and Filter DGE Object
#' 
#' @param count_matrix Matrix. Count matrix with miRNAs as rows and samples as columns
#' @param sample_info Data frame. Sample information
#' @param params List. Analysis parameters
#' @return List containing filtered DGE object and design matrix
#' @export
create_and_filter_dge <- function(count_matrix, sample_info, params) {
  
  # Verify dimensions match
  if (ncol(count_matrix) != nrow(sample_info)) {
    stop(paste("Dimension mismatch: count_matrix has", ncol(count_matrix), 
               "samples but sample_info has", nrow(sample_info), "samples"))
  }
  
  # Verify sample order matches
  if (!all(colnames(count_matrix) == sample_info$Sample)) {
    stop("Sample order mismatch between count_matrix and sample_info")
  }
  
  # Create DGEList object
  dge <- DGEList(counts = count_matrix, samples = sample_info)
  # Simplified group assignment - just use the Group column directly
  dge$samples$group <- factor(dge$samples$Group)
  
  # Display group information
  cat("\nSample groups:\n")
  group_table <- table(dge$samples$group)
  print(group_table)
  
  # Check if we have enough samples per group
  min_samples_per_group <- min(group_table)
  if (min_samples_per_group < 2) {
    warning(paste("Some groups have fewer than 2 samples. Minimum:", min_samples_per_group))
  }
  
  # Filter low-expressed genes
  cat(paste("\nFiltering genes with CPM >", params$MIN_CPM, "in at least", params$MIN_SAMPLES, "samples...\n"))
  keep <- rowSums(cpm(dge) > params$MIN_CPM) >= params$MIN_SAMPLES
  dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
  dge_filtered <- calcNormFactors(dge_filtered)
  
  cat(paste("Retained", nrow(dge_filtered), "out of", nrow(dge), "miRNAs\n"))
  cat(paste("Filtering rate:", round((1 - nrow(dge_filtered)/nrow(dge)) * 100, 1), "%\n"))
  
  # Create design matrix - simplified for group-only design
  design <- model.matrix(~ 0 + group, data = dge_filtered$samples)
  colnames(design) <- levels(dge_filtered$samples$group)
  
  # Verify design matrix dimensions
  cat(paste("Design matrix dimensions:", nrow(design), "x", ncol(design), "\n"))
  cat(paste("DGE object dimensions:", nrow(dge_filtered), "x", ncol(dge_filtered), "\n"))
  
  if (nrow(design) != ncol(dge_filtered)) {
    stop(paste("Design matrix has", nrow(design), "rows but DGE object has", ncol(dge_filtered), "columns"))
  }
  
  return(list(
    dge_original = dge,
    dge_filtered = dge_filtered,
    design = design
  ))
}

#' Print Data Structure Diagnostic
#' 
#' @param dge_filtered DGEList object (filtered)
#' @return NULL (prints diagnostic information)
#' @export
print_data_diagnostic <- function(dge_filtered) {
  cat("\n=== DATA STRUCTURE DIAGNOSTIC ===\n")
  cat("Sample groups:\n")
  print(table(dge_filtered$samples$group))
  cat("\nGroup distribution:\n")
  print(table(dge_filtered$samples$Group))
  cat("\nSample names:\n")
  print(dge_filtered$samples$Sample)
  cat("\nDGE object dimensions:", nrow(dge_filtered), "x", ncol(dge_filtered), "\n")
  
  return(invisible(NULL))
}