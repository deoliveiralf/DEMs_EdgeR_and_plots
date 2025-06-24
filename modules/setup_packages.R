# =============================================================================
# PACKAGE SETUP MODULE
# =============================================================================
# This module handles loading and installation of required packages

#' Load Required Packages
#' 
#' Installs and loads all packages required for the miRNA analysis
#' 
#' @return NULL (loads packages into environment)
#' @export
load_required_packages <- function() {
  
  # Define required packages
  required_packages <- c(
    "edgeR", "ggplot2", "dplyr", "pheatmap", "RColorBrewer", 
    "gridExtra", "knitr", "ggrepel", "VennDiagram", 
    "ComplexHeatmap", "tidyr", "viridis", "circlize"
  )
  
  # Function to install and load packages
  install_and_load <- function(packages) {
    for (pkg in packages) {
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat(paste("Installing package:", pkg, "\n"))
        
        # Try BiocManager first for Bioconductor packages
        if (pkg %in% c("edgeR", "ComplexHeatmap")) {
          if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
          }
          BiocManager::install(pkg, ask = FALSE, update = FALSE)
        } else {
          install.packages(pkg, dependencies = TRUE)
        }
        
        # Load the package
        library(pkg, character.only = TRUE)
      }
    }
  }
  
  # Install and load packages
  cat("Loading required packages...\n")
  install_and_load(required_packages)
  
  # Verify all packages are loaded
  loaded_packages <- sapply(required_packages, function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
  })
  
  if (all(loaded_packages)) {
    cat("All required packages loaded successfully!\n")
  } else {
    failed_packages <- names(loaded_packages)[!loaded_packages]
    warning(paste("Failed to load packages:", paste(failed_packages, collapse = ", ")))
  }
  
  return(invisible(NULL))
}

#' Get Session Information
#' 
#' Captures and saves session information for reproducibility
#' 
#' @param filename Character. Output filename for session info
#' @return NULL (saves session info to file)
#' @export
save_session_info <- function(filename = "session_info.txt") {
  sessionInfo_output <- capture.output(sessionInfo())
  writeLines(sessionInfo_output, filename)
  cat(paste("Session information saved to", filename, "\n"))
  return(invisible(NULL))
}