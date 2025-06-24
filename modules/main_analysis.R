# =============================================================================
# MAIN miRNA DIFFERENTIAL EXPRESSION ANALYSIS SCRIPT - UPDATED FOR GROUPS ONLY
# =============================================================================
# This is the main script that orchestrates the entire analysis
# by calling functions from modularized files

# Clear environment
rm(list = ls())

# Source all module files
source("modules/setup_packages.R")
source("modules/data_handling.R")
source("modules/normalization_filtering.R")
source("modules/differential_expression.R")
source("modules/heatmaps.R")
source("modules/volcano_plots.R")
source("modules/ma_plots.R")
source("modules/venn_diagrams.R")
source("modules/summary_reports.R")
source("modules/diagnostic_module.R")

# =============================================================================
# ANALYSIS PARAMETERS - MODIFY THESE AS NEEDED
# =============================================================================

# Set working directory - MODIFY THIS PATH
setwd("path/to/your/directory")

# Analysis parameters
params <- list(
  # Statistical thresholds
  STAT_THRESHOLD_TYPE = "pvalue",  # "FDR" or "pvalue"
  FDR_THRESHOLD = 0.05,
  PVALUE_THRESHOLD = 0.05,
  LOGFC_THRESHOLD = 1,
  
  # Filtering parameters
  MIN_CPM = 5,
  MIN_SAMPLES = 3,
  
  # Visualization parameters
  TOP_N_HEATMAP = 500,
  TOP_N_VOLCANO_LABELS = 10
)

# Set the threshold based on choice
if (params$STAT_THRESHOLD_TYPE == "FDR") {
  params$stat_threshold <- params$FDR_THRESHOLD
  params$stat_column <- "FDR"
  params$stat_label <- "FDR"
} else {
  params$stat_threshold <- params$PVALUE_THRESHOLD
  params$stat_column <- "PValue"
  params$stat_label <- "P-value"
}

# =============================================================================
# MAIN ANALYSIS WORKFLOW
# =============================================================================

cat("Starting miRNA Differential Expression Analysis...\n")
cat(paste("Using", params$stat_label, "threshold of", params$stat_threshold, 
          "and logFC threshold of", params$LOGFC_THRESHOLD, "\n"))

# 1. Setup and load packages
cat("1. Loading required packages...\n")
load_required_packages()

# 2. Load and validate data
cat("2. Loading and validating data...\n")
data_files <- list(
  sample_info = "sample_info.csv",  # Updated to use your simplified file
  count_matrix = "my_merged_results_reorganized.tsv"
)

raw_data <- load_and_validate_data(data_files)

# 3. Create DGE object and filter
cat("3. Creating DGE object and filtering low-expressed genes...\n")
dge_data <- create_and_filter_dge(raw_data$count_matrix, raw_data$sample_info, params)

# 4. Normalization and dispersion estimation
cat("4. Normalizing data and estimating dispersions...\n")
normalized_data <- normalize_and_estimate_dispersions(dge_data$dge_filtered, dge_data$design)

# Save normalized counts
write.table(normalized_data$norm_counts, 
            file = "normalized_counts_CPM.txt",
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# 5. Create BCV plots
cat("5. Creating BCV plots...\n")
create_bcv_plots(normalized_data$dge_filtered, "BCV_plots_comprehensive.png")

# 6. Differential expression analysis
cat("6. Performing differential expression analysis...\n")

# DEBUG: Check what group levels are actually available
cat("\n=== DEBUGGING GROUP LEVELS ===\n")
cat("Available group levels in design matrix:\n")
print(colnames(normalized_data$design))
cat("\nGroup factor levels:\n")
print(levels(normalized_data$dge_filtered$samples$group))
cat("\nDesign matrix structure:\n")
print(normalized_data$design)
cat("===============================\n\n")

# Get the actual group names from the design matrix
available_groups <- colnames(normalized_data$design)
cat("Available groups for contrasts:", paste(available_groups, collapse = ", "), "\n")

# Create contrast matrix based on available groups
cat("Available groups:", paste(available_groups, collapse = ", "), "\n")

# Check which contrasts can be created with available groups
possible_contrasts <- list()

# Check Young Leaves vs Mature Leaves
if (all(c("Young_Leaves", "Mature_Leaves") %in% available_groups)) {
  possible_contrasts[["Young_Leaves_vs_Mature_Leaves"]] <- "Young_Leaves - Mature_Leaves"
  cat("??? Can create: Young_Leaves vs Mature_Leaves contrast\n")
} else {
  cat("??? Cannot create Young_Leaves vs Mature_Leaves contrast (missing:", 
      paste(setdiff(c("Young_Leaves", "Mature_Leaves"), available_groups), collapse = ", "), ")\n")
}

# Check Young Roots vs Mature Roots
if (all(c("Young_Roots", "Mature_Roots") %in% available_groups)) {
  possible_contrasts[["Young_Roots_vs_Mature_Roots"]] <- "Young_Roots - Mature_Roots"
  cat("??? Can create: Young_Roots vs Mature_Roots contrast\n")
} else {
  cat("??? Cannot create Young_Roots vs Mature_Roots contrast (missing:", 
      paste(setdiff(c("Young_Roots", "Mature_Roots"), available_groups), collapse = ", "), ")\n")
}

# Check Fifth Internode vs Third Internode
if (all(c("Fifth_Internode", "Third_Internode") %in% available_groups)) {
  possible_contrasts[["Fifth_Internode_vs_Third_Internode"]] <- "Fifth_Internode - Third_Internode"
  cat("??? Can create: Fifth_Internode vs Third_Internode contrast\n")
} else {
  cat("??? Cannot create Fifth_Internode vs Third_Internode contrast (missing:", 
      paste(setdiff(c("Fifth_Internode", "Third_Internode"), available_groups), collapse = ", "), ")\n")
}

# Stop if no contrasts can be created
if (length(possible_contrasts) == 0) {
  stop("ERROR: None of the requested contrasts can be created with available groups.")
}

cat("\nCreating", length(possible_contrasts), "possible contrasts...\n")

# Create contrast matrix with only possible contrasts
if ("Young_Leaves_vs_Mature_Leaves" %in% names(possible_contrasts) && 
    "Young_Roots_vs_Mature_Roots" %in% names(possible_contrasts) && 
    "Fifth_Internode_vs_Third_Internode" %in% names(possible_contrasts)) {
  # All three contrasts possible
  contrast_matrix <- makeContrasts(
    Young_Leaves_vs_Mature_Leaves = Young_Leaves - Mature_Leaves,        # Young Leaves vs Mature Leaves
    Young_Roots_vs_Mature_Roots = Young_Roots - Mature_Roots,        # Young Roots vs Mature Roots
    Fifth_Internode_vs_Third_Internode = Fifth_Internode - Third_Internode,        # Fifth Internode vs Third Internode
    levels = normalized_data$design
  )
} else if ("Young_Roots_vs_Mature_Roots" %in% names(possible_contrasts) && 
           "Fifth_Internode_vs_Third_Internode" %in% names(possible_contrasts)) {
  # Only Young_Roots_vs_Mature_Roots and Fifth_Internode_vs_Third_Internode possible (missing Young_Leaves)
  contrast_matrix <- makeContrasts(
    Young_Roots_vs_Mature_Roots = Young_Roots - Mature_Roots,        # Young Roots vs Mature Roots
    Fifth_Internode_vs_Third_Internode = Fifth_Internode - Third_Internode,        # Fifth Internode vs Third Internode
    levels = normalized_data$design
  )
} else if ("Young_Leaves_vs_Mature_Leaves" %in% names(possible_contrasts) && 
           "Fifth_Internode_vs_Third_Internode" %in% names(possible_contrasts)) {
  # Only Young_Leaves_vs_Mature_Leaves and Fifth_Internode_vs_Third_Internode possible
  contrast_matrix <- makeContrasts(
    Young_Leaves_vs_Mature_Leaves = Young_Leaves - Mature_Leaves,        # Young Leaves vs Mature Leaves
    Fifth_Internode_vs_Third_Internode = Fifth_Internode - Third_Internode,        # Fifth Internode vs Third Internode
    levels = normalized_data$design
  )
} else if ("Young_Leaves_vs_Mature_Leaves" %in% names(possible_contrasts) && 
           "Young_Roots_vs_Mature_Roots" %in% names(possible_contrasts)) {
  # Only Young_Leaves_vs_Mature_Leaves and Young_Roots_vs_Mature_Roots possible
  contrast_matrix <- makeContrasts(
    Young_Leaves_vs_Mature_Leaves = Young_Leaves - Mature_Leaves,        # Young Leaves vs Mature Leaves
    Young_Roots_vs_Mature_Roots = Young_Roots - Mature_Roots,        # Young Roots vs Mature Roots
    levels = normalized_data$design
  )
} else if ("Young_Roots_vs_Mature_Roots" %in% names(possible_contrasts)) {
  # Only Young_Roots_vs_Mature_Roots possible
  contrast_matrix <- makeContrasts(
    Young_Roots_vs_Mature_Roots = Young_Roots - Mature_Roots,        # Young Roots vs Mature Roots
    levels = normalized_data$design
  )
} else if ("Fifth_Internode_vs_Third_Internode" %in% names(possible_contrasts)) {
  # Only Fifth_Internode_vs_Third_Internode possible
  contrast_matrix <- makeContrasts(
    Fifth_Internode_vs_Third_Internode = Fifth_Internode - Third_Internode,        # Fifth Internode vs Third Internode
    levels = normalized_data$design
  )
} else if ("Young_Leaves_vs_Mature_Leaves" %in% names(possible_contrasts)) {
  # Only Young_Leaves_vs_Mature_Leaves possible
  contrast_matrix <- makeContrasts(
    Young_Leaves_vs_Mature_Leaves = Young_Leaves - Mature_Leaves,        # Young Leaves vs Mature Leaves
    levels = normalized_data$design
  )
}

cat("Created contrasts:", paste(colnames(contrast_matrix), collapse = ", "), "\n")

# Print the contrast matrix for verification
cat("\nContrast matrix:\n")
print(contrast_matrix)

de_results <- perform_comprehensive_DE_analysis(normalized_data$fit, contrast_matrix, params)

# 7. Create visualizations (with error handling)
cat("7. Creating visualizations...\n")

# Global expression heatmap (with error handling)
cat("   - Creating global expression heatmap...\n")
tryCatch({
  global_heatmap <- create_global_expression_heatmap_fixed(
    normalized_data$dge_filtered, 
    top_n = params$TOP_N_HEATMAP
  )
  cat("   ??? Global heatmap created successfully\n")
}, error = function(e) {
  cat("   ??? Error creating global heatmap:", e$message, "\n")
  cat("   - Skipping global heatmap\n")
})

# DEM-specific heatmaps (with error handling)
cat("   - Creating DEM-specific heatmaps...\n")
tryCatch({
  dem_heatmaps <- create_dem_heatmaps_fixed(
    normalized_data$dge_filtered, 
    de_results$results, 
    params$stat_column, 
    params$stat_threshold, 
    params$LOGFC_THRESHOLD
  )
  cat("   ??? DEM heatmaps created successfully\n")
}, error = function(e) {
  cat("   ??? Error creating DEM heatmaps:", e$message, "\n")
  cat("   - Skipping DEM heatmaps\n")
  dem_heatmaps <- NULL
})

# Volcano plots (with error handling)
cat("   - Creating volcano plots...\n")
tryCatch({
  volcano_results <- create_all_volcano_plots(de_results$results, params)
  cat("   ??? Volcano plots created successfully\n")
}, error = function(e) {
  cat("   ??? Error creating volcano plots:", e$message, "\n")
  cat("   - Skipping volcano plots\n")
  volcano_results <- NULL
})

# MA plots (with error handling)
cat("   - Creating MA plots...\n")
tryCatch({
  create_ma_plots(normalized_data$fit, contrast_matrix, params)
  cat("   ??? MA plots created successfully\n")
}, error = function(e) {
  cat("   ??? Error creating MA plots:", e$message, "\n")
  cat("   - Skipping MA plots\n")
})

# Venn diagram (with error handling)
cat("   - Creating Venn diagram...\n")
tryCatch({
  venn_results <- create_venn_diagram(de_results$results, params)
  cat("   ??? Venn diagram created successfully\n")
}, error = function(e) {
  cat("   ??? Error creating Venn diagram:", e$message, "\n")
  cat("   - Skipping Venn diagram\n")
  venn_results <- NULL
})

# 8. Generate summary report (with error handling)
cat("8. Generating summary report...\n")

tryCatch({
  # Generate comprehensive summary
  generate_comprehensive_summary(
    raw_data = raw_data,
    filtered_data = normalized_data,
    de_results = de_results,
    dem_heatmaps = dem_heatmaps,
    venn_results = venn_results,
    params = params
  )
  cat("   ??? Summary report generated successfully\n")
}, error = function(e) {
  cat("   ??? Error generating summary report:", e$message, "\n")
  cat("   - Creating basic summary instead\n")
  
  # Create a basic summary table
  basic_summary <- create_de_summary_table(de_results$results, params)
  write.csv(basic_summary, "basic_DE_summary.csv", row.names = FALSE)
  cat("   ??? Basic summary saved as 'basic_DE_summary.csv'\n")
})

cat("\n", rep("=", 80), "\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("Check the generated files in your working directory.\n")
cat("Available groups in your dataset:", paste(levels(normalized_data$dge_filtered$samples$group), collapse = ", "), "\n")
cat(rep("=", 80), "\n")
