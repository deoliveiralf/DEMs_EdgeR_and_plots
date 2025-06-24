# =============================================================================
# HEATMAP FUNCTIONS MODULE - CLEANED VERSION (Only ComplexHeatmap outputs)
# =============================================================================
# This module contains all heatmap-related functions for miRNA analysis

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

# Function to calculate z-scores
calculate_zscore <- function(logcpm_data) {
  # Calculate z-scores: (x - mean) / sd for each gene (row)
  zscore_data <- t(apply(logcpm_data, 1, function(x) {
    if (sd(x) == 0) {
      rep(0, length(x))  # Handle constant values
    } else {
      (x - mean(x)) / sd(x)
    }
  }))
  colnames(zscore_data) <- colnames(logcpm_data)
  return(zscore_data)
}

# Function to identify and annotate gene clusters
identify_gene_clusters <- function(logcpm_data, n_clusters = 4) {
  # Perform hierarchical clustering on genes (rows)
  gene_dist <- dist(logcpm_data, method = "euclidean")
  gene_hclust <- hclust(gene_dist, method = "complete")
  
  # Cut tree to get clusters
  gene_clusters <- cutree(gene_hclust, k = n_clusters)
  
  # Create cluster annotation data frame
  cluster_annotation <- data.frame(
    Cluster = paste("Cluster", gene_clusters),
    row.names = names(gene_clusters)
  )
  
  # Create colors for clusters
  cluster_colors <- create_safe_color_palette(cluster_annotation$Cluster, "Set1")
  
  return(list(
    annotation = cluster_annotation,
    colors = list(Cluster = cluster_colors),
    clusters = gene_clusters,
    dendrogram = gene_hclust
  ))
}

# Function to create global expression heatmap for all miRNAs - ENHANCED VERSION
create_global_expression_heatmap_fixed <- function(dge_obj, top_n = 500, 
                                                   filename_scaled = "global_expression_heatmap_scaled.png",
                                                   filename_zscore = "global_expression_heatmap_zscore.png",
                                                   n_clusters = 4) {
  
  # Get normalized log-CPM values
  logcpm <- cpm(dge_obj, log = TRUE, normalized.lib.sizes = TRUE)
  
  # Select most variable genes if dataset is large
  if (nrow(logcpm) > top_n) {
    var_genes <- apply(logcpm, 1, var)
    top_var_genes <- names(sort(var_genes, decreasing = TRUE))[1:top_n]
    logcpm_subset <- logcpm[top_var_genes, ]
    cat(paste("Using top", top_n, "most variable miRNAs for global heatmap\n"))
  } else {
    logcpm_subset <- logcpm
    cat(paste("Using all", nrow(logcpm), "miRNAs for global heatmap\n"))
  }
  
  # Calculate z-scores
  zscore_subset <- calculate_zscore(logcpm_subset)
  cat("Z-scores calculated for expression data\n")
  
  # Identify gene clusters for both scaled and z-score data
  cluster_info_scaled <- identify_gene_clusters(t(scale(t(logcpm_subset))), n_clusters)
  cluster_info_zscore <- identify_gene_clusters(zscore_subset, n_clusters)
  
  cat(paste("Identified", n_clusters, "gene clusters\n"))
  
  # Print cluster summaries
  cat("Gene cluster distribution (scaled data):\n")
  print(table(cluster_info_scaled$clusters))
  cat("Gene cluster distribution (z-score data):\n")
  print(table(cluster_info_zscore$clusters))
  
  # Prepare sample annotation data - ONLY GROUP
  annotation_col <- data.frame(
    Group = dge_obj$samples$group,
    row.names = colnames(logcpm_subset)
  )
  
  # Create safe color palettes for samples - ONLY GROUP
  group_colors <- create_safe_color_palette(dge_obj$samples$group, "Set3")
  
  sample_annotation_colors <- list(
    Group = group_colors
  )
  
  # ===============================
  # CREATE SCALED HEATMAP
  # ===============================
  
  png(filename_scaled, width = 4000, height = 3000, res = 300)
  
  tryCatch({
    # Prepare annotations for ComplexHeatmap - samples and gene clusters
    col_ha <- HeatmapAnnotation(
      Group = annotation_col$Group,
      col = sample_annotation_colors,
      annotation_name_gp = gpar(fontsize = 10),
      annotation_legend_param = list(
        Group = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8))
      )
    )
    
    row_ha <- rowAnnotation(
      Cluster = cluster_info_scaled$annotation$Cluster,
      col = cluster_info_scaled$colors,
      annotation_name_gp = gpar(fontsize = 10),
      annotation_legend_param = list(
        Cluster = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8))
      )
    )
    
    # Scale the data
    logcpm_scaled <- t(scale(t(logcpm_subset)))
    
    ht <- Heatmap(
      logcpm_scaled,
      name = "Scaled\nlog-CPM",
      col = viridis::viridis(100),
      top_annotation = col_ha,
      left_annotation = row_ha,
      clustering_distance_rows = "pearson",
      clustering_distance_columns = "pearson",
      clustering_method_rows = "complete",
      clustering_method_columns = "complete",
      show_row_names = ifelse(nrow(logcpm_subset) <= 50, TRUE, FALSE),
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 8),
      column_title = paste("Global miRNA Expression Heatmap - Scaled (Top", nrow(logcpm_subset), "most variable)"),
      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
      rect_gp = gpar(col = "white", lwd = 0.5),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8)
      )
    )
    
    draw(ht)
  }, error = function(e) {
    cat(paste("Error creating scaled ComplexHeatmap:", e$message, "\n"))
  })
  
  dev.off()
  
  cat(paste("Scaled global expression heatmap saved as:", filename_scaled, "\n"))
  
  # ===============================
  # CREATE Z-SCORE HEATMAP
  # ===============================
  
  png(filename_zscore, width = 4000, height = 3000, res = 300)
  
  tryCatch({
    # Prepare annotations for ComplexHeatmap - samples and gene clusters for z-score
    col_ha_z <- HeatmapAnnotation(
      Group = annotation_col$Group,
      col = sample_annotation_colors,
      annotation_name_gp = gpar(fontsize = 10),
      annotation_legend_param = list(
        Group = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8))
      )
    )
    
    row_ha_z <- rowAnnotation(
      Cluster = cluster_info_zscore$annotation$Cluster,
      col = cluster_info_zscore$colors,
      annotation_name_gp = gpar(fontsize = 10),
      annotation_legend_param = list(
        Cluster = list(title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 8))
      )
    )
    
    ht_z <- Heatmap(
      zscore_subset,
      name = "Z-Score",
      col = colorRampPalette(c("blue", "white", "red"))(100),
      top_annotation = col_ha_z,
      left_annotation = row_ha_z,
      clustering_distance_rows = "pearson",
      clustering_distance_columns = "pearson",
      clustering_method_rows = "complete",
      clustering_method_columns = "complete",
      show_row_names = ifelse(nrow(logcpm_subset) <= 50, TRUE, FALSE),
      show_column_names = TRUE,
      row_names_gp = gpar(fontsize = 6),
      column_names_gp = gpar(fontsize = 8),
      column_title = paste("Global miRNA Expression Heatmap - Z-Score (Top", nrow(logcpm_subset), "most variable)"),
      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
      rect_gp = gpar(col = "white", lwd = 0.5),
      heatmap_legend_param = list(
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8)
      )
    )
    
    draw(ht_z)
  }, error = function(e) {
    cat(paste("Error creating z-score ComplexHeatmap:", e$message, "\n"))
  })
  
  dev.off()
  
  cat(paste("Z-score global expression heatmap saved as:", filename_zscore, "\n"))
  
  return(list(
    data_scaled = logcpm_scaled, 
    data_zscore = zscore_subset,
    annotation = annotation_col,
    clusters_scaled = cluster_info_scaled,
    clusters_zscore = cluster_info_zscore
  ))
}

# FIXED Function to extract group names from comparison string
extract_comparison_groups <- function(comparison_name) {
  cat(paste("Parsing comparison:", comparison_name, "\n"))
  
  # Handle different comparison naming patterns
  # Examples: "Group_F_J_vs_Group_F_M", "Group_RJ_vs_RM", "FJ_vs_FM"
  
  # Split by "_vs_" first
  parts <- strsplit(comparison_name, "_vs_")[[1]]
  
  if (length(parts) != 2) {
    stop(paste("Cannot parse comparison name:", comparison_name, "- should contain exactly one '_vs_'"))
  }
  
  group1 <- parts[1]
  group2 <- parts[2]
  
  # Remove "Group_" prefix if present in both parts
  if (grepl("^Group_", group1)) {
    group1 <- gsub("^Group_", "", group1)
  }
  if (grepl("^Group_", group2)) {
    group2 <- gsub("^Group_", "", group2)
  }
  
  cat(paste("Extracted groups:", group1, "and", group2, "\n"))
  
  return(c(group1, group2))
}

# CLEANED Function to create DEM-specific heatmaps for each comparison (ComplexHeatmap only)
create_dem_heatmaps_fixed <- function(dge_obj, de_results, stat_column = "FDR", 
                                      stat_threshold = 0.05, logFC_threshold = 1,
                                      debug = TRUE) {
  
  # Get normalized log-CPM values
  logcpm <- cpm(dge_obj, log = TRUE, normalized.lib.sizes = TRUE)
  
  heatmap_summaries <- list()
  
  # Print available groups for debugging
  if (debug) {
    cat("Available groups in dge_obj$samples$group:", unique(dge_obj$samples$group), "\n")
    cat("Sample group assignments:\n")
    print(table(dge_obj$samples$group))
  }
  
  for (comparison in names(de_results)) {
    cat(paste("\n=== Processing comparison:", comparison, "===\n"))
    
    # Extract the two groups being compared
    tryCatch({
      comparison_groups <- extract_comparison_groups(comparison)
      group1 <- comparison_groups[1]
      group2 <- comparison_groups[2]
      
      cat(paste("Looking for samples with groups:", group1, "and", group2, "\n"))
    }, error = function(e) {
      cat(paste("Error extracting groups from comparison", comparison, ":", e$message, "\n"))
      cat("Skipping this comparison\n")
      next
    })
    
    # Check if the extracted groups exist in the data
    available_groups <- unique(dge_obj$samples$group)
    groups_found <- c(group1, group2) %in% available_groups
    
    if (!all(groups_found)) {
      cat(paste("Groups not found in data. Available groups:", paste(available_groups, collapse = ", "), "\n"))
      cat(paste("Trying alternative group matching...\n"))
      
      # Try to find groups that contain the extracted group names
      potential_group1 <- available_groups[grepl(group1, available_groups, fixed = TRUE)]
      potential_group2 <- available_groups[grepl(group2, available_groups, fixed = TRUE)]
      
      if (length(potential_group1) == 1 && length(potential_group2) == 1) {
        group1 <- potential_group1
        group2 <- potential_group2
        cat(paste("Found matching groups:", group1, "and", group2, "\n"))
      } else {
        cat("Could not find matching groups - skipping this comparison\n")
        heatmap_summaries[[comparison]] <- list(n_dems = 0, dem_genes = character(0), error = "Groups not found")
        next
      }
    }
    
    # Filter samples to include only those from the two groups being compared
    samples_in_comparison <- dge_obj$samples$group %in% c(group1, group2)
    
    if (sum(samples_in_comparison) == 0) {
      cat(paste("No samples found for groups", group1, "and", group2, "- skipping\n"))
      heatmap_summaries[[comparison]] <- list(n_dems = 0, dem_genes = character(0), error = "No samples found")
      next
    }
    
    # Subset expression data and sample information
    logcpm_subset <- logcpm[, samples_in_comparison]
    samples_subset <- dge_obj$samples[samples_in_comparison, ]
    
    cat(paste("Using", sum(samples_in_comparison), "samples for comparison:\n"))
    print(table(samples_subset$group))
    
    # Prepare sample annotations for the subset - ONLY GROUP
    annotation_col <- data.frame(
      Group = samples_subset$group,
      row.names = colnames(logcpm_subset)
    )
    
    # Create safe color palettes for the subset - ONLY GROUP
    group_colors <- create_safe_color_palette(samples_subset$group, "Set3")
    
    annotation_colors <- list(
      Group = group_colors
    )
    
    # Get DEMs for this comparison
    results <- de_results[[comparison]]
    
    # Debug: show some statistics
    if (debug) {
      cat(paste("Total genes in results:", nrow(results), "\n"))
      cat(paste("Genes with |logFC| >", logFC_threshold, ":", sum(abs(results$logFC) > logFC_threshold, na.rm = TRUE), "\n"))
      cat(paste("Genes with", stat_column, "<", stat_threshold, ":", sum(results[[stat_column]] < stat_threshold, na.rm = TRUE), "\n"))
    }
    
    # Apply thresholds
    dem_genes <- rownames(results[abs(results$logFC) > logFC_threshold & 
                                    results[[stat_column]] < stat_threshold, ])
    
    if (length(dem_genes) == 0) {
      cat(paste("No DEMs found for", comparison, "with current thresholds - skipping heatmap\n"))
      heatmap_summaries[[comparison]] <- list(n_dems = 0, dem_genes = character(0))
      next
    }
    
    cat(paste("Found", length(dem_genes), "DEMs for", comparison, "\n"))
    
    if (length(dem_genes) == 1) {
      cat(paste("Only 1 DEM found for", comparison, "- creating simple plot\n"))
      # Create a simple barplot for single gene using only comparison samples
      single_gene_data <- logcpm_subset[dem_genes, ]
      filename_single <- paste0("single_DEM_", gsub("[^A-Za-z0-9]", "_", comparison), "_", stat_column, "_", stat_threshold, ".png")
      
      png(filename_single, width = 2400, height = 1800, res = 300)
      par(mar = c(8, 4, 4, 2))
      barplot(single_gene_data, 
              main = paste("Expression of", dem_genes, "\nin", group1, "vs", group2),
              ylab = "log-CPM", las = 2, 
              col = group_colors[samples_subset$group])
      legend("topright", legend = names(group_colors), fill = group_colors, cex = 0.8)
      dev.off()
      
      heatmap_summaries[[comparison]] <- list(n_dems = 1, dem_genes = dem_genes)
      next
    }
    
    # Extract expression data for DEMs using only comparison samples
    dem_logcpm <- logcpm_subset[dem_genes, , drop = FALSE]
    
    # Identify gene clusters for DEM heatmap
    dem_cluster_info <- identify_gene_clusters(t(scale(t(dem_logcpm))), min(4, max(2, length(dem_genes) %/% 5)))
    
    cat(paste("Identified", length(unique(dem_cluster_info$clusters)), "gene clusters for DEMs\n"))
    cat("DEM cluster distribution:\n")
    print(table(dem_cluster_info$clusters))
    
    # Clean up comparison name for title
    clean_comparison <- paste(group1, "vs", group2)
    
    # Create filename - only ComplexHeatmap version
    filename <- paste0("DEM_heatmap_", gsub("[^A-Za-z0-9]", "_", comparison), "_", stat_column, "_", stat_threshold, "_pairwise.png")
    
    # Determine heatmap dimensions based on number of DEMs
    if (length(dem_genes) <= 20) {
      height <- 2400
      fontsize_row <- 8
    } else if (length(dem_genes) <= 50) {
      height <- 3000
      fontsize_row <- 7
    } else {
      height <- 4000
      fontsize_row <- 6
    }
    
    # Adjust width based on number of samples
    width <- max(2000, 200 * ncol(dem_logcpm) + 1000)
    
    # Create ComplexHeatmap version only
    png(filename, width = width, height = height, res = 300)
    
    tryCatch({
      # Prepare annotations for ComplexHeatmap - samples and gene clusters
      col_ha <- HeatmapAnnotation(
        Group = annotation_col$Group,
        col = annotation_colors,
        annotation_name_gp = gpar(fontsize = 8)
      )
      
      row_ha <- rowAnnotation(
        Cluster = dem_cluster_info$annotation$Cluster,
        col = dem_cluster_info$colors,
        annotation_name_gp = gpar(fontsize = 8),
        annotation_legend_param = list(
          Cluster = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
        )
      )
      
      # Scale the data
      dem_logcpm_scaled <- t(scale(t(dem_logcpm)))
      
      ht <- Heatmap(
        dem_logcpm_scaled,
        name = "Scaled\nlog-CPM",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        top_annotation = col_ha,
        left_annotation = row_ha,
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson",
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = fontsize_row),
        column_names_gp = gpar(fontsize = 8),
        column_title = paste("DEMs in", clean_comparison, "(", length(dem_genes), "DEMs,", ncol(dem_logcpm), "samples)"),
        column_title_gp = gpar(fontsize = 11, fontface = "bold"),
        rect_gp = gpar(col = "white", lwd = 0.5)
      )
      
      draw(ht)
    }, error = function(e) {
      cat(paste("Error creating ComplexHeatmap for", comparison, ":", e$message, "\n"))
    })
    
    dev.off()
    
    cat(paste("DEM heatmap saved as:", filename, "\n"))
    
    # Store summary - simplified without regulation counts
    regulation_info <- results[dem_genes, "regulation"]
    heatmap_summaries[[comparison]] <- list(
      n_dems = length(dem_genes),
      dem_genes = dem_genes,
      n_up = sum(regulation_info == "Up-regulated"),
      n_down = sum(regulation_info == "Down-regulated"),
      n_samples = ncol(dem_logcpm),
      groups_compared = paste(group1, "vs", group2),
      group1 = group1,
      group2 = group2
    )
  }
  
  # Create summary table
  summary_df <- data.frame(
    Comparison = names(heatmap_summaries),
    Groups_Compared = sapply(heatmap_summaries, function(x) ifelse(is.null(x$groups_compared), "", x$groups_compared)),
    Samples_Used = sapply(heatmap_summaries, function(x) ifelse(is.null(x$n_samples), 0, x$n_samples)),
    Total_DEMs = sapply(heatmap_summaries, function(x) x$n_dems),
    Up_regulated = sapply(heatmap_summaries, function(x) ifelse(is.null(x$n_up), 0, x$n_up)),
    Down_regulated = sapply(heatmap_summaries, function(x) ifelse(is.null(x$n_down), 0, x$n_down)),
    stringsAsFactors = FALSE
  )
  
  write.csv(summary_df, paste0("DEM_heatmap_summary_pairwise_", stat_column, "_", stat_threshold, ".csv"), row.names = FALSE)
  
  cat("\n=== PAIRWISE DEM HEATMAP SUMMARY ===\n")
  print(summary_df)
  
  return(heatmap_summaries)
}

# Helper function to diagnose group matching issues
diagnose_group_matching <- function(dge_obj, comparison_names) {
  cat("=== GROUP MATCHING DIAGNOSIS ===\n")
  cat("Available groups in dge_obj:\n")
  print(unique(dge_obj$samples$group))
  cat("\nSample distribution:\n")
  print(table(dge_obj$samples$group))
  
  cat("\nComparison names to parse:\n")
  for (comp in comparison_names) {
    cat(paste("Comparison:", comp, "\n"))
    tryCatch({
      groups <- extract_comparison_groups(comp)
      cat(paste("  Extracted groups:", groups[1], "and", groups[2], "\n"))
      
      # Check if groups exist
      group_exists <- groups %in% unique(dge_obj$samples$group)
      cat(paste("  Group", groups[1], "exists:", group_exists[1], "\n"))
      cat(paste("  Group", groups[2], "exists:", group_exists[2], "\n"))
      
    }, error = function(e) {
      cat(paste("  Error:", e$message, "\n"))
    })
    cat("\n")
  }
}