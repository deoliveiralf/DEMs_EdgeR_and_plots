# =============================================================================
# VENN DIAGRAM MODULE
# =============================================================================
# This module handles Venn diagram creation for overlapping DEMs

# Function to create Venn diagram for overlapping DEMs
create_venn_diagram <- function(de_results, params) {  # Fixed: accepts params object
  
  # Extract parameters from params object
  stat_column <- params$stat_column
  stat_threshold <- params$stat_threshold
  logFC_threshold <- params$LOGFC_THRESHOLD
  
  # Extract significantly DE miRNAs for each comparison
  de_lists <- list()
  comparison_names <- names(de_results)
  
  for (comp in comparison_names) {
    results <- de_results[[comp]]
    significant_mirnas <- rownames(results[abs(results$logFC) > logFC_threshold & 
                                             results[[stat_column]] < stat_threshold, ])
    
    # Clean up comparison name for display
    clean_name <- gsub("Group_", "", comp)
    clean_name <- gsub("_vs_", " vs ", clean_name)
    de_lists[[clean_name]] <- significant_mirnas
  }
  
  # Only proceed if we have exactly 3 comparisons and at least one has DEMs
  if (length(de_lists) == 3 && any(sapply(de_lists, length) > 0)) {
    
    tryCatch({
      venn_plot <- venn.diagram(
        x = de_lists,
        category.names = names(de_lists),
        filename = NULL,
        
        # Output features
        height = 3000,
        width = 3000,
        resolution = 300,
        
        # Circles
        lwd = 2,
        lty = 'blank',
        fill = c('#E31A1C', '#1F78B4', '#33A02C'),
        
        # Numbers
        cex = 1.5,
        fontface = "bold",
        fontfamily = "sans",
        
        # Set names
        cat.cex = 1.2,
        cat.fontface = "bold",
        cat.default.pos = "outer",
        cat.pos = c(-27, 27, 135),
        cat.dist = c(0.055, 0.055, 0.085),
        cat.fontfamily = "sans",
        rotation = 1
      )
      
      # Save the plot
      png("venn_diagram_DEMs.png", width = 3000, height = 3000, res = 300)
      grid.draw(venn_plot)
      dev.off()
      
      # Create summary of overlaps
      all_de_mirnas <- unique(unlist(de_lists))
      overlap_summary <- data.frame(
        Category = c("Total unique DEMs", names(de_lists), "Common to all three"),
        Count = c(length(all_de_mirnas), sapply(de_lists, length), 
                  length(Reduce(intersect, de_lists))),
        stringsAsFactors = FALSE
      )
      
      # Calculate pairwise overlaps
      pairs <- combn(names(de_lists), 2, simplify = FALSE)
      for (pair in pairs) {
        overlap_count <- length(intersect(de_lists[[pair[1]]], de_lists[[pair[2]]]))
        overlap_summary <- rbind(overlap_summary, 
                                 data.frame(Category = paste(pair[1], "???", pair[2]), 
                                            Count = overlap_count,
                                            stringsAsFactors = FALSE))
      }
      
      write.csv(overlap_summary, "venn_diagram_summary.csv", row.names = FALSE)
      
      cat("\n=== VENN DIAGRAM SUMMARY ===\n")
      print(overlap_summary)
      
      return(list(venn_plot = venn_plot, summary = overlap_summary, de_lists = de_lists))
      
    }, error = function(e) {
      warning(paste("Could not create Venn diagram:", e$message))
      return(NULL)
    })
  } else {
    cat("Venn diagram not created: Need exactly 3 comparisons with at least one having DEMs\n")
    return(NULL)
  }
}