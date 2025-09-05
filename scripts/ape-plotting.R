#!/usr/bin/env Rscript
# Phylogenetic tree plotting script using ape package
# This script reads IQ-TREE output and creates publication-ready plots

# Load required libraries
library(ape)
library(phytools)
library(stringr)
library(data.table)
library(dplyr)

# Configuration - EDIT THESE PATHS AS NEEDED
# ==========================================

# Output directory for plots
output_dir <- "plots/"

# IQ-TREE results (treefile format)  
tree_files <- c(
  "iqtree_results/vgsc_focal.fasta.treefile",
  "iqtree_results/vgsc_upstream.fasta.treefile", 
  "iqtree_results/vgsc_downstream.fasta.treefile"
)

# Corresponding metadata files
metadata_files <- c(
  "fastas/vgsc_focal.metadata.tsv",
  "fastas/vgsc_upstream.metadata.tsv",
  "fastas/vgsc_downstream.metadata.tsv" 
)

# Plot titles
plot_titles <- c("Focal haplotype", "Upstream", "Downstream")

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to process a single tree and metadata file
process_tree <- function(tree_file, metadata_file, plot_title) {
  
  # Check if files exist
  if (!file.exists(tree_file)) {
    warning(paste("Tree file not found:", tree_file))
    return(NULL)
  }
  
  if (!file.exists(metadata_file)) {
    warning(paste("Metadata file not found:", metadata_file))
    return(NULL)
  }
  
  # Read tree and metadata
  cat(paste("Processing:", tree_file, "\n"))
  phy <- read.tree(tree_file)
  meta <- fread(metadata_file) %>% as.data.frame()
  
  # Apply midpoint rooting
  phy_p <- midpoint.root(phy)
  
  # Match metadata to tree tips
  # Handle both 'hap' and 'sample_id' column names for compatibility
  tip_col <- if("hap" %in% colnames(meta)) "hap" else "sample_id"
  
  if (!tip_col %in% colnames(meta)) {
    warning("Could not find tip identifier column ('hap' or 'sample_id') in metadata")
    return(NULL)
  }
  
  meta <- meta %>% arrange(factor(.data[[tip_col]], levels = phy_p$tip.label))
  
  # Set branch length constraints for visualization
  phy_p$edge.length[phy_p$edge.length > 0.005] <- 0.005
  phy_p$edge.length[phy_p$edge.length == 0] <- 5e-5
  
  # Add metadata to tree object for easy access
  phy_p$metadata <- meta
  phy_p$plot_title <- plot_title
  
  return(phy_p)
}


# Function to assign colors based on metadata column
assign_tip_colors <- function(meta, color_var) {
  
  # Check if the column exists
  if (!color_var %in% colnames(meta)) {
    warning(paste("Column", color_var, "not found in metadata. Using default colors."))
    return(rep("grey", nrow(meta)))
  }
  
  if (color_var == 'taxon' || color_var == 'aim_species') {
    # Color by taxonomic groups
    tip_colors <- case_when(
      meta[[color_var]] == 'gambiae' ~ 'indianred',
      meta[[color_var]] == 'coluzzii' ~ 'dodgerblue', 
      meta[[color_var]] == 'arabiensis' ~ 'aquamarine',
      meta[[color_var]] == 'melas' ~ 'cornsilk',
      meta[[color_var]] == 'merus' ~ 'cornsilk4',
      meta[[color_var]] == 'quadriannulatus' ~ 'darkolivegreen',
      meta[[color_var]] == 'fontenillei' ~ 'orange',
      TRUE ~ 'grey'
    )
  } else if (color_var == 'karyotype') {
    # Color by karyotype
    tip_colors <- case_when(
      meta[[color_var]] == '2l+a' ~ 'bisque2',
      meta[[color_var]] == '2la' ~ 'dodgerblue',
      TRUE ~ 'grey'
    )
  } else if (color_var == 'country') {
    # Use rainbow colors for countries
    unique_vals <- unique(meta[[color_var]])
    color_palette <- rainbow(length(unique_vals))
    names(color_palette) <- unique_vals
    tip_colors <- color_palette[meta[[color_var]]]
  } else {
    # Generic coloring for other variables
    unique_vals <- unique(meta[[color_var]])
    if (length(unique_vals) <= 12) {
      color_palette <- c("indianred", "dodgerblue", "forestgreen", "orange", 
                        "purple", "brown", "pink", "grey", "yellow", 
                        "cyan", "magenta", "darkgreen")[1:length(unique_vals)]
    } else {
      color_palette <- rainbow(length(unique_vals))
    }
    names(color_palette) <- unique_vals
    tip_colors <- color_palette[meta[[color_var]]]
  }
  
  return(tip_colors)
}

# Function to plot a phylogenetic tree
plot_tree <- function(tree_obj, color_var = "taxon", tree_type = "unrooted") {
  
  if (is.null(tree_obj) || is.null(tree_obj$metadata)) {
    warning("Invalid tree object or missing metadata")
    return(NULL)
  }
  
  # Get tip colors
  tip_colors <- assign_tip_colors(tree_obj$metadata, color_var)
  
  # Plot tree based on type
  if (tree_type == "unrooted") {
    plot.phylo(tree_obj, type = "unrooted",
               use.edge.length = TRUE, 
               show.tip.label = FALSE,  # Hide tip labels for cleaner look
               show.node.label = FALSE,
               tip.color = tip_colors,
               edge.color = "slategray3",
               font = 1, 
               edge.width = 1.5,
               main = paste(tree_obj$plot_title, "-", color_var))
  } else if (tree_type == "phylogram") {
    plot.phylo(tree_obj, type = "phylogram",
               use.edge.length = TRUE,
               show.tip.label = TRUE,
               show.node.label = FALSE, 
               tip.color = tip_colors,
               edge.color = "slategray3",
               font = 0.8,
               edge.width = 1,
               main = paste(tree_obj$plot_title, "-", color_var))
  }
}

# MAIN EXECUTION
# ===============

cat("Starting phylogenetic tree plotting script...\n")

# Variables for coloring (edit as needed)
color_variables <- c("taxon", "country")  # Add more variables as needed

# Process all trees
processed_trees <- list()
for (i in 1:length(tree_files)) {
  tree_obj <- process_tree(tree_files[i], metadata_files[i], plot_titles[i])
  if (!is.null(tree_obj)) {
    processed_trees[[i]] <- tree_obj
  }
}

# Remove NULL entries
processed_trees <- processed_trees[!sapply(processed_trees, is.null)]

if (length(processed_trees) == 0) {
  stop("No valid tree files were processed. Check file paths and formats.")
}

cat(paste("Successfully processed", length(processed_trees), "trees.\n"))

# Create plots for each tree and color variable
for (color_var in color_variables) {
  
  cat(paste("\nCreating plots colored by:", color_var, "\n"))
  
  # Standard size plots
  pdf_file <- file.path(output_dir, paste0("phylo_trees_", color_var, ".pdf"))
  pdf(pdf_file, height = 8, width = 10)
  
  par(mfrow = c(ceiling(length(processed_trees)/2), 2))  # Arrange plots in grid
  
  for (tree_obj in processed_trees) {
    plot_tree(tree_obj, color_var = color_var, tree_type = "unrooted")
  }
  
  dev.off()
  cat(paste("Saved:", pdf_file, "\n"))
  
  # Large detailed plots (one per page)
  pdf_file_large <- file.path(output_dir, paste0("phylo_trees_", color_var, "_large.pdf"))
  pdf(pdf_file_large, height = 12, width = 16)
  
  for (tree_obj in processed_trees) {
    par(mfrow = c(1, 1))  # One plot per page
    plot_tree(tree_obj, color_var = color_var, tree_type = "phylogram")
  }
  
  dev.off()
  cat(paste("Saved:", pdf_file_large, "\n"))
}

cat("\n=== Plotting completed successfully! ===\n")
cat("Output files saved in:", output_dir, "\n")