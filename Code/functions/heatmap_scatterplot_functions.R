# ============================================================================
# Heatmap Functions for Cell-Hub
# ============================================================================
# This file contains helper functions for generating and customizing heatmaps
# Author: Gaspard
# Last updated: 2025-11-10
# ============================================================================



#' Get Scaled Genes from Seurat Object
#' 
#' @param seurat_object Seurat object
#' @param assay Assay name
#' @return Character vector of gene names in scale.data
get_scaled_genes <- function(seurat_object, assay) {
  tryCatch({
    rownames(LayerData(seurat_object, assay = assay, layer = "scale.data"))
  }, error = function(e) {
    character(0)
  })
}

#' Apply Theme to Heatmap
#' 
#' @param plot ggplot object
#' @param angle Angle for X-axis labels
#' @param show_gene_names Show gene names on Y-axis
#' @return ggplot object with applied theme
apply_heatmap_theme <- function(plot, angle = 45, show_gene_names = TRUE) {
  
  # Build theme with X-axis rotation and Y-axis gene names
  theme_elements <- list()
  
  # X-axis (cluster labels) - rotated
  if (angle > 0) {
    theme_elements$axis.text.x <- element_text(angle = angle, hjust = 1)
  }
  
  # Y-axis (gene names) - make them visible
  if (show_gene_names) {
    theme_elements$axis.text.y <- element_text(size = 10, face = "bold", hjust = 1)
  }
  
  plot <- plot + do.call(theme, theme_elements)
  
  return(plot)
}



#' Validate and Parse Gene List from Text Input
#' 
#' @param gene_text Text input containing comma-separated gene names
#' @param seurat_object Seurat object to validate genes against
#' @param assay Assay to check genes in
#' @return List with valid_genes and missing_genes
parse_and_validate_genes <- function(gene_text, seurat_object, assay = "RNA") {
  
  # Parse gene text
  if (nchar(trimws(gene_text)) == 0) {
    stop("Please enter at least one gene name.")
  }
  
  selected_genes <- trimws(unlist(strsplit(gene_text, ",")))
  selected_genes <- selected_genes[nchar(selected_genes) > 0]
  
  if (length(selected_genes) == 0) {
    stop("No valid gene names provided.")
  }
  
  # Check which genes exist
  available_genes <- rownames(seurat_object[[assay]])
  valid_genes_idx <- selected_genes %in% available_genes
  
  valid_genes <- selected_genes[valid_genes_idx]
  missing_genes <- selected_genes[!valid_genes_idx]
  
  if (length(valid_genes) == 0) {
    stop(paste0("None of the genes found in ", assay, " assay. Missing: ",
                paste(missing_genes, collapse = ", ")))
  }
  
  return(list(
    valid_genes = valid_genes,
    missing_genes = missing_genes
  ))
}


#' Get top marker genes per cluster
#' 
#' @param seurat_object Seurat object with clustering
#' @param n_genes Number of top genes per cluster (default: 10)
#' @param assay Assay to use
#' @param min_pct Minimum percentage of cells expressing gene
#' @param logfc_threshold Log fold change threshold
#' @param clusters_to_include Optional vector of cluster names to restrict results
#' @return List with: genes (unique character vector), markers_table (data.frame per cluster)
get_top_markers <- function(seurat_object,
                            n_genes = 10,
                            assay = "RNA",
                            min_pct = 0.25,
                            logfc_threshold = 1,
                            clusters_to_include = NULL) {
  idents_vec <- as.character(Idents(seurat_object))
  valid_cells <- colnames(seurat_object)[!is.na(idents_vec)]
  if (length(valid_cells) < ncol(seurat_object)) {
    message("[get_top_markers] Dropping ", ncol(seurat_object) - length(valid_cells),
            " cells with NA idents")
    seurat_object <- subset(seurat_object, cells = valid_cells)
  }
  Idents(seurat_object) <- droplevels(Idents(seurat_object))
  # FindAllMarkers requires counts layer — integrated assay has none, use RNA
  de_assay <- if (assay == "integrated") "RNA" else assay
  DefaultAssay(seurat_object) <- de_assay
  seurat_object <- JoinLayers(seurat_object)
  # Ensure counts layer is dgCMatrix — FindAllMarkers requires sparse matrix format
  counts_layer <- LayerData(seurat_object, assay = de_assay, layer = "counts")
  if (!inherits(counts_layer, "dgCMatrix")) {
    LayerData(seurat_object, assay = de_assay, layer = "counts") <- as(counts_layer, "CsparseMatrix")
  }
  markers <- FindAllMarkers(seurat_object, min.pct = min_pct, assay = de_assay)
  if (is.null(markers) || nrow(markers) == 0) {
    stop("No DE genes identified. Try lowering min_pct or logfc_threshold.")
  }
  if (!is.null(clusters_to_include) && length(clusters_to_include) > 0) {
    markers <- markers %>%
      dplyr::filter(cluster %in% clusters_to_include)
    message(paste("Filtered markers to", length(unique(markers$cluster)),
                  "selected clusters:", paste(clusters_to_include, collapse = ", ")))
  }
  top_markers_df <- markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > logfc_threshold) %>%
    slice_head(n = n_genes) %>%
    ungroup()
  top_genes <- unique(top_markers_df$gene)
  return(list(
    genes = top_genes,
    markers_table = as.data.frame(top_markers_df)
  ))
}
#' Generate Average Expression Heatmap
#' 
#' @param seurat_object Seurat object
#' @param genes Character vector of gene names
#' @param assay Assay to use (default: "RNA")
#' @param clusters Optional vector of cluster names to include
#' @param group_by Column name to group by (default: NULL uses Idents)
#' @param scale_rows Scale expression values by row (gene) (default: TRUE)
#' @param color_palette Color palette to use (default: viridis)
#' @return ggplot object with average expression heatmap
generate_average_expression_heatmap <- function(seurat_object,
                                                genes,
                                                assay = "RNA",
                                                clusters = NULL,
                                                group_by = NULL,
                                                scale_rows = TRUE,
                                                color_palette = "viridis") {
  DefaultAssay(seurat_object) <- assay
  if (is.null(group_by)) {
    # Drop NA idents — Seurat v5 subsets can leave NA levels that AverageExpression
    # groups into a spurious "NA" column
    idents_vec <- as.character(Idents(seurat_object))
    valid_cells <- colnames(seurat_object)[!is.na(idents_vec)]
    if (length(valid_cells) < ncol(seurat_object)) {
      message("[Heatmap] Dropping ", ncol(seurat_object) - length(valid_cells),
              " cells with NA idents")
      seurat_object <- subset(seurat_object, cells = valid_cells)
    }
    Idents(seurat_object) <- droplevels(Idents(seurat_object))
    grouping_factor <- Idents(seurat_object)
    group_values <- levels(grouping_factor)
  } else {
    if (!group_by %in% colnames(seurat_object@meta.data)) {
      stop(paste("Column", group_by, "not found in metadata"))
    }
    grouping_factor <- seurat_object@meta.data[[group_by]]
    group_values <- unique(as.character(grouping_factor))
    group_values <- group_values[!is.na(group_values)]
  }
  if (!is.null(clusters)) {
    valid_groups <- clusters %in% group_values
    if (sum(valid_groups) == 0) stop("None of the specified groups is valid.")
    cells_to_keep <- colnames(seurat_object)[as.character(grouping_factor) %in% clusters]
    seurat_subset <- subset(seurat_object, cells = cells_to_keep)
    group_order <- clusters
  } else {
    seurat_subset <- seurat_object
    group_order <- if (is.null(group_by)) levels(Idents(seurat_subset)) else sort(unique(as.character(seurat_subset@meta.data[[group_by]])))
    group_order <- group_order[!is.na(group_order)]
  }
  if (is.null(group_by)) {
    avg_exp <- AverageExpression(
      seurat_subset,
      features = genes,
      assays = assay,
      group.by = "ident"
    )
  } else {
    seurat_subset$temp_group_avg <- as.character(seurat_subset@meta.data[[group_by]])
    avg_exp <- AverageExpression(
      seurat_subset,
      features = genes,
      assays = assay,
      group.by = "temp_group_avg"
    )
    seurat_subset$temp_group_avg <- NULL
  }
  exp_matrix <- avg_exp[[assay]]
  # Remove NA column if it slipped through
  exp_matrix <- exp_matrix[, !colnames(exp_matrix) %in% c("NA", NA), drop = FALSE]
  exp_matrix <- exp_matrix[genes[genes %in% rownames(exp_matrix)], , drop = FALSE]
  if (scale_rows) exp_matrix <- t(scale(t(exp_matrix)))
  exp_df <- as.data.frame(exp_matrix)
  exp_df$Gene <- rownames(exp_df)
  exp_long <- tidyr::pivot_longer(exp_df, cols = -Gene, names_to = "Group", values_to = "Expression")
  gene_levels <- rev(unique(genes[genes %in% unique(exp_long$Gene)]))
  group_levels <- unique(group_order[group_order %in% unique(exp_long$Group)])
  exp_long$Gene <- factor(exp_long$Gene, levels = gene_levels)
  exp_long$Group <- factor(exp_long$Group, levels = group_levels)
  p <- ggplot(exp_long, aes(x = Group, y = Gene, fill = Expression)) +
    geom_tile(color = "white", linewidth = 0.5) +
    labs(
      title = paste("Average Expression Heatmap -", assay, "assay"),
      x = if (is.null(group_by)) "Cluster" else tools::toTitleCase(group_by),
      y = "Gene",
      fill = if (scale_rows) "Scaled\nExpression" else "Average\nExpression"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.text.y = element_text(size = 10, face = "bold", hjust = 1),
      panel.grid = element_blank(),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  if (color_palette == "viridis") {
    p <- p + scale_fill_viridis_c(option = "viridis", name = if (scale_rows) "Scaled\nExpression" else "Average\nExpression")
  } else if (color_palette == "magma") {
    p <- p + scale_fill_viridis_c(option = "magma", name = if (scale_rows) "Scaled\nExpression" else "Average\nExpression")
  } else if (color_palette == "inferno") {
    p <- p + scale_fill_viridis_c(option = "inferno", name = if (scale_rows) "Scaled\nExpression" else "Average\nExpression")
  } else if (color_palette == "plasma") {
    p <- p + scale_fill_viridis_c(option = "plasma", name = if (scale_rows) "Scaled\nExpression" else "Average\nExpression")
  } else if (color_palette == "RdYlBu") {
    p <- p + scale_fill_gradient2(
      low = "blue", mid = "yellow", high = "red",
      midpoint = if (scale_rows) 0 else median(exp_long$Expression, na.rm = TRUE),
      name = if (scale_rows) "Scaled\nExpression" else "Average\nExpression"
    )
  } else if (color_palette == "RdBu") {
    p <- p + scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = if (scale_rows) 0 else median(exp_long$Expression, na.rm = TRUE),
      name = if (scale_rows) "Scaled\nExpression" else "Average\nExpression"
    )
  } else if (color_palette == "BlueRed") {
    p <- p + scale_fill_gradient2(
      low = "#0571b0", mid = "#f7f7f7", high = "#ca0020",
      midpoint = if (scale_rows) 0 else median(exp_long$Expression, na.rm = TRUE),
      name = if (scale_rows) "Scaled\nExpression" else "Average\nExpression"
    )
  } else if (color_palette == "YellowRed") {
    p <- p + scale_fill_gradient(
      low = "yellow", high = "red",
      name = if (scale_rows) "Scaled\nExpression" else "Average\nExpression"
    )
  }
  return(p)
}


#' Generate Average Expression Heatmap with Split
#' 
#' @param seurat_object Seurat object
#' @param genes Character vector of gene names
#' @param clusters Vector of cluster names in desired order
#' @param assay Assay to use (default: "RNA")
#' @param split_by Optional metadata column to split by (creates subgroups)
#' @param scale_rows Scale expression values by row (gene) (default: TRUE)
#' @param color_palette Color palette to use (default: viridis)
#' @return ggplot object with heatmap
#' Generate Split Heatmap - USE DASH not underscore
#' Generate split heatmap with average expression
#'
#' Creates heatmap showing average gene expression across clusters with optional splitting
#' Supports dataset filtering and flexible grouping (by dataset or by cluster)
#'
#' @param seurat_object Seurat object
#' @param genes Character vector of gene names
#' @param clusters Character vector of cluster identities to include
#' @param assay Assay name (default: "RNA")
#' @param split_by Metadata column name for splitting (default: NULL for no split)
#' @param scale_rows Scale expression by row/gene (default: TRUE)
#' @param color_palette Color palette name (default: "viridis")
#' @param datasets_to_include Character vector of dataset names to include (default: NULL for all)
#' @param group_by_dataset Logical, if TRUE group by dataset, if FALSE group by cluster (default: TRUE)
#'
#' @return ggplot object
#' @export
generate_split_heatmap <- function(seurat_object,
                                   genes,
                                   clusters,
                                   assay = "RNA",
                                   split_by = NULL,
                                   scale_rows = TRUE,
                                   color_palette = "viridis",
                                   datasets_to_include = NULL,
                                   group_by_dataset = TRUE) {
  DefaultAssay(seurat_object) <- assay
  idents_vec <- as.character(Idents(seurat_object))
  valid_cells <- colnames(seurat_object)[!is.na(idents_vec)]
  if (length(valid_cells) < ncol(seurat_object)) {
    message("[Heatmap] Dropping ", ncol(seurat_object) - length(valid_cells), " cells with NA idents")
    seurat_object <- subset(seurat_object, cells = valid_cells)
  }
  Idents(seurat_object) <- droplevels(Idents(seurat_object))
  if (!is.null(clusters) && length(clusters) > 0) {
    clusters <- clusters[clusters %in% levels(Idents(seurat_object))]
    if (length(clusters) == 0) stop("None of the specified clusters found in object.")
    seurat_subset <- subset(seurat_object, idents = clusters)
    Idents(seurat_subset) <- droplevels(factor(Idents(seurat_subset), levels = clusters))
  } else {
    seurat_subset <- seurat_object
    clusters <- levels(Idents(seurat_subset))
  }
  if (is.null(split_by) || split_by == "None") {
    avg_exp <- AverageExpression(
      seurat_subset,
      features = genes,
      assays = assay,
      group.by = "ident"
    )
    exp_matrix <- avg_exp[[assay]]
    exp_matrix <- exp_matrix[, !colnames(exp_matrix) %in% c("NA", NA), drop = FALSE]
    # AverageExpression replaces underscores with dashes in column names —
    # build a mapping from sanitized names back to original cluster names
    sanitized_clusters <- gsub("_", "-", clusters)
    name_map <- setNames(clusters, sanitized_clusters)
    colnames(exp_matrix) <- ifelse(
      colnames(exp_matrix) %in% names(name_map),
      name_map[colnames(exp_matrix)],
      colnames(exp_matrix)
    )
    exp_matrix <- exp_matrix[genes[genes %in% rownames(exp_matrix)], , drop = FALSE]
    if (scale_rows) exp_matrix <- t(scale(t(exp_matrix)))
    exp_df <- as.data.frame(exp_matrix)
    exp_df$Gene <- rownames(exp_df)
    exp_long <- tidyr::pivot_longer(exp_df, cols = -Gene, names_to = "Group", values_to = "Expression")
    gene_levels <- rev(unique(genes[genes %in% unique(exp_long$Gene)]))
    group_levels <- unique(clusters[clusters %in% unique(exp_long$Group)])
    exp_long$Gene <- factor(exp_long$Gene, levels = gene_levels)
    exp_long$Group <- factor(exp_long$Group, levels = group_levels)
    x_label <- "Cluster"
    subtitle <- NULL
  } else {
    if (!split_by %in% colnames(seurat_subset@meta.data)) {
      stop(paste("Column", split_by, "not found in metadata"))
    }
    valid_cells <- !is.na(seurat_subset@meta.data[[split_by]])
    if (sum(!valid_cells) > 0) {
      seurat_subset <- subset(seurat_subset, cells = colnames(seurat_subset)[valid_cells])
    }
    if (!is.null(datasets_to_include) && length(datasets_to_include) > 0) {
      cells_to_keep <- seurat_subset@meta.data[[split_by]] %in% datasets_to_include
      if (sum(cells_to_keep) == 0) stop("No cells found for selected datasets")
      seurat_subset <- subset(seurat_subset, cells = colnames(seurat_subset)[cells_to_keep])
      message(paste("Filtered to datasets:", paste(datasets_to_include, collapse = ", ")))
    }
    current_idents <- as.character(Idents(seurat_subset))
    split_vals <- as.character(seurat_subset@meta.data[[split_by]])
    temp_groups <- paste0(current_idents, "-", split_vals)
    names(temp_groups) <- colnames(seurat_subset)
    seurat_subset$temp_group_heatmap <- temp_groups
    avg_exp <- AverageExpression(
      seurat_subset,
      features = genes,
      assays = assay,
      group.by = "temp_group_heatmap"
    )
    seurat_subset$temp_group_heatmap <- NULL
    exp_matrix <- avg_exp[[assay]]
    exp_matrix <- exp_matrix[, !colnames(exp_matrix) %in% c("NA", NA), drop = FALSE]
    exp_matrix <- exp_matrix[genes[genes %in% rownames(exp_matrix)], , drop = FALSE]
    if (scale_rows) exp_matrix <- t(scale(t(exp_matrix)))
    exp_df <- as.data.frame(exp_matrix)
    exp_df$Gene <- rownames(exp_df)
    exp_long <- tidyr::pivot_longer(exp_df, cols = -Gene, names_to = "Group", values_to = "Expression")
    gene_levels <- rev(unique(genes[genes %in% unique(exp_long$Gene)]))
    exp_long$Gene <- factor(exp_long$Gene, levels = gene_levels)
    actual_groups <- unique(as.character(exp_long$Group))
    parsed_groups <- data.frame(
      full_name = actual_groups,
      cluster = sub("^([^-]+)-.*", "\\1", actual_groups),
      split_value = sub("^[^-]+-(.*)$", "\\1", actual_groups),
      stringsAsFactors = FALSE
    )
    clusters_in_data <- unique(parsed_groups$cluster)
    split_values_in_data <- unique(parsed_groups$split_value)
    group_order <- character()
    if (group_by_dataset) {
      for (split_val in sort(split_values_in_data)) {
        for (clust in clusters) {
          matching_group <- parsed_groups$full_name[
            parsed_groups$cluster == clust & parsed_groups$split_value == split_val
          ]
          if (length(matching_group) > 0) group_order <- c(group_order, matching_group[1])
        }
      }
      grouping_label <- paste0("grouped by ", split_by)
    } else {
      for (clust in clusters) {
        for (split_val in sort(split_values_in_data)) {
          matching_group <- parsed_groups$full_name[
            parsed_groups$cluster == clust & parsed_groups$split_value == split_val
          ]
          if (length(matching_group) > 0) group_order <- c(group_order, matching_group[1])
        }
      }
      grouping_label <- "grouped by cluster"
    }
    exp_long$Group <- factor(exp_long$Group, levels = group_order)
    x_label <- paste("Cluster x", split_by)
    subtitle <- paste("Split by:", split_by, paste0("(", grouping_label, ")"))
  }
  p <- ggplot(exp_long, aes(x = Group, y = Gene, fill = Expression)) +
    geom_tile(color = "white", linewidth = 0.5) +
    labs(
      title = paste("Heatmap -", assay, "assay"),
      subtitle = subtitle,
      x = x_label,
      y = "Gene",
      fill = if (scale_rows) "Scaled\nExpression" else "Average\nExpression"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 10, face = "bold", hjust = 1),
      panel.grid = element_blank(),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  if (color_palette == "viridis") {
    p <- p + scale_fill_viridis_c(option = "viridis")
  } else if (color_palette == "magma") {
    p <- p + scale_fill_viridis_c(option = "magma")
  } else if (color_palette == "inferno") {
    p <- p + scale_fill_viridis_c(option = "inferno")
  } else if (color_palette == "plasma") {
    p <- p + scale_fill_viridis_c(option = "plasma")
  } else if (color_palette == "RdYlBu") {
    p <- p + scale_fill_gradient2(low = "blue", mid = "yellow", high = "red",
                                  midpoint = if (scale_rows) 0 else median(exp_long$Expression, na.rm = TRUE))
  } else if (color_palette == "RdBu") {
    p <- p + scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                  midpoint = if (scale_rows) 0 else median(exp_long$Expression, na.rm = TRUE))
  } else if (color_palette == "BlueRed") {
    p <- p + scale_fill_gradient2(low = "#0571b0", mid = "#f7f7f7", high = "#ca0020",
                                  midpoint = if (scale_rows) 0 else median(exp_long$Expression, na.rm = TRUE))
  } else if (color_palette == "YellowRed") {
    p <- p + scale_fill_gradient(low = "yellow", high = "red")
  }
  return(p)
}


#' Generate Feature Scatter Plot with Cluster Filtering and Metadata Coloring
#' 
#' @param seurat_object Seurat object
#' @param gene1 First gene name
#' @param gene2 Second gene name
#' @param clusters Optional vector of cluster names to include (NULL = all clusters)
#' @param assay Assay to use (default: "RNA")
#' @param color_by Column to color by: "cluster" or metadata column name (default: "cluster")
#' @param split_by Optional metadata column to split plot into facets (default: NULL)
#' @param threshold_gene1 Minimum expression threshold for gene1 (default: 0)
#' @param threshold_gene2 Minimum expression threshold for gene2 (default: 0)
#' @return ggplot object with scatter plot
generate_feature_scatter <- function(seurat_object,
                                     gene1,
                                     gene2,
                                     clusters = NULL,
                                     assay = "RNA",
                                     color_by = "cluster",
                                     split_by = NULL,
                                     threshold_gene1 = 0,
                                     threshold_gene2 = 0) {
  
  # Set and validate assay
  if (!(assay %in% names(seurat_object@assays))) {
    warning(paste("Assay", assay, "not found. Using RNA instead."))
    assay <- "RNA"
  }
  DefaultAssay(seurat_object) <- assay
  
  # CRITICAL: Join layers for Assay5 before FetchData
  if (inherits(seurat_object[[assay]], "Assay5")) {
    message("Joining layers for Assay5")
    seurat_object <- JoinLayers(seurat_object, assay = assay)
  }
  
  # Validate genes exist
  available_genes <- rownames(seurat_object[[assay]])
  
  if (!(gene1 %in% available_genes)) {
    stop(paste("Gene", gene1, "not found in", assay, "assay"))
  }
  
  if (!(gene2 %in% available_genes)) {
    stop(paste("Gene", gene2, "not found in", assay, "assay"))
  }
  
  # Filter by clusters if specified
  if (!is.null(clusters) && length(clusters) > 0) {
    available_clusters <- levels(Idents(seurat_object))
    valid_clusters <- clusters[clusters %in% available_clusters]
    
    if (length(valid_clusters) == 0) {
      stop("None of the specified clusters found in the dataset.")
    }
    
    seurat_subset <- subset(seurat_object, idents = valid_clusters)
  } else {
    seurat_subset <- seurat_object
  }
  
  # Validate color_by and split_by
  if (color_by != "cluster" && !(color_by %in% colnames(seurat_subset@meta.data))) {
    warning(paste("Column", color_by, "not found in metadata. Using cluster instead."))
    color_by <- "cluster"
  }
  
  if (!is.null(split_by) && split_by != "None") {
    if (!(split_by %in% colnames(seurat_subset@meta.data))) {
      warning(paste("Column", split_by, "not found in metadata. Split disabled."))
      split_by <- NULL
    }
  } else {
    split_by <- NULL
  }
  
  # Extract expression data - use GetAssayData for Assay5 compatibility
  message(paste("Extracting expression data for", gene1, "and", gene2))
  data_layer <- GetAssayData(seurat_subset, assay = assay, layer = "data")
  
  gene1_exp <- data_layer[gene1, ]
  gene2_exp <- data_layer[gene2, ]
  
  # Create data frame
  plot_data <- data.frame(
    Gene1 = as.numeric(gene1_exp),
    Gene2 = as.numeric(gene2_exp),
    row.names = colnames(seurat_subset)
  )
  
  # Apply thresholds - filter cells below threshold
  if (threshold_gene1 > 0 || threshold_gene2 > 0) {
    cells_above_threshold <- (plot_data$Gene1 >= threshold_gene1) & (plot_data$Gene2 >= threshold_gene2)
    plot_data <- plot_data[cells_above_threshold, , drop = FALSE]
    
    # Also subset seurat object to keep metadata aligned
    cells_to_keep <- rownames(plot_data)
    seurat_subset <- subset(seurat_subset, cells = cells_to_keep)
    
    message(paste("After thresholding:", nrow(plot_data), "cells retained"))
  }
  
  # Add color column
  if (color_by == "cluster") {
    plot_data$Color <- as.character(Idents(seurat_subset))
  } else {
    plot_data$Color <- seurat_subset@meta.data[[color_by]]
  }
  
  # Add split column if needed
  if (!is.null(split_by)) {
    plot_data$Split <- seurat_subset@meta.data[[split_by]]
  }
  
  # Create plot
  if (!is.null(split_by)) {
    # With facets
    plot <- ggplot(plot_data, aes(x = Gene1, y = Gene2, color = Color)) +
      geom_point(size = 0.5, alpha = 0.5) +
      facet_wrap(~ Split, ncol = 3) +
      labs(
        x = gene1,
        y = gene2,
        color = if (color_by == "cluster") "Cluster" else color_by
      ) +
      theme_bw() +
      theme(
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 11)
      )
  } else {
    # No facets
    plot <- ggplot(plot_data, aes(x = Gene1, y = Gene2, color = Color)) +
      geom_point(size = 0.5, alpha = 0.5) +
      labs(
        x = gene1,
        y = gene2,
        color = if (color_by == "cluster") "Cluster" else color_by
      ) +
      theme_bw()
  }
  
  # Customize title
  title_parts <- c(paste(gene1, "vs", gene2))
  if (!is.null(clusters) && length(clusters) > 0) {
    title_parts <- c(title_parts, paste("(Clusters:", paste(clusters, collapse = ", "), ")"))
  }
  if (color_by != "cluster") {
    title_parts <- c(title_parts, paste("- Colored by", color_by))
  }
  if (!is.null(split_by)) {
    title_parts <- c(title_parts, paste("- Split by", split_by))
  }
  if (threshold_gene1 > 0 || threshold_gene2 > 0) {
    title_parts <- c(title_parts, paste0("(Thresholds: ", gene1, "≥", threshold_gene1, ", ", gene2, "≥", threshold_gene2, ")"))
  }
  
  plot <- plot + labs(title = paste(title_parts, collapse = " "))
  
  return(plot)
}



# ===============================================================================
# HEATMAP WITH DOTPLOT-STYLE SPLIT/GROUP SYSTEM
# ===============================================================================
# This implements the EXACT same logic as DotPlot for consistency
# Based on working code from multiple_datasets_server_newarch.R lines 1801-1876
# ===============================================================================

# ===============================================================================
# FILE 1: heatmap_functions_newarch.R
# ===============================================================================
# ADD this new function (don't replace generate_split_heatmap):

#' Generate heatmap with DotPlot-style split/group logic
#'
#' Implements the same split/group system as DotPlot for consistency
#'
#' @param seurat_object Seurat object
#' @param genes Character vector of gene names
#' @param assay Assay name (default: "RNA")
#' @param group_by Primary grouping column (like group_by_select in DotPlot)
#' @param metadata_to_compare Optional secondary column for split (like in DotPlot)
#' @param comparison_mode "split" (combine) or "subset" (filter)
#' @param cluster_selection Vector of values to show (from cluster_order input)
#' @param scale_rows Scale by gene (default: TRUE)
#' @param color_palette Color palette (default: "viridis")
#'
#' @return ggplot object
#' @export
generate_heatmap_with_split <- function(seurat_object,
                                        genes,
                                        assay = "RNA",
                                        group_by = "seurat_clusters",
                                        metadata_to_compare = NULL,
                                        comparison_mode = "split",
                                        cluster_selection = NULL,
                                        scale_rows = TRUE,
                                        color_palette = "viridis") {
  
  DefaultAssay(seurat_object) <- assay
  
  # LOGIC 1: group_by="dataset" with metadata_to_compare
  if (group_by == "dataset") {
    if (!is.null(metadata_to_compare) && metadata_to_compare != "") {
      if (comparison_mode == "split") {
        # MODE SPLIT: Create combined categories (dataset_cluster)
        # Example: WT_IIa, WT_IIb, mdx_IIa, mdx_IIb
        seurat_object$heatmap_category <- paste(seurat_object$dataset,
                                                seurat_object@meta.data[[metadata_to_compare]],
                                                sep = "_")
        
        # Filter to valid combinations
        if (!is.null(cluster_selection) && length(cluster_selection) > 0) {
          valid_combinations <- unlist(lapply(unique(seurat_object$dataset), function(ds) {
            paste(ds, cluster_selection, sep = "_")
          }))
          
          cells_to_keep <- seurat_object$heatmap_category %in% valid_combinations
          seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
          seurat_object$heatmap_category <- factor(seurat_object$heatmap_category,
                                                   levels = valid_combinations)
        }
        
        group_column <- "heatmap_category"
        
      } else {
        # MODE SUBSET: Compare datasets for selected clusters only
        # Example: Compare WT vs mdx for IIa only
        cells_subset <- seurat_object@meta.data[[metadata_to_compare]] %in% cluster_selection
        seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_subset])
        group_column <- "dataset"
      }
    } else {
      # Simple grouping by dataset
      group_column <- "dataset"
    }
    
  } else {
    # LOGIC 2: group_by != "dataset" (standard grouping by any metadata column)
    if (!is.null(cluster_selection) && length(cluster_selection) > 0) {
      cells_to_keep <- seurat_object@meta.data[[group_by]] %in% cluster_selection
      seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
      seurat_object@meta.data[[group_by]] <- factor(
        seurat_object@meta.data[[group_by]],
        levels = cluster_selection
      )
    }
    group_column <- group_by
  }
  
  # Compute average expression
  avg_exp <- AverageExpression(
    seurat_object,
    features = genes,
    assays = assay,
    group.by = group_column
  )
  
  exp_matrix <- avg_exp[[assay]]
  exp_matrix <- exp_matrix[genes[genes %in% rownames(exp_matrix)], , drop = FALSE]
  
  # Scale by row if requested
  if (scale_rows) {
    exp_matrix <- t(scale(t(exp_matrix)))
  }
  
  # Convert to long format for ggplot
  exp_df <- as.data.frame(exp_matrix)
  exp_df$Gene <- rownames(exp_df)
  exp_long <- tidyr::pivot_longer(exp_df, cols = -Gene, names_to = "Group", values_to = "Expression")
  
  # Factor genes and groups - CRITICAL: ensure unique levels
  gene_levels <- rev(unique(genes[genes %in% unique(exp_long$Gene)]))
  group_levels <- unique(colnames(exp_matrix))
  exp_long$Gene <- factor(exp_long$Gene, levels = gene_levels)
  exp_long$Group <- factor(exp_long$Group, levels = group_levels)
  
  # Create plot
  p <- ggplot(exp_long, aes(x = Group, y = Gene, fill = Expression)) +
    geom_tile(color = "white", linewidth = 0.5) +
    labs(
      title = paste("Heatmap -", assay, "assay"),
      x = tools::toTitleCase(gsub("_", " ", group_by)),
      y = "Gene",
      fill = if (scale_rows) "Scaled\nExpression" else "Average\nExpression"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      axis.text.y = element_text(size = 10, face = "bold", hjust = 1),
      panel.grid = element_blank(),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  
  # Apply color palette
  if (color_palette == "viridis") {
    p <- p + scale_fill_viridis_c(option = "viridis")
  } else if (color_palette == "magma") {
    p <- p + scale_fill_viridis_c(option = "magma")
  } else if (color_palette == "inferno") {
    p <- p + scale_fill_viridis_c(option = "inferno")
  } else if (color_palette == "plasma") {
    p <- p + scale_fill_viridis_c(option = "plasma")
  } else if (color_palette == "RdYlBu") {
    p <- p + scale_fill_gradient2(low = "blue", mid = "yellow", high = "red", 
                                  midpoint = if (scale_rows) 0 else median(exp_long$Expression, na.rm = TRUE))
  } else if (color_palette == "RdBu") {
    p <- p + scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                                  midpoint = if (scale_rows) 0 else median(exp_long$Expression, na.rm = TRUE))
  } else if (color_palette == "BlueRed") {
    p <- p + scale_fill_gradient2(low = "#0571b0", mid = "#f7f7f7", high = "#ca0020", 
                                  midpoint = if (scale_rows) 0 else median(exp_long$Expression, na.rm = TRUE))
  } else if (color_palette == "YellowRed") {
    p <- p + scale_fill_gradient(low = "yellow", high = "red")
  }
  
  # Clean up temporary columns
  if ("heatmap_category" %in% colnames(seurat_object@meta.data)) {
    seurat_object$heatmap_category <- NULL
  }
  
  return(p)
}

