# ============================================================================
# Heatmap Functions for Cell-Hub
# ============================================================================
# This file contains helper functions for generating and customizing heatmaps
# Author: Gaspard
# Last updated: 2025-11-10
# ============================================================================




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
    message("[get_top_markers] Dropping ", ncol(seurat_object) - length(valid_cells), " cells with NA idents")
    seurat_object <- subset(seurat_object, cells = valid_cells)
  }
  Idents(seurat_object) <- droplevels(Idents(seurat_object))
  # DE reads RNA (the integrated assay holds ~2,000 corrected genes and no counts); layers joined and, for SCTransform-only
  # objects whose RNA assay has no data layer, a log-normalized copy is made (FindAllMarkers otherwise returns 0 rows silently)
  de_assay <- if (assay == "integrated") "RNA" else assay
  seurat_object <- prepare_rna_assay_for_de(seurat_object, de_assay)$obj
  # FindAllMarkers requires a sparse counts layer
  counts_layer <- LayerData(seurat_object, assay = de_assay, layer = "counts")
  if (!inherits(counts_layer, "dgCMatrix")) LayerData(seurat_object, assay = de_assay, layer = "counts") <- as(counts_layer, "CsparseMatrix")
  markers <- FindAllMarkers(seurat_object, min.pct = min_pct, assay = de_assay, verbose = FALSE)
  if (is.null(markers) || nrow(markers) == 0) stop("No DE genes identified. Try lowering min_pct or logfc_threshold.")
  if (!is.null(clusters_to_include) && length(clusters_to_include) > 0) {
    markers <- markers %>% dplyr::filter(cluster %in% clusters_to_include)
    message(paste("Filtered markers to", length(unique(markers$cluster)), "selected clusters:", paste(clusters_to_include, collapse = ", ")))
  }
  top_markers_df <- markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > logfc_threshold) %>%
    slice_head(n = n_genes) %>%
    ungroup()
  if (nrow(top_markers_df) == 0) stop(sprintf("No marker with log2FC > %g in the selected clusters", logfc_threshold))
  return(list(genes = unique(top_markers_df$gene), markers_table = as.data.frame(top_markers_df)))
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
  if (!(assay %in% names(seurat_object@assays))) {
    warning(paste("Assay", assay, "not found. Using RNA instead."))
    assay <- "RNA"
  }
  DefaultAssay(seurat_object) <- assay
  available_genes <- rownames(seurat_object[[assay]])
  if (!(gene1 %in% available_genes)) stop(paste("Gene", gene1, "not found in", assay, "assay"))
  if (!(gene2 %in% available_genes)) stop(paste("Gene", gene2, "not found in", assay, "assay"))
  # Cells of the selected clusters (identities), chosen before any data is read
  idents <- as.character(Idents(seurat_object))
  cells <- colnames(seurat_object)
  if (!is.null(clusters) && length(clusters) > 0) {
    valid_clusters <- clusters[clusters %in% unique(idents)]
    if (length(valid_clusters) == 0) stop("None of the specified clusters found in the dataset.")
    cells <- cells[idents %in% valid_clusters]
  }
  if (color_by != "cluster" && !(color_by %in% colnames(seurat_object@meta.data))) {
    warning(paste("Column", color_by, "not found in metadata. Using cluster instead."))
    color_by <- "cluster"
  }
  if (is.null(split_by) || split_by == "None" || !(split_by %in% colnames(seurat_object@meta.data))) split_by <- NULL
  # FetchData reads every layer of a split Seurat v5 assay: no JoinLayers copy of the whole assay is needed for two genes
  expr <- FetchData(seurat_object, vars = c(gene1, gene2), cells = cells, layer = "data")
  plot_data <- data.frame(Gene1 = as.numeric(expr[[1]]), Gene2 = as.numeric(expr[[2]]), row.names = rownames(expr))
  if (threshold_gene1 > 0 || threshold_gene2 > 0) {
    keep <- plot_data$Gene1 >= threshold_gene1 & plot_data$Gene2 >= threshold_gene2
    if (!any(keep)) stop(sprintf("No cell has %s >= %g and %s >= %g", gene1, threshold_gene1, gene2, threshold_gene2))
    plot_data <- plot_data[keep, , drop = FALSE]
    message(paste("After thresholding:", nrow(plot_data), "cells retained"))
  }
  # Colour and panel values matched by barcode
  plot_data$Color <- if (color_by == "cluster") idents[match(rownames(plot_data), colnames(seurat_object))] else seurat_object@meta.data[rownames(plot_data), color_by]
  if (!is.null(split_by)) plot_data$Split <- seurat_object@meta.data[rownames(plot_data), split_by]
  plot <- ggplot(plot_data, aes(x = Gene1, y = Gene2, color = Color)) +
    geom_point(size = 0.5, alpha = 0.5) +
    labs(x = gene1, y = gene2, color = if (color_by == "cluster") "Cluster" else color_by) +
    theme_bw()
  if (!is.null(split_by)) {
    plot <- plot + facet_wrap(~ Split, ncol = 3) +
      theme(strip.background = element_rect(fill = "white"), strip.text = element_text(face = "bold", size = 11))
  }
  title_parts <- c(paste(gene1, "vs", gene2))
  if (!is.null(clusters) && length(clusters) > 0) title_parts <- c(title_parts, paste("(Clusters:", paste(clusters, collapse = ", "), ")"))
  if (color_by != "cluster") title_parts <- c(title_parts, paste("- Colored by", color_by))
  if (!is.null(split_by)) title_parts <- c(title_parts, paste("- Split by", split_by))
  if (threshold_gene1 > 0 || threshold_gene2 > 0) title_parts <- c(title_parts, paste0("(Thresholds: ", gene1, "≥", threshold_gene1, ", ", gene2, "≥", threshold_gene2, ")"))
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


#' Top differentially expressed genes between two groups of a metadata column, within each selected cluster
#'
#' Ranking for display: cell-level Wilcoxon test (FindMarkers, presto backend when installed). Cells of one animal are
#' not independent replicates, so the adjusted p-values only rank genes; sample-level inference uses the pseudobulk DE.
#' @param seurat_object Seurat object with active identities set to the clusters
#' @param split_by Metadata column holding the conditions (e.g. genotype)
#' @param group1,group2 Values of split_by pooled into group 1 (numerator of log2FC) and group 2
#' @param clusters Clusters tested one by one (NULL = all)
#' @param n_genes Number of genes kept per cluster
#' @param direction "both" (largest |log2FC|), "up" (higher in group 1) or "down" (higher in group 2)
#' @return List with genes (unique, ordered by cluster then log2FC), markers_table (data.frame) and skipped (clusters not tested)
get_top_condition_de_genes <- function(seurat_object, split_by, group1, group2, clusters = NULL, n_genes = 10, direction = "both",
                                       assay = "RNA", min_pct = 0.1, logfc_threshold = 0.25, max_padj = 0.05, min_cells = 3) {
  if (!split_by %in% colnames(seurat_object@meta.data)) stop("Column ", split_by, " not found in metadata")
  if (length(group1) == 0 || length(group2) == 0) stop("Select at least one value for each group")
  shared <- intersect(group1, group2)
  if (length(shared) > 0) stop("Value(s) in both groups: ", paste(shared, collapse = ", "))
  if (!direction %in% c("both", "up", "down")) stop("direction must be 'both', 'up' or 'down'")
  # RNA log-normalized data with joined layers (FindMarkers reads one matrix)
  obj <- prepare_rna_assay_for_de(seurat_object, assay)$obj
  idents <- as.character(Idents(obj))
  if (is.null(clusters) || length(clusters) == 0) clusters <- levels(droplevels(Idents(obj)))
  split_vals <- trimws(as.character(obj@meta.data[[split_by]]))
  condition <- ifelse(split_vals %in% group1, "Group1", ifelse(split_vals %in% group2, "Group2", NA_character_))
  tables <- list(); skipped <- character()
  for (cl in clusters) {
    in_cl <- idents == cl & !is.na(condition)
    n1 <- sum(in_cl & condition == "Group1"); n2 <- sum(in_cl & condition == "Group2")
    if (n1 < min_cells || n2 < min_cells) { skipped <- c(skipped, sprintf("%s: %d vs %d cells", cl, n1, n2)); next }
    sub <- subset(obj, cells = colnames(obj)[in_cl])
    # Labels matched by barcode, so the result does not depend on the cell order returned by subset()
    sub$heatmap_condition <- condition[match(colnames(sub), colnames(obj))]
    Idents(sub) <- factor(sub$heatmap_condition, levels = c("Group1", "Group2"))
    de <- FindMarkers(sub, ident.1 = "Group1", ident.2 = "Group2", assay = assay, min.pct = min_pct, logfc.threshold = logfc_threshold, verbose = FALSE)
    de <- de[!is.na(de$p_val_adj) & de$p_val_adj < max_padj, , drop = FALSE]
    if (direction == "up") de <- de[de$avg_log2FC > 0, , drop = FALSE]
    if (direction == "down") de <- de[de$avg_log2FC < 0, , drop = FALSE]
    if (nrow(de) == 0) { skipped <- c(skipped, sprintf("%s: no gene with adj. p < %g", cl, max_padj)); next }
    # Strongest effects first, then displayed from most up to most down within the cluster
    de <- head(de[order(-abs(de$avg_log2FC), de$p_val_adj), , drop = FALSE], n_genes)
    de <- de[order(-de$avg_log2FC), , drop = FALSE]
    tables[[cl]] <- data.frame(cluster = cl, gene = rownames(de), avg_log2FC = de$avg_log2FC, pct.1 = de$pct.1, pct.2 = de$pct.2,
                               p_val = de$p_val, p_val_adj = de$p_val_adj, n_cells_group1 = n1, n_cells_group2 = n2, stringsAsFactors = FALSE)
  }
  if (length(tables) == 0) stop("No differentially expressed gene found (", paste(skipped, collapse = "; "), ")")
  markers_table <- do.call(rbind, tables)
  rownames(markers_table) <- NULL
  markers_table$group1 <- paste(group1, collapse = " + ")
  markers_table$group2 <- paste(group2, collapse = " + ")
  return(list(genes = unique(markers_table$gene), markers_table = markers_table, skipped = skipped))
}

#' Average expression per group computed directly from the log-normalized data layer
#'
#' Same convention as Seurat AverageExpression on the data layer (mean in linear space) returned on the log scale,
#' log1p(mean(expm1(x))), without Seurat renaming the groups (underscores, leading digits).
#' @param seurat_object Seurat object whose assay has a single (joined) data layer
#' @param groups Factor of length ncol(seurat_object); NA cells are ignored; its levels give the column order
#' @return Genes x groups numeric matrix
compute_group_average_expression <- function(seurat_object, genes, groups, assay = "RNA") {
  if (length(groups) != ncol(seurat_object)) stop("groups must have one value per cell")
  # Split data layers (merged, unjoined object): LayerData() would return the first layer only and silently bias every average
  if (inherits(seurat_object[[assay]], "Assay5") && length(SeuratObject::Layers(seurat_object[[assay]], search = "data")) > 1) {
    seurat_object <- SeuratObject::JoinLayers(seurat_object, assay = assay)
  }
  data <- SeuratObject::LayerData(seurat_object, assay = assay, layer = "data")
  if (ncol(data) != length(groups)) stop("The data layer of the ", assay, " assay does not cover every cell")
  data <- data[, colnames(seurat_object), drop = FALSE]
  genes <- unique(genes)
  missing_genes <- setdiff(genes, rownames(data))
  if (length(missing_genes) > 0) stop("Genes not found in the ", assay, " assay: ", paste(missing_genes, collapse = ", "))
  data <- data[genes, , drop = FALSE]
  if (inherits(data, "dgCMatrix")) data@x <- expm1(data@x) else data <- expm1(as.matrix(data))
  groups <- droplevels(as.factor(groups))
  keep <- which(!is.na(groups))
  indicator <- Matrix::sparseMatrix(i = keep, j = as.integer(groups[keep]), x = 1, dims = c(length(groups), nlevels(groups)), dimnames = list(NULL, levels(groups)))
  sums <- as.matrix(data %*% indicator)
  avg <- sweep(sums, 2, Matrix::colSums(indicator), "/")
  rownames(avg) <- genes
  return(log1p(avg))
}

#' Heatmap of a genes x groups matrix, with an optional second axis drawn as facets
#'
#' @param avg_matrix Genes x groups matrix (log-scale average expression); row order is kept from top to bottom
#' @param col_info data.frame with id (column names of avg_matrix, in display order), x (axis label) and facet (panel, NA if none)
#' @param scale_mode "zscore" (per gene across columns) or "none" (average log-normalized expression)
#' @return ggplot object
plot_grouped_expression_heatmap <- function(avg_matrix, col_info, scale_mode = "zscore", color_palette = "viridis",
                                            title = NULL, subtitle = NULL, x_title = "Cluster") {
  mat <- avg_matrix[, col_info$id, drop = FALSE]
  if (scale_mode == "zscore") {
    mat <- t(scale(t(mat)))
    # Genes with identical values in every column have no z-score
    mat[is.na(mat)] <- 0
  }
  long <- data.frame(Gene = rep(rownames(mat), times = ncol(mat)), id = rep(colnames(mat), each = nrow(mat)), Expression = as.vector(mat), stringsAsFactors = FALSE)
  long$x <- col_info$x[match(long$id, col_info$id)]
  long$facet <- col_info$facet[match(long$id, col_info$id)]
  long$Gene <- factor(long$Gene, levels = rev(unique(rownames(mat))))
  long$x <- factor(long$x, levels = unique(col_info$x))
  has_facet <- !all(is.na(col_info$facet))
  if (has_facet) long$facet <- factor(long$facet, levels = unique(col_info$facet))
  midpoint <- if (scale_mode == "zscore") 0 else stats::median(long$Expression, na.rm = TRUE)
  fill_scale <- switch(color_palette,
                       magma = scale_fill_viridis_c(option = "magma"),
                       inferno = scale_fill_viridis_c(option = "inferno"),
                       plasma = scale_fill_viridis_c(option = "plasma"),
                       RdYlBu = scale_fill_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = midpoint),
                       RdBu = scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = midpoint),
                       BlueRed = scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = midpoint),
                       YellowRed = scale_fill_gradient(low = "#FFFFCC", high = "#BD0026"),
                       scale_fill_viridis_c(option = "viridis"))
  p <- ggplot(long, aes(x = x, y = Gene, fill = Expression)) +
    geom_tile(color = "white", linewidth = 0.5) + fill_scale +
    labs(title = title, subtitle = subtitle, x = x_title, y = "Gene",
         fill = if (scale_mode == "zscore") "Row z-score" else "Average\nlog expression") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10), axis.text.y = element_text(size = 10, face = "bold", hjust = 1),
          panel.grid = element_blank(), axis.title = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          strip.text = element_text(face = "bold", size = 11), panel.spacing.x = unit(4, "pt"))
  # Second comparison axis: one panel per level, panel width proportional to its number of columns
  if (has_facet) p <- p + facet_grid(. ~ facet, scales = "free_x", space = "free_x")
  return(p)
}