# ==============================================================================
# trajectory_functions.R
# Pure helper functions of the trajectory module (Monocle 3), no Shiny reactive code.
# Conversion of a Seurat object to a cell_data_set, principal graph learning (seed 42), root selection and pseudotime,
# genes varying along the trajectory (graph_test) or correlated with pseudotime, and the plots and tables of the module.
# Reactive logic lives in trajectory_server.R.
# ==============================================================================

############################## Object checks ##############################

# Cluster column of a cell_data_set: annotated names (ClusterIdents) first, then Seurat cluster numbers
get_cds_cluster_column <- function(cds) {
  cols <- colnames(colData(cds))
  if ("ClusterIdents" %in% cols) return("ClusterIdents")
  if ("seurat_clusters" %in% cols) return("seurat_clusters")
  stop("No cluster column (ClusterIdents or seurat_clusters) found in the object")
}

# Sorted cluster names of a cell_data_set, or NULL when it has no cluster column
get_cds_clusters <- function(cds) {
  col <- tryCatch(get_cds_cluster_column(cds), error = function(e) NULL)
  if (is.null(col)) return(NULL)
  return(sort(unique(as.character(colData(cds)[[col]]))))
}

# TRUE when the principal graph has been learned
has_trajectory_graph <- function(cds) {
  return(!is.null(cds@principal_graph) && length(cds@principal_graph) > 0)
}

# TRUE when pseudotime has been computed for at least one cell
has_valid_pseudotime <- function(cds) {
  pt <- tryCatch(pseudotime(cds), error = function(e) NULL)
  return(!is.null(pt) && !all(is.na(pt)))
}

# Checks required before Monocle conversion; returns the object with split counts layers joined (Seurat v5 rules 1-2)
prepare_seurat_for_monocle <- function(seurat_obj) {
  if (!inherits(seurat_obj, "Seurat")) stop("The object is not a Seurat object")
  if (!"RNA" %in% names(seurat_obj@assays)) stop("No RNA assay found: Monocle conversion requires RNA counts")
  # The conversion transfers the reduction named exactly 'umap'
  if (!"umap" %in% names(seurat_obj@reductions)) stop("No 'umap' reduction found: run UMAP before trajectory analysis")
  if (inherits(seurat_obj[["RNA"]], "Assay5") && length(SeuratObject::Layers(seurat_obj[["RNA"]], search = "counts")) > 1) {
    seurat_obj <- SeuratObject::JoinLayers(seurat_obj, assay = "RNA")
  }
  return(seurat_obj)
}

############################## Conversion and graph ##############################

# One partition for all cells and the Seurat clusters as Monocle clusters (no re-clustering), so the trajectory
# connects the existing clusters
set_monocle_clusters <- function(cds, cluster_column) {
  cell_names <- colnames(cds)
  partitions <- factor(rep("1", ncol(cds)))
  names(partitions) <- cell_names
  clusters <- factor(colData(cds)[[cluster_column]])
  names(clusters) <- cell_names
  cds@clusters <- S4Vectors::SimpleList()
  cds@clusters[["UMAP"]] <- S4Vectors::SimpleList(partitions = partitions, clusters = clusters)
  return(cds)
}

# Seurat object -> Monocle 3 cell_data_set: RNA counts, metadata (ClusterIdents added from the identities when absent),
# PCA and UMAP embeddings transferred, size factors estimated, Seurat clusters kept
seurat_to_monocle <- function(seurat_obj) {
  seurat_obj <- prepare_seurat_for_monocle(seurat_obj)
  counts <- SeuratObject::LayerData(seurat_obj, assay = "RNA", layer = "counts")
  cell_metadata <- seurat_obj@meta.data
  if (!"ClusterIdents" %in% colnames(cell_metadata)) {
    if (length(Idents(seurat_obj)) == ncol(seurat_obj)) {
      cell_metadata$ClusterIdents <- as.character(Idents(seurat_obj))
    } else if ("seurat_clusters" %in% colnames(cell_metadata)) {
      cell_metadata$ClusterIdents <- as.character(cell_metadata$seurat_clusters)
    } else {
      stop("No cluster information found in the Seurat object")
    }
  }
  gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  cds <- new_cell_data_set(expression_data = counts, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
  if ("pca" %in% names(seurat_obj@reductions)) {
    pca <- Embeddings(seurat_obj, "pca")
    if (nrow(pca) == ncol(cds)) reducedDims(cds)[["PCA"]] <- pca
  }
  umap <- Embeddings(seurat_obj, "umap")
  if (nrow(umap) != ncol(cds)) stop("UMAP dimensions do not match the number of cells")
  reducedDims(cds)[["UMAP"]] <- umap
  cds <- estimate_size_factors(cds)
  cds <- set_monocle_clusters(cds, "ClusterIdents")
  return(cds)
}

# Principal graph on the UMAP with one partition and no closed loops (seed 42: the graph initialization is random)
learn_trajectory_graph <- function(cds, seed = 42) {
  cds <- set_monocle_clusters(cds, get_cds_cluster_column(cds))
  set.seed(seed)
  cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)
  return(cds)
}

############################## Root and pseudotime ##############################

# Root cell of a cluster: the cell closest to the cluster centre on the UMAP
find_root_cell <- function(cds, root_cluster) {
  cluster_column <- get_cds_cluster_column(cds)
  in_cluster <- which(as.character(colData(cds)[[cluster_column]]) == as.character(root_cluster))
  if (length(in_cluster) == 0) stop("No cells found in the selected cluster")
  coords <- reducedDims(cds)[["UMAP"]][in_cluster, , drop = FALSE]
  distances <- sqrt(rowSums(sweep(coords, 2, colMeans(coords), "-")^2))
  return(list(cell = colnames(cds)[in_cluster][which.min(distances)], n_cells = length(in_cluster)))
}

# Pseudotime from a root cell, also stored in colData$pseudotime
order_cells_from_root <- function(cds, root_cell) {
  cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cell)
  pt <- pseudotime(cds)
  if (is.null(pt) || all(is.na(pt))) stop("Pseudotime could not be calculated from this root")
  colData(cds)$pseudotime <- pt
  return(cds)
}

# History entry of one pseudotime calculation (the picker lets the user switch between roots)
new_pseudotime_entry <- function(cds, root_cluster, root_cell, label = paste("Root: Cluster", root_cluster)) {
  id <- paste0("Root_Cluster_", gsub("[^A-Za-z0-9]+", "_", root_cluster), "_", format(Sys.time(), "%H%M%S"))
  return(list(id = id, label = label, root_cluster = root_cluster, root_cell = root_cell, monocle_object = cds, timestamp = Sys.time()))
}

# Picker choices of the pseudotime history: entry ids labelled by their root
pseudotime_choices <- function(history) {
  if (length(history) == 0) return(character(0))
  return(stats::setNames(names(history), vapply(history, function(x) x$label, character(1))))
}

# File name part describing the selected pseudotime root
pseudotime_root_tag <- function(entry) {
  if (is.null(entry)) return(NULL)
  return(gsub("[^A-Za-z0-9]+", "_", as.character(entry$root_cluster)))
}

############################## Trajectory plots ##############################

# Trajectory graph on the UMAP, cells coloured by cluster
plot_trajectory_clusters <- function(cds, title = "Trajectory on Clusters") {
  return(plot_cells(cds, color_cells_by = get_cds_cluster_column(cds), cell_size = 0.5, label_groups_by_cluster = TRUE,
                    label_leaves = FALSE, label_branch_points = TRUE, graph_label_size = 1.5) +
           ggplot2::ggtitle(title) + ggplot2::theme_minimal())
}

# Trajectory graph on the UMAP, cells coloured by pseudotime
plot_trajectory_pseudotime <- function(cds, title = "Pseudotime Trajectory") {
  return(plot_cells(cds, color_cells_by = "pseudotime", cell_size = 0.5, label_cell_groups = FALSE, label_leaves = FALSE,
                    label_branch_points = TRUE, graph_label_size = 1.5) +
           ggplot2::ggtitle(title) + ggplot2::theme_minimal() + viridis::scale_color_viridis(name = "Pseudotime", option = "plasma"))
}

############################## Genes along the trajectory ##############################

# Genes with structure along the principal graph (Moran's I, graph_test), independent of the root.
# Returns list(significant, n_tested, min_q); significant is sorted by Moran's I and may be empty.
run_trajectory_graph_test <- function(cds, q_threshold = 0.05, max_cores = 4) {
  if (!has_trajectory_graph(cds)) stop("Trajectory not calculated: construct the graph first")
  cores <- max(1, min(max_cores, parallel::detectCores() - 1, na.rm = TRUE))
  markers <- graph_test(cds, neighbor_graph = "principal_graph", cores = cores)
  if (is.null(markers) || nrow(markers) == 0) stop("graph_test returned no result")
  markers <- as.data.frame(markers)
  significant <- markers[!is.na(markers$q_value) & markers$q_value < q_threshold, , drop = FALSE]
  significant <- significant[order(-significant$morans_I), , drop = FALSE]
  return(list(significant = significant, n_tested = nrow(markers), min_q = suppressWarnings(min(markers$q_value, na.rm = TRUE))))
}

## Group: Spearman correlation of each gene with pseudotime (genes detected in at least min_cells cells), Benjamini-Hochberg.
## Group: Expression = counts divided by the cell size factors (no log: Spearman uses ranks). Cells with a non-finite pseudotime
## Group: (NA, or Inf when unreachable from the root) are excluded. Returns list(significant, n_tested);
## Group: significant is sorted by |correlation| and may be empty.
run_pseudotime_correlation <- function(cds, q_threshold = 0.05, min_cells = 10) {
  if (!has_valid_pseudotime(cds)) stop("Pseudotime not calculated: select a root first")
  if (is.null(SummarizedExperiment::colData(cds)$Size_Factor)) cds <- monocle3::estimate_size_factors(cds)
  pt <- pseudotime(cds)
  valid <- is.finite(pt)
  if (sum(valid) < min_cells) stop("Too few cells with a finite pseudotime")
  expr <- monocle3::normalized_counts(cds, norm_method = "size_only")[, valid, drop = FALSE]
  pt <- pt[valid]
  keep <- Matrix::rowSums(expr > 0) >= min_cells
  expr <- expr[keep, , drop = FALSE]
  if (nrow(expr) == 0) stop(sprintf("No gene detected in at least %d cells with pseudotime", min_cells))
  results <- do.call(rbind, lapply(seq_len(nrow(expr)), function(i) {
    test <- suppressWarnings(stats::cor.test(as.numeric(expr[i, ]), pt, method = "spearman", exact = FALSE))
    data.frame(gene = rownames(expr)[i], correlation = as.numeric(test$estimate), p_value = test$p.value, stringsAsFactors = FALSE)
  }))
  results <- results[!is.na(results$correlation), , drop = FALSE]
  results$q_value <- stats::p.adjust(results$p_value, method = "BH")
  significant <- results[results$q_value < q_threshold, , drop = FALSE]
  significant <- significant[order(-abs(significant$correlation)), , drop = FALSE]
  rownames(significant) <- significant$gene
  return(list(significant = significant, n_tested = nrow(results)))
}

# Display table of graph_test results (p-values kept numeric, formatted by DT)
graph_test_datatable <- function(results) {
  df <- as.data.frame(results)
  if (!"gene_short_name" %in% colnames(df)) df$gene_short_name <- rownames(df)
  cols <- intersect(c("gene_short_name", "morans_I", "q_value", "p_value", "status"), colnames(df))
  df <- df[, cols, drop = FALSE]
  colnames(df) <- c(gene_short_name = "Gene", morans_I = "Moran's I", q_value = "Q-value", p_value = "P-value", status = "Status")[cols]
  return(DT::datatable(df, rownames = FALSE, filter = "top", caption = "Genes with structure along the trajectory (independent of the root)",
                       options = list(pageLength = 25, scrollX = TRUE, order = list(list(1, "desc")))) |>
           DT::formatRound("Moran's I", digits = 4) |>
           DT::formatSignif(intersect(c("Q-value", "P-value"), colnames(df)), digits = 3) |>
           DT::formatStyle("Moran's I", color = DT::styleInterval(0.3, c("black", "darkred")), fontWeight = DT::styleInterval(0.3, c("normal", "bold"))) |>
           DT::formatStyle("Q-value", backgroundColor = DT::styleInterval(c(0.01, 0.05), c("lightgreen", "lightyellow", "white"))))
}

# Display table of the pseudotime correlation (red = increases, blue = decreases along pseudotime)
correlation_datatable <- function(results, caption = NULL) {
  df <- results[, c("gene", "correlation", "p_value", "q_value"), drop = FALSE]
  colnames(df) <- c("Gene", "Correlation", "P-value", "Q-value")
  df$Direction <- ifelse(df$Correlation > 0, "Up", "Down")
  df <- df[, c("Gene", "Correlation", "Direction", "P-value", "Q-value"), drop = FALSE]
  return(DT::datatable(df, rownames = FALSE, filter = "top", caption = caption,
                       options = list(pageLength = 25, scrollX = TRUE, order = list(list(1, "desc")))) |>
           DT::formatRound("Correlation", digits = 3) |>
           DT::formatSignif(c("P-value", "Q-value"), digits = 3) |>
           DT::formatStyle("Correlation", color = DT::styleInterval(c(-0.3, 0.3), c("blue", "gray", "red")),
                           fontWeight = DT::styleInterval(c(-0.5, 0.5), c("bold", "normal", "bold"))) |>
           DT::formatStyle("Direction", backgroundColor = DT::styleEqual(c("Up", "Down"), c("#ffcccc", "#ccccff"))) |>
           DT::formatStyle("Q-value", backgroundColor = DT::styleInterval(c(0.01, 0.05), c("lightgreen", "lightyellow", "white"))))
}

############################## Gene plots ##############################

# Expression of one gene on the trajectory UMAP
plot_gene_on_trajectory <- function(cds, gene, show_graph = TRUE, cell_size = 0.5) {
  if (!gene %in% rownames(cds)) stop("Selected gene not found in the dataset")
  return(plot_cells(cds, genes = gene, show_trajectory_graph = show_graph, label_cell_groups = FALSE, label_leaves = FALSE,
                    label_branch_points = show_graph, cell_size = cell_size, cell_stroke = 0.1,
                    trajectory_graph_color = "black", trajectory_graph_segment_size = 0.75) +
           ggplot2::ggtitle(paste("Expression of", gene, "along trajectory")) +
           ggplot2::theme_minimal(base_size = 14) +
           ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16), legend.position = "right") +
           viridis::scale_color_viridis(name = "Expression", option = "magma"))
}

## Group: Cells of the selected clusters with a finite pseudotime, pseudotime carried over (colData$Pseudotime and the principal graph slot)
subset_cds_by_clusters <- function(cds, clusters) {
  pt <- pseudotime(cds)
  if (is.null(pt) || all(is.na(pt))) stop("Pseudotime not calculated: select a root first")
  # Cells unreachable from the root (Inf) would break the spline trend of plot_genes_in_pseudotime
  keep <- is.finite(pt)
  if (length(clusters) > 0) keep <- keep & as.character(colData(cds)[[get_cds_cluster_column(cds)]]) %in% clusters
  if (!any(keep)) stop("No cells with a finite pseudotime in the selected clusters")
  if (all(keep)) return(cds)
  sub <- cds[, keep]
  colData(sub)$Pseudotime <- pt[keep]
  if ("UMAP" %in% names(cds@principal_graph_aux) && "pseudotime" %in% names(cds@principal_graph_aux[["UMAP"]])) {
    sub@principal_graph_aux[["UMAP"]]$pseudotime <- pt[keep]
  }
  return(sub)
}

# Expression of several genes along pseudotime, cells coloured by cluster, with a spline trend (df = 3) per gene
plot_genes_along_pseudotime <- function(cds, genes, min_expr = 0.1, cell_size = 0.75, trend_size = 1) {
  if (!has_valid_pseudotime(cds)) stop("Pseudotime not calculated: select a root first")
  if (!"Pseudotime" %in% colnames(colData(cds))) colData(cds)$Pseudotime <- pseudotime(cds)
  genes <- genes[genes %in% rownames(cds)]
  if (length(genes) == 0) stop("Selected genes not found in the dataset")
  p <- plot_genes_in_pseudotime(cds[genes, ], color_cells_by = get_cds_cluster_column(cds), min_expr = min_expr, cell_size = cell_size,
                                nrow = NULL, ncol = min(3, length(genes)), panel_order = genes,
                                trend_formula = "~ splines::ns(Pseudotime, df=3)", label_by_short_name = TRUE)
  # Width of the trend lines (geom_line layers added by Monocle)
  for (i in seq_along(p$layers)) if (inherits(p$layers[[i]]$geom, "GeomLine")) p$layers[[i]]$aes_params$linewidth <- trend_size
  p <- p + ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_blank(), legend.position = "right",
                   legend.title = ggplot2::element_text(size = 12, face = "bold"), legend.text = ggplot2::element_text(size = 10),
                   strip.text = ggplot2::element_text(size = 13, face = "bold"), strip.background = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(1, "lines"), panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
                   panel.grid.minor = ggplot2::element_blank(), panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
                   axis.text = ggplot2::element_text(size = 11, color = "black"), axis.title = ggplot2::element_text(size = 12, face = "bold"),
                   axis.line = ggplot2::element_line(color = "black", linewidth = 0.5), panel.background = ggplot2::element_rect(fill = "white", color = NA)) +
    ggplot2::labs(x = "Pseudotime", y = "Expression") +
    ggplot2::scale_y_continuous(labels = function(x) ifelse(x == floor(x), as.character(as.integer(x)), as.character(x)))
  if (length(genes) == 1) p <- p + ggplot2::theme(aspect.ratio = 0.6)
  return(list(plot = p, genes = genes))
}

# Display options applied on top of a stored trajectory or pseudotime ggplot (the stored plot is not modified)
style_trajectory_plot <- function(plot, hide_grid = FALSE, hide_axes = FALSE, hide_title = FALSE) {
  stopifnot(inherits(plot, "ggplot"))
  styled <- plot
  if (isTRUE(hide_grid)) styled <- styled + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())
  # Axis titles, numbers, ticks and lines removed together (UMAP coordinates carry no biological meaning)
  if (isTRUE(hide_axes)) styled <- styled + ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
                                                           axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank())
  if (isTRUE(hide_title)) styled <- styled + ggplot2::labs(title = NULL)
  return(styled)
}