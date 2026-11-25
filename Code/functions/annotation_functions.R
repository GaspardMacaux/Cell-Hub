# annotation_functions.R
# Automated cell-type annotation: SingleR with celldex references, and label transfer from a user-provided annotated
# Seurat reference (Seurat FindTransferAnchors / TransferData). Pure functions, no Shiny reactive code.

############################## Shared utilities ##############################

# Per-cluster summary of per-cell labels: majority label, agreement, unassigned fraction and confidence flag
# cell_df columns: cell, cluster, label, pruned_label (NA = unreliable), delta_next, max_score
summarize_annotation_by_cluster <- function(cell_df, min_agreement = 0.5, max_unassigned = 0.3) {
  rows <- lapply(split(cell_df, cell_df$cluster), function(d) {
    lab <- d$pruned_label
    tab <- sort(table(lab[!is.na(lab)]), decreasing = TRUE)
    n <- nrow(d)
    top <- if (length(tab) > 0) names(tab)[1] else NA_character_
    pct_top <- if (length(tab) > 0) as.numeric(tab[1]) / n else 0
    pct_unassigned <- mean(is.na(lab))
    data.frame(cluster = as.character(d$cluster[1]), n_cells = n, predicted_label = top,
               agreement_pct = round(100 * pct_top, 1),
               second_label = if (length(tab) > 1) names(tab)[2] else NA_character_,
               second_pct = if (length(tab) > 1) round(100 * as.numeric(tab[2]) / n, 1) else 0,
               unassigned_pct = round(100 * pct_unassigned, 1),
               median_score = round(stats::median(d$max_score, na.rm = TRUE), 3),
               median_delta = round(stats::median(d$delta_next, na.rm = TRUE), 3),
               confidence = if (is.na(top) || pct_top < min_agreement || pct_unassigned > max_unassigned) "Low" else "High",
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  # Numeric-aware cluster order (0, 1, 2, ..., 10 rather than 0, 1, 10, 2)
  out <- out[order(suppressWarnings(as.numeric(out$cluster)), out$cluster), ]
  rownames(out) <- NULL
  out
}

# Store annotation results as metadata columns <prefix>_label, _pruned_label, _score, _delta, _cluster_label
# (existing cluster columns are never modified)
add_annotation_metadata <- function(seurat_obj, result, prefix) {
  cd <- result$cells
  idx <- match(colnames(seurat_obj), cd$cell)
  pruned <- cd$pruned_label[idx]
  pruned[is.na(pruned)] <- "Unassigned"
  cluster_label <- result$clusters$predicted_label[match(cd$cluster, result$clusters$cluster)][idx]
  cluster_label[is.na(cluster_label)] <- "Unassigned"
  seurat_obj@meta.data[[paste0(prefix, "_label")]] <- cd$label[idx]
  seurat_obj@meta.data[[paste0(prefix, "_pruned_label")]] <- pruned
  seurat_obj@meta.data[[paste0(prefix, "_score")]] <- cd$max_score[idx]
  seurat_obj@meta.data[[paste0(prefix, "_delta")]] <- cd$delta_next[idx]
  seurat_obj@meta.data[[paste0(prefix, "_cluster_label")]] <- cluster_label
  seurat_obj
}

# New name for each cluster from a cluster summary: low-confidence clusters can keep their current name,
# and clusters sharing a label can be numbered (Neurons_1, Neurons_2) instead of being merged
build_annotation_cluster_names <- function(cluster_summary, keep_low_confidence = TRUE, make_unique = TRUE) {
  keep <- is.na(cluster_summary$predicted_label) | (keep_low_confidence & cluster_summary$confidence == "Low")
  new_names <- ifelse(keep, as.character(cluster_summary$cluster), cluster_summary$predicted_label)
  if (make_unique) {
    assigned <- !keep
    # Duplicates counted over all new names, so an assigned label equal to a kept cluster name is numbered too
    counts <- table(new_names)
    shared <- assigned & new_names %in% names(counts)[counts > 1]
    rank_in_label <- stats::ave(seq_along(new_names), ifelse(assigned, new_names, paste0("\r", new_names)), FUN = seq_along)
    new_names[shared] <- paste0(new_names[shared], "_", rank_in_label[shared])
    # A generated name may still collide with an existing kept name (e.g. a cluster already called "Neurons_1")
    new_names <- make.unique(new_names, sep = "_")
  }
  stats::setNames(new_names, cluster_summary$cluster)
}

# Log-normalized expression matrix of an assay; Seurat v5 split layers are joined first (GetAssayData fails with multiple layers)
extract_normalized_matrix <- function(seurat_obj, assay = DefaultAssay(seurat_obj)) {
  if (!assay %in% names(seurat_obj@assays)) stop(sprintf("Assay '%s' not found", assay))
  if (inherits(seurat_obj[[assay]], "Assay5")) seurat_obj <- JoinLayers(seurat_obj, assay = assay)
  if (!"data" %in% SeuratObject::Layers(seurat_obj[[assay]])) stop(sprintf("No normalized 'data' layer in assay '%s': run normalization first", assay))
  LayerData(seurat_obj, assay = assay, layer = "data")
}

############################## SingleR (celldex references) ##############################

# celldex references available per species (names = labels shown in the UI, values = celldex function names)
get_singler_references <- function(species) {
  refs <- list(
    human = c("Human Primary Cell Atlas (broad)" = "HumanPrimaryCellAtlasData",
              "Blueprint/ENCODE (immune, stromal)" = "BlueprintEncodeData",
              "Monaco (human immune, PBMC)" = "MonacoImmuneData",
              "DICE (human immune)" = "DatabaseImmuneCellExpressionData",
              "Novershtern (hematopoietic)" = "NovershternHematopoieticData"),
    mouse = c("Mouse RNA-seq (broad)" = "MouseRNAseqData",
              "ImmGen (mouse immune)" = "ImmGenData")
  )
  if (is.null(species) || !species %in% names(refs)) return(character(0))
  refs[[species]]
}

# Load a celldex reference, cached as RDS so that it is downloaded only once (or pre-cached in the Docker image)
load_singler_reference <- function(ref_name, cache_dir) {
  allowed <- unname(c(get_singler_references("human"), get_singler_references("mouse")))
  if (!ref_name %in% allowed) stop(sprintf("Unknown SingleR reference: %s", ref_name))
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(cache_dir, paste0(ref_name, ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))
  if (!requireNamespace("celldex", quietly = TRUE)) stop("Package 'celldex' is not installed and the reference is not cached")
  # celldex returns gene symbols by default (ensembl = FALSE); the first call downloads the data
  ref <- getExportedValue("celldex", ref_name)()
  saveRDS(ref, cache_file)
  ref
}

# SingleR on a log-normalized genes x cells matrix
# mode = "cell": one label per cell then a majority vote per cluster; mode = "cluster": one label per cluster (aggregated, faster)
singler_annotate_matrix <- function(test_matrix, reference, label_column = "label.main", clusters = NULL, mode = c("cell", "cluster"),
                                    min_common_genes = 500, min_agreement = 0.5, max_unassigned = 0.3) {
  mode <- match.arg(mode)
  if (!requireNamespace("SingleR", quietly = TRUE)) stop("Package 'SingleR' is not installed")
  if (!label_column %in% colnames(SummarizedExperiment::colData(reference))) stop(sprintf("Label column '%s' not found in the reference", label_column))
  common <- intersect(rownames(test_matrix), rownames(reference))
  if (length(common) < min_common_genes) {
    stop(sprintf("Only %d genes shared with the reference (minimum %d): check that species and gene symbols match", length(common), min_common_genes))
  }
  test_matrix <- test_matrix[common, , drop = FALSE]
  reference <- reference[common, ]
  labels <- as.character(SummarizedExperiment::colData(reference)[[label_column]])
  clusters <- if (is.null(clusters)) rep("all", ncol(test_matrix)) else as.character(clusters)
  if (length(clusters) != ncol(test_matrix)) stop("clusters must have one value per cell")
  set.seed(42)
  if (mode == "cell") {
    pred <- SingleR::SingleR(test = test_matrix, ref = reference, labels = labels, de.method = "classic")
    idx <- seq_len(ncol(test_matrix))
  } else {
    pred <- SingleR::SingleR(test = test_matrix, ref = reference, labels = labels, clusters = clusters, de.method = "classic")
    idx <- match(clusters, rownames(pred))
  }
  max_score <- apply(pred$scores, 1, max)
  cell_df <- data.frame(cell = colnames(test_matrix), cluster = clusters,
                        label = pred$labels[idx], pruned_label = pred$pruned.labels[idx],
                        delta_next = pred$delta.next[idx], max_score = max_score[idx], stringsAsFactors = FALSE)
  list(cells = cell_df,
       clusters = summarize_annotation_by_cluster(cell_df, min_agreement, max_unassigned),
       params = list(method = "SingleR", label_column = label_column, mode = mode, n_common_genes = length(common),
                     min_agreement = min_agreement, max_unassigned = max_unassigned))
}

# SingleR on a Seurat object (log-normalized data of the chosen assay)
run_singler_annotation <- function(seurat_obj, reference, label_column = "label.main", clusters = NULL, mode = "cell",
                                   assay = DefaultAssay(seurat_obj), min_agreement = 0.5, max_unassigned = 0.3) {
  test_matrix <- extract_normalized_matrix(seurat_obj, assay)
  result <- singler_annotate_matrix(test_matrix, reference, label_column = label_column, clusters = clusters, mode = mode,
                                    min_agreement = min_agreement, max_unassigned = max_unassigned)
  result$params$assay <- assay
  result
}

############################## Label transfer from a user reference ##############################

# Read an annotated reference: a Seurat object, or a Cell-Hub list with a $seurat element
read_reference_object <- function(path) {
  obj <- readRDS(path)
  if (is.list(obj) && !inherits(obj, "Seurat") && inherits(obj$seurat, "Seurat")) obj <- obj$seurat
  if (!inherits(obj, "Seurat")) stop("The reference file must contain a Seurat object")
  if (!"RNA" %in% names(obj@assays)) stop("The reference must contain an 'RNA' assay")
  obj
}

# Categorical metadata columns usable as labels (between 2 and max_levels distinct values)
get_reference_label_columns <- function(reference, max_levels = 300) {
  meta <- reference@meta.data
  cols <- colnames(meta)[vapply(meta, function(x) {
    if (!(is.character(x) || is.factor(x))) return(FALSE)
    n <- length(unique(x[!is.na(x)]))
    n >= 2 && n <= max_levels
  }, logical(1))]
  setdiff(cols, c("orig.ident"))
}

# Downsample the reference per label (seed 42), keep only the RNA assay, join layers, normalize and find variable genes if needed
prepare_reference_object <- function(reference, label_column, max_cells_per_label = 500) {
  if (!label_column %in% colnames(reference@meta.data)) stop(sprintf("Label column '%s' not found in the reference", label_column))
  DefaultAssay(reference) <- "RNA"
  labels <- as.character(reference@meta.data[[label_column]])
  keep <- !is.na(labels) & nzchar(labels)
  if (max_cells_per_label > 0) {
    set.seed(42)
    kept_idx <- unlist(lapply(split(which(keep), labels[keep]), function(i) if (length(i) > max_cells_per_label) sample(i, max_cells_per_label) else i))
    keep <- seq_len(ncol(reference)) %in% kept_idx
  }
  reference <- subset(reference, cells = colnames(reference)[keep])
  reference <- DietSeurat(reference, assays = "RNA")
  if (inherits(reference[["RNA"]], "Assay5")) reference <- JoinLayers(reference, assay = "RNA")
  if (!"data" %in% SeuratObject::Layers(reference[["RNA"]])) reference <- NormalizeData(reference, verbose = FALSE)
  if (length(VariableFeatures(reference)) == 0) reference <- FindVariableFeatures(reference, nfeatures = 2000, verbose = FALSE)
  reference
}

# Seurat label transfer (PCA projection anchors) from a prepared reference to the query's RNA assay
# Per-cell confidence = prediction.score.max; cells below min_score are reported as unassigned
run_label_transfer <- function(query, reference, label_column, clusters = NULL, dims = 30, min_score = 0.5,
                               min_agreement = 0.5, max_unassigned = 0.3, min_common_genes = 500) {
  if (!"RNA" %in% names(query@assays)) stop("The dataset must contain an 'RNA' assay")
  DefaultAssay(query) <- "RNA"
  if (inherits(query[["RNA"]], "Assay5")) query <- JoinLayers(query, assay = "RNA")
  if (!"data" %in% SeuratObject::Layers(query[["RNA"]])) query <- NormalizeData(query, verbose = FALSE)
  common <- intersect(rownames(query), rownames(reference))
  if (length(common) < min_common_genes) {
    stop(sprintf("Only %d genes shared with the reference (minimum %d): check that species and gene symbols match", length(common), min_common_genes))
  }
  features <- intersect(VariableFeatures(reference), common)
  if (length(features) < 200) stop(sprintf("Only %d reference variable genes found in the dataset (minimum 200)", length(features)))
  dims <- seq_len(min(dims, length(features) - 1, ncol(reference) - 1))
  clusters <- if (is.null(clusters)) rep("all", ncol(query)) else as.character(clusters)
  set.seed(42)
  anchors <- FindTransferAnchors(reference = reference, query = query, normalization.method = "LogNormalize",
                                 reference.assay = "RNA", query.assay = "RNA", features = features, dims = dims, verbose = FALSE)
  n_anchors <- nrow(anchors@anchors)
  # k.weight must stay below the number of anchors, otherwise TransferData stops
  k_weight <- max(5, min(50, n_anchors - 1))
  refdata <- stats::setNames(as.character(reference@meta.data[[label_column]]), colnames(reference))
  pred <- TransferData(anchorset = anchors, refdata = refdata, dims = dims, k.weight = k_weight, verbose = FALSE)
  score_cols <- setdiff(grep("^prediction\\.score\\.", colnames(pred), value = TRUE), "prediction.score.max")
  scores <- as.matrix(pred[, score_cols, drop = FALSE])
  second_best <- if (ncol(scores) > 1) apply(scores, 1, function(s) sort(s, decreasing = TRUE)[2]) else rep(0, nrow(scores))
  label <- as.character(pred$predicted.id)
  max_score <- pred$prediction.score.max
  cell_df <- data.frame(cell = rownames(pred), cluster = clusters[match(rownames(pred), colnames(query))],
                        label = label, pruned_label = ifelse(max_score >= min_score, label, NA_character_),
                        delta_next = max_score - second_best, max_score = max_score, stringsAsFactors = FALSE)
  list(cells = cell_df,
       clusters = summarize_annotation_by_cluster(cell_df, min_agreement, max_unassigned),
       params = list(method = "Seurat label transfer", label_column = label_column, n_common_genes = length(common),
                     n_features = length(features), dims = max(dims), n_anchors = n_anchors, k_weight = k_weight,
                     n_reference_cells = ncol(reference), min_score = min_score,
                     min_agreement = min_agreement, max_unassigned = max_unassigned))
}

# Expression assay for annotation: RNA holds all genes, an "integrated" assay only the corrected variable genes
annotation_assay_merge <- function(obj) {
  if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
}

