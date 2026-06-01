# workspace_management_newarch.R



cleanTempDirectories <- function(prefix_patterns = c("single_dataset", "spatial", "multiome")) {
  # Clean temporary directories created by the app
  
  temp_dirs <- list.dirs(tempdir(), full.names = TRUE, recursive = FALSE)
  
  # Filter directories matching our patterns
  app_temp_dirs <- character(0)
  for (pattern in prefix_patterns) {
    matching_dirs <- temp_dirs[grepl(pattern, basename(temp_dirs))]
    app_temp_dirs <- c(app_temp_dirs, matching_dirs)
  }
  
  # Remove directories
  removed_count <- 0
  for (dir in app_temp_dirs) {
    if (dir.exists(dir)) {
      unlink(dir, recursive = TRUE)
      message(paste("Removed temp directory:", basename(dir)))
      removed_count <- removed_count + 1
    }
  }
  
  message(paste("Removed", removed_count, "temporary directories"))
}

cleanWorkspace <- function(module = "all", custom_patterns = NULL) {
  # General workspace cleaning function
  # Args:
  #   module: "single", "spatial", "multiome", or "all"
  #   custom_patterns: Additional patterns to clean
  
  message(paste("Cleaning workspace for module:", module))
  
  patterns <- switch(module,
                     "single" = c("single_dataset", "_plot$", "_data$", "gene_list", "subset_", "normalized_", "scaled_"),
                     "spatial" = c("spatial_", "tissue_", "spot_"),
                     "multiome" = c("multiome_", "atac_", "rna_", "peak_"),
                     # ✅ ADD: Support for multiple datasets
                     "multiple" = c("multiple_datasets_", "_merge$", "_plot_merge$", "clustering_plot_merge", 
                                    "feature_plot_merge", "vln_plot_merge", "dot_plot_merge", "ridge_plot_merge", 
                                    "heatmap_plot_multidataset", "seurat_objects", "data_loaded", "merged_gene_tables"),
                     "all" = c("_plot$", "_data$", "_seurat$", "_object$", "gene_list", "subset_", "normalized_", 
                               "scaled_", "feature_", "neighbors_", "clustering_", "spatial_", "multiome_", 
                               "multiple_datasets_", "_merge$")
  )
  
  if (!is.null(custom_patterns)) {
    patterns <- c(patterns, custom_patterns)
  }
  
  
  # Clean temp directories
  temp_patterns <- switch(module,
                          "multiple" = "multiple_datasets",
                          "single" = "single_dataset", 
                          "spatial" = "spatial",
                          "multiome" = "multiome",
                          c("single_dataset", "multiple_datasets", "spatial", "multiome")
  )
  
  cleanTempDirectories(prefix_patterns = temp_patterns)
  
  gc()
  message(paste("Workspace cleaning completed for module:", module))
}

createTempDirectory <- function(prefix = "analysis") {
  # Create a timestamped temporary directory
  
  temp_dir <- file.path(tempdir(), paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  message(paste("Created temp directory:", temp_dir))
  
  return(temp_dir)
}


#' Collect cluster composition statistics from a Seurat object
#'
#' Extracts cluster names, nuclei counts per cluster, and per-dataset breakdown
#' for multi-dataset objects. Returns rows formatted for the parameters export table.
#'
#' @param seurat_obj Seurat object after clustering
#' @param module Character string: "single" or "multiple"
#' @return List of parameter rows (same format as collect_analysis_params internals)
collect_cluster_stats <- function(seurat_obj, module = "single") {
  params <- list()
  
  if (is.null(seurat_obj)) return(params)
  
  # Resolve cluster identities
  if (module == "multiple" && "ClusterIdents" %in% colnames(seurat_obj@meta.data)) {
    cluster_vec <- seurat_obj@meta.data[["ClusterIdents"]]
  } else {
    cluster_vec <- as.character(Idents(seurat_obj))
  }
  
  cluster_table <- sort(table(cluster_vec), decreasing = TRUE)
  total_nuclei  <- ncol(seurat_obj)
  n_clusters    <- length(cluster_table)
  cluster_names <- names(cluster_table)
  meta          <- seurat_obj@meta.data
  has_genes     <- "nFeature_RNA" %in% colnames(meta)
  has_mito      <- "percent.mt"   %in% colnames(meta)
  
  # ── Total nuclei and cluster count ──────────────────────────────────────────
  params[["total_nuclei"]] <- list(
    category    = "Cluster Composition",
    parameter   = "Total nuclei (post-clustering)",
    value       = total_nuclei,
    description = "Total number of nuclei in the object at time of export"
  )
  params[["n_clusters"]] <- list(
    category    = "Cluster Composition",
    parameter   = "Number of clusters",
    value       = n_clusters,
    description = "Number of distinct clusters identified"
  )
  params[["cluster_names_all"]] <- list(
    category    = "Cluster Composition",
    parameter   = "Cluster names (all)",
    value       = paste(cluster_names, collapse = " | "),
    description = "All cluster names in descending order of size"
  )
  
  # ── Per-cluster count + median QC metrics ────────────────────────────────────
  for (cl_name in cluster_names) {
    safe_key     <- paste0("cluster_", gsub("[^a-zA-Z0-9]", "_", cl_name))
    cells_in_cl  <- which(cluster_vec == cl_name)
    n_in_cluster <- length(cells_in_cl)
    pct_in_cluster <- round(n_in_cluster / total_nuclei * 100, 1)
    
    med_genes <- if (has_genes) round(median(meta$nFeature_RNA[cells_in_cl], na.rm = TRUE), 0) else "N/A"
    med_mito  <- if (has_mito)  round(median(meta$percent.mt[cells_in_cl],   na.rm = TRUE), 2) else "N/A"
    
    params[[safe_key]] <- list(
      category    = "Cluster Composition",
      parameter   = paste0("Cluster: ", cl_name),
      value       = n_in_cluster,
      description = paste0(pct_in_cluster, "% of total | ",
                           "median genes: ", med_genes, " | ",
                           "median %mito: ", med_mito)
    )
  }
  
  # ── Per-dataset cluster breakdown (multiple only) ────────────────────────────
  if (module == "multiple" && "dataset" %in% colnames(meta)) {
    dataset_vec   <- meta[["dataset"]]
    cross_tab     <- table(cluster_vec, dataset_vec)
    dataset_names <- colnames(cross_tab)
    
    for (ds_name in dataset_names) {
      n_ds    <- sum(dataset_vec == ds_name)
      safe_ds <- paste0("dataset_total_", gsub("[^a-zA-Z0-9]", "_", ds_name))
      params[[safe_ds]] <- list(
        category    = paste0("Dataset Composition - ", ds_name),
        parameter   = paste0("[", ds_name, "] Total nuclei"),
        value       = n_ds,
        description = paste0("Total nuclei in dataset '", ds_name,
                             "' (", round(n_ds / total_nuclei * 100, 1), "% of merged object)")
      )
    }
    
    for (cl_name in cluster_names) {
      for (ds_name in dataset_names) {
        n_cell   <- as.integer(cross_tab[cl_name, ds_name])
        n_clust  <- as.integer(cluster_table[[cl_name]])
        safe_key <- paste0("ds_", gsub("[^a-zA-Z0-9]", "_", ds_name),
                           "_cl_", gsub("[^a-zA-Z0-9]", "_", cl_name))
        params[[safe_key]] <- list(
          category    = paste0("Dataset Composition - ", ds_name),
          parameter   = paste0("[", ds_name, "] Cluster: ", cl_name),
          value       = n_cell,
          description = paste0(round(n_cell / n_clust * 100, 1),
                               "% of cluster '", cl_name, "' comes from ", ds_name)
        )
      }
    }
  }
  
  return(params)
}

#' Collect analysis parameters for export as CSV
#'
#' Builds a structured data.frame summarizing all key analytical parameters
#' used during the session. Captures QC thresholds, pre/post-QC stats,
#' normalization, dimensionality reduction, clustering, integration settings,
#' and cluster composition.
#'
#' @param input Shiny input object
#' @param module Character string: "single" or "multiple"
#' @param qc_stats List with pre/post QC metrics (NULL if QC not yet run)
#' @param seurat_obj Seurat object for cluster composition (NULL if not yet clustered)
#' @return data.frame with columns: category, parameter, value, description
collect_analysis_params <- function(input, module = "single", 
                                    qc_stats = NULL, seurat_obj = NULL) {
  
  params <- list()
  
  algo_labels <- c(
    "1" = "Original Louvain",
    "2" = "Louvain with Multilevel Refinement",
    "3" = "SLM Algorithm"
  )
  
  # ── QC THRESHOLDS ────────────────────────────────────────────────────────────
  params[["qc_nfeature_min"]] <- list(
    category    = "QC Thresholds",
    parameter   = "Min nFeature_RNA",
    value       = input$nFeature_range[1],
    description = "Minimum number of unique genes per nucleus (lower QC cutoff)"
  )
  params[["qc_nfeature_max"]] <- list(
    category    = "QC Thresholds",
    parameter   = "Max nFeature_RNA",
    value       = input$nFeature_range[2],
    description = "Maximum number of unique genes per nucleus (upper QC cutoff)"
  )
  params[["qc_mt_max"]] <- list(
    category    = "QC Thresholds",
    parameter   = "Max mitochondrial %",
    value       = input$percent.mt_max,
    description = "Maximum percentage of mitochondrial reads per nucleus"
  )
  
  # ── PRE/POST QC STATISTICS ───────────────────────────────────────────────────
  if (!is.null(qc_stats)) {
    
    if (module == "single") {
      params[["qc_n_before"]] <- list(
        category    = "QC Statistics",
        parameter   = "Nuclei before QC",
        value       = qc_stats$n_before,
        description = "Total number of nuclei before QC filtering"
      )
      params[["qc_n_after"]] <- list(
        category    = "QC Statistics",
        parameter   = "Nuclei after QC",
        value       = qc_stats$n_after,
        description = "Total number of nuclei retained after QC filtering"
      )
      params[["qc_pct_retained"]] <- list(
        category    = "QC Statistics",
        parameter   = "% nuclei retained",
        value       = round(qc_stats$n_after / qc_stats$n_before * 100, 1),
        description = "Percentage of nuclei passing QC filters"
      )
      params[["qc_median_genes_before"]] <- list(
        category    = "QC Statistics",
        parameter   = "Median genes/nucleus (before QC)",
        value       = qc_stats$median_genes_before,
        description = "Median nFeature_RNA before QC filtering"
      )
      params[["qc_median_genes_after"]] <- list(
        category    = "QC Statistics",
        parameter   = "Median genes/nucleus (after QC)",
        value       = qc_stats$median_genes_after,
        description = "Median nFeature_RNA after QC filtering"
      )
      params[["qc_median_counts_before"]] <- list(
        category    = "QC Statistics",
        parameter   = "Median UMI counts/nucleus (before QC)",
        value       = qc_stats$median_counts_before,
        description = "Median nCount_RNA before QC filtering"
      )
      params[["qc_median_counts_after"]] <- list(
        category    = "QC Statistics",
        parameter   = "Median UMI counts/nucleus (after QC)",
        value       = qc_stats$median_counts_after,
        description = "Median nCount_RNA after QC filtering"
      )
      params[["qc_median_mt_before"]] <- list(
        category    = "QC Statistics",
        parameter   = "Median %mito (before QC)",
        value       = qc_stats$median_mt_before,
        description = "Median mitochondrial percentage before QC filtering"
      )
      params[["qc_median_mt_after"]] <- list(
        category    = "QC Statistics",
        parameter   = "Median %mito (after QC)",
        value       = qc_stats$median_mt_after,
        description = "Median mitochondrial percentage after QC filtering"
      )
      
    } else if (module == "multiple") {
      dataset_names <- qc_stats$dataset_names
      per_dataset   <- qc_stats$per_dataset
      
      for (ds_name in dataset_names) {
        ds <- per_dataset[[ds_name]]
        if (is.null(ds)) next
        safe_key <- gsub("[^a-zA-Z0-9]", "_", ds_name)
        
        params[[paste0(safe_key, "_n_before")]] <- list(
          category    = paste0("QC Statistics - ", ds_name),
          parameter   = paste0("[", ds_name, "] Nuclei before QC"),
          value       = ds$n_before,
          description = paste("Total nuclei before QC filtering for dataset:", ds_name)
        )
        params[[paste0(safe_key, "_n_after")]] <- list(
          category    = paste0("QC Statistics - ", ds_name),
          parameter   = paste0("[", ds_name, "] Nuclei after QC"),
          value       = ds$n_after,
          description = paste("Nuclei retained after QC filtering for dataset:", ds_name)
        )
        params[[paste0(safe_key, "_pct")]] <- list(
          category    = paste0("QC Statistics - ", ds_name),
          parameter   = paste0("[", ds_name, "] % nuclei retained"),
          value       = round(ds$n_after / ds$n_before * 100, 1),
          description = paste("Percentage of nuclei passing QC for dataset:", ds_name)
        )
        params[[paste0(safe_key, "_median_genes_before")]] <- list(
          category    = paste0("QC Statistics - ", ds_name),
          parameter   = paste0("[", ds_name, "] Median genes/nucleus (before QC)"),
          value       = ifelse(is.na(ds$median_genes_before), "N/A", ds$median_genes_before),
          description = paste("Median nFeature_RNA before QC for dataset:", ds_name)
        )
        params[[paste0(safe_key, "_median_genes_after")]] <- list(
          category    = paste0("QC Statistics - ", ds_name),
          parameter   = paste0("[", ds_name, "] Median genes/nucleus (after QC)"),
          value       = ds$median_genes_after,
          description = paste("Median nFeature_RNA after QC for dataset:", ds_name)
        )
        params[[paste0(safe_key, "_median_counts_before")]] <- list(
          category    = paste0("QC Statistics - ", ds_name),
          parameter   = paste0("[", ds_name, "] Median UMI counts/nucleus (before QC)"),
          value       = ifelse(is.na(ds$median_counts_before), "N/A", ds$median_counts_before),
          description = paste("Median nCount_RNA before QC for dataset:", ds_name)
        )
        params[[paste0(safe_key, "_median_counts_after")]] <- list(
          category    = paste0("QC Statistics - ", ds_name),
          parameter   = paste0("[", ds_name, "] Median UMI counts/nucleus (after QC)"),
          value       = ds$median_counts_after,
          description = paste("Median nCount_RNA after QC for dataset:", ds_name)
        )
        params[[paste0(safe_key, "_median_mt_before")]] <- list(
          category    = paste0("QC Statistics - ", ds_name),
          parameter   = paste0("[", ds_name, "] Median %mito (before QC)"),
          value       = ifelse(is.na(ds$median_mt_before), "N/A", ds$median_mt_before),
          description = paste("Median %mito before QC for dataset:", ds_name)
        )
        params[[paste0(safe_key, "_median_mt_after")]] <- list(
          category    = paste0("QC Statistics - ", ds_name),
          parameter   = paste0("[", ds_name, "] Median %mito (after QC)"),
          value       = ds$median_mt_after,
          description = paste("Median %mito after QC for dataset:", ds_name)
        )
      }
    }
  }
  
  # ── NORMALIZATION & FEATURE SELECTION ────────────────────────────────────────
  if (module == "single") {
    params[["norm_method"]] <- list(
      category    = "Normalization",
      parameter   = "Normalization method",
      value       = input$normalization_method_single,
      description = "Normalization method applied to raw counts"
    )
    params[["scale_factor"]] <- list(
      category    = "Normalization",
      parameter   = "Scale factor",
      value       = input$scale_factor,
      description = "Scale factor used in LogNormalize (ignored for SCTransform)"
    )
    params[["num_var_features"]] <- list(
      category    = "Normalization",
      parameter   = "Variable features (n)",
      value       = input$num_var_features,
      description = "Number of highly variable features selected for PCA"
    )
  } else if (module == "multiple") {
    norm_val <- tryCatch(input$normalization_method_merge, error = function(e) "N/A")
    params[["norm_method"]] <- list(
      category    = "Normalization",
      parameter   = "Normalization method",
      value       = norm_val,
      description = "Normalization method applied to each dataset before integration"
    )
    params[["scale_factor"]] <- list(
      category    = "Normalization",
      parameter   = "Scale factor",
      value       = input$scale_factor,
      description = "Scale factor used in LogNormalize (ignored for SCTransform)"
    )
    params[["num_var_features"]] <- list(
      category    = "Normalization",
      parameter   = "Variable features (n)",
      value       = input$num_var_features,
      description = "Number of highly variable features selected per dataset"
    )
  }
  
  # ── DIMENSIONALITY REDUCTION & CLUSTERING ────────────────────────────────────
  if (module == "single") {
    params[["n_pcs"]] <- list(
      category    = "Dimensionality Reduction",
      parameter   = "Number of PCs (FindNeighbors / UMAP)",
      value       = input$dimension_1,
      description = "Number of principal components used for FindNeighbors and RunUMAP"
    )
    params[["resolution"]] <- list(
      category    = "Clustering",
      parameter   = "Clustering resolution",
      value       = input$resolution_step1,
      description = "Resolution parameter for FindClusters (higher = more clusters)"
    )
    params[["algorithm"]] <- list(
      category    = "Clustering",
      parameter   = "Clustering algorithm",
      value       = algo_labels[as.character(input$algorithm_select)],
      description = "Graph-based clustering algorithm used in FindClusters"
    )
  } else if (module == "multiple") {
    params[["n_pcs"]] <- list(
      category    = "Dimensionality Reduction",
      parameter   = "Number of PCs (FindNeighbors / UMAP)",
      value       = input$dimension_2,
      description = "Number of principal components used for FindNeighbors and RunUMAP"
    )
    params[["resolution"]] <- list(
      category    = "Clustering",
      parameter   = "Clustering resolution",
      value       = input$resolution_step2,
      description = "Resolution parameter for FindClusters (higher = more clusters)"
    )
    params[["algorithm"]] <- list(
      category    = "Clustering",
      parameter   = "Clustering algorithm",
      value       = tryCatch(algo_labels[as.character(input$algorithm_select)], 
                             error = function(e) "N/A"),
      description = "Graph-based clustering algorithm used in FindClusters"
    )
  }
  
  # ── INTEGRATION (multiple only) ───────────────────────────────────────────────
  if (module == "multiple" && !is.null(seurat_obj)) {
    
    integ_info <- extract_integration_info(seurat_obj)
    
    params[["n_datasets"]] <- list(
      category    = "Integration",
      parameter   = "Number of datasets integrated",
      value       = length(integ_info$dataset_names),
      description = "Total number of datasets merged or integrated"
    )
    params[["dataset_names_list"]] <- list(
      category    = "Integration",
      parameter   = "Dataset names",
      value       = paste(integ_info$dataset_names, collapse = " | "),
      description = "Names of all integrated datasets"
    )
    params[["integration_method"]] <- list(
      category    = "Integration",
      parameter   = "Integration method",
      value       = integ_info$integration_method,
      description = "Method used to combine datasets (simple_merge, standard, or harmony)"
    )
    if (!is.null(integ_info$harmony_vars)) {
      params[["harmony_vars"]] <- list(
        category    = "Integration",
        parameter   = "Harmony - Batch correction variable",
        value       = integ_info$harmony_vars,
        description = "Metadata variable used for Harmony batch correction"
      )
    }
    
    # Also override normalization from object commands if not found via input
    params[["norm_method"]]$value <- integ_info$normalization_method
  }
  
  # ── DOUBLETFINDER (single only) ───────────────────────────────────────────────
  if (module == "single") {
    params[["doublet_rate"]] <- list(
      category    = "DoubletFinder",
      parameter   = "Expected doublet rate (%)",
      value       = input$doublet_rate,
      description = "Expected doublet rate (~1% per 1000 cells loaded)"
    )
    params[["doublet_pcs"]] <- list(
      category    = "DoubletFinder",
      parameter   = "PCs used",
      value       = input$pc_use,
      description = "Number of PCs used in DoubletFinder paramSweep"
    )
    params[["doublet_pN"]] <- list(
      category    = "DoubletFinder",
      parameter   = "pN",
      value       = input$pN_value,
      description = "Proportion of artificial doublets generated (pN parameter)"
    )
    params[["doublet_pK"]] <- list(
      category    = "DoubletFinder",
      parameter   = "pK",
      value       = input$pK_value,
      description = "Neighborhood size for doublet scoring (dataset-specific, optimize with pK sweep)"
    )
  }
  
  # ── CLUSTER COMPOSITION (from Seurat object at export time) ──────────────────
  if (!is.null(seurat_obj)) {
    cluster_params <- tryCatch(
      collect_cluster_stats(seurat_obj, module = module),
      error = function(e) {
        message("Warning: could not extract cluster stats: ", e$message)
        list()
      }
    )
    params <- c(params, cluster_params)
  }
  
  # ── REPRODUCIBILITY ───────────────────────────────────────────────────────────
  params[["seed"]] <- list(
    category    = "Reproducibility",
    parameter   = "Random seed",
    value       = 42,
    description = "Fixed seed used in RunPCA, RunUMAP, FindNeighbors, FindClusters"
  )
  params[["export_date"]] <- list(
    category    = "Reproducibility",
    parameter   = "Export timestamp",
    value       = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    description = "Date and time of parameter export"
  )
  
  # Convert to data.frame - coerce all values to character to handle
  # mixed types (NA, NULL, length-0 vectors) without rbind failures
  params_df <- do.call(rbind, lapply(params, function(x) {
    as.data.frame(
      lapply(x, function(v) {
        if (is.null(v) || length(v) == 0) return(NA_character_)
        as.character(v[1])  # [1] to guard against accidental length > 1
      }),
      stringsAsFactors = FALSE
    )
  }))
  rownames(params_df) <- NULL
  params_df <- params_df[order(params_df$category, params_df$parameter), ]
  
  return(params_df)
}

#' Extract integration and normalization info from a Seurat object
#'
#' Reads metadata columns and Seurat command history to recover
#' integration method, normalization method, and dataset names
#' without relying on reactive tracking.
#'
#' @param seurat_obj Seurat object after integration
#' @return List with: integration_method, normalization_method, 
#'         harmony_vars, dataset_names
extract_integration_info <- function(seurat_obj) {
  
  result <- list(
    integration_method   = "unknown",
    normalization_method = "unknown",
    harmony_vars         = NULL,
    dataset_names        = character(0)
  )
  
  if (is.null(seurat_obj)) return(result)
  
  # ── Integration method from metadata column ──────────────────────────────────
  # harmony stores "harmony" in $integration_method
  # simple merge stores "simple_merge" in $merge_method
  if ("integration_method" %in% colnames(seurat_obj@meta.data)) {
    result$integration_method <- unique(seurat_obj@meta.data$integration_method)[1]
  } else if ("merge_method" %in% colnames(seurat_obj@meta.data)) {
    result$integration_method <- unique(seurat_obj@meta.data$merge_method)[1]
  }
  
  # ── Dataset names from $dataset column ──────────────────────────────────────
  if ("dataset" %in% colnames(seurat_obj@meta.data)) {
    result$dataset_names <- unique(as.character(seurat_obj@meta.data$dataset))
  } else if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
    result$dataset_names <- unique(as.character(seurat_obj@meta.data$orig.ident))
  }
  
  # ── Normalization method from Seurat command history (@commands slot) ────────
  # Seurat automatically logs all commands with their parameters
  if (length(seurat_obj@commands) > 0) {
    cmd_names <- names(seurat_obj@commands)
    
    # NormalizeData command stores normalization.method param
    normalize_cmd <- grep("^NormalizeData", cmd_names, value = TRUE)
    if (length(normalize_cmd) > 0) {
      cmd <- seurat_obj@commands[[normalize_cmd[1]]]
      norm_method <- tryCatch(
        cmd@params$normalization.method,
        error = function(e) NULL
      )
      if (!is.null(norm_method)) {
        result$normalization_method <- norm_method
      } else {
        result$normalization_method <- "LogNormalize"  # Seurat default
      }
    }
    
    # SCTransform command = SCTransform normalization
    sct_cmd <- grep("^SCTransform", cmd_names, value = TRUE)
    if (length(sct_cmd) > 0) {
      result$normalization_method <- "SCTransform"
    }
    
    # RunHarmony command stores group.by.vars param
    harmony_cmd <- grep("^RunHarmony", cmd_names, value = TRUE)
    if (length(harmony_cmd) > 0) {
      cmd <- seurat_obj@commands[[harmony_cmd[1]]]
      harmony_vars <- tryCatch(
        cmd@params$group.by.vars,
        error = function(e) NULL
      )
      if (!is.null(harmony_vars)) {
        result$harmony_vars <- paste(harmony_vars, collapse = ", ")
      }
    }
  }
  
  return(result)
}

# ============================================================================
# find_stored_colors()
#
# Robustly scans an S4 object's misc and options slots for a named hex-color
# vector that plausibly represents cluster colors.
#
# Strategy (in order of confidence):
#   1. Look for a slot named exactly "cluster_colors" in @misc — highest trust
#   2. Scan all other slots for named hex vectors with high overlap with the
#      provided cluster names (overlap >= 0.6 both ways to avoid false positives)
#   3. Among candidates, pick the one with the highest overlap score
#
# Args:
#   obj           : S4 object (Seurat, CellChat, or similar)
#   cluster_names : character vector of expected cluster/cell-type names.
#                   If NULL, returns best candidate without overlap filtering.
#   min_overlap   : minimum fraction of cluster_names that must be found in the
#                   candidate AND fraction of candidate names that must be in
#                   cluster_names (bidirectional check). Default 0.6.
#
# Returns:
#   A named character vector of hex colors, or NULL if nothing credible found.
# ============================================================================
find_stored_colors <- function(obj, cluster_names = NULL, min_overlap = 0.6) {
  hex_pattern <- "^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$"
  
  candidates <- list()
  
  available_slots <- tryCatch(slotNames(obj), error = function(e) character(0))
  
  if ("misc" %in% available_slots) {
    tryCatch({
      misc_val <- slot(obj, "misc")
      if (is.list(misc_val) && length(misc_val) > 0) {
        for (nm in names(misc_val)) {
          candidates[[paste0("misc::", nm)]] <- misc_val[[nm]]
        }
      }
    }, error = function(e) NULL)
  }
  
  if ("options" %in% available_slots) {
    tryCatch({
      opts <- slot(obj, "options")
      if (is.list(opts) && length(opts) > 0) {
        for (nm in names(opts)) {
          candidates[[paste0("options::", nm)]] <- opts[[nm]]
        }
      }
    }, error = function(e) NULL)
  }
  
  if (length(candidates) == 0) return(NULL)
  
  # Fast path 1: canonical Cell-Hub slot in Seurat @misc
  canonical_misc <- "misc::cluster_colors"
  if (canonical_misc %in% names(candidates)) {
    v <- candidates[[canonical_misc]]
    if (is.character(v) && !is.null(names(v)) && all(grepl(hex_pattern, v))) {
      if (is.null(cluster_names) || any(cluster_names %in% names(v))) {
        return(v)
      }
    }
  }
  
  # Fast path 2: canonical Cell-Hub slot in CellChat @options
  canonical_options <- "options::cluster_colors"
  if (canonical_options %in% names(candidates)) {
    v <- candidates[[canonical_options]]
    if (is.character(v) && !is.null(names(v)) && all(grepl(hex_pattern, v))) {
      if (is.null(cluster_names) || any(cluster_names %in% names(v))) {
        return(v)
      }
    }
  }
  
  # General scan with bidirectional overlap scoring
  is_hex_vector <- function(x) {
    is.character(x) && !is.null(names(x)) &&
      length(x) >= 2 && all(grepl(hex_pattern, x))
  }
  
  overlap_score <- function(v_names, ref_names) {
    if (is.null(ref_names)) return(1)
    n_common        <- sum(ref_names %in% v_names)
    forward_overlap <- n_common / length(ref_names)
    reverse_overlap <- n_common / length(v_names)
    min(forward_overlap, reverse_overlap)
  }
  
  scored <- list()
  for (key in names(candidates)) {
    if (key %in% c(canonical_misc, canonical_options)) next
    v <- candidates[[key]]
    if (!is_hex_vector(v)) next
    score <- overlap_score(names(v), cluster_names)
    if (score >= min_overlap) {
      scored[[key]] <- list(colors = v, score = score)
    }
  }
  
  if (length(scored) == 0) return(NULL)
  
  best_key <- names(which.max(sapply(scored, function(x) x$score)))
  message(sprintf("[find_stored_colors] Using '%s' (overlap: %.2f)",
                  best_key, scored[[best_key]]$score))
  return(scored[[best_key]]$colors)
}