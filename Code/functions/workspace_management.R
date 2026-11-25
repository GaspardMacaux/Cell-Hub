# workspace_management.R



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
#' Builds a structured data.frame of the parameters actually applied during the session, never the values currently
#' shown in the interface: QC thresholds and pre/post-QC statistics recorded when QC was run (qc_stats, saved with the
#' session), normalization, PCA, neighbour graph, UMAP and clustering read from the Seurat command history, integration,
#' DoubletFinder settings and cluster composition read from the object, and package versions.
#'
#' @param module Character string: "single" or "multiple"
#' @param qc_stats List recorded when QC was run: thresholds, n/medians before and after (single) or per_dataset
#'   (multiple), and for integrated objects the per-dataset normalization settings (NULL if QC not run in this session)
#' @param seurat_obj Seurat object (NULL if not yet loaded)
#' @return data.frame with columns: category, parameter, value, description
collect_analysis_params <- function(module = "single", qc_stats = NULL, seurat_obj = NULL) {
  
  params <- list()
  
  algo_labels <- c(
    "1" = "Original Louvain",
    "2" = "Louvain with Multilevel Refinement",
    "3" = "SLM Algorithm"
  )
  
  # ── QC THRESHOLDS (values applied, recorded when the QC filter was run; never read from the interface) ──
  thr <- qc_stats$thresholds
  thr_value <- function(key) {
    if (is.null(thr)) return("not recorded (QC not applied in this session)")
    v <- thr[[key]]
    if (is.null(v) || length(v) == 0 || all(is.na(v))) return("none")
    if (length(unique(v)) > 1) return(paste0("differs by dataset (", paste(v, collapse = " / "), ")"))
    v[1]
  }
  params[["qc_nfeature_min"]] <- list(
    category    = "QC Thresholds",
    parameter   = "Min nFeature_RNA",
    value       = thr_value("min_features"),
    description = "Cells kept with nFeature_RNA strictly above this value"
  )
  params[["qc_nfeature_max"]] <- list(
    category    = "QC Thresholds",
    parameter   = "Max nFeature_RNA",
    value       = thr_value("max_features"),
    description = "Cells kept with nFeature_RNA strictly below this value"
  )
  params[["qc_mt_max"]] <- list(
    category    = "QC Thresholds",
    parameter   = "Max mitochondrial %",
    value       = thr_value("max_mt_percent"),
    description = "Cells kept with a mitochondrial percentage strictly below this value"
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
  
  # ── NORMALIZATION, DIMENSIONALITY REDUCTION AND CLUSTERING (values applied, from the Seurat command history) ──
  # merge() drops the command history: for integrated objects the per-dataset values recorded at integration
  # (qc_stats) are used when the merged object holds no NormalizeData / SCTransform entry
  applied <- if (!is.null(seurat_obj)) extract_applied_params(seurat_obj) else extract_applied_params(NULL)
  pick <- function(object_value, session_value) {
    if (!identical(as.character(object_value), "N/A")) return(object_value)
    session_value <- setdiff(unique(as.character(session_value)), "N/A")
    if (length(session_value) > 0) return(paste(session_value, collapse = " / "))
    "N/A"
  }
  norm_value <- pick(applied$normalization_method, qc_stats$normalization_method)
  if (identical(norm_value, "N/A") && module == "multiple" && !is.null(seurat_obj)) {
    integ_norm <- extract_integration_info(seurat_obj)$normalization_method
    if (!identical(integ_norm, "unknown")) norm_value <- integ_norm
  }
  params[["norm_method"]] <- list(
    category    = "Normalization",
    parameter   = "Normalization method",
    value       = norm_value,
    description = if (module == "multiple") "Normalization of the datasets (per dataset before integration, or of the merged object)" else "Normalization method applied to raw counts"
  )
  params[["scale_factor"]] <- list(
    category    = "Normalization",
    parameter   = "Scale factor",
    value       = pick(applied$scale_factor, qc_stats$scale_factor),
    description = "Scale factor used in LogNormalize (N/A for SCTransform)"
  )
  params[["num_var_features"]] <- list(
    category    = "Normalization",
    parameter   = "Variable features (n)",
    value       = pick(applied$n_variable_features, qc_stats$n_variable_features),
    description = "Number of highly variable features selected for PCA"
  )
  params[["n_pcs_computed"]] <- list(
    category    = "Dimensionality Reduction",
    parameter   = "PCs computed (RunPCA npcs)",
    value       = applied$npcs,
    description = "Number of principal components computed"
  )
  params[["n_pcs"]] <- list(
    category    = "Dimensionality Reduction",
    parameter   = "Number of PCs (FindNeighbors / UMAP)",
    value       = applied$dims,
    description = paste0("Dimensions of the '", applied$reduction, "' reduction used by the last FindNeighbors call")
  )
  params[["k_param"]] <- list(
    category    = "Dimensionality Reduction",
    parameter   = "k neighbors (FindNeighbors k.param)",
    value       = applied$k_param,
    description = "Number of nearest neighbours of the shared nearest-neighbour graph"
  )
  params[["umap_n_neighbors"]] <- list(
    category    = "Dimensionality Reduction",
    parameter   = "UMAP n.neighbors",
    value       = applied$umap_n_neighbors,
    description = "Local neighbourhood size of the UMAP embedding"
  )
  params[["umap_min_dist"]] <- list(
    category    = "Dimensionality Reduction",
    parameter   = "UMAP min.dist",
    value       = applied$umap_min_dist,
    description = "Minimum distance between points in the UMAP embedding"
  )
  params[["resolution"]] <- list(
    category    = "Clustering",
    parameter   = "Clustering resolution",
    value       = applied$resolution,
    description = "Resolution of the last FindClusters call (higher = more clusters)"
  )
  params[["algorithm"]] <- list(
    category    = "Clustering",
    parameter   = "Clustering algorithm",
    value       = if (as.character(applied$algorithm) %in% names(algo_labels)) algo_labels[[as.character(applied$algorithm)]] else applied$algorithm,
    description = paste0("Graph-based clustering algorithm of the last FindClusters call (graph: ", applied$cluster_graph, ")")
  )
  
  # ── INTEGRATION (multiple only) ───────────────────────────────────────────────
  if (module == "multiple" && !is.null(seurat_obj)) {
    
    integ_info <- extract_integration_info(seurat_obj)
    # The method recorded by the integration step of this session takes precedence over what is read from the object
    if (!is.null(qc_stats$integration_method) && nzchar(as.character(qc_stats$integration_method)[1])) integ_info$integration_method <- as.character(qc_stats$integration_method)[1]
    if (!is.null(qc_stats$harmony_vars) && length(qc_stats$harmony_vars) > 0) integ_info$harmony_vars <- paste(qc_stats$harmony_vars, collapse = ", ")
    
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
      description = "Method used to combine the datasets (Seurat CCA anchors, Harmony, fastMNN or merge without correction)"
    )
    if (!is.null(integ_info$harmony_vars)) {
      params[["harmony_vars"]] <- list(
        category    = "Integration",
        parameter   = "Harmony - Batch correction variable",
        value       = integ_info$harmony_vars,
        description = "Metadata variable used for Harmony batch correction"
      )
    }
    # Settings of the integration step stored in the object (CCA anchor parameters, Harmony settings)
    integ_params <- seurat_obj@misc$integration_params
    for (nm in setdiff(names(integ_params), "method")) {
      params[[paste0("integ_", nm)]] <- list(
        category    = "Integration",
        parameter   = paste("Integration -", nm),
        value       = paste(as.character(unlist(integ_params[[nm]])), collapse = ", "),
        description = "Value applied (stored in the object at integration)"
      )
    }
    
  }
  
  # ── DOUBLETFINDER (single only) ───────────────────────────────────────────────
  # Parameters are read from the object (values actually applied), never from UI inputs,
  # so the export reflects what was run — including the automatically optimized pK
  if (module == "single") {
    df_params <- if (!is.null(seurat_obj)) seurat_obj@misc$doubletfinder_params else NULL
    if (is.null(df_params)) {
      params[["doublet_status"]] <- list(
        category    = "DoubletFinder",
        parameter   = "Status",
        value       = "Not run",
        description = "DoubletFinder was not run on this object"
      )
    } else {
      params[["doublet_rate"]] <- list(
        category    = "DoubletFinder",
        parameter   = "Expected doublet rate (%)",
        value       = df_params$doublet_rate,
        description = "Expected doublet rate used to compute nExp (~1% per 1000 cells loaded)"
      )
      params[["doublet_pcs"]] <- list(
        category    = "DoubletFinder",
        parameter   = "PCs used",
        value       = df_params$n_pcs,
        description = "Number of principal components used by paramSweep and doubletFinder"
      )
      params[["doublet_pN"]] <- list(
        category    = "DoubletFinder",
        parameter   = "pN",
        value       = df_params$pN,
        description = "Proportion of artificial doublets generated (pN parameter)"
      )
      params[["doublet_pK"]] <- list(
        category    = "DoubletFinder",
        parameter   = "pK",
        value       = df_params$pK,
        description = paste0("PC neighborhood size applied (source: ", df_params$pK_source, ")")
      )
      params[["doublet_nExp"]] <- list(
        category    = "DoubletFinder",
        parameter   = "nExp",
        value       = df_params$nExp,
        description = "Number of doublets expected (doublet rate x cells at detection time)"
      )
      params[["doublet_detected"]] <- list(
        category    = "DoubletFinder",
        parameter   = "Doublets detected",
        value       = df_params$n_doublets,
        description = paste0("Cells classified as doublets out of ", df_params$n_cells, " at detection time")
      )
      params[["doublet_version"]] <- list(
        category    = "DoubletFinder",
        parameter   = "DoubletFinder version",
        value       = df_params$package_version,
        description = "Installed DoubletFinder package version (seed 42 for paramSweep and doubletFinder)"
      )
    }
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
    value       = if (!identical(as.character(applied$cluster_seed), "N/A")) applied$cluster_seed else 42,
    description = "Fixed seed used in RunPCA, RunUMAP, FindNeighbors, FindClusters (value read from the last FindClusters call)"
  )
  params[["r_version"]] <- list(
    category    = "Reproducibility",
    parameter   = "R version",
    value       = paste(R.version$major, R.version$minor, sep = "."),
    description = "R version of the session that produced the export"
  )
  params[["seurat_version"]] <- list(
    category    = "Reproducibility",
    parameter   = "Seurat / SeuratObject versions",
    value       = paste(as.character(utils::packageVersion("Seurat")), as.character(utils::packageVersion("SeuratObject")), sep = " / "),
    description = "Package versions used for the analysis"
  )
  # Versions of the other analysis packages installed in the session (absent packages are left out)
  key_packages <- c("harmony", "batchelor", "DoubletFinder", "SingleR", "DESeq2", "fgsea", "CellChat", "monocle3", "presto")
  installed <- key_packages[vapply(key_packages, requireNamespace, logical(1), quietly = TRUE)]
  params[["package_versions"]] <- list(
    category    = "Reproducibility",
    parameter   = "Other package versions",
    value       = paste(vapply(installed, function(p) paste(p, as.character(utils::packageVersion(p))), character(1)), collapse = "; "),
    description = "Versions of the analysis packages installed in the session"
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

# Parameters actually applied, read from the Seurat command history (Seurat records every call of NormalizeData,
# SCTransform, FindVariableFeatures, RunPCA, FindNeighbors, FindClusters and RunUMAP with its parameters and time stamp).
# The most recent call of each kind is used; "N/A" when the step was not run on this object. merge() drops the history,
# so integrated objects only hold the calls made after the merge.
extract_applied_params <- function(seurat_obj) {
  out <- list(normalization_method = "N/A", scale_factor = "N/A", n_variable_features = "N/A", npcs = "N/A",
              k_param = "N/A", dims = "N/A", reduction = "N/A", umap_n_neighbors = "N/A", umap_min_dist = "N/A",
              resolution = "N/A", algorithm = "N/A", cluster_graph = "N/A", cluster_seed = "N/A")
  cmds <- if (is.null(seurat_obj)) list() else tryCatch(seurat_obj@commands, error = function(e) list())
  if (length(cmds) == 0) return(out)
  # Most recent command whose name starts with prefix (re-running a step replaces its entry in place, so the list order is not chronological)
  last_cmd <- function(prefix) {
    hit <- grep(paste0("^", prefix, "(\\.|$)"), names(cmds), value = TRUE)
    if (length(hit) == 0) return(NULL)
    stamps <- vapply(hit, function(h) as.numeric(tryCatch(cmds[[h]]@time.stamp, error = function(e) NA)), numeric(1))
    cmds[[hit[if (all(is.na(stamps))) length(hit) else which.max(stamps)]]]
  }
  stamp <- function(cmd) if (is.null(cmd)) -Inf else as.numeric(cmd@time.stamp)
  val <- function(v) if (is.null(v) || length(v) == 0) "N/A" else paste(as.character(v), collapse = ", ")
  # Ranges of dimensions written compactly (1:30 -> "1-30")
  dims_text <- function(d) {
    if (is.null(d) || length(d) == 0) return("N/A")
    d <- as.integer(d)
    if (length(d) > 1 && all(diff(d) == 1)) paste0(min(d), "-", max(d)) else paste(d, collapse = ", ")
  }
  sct <- last_cmd("SCTransform"); nd <- last_cmd("NormalizeData"); fv <- last_cmd("FindVariableFeatures")
  if (!is.null(sct) && stamp(sct) >= stamp(nd)) {
    out$normalization_method <- "SCTransform"
    out$n_variable_features <- val(sct@params$variable.features.n)
  } else if (!is.null(nd)) {
    out$normalization_method <- val(nd@params$normalization.method)
    out$scale_factor <- val(nd@params$scale.factor)
    if (!is.null(fv)) out$n_variable_features <- val(fv@params$nfeatures)
  }
  pc <- last_cmd("RunPCA")
  if (!is.null(pc)) out$npcs <- val(pc@params$npcs)
  fn <- last_cmd("FindNeighbors")
  if (!is.null(fn)) {
    out$k_param <- val(fn@params$k.param)
    out$dims <- dims_text(fn@params$dims)
    out$reduction <- val(fn@params$reduction)
  }
  um <- last_cmd("RunUMAP")
  if (!is.null(um)) {
    out$umap_n_neighbors <- val(um@params$n.neighbors)
    out$umap_min_dist <- val(um@params$min.dist)
  }
  fc <- last_cmd("FindClusters")
  if (!is.null(fc)) {
    out$resolution <- val(fc@params$resolution)
    out$algorithm <- val(fc@params$algorithm)
    out$cluster_graph <- val(fc@params$graph.name)
    out$cluster_seed <- val(fc@params$random.seed)
  }
  return(out)
}
#' Extract integration and normalization info from a Seurat object
#'
#' Reads metadata columns, reductions and the Seurat command history to recover the integration method, the
#' normalization method and the dataset names without relying on reactive tracking. Objects reloaded from a file keep
#' their metadata and reductions but may have lost part of the command history (merge() does not keep it), so several
#' sources are checked in turn.
#'
#' @param seurat_obj Seurat object after integration
#' @return List with: integration_method, normalization_method, harmony_vars, dataset_names
extract_integration_info <- function(seurat_obj) {
  result <- list(integration_method = "unknown", normalization_method = "unknown", harmony_vars = NULL, dataset_names = character(0))
  if (is.null(seurat_obj)) return(result)
  meta <- seurat_obj@meta.data
  cmd_names <- tryCatch(names(seurat_obj@commands), error = function(e) character(0))
  reductions <- names(seurat_obj@reductions)
  # ── Integration method: metadata column, then command history, then integrated reductions ──
  if ("integration_method" %in% colnames(meta)) {
    result$integration_method <- as.character(unique(meta$integration_method)[1])
  } else if ("merge_method" %in% colnames(meta)) {
    result$integration_method <- as.character(unique(meta$merge_method)[1])
  } else if (any(grepl("^RunHarmony", cmd_names)) || "harmony" %in% reductions) {
    result$integration_method <- "harmony"
  } else if (any(grepl("^IntegrateLayers", cmd_names))) {
    cmd <- seurat_obj@commands[[tail(grep("^IntegrateLayers", cmd_names, value = TRUE), 1)]]
    method <- tryCatch(cmd@params$method, error = function(e) NULL)
    result$integration_method <- if (is.function(method) || is.null(method)) "IntegrateLayers" else as.character(method)[1]
  } else if (any(grepl("^FindIntegrationAnchors|^IntegrateData", cmd_names)) || "integrated" %in% names(seurat_obj@assays)) {
    result$integration_method <- "Seurat CCA anchors"
  } else if (any(grepl("mnn", reductions, ignore.case = TRUE))) {
    result$integration_method <- "fastMNN"
  } else if (any(grepl("cca", reductions, ignore.case = TRUE))) {
    result$integration_method <- "Seurat CCA anchors"
  }
  # ── Dataset names from the dataset column ──
  if ("dataset" %in% colnames(meta)) {
    result$dataset_names <- unique(as.character(meta$dataset))
  } else if ("orig.ident" %in% colnames(meta)) {
    result$dataset_names <- unique(as.character(meta$orig.ident))
  }
  # ── Normalization: SCTransform assay or command, NormalizeData command, otherwise log-normalized data layer ──
  if ("SCT" %in% names(seurat_obj@assays) || any(grepl("^SCTransform", cmd_names))) {
    result$normalization_method <- "SCTransform"
  } else if (any(grepl("^NormalizeData", cmd_names))) {
    cmd <- seurat_obj@commands[[tail(grep("^NormalizeData", cmd_names, value = TRUE), 1)]]
    norm_method <- tryCatch(cmd@params$normalization.method, error = function(e) NULL)
    result$normalization_method <- if (!is.null(norm_method)) as.character(norm_method) else "LogNormalize"
  } else if (length(tryCatch(SeuratObject::Layers(seurat_obj[["RNA"]], search = "data"), error = function(e) character(0))) > 0) {
    result$normalization_method <- "LogNormalize (data layer present; method not recorded)"
  }
  # ── Harmony variables from the RunHarmony command ──
  harmony_cmd <- grep("^RunHarmony", cmd_names, value = TRUE)
  if (length(harmony_cmd) > 0) {
    harmony_vars <- tryCatch(seurat_obj@commands[[tail(harmony_cmd, 1)]]@params$group.by.vars, error = function(e) NULL)
    if (!is.null(harmony_vars)) result$harmony_vars <- paste(harmony_vars, collapse = ", ")
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


# Session state stored inside the Seurat object, so one .rds file restores the whole analysis.
# Only lightweight tabular results are kept: plots are never stored (ggplot objects embed full data copies)
pack_session_state <- function(seurat_obj, gene_tables = list(), qc_stats = NULL, module = "single") {
  seurat_obj@misc$cellhub_session <- list(
    format_version = 1L,
    module         = module,
    saved_at       = Sys.time(),
    seurat_version = as.character(utils::packageVersion("Seurat")),
    qc_stats       = qc_stats,
    gene_tables    = gene_tables
  )
  return(seurat_obj)
}

# Read back the session state written by pack_session_state(); found = FALSE for plain Seurat objects
unpack_session_state <- function(seurat_obj) {
  state <- tryCatch(seurat_obj@misc$cellhub_session, error = function(e) NULL)
  if (is.null(state)) return(list(found = FALSE, gene_tables = list(), qc_stats = NULL, saved_at = NULL, module = NULL))
  gene_tables <- if (is.list(state$gene_tables)) state$gene_tables else list()
  return(list(found = TRUE, gene_tables = gene_tables, qc_stats = state$qc_stats,
              saved_at = state$saved_at, module = state$module))
}

# Human-readable summary of what a Seurat object contains (Item / Value table)
summarize_seurat_object <- function(seurat_obj) {
  if (is.null(seurat_obj)) return(data.frame(Item = character(0), Value = character(0)))
  meta <- seurat_obj@meta.data
  cluster_col <- if ("ClusterIdents" %in% colnames(meta)) "ClusterIdents" else if ("seurat_clusters" %in% colnames(meta)) "seurat_clusters" else NULL
  n_clusters <- if (!is.null(cluster_col)) length(unique(meta[[cluster_col]])) else length(levels(Idents(seurat_obj)))
  # Layers() only exists for Seurat v5 Assay5 objects
  rna_layers <- tryCatch(paste(Layers(seurat_obj[["RNA"]]), collapse = ", "), error = function(e) "N/A")
  df_params <- seurat_obj@misc$doubletfinder_params
  items <- c(
    "Cells / nuclei"         = format(ncol(seurat_obj), big.mark = ","),
    "Genes (features)"       = format(nrow(seurat_obj), big.mark = ","),
    "Assays"                 = paste(names(seurat_obj@assays), collapse = ", "),
    "Default assay"          = DefaultAssay(seurat_obj),
    "RNA layers"             = rna_layers,
    "Dimensional reductions" = if (length(seurat_obj@reductions) > 0) paste(names(seurat_obj@reductions), collapse = ", ") else "none",
    "Cluster column"         = if (!is.null(cluster_col)) cluster_col else "Idents()",
    "Clusters"               = as.character(n_clusters),
    "Metadata columns"       = as.character(ncol(meta)),
    "Custom cluster colors"  = if (!is.null(seurat_obj@misc$cluster_colors)) "stored" else "none",
    "DoubletFinder"          = if (!is.null(df_params)) sprintf("run (pK = %s, %s doublets)", df_params$pK, df_params$n_doublets) else "not run",
    "Object size in memory"  = format(object.size(seurat_obj), units = "auto")
  )
  return(data.frame(Item = names(items), Value = unname(items), stringsAsFactors = FALSE))
}

# Human-readable labels of the stored-table types (Single, Multiple and enrichment modules)
stored_table_type_labels <- function() {
  return(c(single_cluster    = "Cluster vs all",
           cluster_group     = "Differential expression",
           metadata_group    = "Condition comparison",
           exclusive_markers = "Exclusive biomarkers",
           enrichment        = "Pathway enrichment"))
}

# Stored tables usable as differential expression input (Venn diagrams, volcano plots): enrichment results are left out
de_tables_only <- function(gene_tables) {
  if (length(gene_tables) == 0) return(list())
  return(gene_tables[vapply(gene_tables, function(e) !identical(e$type, "enrichment"), logical(1))])
}

# Stored analyses for the parameter export (CSV and PDF report), same columns as collect_analysis_params():
# - category "Stored Analyses": one row per stored table (what was compared, result size, settings);
# - category "Pathway Enrichment - <description>": every parameter of each stored enrichment run (gene sets and their
#   version, universe, thresholds, seed), so the report documents the annotation release used.
collect_stored_analyses <- function(gene_tables) {
  empty <- data.frame(category = character(0), parameter = character(0), value = character(0), description = character(0), stringsAsFactors = FALSE)
  if (length(gene_tables) == 0) return(empty)
  labels <- stored_table_type_labels()
  # Settings as "name = value" pairs (vectors collapsed with commas)
  flatten <- function(p) {
    if (length(p) == 0) return("")
    paste(vapply(names(p), function(k) paste0(k, " = ", paste(as.character(unlist(p[[k]])), collapse = ", ")), character(1)), collapse = "; ")
  }
  rows <- list()
  for (id in names(gene_tables)) {
    entry <- gene_tables[[id]]
    type <- entry$type %||% "unknown"
    label <- if (type %in% names(labels)) labels[[type]] else type
    tbl <- entry$data
    desc <- entry$description %||% id
    if (identical(type, "enrichment")) {
      n_sig <- entry$n_significant %||% c(ORA = NA, GSEA = NA)
      value <- sprintf("%s: ORA %s, GSEA %s significant sets (adjusted p <= %s)", label, n_sig[["ORA"]], n_sig[["GSEA"]], entry$max_p_adj %||% 0.05)
      rows[[length(rows) + 1]] <- data.frame(category = "Stored Analyses", parameter = desc, value = value,
                                             description = "Parameters in the Pathway Enrichment section", stringsAsFactors = FALSE)
      p <- entry$parameters
      if (length(p) > 0) {
        rows[[length(rows) + 1]] <- data.frame(category = paste0("Pathway Enrichment - ", desc), parameter = names(p),
                                               value = vapply(p, function(v) paste(as.character(unlist(v)), collapse = ", "), character(1)),
                                               description = "", stringsAsFactors = FALSE)
      }
    } else {
      n_sig <- if (!is.null(tbl) && "p_val_adj" %in% colnames(tbl)) sum(suppressWarnings(as.numeric(tbl$p_val_adj)) < 0.05, na.rm = TRUE) else NA
      value <- sprintf("%s: %d genes (%s with adjusted p < 0.05)", label, if (!is.null(tbl)) nrow(tbl) else 0L, n_sig)
      rows[[length(rows) + 1]] <- data.frame(category = "Stored Analyses", parameter = desc, value = value,
                                             description = flatten(entry$parameters), stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  return(out)
}

# One row per stored result table, labelled with its human-readable description (ID kept for deletion)
summarize_gene_tables <- function(gene_tables) {
  empty <- data.frame(ID = character(0), Table = character(0), Type = character(0), Genes = integer(0),
                      `Significant (adj. p < 0.05)` = integer(0), Created = character(0), Size = character(0),
                      check.names = FALSE, stringsAsFactors = FALSE)
  if (length(gene_tables) == 0) return(empty)
  type_labels <- stored_table_type_labels()
  rows <- lapply(names(gene_tables), function(id) {
    entry <- gene_tables[[id]]
    tbl   <- entry$data
    type  <- if (!is.null(entry$type)) entry$type else "unknown"
    # Enrichment tables hold gene sets (p.adjust column), DE tables hold genes (p_val_adj column)
    n_sig <- if (identical(type, "enrichment") && !is.null(tbl) && "p.adjust" %in% colnames(tbl)) {
      sum(tbl$p.adjust < 0.05, na.rm = TRUE)
    } else if (!is.null(tbl) && "p_val_adj" %in% colnames(tbl)) {
      sum(suppressWarnings(as.numeric(tbl$p_val_adj)) < 0.05, na.rm = TRUE)
    } else NA_integer_
    data.frame(
      ID      = id,
      Table   = if (!is.null(entry$description)) entry$description else id,
      Type    = if (type %in% names(type_labels)) type_labels[[type]] else type,
      Genes   = if (!is.null(tbl)) nrow(tbl) else 0L,
      `Significant (adj. p < 0.05)` = n_sig,
      Created = if (!is.null(entry$timestamp)) format(entry$timestamp, "%Y-%m-%d %H:%M") else "",
      Size    = format(object.size(entry), units = "auto"),
      check.names = FALSE, stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$Created, decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  return(out)
}


get_cluster_colors_merge <- function(seurat_object) {
  if (is.null(seurat_object)) return(NULL)
  colors <- find_stored_colors(seurat_object, cluster_names = levels(Idents(seurat_object)))
  if (!is.null(colors)) {
    message("Using stored cluster colors for multiple datasets")
    return(colors)
  }
  message("Using Seurat default colors for multiple datasets")
  return(NULL)
}