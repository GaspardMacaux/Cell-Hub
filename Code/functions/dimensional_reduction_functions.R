# dimensional_reduction_functions_newarch.R

############################## UMAP Reproducibility Functions ##############################

runUMAP_reproducible <- function(object, 
                                 dims, 
                                 seed = 42, 
                                 n.neighbors = 30L, 
                                 min.dist = 0.3, 
                                 metric = "cosine", 
                                 verbose = FALSE, 
                                 ...) {
  # Run UMAP with reproducible seed settings
  # Args:
  #   object: Seurat object with PCA reduction computed
  #   dims: Integer vector of PCA dimensions to use (e.g., 1:30)
  #   seed: Integer seed for reproducibility (default: 42)
  #   n.neighbors: Number of neighbors for UMAP (default: 30L)
  #   min.dist: Minimum distance parameter for UMAP (default: 0.3)
  #   metric: Distance metric to use (default: "cosine")
  #   verbose: Whether to print progress messages (default: FALSE)
  #   ...: Additional parameters passed to Seurat::RunUMAP
  # Returns:
  #   Seurat object with UMAP reduction added and parameters stored in @misc
  
  tryCatch({
    # Input validation
    if (!inherits(object, "Seurat")) {
      stop("Input must be a Seurat object")
    }
    
    if (!"pca" %in% names(object@reductions)) {
      stop("PCA reduction not found. Run PCA first using RunPCA()")
    }
    
    if (!is.numeric(dims) || length(dims) == 0) {
      stop("'dims' must be a non-empty numeric vector")
    }
    
    # Check available PCs and adjust dims if necessary
    available_pcs <- ncol(object@reductions$pca@cell.embeddings)
    if (max(dims) > available_pcs) {
      warning(paste0(
        "Requested dims (max=", max(dims), ") exceeds available PCs (", 
        available_pcs, "). Adjusting to use available dimensions."
      ))
      dims <- dims[dims <= available_pcs]
    }
    
    if (!is.numeric(seed) || length(seed) != 1) {
      stop("'seed' must be a single integer value")
    }
    
    # Set R's global seed for reproducibility
    set.seed(seed)
    
    if (verbose) {
      message(paste0(
        "Running UMAP with: ",
        "dims=1:", max(dims), 
        " | seed=", seed,
        " | n.neighbors=", n.neighbors,
        " | min.dist=", min.dist,
        " | metric=", metric
      ))
    }
    
    # Run UMAP with explicit parameters for reproducibility
    object <- RunUMAP(
      object = object,
      dims = dims,
      seed.use = seed,
      n.neighbors = n.neighbors,
      min.dist = min.dist,
      metric = metric,
      verbose = verbose,
      ...
    )
    
    # Store metadata about UMAP computation for traceability
    object@misc$umap_params <- list(
      dims = dims,
      seed = seed,
      n.neighbors = n.neighbors,
      min.dist = min.dist,
      metric = metric,
      computation_date = Sys.time()
    )
    
    if (verbose) {
      message("UMAP computation completed successfully")
    }
    
    return(object)
    
  }, error = function(e) {
    message(paste("ERROR in runUMAP_reproducible:", e$message))
    stop(e)
  })
}


hasValidUMAP <- function(object, min_cells = 10) {
  # Check if UMAP reduction exists and is valid
  # Args:
  #   object: Seurat object to check
  #   min_cells: Minimum number of cells expected (default: 10)
  # Returns:
  #   Logical TRUE if valid UMAP exists, FALSE otherwise
  
  tryCatch({
    if (!inherits(object, "Seurat")) {
      warning("Input is not a Seurat object")
      return(FALSE)
    }
    
    # Check if UMAP reduction exists
    if (!"umap" %in% names(object@reductions)) {
      return(FALSE)
    }
    
    # Check UMAP dimensions
    umap_coords <- object@reductions$umap@cell.embeddings
    
    if (is.null(umap_coords) || nrow(umap_coords) < min_cells) {
      warning(paste0(
        "UMAP exists but has insufficient cells: ",
        nrow(umap_coords), " < ", min_cells
      ))
      return(FALSE)
    }
    
    # Check for proper dimensions (should be 2D)
    if (ncol(umap_coords) != 2) {
      warning(paste0(
        "UMAP has unexpected dimensions: ",
        ncol(umap_coords), " (expected 2)"
      ))
      return(FALSE)
    }
    
    # Check for NA values
    if (any(is.na(umap_coords))) {
      warning("UMAP contains NA values")
      return(FALSE)
    }
    
    return(TRUE)
    
  }, error = function(e) {
    message(paste("ERROR in hasValidUMAP:", e$message))
    return(FALSE)
  })
}


ensureUMAP <- function(object, 
                       dims = 1:30, 
                       seed = 42, 
                       force = FALSE, 
                       ...) {
  # Ensure UMAP exists and is valid, computing if necessary
  # Args:
  #   object: Seurat object
  #   dims: Dimensions to use for UMAP (default: 1:30)
  #   seed: Seed for reproducibility (default: 42)
  #   force: Logical, if TRUE recomputes UMAP even if it exists (default: FALSE)
  #   ...: Additional parameters passed to runUMAP_reproducible
  # Returns:
  #   Seurat object with valid UMAP reduction
  
  tryCatch({
    if (!inherits(object, "Seurat")) {
      stop("Input must be a Seurat object")
    }
    
    # Check if recomputation is needed
    needs_computation <- force || !hasValidUMAP(object)
    
    if (needs_computation) {
      message("Computing UMAP...")
      object <- runUMAP_reproducible(
        object = object,
        dims = dims,
        seed = seed,
        verbose = TRUE,
        ...
      )
    } else {
      message("Valid UMAP already exists. Use force=TRUE to recompute.")
    }
    
    return(object)
    
  }, error = function(e) {
    message(paste("ERROR in ensureUMAP:", e$message))
    stop(e)
  })
}


getUMAPParams <- function(object) {
  # Extract UMAP parameters from Seurat object
  # Args:
  #   object: Seurat object
  # Returns:
  #   List of UMAP parameters or NULL if not found
  
  tryCatch({
    if (!inherits(object, "Seurat")) {
      stop("Input must be a Seurat object")
    }
    
    if (!"umap_params" %in% names(object@misc)) {
      warning(paste(
        "UMAP parameters not found in object@misc.",
        "Object may not have been created with runUMAP_reproducible()"
      ))
      return(NULL)
    }
    
    return(object@misc$umap_params)
    
  }, error = function(e) {
    message(paste("ERROR in getUMAPParams:", e$message))
    return(NULL)
  })
}


compareUMAPParams <- function(object1, object2) {
  # Compare UMAP parameters between two Seurat objects
  # Args:
  #   object1: First Seurat object
  #   object2: Second Seurat object
  # Returns:
  #   Logical TRUE if parameters match, FALSE otherwise
  
  tryCatch({
    params1 <- getUMAPParams(object1)
    params2 <- getUMAPParams(object2)
    
    if (is.null(params1) || is.null(params2)) {
      warning("Cannot compare: UMAP parameters not found in one or both objects")
      return(FALSE)
    }
    
    # Remove computation date for comparison
    params1$computation_date <- NULL
    params2$computation_date <- NULL
    
    # Compare parameters
    differences <- list()
    all_params <- unique(c(names(params1), names(params2)))
    
    for (param in all_params) {
      val1 <- params1[[param]]
      val2 <- params2[[param]]
      
      if (!identical(val1, val2)) {
        differences[[param]] <- list(object1 = val1, object2 = val2)
      }
    }
    
    if (length(differences) > 0) {
      message("UMAP parameters differ:")
      print(differences)
      return(FALSE)
    } else {
      message("UMAP parameters match between objects")
      return(TRUE)
    }
    
  }, error = function(e) {
    message(paste("ERROR in compareUMAPParams:", e$message))
    return(FALSE)
  })
}


############################## PCA Reproducibility Functions ##############################

runPCA_reproducible <- function(object,
                                npcs = 50,
                                seed = 42,
                                features = NULL,
                                verbose = FALSE,
                                ...) {
  # Run PCA with reproducible seed settings
  # Args:
  #   object: Seurat object with normalized data
  #   npcs: Number of principal components to compute (default: 50)
  #   seed: Integer seed for reproducibility (default: 42)
  #   features: Features to use for PCA (default: NULL uses variable features)
  #   verbose: Whether to print progress messages (default: FALSE)
  #   ...: Additional parameters passed to Seurat::RunPCA
  # Returns:
  #   Seurat object with PCA reduction added and parameters stored in @misc
  
  tryCatch({
    # Input validation
    if (!inherits(object, "Seurat")) {
      stop("Input must be a Seurat object")
    }
    
    if (!is.numeric(npcs) || length(npcs) != 1 || npcs < 1) {
      stop("'npcs' must be a positive integer")
    }
    
    if (!is.numeric(seed) || length(seed) != 1) {
      stop("'seed' must be a single integer value")
    }
    
    # Set R's global seed for reproducibility
    set.seed(seed)
    
    if (verbose) {
      message(paste0(
        "Running PCA with: ",
        "npcs=", npcs,
        " | seed=", seed
      ))
    }
    
    # Run PCA with explicit parameters
    object <- RunPCA(
      object = object,
      npcs = npcs,
      features = features,
      seed.use = seed,
      verbose = verbose,
      ...
    )
    
    # Store metadata about PCA computation for traceability
    object@misc$pca_params <- list(
      npcs = npcs,
      seed = seed,
      features = if (!is.null(features)) length(features) else "variable_features",
      computation_date = Sys.time()
    )
    
    if (verbose) {
      message("PCA computation completed successfully")
    }
    
    return(object)
    
  }, error = function(e) {
    message(paste("ERROR in runPCA_reproducible:", e$message))
    stop(e)
  })
}


hasValidPCA <- function(object, min_pcs = 10) {
  # Check if PCA reduction exists and is valid
  # Args:
  #   object: Seurat object to check
  #   min_pcs: Minimum number of PCs expected (default: 10)
  # Returns:
  #   Logical TRUE if valid PCA exists, FALSE otherwise
  
  tryCatch({
    if (!inherits(object, "Seurat")) {
      warning("Input is not a Seurat object")
      return(FALSE)
    }
    
    # Check if PCA reduction exists
    if (!"pca" %in% names(object@reductions)) {
      return(FALSE)
    }
    
    # Check PCA dimensions
    pca_embeddings <- object@reductions$pca@cell.embeddings
    
    if (is.null(pca_embeddings)) {
      return(FALSE)
    }
    
    n_pcs <- ncol(pca_embeddings)
    if (n_pcs < min_pcs) {
      warning(paste0(
        "PCA exists but has insufficient PCs: ",
        n_pcs, " < ", min_pcs
      ))
      return(FALSE)
    }
    
    # Check for NA values
    if (any(is.na(pca_embeddings))) {
      warning("PCA contains NA values")
      return(FALSE)
    }
    
    return(TRUE)
    
  }, error = function(e) {
    message(paste("ERROR in hasValidPCA:", e$message))
    return(FALSE)
  })
}


############################## Neighbors Computation Functions ##############################

findNeighbors_reproducible <- function(object,
                                       dims,
                                       k.param = 20,
                                       seed = 42,
                                       verbose = FALSE,
                                       ...) {
  # Find neighbors with reproducible seed settings
  # Args:
  #   object: Seurat object with PCA reduction
  #   dims: Integer vector of PCA dimensions to use (e.g., 1:30)
  #   k.param: Number of neighbors to find (default: 20)
  #   seed: Integer seed for reproducibility (default: 42)
  #   verbose: Whether to print progress messages (default: FALSE)
  #   ...: Additional parameters passed to Seurat::FindNeighbors
  # Returns:
  #   Seurat object with neighbor graphs added
  
  tryCatch({
    # Input validation
    if (!inherits(object, "Seurat")) {
      stop("Input must be a Seurat object")
    }
    
    if (!"pca" %in% names(object@reductions)) {
      stop("PCA reduction not found. Run PCA first using RunPCA()")
    }
    
    if (!is.numeric(dims) || length(dims) == 0) {
      stop("'dims' must be a non-empty numeric vector")
    }
    
    # Check available PCs and adjust dims if necessary
    available_pcs <- ncol(object@reductions$pca@cell.embeddings)
    if (max(dims) > available_pcs) {
      warning(paste0(
        "Requested dims (max=", max(dims), ") exceeds available PCs (", 
        available_pcs, "). Adjusting to use available dimensions."
      ))
      dims <- dims[dims <= available_pcs]
    }
    
    # Set seed for reproducibility
    set.seed(seed)
    
    if (verbose) {
      message(paste0(
        "Finding neighbors with: ",
        "dims=1:", max(dims),
        " | k.param=", k.param,
        " | seed=", seed
      ))
    }
    
    # Find neighbors
    object <- FindNeighbors(
      object = object,
      dims = dims,
      k.param = k.param,
      verbose = verbose,
      ...
    )
    
    # Store metadata
    object@misc$neighbors_params <- list(
      dims = dims,
      k.param = k.param,
      seed = seed,
      computation_date = Sys.time()
    )
    
    if (verbose) {
      message("Neighbor computation completed successfully")
    }
    
    return(object)
    
  }, error = function(e) {
    message(paste("ERROR in findNeighbors_reproducible:", e$message))
    stop(e)
  })
}


############################## Clustering Reproducibility Functions ##############################

findClusters_reproducible <- function(object,
                                      resolution = 0.5,
                                      algorithm = 1,
                                      seed = 42,
                                      verbose = FALSE,
                                      ...) {
  # Find clusters with reproducible seed settings
  # Args:
  #   object: Seurat object with neighbor graph computed
  #   resolution: Clustering resolution parameter (default: 0.5)
  #   algorithm: Clustering algorithm (1=Louvain, 2=Louvain multilevel, 3=SLM, 4=Leiden)
  #   seed: Integer seed for reproducibility (default: 42)
  #   verbose: Whether to print progress messages (default: FALSE)
  #   ...: Additional parameters passed to Seurat::FindClusters
  # Returns:
  #   Seurat object with clusters added to metadata
  
  tryCatch({
    if (!inherits(object, "Seurat")) {
      stop("Input must be a Seurat object")
    }
    
    # Accept any assay-prefixed nn/snn graph (RNA_nn, SCT_snn, etc.)
    has_graph <- any(grepl("_nn$|_snn$", names(object@graphs)))
    if (!has_graph) {
      stop("No neighbor graph found. Run FindNeighbors first.")
    }
    
    set.seed(seed)
    
    if (verbose) {
      message(paste0(
        "Finding clusters with: ",
        "resolution=", resolution,
        " | algorithm=", algorithm,
        " | seed=", seed
      ))
    }
    
    object <- FindClusters(
      object      = object,
      resolution  = resolution,
      algorithm   = algorithm,
      random.seed = seed,
      verbose     = verbose,
      ...
    )
    
    object@misc$clustering_params <- list(
      resolution       = resolution,
      algorithm        = algorithm,
      seed             = seed,
      n_clusters       = length(unique(Idents(object))),
      computation_date = Sys.time()
    )
    
    if (verbose) {
      message(paste0(
        "Clustering completed successfully. Found ",
        length(unique(Idents(object))), " clusters."
      ))
    }
    
    return(object)
    
  }, error = function(e) {
    message(paste("ERROR in findClusters_reproducible:", e$message))
    stop(e)
  })
}

############################## Complete Pipeline Functions ##############################

runDimensionalReductionPipeline <- function(object,
                                            npcs = 50,
                                            umap_dims = 1:30,
                                            neighbors_dims = 1:30,
                                            resolution = 0.5,
                                            seed = 42,
                                            verbose = TRUE) {
  # Run complete dimensional reduction pipeline with reproducible settings
  # Args:
  #   object: Seurat object with normalized data
  #   npcs: Number of principal components to compute (default: 50)
  #   umap_dims: Dimensions to use for UMAP (default: 1:30)
  #   neighbors_dims: Dimensions to use for neighbor finding (default: 1:30)
  #   resolution: Clustering resolution (default: 0.5)
  #   seed: Integer seed for reproducibility (default: 42)
  #   verbose: Whether to print progress messages (default: TRUE)
  # Returns:
  #   Seurat object with PCA, UMAP, neighbors, and clusters computed
  
  tryCatch({
    if (verbose) {
      message("=== Starting Dimensional Reduction Pipeline ===")
    }
    
    # Step 1: PCA
    if (verbose) message("Step 1/5: Computing PCA...")
    object <- runPCA_reproducible(
      object = object,
      npcs = npcs,
      seed = seed,
      verbose = verbose
    )
    
    # Step 2: Find Neighbors
    if (verbose) message("Step 2/5: Finding neighbors...")
    object <- findNeighbors_reproducible(
      object = object,
      dims = neighbors_dims,
      seed = seed,
      verbose = verbose
    )
    
    # Step 3: UMAP
    if (verbose) message("Step 3/5: Computing UMAP...")
    object <- runUMAP_reproducible(
      object = object,
      dims = umap_dims,
      seed = seed,
      verbose = verbose
    )
    
    # Step 4: Clustering
    if (verbose) message("Step 4/5: Finding clusters...")
    object <- findClusters_reproducible(
      object = object,
      resolution = resolution,
      seed = seed,
      verbose = verbose
    )
    
    # Step 5: Store pipeline metadata
    if (verbose) message("Step 5/5: Storing pipeline metadata...")
    object@misc$pipeline_params <- list(
      npcs = npcs,
      umap_dims = umap_dims,
      neighbors_dims = neighbors_dims,
      resolution = resolution,
      seed = seed,
      computation_date = Sys.time()
    )
    
    if (verbose) {
      message("=== Dimensional Reduction Pipeline Completed Successfully ===")
      message(paste("Total cells:", ncol(object)))
      message(paste("PCs computed:", npcs))
      message(paste("Clusters found:", length(unique(Idents(object)))))
    }
    
    return(object)
    
  }, error = function(e) {
    message(paste("ERROR in runDimensionalReductionPipeline:", e$message))
    stop(e)
  })
}



runUMAP3D_reproducible <- function(object, 
                                   dims, 
                                   seed = 42, 
                                   n.neighbors = 30L, 
                                   min.dist = 0.3, 
                                   metric = "cosine",
                                   verbose = FALSE, 
                                   ...) {
  # Run 3D UMAP with reproducible seed settings
  # Args:
  #   object: Seurat object with PCA reduction computed
  #   dims: Integer vector of PCA dimensions to use (e.g., 1:30)
  #   seed: Integer seed for reproducibility (default: 42)
  #   n.neighbors: Number of neighbors for UMAP (default: 30L)
  #   min.dist: Minimum distance parameter for UMAP (default: 0.3)
  #   metric: Distance metric to use (default: "cosine")
  #   verbose: Whether to print progress messages (default: FALSE)
  #   ...: Additional parameters passed to Seurat::RunUMAP
  # Returns:
  #   Seurat object with 3D UMAP reduction added in "umap3d" slot
  
  tryCatch({
    # Input validation
    if (!inherits(object, "Seurat")) {
      stop("Input must be a Seurat object")
    }
    
    if (!"pca" %in% names(object@reductions)) {
      stop("PCA reduction not found. Run PCA first using RunPCA()")
    }
    
    if (!is.numeric(dims) || length(dims) == 0) {
      stop("'dims' must be a non-empty numeric vector")
    }
    
    # Check available PCs and adjust dims if necessary
    available_pcs <- ncol(object@reductions$pca@cell.embeddings)
    if (max(dims) > available_pcs) {
      warning(paste0(
        "Requested dims (max=", max(dims), ") exceeds available PCs (", 
        available_pcs, "). Adjusting to use available dimensions."
      ))
      dims <- dims[dims <= available_pcs]
    }
    
    if (!is.numeric(seed) || length(seed) != 1) {
      stop("'seed' must be a single integer value")
    }
    
    # Set R's global seed for reproducibility
    set.seed(seed)
    
    if (verbose) {
      message(paste0(
        "Running 3D UMAP with: ",
        "dims=1:", max(dims), 
        " | seed=", seed,
        " | n.neighbors=", n.neighbors,
        " | min.dist=", min.dist,
        " | metric=", metric
      ))
    }
    
    # Run 3D UMAP with explicit parameters for reproducibility
    object <- RunUMAP(
      object = object,
      dims = dims,
      seed.use = seed,
      n.neighbors = n.neighbors,
      min.dist = min.dist,
      metric = metric,
      n.components = 3L,
      reduction.name = "umap3d",
      verbose = verbose,
      ...
    )
    
    # Store metadata about 3D UMAP computation for traceability
    object@misc$umap3d_params <- list(
      dims = dims,
      seed = seed,
      n.neighbors = n.neighbors,
      min.dist = min.dist,
      metric = metric,
      n.components = 3L,
      computation_date = Sys.time()
    )
    
    if (verbose) {
      message("3D UMAP computation completed successfully")
    }
    
    return(object)
    
  }, error = function(e) {
    message(paste("ERROR in runUMAP3D_reproducible:", e$message))
    stop(e)
  })
}


hasValid3DUMAP <- function(object, min_cells = 10) {
  # Check if 3D UMAP reduction exists and is valid
  # Args:
  #   object: Seurat object to check
  #   min_cells: Minimum number of cells expected (default: 10)
  # Returns:
  #   Logical TRUE if valid 3D UMAP exists, FALSE otherwise
  
  tryCatch({
    if (!inherits(object, "Seurat")) {
      warning("Input is not a Seurat object")
      return(FALSE)
    }
    
    # Check if 3D UMAP reduction exists
    if (!"umap3d" %in% names(object@reductions)) {
      return(FALSE)
    }
    
    # Check UMAP dimensions
    umap_coords <- object@reductions$umap3d@cell.embeddings
    
    if (is.null(umap_coords) || nrow(umap_coords) < min_cells) {
      warning(paste0(
        "3D UMAP exists but has insufficient cells: ",
        nrow(umap_coords), " < ", min_cells
      ))
      return(FALSE)
    }
    
    # Check for proper dimensions (must be 3D)
    if (ncol(umap_coords) != 3) {
      warning(paste0(
        "3D UMAP has unexpected dimensions: ",
        ncol(umap_coords), " (expected 3)"
      ))
      return(FALSE)
    }
    
    # Check for NA values
    if (any(is.na(umap_coords))) {
      warning("3D UMAP contains NA values")
      return(FALSE)
    }
    
    return(TRUE)
    
  }, error = function(e) {
    message(paste("ERROR in hasValid3DUMAP:", e$message))
    return(FALSE)
  })
}
