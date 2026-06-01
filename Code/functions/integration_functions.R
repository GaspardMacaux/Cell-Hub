# integration_functions_newarch.R


############################## Seurat v5 Layer Management ##############################

prepareSeuratV5ForIntegration <- function(seurat_object) {
  # Prepare Seurat v5 object for integration by managing layers properly
  # Args:
  #   seurat_object: Seurat object to prepare
  # Returns:
  #   Prepared Seurat object with proper layer structure
  
  if (packageVersion("Seurat") >= "5.0.0") {
    message("Preparing Seurat v5 object for integration...")
    
    tryCatch({
      # Join layers in RNA assay if multiple layers exist
      if ("RNA" %in% names(seurat_object@assays)) {
        rna_layers <- names(seurat_object[["RNA"]]@layers)
        if (length(rna_layers) > 1) {
          message(paste("Found", length(rna_layers), "RNA layers, joining..."))
          seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
        }
      }
      
      # Ensure proper layer structure
      if ("RNA" %in% names(seurat_object@assays)) {
        existing <- names(seurat_object[["RNA"]]@layers)
        # Fallback counts from data only — never from scale.data or unknown layer order
        if (!"counts" %in% existing && "data" %in% existing) {
          message("Adding counts layer from data...")
          seurat_object[["RNA"]]@layers$counts <- seurat_object[["RNA"]]@layers[["data"]]
        }
        # Fallback data from counts only — never from scale.data or unknown layer order
        if (!"data" %in% existing && "counts" %in% existing) {
          message("Adding data layer from counts...")
          seurat_object[["RNA"]]@layers$data <- seurat_object[["RNA"]]@layers[["counts"]]
        }
      }
      
      message("Seurat v5 object prepared successfully")
      
    }, error = function(e) {
      message(paste("Warning: Could not prepare Seurat v5 object:", e$message))
    })
  }
  
  return(seurat_object)
}

cleanIntegratedSeuratV5 <- function(seurat_object) {
  # Clean integrated Seurat v5 object to resolve layer issues
  # Args:
  #   seurat_object: Integrated Seurat object
  # Returns:
  #   Cleaned Seurat object
  
  if (packageVersion("Seurat") >= "5.0.0") {
    message("Cleaning integrated Seurat v5 object...")
    
    tryCatch({
      # Join layers in all assays
      for (assay_name in names(seurat_object@assays)) {
        assay_layers <- names(seurat_object[[assay_name]]@layers)
        if (length(assay_layers) > 1) {
          message(paste("Joining layers in", assay_name, "assay"))
          seurat_object[[assay_name]] <- JoinLayers(seurat_object[[assay_name]])
        }
      }
      
      message("Seurat v5 object cleaned successfully")
      
    }, error = function(e) {
      message(paste("Warning: Could not clean Seurat v5 object:", e$message))
    })
  }
  
  return(seurat_object)
}



############################## Data Preprocessing Functions ##############################

preprocessRawDataset <- function(file_path, dataset_type, species,
                                 dataset_name, qc_params,
                                 normalization_method = "LogNormalize",
                                 log_function = message) {
  # Preprocess a single raw dataset for integration
  # Args:
  #   file_path: Path to the uploaded file (ZIP, H5AD, or RDS)
  #   dataset_type: "snRNA_merge", "multiome_merge", or "seurat_object_merge"
  #   species: Species for mitochondrial pattern
  #   dataset_name: Name for the dataset
  #   qc_params: List with min_features_merge, max_features_merge, max_mt_percent_merge
  #   normalization_method: "LogNormalize", "SCTransform", "CLR", or "RC"
  #   log_function: Function to call for logging
  # Returns:
  #   Processed Seurat object ready for integration
  log_function(paste("Starting preprocessing for dataset:", dataset_name), "info")
  mt_pattern <- switch(species,
                       "mouse" = "^mt-",
                       "human" = "^MT-",
                       "rat" = "^Mt-",
                       "^mt-")
  log_function(paste("Mitochondrial pattern:", mt_pattern), "info")
  log_function("QC Parameters:", "info")
  log_function(paste("  - Min features:", qc_params$min_features_merge), "info")
  log_function(paste("  - Max features:", qc_params$max_features_merge), "info")
  log_function(paste("  - Max MT%:", qc_params$max_mt_percent_merge), "info")
  log_function(paste("  - Normalization method:", normalization_method), "info")
  seurat_object <- NULL
  file_extension <- tolower(tools::file_ext(file_path))
  tryCatch({
    if (dataset_type == "seurat_object_merge") {
      # Processed object path (RDS or H5AD)
      log_function("Loading pre-processed object", "info")
      seurat_object <- loadSeuratObject(
        rds_path = file_path,
        add_dataset_column = FALSE,
        file_format = tolower(tools::file_ext(file_path)),
        log_function = log_function
      )
      if (!inherits(seurat_object, "Seurat")) {
        stop("Loaded file is not a valid Seurat object")
      }
      DefaultAssay(seurat_object) <- "RNA"
      cells_before <- ncol(seurat_object)
      genes_before <- nrow(seurat_object)
      log_function(paste("Cells before QC:", format(cells_before, big.mark = ",")), "info")
      log_function(paste("Genes:", format(genes_before, big.mark = ",")), "info")
      if ("dataset" %in% colnames(seurat_object@meta.data)) {
        original_datasets <- unique(seurat_object@meta.data$dataset)
        log_function(paste("✓ Object contains datasets:",
                           paste(original_datasets, collapse = ", ")), "info")
        log_function(paste("  Total:", length(original_datasets), "original datasets"), "info")
      } else {
        log_function("⚠ No 'dataset' column found in object", "warning")
        log_function("  Creating one based on project name", "info")
        seurat_object@meta.data$dataset <- seurat_object@project.name
      }
      seurat_object@meta.data$integration_group <- dataset_name
      log_function(paste("✓ Created 'integration_group' column with value:", dataset_name), "success")
      log_function("✓ Preserved 'dataset' column with original sample names", "success")
      if (!"integration_level" %in% colnames(seurat_object@meta.data)) {
        seurat_object@meta.data$integration_level <- 1
      }
      seurat_object@meta.data$integration_level <- seurat_object@meta.data$integration_level + 1
      log_function(paste("✓ Integration level:",
                         unique(seurat_object@meta.data$integration_level)), "info")
      seurat_object <- validateAndCompleteSeurat(seurat_object, log_function)
    } else if (file_extension == "h5ad") {
      # Raw H5AD path — load counts via schard, then run full preprocessing
      log_function("Processing raw H5AD (AnnData) data", "info")
      seurat_object <- loadFromH5AD(
        h5ad_path = file_path,
        dataset_type = gsub("_merge$", "", dataset_type),
        mt_pattern = mt_pattern,
        log_function = log_function
      )
      # Strip imported reductions — user will re-compute via integration pipeline
      if (length(seurat_object@reductions) > 0) {
        log_function(paste("Clearing imported reductions:",
                           paste(names(seurat_object@reductions), collapse = ", ")), "info")
        for (red_name in names(seurat_object@reductions)) {
          seurat_object[[red_name]] <- NULL
        }
      }
      seurat_object@project.name <- dataset_name
      cells_before <- ncol(seurat_object)
      genes_before <- nrow(seurat_object)
      log_function(paste("Cells loaded:", format(cells_before, big.mark = ",")), "info")
      log_function(paste("Genes:", format(genes_before, big.mark = ",")), "info")
      seurat_object@meta.data$integration_group <- dataset_name
      # Normalization + variable features + scaling + PCA
      log_function("Normalizing and preprocessing...", "info")
      if (normalization_method == "SCTransform") {
        log_function("Using SCTransform for normalization", "info")
        seurat_object <- SCTransform(seurat_object, variable.features.n = 4000, verbose = FALSE)
        log_function("SCTransform completed", "success")
        DefaultAssay(seurat_object) <- "SCT"
        seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object),
                                npcs = 50, verbose = FALSE)
        log_function(paste("Variable features found:", length(VariableFeatures(seurat_object))), "info")
        log_function(paste("Default assay set to:", DefaultAssay(seurat_object)), "info")
      } else {
        log_function(paste0("Using ", normalization_method, " for normalization"), "info")
        seurat_object <- NormalizeData(seurat_object, normalization.method = normalization_method,
                                       verbose = FALSE)
        seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst",
                                              nfeatures = 4000, verbose = FALSE)
        seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object), verbose = FALSE)
        seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object),
                                verbose = FALSE)
        DefaultAssay(seurat_object) <- "RNA"
        log_function(paste("Variable features found:", length(VariableFeatures(seurat_object))), "info")
      }
      log_function("Preprocessing completed", "success")
    } else {
      # Raw 10X ZIP path (existing logic)
      log_function("Processing raw 10X data from ZIP", "info")
      temp_dir <- createTempDirectory(paste0("dataset_", dataset_name))
      log_function("Extracting ZIP file...", "info")
      unzip(file_path, exdir = temp_dir)
      log_function("Searching for 10X files in extracted directory...", "info")
      valid_folder <- find10XFolder(temp_dir)
      if (is.null(valid_folder)) {
        stop("No complete set of 10X files found in ZIP. Please ensure the ZIP contains barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz files.")
      }
      log_function(paste("Found 10X files in:", valid_folder), "info")
      log_function("Reading 10X data...", "info")
      data_10x <- Read10X(valid_folder)
      if (is.list(data_10x) && "Gene Expression" %in% names(data_10x)) {
        cells_before <- ncol(data_10x$`Gene Expression`)
        genes_before <- nrow(data_10x$`Gene Expression`)
      } else {
        cells_before <- ncol(data_10x)
        genes_before <- nrow(data_10x)
      }
      log_function(paste("Raw matrix dimensions:", format(genes_before, big.mark = ","),
                         "genes x", format(cells_before, big.mark = ","), "barcodes"), "info")
      log_function("Creating Seurat object with initial filtering...", "info")
      if (dataset_type == "snRNA_merge") {
        seurat_object <- CreateSeuratObject(counts = data_10x, project = dataset_name,
                                            min.cells = 3, min.features = qc_params$min_features_merge)
      } else if (dataset_type == "multiome_merge") {
        if (is.list(data_10x) && "Gene Expression" %in% names(data_10x)) {
          log_function("Extracting Gene Expression assay from multiome data", "info")
          seurat_object <- CreateSeuratObject(counts = data_10x$`Gene Expression`, project = dataset_name,
                                              min.cells = 3, min.features = qc_params$min_features_merge)
        } else {
          stop("Could not find Gene Expression data in multiome file")
        }
      }
      log_function(paste("Cells after initial filtering:", format(ncol(seurat_object), big.mark = ",")), "info")
      unlink(temp_dir, recursive = TRUE)
      seurat_object@meta.data$integration_group <- dataset_name
      # Normalization + variable features + scaling + PCA
      log_function("Normalizing and preprocessing...", "info")
      if (normalization_method == "SCTransform") {
        log_function("Using SCTransform for normalization", "info")
        seurat_object <- SCTransform(seurat_object, variable.features.n = 4000, verbose = FALSE)
        log_function("SCTransform completed", "success")
        DefaultAssay(seurat_object) <- "SCT"
        seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object),
                                npcs = 50, verbose = FALSE)
        log_function(paste("Variable features found:", length(VariableFeatures(seurat_object))), "info")
        log_function(paste("Default assay set to:", DefaultAssay(seurat_object)), "info")
      } else {
        log_function(paste0("Using ", normalization_method, " for normalization"), "info")
        seurat_object <- NormalizeData(seurat_object, normalization.method = normalization_method,
                                       verbose = FALSE)
        seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst",
                                              nfeatures = 4000, verbose = FALSE)
        seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object), verbose = FALSE)
        seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object),
                                verbose = FALSE)
        DefaultAssay(seurat_object) <- "RNA"
        log_function(paste("Variable features found:", length(VariableFeatures(seurat_object))), "info")
      }
      log_function("Preprocessing completed", "success")
    }
    seurat_object <- applyQCAndPreprocessing(
      seurat_object = seurat_object,
      mt_pattern = mt_pattern,
      qc_params = qc_params,
      dataset_name = dataset_name,
      is_preprocessed = TRUE,
      log_function = log_function
    )
    log_function(paste("Dataset", dataset_name, "preprocessing completed successfully"), "success")
    return(seurat_object)
  }, error = function(e) {
    if (exists("temp_dir") && dir.exists(temp_dir)) {
      unlink(temp_dir, recursive = TRUE)
    }
    log_function(paste("Error preprocessing dataset:", conditionMessage(e)[1]), "error")
    stop(paste("Error preprocessing dataset", dataset_name, ":", conditionMessage(e)[1]))
  })
}

applyQCAndPreprocessing <- function(seurat_object, mt_pattern, qc_params, 
                                    dataset_name, is_preprocessed = FALSE,
                                    log_function = message) {
  # Apply QC filtering and standard preprocessing steps with v5 compatibility
  # Args:
  #   seurat_object: Seurat object to process
  #   mt_pattern: Mitochondrial gene pattern
  #   qc_params: List with min_features_merge, max_features_merge, max_mt_percent_merge
  #   dataset_name: Dataset identifier
  #   is_preprocessed: Whether object is already processed
  #   log_function: Function for logging
  # Returns:
  #   QC-filtered and preprocessed Seurat object
  
  log_function("Applying QC and preprocessing steps", "info")
  
  # Prepare for Seurat v5 if needed
  seurat_object <- prepareSeuratV5ForIntegration(seurat_object)
  
  # Add mitochondrial percentage if not present
  if (!"percent.mt" %in% colnames(seurat_object@meta.data)) {
    log_function("Calculating mitochondrial percentage...", "info")
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = mt_pattern)
  }
  
  # Log pre-filter statistics
  pre_filter_cells <- ncol(seurat_object)
  log_function(paste("Pre-QC cell count:", format(pre_filter_cells, big.mark = ",")), "info")
  
  # Apply QC filters using UI parameters
  log_function("Applying QC filters...", "info")
  seurat_object <- subset(seurat_object,
                          subset = nFeature_RNA > qc_params$min_features_merge &
                            nFeature_RNA < qc_params$max_features_merge &
                            percent.mt < qc_params$max_mt_percent_merge)
  
  post_filter_cells <- ncol(seurat_object)
  pct_retained <- round(post_filter_cells/pre_filter_cells * 100, 2)
  
  log_function(paste("Post-QC cell count:", format(post_filter_cells, big.mark = ",")), "success")
  log_function(paste("Retained:", pct_retained, "% of cells"), "success")
  log_function(paste("Removed:", format(pre_filter_cells - post_filter_cells, big.mark = ","), "cells"), "info")
  
  # Validate that cells remain after filtering
  if (post_filter_cells == 0) {
    log_function("ERROR: No cells remain after QC filtering!", "error")
    log_function("Consider relaxing QC parameters", "warning")
    stop(paste("No cells remain after QC filtering for dataset", dataset_name, 
               "- consider relaxing QC parameters"))
  }
  
  # Apply standard preprocessing if not already done
  if (!is_preprocessed) {
    log_function("Applying normalization, scaling, and PCA", "info")
    
    # Use layer-aware functions for Seurat v5
    log_function("  - Normalizing data...", "info")
    seurat_object <- NormalizeData(seurat_object, verbose = FALSE)
    
    log_function("  - Finding variable features...", "info")
    seurat_object <- FindVariableFeatures(seurat_object,
                                          selection.method = "vst",
                                          nfeatures = 4000,
                                          verbose = FALSE)
    
    log_function("  - Scaling data...", "info")
    seurat_object <- ScaleData(seurat_object,
                               features = rownames(seurat_object),
                               verbose = FALSE)
    
    log_function("  - Running PCA...", "info")
    seurat_object <- RunPCA(seurat_object,
                            features = VariableFeatures(seurat_object),
                            verbose = FALSE)
    
    log_function("Preprocessing completed", "success")
  } else {
    log_function("Object already preprocessed, skipping normalization steps", "info")
  }
  
  # Add dataset metadata
  seurat_object$dataset <- dataset_name
  seurat_object$orig.ident <- dataset_name
  
  # Log final QC metrics
  if ("nFeature_RNA" %in% colnames(seurat_object@meta.data)) {
    log_function(paste("Median features/cell:", round(median(seurat_object$nFeature_RNA))), "info")
  }
  if ("nCount_RNA" %in% colnames(seurat_object@meta.data)) {
    log_function(paste("Median UMI/cell:", round(median(seurat_object$nCount_RNA))), "info")
  }
  if ("percent.mt" %in% colnames(seurat_object@meta.data)) {
    log_function(paste("Median MT%:", round(median(seurat_object$percent.mt), 2), "%"), "info")
  }
  
  return(seurat_object)
}

############################## Batch Processing Function ##############################

processDatasetsForIntegration <- function(file_inputs, dataset_types, dataset_names, 
                                          species, qc_params) {
  # Process multiple datasets for integration
  # Args:
  #   file_inputs: List of file paths
  #   dataset_types: Vector of dataset types
  #   dataset_names: Vector of dataset names
  #   species: Species identifier
  #   qc_params: QC parameters from UI
  # Returns:
  #   List of processed Seurat objects
  
  if (length(file_inputs) != length(dataset_types) || 
      length(file_inputs) != length(dataset_names)) {
    stop("Mismatch in number of files, types, and names")
  }
  
  message(paste("Processing", length(file_inputs), "datasets for integration"))
  
  seurat_list <- list()
  
  for (i in 1:length(file_inputs)) {
    message(paste("Processing dataset", i, "of", length(file_inputs)))
    
    processed_object <- preprocessRawDataset(
      file_path = file_inputs[[i]],
      dataset_type = dataset_types[i],
      species = species,
      dataset_name = dataset_names[i],
      qc_params = qc_params
    )
    
    seurat_list[[i]] <- processed_object
  }
  
  # Validate all objects
  validateSeuratList(seurat_list)
  
  message("All datasets processed successfully")
  return(seurat_list)
}

############################## Integration Functions ##############################

performDataIntegration <- function(seurat_list, integration_method = "standard", 
                                   harmony_vars = NULL, harmony_dims = NULL,
                                   mnn_dims = NULL, mnn_k = NULL) {
  # Perform data integration using specified method
  # Args:
  #   seurat_list: List of preprocessed Seurat objects
  #   integration_method: "standard", "harmony", "fastmnn", or "simple"
  #   harmony_vars: Variables to correct for with Harmony (optional)
  #   harmony_dims: Number of dimensions for Harmony (optional)
  #   mnn_dims: Number of corrected dimensions for fastMNN (optional, default 50)
  #   mnn_k: Number of MNN neighbors for fastMNN (optional, default 20)
  # Returns:
  #   Integrated Seurat object
  
  if (length(seurat_list) < 2) {
    stop("At least two datasets are required for integration")
  }
  
  message(paste("Starting", integration_method, "integration with", 
                length(seurat_list), "datasets"))
  
  validateSeuratList(seurat_list)
  
  if (integration_method == "simple") {
    return(performSimpleMerge(seurat_list))
    
  } else if (integration_method == "harmony") {
    if (is.null(harmony_vars)) harmony_vars <- "dataset"
    if (is.null(harmony_dims)) harmony_dims <- 30
    return(performHarmonyIntegration(seurat_list, 
                                     harmony_vars = harmony_vars,
                                     harmony_dims = harmony_dims))
    
  } else if (integration_method == "fastmnn") {
    if (is.null(mnn_dims)) mnn_dims <- 50
    if (is.null(mnn_k)) mnn_k <- 20
    return(performFastMNNIntegration(seurat_list,
                                     mnn_dims = mnn_dims,
                                     mnn_k = mnn_k))
    
  } else if (integration_method == "standard") {
    return(performStandardIntegration(seurat_list))
    
  } else {
    stop(paste("Unknown integration method:", integration_method))
  }
}


performStandardIntegration <- function(seurat_list) {
  # Perform standard Seurat integration workflow with v5 compatibility
  # Args:
  #   seurat_list: List of preprocessed Seurat objects
  # Returns:
  #   Integrated Seurat object with "integrated" assay
  
  message("Starting standard Seurat integration")
  
  # Prepare objects for Seurat v5 integration
  seurat_list <- lapply(seurat_list, prepareSeuratV5ForIntegration)
  
  # Calculate dataset statistics for parameter adjustment
  dataset_sizes <- sapply(seurat_list, function(obj) ncol(obj))
  min_cells <- min(dataset_sizes)
  total_cells <- sum(dataset_sizes)
  
  message(paste("Dataset sizes:", paste(dataset_sizes, collapse = ", ")))
  message(paste("Minimum cells:", min_cells, "Total cells:", total_cells))
  
  # Adjust integration parameters based on dataset sizes
  integration_params <- calculateIntegrationParameters(min_cells)
  
  message(paste("Using integration parameters:",
                "k.filter =", integration_params$k_filter,
                "k.score =", integration_params$k_score,
                "k.anchor =", integration_params$k_anchor,
                "k.weight =", integration_params$k_weight))
  
  # Select integration features
  message("Selecting integration features")
  features <- SelectIntegrationFeatures(
    object.list = seurat_list,
    nfeatures = 2000
  )
  
  # Find integration anchors
  message("Finding integration anchors")
  anchors <- FindIntegrationAnchors(
    object.list = seurat_list,
    dims = 1:30,
    k.filter = integration_params$k_filter,
    k.score = integration_params$k_score,
    k.anchor = integration_params$k_anchor,
    anchor.features = features
  )
  
  # Integrate data
  message("Integrating datasets")
  integrated_object <- IntegrateData(
    anchorset = anchors,
    dims = 1:30,
    features.to.integrate = features,
    k.weight = integration_params$k_weight
  )
  
  # Clean integrated object for v5 compatibility
  integrated_object <- cleanIntegratedSeuratV5(integrated_object)
  
  # Set default assay and clean up metadata
  DefaultAssay(integrated_object) <- "integrated"
  integrated_object <- cleanIntegratedMetadata(integrated_object)
  
  message("Standard integration completed successfully")
  return(integrated_object)
}


performSimpleMerge <- function(seurat_list, dataset_names = NULL) {
  # Perform simple merge without integration, preserving original clusters
  # Args:
  #   seurat_list: List of preprocessed Seurat objects (can be named list)
  #   dataset_names: Optional. If NULL, will auto-detect from objects or list names
  # Returns:
  #   Merged Seurat object with preserved cluster and dataset information + reductions
  
  message("Starting simple merge without integration")
  
  if (is.null(dataset_names)) {
    message("Auto-detecting dataset names...")
    
    if (!is.null(names(seurat_list)) && all(names(seurat_list) != "")) {
      dataset_names <- names(seurat_list)
      message(paste("✓ Using list names:", paste(dataset_names, collapse = ", ")))
    } 
    else {
      dataset_names <- sapply(seurat_list, function(obj) {
        if (!is.null(obj@project.name) && obj@project.name != "") {
          return(obj@project.name)
        } else {
          return(NULL)
        }
      })
      
      if (!any(is.null(dataset_names))) {
        message(paste("✓ Using project names:", paste(dataset_names, collapse = ", ")))
      } else {
        dataset_names <- paste0("Dataset_", 1:length(seurat_list))
        message(paste("⚠ No names found, using default:", paste(dataset_names, collapse = ", ")))
      }
    }
  } else {
    message(paste("✓ Using provided names:", paste(dataset_names, collapse = ", ")))
  }
  
  if (length(dataset_names) != length(seurat_list)) {
    stop(paste("Number of dataset names (", length(dataset_names), 
               ") must match number of Seurat objects (", length(seurat_list), ")"))
  }
  
  for (i in 1:length(seurat_list)) {
    dataset_name <- dataset_names[i]
    
    message(paste("\nProcessing", dataset_name, "..."))
    
    if ("integration_group" %in% colnames(seurat_list[[i]]@meta.data)) {
      message("  ✓ Pre-integrated object detected (has integration_group column)")
      
      if ("dataset" %in% colnames(seurat_list[[i]]@meta.data)) {
        original_samples <- unique(seurat_list[[i]]@meta.data$dataset)
        message(paste("  ✓ Preserving original samples:", 
                      paste(head(original_samples, 5), collapse = ", "),
                      if(length(original_samples) > 5) "..." else ""))
      }
      
    } else {
      message("  ✓ First-time integration (raw data)")
      seurat_list[[i]]$dataset <- dataset_name
      seurat_list[[i]]$integration_group <- dataset_name
    }
    
    if (!"dataset_origin" %in% colnames(seurat_list[[i]]@meta.data)) {
      seurat_list[[i]]$dataset_origin <- dataset_name
    }
    
    existing_meta <- colnames(seurat_list[[i]]@meta.data)
    message(paste("  Metadata columns:", length(existing_meta)))
    
    cluster_names <- NULL
    
    if ("ClusterIdents" %in% existing_meta) {
      cluster_names <- as.character(seurat_list[[i]]$ClusterIdents)
      message("  ✓ Using ClusterIdents column for cluster names")
    }
    else if ("annotated_clusters" %in% existing_meta) {
      cluster_names <- as.character(seurat_list[[i]]$annotated_clusters)
      message("  ✓ Using annotated_clusters column for cluster names")
    }
    else if (!all(grepl("^[0-9]+$", as.character(Idents(seurat_list[[i]]))))) {
      cluster_names <- as.character(Idents(seurat_list[[i]]))
      message("  ✓ Using Idents for cluster names")
    }
    else if ("seurat_clusters" %in% existing_meta) {
      cluster_names <- paste0("C", seurat_list[[i]]$seurat_clusters)
      message("  ⚠ Using numbered clusters (no annotation found)")
    }
    else {
      cluster_names <- as.character(Idents(seurat_list[[i]]))
      message("  ⚠ Using default identities")
    }
    
    seurat_list[[i]]$original_clusters <- paste0(dataset_name, "_", cluster_names)
    seurat_list[[i]]$cluster_name_only <- cluster_names
    
    n_unique_clusters <- length(unique(cluster_names))
    message(paste("  ✓ Preserved", n_unique_clusters, "unique cluster names"))
    message(paste("    Example clusters:", paste(head(unique(cluster_names), 3), collapse = ", ")))
    
    if ("seurat_clusters" %in% existing_meta) {
      seurat_list[[i]]$original_seurat_clusters <- as.character(seurat_list[[i]]$seurat_clusters)
    }
    
    if ("cell_type" %in% existing_meta) {
      seurat_list[[i]]$original_cell_type <- as.character(seurat_list[[i]]$cell_type)
      message("  ✓ Preserved cell_type column")
    }
    
    if ("ClusterIdents" %in% existing_meta) {
      seurat_list[[i]]$original_ClusterIdents <- as.character(seurat_list[[i]]$ClusterIdents)
      message("  ✓ Preserved ClusterIdents column")
    }
    
    custom_cols <- setdiff(existing_meta, c("nCount_RNA", "nFeature_RNA", "orig.ident"))
    if (length(custom_cols) > 0) {
      message(paste("  ✓ Preserving", length(custom_cols), "additional metadata columns"))
    }
  }
  
  message("\nMerging objects...")
  merged_object <- merge(
    seurat_list[[1]], 
    y = seurat_list[-1],
    add.cell.ids = dataset_names,
    project = "Merged_Datasets",
    merge.data = TRUE
  )
  
  message("\n========== CHECKING SEURAT VERSION ==========")
  seurat_version <- packageVersion("Seurat")
  message(paste("Seurat version:", seurat_version))
  
  if (seurat_version >= "5.0.0") {
    message("Seurat v5 detected, checking for layers...")
    
    has_layers <- tryCatch({
      test_assay <- merged_object[["RNA"]]
      "layers" %in% slotNames(test_assay)
    }, error = function(e) {
      FALSE
    })
    
    if (has_layers) {
      message("✓ Layers slot found, joining layers...")
      tryCatch({
        for (assay_name in names(merged_object@assays)) {
          assay_layers <- names(merged_object[[assay_name]]@layers)
          if (length(assay_layers) > 1) {
            message(paste("  Joining", length(assay_layers), "layers in", assay_name, "assay"))
            merged_object[[assay_name]] <- JoinLayers(merged_object[[assay_name]])
            message(paste("  ✓", assay_name, "layers joined"))
          } else {
            message(paste("  ✓", assay_name, "has single layer, no joining needed"))
          }
        }
        message("✓ All layers joined successfully")
      }, error = function(e) {
        warning(paste("Could not join layers:", e$message, 
                      "\nContinuing without joining - may cause issues with some functions"))
      })
    } else {
      message("✓ No layers slot found (object created with older Seurat version)")
      message("  This is normal - continuing without layer joining")
    }
  } else {
    message("✓ Seurat v4 or earlier detected - no layer joining needed")
  }
  message("=============================================\n")
  
  if ("cluster_name_only" %in% colnames(merged_object@meta.data)) {
    Idents(merged_object) <- "cluster_name_only"
    message("✓ Set identities to cluster names")
  } else if ("original_clusters" %in% colnames(merged_object@meta.data)) {
    Idents(merged_object) <- "original_clusters"
    message("✓ Set identities to original clusters")
  } else {
    Idents(merged_object) <- "dataset"
    message("✓ Set identities to dataset")
  }
  
  merged_object$merge_method <- "simple_merge"
  
  message("\n========== MERGE SUMMARY ==========")
  message(paste("Total cells:", ncol(merged_object)))
  message(paste("Total features:", nrow(merged_object)))
  message(paste("\nDatasets merged:", length(unique(merged_object$dataset))))
  
  dataset_counts <- table(merged_object$dataset)
  for (ds in names(dataset_counts)) {
    message(paste("  -", ds, ":", dataset_counts[ds], "cells"))
  }
  
  if ("cluster_name_only" %in% colnames(merged_object@meta.data)) {
    unique_clusters <- unique(merged_object$cluster_name_only)
    message(paste("\nUnique cluster names:", length(unique_clusters)))
    message(paste("Clusters:", paste(head(unique_clusters, 10), collapse = ", "), 
                  if(length(unique_clusters) > 10) "..." else ""))
  }
  
  message("===================================\n")
  
  message("\n========== COMPUTING REDUCTIONS ==========")
  message("Running normalization pipeline...")
  
  DefaultAssay(merged_object) <- "RNA"
  
  message("  [1/5] Normalizing...")
  merged_object <- NormalizeData(merged_object, verbose = FALSE)
  
  message("  [2/5] Finding variable features...")
  merged_object <- FindVariableFeatures(merged_object, 
                                        selection.method = "vst",
                                        nfeatures = 2000,
                                        verbose = FALSE)
  
  message("  [3/5] Scaling data...")
  merged_object <- ScaleData(merged_object, verbose = FALSE)
  
  message("  [4/5] Computing PCA...")
  merged_object <- runPCA_reproducible(
    object = merged_object,
    npcs = 30,
    seed = 42,
    verbose = FALSE
  )
  
  message("  [5/5] Computing UMAP...")
  merged_object <- runUMAP_reproducible(
    object = merged_object,
    dims = 1:20,
    seed = 42,
    verbose = FALSE
  )
  
  message("✓ Reductions computed successfully")
  message("==========================================\n")
  
  if (!("umap" %in% names(merged_object@reductions))) {
    stop("ERROR: UMAP was not created. Object cannot be plotted.")
  }
  
  message("✓ Merged object ready for visualization!")
  message("\nPlotting options:")
  message("  DimPlot(merged, reduction = 'umap', group.by = 'dataset')")
  message("  DimPlot(merged, reduction = 'umap', group.by = 'cluster_name_only')")
  message("  DotPlot(merged, features = c('gene1', 'gene2'))")
  
  return(merged_object)
}



performHarmonyIntegration <- function(seurat_list, harmony_vars = "dataset", 
                                      harmony_dims = 30) {
  n_datasets <- length(seurat_list)
  message(paste("Starting Harmony with", n_datasets, "datasets"))
  
  gc()
  
  message("Step 1/3: Merging datasets")
  merged_object <- performSimpleMerge(seurat_list)
  
  rm(seurat_list)
  gc()
  
  message("Step 2/3: Normalizing and preprocessing")
  
  merged_object <- NormalizeData(merged_object, verbose = FALSE)
  
  n_cells <- ncol(merged_object)
  n_features <- if (n_cells > 150000) {
    2000
  } else if (n_cells > 100000) {
    2500
  } else {
    3000
  }
  
  message(paste("Using", n_features, "variable features (optimized for 48 GB RAM)"))
  
  merged_object <- FindVariableFeatures(
    merged_object,
    selection.method = "vst",
    nfeatures = n_features,
    verbose = FALSE
  )
  
  message("Scaling ONLY variable features (memory optimization)")
  merged_object <- ScaleData(
    merged_object,
    features = VariableFeatures(merged_object),
    verbose = FALSE
  )
  
  gc()
  
  message("Running PCA")
  merged_object <- RunPCA(
    merged_object,
    features = VariableFeatures(merged_object),
    npcs = 50,
    verbose = FALSE
  )
  
  gc()
  
  message("Step 3/3: Running Harmony batch correction")
  
  available_vars <- colnames(merged_object@meta.data)
  valid_vars <- harmony_vars[harmony_vars %in% available_vars]
  
  if (length(valid_vars) == 0) {
    warning("Harmony variables not found, using 'dataset'")
    valid_vars <- "dataset"
  }
  
  vars_to_remove <- c()
  for (var in valid_vars) {
    unique_levels <- unique(merged_object@meta.data[[var]])
    unique_levels <- unique_levels[!is.na(unique_levels)]
    n_levels <- length(unique_levels)
    
    if (n_levels < 2) {
      message(paste("⚠ WARNING: Variable", var, "has only", n_levels, "level(s)"))
      message(paste("    Values:", paste(unique_levels, collapse = ", ")))
      vars_to_remove <- c(vars_to_remove, var)
    } else {
      message(paste("✓ Variable", var, "has", n_levels, "levels:", 
                    paste(head(unique_levels, 5), collapse = ", "),
                    if(n_levels > 5) "..." else ""))
    }
  }
  
  if (length(vars_to_remove) > 0) {
    valid_vars <- setdiff(valid_vars, vars_to_remove)
    message(paste("Removed", length(vars_to_remove), "variable(s) with insufficient levels"))
  }
  
  if (length(valid_vars) == 0) {
    stop(paste(
      "ERROR: No valid batch correction variables found.",
      "\nAll specified variables have less than 2 levels.",
      "\nHarmony requires at least one variable with 2+ levels to correct batch effects.",
      "\n\nSuggestions:",
      "\n  1. Check that your datasets are actually different (not the same dataset loaded 4 times)",
      "\n  2. Verify metadata columns have multiple values across datasets",
      "\n  3. Try 'Simple Merge' if no batch correction is needed",
      "\n\nAvailable metadata columns:", paste(available_vars, collapse = ", ")
    ))
  }
  
  dims_use <- min(harmony_dims, 30)
  max_iter <- 20
  
  message(paste("Harmony parameters: dims =", dims_use, 
                ", max.iter =", max_iter,
                ", group.by.vars =", paste(valid_vars, collapse = ", ")))
  
  if (.Platform$OS.type == "windows") {
    message(paste("RAM before Harmony:", round(memory.size()/1024, 1), "GB"))
  }
  
  merged_object <- RunHarmony(
    merged_object,
    group.by.vars = valid_vars,
    dims.use = 1:dims_use,
    max.iter.harmony = max_iter,
    verbose = TRUE,
    plot_convergence = FALSE
  )
  
  gc()
  
  if (.Platform$OS.type == "windows") {
    message(paste("RAM after Harmony:", round(memory.size()/1024, 1), "GB"))
  }
  
  merged_object$integration_method <- "harmony"
  DefaultAssay(merged_object) <- "RNA"
  
  message(paste("✓ Harmony completed successfully:",
                ncol(merged_object), "cells integrated"))
  
  return(merged_object)
}
############################## Helper Functions ##############################

calculateIntegrationParameters <- function(min_cells) {
  # Calculate appropriate integration parameters based on dataset size
  # Args:
  #   min_cells: Minimum number of cells across datasets
  # Returns:
  #   List of integration parameters
  
  # Conservative parameter calculation to avoid errors
  k_filter <- min(30, max(5, floor(min_cells / 15)))
  k_score <- min(20, max(5, floor(min_cells / 20)))
  k_anchor <- min(5, max(2, floor(min_cells / 50)))
  k_weight <- min(30, max(5, floor(min_cells / 15)))
  
  return(list(
    k_filter = k_filter,
    k_score = k_score,
    k_anchor = k_anchor,
    k_weight = k_weight
  ))
}

cleanIntegratedMetadata <- function(seurat_object) {
  # Clean and standardize metadata for integrated object
  # Args:
  #   seurat_object: Integrated Seurat object
  # Returns:
  #   Seurat object with cleaned metadata
  
  # Standardize dataset column name
  metadata <- seurat_object@meta.data
  if ("orig.ident" %in% colnames(metadata) && !"dataset" %in% colnames(metadata)) {
    colnames(metadata)[colnames(metadata) == "orig.ident"] <- "dataset"
    seurat_object <- AddMetaData(seurat_object, metadata = metadata)
  }
  
  return(seurat_object)
}

validateSeuratList <- function(seurat_list) {
  # Validate that all objects in list are proper Seurat objects
  # Args:
  #   seurat_list: List of objects to validate
  # Returns:
  #   TRUE if all valid, throws error otherwise
  
  if (length(seurat_list) == 0) {
    stop("Empty Seurat object list provided")
  }
  
  for (i in 1:length(seurat_list)) {
    if (is.null(seurat_list[[i]])) {
      stop(paste("Seurat object", i, "is NULL"))
    }
    
    if (!inherits(seurat_list[[i]], "Seurat")) {
      stop(paste("Object", i, "is not a valid Seurat object"))
    }
    
    if (ncol(seurat_list[[i]]) == 0) {
      stop(paste("Seurat object", i, "has no cells"))
    }
  }
  
  message("All Seurat objects validated successfully")
  return(TRUE)
}

extractQCParams <- function(input) {
  # Extract QC parameters from Shiny input
  # Args:
  #   input: Shiny input object
  # Returns:
  #   List of QC parameters
  
  return(list(
    min_features_merge = input$min_features_merge,
    max_features_merge = input$max_features_merge,
    max_mt_percent_merge = input$max_mt_percent_merge
  ))
}

############################## Metadata Management Functions ##############################

addDatasetMetadata <- function(seurat_object, metadata_fields, metadata_values, dataset_name) {
  # Add custom metadata to a specific dataset
  # Args:
  #   seurat_object: Integrated Seurat object
  #   metadata_fields: Vector of metadata field names
  #   metadata_values: Vector of metadata values for this dataset
  #   dataset_name: Name of the dataset to update
  # Returns:
  #   Seurat object with updated metadata
  
  if (length(metadata_fields) != length(metadata_values)) {
    stop("Metadata fields and values must have the same length")
  }
  
  # Find cells belonging to this dataset
  dataset_cells <- which(seurat_object@meta.data$dataset == dataset_name)
  
  if (length(dataset_cells) == 0) {
    warning(paste("No cells found for dataset:", dataset_name))
    return(seurat_object)
  }
  
  # Add metadata for each field
  for (i in 1:length(metadata_fields)) {
    field_name <- metadata_fields[i]
    field_value <- metadata_values[i]
    
    # Initialize field if it doesn't exist
    if (!field_name %in% colnames(seurat_object@meta.data)) {
      seurat_object@meta.data[[field_name]] <- NA
    }
    
    # Set values for this dataset
    seurat_object@meta.data[dataset_cells, field_name] <- field_value
  }
  
  message(paste("Added metadata for dataset:", dataset_name))
  return(seurat_object)
}

processMetadataFromUI <- function(seurat_object, input, num_fields, split_column = "dataset") {
  # Process metadata inputs from UI and add to Seurat object
  # Args:
  #   seurat_object: Integrated Seurat object
  #   input: Shiny input object
  #   num_fields: Number of metadata fields to process
  #   split_column: Column to use for splitting metadata values
  # Returns:
  #   Seurat object with updated metadata
  
  message("=== processMetadataFromUI DEBUG ===")
  message(paste("Split column requested:", split_column))
  message(paste("Available columns:", paste(head(colnames(seurat_object@meta.data), 20), collapse = ", ")))
  
  # CRITICAL FIX: Verify split_column exists
  if (!split_column %in% colnames(seurat_object@meta.data)) {
    warning(paste("Split column", split_column, "not found in metadata"))
    
    # Try to find alternative columns
    possible_cols <- c("dataset", "orig.ident", "dataset_origin", "integration_group")
    found_col <- NULL
    for (col in possible_cols) {
      if (col %in% colnames(seurat_object@meta.data)) {
        found_col <- col
        break
      }
    }
    
    if (!is.null(found_col)) {
      message(paste("Using alternative column:", found_col))
      split_column <- found_col
    } else {
      # No split column found - apply metadata globally
      message("No split column found - applying metadata globally to all cells")
      split_values <- "all_cells"
      
      for (j in 1:num_fields) {
        metadata_field_name <- input[[paste0("metadata_name_", j)]]
        
        if (is.null(metadata_field_name) || metadata_field_name == "") {
          next
        }
        
        # Get value from first available input
        metadata_field_value <- input[[paste0("metadata_value_all_cells_", j)]]
        if (is.null(metadata_field_value)) {
          metadata_field_value <- input[[paste0("metadata_value_1_", j)]]
        }
        
        if (!is.null(metadata_field_value) && metadata_field_value != "") {
          seurat_object@meta.data[[metadata_field_name]] <- metadata_field_value
          message(paste("  Applied", metadata_field_name, "=", metadata_field_value, 
                        "to all", ncol(seurat_object), "cells"))
        }
      }
      
      message("==================================")
      return(seurat_object)
    }
  }
  
  # Get unique values from the split column
  split_values <- unique(seurat_object@meta.data[[split_column]])
  split_values <- split_values[!is.na(split_values)]
  
  if (length(split_values) == 0) {
    stop(paste("No valid values found in split column:", split_column))
  }
  
  message(paste("Processing metadata split by column:", split_column))
  message(paste("Split values:", paste(split_values, collapse = ", ")))
  
  for (j in 1:num_fields) {
    metadata_field_name <- input[[paste0("metadata_name_", j)]]
    
    if (is.null(metadata_field_name) || metadata_field_name == "") {
      next  # Skip empty fields
    }
    
    message(paste("Processing field:", metadata_field_name))
    
    # Initialize field in metadata
    if (!metadata_field_name %in% colnames(seurat_object@meta.data)) {
      seurat_object@meta.data[[metadata_field_name]] <- NA
    }
    
    # Process each split value
    for (split_value in split_values) {
      # Try different input name formats
      metadata_field_value <- input[[paste0("metadata_value_", split_value, "_", j)]]
      if (is.null(metadata_field_value)) {
        metadata_field_value <- input[[paste0("metadata_value_", j, "_", split_value)]]
      }
      
      message(paste("  Looking for input:", paste0("metadata_value_", split_value, "_", j)))
      message(paste("  Value found:", ifelse(is.null(metadata_field_value), "NULL", metadata_field_value)))
      
      if (!is.null(metadata_field_value) && metadata_field_value != "") {
        # Find cells with this split value
        split_cells <- which(seurat_object@meta.data[[split_column]] == split_value)
        
        if (length(split_cells) > 0) {
          seurat_object@meta.data[split_cells, metadata_field_name] <- metadata_field_value
          
          message(paste("  ✓ Assigned", metadata_field_value, "to", length(split_cells), 
                        "cells where", split_column, "=", split_value))
        } else {
          warning(paste("  ⚠ No cells found where", split_column, "=", split_value))
        }
      } else {
        message(paste("  ⚠ No value provided for", split_column, "=", split_value))
      }
    }
  }
  
  message("Processed metadata from UI successfully")
  message("==================================")
  return(seurat_object)
}
############################## Memory Management Functions ##############################

cleanupIntegrationMemory <- function() {
  # Clean up memory after integration process
  # Returns:
  #   Nothing (performs garbage collection)
  
  # Force garbage collection
  gc()
  
  # Clean temporary directories related to integration
  temp_dirs <- list.dirs(tempdir(), full.names = TRUE, recursive = FALSE)
  integration_dirs <- temp_dirs[grepl("dataset_|integration_", basename(temp_dirs))]
  
  for (dir in integration_dirs) {
    if (dir.exists(dir)) {
      unlink(dir, recursive = TRUE)
      message(paste("Cleaned temporary directory:", basename(dir)))
    }
  }
  
  message("Memory cleanup completed")
}

monitorIntegrationProgress <- function(current_step, total_steps, step_name) {
  # Monitor and log integration progress
  # Args:
  #   current_step: Current step number
  #   total_steps: Total number of steps
  #   step_name: Name of current step
  # Returns:
  #   Progress percentage
  
  progress_pct <- round((current_step / total_steps) * 100, 1)
  message(paste("Integration Progress:", progress_pct, "% -", step_name))
  
  return(progress_pct)
}

performFastMNNIntegration <- function(seurat_list, mnn_dims = 50, mnn_k = 20) {
  # Perform fastMNN integration via batchelor on log-normalized counts
  # Args:
  #   seurat_list: List of preprocessed Seurat objects
  #   mnn_dims: Number of dimensions for the corrected embedding (default 50)
  #   mnn_k: Number of nearest neighbors for MNN matching (default 20)
  # Returns:
  #   Merged Seurat object with "integrated.mnn" reduction
  # Note: fastMNN keeps LogNormalize + FindVariableFeatures but skips ScaleData
  #       and RunPCA upstream — fastMNN runs its own multiBatchPCA internally
  
  n_datasets <- length(seurat_list)
  message(paste("Starting fastMNN with", n_datasets, "datasets"))
  
  if (!requireNamespace("batchelor", quietly = TRUE)) {
    stop("Package 'batchelor' is required for fastMNN integration")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required for fastMNN integration")
  }
  
  gc()
  
  message("Step 1/4: Merging datasets")
  merged_object <- performSimpleMerge(seurat_list)
  
  rm(seurat_list)
  gc()
  
  message("Step 2/4: Log-normalizing merged object")
  merged_object <- NormalizeData(merged_object, verbose = FALSE)
  
  n_cells <- ncol(merged_object)
  n_features <- if (n_cells > 150000) {
    2000
  } else if (n_cells > 100000) {
    2500
  } else {
    3000
  }
  
  message(paste("Step 3/4: Selecting", n_features, "variable features"))
  merged_object <- FindVariableFeatures(
    merged_object,
    selection.method = "vst",
    nfeatures = n_features,
    verbose = FALSE
  )
  
  # Validate that 'dataset' metadata exists and has at least 2 levels
  if (!"dataset" %in% colnames(merged_object@meta.data)) {
    stop("fastMNN requires a 'dataset' metadata column for batch assignment")
  }
  batch_vec <- as.character(merged_object$dataset)
  unique_batches <- unique(batch_vec[!is.na(batch_vec)])
  if (length(unique_batches) < 2) {
    stop(paste(
      "ERROR: fastMNN requires at least 2 batches in 'dataset' column.",
      "\nFound only", length(unique_batches), "level:", paste(unique_batches, collapse = ", "),
      "\n\nSuggestions:",
      "\n  1. Check that your datasets are actually different",
      "\n  2. Try 'Simple Merge' if no batch correction is needed"
    ))
  }
  
  hvg <- VariableFeatures(merged_object)
  log_mat <- GetAssayData(merged_object, assay = "RNA", layer = "data")
  
  # Subset to HVGs (fastMNN's subset.row argument supports gene names)
  log_mat <- log_mat[hvg, , drop = FALSE]
  
  message(paste("Step 4/4: Running fastMNN (d =", mnn_dims, ", k =", mnn_k, 
                ", batches =", length(unique_batches), ")"))
  
  # Seed for reproducibility — fastMNN uses IRLBA for SVD which has random init
  set.seed(42)
  
  mnn_out <- batchelor::fastMNN(
    log_mat,
    batch = batch_vec,
    d = mnn_dims,
    k = mnn_k,
    BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE),
    BPPARAM = BiocParallel::SerialParam()
  )
  
  # Extract corrected embedding and inject as a Seurat reduction
  mnn_embeddings <- SingleCellExperiment::reducedDim(mnn_out, "corrected")
  rownames(mnn_embeddings) <- colnames(merged_object)
  colnames(mnn_embeddings) <- paste0("MNN_", seq_len(ncol(mnn_embeddings)))
  
  merged_object[["integrated.mnn"]] <- CreateDimReducObject(
    embeddings = mnn_embeddings,
    key = "MNN_",
    assay = DefaultAssay(merged_object)
  )
  
  rm(mnn_out, log_mat)
  gc()
  
  merged_object$integration_method <- "fastmnn"
  DefaultAssay(merged_object) <- "RNA"
  
  message(paste("✓ fastMNN completed successfully:",
                ncol(merged_object), "cells integrated,",
                ncol(mnn_embeddings), "MNN dimensions"))
  
  return(merged_object)
}
