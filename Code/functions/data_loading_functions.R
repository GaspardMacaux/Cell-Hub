# data_loading_functions_newarch.R

############################## Utility Functions ##############################

getDatasetFileName <- function(file_inputs = NULL, default_name = "Dataset") {
  # Extract dataset file name from various input sources
  # Args:
  #   file_inputs: List of file inputs to check (e.g., list(input$file, input$load_seurat_file))
  #   default_name: Default name if no file found
  
  tryCatch({
    file_name <- ""
    
    # If specific inputs provided, check them
    if (!is.null(file_inputs)) {
      for (file_input in file_inputs) {
        if (!is.null(file_input) && 
            !is.null(file_input$name) && 
            nchar(file_input$name) > 0) {
          file_name <- tools::file_path_sans_ext(basename(file_input$name))
          break  # Use first valid file found
        }
      }
    }
    
    # Return file name or default
    if (nchar(file_name) == 0) {
      return(default_name)
    } else {
      return(file_name)
    }
    
  }, error = function(e) {
    message(paste("Error extracting file name:", e$message))
    return(default_name)
  })
}

getClusters <- function(seurat_object) {
  # Get clusters from Seurat object
  # Args:
  #   seurat_object: Seurat object
  # Returns:
  #   Vector of cluster levels
  
  if (is.null(seurat_object)) {
    return(NULL)
  }
  
  tryCatch({
    return(levels(Idents(seurat_object)))
  }, error = function(e) {
    message(paste("Error getting clusters:", e$message))
    return(NULL)
  })
}

find10XFolder <- function(base_dir) {
  message(paste("Searching for 10X files in:", base_dir))
  folders <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
  for (folder in folders) {
    compressed_files <- list.files(
      folder,
      pattern = "barcodes\\.tsv\\.gz$|features\\.tsv\\.gz$|matrix\\.mtx\\.gz$|genes\\.tsv\\.gz$",
      full.names = TRUE
    )
    uncompressed_files <- list.files(
      folder,
      pattern = "barcodes\\.tsv$|features\\.tsv$|matrix\\.mtx$|genes\\.tsv$",
      full.names = TRUE
    )
    if (length(compressed_files) >= 3 || length(uncompressed_files) >= 3) {
      message(paste("Found 10X files in:", folder))
      return(folder)
    }
  }
  return(NULL)
}

compressFiles <- function(files) {
  gz_files <- sapply(files, function(file) {
    gz_file <- paste0(file, ".gz")
    if (!file.exists(gz_file)) {
      message(paste("Compressing file:", file))
      R.utils::gzip(file, destname = gz_file, overwrite = TRUE)
      message(paste("Compressed", file, "to", gz_file))
    } else {
      message(gz_file, "already exists.")
    }
    return(gz_file)
  })
  return(gz_files)
}

extractAndProcessZip <- function(zip_path, dataset_name = NULL) {
  # Extract and process ZIP file with handling for nested ZIP files
  # Args:
  #   zip_path: Path to ZIP file
  #   dataset_name: Name for dataset (optional)
  # Returns:
  #   Path to folder containing 10X files or NULL
  
  # Create unique identifier for this dataset
  if (is.null(dataset_name)) {
    dataset_name <- format(Sys.time(), "%Y%m%d_%H%M%S")
  }
  
  # Create specific temporary folder for this dataset
  temp_base_dir <- file.path(tempdir(), "single_dataset")
  dir.create(temp_base_dir, showWarnings = FALSE, recursive = TRUE)
  
  dataset_dir <- file.path(temp_base_dir, paste0("dataset_", dataset_name))
  dir.create(dataset_dir, showWarnings = FALSE, recursive = TRUE)
  
  message(paste("Extracting zip file to unique directory:", dataset_dir))
  
  # Extract ZIP file
  tryCatch({
    unzip(zip_path, exdir = dataset_dir)
    
    # Search for nested ZIP files
    nested_zips <- list.files(dataset_dir, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
    
    while (length(nested_zips) > 0) {
      for (nested_zip in nested_zips) {
        message(paste("Extracting nested zip file:", nested_zip))
        nested_dir <- file.path(dirname(nested_zip), tools::file_path_sans_ext(basename(nested_zip)))
        dir.create(nested_dir, showWarnings = FALSE)
        unzip(nested_zip, exdir = nested_dir)
        file.remove(nested_zip)  # Remove ZIP after extraction
      }
      # Check if there are remaining nested ZIP files
      nested_zips <- list.files(dataset_dir, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
    }
    
    # Search for folder containing 10X files
    valid_folder <- find10XFolder(dataset_dir)
    
    if (is.null(valid_folder)) {
      stop("No complete set of 10X files found in the uploaded ZIP")
    }
    
    return(valid_folder)
    
  }, error = function(e) {
    message(paste("Error extracting zip:", e$message))
    return(NULL)
  })
}

validateAndCompleteSeurat <- function(seurat_object, log_function = message) {
  # Validate and complete preprocessing steps for a Seurat object
  # Compatible with both Seurat v4 and v5
  
  log_function("Validating Seurat object preprocessing...", "info")
  
  # Determine Seurat version
  seurat_v5 <- packageVersion("Seurat") >= "5.0.0"
  
  # Check 1: Normalized data layer
  has_data <- FALSE
  if (seurat_v5) {
    # Seurat v5: Check layers
    has_data <- "data" %in% names(seurat_object[["RNA"]]@layers)
  } else {
    # Seurat v4: Check slot
    has_data <- length(seurat_object[["RNA"]]@data) > 0
  }
  
  if (!has_data) {
    log_function("Missing normalized data - normalizing...", "warning")
    seurat_object <- NormalizeData(seurat_object, verbose = FALSE)
    log_function("Normalization completed", "success")
  } else {
    log_function("Normalized data found", "success")
  }
  
  # Check 2: Variable features
  has_variable_features <- length(VariableFeatures(seurat_object)) > 0
  
  if (!has_variable_features) {
    log_function("Missing variable features - calculating...", "warning")
    seurat_object <- FindVariableFeatures(seurat_object,
                                          selection.method = "vst",
                                          nfeatures = 4000,
                                          verbose = FALSE)
    log_function(paste("Found", length(VariableFeatures(seurat_object)), "variable features"), "success")
  } else {
    log_function(paste("Found", length(VariableFeatures(seurat_object)), "variable features"), "success")
  }
  
  # Check 3: Scaled data - FIXED FOR SEURAT V5
  has_scaled <- FALSE
  if (seurat_v5) {
    # Seurat v5: Check if scale.data layer exists
    tryCatch({
      has_scaled <- "scale.data" %in% names(seurat_object[["RNA"]]@layers)
    }, error = function(e) {
      # If error, assume no scaled data
      has_scaled <<- FALSE
    })
  } else {
    # Seurat v4: Check slot
    has_scaled <- length(seurat_object[["RNA"]]@scale.data) > 0
  }
  
  if (!has_scaled) {
    log_function("Missing scaled data - scaling...", "warning")
    seurat_object <- ScaleData(seurat_object,
                               features = VariableFeatures(seurat_object),
                               verbose = FALSE)
    log_function("Scaling completed", "success")
  } else {
    log_function("Scaled data found", "success")
  }
  
  # Check 4: PCA
  has_pca <- "pca" %in% names(seurat_object@reductions)
  if (!has_pca) {
    log_function("Missing PCA - running...", "warning")
    seurat_object <- runPCA_reproducible(seurat_object,
                                         features = VariableFeatures(seurat_object),
                                         verbose = FALSE)
    log_function("PCA completed", "success")
  } else {
    log_function("PCA found", "success")
  }
  
  log_function("Object validation and completion finished", "success")
  
  return(seurat_object)
}

############################## Main Loading Functions ##############################

loadRawData <- function(file_path, dataset_type, species, log_function = message) {
  # Load raw data from ZIP (10X), H5 (10X HDF5), or H5AD (AnnData) format
  # Args:
  #   file_path: Path to ZIP, H5, or H5AD file
  #   dataset_type: Type of dataset ("snRNA" or "multiome")
  #   species: Species for mitochondrial pattern
  #   log_function: Function for logging (default: message)
  # Returns:
  #   Seurat object
  file_extension <- tolower(tools::file_ext(file_path))
  mt_pattern <- switch(species,
                       "human" = "^MT-",
                       "mouse" = "^mt-|^Mt-",
                       "rat"   = "^Mt-",
                       "^MT-")
  log_function(paste("Mitochondrial pattern:", mt_pattern), "info")
  seuratObj <- NULL
  if (file_extension == "zip") {
    log_function("Processing ZIP file...", "info")
    seuratObj <- loadFromZip(file_path, dataset_type, mt_pattern, log_function)
  } else if (file_extension %in% c("h5", "hdf5")) {
    log_function("Processing H5 file...", "info")
    seuratObj <- loadFromH5(file_path, dataset_type, mt_pattern, log_function)
  } else if (file_extension == "h5ad") {
    log_function("Processing H5AD (AnnData) file...", "info")
    seuratObj <- loadFromH5AD(h5ad_path = file_path,
                              dataset_type = dataset_type,
                              mt_pattern = mt_pattern,
                              log_function = log_function)
  } else {
    stop("Unsupported file format. Provide ZIP, H5, HDF5, or H5AD file.")
  }
  if (is.null(seuratObj)) {
    stop("Failed to create Seurat object from raw data")
  }
  seuratObj@meta.data$dataset_type <- dataset_type
  seuratObj@meta.data$species <- species
  seuratObj@meta.data$source_format <- file_extension
  return(seuratObj)
}



loadFromZip <- function(zip_path, dataset_type, mt_pattern, log_function = message) {
  # Load data from ZIP file containing 10X files
  
  # Extract ZIP and find 10X files
  log_function("Extracting ZIP file...", "info")
  valid_folder <- extractAndProcessZip(zip_path)
  
  if (is.null(valid_folder)) {
    stop("Could not extract and process ZIP file")
  }
  
  # Load 10X data
  log_function("Reading 10X data...", "info")
  data_10x <- Read10X(valid_folder)
  
  # Detect data characteristics
  n_barcodes <- ifelse(is.list(data_10x), ncol(data_10x[[1]]), ncol(data_10x))
  n_features <- ifelse(is.list(data_10x), nrow(data_10x[[1]]), nrow(data_10x))
  
  log_function(paste("Matrix dimensions:", format(n_features, big.mark = ","), 
                     "genes x", format(n_barcodes, big.mark = ","), "barcodes"), "info")
  
  # Create Seurat object based on dataset type
  log_function("Creating Seurat object...", "info")
  
  if (dataset_type == "snRNA") {
    seuratObj <- CreateSeuratObject(
      counts = data_10x,
      project = "RawData"
    )
    
  } else if (dataset_type == "multiome") {
    if (is.list(data_10x) && !is.matrix(data_10x)) {
      # Multimodal data
      log_function("Multiome data detected - identifying RNA assay", "info")
      
      rna_assay <- NULL
      if ("Gene Expression" %in% names(data_10x)) {
        rna_assay <- "Gene Expression"
      } else {
        rna_candidates <- grep("RNA|Gene|Expression", names(data_10x), ignore.case = TRUE, value = TRUE)
        if (length(rna_candidates) > 0) {
          rna_assay <- rna_candidates[1]
        }
      }
      
      if (is.null(rna_assay)) {
        stop("Could not identify RNA assay in multimodal data")
      }
      
      log_function(paste("Using RNA assay:", rna_assay), "info")
      
      seuratObj <- CreateSeuratObject(
        counts = data_10x[[rna_assay]],
        project = "RawData"
      )
      
      # Add other assays
      log_function("Adding additional assays...", "info")
      for (assay_name in names(data_10x)) {
        if (assay_name != rna_assay) {
          log_function(paste("  - Adding assay:", assay_name), "debug")
          seuratObj[[assay_name]] <- CreateAssayObject(counts = data_10x[[assay_name]])
        }
      }
      
    } else {
      # Single assay
      seuratObj <- CreateSeuratObject(
        counts = data_10x,
        project = "RawData"
      )
    }
  }
  
  # Add mitochondrial percentage
  log_function("Calculating mitochondrial percentage...", "info")
  seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = mt_pattern)
  
  # Log results
  log_function(paste("Cells loaded:", format(ncol(seuratObj), big.mark = ",")), "success")
  log_function(paste("Genes loaded:", format(nrow(seuratObj), big.mark = ",")), "success")
  
  return(seuratObj)
}

loadFromH5 <- function(h5_path, dataset_type, mt_pattern, log_function = message) {
  # Load data from H5 file
  
  log_function("Reading H5 file structure...", "info")
  h5_data <- Read10X_h5(h5_path)
  
  # Detect file characteristics
  n_barcodes <- ifelse(is.list(h5_data), ncol(h5_data[[1]]), ncol(h5_data))
  n_features <- ifelse(is.list(h5_data), nrow(h5_data[[1]]), nrow(h5_data))
  
  log_function(paste("Matrix dimensions:", format(n_features, big.mark = ","), 
                     "genes x", format(n_barcodes, big.mark = ","), "barcodes"), "info")
  
  # Decide on filtering strategy
  is_raw_file <- n_barcodes > 500000
  
  if (is_raw_file) {
    log_function("Detected RAW file - applying initial filtering", "warning")
    log_function("Filtering: min.cells=3, min.features=200", "info")
    min_cells_filter <- 3
    min_features_filter <- 200
  } else {
    log_function("Detected FILTERED file - loading without initial filtering", "info")
    min_cells_filter <- 0
    min_features_filter <- 0
  }
  
  # Create Seurat object based on dataset type
  log_function("Creating Seurat object...", "info")
  
  if (dataset_type == "snRNA") {
    seuratObj <- CreateSeuratObject(
      counts = h5_data,
      project = "RawData",
      min.cells = min_cells_filter,
      min.features = min_features_filter
    )
    
  } else if (dataset_type == "multiome") {
    if (is.list(h5_data) && !is.matrix(h5_data)) {
      # Multimodal H5
      log_function("Multiome data detected - identifying RNA assay", "info")
      
      rna_assay <- NULL
      if ("Gene Expression" %in% names(h5_data)) {
        rna_assay <- "Gene Expression"
      } else {
        rna_candidates <- grep("RNA|Gene|Expression", names(h5_data), ignore.case = TRUE, value = TRUE)
        if (length(rna_candidates) > 0) {
          rna_assay <- rna_candidates[1]
        }
      }
      
      if (is.null(rna_assay)) {
        stop("Could not identify RNA assay in H5 file")
      }
      
      log_function(paste("Using RNA assay:", rna_assay), "info")
      rna_data <- h5_data[[rna_assay]]
      
      seuratObj <- CreateSeuratObject(
        counts = rna_data,
        project = "RawData",
        min.cells = min_cells_filter,
        min.features = min_features_filter
      )
      
      # Add other assays
      log_function("Adding additional assays...", "info")
      for (assay_name in names(h5_data)) {
        if (assay_name != rna_assay) {
          log_function(paste("  - Adding assay:", assay_name), "debug")
          seuratObj[[assay_name]] <- CreateAssayObject(counts = h5_data[[assay_name]])
        }
      }
      
    } else {
      # Single assay H5
      seuratObj <- CreateSeuratObject(
        counts = h5_data,
        project = "RawData",
        min.cells = min_cells_filter,
        min.features = min_features_filter
      )
    }
  }
  
  # Add mitochondrial percentage
  log_function("Calculating mitochondrial percentage...", "info")
  seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = mt_pattern)
  
  # Log filtering results
  cells_kept <- ncol(seuratObj)
  genes_kept <- nrow(seuratObj)
  pct_cells_kept <- round(100 * cells_kept / n_barcodes, 2)
  
  log_function(paste("Cells retained:", format(cells_kept, big.mark = ","), 
                     paste0("(", pct_cells_kept, "%)")), "success")
  log_function(paste("Genes retained:", format(genes_kept, big.mark = ",")), "success")
  
  if (is_raw_file && pct_cells_kept < 5) {
    log_function(paste("Only", pct_cells_kept, "% of barcodes retained - this is normal for RAW files"), "info")
    log_function("Empty droplets successfully removed", "success")
  }
  
  return(seuratObj)
}


# Load and restore colors after loading Seurat object for single dataset
loadSeuratObject <- function(rds_path, add_dataset_column = FALSE, dataset_name = NULL,
                             clean_before = FALSE, module_type = "single",
                             file_format = NULL,
                             log_function = message) {
  # Load processed Seurat object from RDS or H5AD file
  # Handles: Seurat v4 objects (auto-converts to v5), combined Seurat+CellChat,
  #          H5AD via schard, multi-layer v5 objects (auto-joins)
  # Args:
  #   rds_path: Path to RDS or H5AD file (Shiny temp path — may lack extension)
  #   add_dataset_column: Whether to add/verify dataset column in metadata
  #   dataset_name: Name for dataset column (if NULL, extracted from filename)
  #   clean_before: Whether to clean workspace before loading
  #   module_type: Module type for cleaning ("single", "multiple", etc.)
  #   file_format: File format override ("rds", "h5ad"). If NULL, inferred from rds_path extension.
  #   log_function: Function for logging (default: message)
  # Returns:
  #   Seurat v5 object with Assay5 class assays, joined layers
  if (clean_before) {
    log_function("Cleaning workspace...", "info")
    cleanWorkspace(module_type)
  }
  # Determine file format: use explicit override if provided, else infer from path
  if (is.null(file_format) || !nzchar(file_format)) {
    file_format <- tolower(tools::file_ext(rds_path))
  }
  if (file_format == "h5ad") {
    log_function("Detected H5AD (AnnData) file - loading via schard", "info")
    inspectH5ADFile(rds_path, log_function)
    seuratObj <- loadFromH5AD(h5ad_path = rds_path,
                              dataset_type = "snRNA",
                              mt_pattern = "^MT-|^mt-|^Mt-",
                              log_function = log_function)
    log_function("✓ H5AD loaded as Seurat object", "success")
  } else {
    log_function("Loading Seurat object from RDS...", "info")
    loaded_data <- readRDS(rds_path)
    if (is.list(loaded_data) && !is.null(loaded_data$seurat) && !is.null(loaded_data$cellchat)) {
      log_function("✓ Detected combined Seurat + CellChat file", "success")
      log_function(paste("  Found", length(loaded_data$cellchat), "CellChat object(s)"), "info")
      log_function("  Extracting Seurat object only...", "info")
      seuratObj <- loaded_data$seurat
      if (!is.null(loaded_data$cellchat) && length(loaded_data$cellchat) > 0) {
        seuratObj@misc$cellchat_available <- TRUE
        seuratObj@misc$cellchat_names <- names(loaded_data$cellchat)
        log_function(paste("  Note: CellChat objects available:",
                           paste(names(loaded_data$cellchat), collapse = ", ")), "info")
        log_function("  (Use CellChat module to load full analysis)", "info")
      }
    } else {
      seuratObj <- loaded_data
    }
  }
  if (!inherits(seuratObj, "Seurat")) {
    log_function("File is not a valid Seurat object", "error")
    stop("File is not a valid Seurat object")
  }
  log_function("✓ Valid Seurat object detected", "success")
  # Convert Seurat v4 Assay objects to v5 Assay5 if needed
  # UpdateSeuratObject() updates the outer Seurat structure but does NOT convert
  # Assay → Assay5. We must do that explicitly with as(assay, "Assay5").
  v4_assays_found <- FALSE
  for (assay_name in names(seuratObj@assays)) {
    if (inherits(seuratObj[[assay_name]], "Assay") &&
        !inherits(seuratObj[[assay_name]], "Assay5")) {
      v4_assays_found <- TRUE
      break
    }
  }
  if (v4_assays_found) {
    log_function("Detected Seurat v4 assay format — converting to v5...", "info")
    # First update the outer object structure (slots, keys, dimreduc names)
    tryCatch({
      seuratObj <- UpdateSeuratObject(seuratObj)
      log_function("✓ Object structure updated", "success")
    }, error = function(e) {
      log_function(paste("UpdateSeuratObject warning:", conditionMessage(e)[1]), "warning")
    })
    # Then force-convert each Assay → Assay5
    for (assay_name in names(seuratObj@assays)) {
      if (inherits(seuratObj[[assay_name]], "Assay") &&
          !inherits(seuratObj[[assay_name]], "Assay5")) {
        tryCatch({
          log_function(paste("  Converting assay:", assay_name), "info")
          seuratObj[[assay_name]] <- as(seuratObj[[assay_name]], "Assay5")
          log_function(paste("  ✓", assay_name, "converted to Assay5"), "success")
        }, error = function(e) {
          log_function(paste("  ✗ Failed to convert", assay_name, ":",
                             conditionMessage(e)[1]), "warning")
        })
      }
    }
  }
  # Ensure a usable default assay
  if ("RNA" %in% names(seuratObj@assays)) {
    DefaultAssay(seuratObj) <- "RNA"
  } else {
    first_assay <- names(seuratObj@assays)[1]
    DefaultAssay(seuratObj) <- first_assay
    log_function(paste("No RNA assay found — using", first_assay, "as default"), "info")
  }
  # Check and preserve cluster colors
  if (!is.null(seuratObj@misc$cluster_colors)) {
    log_function("Found existing cluster colors in Seurat object", "info")
    cluster_colors <- seuratObj@misc$cluster_colors
    current_clusters <- as.character(unique(Idents(seuratObj)))
    stored_clusters <- names(cluster_colors)
    if (all(current_clusters %in% stored_clusters)) {
      log_function("Cluster colors are compatible with current clusters", "success")
    } else {
      log_function("Cluster colors don't match current clusters, will generate new ones", "warning")
      seuratObj@misc$cluster_colors <- NULL
    }
  }
  # Join multi-fragment layers (Seurat v5 only, Assay5 only)
  if (packageVersion("Seurat") >= "5.0.0") {
    log_function("Checking for layers to join...", "info")
    for (assay_name in names(seuratObj@assays)) {
      assay_obj <- seuratObj[[assay_name]]
      if (!inherits(assay_obj, "Assay5")) {
        log_function(paste("Skipping", assay_name, "— still v4 Assay class (conversion may have failed)"), "info")
        next
      }
      tryCatch({
        n_layers <- length(names(assay_obj@layers))
        if (n_layers > 1) {
          log_function(paste("Found", n_layers, "layers in", assay_name, ", joining..."), "info")
          seuratObj[[assay_name]] <- JoinLayers(seuratObj[[assay_name]])
          log_function(paste("✓", assay_name, "layers joined"), "success")
        } else {
          log_function(paste(assay_name, ":", n_layers, "layer(s), no join needed"), "info")
        }
      }, error = function(e) {
        log_function(paste("Could not join layers for", assay_name, ":",
                           conditionMessage(e)[1]), "warning")
      })
    }
  }
  # Handle dataset column for multiple datasets workflow
  if (add_dataset_column) {
    if (!"dataset" %in% colnames(seuratObj@meta.data)) {
      log_function("Adding 'dataset' column to metadata", "info")
      if (is.null(dataset_name)) {
        dataset_name <- tools::file_path_sans_ext(basename(rds_path))
        if (is.null(dataset_name) || dataset_name == "") {
          dataset_name <- "Dataset"
        }
      }
      seuratObj$dataset <- dataset_name
      log_function(paste("Set dataset name to:", dataset_name), "success")
    } else {
      datasets <- unique(seuratObj$dataset)
      log_function(paste("Found existing datasets:", paste(datasets, collapse = ", ")), "info")
    }
  }
  log_function(paste("Loaded Seurat object:", format(ncol(seuratObj), big.mark = ","),
                     "cells,", format(nrow(seuratObj), big.mark = ","), "genes"), "success")
  log_function(paste("Available assays:", paste(names(seuratObj@assays), collapse = ", ")), "info")
  if (length(seuratObj@reductions) > 0) {
    log_function(paste("Available reductions:", paste(names(seuratObj@reductions), collapse = ", ")), "info")
  }
  return(seuratObj)
}

generateInitialPlot <- function(seurat_obj, remove_axes = FALSE, remove_legend = FALSE) {
  # Generate initial clustering plot with color restoration
  # Args:
  #   seurat_obj: Seurat object
  #   remove_axes: Remove axes from plot
  #   remove_legend: Remove legend from plot
  # Returns:
  #   List with plot and updated seurat object
  
  if (is.null(seurat_obj) || !"umap" %in% names(seurat_obj@reductions)) {
    return(NULL)
  }
  
  # Handle cluster colors
  if (!is.null(seurat_obj@misc$cluster_colors)) {
    message("Restoring cluster colors from Seurat object.")
    colors <- seurat_obj@misc$cluster_colors
    # Create plot with custom colors
    plot <- DimPlot(seurat_obj, group.by = "ident", label = FALSE) +
      ggtitle("") +
      scale_color_manual(values = colors)
  } else {
    message("No cluster colors found. Using Seurat default colors.")
    
    plot <- DimPlot(seurat_obj, group.by = "ident", label = FALSE) + 
      ggtitle("")
  }
  
  if (remove_axes) { plot <- plot + NoAxes() }
  if (remove_legend) { plot <- plot + NoLegend() }
  
  return(list(plot = plot, seurat_obj = seurat_obj))
}

############################## Logging Functions ##############################

createLogSystem <- function() {
  # Create a reactive log system for a module
  # Returns:
  #   List with log_value (reactiveVal) and add_log (function)
  
  log_value <- reactiveVal("")
  
  add_log <- function(message, type = "info") {
    timestamp <- format(Sys.time(), "%H:%M:%S")
    
    # Add color prefix based on type
    prefix <- switch(type,
                     "info" = "[INFO]",
                     "warning" = "[⚠️ WARN]",
                     "error" = "[❌ ERROR]",
                     "success" = "[✓ SUCCESS]",
                     "debug" = "[🔍 DEBUG]",
                     "[INFO]"
    )
    
    new_log <- paste0(timestamp, " ", prefix, " ", message)
    current_logs <- log_value()
    
    if (current_logs == "") {
      log_value(new_log)
    } else {
      log_value(paste(current_logs, new_log, sep = "\n"))
    }
  }
  
  clear_logs <- function() {
    log_value("")
  }
  
  return(list(
    log_value = log_value,
    add_log = add_log,
    clear_logs = clear_logs
  ))
}

logFileInfo <- function(file_name, file_path, log_function) {
  # Log basic file information
  # Args:
  #   file_name: Name of the uploaded file
  #   file_path: Path to the file
  #   log_function: Function to call for logging
  
  file_extension <- tolower(tools::file_ext(file_name))
  file_size_mb <- round(file.size(file_path) / 1024^2, 1)
  
  log_function(paste("File detected:", file_name), "info")
  log_function(paste("File extension:", file_extension), "info")
  log_function(paste("File size:", file_size_mb, "MB"), "info")
  
  return(file_extension)
}

inspectH5File <- function(h5_path, log_function) {
  # Inspect H5 file and log diagnostics without loading full data
  # Args:
  #   h5_path: Path to H5 file
  #   log_function: Function to call for logging
  # Returns:
  #   List with n_barcodes, assays, is_raw
  
  log_function("Inspecting H5 file structure...", "info")
  
  h5_info <- tryCatch({
    h5_conn <- hdf5r::H5File$new(h5_path, mode = "r")
    
    # Get barcode count
    n_barcodes <- NA
    if (h5_conn$exists("matrix/barcodes")) {
      n_barcodes <- h5_conn[["matrix/barcodes"]]$dims
    }
    
    # Get available assays
    assays <- character(0)
    if (h5_conn$exists("matrix/features")) {
      features <- h5_conn[["matrix/features"]]
      if ("feature_type" %in% names(features)) {
        assays <- unique(as.character(features[["feature_type"]][]))
      }
    }
    
    h5_conn$close()
    
    list(n_barcodes = n_barcodes, assays = assays, success = TRUE)
  }, error = function(e) {
    log_function(paste("Could not inspect H5 structure:", e$message), "warning")
    list(n_barcodes = NA, assays = character(0), success = FALSE)
  })
  
  # Log diagnostics
  if (!is.na(h5_info$n_barcodes)) {
    log_function(paste("Barcodes detected:", format(h5_info$n_barcodes, big.mark = ",")), "info")
    
    is_raw <- h5_info$n_barcodes > 500000
    h5_info$is_raw <- is_raw
    
    if (is_raw) {
      log_function("⚠️ WARNING: File contains >500k barcodes", "warning")
      log_function("This looks like a RAW file (raw_feature_bc_matrix.h5)", "warning")
      log_function("Recommendation: Use filtered_feature_bc_matrix.h5 for faster loading", "warning")
      log_function("Initial filtering will be applied to remove empty droplets", "info")
    } else if (h5_info$n_barcodes < 100000) {
      log_function("✓ File appears to be pre-filtered (filtered_feature_bc_matrix.h5)", "success")
    }
  } else {
    h5_info$is_raw <- FALSE
  }
  
  if (length(h5_info$assays) > 0) {
    log_function(paste("Available assays:", paste(h5_info$assays, collapse = ", ")), "info")
  }
  
  return(h5_info)
}

logSeuratObjectStats <- function(seurat_obj, log_function) {
  # Log comprehensive statistics about a Seurat object
  # Args:
  #   seurat_obj: Seurat object to analyze
  #   log_function: Function to call for logging
  
  log_function("=== LOADING COMPLETED ===", "success")
  log_function(paste("Final cells:", format(ncol(seurat_obj), big.mark = ",")), "success")
  log_function(paste("Final genes:", format(nrow(seurat_obj), big.mark = ",")), "success")
  log_function(paste("Available assays:", paste(names(seurat_obj@assays), collapse = ", ")), "info")
  
  # QC metrics
  if ("nFeature_RNA" %in% colnames(seurat_obj@meta.data)) {
    log_function(paste("Median nFeature_RNA:", round(median(seurat_obj$nFeature_RNA))), "info")
  }
  
  if ("nCount_RNA" %in% colnames(seurat_obj@meta.data)) {
    log_function(paste("Median nCount_RNA:", round(median(seurat_obj$nCount_RNA))), "info")
  }
  
  if ("percent.mt" %in% colnames(seurat_obj@meta.data)) {
    log_function(paste("Median percent.mt:", round(median(seurat_obj$percent.mt), 2), "%"), "info")
  }
}

renderLogsUI <- function(output_id, title = "Data Loading Logs") {
  # Create a standardized logs UI panel
  # Args:
  #   output_id: ID for the verbatimTextOutput
  #   title: Title for the logs panel
  # Returns:
  #   Shiny tags div with logs display
  
  tags$div(
    style = "margin-top: 20px; background-color: #f8f9fa; padding: 20px; 
           border-radius: 8px; border: 2px solid #dee2e6;",
    tags$h4(
      icon("terminal"), paste(" ", title), 
      style = "color: #495057; margin-top: 0; margin-bottom: 15px; font-weight: bold;"
    ),
    tags$pre(
      style = "background-color: #1e272e; color: #00ff00; padding: 15px; 
             border-radius: 6px; font-family: 'Courier New', monospace; 
             font-size: 11px; max-height: 400px; overflow-y: auto; 
             margin: 0; border-left: 4px solid #667eea; white-space: pre-wrap;",
      verbatimTextOutput(output_id)
    ),
    tags$div(
      style = "margin-top: 10px; padding: 10px; background-color: #fff3cd; 
             border-radius: 5px; border-left: 3px solid #ffc107;",
      icon("info-circle"), tags$strong(" Tip: "),
      "These logs help diagnose loading issues. If you see warnings about RAW files, 
      consider using filtered_feature_bc_matrix.h5 for faster loading."
    )
  )
}


inspectH5ADFile <- function(h5ad_path, log_function = message) {
  # Inspect H5AD (AnnData) file structure and log diagnostics without loading
  # Args:
  #   h5ad_path: Path to .h5ad file
  #   log_function: Logging function
  # Returns:
  #   List with structural metadata
  log_function("Inspecting H5AD file structure...", "info")
  result <- list(
    n_cells = NA, n_genes = NA,
    has_raw = FALSE, has_layers_counts = FALSE,
    obsm_keys = character(0), obs_columns = character(0),
    encoding_version = NA_character_, success = FALSE
  )
  tryCatch({
    h5 <- hdf5r::H5File$new(h5ad_path, mode = "r")
    on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
    # Encoding version (anndata >= 0.8)
    if (h5$attr_exists("encoding-version")) {
      result$encoding_version <- as.character(hdf5r::h5attr(h5, "encoding-version"))
    }
    # X dimensions: dataset (dense) or group (sparse with shape attr)
    if (h5$exists("X")) {
      x_obj <- h5[["X"]]
      if (inherits(x_obj, "H5D")) {
        dims <- x_obj$dims
        result$n_cells <- dims[1]; result$n_genes <- dims[2]
      } else if (inherits(x_obj, "H5Group") && x_obj$attr_exists("shape")) {
        shape <- hdf5r::h5attr(x_obj, "shape")
        result$n_cells <- shape[1]; result$n_genes <- shape[2]
      }
    }
    # Raw counts presence
    result$has_raw <- h5$exists("raw") || h5$exists("raw.X")
    result$has_layers_counts <- h5$exists("layers/counts") ||
      h5$exists("layers/raw_counts") ||
      h5$exists("layers/raw")
    # Embeddings
    if (h5$exists("obsm")) result$obsm_keys <- names(h5[["obsm"]])
    # Metadata columns
    if (h5$exists("obs")) result$obs_columns <- names(h5[["obs"]])
    result$success <- TRUE
  }, error = function(e) {
    log_function(paste("H5AD inspection failed:", conditionMessage(e)[1]), "warning")
  })
  # Log findings
  if (result$success) {
    if (!is.na(result$encoding_version) && nzchar(result$encoding_version)) {
      log_function(paste("AnnData encoding version:", result$encoding_version), "info")
    }
    if (!is.na(result$n_cells)) {
      log_function(paste("Dimensions:", format(result$n_cells, big.mark = ","), "cells x",
                         format(result$n_genes, big.mark = ","), "genes"), "info")
    }
    log_function(paste("raw slot present:", result$has_raw), "info")
    log_function(paste("layers/counts present:", result$has_layers_counts), "info")
    if (length(result$obsm_keys) > 0) {
      log_function(paste("Embeddings (obsm):", paste(result$obsm_keys, collapse = ", ")), "info")
    }
    if (length(result$obs_columns) > 0) {
      n_show <- min(10, length(result$obs_columns))
      cols_show <- paste(result$obs_columns[seq_len(n_show)], collapse = ", ")
      if (length(result$obs_columns) > n_show) {
        cols_show <- paste0(cols_show, ", ... (+", length(result$obs_columns) - n_show, " more)")
      }
      log_function(paste("Metadata columns:", cols_show), "info")
    }
    if (!result$has_raw && !result$has_layers_counts) {
      log_function("⚠️ No raw counts detected - X will be used as counts", "warning")
      log_function("If X is normalized, downstream NormalizeData would corrupt the data", "warning")
    }
  }
  return(result)
}

cleanH5ADFactors <- function(seurat_obj, log_function = message) {
  # Apply droplevels() on all factor columns in meta.data (CellChat compatibility, rule #5)
  # Args:
  #   seurat_obj: Seurat object
  #   log_function: Logging function
  # Returns:
  #   Seurat object with cleaned factor levels
  factor_cols <- vapply(seurat_obj@meta.data, is.factor, logical(1))
  if (any(factor_cols)) {
    n_factors <- sum(factor_cols)
    log_function(paste("Cleaning", n_factors, "factor column(s) with droplevels()"), "info")
    for (col_name in names(seurat_obj@meta.data)[factor_cols]) {
      seurat_obj@meta.data[[col_name]] <- droplevels(seurat_obj@meta.data[[col_name]])
    }
  }
  return(seurat_obj)
}

loadFromH5AD <- function(h5ad_path, dataset_type = "snRNA",
                         mt_pattern = "^MT-|^mt-|^Mt-",
                         log_function = message) {
  # Load Seurat object from H5AD (AnnData) file via schard
  # Handles Ensembl ID → gene symbol remapping before Seurat object creation
  # Args:
  #   h5ad_path: Path to .h5ad file
  #   dataset_type: "snRNA" or "multiome" (informational only)
  #   mt_pattern: Mitochondrial gene regex (only used if percent.mt absent in obs)
  #   log_function: Logging function
  # Returns:
  #   Seurat object with gene symbols, counts, embeddings, and metadata
  if (!requireNamespace("schard", quietly = TRUE)) {
    stop("Package 'schard' is required to read .h5ad files. ",
         "Install via: remotes::install_github('cellgeni/schard')")
  }
  # ── Step 1: Pre-read h5ad var to detect Ensembl IDs and find gene symbol column ──
  gene_symbols <- NULL
  tryCatch({
    h5 <- hdf5r::H5File$new(h5ad_path, mode = "r")
    on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
    if (h5$exists("var")) {
      var_group <- h5[["var"]]
      var_cols <- names(var_group)
      # Read current var index (will become rownames in schard output)
      var_index <- NULL
      if ("_index" %in% var_cols) {
        idx_obj <- var_group[["_index"]]
        var_index <- if (inherits(idx_obj, "H5Group") && idx_obj$exists("categories") && idx_obj$exists("codes")) {
          idx_obj[["categories"]]$read()[idx_obj[["codes"]]$read() + 1L]
        } else {
          idx_obj$read()
        }
      }
      # Check if index looks like Ensembl IDs
      ensembl_pattern <- "^ENS[A-Z]*G[0-9]+"
      if (!is.null(var_index) && mean(grepl(ensembl_pattern, var_index)) > 0.5) {
        log_function(paste0("Detected Ensembl IDs in var index (", length(var_index), " genes)"), "info")
        # Search for gene symbol column
        symbol_candidates <- c("feature_name", "gene_short_name", "gene_symbols",
                               "gene_symbol", "gene_name", "symbol", "name")
        var_cols_lower <- tolower(var_cols)
        for (cand in symbol_candidates) {
          match_idx <- which(var_cols_lower == tolower(cand))
          if (length(match_idx) > 0) {
            col_name <- var_cols[match_idx[1]]
            col_obj <- var_group[[col_name]]
            values <- if (inherits(col_obj, "H5Group") && col_obj$exists("categories") && col_obj$exists("codes")) {
              col_obj[["categories"]]$read()[col_obj[["codes"]]$read() + 1L]
            } else {
              col_obj$read()
            }
            # Verify it contains actual gene symbols, not more Ensembl IDs
            if (mean(grepl(ensembl_pattern, values)) < 0.5 && mean(nchar(values) > 0) > 0.9) {
              log_function(paste("Found gene symbol column:", col_name), "info")
              # Handle NA/empty
              bad <- is.na(values) | values == ""
              values[bad] <- var_index[bad]
              # Handle duplicates
              n_dupes <- sum(duplicated(values))
              if (n_dupes > 0) {
                log_function(paste("Resolving", n_dupes, "duplicate gene symbols"), "info")
                values <- make.unique(values)
              }
              gene_symbols <- values
              log_function(paste("Examples:", paste(head(gene_symbols, 5), collapse = ", ")), "info")
              break
            }
          }
        }
        if (is.null(gene_symbols)) {
          log_function("No gene symbol column found — keeping Ensembl IDs", "warning")
          log_function(paste("Available var columns:", paste(var_cols, collapse = ", ")), "info")
        }
      } else {
        log_function("Gene names already appear to be symbols — no remapping needed", "info")
      }
    }
  }, error = function(e) {
    log_function(paste("Var pre-read failed:", conditionMessage(e)[1]), "warning")
  })
  # ── Step 2: Load via schard ──
  log_function("Reading H5AD via schard (this may take a while for large files)...", "info")
  seuratObj <- tryCatch({
    schard::h5ad2seurat(h5ad_path, use.raw = TRUE)
  }, error = function(e) {
    log_function(paste("schard with use.raw=TRUE failed:", conditionMessage(e)[1]), "warning")
    log_function("Falling back to use.raw=FALSE (X matrix as counts)", "info")
    tryCatch({
      schard::h5ad2seurat(h5ad_path, use.raw = FALSE)
    }, error = function(e2) {
      stop(paste("schard could not parse H5AD:", conditionMessage(e2)[1]))
    })
  })
  if (!inherits(seuratObj, "Seurat")) {
    stop("schard returned a non-Seurat object")
  }
  # ── Step 3: Apply gene symbol remapping if detected ──
  if (!is.null(gene_symbols) && length(gene_symbols) == nrow(seuratObj)) {
    log_function("Applying Ensembl → gene symbol remapping...", "info")
    # Extract counts, rename, rebuild — one shot, before Seurat caches anything
    assay_name <- names(seuratObj@assays)[1]
    assay_obj <- seuratObj[[assay_name]]
    if (inherits(assay_obj, "Assay")) {
      # v4 Assay: direct Dimnames access on sparse matrices
      if (length(assay_obj@counts) > 0) assay_obj@counts@Dimnames[[1]] <- gene_symbols
      if (length(assay_obj@data) > 0) assay_obj@data@Dimnames[[1]] <- gene_symbols
      if (nrow(assay_obj@meta.features) > 0) rownames(assay_obj@meta.features) <- gene_symbols
    } else if (inherits(assay_obj, "Assay5")) {
      # v5 Assay5: rename each layer + features metadata
      for (ln in names(assay_obj@layers)) {
        assay_obj@layers[[ln]]@Dimnames[[1]] <- gene_symbols
      }
      if (nrow(assay_obj@meta.data) > 0) rownames(assay_obj@meta.data) <- gene_symbols
    }
    seuratObj[[assay_name]] <- assay_obj
    # Verify
    actual_genes <- head(rownames(seuratObj), 5)
    if (any(grepl("^ENS[A-Z]*G[0-9]+", actual_genes))) {
      log_function("⚠️ Remapping may not have persisted — Ensembl IDs still detected", "warning")
    } else {
      log_function(paste0("✓ Remapped → gene symbols (", paste(actual_genes, collapse = ", "), ")"), "success")
    }
  } else if (!is.null(gene_symbols)) {
    log_function(paste("Gene count mismatch: mapping has", length(gene_symbols),
                       "but object has", nrow(seuratObj), "— skipping remap"), "warning")
  }
  # ── Step 4: Verify counts ──
  counts_layer <- tryCatch(
    LayerData(seuratObj, assay = names(seuratObj@assays)[1], layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(counts_layer) || length(counts_layer) == 0) {
    stop("No counts matrix found in H5AD - cannot proceed")
  }
  sample_values <- if (inherits(counts_layer, "dgCMatrix")) {
    head(counts_layer@x, 1000)
  } else {
    head(as.vector(counts_layer), 1000)
  }
  is_integer_like <- length(sample_values) > 0 && all(abs(sample_values - round(sample_values)) < 1e-9)
  if (!is_integer_like) {
    log_function("⚠️ counts matrix contains non-integer values", "warning")
    log_function("This H5AD likely stores normalized data - downstream QC may produce wrong results", "warning")
  } else {
    log_function("✓ Counts matrix appears integer-valued", "success")
  }
  # ── Step 5: Metadata ──
  seuratObj@project.name <- tools::file_path_sans_ext(basename(h5ad_path))
  if (!"percent.mt" %in% colnames(seuratObj@meta.data)) {
    log_function("Computing percent.mt (not present in metadata)...", "info")
    seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = mt_pattern)
  } else {
    log_function("✓ percent.mt already present in metadata", "info")
  }
  seuratObj <- cleanH5ADFactors(seuratObj, log_function)
  # ── Step 6: Auto-detect clustering column ──
  cluster_candidates <- c("leiden", "louvain", "cell_type", "celltype", "CellType",
                          "cell_type_ontology_term_id", "cluster", "clusters",
                          "seurat_clusters", "annotation", "cell_annotation",
                          "author_cell_type")
  found_cluster_col <- NULL
  for (cand in cluster_candidates) {
    if (cand %in% colnames(seuratObj@meta.data)) {
      found_cluster_col <- cand
      break
    }
  }
  if (!is.null(found_cluster_col)) {
    log_function(paste("Setting active identity from column:", found_cluster_col), "info")
    Idents(seuratObj) <- found_cluster_col
    n_clusters <- length(levels(Idents(seuratObj)))
    log_function(paste("Found", n_clusters, "clusters"), "success")
  } else {
    log_function("No standard clustering column found in metadata - Idents unchanged", "warning")
    log_function(paste("Available columns:", paste(colnames(seuratObj@meta.data), collapse = ", ")), "info")
  }
  # ── Step 7: Log results ──
  log_function(paste("Cells loaded:", format(ncol(seuratObj), big.mark = ",")), "success")
  log_function(paste("Genes loaded:", format(nrow(seuratObj), big.mark = ",")), "success")
  if (length(seuratObj@reductions) > 0) {
    log_function(paste("Imported reductions:", paste(names(seuratObj@reductions), collapse = ", ")), "info")
  }
  meta_cols <- colnames(seuratObj@meta.data)
  if (length(meta_cols) > 0) {
    n_show <- min(15, length(meta_cols))
    cols_show <- paste(meta_cols[seq_len(n_show)], collapse = ", ")
    if (length(meta_cols) > n_show) cols_show <- paste0(cols_show, ", ...")
    log_function(paste("Metadata columns:", cols_show), "info")
  }
  return(seuratObj)
}


slimSeuratObject <- function(seurat_obj, keep_counts = TRUE, keep_scale_data = FALSE,
                             log_function = message) {
  # Reduce memory footprint by removing redundant assays and layers
  # Keeps: RNA assay, all reductions, all graphs, full metadata
  # Drops: extra assays (SCT, integrated), scale.data by default
  # Args:
  #   seurat_obj: Seurat object
  #   keep_counts: Keep counts layer (needed for DEG, CellChat)
  #   keep_scale_data: Keep scale.data layer (needed for heatmaps — can be recomputed)
  #   log_function: Logging function
  # Returns:
  #   Slimmed Seurat object
  mem_before <- as.numeric(object.size(seurat_obj)) / 1024^2
  # Determine which layers to keep
  layers_to_keep <- "data"
  if (keep_counts) layers_to_keep <- c(layers_to_keep, "counts")
  if (keep_scale_data) layers_to_keep <- c(layers_to_keep, "scale.data")
  # Determine which assays to keep: RNA + default assay (may differ after SCTransform)
  default_assay <- DefaultAssay(seurat_obj)
  assays_to_keep <- unique(c("RNA", default_assay))
  all_assays <- names(seurat_obj@assays)
  dropped_assays <- setdiff(all_assays, assays_to_keep)
  if (length(dropped_assays) > 0) {
    log_function(paste("Dropping assays:", paste(dropped_assays, collapse = ", ")), "info")
  }
  # Keep all existing reductions and graphs
  existing_reductions <- names(seurat_obj@reductions)
  existing_graphs <- names(seurat_obj@graphs)
  # Check if all assays are v5 (Assay5) — DietSeurat layers param requires v5
  all_v5 <- all(vapply(names(seurat_obj@assays), function(a) {
    inherits(seurat_obj[[a]], "Assay5")
  }, logical(1)))
  if (all_v5) {
    seurat_obj <- DietSeurat(
      seurat_obj,
      layers = layers_to_keep,
      assays = intersect(assays_to_keep, all_assays),
      dimreducs = existing_reductions,
      graphs = existing_graphs
    )
  } else {
    log_function("Seurat v4 assay detected — skipping layer-level optimization", "info")
    seurat_obj <- DietSeurat(
      seurat_obj,
      assays = intersect(assays_to_keep, all_assays),
      dimreducs = existing_reductions,
      graphs = existing_graphs
    )
  }
  # Ensure default assay is set correctly after Diet
  if (default_assay %in% names(seurat_obj@assays)) {
    DefaultAssay(seurat_obj) <- default_assay
  } else {
    DefaultAssay(seurat_obj) <- names(seurat_obj@assays)[1]
  }
  gc()
  mem_after <- as.numeric(object.size(seurat_obj)) / 1024^2
  mem_saved <- mem_before - mem_after
  log_function(paste0("Memory optimization: ", round(mem_before), " MB → ",
                      round(mem_after), " MB (saved ", round(mem_saved), " MB)"), "success")
  return(seurat_obj)
}

logMemoryUsage <- function(log_function = message, seurat_obj = NULL, label = NULL) {
  # Log memory usage — works in Docker, local Linux, macOS, Windows
  # Uses only base R: gc(), file.exists(), readLines(), object.size()
  # Args:
  #   log_function: Logging function (e.g. log_system$add_log)
  #   seurat_obj: Optional Seurat object to measure
  #   label: Optional label for the object
  log_function("--- Memory Usage ---", "info")
  # R heap via gc() — works everywhere
  tryCatch({
    gc_result <- gc(verbose = FALSE, reset = FALSE)
    r_heap_mb <- round(sum(gc_result[, 2]), 1)
    log_function(paste0("R heap: ", format(r_heap_mb, big.mark = ","), " MB"), "info")
  }, error = function(e) NULL)
  # System RAM via /proc/meminfo — works on Linux (Docker + local Ubuntu)
  # Silently skipped on macOS/Windows
  tryCatch({
    if (file.exists("/proc/meminfo")) {
      meminfo <- readLines("/proc/meminfo", warn = FALSE)
      total_kb <- as.numeric(gsub("[^0-9]", "", grep("^MemTotal:", meminfo, value = TRUE)))
      avail_kb <- as.numeric(gsub("[^0-9]", "", grep("^MemAvailable:", meminfo, value = TRUE)))
      used_gb <- round((total_kb - avail_kb) / 1024^2, 1)
      total_gb <- round(total_kb / 1024^2, 1)
      avail_gb <- round(avail_kb / 1024^2, 1)
      pct <- round(100 * (total_kb - avail_kb) / total_kb, 0)
      log_function(paste0("System RAM: ", used_gb, " / ", total_gb, " GB (",
                          pct, "% used, ", avail_gb, " GB free)"), "info")
      if (pct > 85) {
        log_function("⚠️ RAM usage is high — consider enabling memory optimization", "warning")
      }
    }
  }, error = function(e) NULL)
  # Seurat object size via object.size() — works everywhere
  if (!is.null(seurat_obj) && inherits(seurat_obj, "Seurat")) {
    obj_mb <- round(as.numeric(object.size(seurat_obj)) / 1024^2, 1)
    obj_label <- if (!is.null(label)) label else "Seurat object"
    n_assays <- length(names(seurat_obj@assays))
    assay_names <- paste(names(seurat_obj@assays), collapse = ", ")
    n_layers <- sum(vapply(names(seurat_obj@assays), function(a) {
      assay <- seurat_obj[[a]]
      if (inherits(assay, "Assay5")) {
        length(names(assay@layers))
      } else {
        # Seurat v4 Assay: count non-empty slots as layers
        sum(c(length(assay@counts) > 0, length(assay@data) > 0, length(assay@scale.data) > 0))
      }
    }, integer(1)))

    log_function(paste0(obj_label, ": ~", format(obj_mb, big.mark = ","), " MB (",
                        n_assays, " assay(s): ", assay_names, ", ",
                        n_layers, " layer(s) total)"), "info")
  }
  log_function("--------------------", "info")
}


remapEnsemblToSymbols <- function(seurat_obj, h5ad_path, log_function = message) {
  # Detect Ensembl IDs as rownames and remap to gene symbols from h5ad var metadata
  # Common in public datasets (CELLxGENE, HCA, GEO)
  # Args:
  #   seurat_obj: Seurat object with potential Ensembl ID rownames
  #   h5ad_path: Path to original h5ad file (to read var columns)
  #   log_function: Logging function
  # Returns:
  #   Seurat object with gene symbols as rownames (if mapping found)
  current_genes <- rownames(seurat_obj)
  # Detect Ensembl IDs: ENSG (human), ENSMUSG (mouse), ENSRNOG (rat)
  ensembl_pattern <- "^ENS[A-Z]*G[0-9]+"
  pct_ensembl <- mean(grepl(ensembl_pattern, current_genes))
  if (pct_ensembl < 0.5) {
    log_function("Gene names appear to be symbols (not Ensembl IDs) — no remapping needed", "info")
    return(seurat_obj)
  }
  log_function(paste0("Detected Ensembl IDs as gene names (", round(pct_ensembl * 100), "%)"), "info")
  log_function("Searching for gene symbol column in H5AD var metadata...", "info")
  # Read var columns from h5ad to find gene symbols
  symbol_col <- NULL
  symbol_values <- NULL
  var_cols <- character(0)
  tryCatch({
    h5 <- hdf5r::H5File$new(h5ad_path, mode = "r")
    on.exit(try(h5$close_all(), silent = TRUE), add = TRUE)
    if (!h5$exists("var")) {
      log_function("No var group found in H5AD — cannot remap", "warning")
      return(seurat_obj)
    }
    var_group <- h5[["var"]]
    var_cols <- names(var_group)
    # Priority order for gene symbol columns (most common naming conventions)
    symbol_candidates <- c("feature_name", "gene_short_name", "gene_symbols",
                           "gene_symbol", "gene_name", "symbol", "name",
                           "gene_ids", "_index", "index")
    var_cols_lower <- tolower(var_cols)
    for (cand in symbol_candidates) {
      match_idx <- which(var_cols_lower == tolower(cand))
      if (length(match_idx) > 0) {
        col_name <- var_cols[match_idx[1]]
        col_obj <- var_group[[col_name]]
        # h5ad categorical columns are stored as groups with codes + categories
        if (inherits(col_obj, "H5Group")) {
          if (col_obj$exists("categories") && col_obj$exists("codes")) {
            categories <- col_obj[["categories"]]$read()
            codes <- col_obj[["codes"]]$read()
            values <- categories[codes + 1L]
          } else {
            next
          }
        } else {
          values <- col_obj$read()
        }
        # Check if this column looks like symbols (not Ensembl IDs themselves)
        if (mean(grepl(ensembl_pattern, values)) < 0.5 && mean(nchar(values) > 0) > 0.9) {
          symbol_col <- col_name
          symbol_values <- values
          break
        }
      }
    }
  }, error = function(e) {
    log_function(paste("Error reading var metadata:", conditionMessage(e)[1]), "warning")
  })
  if (is.null(symbol_col) || is.null(symbol_values)) {
    log_function("No gene symbol column found in var metadata — keeping Ensembl IDs", "warning")
    if (length(var_cols) > 0) {
      log_function(paste("Available var columns:", paste(var_cols, collapse = ", ")), "info")
    }
    return(seurat_obj)
  }
  log_function(paste("Found gene symbol column:", symbol_col), "info")
  if (length(symbol_values) != length(current_genes)) {
    log_function(paste("Length mismatch: var has", length(symbol_values),
                       "entries but object has", length(current_genes), "genes"), "warning")
    return(seurat_obj)
  }
  # Handle NA/empty symbols: fall back to Ensembl ID
  bad_idx <- is.na(symbol_values) | symbol_values == ""
  symbol_values[bad_idx] <- current_genes[bad_idx]
  # Handle duplicate symbols
  n_dupes <- sum(duplicated(symbol_values))
  if (n_dupes > 0) {
    log_function(paste("Resolving", n_dupes, "duplicate gene symbols with make.unique()"), "info")
    symbol_values <- make.unique(symbol_values)
  }
  # Apply renaming at the matrix level for each assay (handles both v4 Assay and v5 Assay5)
  for (assay_name in names(seurat_obj@assays)) {
    tryCatch({
      assay_obj <- seurat_obj[[assay_name]]
      if (inherits(assay_obj, "Assay5")) {
        # Seurat v5: rename rownames in each layer + feature metadata
        for (layer_name in names(assay_obj@layers)) {
          rownames(assay_obj@layers[[layer_name]]) <- symbol_values
        }
        if (nrow(assay_obj@meta.data) > 0) {
          rownames(assay_obj@meta.data) <- symbol_values
        }
        # Update cells (features are stored in @features if it exists)
        if (!is.null(assay_obj@features)) {
          rownames(assay_obj@features) <- symbol_values
        }
      } else {
        # Seurat v4: rename rownames in @counts, @data, @scale.data
        if (length(assay_obj@counts) > 0) {
          rownames(assay_obj@counts) <- symbol_values
        }
        if (length(assay_obj@data) > 0) {
          rownames(assay_obj@data) <- symbol_values
        }
        if (length(assay_obj@scale.data) > 0) {
          # scale.data may have fewer genes (only variable features)
          old_scaled_genes <- rownames(assay_obj@scale.data)
          new_scaled_genes <- symbol_values[match(old_scaled_genes, current_genes)]
          rownames(assay_obj@scale.data) <- new_scaled_genes
        }
      }
      seurat_obj[[assay_name]] <- assay_obj
    }, error = function(e) {
      log_function(paste("Could not rename features in assay", assay_name, ":",
                         conditionMessage(e)[1]), "warning")
    })
  }
  log_function(paste0("✓ Remapped ", length(current_genes), " Ensembl IDs → gene symbols"), "success")
  log_function(paste("Examples:", paste(head(symbol_values, 5), collapse = ", ")), "info")
  return(seurat_obj)
}

