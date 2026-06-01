# spatial_server.R - Server functions for Spatial Transcriptomics

spatial_transcriptomic_server <- function(input, output, session) {
  
  ############################## Reactive Variables ##############################
  
  shinyjs::useShinyjs()
  
  # Main spatial object storage
  spatial_obj <- reactiveVal(NULL)
  sketching_applied <- reactiveVal(FALSE)
  available_assays <- reactiveVal(c("Spatial"))
  # Simple assay selection storage
  selected_assay <- reactiveVal("Spatial")
  
  
  # Processing state flags
  processing_states <- reactiveValues(
    qc_applied = FALSE,
    normalized = FALSE,
    pca_computed = FALSE,
    clustered = FALSE
  )
  
  # Individual state flags (for compatibility)
  spatial_normalized <- reactiveVal(FALSE)
  spatial_qc_applied <- reactiveVal(FALSE)
  spatial_pca_computed <- reactiveVal(FALSE)
  spatial_clustered <- reactiveVal(FALSE)
  
  # Plots for download
  spatial_plots <- reactiveValues()
  
  ############################## Helper Functions ##############################
  
  # Update available assays when object changes
  observe({
    req(spatial_obj())
    obj <- spatial_obj()
    assays <- names(obj@assays)
    available_assays(assays)
    
    # Update assay selection menus
    updateSelectInput(session, "spatial_assay_select", choices = assays, selected = assays[1])
    updateSelectInput(session, "spatial_viz_assay", choices = assays, selected = if("RNA" %in% assays) "RNA" else assays[1])
  })
  
  
  # Update object and state in one function
  update_spatial_object <- function(new_obj, state_updates = list()) {
    spatial_obj(new_obj)
    for (state_name in names(state_updates)) {
      processing_states[[state_name]] <- state_updates[[state_name]]
      # Also update individual flags
      if (state_name == "normalized") spatial_normalized(state_updates[[state_name]])
      if (state_name == "qc_applied") spatial_qc_applied(state_updates[[state_name]])
      if (state_name == "pca_computed") spatial_pca_computed(state_updates[[state_name]])
      if (state_name == "clustered") spatial_clustered(state_updates[[state_name]])
    }
  }
  
  # Get current processing status
  get_processing_status <- function() {
    obj <- spatial_obj()
    if (is.null(obj)) return("No data loaded")
    
    status_parts <- c()
    if (processing_states$qc_applied) status_parts <- c(status_parts, "QC")
    if (processing_states$normalized) status_parts <- c(status_parts, "Normalized")
    if (sketching_applied()) status_parts <- c(status_parts, "Sketched")
    if (processing_states$pca_computed) status_parts <- c(status_parts, "PCA")
    if (processing_states$clustered) status_parts <- c(status_parts, "Clustered")
    
    if (length(status_parts) == 0) return("Raw data")
    return(paste(status_parts, collapse = " → "))
  }
  
  # Add debug information display
  output$spatial_debug_info <- renderText({
    req(spatial_obj())
    
    obj <- spatial_obj()
    current_assay <- input$spatial_viz_assay %||% DefaultAssay(obj)
    
    debug_info <- paste(
      "Current assay:", current_assay,
      "\nAvailable assays:", paste(names(obj@assays), collapse = ", "),
      "\nDefault assay:", DefaultAssay(obj),
      "\nGenes in current assay:", nrow(obj[[current_assay]]),
      "\nHas spatial images:", length(obj@images) > 0
    )
    
    return(debug_info)
  })
  # Get spatial dataset name
  get_spatial_dataset_name <- function() {
    if (!is.null(spatial_obj())) {
      if ("orig.ident" %in% colnames(spatial_obj()@meta.data)) {
        return(unique(spatial_obj()@meta.data$orig.ident)[1])
      }
    }
    return("SpatialDataset")
  }
  
  # Validate spatial data structure
  validate_spatial_data <- function(seurat_obj) {
    required_assays <- c("Spatial")
    has_spatial <- any(required_assays %in% names(seurat_obj@assays))
    has_images <- length(seurat_obj@images) > 0
    
    has_coordinates <- TRUE
    if (has_images) {
      tryCatch({
        if (class(seurat_obj@images[[1]])[1] == "VisiumV2") {
          coords <- GetTissueCoordinates(seurat_obj@images[[1]])
          has_coordinates <- !is.null(coords) && nrow(coords) > 0
        } else {
          has_coordinates <- !is.null(seurat_obj@images[[1]]@coordinates)
        }
      }, error = function(e) {
        message("Could not validate coordinates, assuming they exist")
        has_coordinates <- TRUE
      })
    }
    
    return(list(
      valid = has_spatial,
      has_images = has_images,
      has_coordinates = has_coordinates,
      message = if (!has_spatial) "Missing Spatial assay" else "Valid spatial data"
    ))
  }
  
  # Safely extract ZIP file
  safe_extract_zip <- function(zip_path, extract_dir) {
    if (!file.exists(zip_path)) {
      stop("ZIP file does not exist at the specified path")
    }
    
    file_size <- file.info(zip_path)$size
    if (is.na(file_size) || file_size == 0) {
      stop("ZIP file appears to be empty or corrupted")
    }
    
    tryCatch({
      zip_contents <- unzip(zip_path, list = TRUE)
      if (nrow(zip_contents) == 0) {
        stop("ZIP file contains no files")
      }
      
      extracted_files <- unzip(zip_path, exdir = extract_dir)
      return(list(success = TRUE, contents = zip_contents, extracted_files = extracted_files))
      
    }, error = function(e) {
      stop(paste("Failed to process ZIP file:", e$message))
    })
  }
  
  # Organize spatial files into expected structure
  organize_spatial_files <- function(data_dir) {
    all_files <- list.files(data_dir, full.names = TRUE, recursive = FALSE)
    file_names <- basename(all_files)
    
    spatial_dir <- file.path(data_dir, "spatial")
    if (!dir.exists(spatial_dir)) {
      dir.create(spatial_dir, recursive = TRUE)
    }
    
    # Define file patterns for different spatial files
    file_patterns <- list(
      expression_h5 = list(
        patterns = c("filtered_feature_bc_matrix\\.h5", ".*feature.*matrix.*\\.h5", ".*\\.h5"),
        destination = "filtered_feature_bc_matrix.h5",
        target_dir = data_dir
      ),
      scalefactors = list(
        patterns = c("scalefactors_json.*\\.json", ".*scalefactors.*\\.json", ".*\\.json"),
        destination = "scalefactors_json.json",
        target_dir = spatial_dir
      ),
      tissue_positions_csv = list(
        patterns = c("tissue_positions.*\\.csv", ".*positions.*\\.csv", ".*tissue.*\\.csv"),
        destination = "tissue_positions_list.csv",
        target_dir = spatial_dir
      ),
      tissue_positions_parquet = list(
        patterns = c("tissue_positions.*\\.parquet", ".*positions.*\\.parquet"),
        destination = "tissue_positions.parquet",
        target_dir = spatial_dir
      ),
      hires_image = list(
        patterns = c("tissue_hires_image.*\\.png", ".*hires.*\\.png", ".*high.*\\.png"),
        destination = "tissue_hires_image.png",
        target_dir = spatial_dir
      ),
      lowres_image = list(
        patterns = c("tissue_lowres_image.*\\.png", ".*lowres.*\\.png", ".*low.*\\.png"),
        destination = "tissue_lowres_image.png",
        target_dir = spatial_dir
      )
    )
    
    # Match and move files
    moved_files <- list()
    for (file_type in names(file_patterns)) {
      pattern_info <- file_patterns[[file_type]]
      
      for (pattern in pattern_info$patterns) {
        matching_files <- grep(pattern, file_names, ignore.case = TRUE, value = TRUE)
        
        if (length(matching_files) > 0) {
          source_file <- file.path(data_dir, matching_files[1])
          dest_file <- file.path(pattern_info$target_dir, pattern_info$destination)
          
          if (file.copy(source_file, dest_file, overwrite = TRUE)) {
            moved_files[[file_type]] <- dest_file
            message(paste("Organized:", matching_files[1], "->", pattern_info$destination))
            break
          }
        }
      }
    }
    
    return(moved_files)
  }
  
  # Detect and organize spatial data structure
  detect_spatial_data_structure_flexible <- function(data_dir) {
    # First, try to organize files
    organized_files <- organize_spatial_files(data_dir)
    
    # Check for spatial directory
    spatial_dir <- file.path(data_dir, "spatial")
    if (!dir.exists(spatial_dir)) {
      return(list(valid = FALSE, message = "Could not create or find spatial directory"))
    }
    
    # Check required files
    required_files <- c(
      expression = "filtered_feature_bc_matrix.h5",
      scalefactors = file.path("spatial", "scalefactors_json.json")
    )
    
    missing_files <- c()
    found_files <- c()
    
    for (file_name in required_files) {
      file_path <- file.path(data_dir, file_name)
      if (file.exists(file_path)) {
        found_files <- c(found_files, file_name)
      } else {
        missing_files <- c(missing_files, file_name)
      }
    }
    
    # Check for tissue positions (either format)
    tissue_pos_parquet <- file.path(spatial_dir, "tissue_positions.parquet")
    tissue_pos_csv <- file.path(spatial_dir, "tissue_positions_list.csv")
    
    tissue_positions_found <- FALSE
    tissue_positions_format <- "none"
    
    if (file.exists(tissue_pos_parquet)) {
      tissue_positions_found <- TRUE
      tissue_positions_format <- "parquet"
    } else if (file.exists(tissue_pos_csv)) {
      tissue_positions_found <- TRUE
      tissue_positions_format <- "csv"
    }
    
    if (!tissue_positions_found) {
      missing_files <- c(missing_files, "tissue_positions file")
    }
    
    # Check for images (optional)
    has_images <- file.exists(file.path(spatial_dir, "tissue_hires_image.png")) ||
      file.exists(file.path(spatial_dir, "tissue_lowres_image.png"))
    
    # Determine if structure is valid
    valid <- length(missing_files) == 0
    
    if (valid) {
      message <- paste("Valid spatial data structure created. Found:", paste(found_files, collapse = ", "))
    } else {
      message <- paste("Missing required files:", paste(missing_files, collapse = ", "))
    }
    
    return(list(
      valid = valid,
      format = "H5",
      tissue_positions_format = tissue_positions_format,
      has_images = has_images,
      organized_files = organized_files,
      message = message
    ))
  }
  
  ############################## Load Spatial Dataset ##############################
  
  # Enhanced spatial data loading function
  observeEvent(input$spatial_file, {
    req(input$spatial_file)
    
    message("=== SPATIAL FILE UPLOAD STARTED ===")
    
    tryCatch({
      showModal(modalDialog(
        title = "Loading Spatial Data",
        "Processing your spatial transcriptomics data...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      file_info <- input$spatial_file
      temp_dir <- file.path(tempdir(), paste0("spatial_", format(Sys.time(), "%Y%m%d_%H%M%S")))
      dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Extract ZIP
      message("Extracting ZIP file...")
      extracted_files <- unzip(file_info$datapath, exdir = temp_dir)
      message(paste("Extracted files:", paste(basename(extracted_files), collapse = ", ")))
      
      # Get all files recursively to handle nested structures
      all_files <- list.files(temp_dir, full.names = TRUE, recursive = TRUE)
      all_basenames <- basename(all_files)
      
      message(paste("All extracted files:", paste(all_basenames, collapse = ", ")))
      
      # Enhanced file detection with multiple patterns
      h5_files <- grep("(filtered_feature_bc_matrix|matrix|expression).*\\.h5$", all_files, value = TRUE, ignore.case = TRUE)
      json_files <- grep("scalefactors.*\\.json$", all_files, value = TRUE, ignore.case = TRUE)
      csv_files <- grep("tissue_positions.*\\.csv$", all_files, value = TRUE, ignore.case = TRUE)
      parquet_files <- grep("tissue_positions.*\\.parquet$", all_files, value = TRUE, ignore.case = TRUE)
      hires_images <- grep("tissue_hires_image.*\\.png$", all_files, value = TRUE, ignore.case = TRUE)
      lowres_images <- grep("tissue_lowres_image.*\\.png$", all_files, value = TRUE, ignore.case = TRUE)
      
      # Alternative patterns for different resolutions
      if (length(h5_files) == 0) {
        h5_files <- grep("\\.h5$", all_files, value = TRUE)
        message("Using broader H5 file detection")
      }
      
      if (length(json_files) == 0) {
        json_files <- grep("\\.json$", all_files, value = TRUE)
        message("Using broader JSON file detection")
      }
      
      if (length(csv_files) == 0) {
        csv_files <- grep("positions.*\\.csv$", all_files, value = TRUE, ignore.case = TRUE)
        message("Using broader CSV file detection")
      }
      
      if (length(parquet_files) == 0) {
        parquet_files <- grep("positions.*\\.parquet$", all_files, value = TRUE, ignore.case = TRUE)
        message("Using broader Parquet file detection")
      }
      
      message(paste("Enhanced detection - H5:", length(h5_files), "JSON:", length(json_files), 
                    "CSV:", length(csv_files), "Parquet:", length(parquet_files)))
      
      if (length(h5_files) == 0) {
        message("Available files:", paste(all_basenames, collapse = ", "))
        stop("No H5 matrix files found in uploaded data. Please check file structure.")
      }
      
      # Create spatial directory structure
      spatial_dir <- file.path(temp_dir, "spatial")
      if (!dir.exists(spatial_dir)) {
        dir.create(spatial_dir, recursive = TRUE)
      }
      
      # Copy H5 file to root with standard name
      h5_target <- file.path(temp_dir, "filtered_feature_bc_matrix.h5")
      if (!file.exists(h5_target)) {
        file.copy(h5_files[1], h5_target, overwrite = TRUE)
        message(paste("Copied H5 file:", basename(h5_files[1]), "to root"))
      }
      
      # Handle JSON file
      if (length(json_files) > 0) {
        json_target <- file.path(spatial_dir, "scalefactors_json.json")
        if (!file.exists(json_target)) {
          file.copy(json_files[1], json_target, overwrite = TRUE)
          message(paste("Copied JSON file:", basename(json_files[1])))
        }
      } else {
        # Create minimal scalefactors if missing (common issue with some datasets)
        json_target <- file.path(spatial_dir, "scalefactors_json.json")
        minimal_scalefactors <- list(
          tissue_hires_scalef = 1.0,
          tissue_lowres_scalef = 1.0,
          fiducial_diameter_fullres = 1.0,
          spot_diameter_fullres = 1.0
        )
        write(jsonlite::toJSON(minimal_scalefactors, auto_unbox = TRUE), file = json_target)
        message("Created minimal scalefactors file")
      }
      
      # Handle tissue positions - prefer Parquet over CSV
      positions_copied <- FALSE
      if (length(parquet_files) > 0) {
        parquet_target <- file.path(spatial_dir, "tissue_positions.parquet")
        if (!file.exists(parquet_target)) {
          file.copy(parquet_files[1], parquet_target, overwrite = TRUE)
          message(paste("Copied Parquet positions file:", basename(parquet_files[1])))
          positions_copied <- TRUE
        }
      }
      
      if (!positions_copied && length(csv_files) > 0) {
        csv_target <- file.path(spatial_dir, "tissue_positions_list.csv")
        if (!file.exists(csv_target)) {
          file.copy(csv_files[1], csv_target, overwrite = TRUE)
          message(paste("Copied CSV positions file:", basename(csv_files[1])))
          positions_copied <- TRUE
        }
      }
      
      if (!positions_copied) {
        message("Warning: No tissue positions file found. This may cause issues with spatial visualization.")
      }
      
      # Handle image files
      if (length(hires_images) > 0) {
        hires_target <- file.path(spatial_dir, "tissue_hires_image.png")
        if (!file.exists(hires_target)) {
          file.copy(hires_images[1], hires_target, overwrite = TRUE)
          message("Copied hires image")
        }
      }
      
      if (length(lowres_images) > 0) {
        lowres_target <- file.path(spatial_dir, "tissue_lowres_image.png")
        if (!file.exists(lowres_target)) {
          file.copy(lowres_images[1], lowres_target, overwrite = TRUE)
          message("Copied lowres image")
        }
      }
      
      # Enhanced loading with error handling for different formats
      message("Loading spatial data...")
      
      # Try different loading approaches
      spatial_data <- NULL
      loading_success <- FALSE
      
      # Approach 1: Standard Load10X_Spatial
      if (!loading_success) {
        tryCatch({
          message("Attempting standard Load10X_Spatial...")
          spatial_data <- Load10X_Spatial(
            data.dir = temp_dir,
            filename = "filtered_feature_bc_matrix.h5",
            assay = "Spatial",
            slice = "slice1",
            filter.matrix = TRUE,
            to.upper = FALSE
          )
          loading_success <- TRUE
          message("Standard loading successful")
        }, error = function(e) {
          message(paste("Standard loading failed:", e$message))
        })
      }
      
      # Approach 2: Try without filtering
      if (!loading_success) {
        tryCatch({
          message("Attempting Load10X_Spatial without filtering...")
          spatial_data <- Load10X_Spatial(
            data.dir = temp_dir,
            filename = "filtered_feature_bc_matrix.h5",
            assay = "Spatial",
            slice = "slice1",
            filter.matrix = FALSE,
            to.upper = FALSE
          )
          loading_success <- TRUE
          message("Loading without filtering successful")
        }, error = function(e) {
          message(paste("Loading without filtering failed:", e$message))
        })
      }
      
      # Approach 3: Manual loading for problematic datasets
      if (!loading_success) {
        tryCatch({
          message("Attempting manual H5 loading...")
          # Load expression matrix
          h5_data <- Read10X_h5(h5_target)
          
          # Create Seurat object
          spatial_data <- CreateSeuratObject(
            counts = h5_data,
            assay = "Spatial",
            project = "SpatialData"
          )
          
          # Try to add spatial information if available
          if (positions_copied) {
            # This is a simplified approach - may need adjustment for your specific data
            message("Adding spatial coordinates...")
            # --- REPLACE THIS TODO WITH THE BLOCK BELOW ---
            spatial_dir <- file.path(temp_dir, "spatial")
            
            #Load the image + Visium coordinates
            sp_img <- tryCatch({
              Seurat::Read10X_Image(image.dir = spatial_dir, strip.suffix = TRUE)
            }, error = function(e) {
              Seurat::Read10X_Image(image.dir = spatial_dir, strip.suffix = FALSE)
            })
            
            # Align barcodes (with/without suffix “-1”)
            coord_bc <- rownames(sp_img@coordinates)
            obj_bc   <- colnames(spatial_data)
            common   <- intersect(coord_bc, obj_bc)
            if (length(common) == 0) {
              coord_bc2 <- paste0(coord_bc, "-1")
              rownames(sp_img@coordinates) <- coord_bc2
              common <- intersect(coord_bc2, obj_bc)
            }
            if (length(common) == 0) stop("Aucun barcode en commun entre coords et matrice.")
            
            # Restrict the image to the barcodes present in the object
            sp_img <- sp_img[common, ]
            
            #Attach image/coords to Seurat object
            spatial_data[["slice1"]] <- sp_img
            if ("Spatial" %in% Seurat::Assays(spatial_data)) {
              Seurat::DefaultAssay(spatial_data) <- "Spatial"
            }
            
            # Propagate ‘in_tissue’ as meta for QC
            coords <- spatial_data@images[["slice1"]]@coordinates
            mtch   <- match(colnames(spatial_data), rownames(coords))
            spatial_data$in_tissue <- 0L
            spatial_data$in_tissue[!is.na(mtch)] <- as.integer(coords$tissue[mtch])
          }
          
          loading_success <- TRUE
          message("Manual loading successful")
        }, error = function(e) {
          message(paste("Manual loading failed:", e$message))
        })
      }
      
      if (!loading_success || is.null(spatial_data)) {
        stop("All loading methods failed. Please check your data format and structure.")
      }
      
      message(paste("Spatial data loaded successfully:", ncol(spatial_data), "spots,", nrow(spatial_data), "genes"))
      
      # Check if we have actual data
      if (ncol(spatial_data) == 0) {
        stop("No cells/spots found in the loaded data. This may indicate a format issue or empty dataset.")
      }
      
      # Add mitochondrial genes based on species
      if (input$spatial_species_choice == "mouse") {
        spatial_data[["percent.mt"]] <- PercentageFeatureSet(spatial_data, pattern = "^mt-|^Mt-")
      } else if (input$spatial_species_choice == "human") {
        spatial_data[["percent.mt"]] <- PercentageFeatureSet(spatial_data, pattern = "^MT-")
      } else if (input$spatial_species_choice == "rat") {
        spatial_data[["percent.mt"]] <- PercentageFeatureSet(spatial_data, pattern = "^Mt-")
      }
      
      # Update object
      update_spatial_object(spatial_data, list(
        qc_applied = FALSE, 
        normalized = FALSE, 
        pca_computed = FALSE, 
        clustered = FALSE
      ))
      
      # Success message
      n_spots <- ncol(spatial_data)
      n_genes <- nrow(spatial_data)
      
      showNotification(
        paste0("Spatial data loaded successfully! Spots: ", format(n_spots, big.mark = ","), 
               ", Genes: ", format(n_genes, big.mark = ",")), 
        type = "message", 
        duration = 8
      )
      removeModal()
      
      message("=== SPATIAL LOADING COMPLETED ===")
      
    }, error = function(e) {
      message(paste("=== SPATIAL LOADING ERROR ===", e$message))
      removeModal()
      showNotification(
        paste("Error loading spatial data:", e$message, 
              "\nPlease check your file structure and ensure all required files are included."), 
        type = "error",
        duration = 15
      )
    })
  })
  
  
  
  # Variable réactive pour stocker l'image uploadée
  uploaded_image_data <- reactiveVal(NULL)
  
  observeEvent(input$upload_high_res_image, {
    req(input$upload_high_res_image, spatial_obj())
    
    file_path <- input$upload_high_res_image$datapath
    file_ext <- tolower(tools::file_ext(input$upload_high_res_image$name))
    
    tryCatch({
      message("=== LOADING HIGH-RES IMAGE ===")
      
      if(file_ext == "png") {
        img_data <- png::readPNG(file_path)
      } else if(file_ext %in% c("jpg", "jpeg")) {
        img_data <- jpeg::readJPEG(file_path)
      } else if(file_ext == "tiff") {
        img_data <- tiff::readTIFF(file_path)
      } else {
        stop("Unsupported image format")
      }
      
      message("Image loaded: ", paste(dim(img_data), collapse = " x "))
      
      # S'assurer que l'image est en format 0-1
      if(max(img_data) > 1) {
        img_data <- img_data / 255
      }
      
      current_obj <- spatial_obj()
      
      if(length(current_obj@images) == 0) {
        stop("No spatial images found")
      }
      
      # RÉCUPÉRER L'OBJET IMAGE SPATIAL
      image_name <- names(current_obj@images)[1]
      spatial_image <- current_obj@images[[image_name]]
      
      message("Current image class: ", class(spatial_image))
      message("Current image dimensions: ", paste(dim(spatial_image@image), collapse = " x "))
      
      # REMPLACER L'IMAGE
      spatial_image@image <- img_data
      
      # AJUSTER LES SCALE FACTORS
      old_dims <- dim(current_obj@images[[image_name]]@image)[1:2]
      new_dims <- dim(img_data)[1:2]
      
      scale_ratio <- mean(new_dims / old_dims)
      
      # Créer de nouveaux scale factors
      new_scale_factors <- spatial_image@scale.factors
      
      if("lowres" %in% names(new_scale_factors)) {
        new_scale_factors$lowres <- new_scale_factors$lowres * scale_ratio
      }
      if("hires" %in% names(new_scale_factors)) {
        new_scale_factors$hires <- scale_ratio
      }
      
      spatial_image@scale.factors <- new_scale_factors
      
      # REMETTRE L'OBJET IMAGE MODIFIÉ DANS SEURAT
      current_obj@images[[image_name]] <- spatial_image
      
      # VÉRIFIER LA COORDONNÉE TISSUE
      tissue_coords <- GetTissueCoordinates(spatial_image)
      message("Tissue coordinates range: x[", min(tissue_coords[,1]), ",", max(tissue_coords[,1]), "] y[", min(tissue_coords[,2]), ",", max(tissue_coords[,2]), "]")
      
      spatial_obj(current_obj)
      
      message("=== IMAGE REPLACED SUCCESSFULLY ===")
      showNotification(paste("High-res image loaded:", paste(new_dims, collapse = "x")), type = "message")
      
    }, error = function(e) {
      message("ERROR: ", e$message)
      message("Traceback: ")
      print(traceback())
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  
  
  # Observer pour afficher les informations de l'image actuelle
  output$current_image_info <- renderText({
    req(spatial_obj())
    
    if(length(spatial_obj()@images) == 0) {
      return("No spatial images found")
    }
    
    img <- spatial_obj()@images[[1]]@image
    scale_factors <- spatial_obj()@images[[1]]@scale.factors
    
    info <- paste(
      "Current Image Information:",
      paste("- Dimensions:", paste(dim(img), collapse = " x ")),
      paste("- Size (MB):", round(as.numeric(object.size(img)) / 1024^2, 2)),
      paste("- Scale factors:", paste(names(scale_factors), scale_factors, sep = "=", collapse = ", ")),
      sep = "\n"
    )
    
    return(info)
  })
  
  
  ############################## Dataset Information Display ##############################
  
  output$spatial_spots_count <- renderInfoBox({
    count <- if (!is.null(spatial_obj())) ncol(spatial_obj()) else 0
    infoBox(
      title = "Spots/Cells",
      value = format(count, big.mark = ","),
      icon = icon("circle"),
      color = "fuchsia",  
      fill = TRUE
    )
  })
  
  output$spatial_genes_count <- renderInfoBox({
    count <- if (!is.null(spatial_obj())) nrow(spatial_obj()) else 0
    infoBox(
      title = "Genes",
      value = format(count, big.mark = ","),
      icon = icon("dna"),
      color = "teal",  
      fill = TRUE
    )
  })
  
  
  # Dataset summary table
  output$spatial_dataset_summary <- renderDT({
    req(spatial_obj())
    
    obj <- spatial_obj()
    
    summary_data <- data.frame(
      Metric = c(
        "Total Spots/Cells",
        "Total Genes", 
        "Available Assays",
        "Sketched Data",
        "Median Genes per Spot",
        "Median UMI per Spot", 
        "Mean Mitochondrial %",
        "Has Tissue Image",
        "Processing Status"
      ),
      Value = c(
        format(ncol(obj), big.mark = ","),
        format(nrow(obj), big.mark = ","),
        paste(names(obj@assays), collapse = ", "),
        if (sketching_applied()) "Yes" else "No",
        round(median(obj$nFeature_Spatial), 0),
        format(round(median(obj$nCount_Spatial), 0), big.mark = ","),
        paste0(round(mean(obj$percent.mt, na.rm = TRUE), 2), "%"),
        if (length(obj@images) > 0) "Yes" else "No",
        get_processing_status()
      ),
      stringsAsFactors = FALSE
    )
    
    datatable(summary_data, options = list(pageLength = 10, searching = FALSE, lengthChange = FALSE, info = FALSE), rownames = FALSE)
  })
  
  
  # Load processed Seurat object from RDS file
  observeEvent(input$load_spatial_seurat_file, {
    req(input$load_spatial_seurat_file)
    
    message("=== LOADING SPATIAL SEURAT OBJECT FROM RDS ===")
    
    tryCatch({
      showModal(modalDialog(
        title = "Loading Spatial Object",
        "Reading your saved spatial Seurat object...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      # Read the RDS file
      file_path <- input$load_spatial_seurat_file$datapath
      message(paste("Loading file from:", file_path))
      
      # Load the object
      loaded_obj <- readRDS(file_path)
      message("RDS file loaded successfully")
      
      # Validate that it's a Seurat object
      if (!inherits(loaded_obj, "Seurat")) {
        stop("The loaded file is not a valid Seurat object")
      }
      
      # Check if it has spatial data
      has_spatial <- FALSE
      if ("Spatial" %in% names(loaded_obj@assays)) {
        has_spatial <- TRUE
        message("Found Spatial assay")
      }
      
      # Check for images
      has_images <- length(loaded_obj@images) > 0
      if (has_images) {
        message(paste("Found", length(loaded_obj@images), "spatial image(s)"))
      } else {
        message("No spatial images found")
      }
      
      # Get basic info
      n_spots <- ncol(loaded_obj)
      n_genes <- nrow(loaded_obj)
      
      message(paste("Object contains:", n_spots, "spots and", n_genes, "genes"))
      
      # Check what processing has been done
      has_pca <- "pca" %in% names(loaded_obj@reductions)
      has_umap <- "umap" %in% names(loaded_obj@reductions)
      has_clusters <- "seurat_clusters" %in% colnames(loaded_obj@meta.data)
      has_normalized <- FALSE
      
      # Check if data is normalized by looking at the assay data
      if ("Spatial" %in% names(loaded_obj@assays)) {
        assay_data <- GetAssayData(loaded_obj, assay = "Spatial", slot = "data")
        # If max value is much larger than raw counts would be, it's likely normalized
        max_val <- max(assay_data@x[1:min(1000, length(assay_data@x))])
        has_normalized <- max_val < 50  # Normalized data typically has smaller values
      }
      
      # Add mitochondrial percentage if not present
      if (!"percent.mt" %in% colnames(loaded_obj@meta.data)) {
        message("Adding mitochondrial percentage...")
        if (input$spatial_species_choice == "mouse") {
          loaded_obj[["percent.mt"]] <- PercentageFeatureSet(loaded_obj, pattern = "^mt-|^Mt-")
        } else if (input$spatial_species_choice == "human") {
          loaded_obj[["percent.mt"]] <- PercentageFeatureSet(loaded_obj, pattern = "^MT-")
        } else if (input$spatial_species_choice == "rat") {
          loaded_obj[["percent.mt"]] <- PercentageFeatureSet(loaded_obj, pattern = "^Mt-")
        }
      }
      
      # Update the spatial object and processing states
      update_spatial_object(loaded_obj, list(
        qc_applied = has_normalized,  # Assume QC was done if normalized
        normalized = has_normalized,
        pca_computed = has_pca,
        clustered = has_clusters
      ))
      
      # Update sketching status if sketch assay exists
      if ("sketch" %in% names(loaded_obj@assays)) {
        sketching_applied(TRUE)
        message("Found sketch assay")
      }
      
      # Update available assays
      assay_names <- names(loaded_obj@assays)
      available_assays(assay_names)
      updateSelectInput(session, "spatial_assay_select", choices = assay_names, selected = DefaultAssay(loaded_obj))
      updateSelectInput(session, "spatial_viz_assay", choices = assay_names, selected = DefaultAssay(loaded_obj))
      
      # Success message
      info_parts <- c()
      if (has_spatial) info_parts <- c(info_parts, "Spatial data")
      if (has_images) info_parts <- c(info_parts, "Images")
      if (has_normalized) info_parts <- c(info_parts, "Normalized")
      if (has_pca) info_parts <- c(info_parts, "PCA")
      if (has_umap) info_parts <- c(info_parts, "UMAP")
      if (has_clusters) info_parts <- c(info_parts, "Clusters")
      
      success_msg <- paste0(
        "Spatial object loaded successfully!\n",
        "Spots: ", format(n_spots, big.mark = ","), "\n",
        "Genes: ", format(n_genes, big.mark = ","), "\n",
        if (length(info_parts) > 0) paste("Contains:", paste(info_parts, collapse = ", ")) else ""
      )
      
      removeModal()
      
      message("=== SPATIAL OBJECT LOADING COMPLETED ===")
      
    }, error = function(e) {
      message(paste("=== SPATIAL OBJECT LOADING ERROR ===", e$message))
      removeModal()
      showNotification(
        paste("Error loading spatial object:", e$message), 
        type = "error",
        duration = 10
      )
    })
  })
  
  
  
  ############################## Spatial QC & Normalization ##############################
  
  # Info box showing current spot count
  output$spatial_spots_info <- renderInfoBox({
    count <- if (!is.null(spatial_obj())) ncol(spatial_obj()) else 0
    infoBox(
      title = "Current Spots",
      value = format(count, big.mark = ","),
      icon = icon("circle"),
      color = "blue",
      width = 12
    )
  })
  
  # Generate QC plots
  # Generate QC plots - FIXED VERSION with proper layer handling
  observeEvent(input$spatial_qc_plots, {
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      options(warn = -1)  # Suppress warnings for plotting
      
      # Check available layers in the Spatial assay
      spatial_assay <- obj[["Spatial"]]
      available_layers <- names(spatial_assay@layers)
      message(paste("Available layers in Spatial assay:", paste(available_layers, collapse = ", ")))
      
      # Use counts layer explicitly for QC metrics
      if ("counts" %in% available_layers) {
        message("Using 'counts' layer for QC calculations")
      }
      
      output$spatial_tissue_plot <- renderPlot({
        if (length(obj@images) > 0) {
          tryCatch({
            SpatialFeaturePlot(obj, features = "nCount_Spatial", 
                               pt.size.factor = 1.2, 
                               alpha = c(0.1, 1)) +
              theme(legend.position = "bottom") +
              ggtitle("UMI Count Distribution")
          }, error = function(e) {
            message(paste("Spatial plot error:", e$message))
            # Fallback plot if spatial plot fails
            ggplot(obj@meta.data, aes(x = nCount_Spatial, y = nFeature_Spatial)) +
              geom_point(alpha = 0.6, size = 0.5) +
              labs(title = "Spatial QC Metrics", x = "UMI Count", y = "Gene Count") +
              theme_minimal()
          })
        } else {
          ggplot(obj@meta.data, aes(x = nCount_Spatial, y = nFeature_Spatial)) +
            geom_point(alpha = 0.6, size = 0.5) +
            labs(title = "Spatial QC Metrics", x = "UMI Count", y = "Gene Count") +
            theme_minimal()
        }
      }, height = 400, width = 600)  
      
      # Standard QC violin plots
      output$spatial_vlnplot <- renderPlot({
        VlnPlot(obj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),
                ncol = 3, pt.size = 0) +
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      }, height = 450, width = 900)  # Dimensions pour 3 colonnes
      
      output$spatial_qc_violin <- renderPlot({
        VlnPlot(obj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),
                ncol = 1, pt.size = 0) &
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }, height = 400, width = 600)
      
      # Feature scatter plots
      output$spatial_feature_scatter <- renderPlot({
        plot1 <- FeatureScatter(obj, feature1 = "nCount_Spatial", feature2 = "percent.mt", 
                                pt.size = 0.5) + NoLegend()
        plot2 <- FeatureScatter(obj, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", 
                                pt.size = 0.5) + NoLegend()
        plot1 + plot2
      }, height = 400, width = 800)
      
      output$spatial_scatter1 <- renderPlot({
        FeatureScatter(obj, feature1 = "nCount_Spatial", feature2 = "percent.mt", 
                       pt.size = 0.5)
      }, height = 400, width = 500)
      
      output$spatial_scatter2 <- renderPlot({
        FeatureScatter(obj, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", 
                       pt.size = 0.5)
      }, height = 400, width = 500)
      
      options(warn = 0) 
      showNotification("QC plots generated successfully!", type = "message")
    }, error = function(e) {
      options(warn = 0)
      showNotification(paste("Error generating QC plots:", e$message), type = "error")
    })
  })
  
  analyze_data_distribution <- function(obj) {
    features <- obj$nFeature_Spatial
    counts <- obj$nCount_Spatial
    mt <- obj$percent.mt
    features[is.na(features)] <- 0
    counts[is.na(counts)] <- 0
    mt[is.na(mt)] <- 0
    suggested_filters <- list(
      min_features = max(1, round(quantile(features[features > 0], 0.25, na.rm = TRUE) * 0.5)),
      max_features = min(round(quantile(features[features > 0], 0.95, na.rm = TRUE)), 1000),
      min_counts = max(1, round(quantile(counts[counts > 0], 0.25, na.rm = TRUE) * 0.5)),
      max_counts = min(round(quantile(counts[counts > 0], 0.95, na.rm = TRUE)), 5000),
      max_mt = min(round(quantile(mt, 0.95, na.rm = TRUE)), 50)
    )
    would_pass <- features >= suggested_filters$min_features & 
      features <= suggested_filters$max_features &
      counts >= suggested_filters$min_counts &
      counts <= suggested_filters$max_counts &
      mt <= suggested_filters$max_mt
    spots_passing <- sum(would_pass, na.rm = TRUE)
    return(list(suggested_filters = suggested_filters, spots_passing = spots_passing))
  }
  
  observeEvent(input$apply_spatial_qc, {
    req(spatial_obj())
    tryCatch({
      showModal(modalDialog(
        title = "Applying QC Filters",
        "Analyzing and filtering spatial data...",
        easyClose = FALSE,
        footer = NULL
      ))
      obj <- spatial_obj()
      original_spots <- ncol(obj)
      # Analyze data distribution for adaptive filtering
      distribution_analysis <- analyze_data_distribution(obj)
      empty_spots <- obj$nFeature_Spatial == 0 & obj$nCount_Spatial == 0
      empty_spots[is.na(empty_spots)] <- TRUE
      if (sum(!empty_spots) > 0) {
        obj <- obj[, !empty_spots]
      }
      if (isTRUE(input$filter_tissue_spots)) {
        if (!"in_tissue" %in% colnames(obj@meta.data)) {
          showNotification(
            "Cannot apply tissue filter: 'in_tissue' missing (no spatial coordinates).",
            type = "warning"
          )
        } else {
          keep_tissue <- obj$in_tissue == 1L
          keep_tissue[is.na(keep_tissue)] <- FALSE
          if (sum(keep_tissue) == 0) {
            removeModal()
            showNotification("No tissue spots detected — check coordinates.", type = "error")
            return(NULL)
          }
          obj <- obj[, keep_tissue, drop = FALSE]
        }
      }
      user_filters_would_pass <- sum(
        obj$nFeature_Spatial >= input$spatial_nFeature_range[1] & 
          obj$nFeature_Spatial <= input$spatial_nFeature_range[2] &
          obj$nCount_Spatial >= input$spatial_nCount_range[1] &
          obj$nCount_Spatial <= input$spatial_nCount_range[2] &
          obj$percent.mt <= input$spatial_mt_max,
        na.rm = TRUE
      )
      
      # Use adaptive filters if user filters would remove too many spots
      if (user_filters_would_pass < 1000 && distribution_analysis$spots_passing > user_filters_would_pass) {
        filters <- distribution_analysis$suggested_filters
        filter_type <- "adaptive"
      } else {
        filters <- list(
          min_features = input$spatial_nFeature_range[1],
          max_features = input$spatial_nFeature_range[2],
          min_counts   = input$spatial_nCount_range[1],
          max_counts   = input$spatial_nCount_range[2],
          max_mt       = input$spatial_mt_max
        )
        filter_type <- "user-defined"
      }
      
      # Apply filters
      keep_cells <- obj$nFeature_Spatial >= filters$min_features & 
        obj$nFeature_Spatial <= filters$max_features &
        obj$nCount_Spatial   >= filters$min_counts &
        obj$nCount_Spatial   <= filters$max_counts &
        obj$percent.mt       <= filters$max_mt
      keep_cells[is.na(keep_cells)] <- FALSE
      
      # Fallback to minimal filters if nothing passes
      if (sum(keep_cells) == 0) {
        keep_cells <- obj$nFeature_Spatial > 0 & obj$nCount_Spatial > 0 & obj$percent.mt < 90
        keep_cells[is.na(keep_cells)] <- FALSE
        filter_type <- "minimal"
      }
      
      obj <- obj[, keep_cells]
      
      update_spatial_object(obj, list(qc_applied = TRUE))
      
      filtered_spots <- ncol(obj)
      removal_rate <- round((1 - filtered_spots / original_spots) * 100, 1)
      
      success_message <- paste0(
        "QC completed! Filter type: ", filter_type, "\n",
        "Remaining: ", format(filtered_spots, big.mark = ","), " spots (", 100 - removal_rate, "%)"
      )
      
      showNotification(success_message, type = "message", duration = 10)
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("QC failed:", e$message), type = "error")
    })
  })
  
  
  # Normalize spatial data WITH DEBUG LOGS
  observeEvent(input$normalize_spatial_data, {
    message("=== NORMALIZE BUTTON CLICKED ===")
    
    if (is.null(spatial_obj())) {
      message("ERROR: spatial_obj() is NULL")
      showNotification("No spatial data loaded!", type = "error")
      return()
    }
    
    tryCatch({
      obj <- spatial_obj()
      if (!"Spatial" %in% Seurat::Assays(obj)) {
        message("ERROR: 'Spatial' assay missing")
        showNotification("Missing 'Spatial' assay in object.", type = "error")
        return()
      }
      
      n_spots <- ncol(obj)
      method  <- (input$norm_method %||% "Normalize")
      nvar_in <- (input$spatial_var_features %||% 3000)
      message(sprintf("Starting normalization | spots=%s | method=%s | nvar=%s",
                      format(n_spots, big.mark=","), method, nvar_in))
      
      showModal(modalDialog(
        title = "Normalizing Data",
        div(
          paste("Normalizing", format(n_spots, big.mark = ","), "spots using", method, "..."),
          br(), tags$small("This may take several minutes")
        ),
        easyClose = FALSE, footer = NULL
      ))
      
      # GC + quiet
      gc()
      old_warn <- getOption("warn"); options(warn = -1)
      
      normalization_success <- FALSE
      method_used <- "none"
      
      if (identical(method, "SCT")) {
        message("=== PATH: SCTransform ===")
        # SCTransform on Spatial assay -> creates 'SCT'
        if (requireNamespace("glmGamPoi", quietly = TRUE)) {
          message("glmGamPoi available — using method='glmGamPoi'")
          obj <- Seurat::SCTransform(
            object = obj, assay = "Spatial", new.assay.name = "SCT",
            method = "glmGamPoi", verbose = FALSE
          )
          method_used <- "SCTransform (glmGamPoi)"
        } else {
          message("glmGamPoi NOT available — falling back to default SCTransform")
          obj <- Seurat::SCTransform(
            object = obj, assay = "Spatial", new.assay.name = "SCT",
            verbose = FALSE
          )
          method_used <- "SCTransform (default)"
        }
        normalization_success <- TRUE
        
      } else {
        message("=== PATH: Classic Normalize + VST variable features ===")
        # Classic path; adapt nfeatures if very large
        nfeatures_use <- if (n_spots > 10000) {
          message(sprintf("Large dataset detected (>10k). Using min(1000, %d) variable features.", nvar_in))
          min(1000L, as.integer(nvar_in))
        } else {
          as.integer(nvar_in)
        }
        
        # Normalize + VFs on Spatial, then ScaleData
        message("Running NormalizeData...")
        obj <- Seurat::NormalizeData(obj, assay = "Spatial",
                                     normalization.method = "LogNormalize",
                                     scale.factor = 10000, verbose = FALSE)
        
        message(sprintf("Running FindVariableFeatures (nfeatures=%d)...", nfeatures_use))
        obj <- Seurat::FindVariableFeatures(obj, assay = "Spatial",
                                            selection.method = "vst",
                                            nfeatures = nfeatures_use, verbose = FALSE)
        
        message("Running ScaleData...")
        obj <- Seurat::ScaleData(obj, assay = "Spatial", verbose = FALSE)
        
        normalization_success <- TRUE
        method_used <- sprintf("Standard normalization (nfeatures=%d)", nfeatures_use)
      }
      
      options(warn = old_warn)
      
      if (!normalization_success) stop("All normalization methods failed")
      
      # ---------- Variable features plot ----------
      message("Preparing Variable Features plot...")
      output$spatial_variable_features <- renderPlot({
        tryCatch({
          if (identical(method, "SCT")) {
            # SCT stores VFs in SCT assay
            if ("SCT" %in% Seurat::Assays(obj) && length(Seurat::VariableFeatures(obj, assay="SCT")) > 0) {
              Seurat::VariableFeaturePlot(obj, assay = "SCT") +
                ggtitle(paste("Variable Features —", method_used))
            } else {
              ggplot() + geom_text(aes(x=0.5, y=0.5, label="No variable features (SCT)"), size=6) + theme_void()
            }
          } else {
            # Classic: VFs in Spatial assay
            if (length(Seurat::VariableFeatures(obj, assay="Spatial")) > 0) {
              Seurat::VariableFeaturePlot(obj, assay = "Spatial") +
                ggtitle(paste("Variable Features —", method_used))
            } else {
              ggplot() + geom_text(aes(x=0.5, y=0.5, label="No variable features (Spatial)"), size=6) + theme_void()
            }
          }
        }, error = function(e) {
          ggplot() + geom_text(aes(x=0.5, y=0.5, label="Normalization completed"), size=6) + theme_void()
        })
      })
      
      # ---------- Update object & notify ----------
      message("Updating spatial object state...")
      update_spatial_object(obj, list(normalized = TRUE, norm_method = method))
      
      gc()
      message(paste("=== NORMALIZATION COMPLETED ===", method_used))
      showNotification(paste("Normalization completed using:", method_used),
                       type = "message", duration = 8)
      removeModal()
      
    }, error = function(e) {
      options(warn = 0); removeModal(); gc()
      message(paste("=== NORMALIZATION ERROR ===", e$message))
      showNotification(paste("Normalization failed:", e$message,
                             "\nTry applying stricter QC filters to reduce dataset size"),
                       type = "error")
    })
  })
  
  
  
  # Apply assay selection
  observeEvent(input$apply_assay_selection, {
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      chosen_assay <- input$spatial_analysis_assay
      
      if (chosen_assay %in% names(obj@assays)) {
        selected_assay(chosen_assay)
        showNotification(paste("Analysis assay set to:", chosen_assay), type = "message", duration = 5)
      } else {
        showNotification("Selected assay not found", type = "error")
      }
      
    }, error = function(e) {
      showNotification(paste("Error setting assay:", e$message), type = "error")
    })
  })
  
  # Update available assays when object changes
  observe({
    req(spatial_obj())
    obj <- spatial_obj()
    assays <- names(obj@assays)
    
    # Update assay selection menu
    updateSelectInput(session, "spatial_analysis_assay", choices = assays, selected = selected_assay())
  })
  
  
  
  # Display current assay information - FIXED VERSION
  output$spatial_assay_info <- renderText({
    req(spatial_obj())
    obj <- spatial_obj()
    
    assay_info <- paste("Available assays:", paste(names(obj@assays), collapse = ", "))
    
    if (sketching_applied()) {
      # Get the actual number of cells in the sketch assay
      sketch_ncells <- ncol(obj[["sketch"]])
      sketch_info <- paste("Sketched data available with", format(sketch_ncells, big.mark = ","), "cells")
      assay_info <- paste(assay_info, "\n", sketch_info)
    }
    
    selected_assay_name <- input$spatial_viz_assay
    if (!is.null(selected_assay_name) && selected_assay_name %in% names(obj@assays)) {
      # Get number of features and cells for the selected assay
      n_features <- nrow(obj[[selected_assay_name]])
      n_cells <- ncol(obj[[selected_assay_name]])
      current_info <- paste("Current assay:", selected_assay_name, "with", 
                            format(n_features, big.mark = ","), "features and",
                            format(n_cells, big.mark = ","), "cells")
      assay_info <- paste(assay_info, "\n", current_info)
    }
    
    return(assay_info)
  })
  
  
  
  
  
  # Set Default Assay function
  observeEvent(input$set_default_assay, {
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      selected_assay <- input$spatial_default_assay
      
      if (selected_assay %in% names(obj@assays)) {
        # Set as default assay in the Seurat object
        DefaultAssay(obj) <- selected_assay
        
        # Update our reactive values
        spatial_obj(obj)
        selected_assay(selected_assay)
        
        showNotification(paste("Default assay set to:", selected_assay), type = "message", duration = 5)
        
        # Update all assay selection menus to reflect the change
        updateSelectInput(session, "spatial_viz_assay", selected = selected_assay)
        updateSelectInput(session, "spatial_assay_select", selected = selected_assay)
        
      } else {
        showNotification("Selected assay not found", type = "error")
      }
      
    }, error = function(e) {
      showNotification(paste("Error setting assay:", e$message), type = "error")
    })
  })
  
  # Update available assays when object changes
  observe({
    req(spatial_obj())
    obj <- spatial_obj()
    assays <- names(obj@assays)
    current_assay <- DefaultAssay(obj)
    
    # Update assay selection menu
    updateSelectInput(session, "spatial_default_assay", choices = assays, selected = current_assay)
    
    # Update current assay info
    output$current_assay_info <- renderText({
      paste("Current:", current_assay)
    })
  })
  
  output$spatial_current_assay_box <- renderInfoBox({
    req(spatial_obj())
    obj <- spatial_obj()
    
    current_assay <- selected_assay()
    if (current_assay %in% names(obj@assays)) {
      n_features <- nrow(obj[[current_assay]])
      color <- if (current_assay == "sketch") "purple" else if (current_assay == "SCT") "green" else "blue"
    } else {
      n_features <- 0
      color <- "red"
    }
    
    infoBox(
      title = "Selected Assay",
      value = current_assay,
      subtitle = paste(format(n_features, big.mark = ","), "features"),
      icon = icon("layer-group"),
      color = color,
      width = 12
    )
  })
  
  
  
  
  
  ############################## PCA & Dimensionality Reduction ##############################
  
  # Run PCA on spatial data
  observeEvent(input$spatial_scale_pca, {
    req(spatial_obj())
    
    if (!processing_states$normalized) {
      showNotification("Please normalize the data first!", type = "warning")
      return()
    }
    
    tryCatch({
      showModal(modalDialog(title = "Running PCA", "Computing principal components...", easyClose = FALSE, footer = NULL))
      
      obj <- spatial_obj()
      assay_for_pca <- selected_assay()
      
      # Determine scaling strategy
      if ("sketch" %in% names(obj@assays) && assay_for_pca == "sketch") {
        message("Using pre-scaled sketch assay for PCA")
      } else if ("SCT" %in% names(obj@assays) && assay_for_pca == "SCT") {
        message("Using SCT assay for PCA")
      } else {
        # Scale data for other assays
        var_features <- VariableFeatures(obj, assay = assay_for_pca)
        if (length(var_features) > 0) {
          message(paste("Scaling", assay_for_pca, "assay with", length(var_features), "variable features"))
          obj <- ScaleData(obj, assay = assay_for_pca, features = var_features, verbose = FALSE)
        } else {
          stop(paste("No variable features found for", assay_for_pca, "assay"))
        }
      }
      
      # Calculate optimal number of PCs  
      n_spots <- ncol(obj)
      n_features <- nrow(obj[[assay_for_pca]])
      max_pcs <- if (n_spots > 100000) 30 else if (n_spots > 50000) 40 else 50
      max_pcs <- min(max_pcs, n_spots - 1, n_features - 1)
      
      message(paste("Running PCA with", max_pcs, "components on", assay_for_pca, "assay"))
      obj <- RunPCA(obj, assay = assay_for_pca, npcs = max_pcs, verbose = FALSE)
      
      output$spatial_elbow_plot <- renderPlot({
        ElbowPlot(obj, ndims = max_pcs) + 
          ggtitle(paste("PCA Elbow Plot -", assay_for_pca, "assay")) +
          labs(subtitle = paste(n_spots, "spots,", max_pcs, "components"))
      })
      
      update_spatial_object(obj, list(pca_computed = TRUE))
      showNotification(paste("PCA completed using", assay_for_pca, "assay"), type = "message")
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("PCA failed:", e$message), type = "error")
    })
  })
  
  
  ############################## Clustering ##############################
  
  # Combined Neighbors + UMAP function
  observeEvent(input$spatial_run_neighbors_umap, {
    message("=== NEIGHBORS + UMAP PIPELINE STARTED ===")
    
    if (is.null(spatial_obj())) {
      showNotification("No spatial data loaded!", type = "error")
      return()
    }
    
    if (!processing_states$pca_computed) {
      showNotification("Please run PCA first!", type = "warning")
      return()
    }
    
    tryCatch({
      obj <- spatial_obj()
      
      # Check PCA availability
      if (!"pca" %in% names(obj@reductions)) {
        stop("PCA reduction not found in object")
      }
      
      # Get parameters from UI inputs
      cluster_dims <- input$spatial_cluster_dims %||% 15
      umap_dims <- input$spatial_umap_dims %||% 30
      
      available_pcs <- ncol(obj@reductions$pca@cell.embeddings)
      max_dims_neighbors <- min(cluster_dims, available_pcs)
      max_dims_umap <- min(umap_dims, available_pcs)
      
      showModal(modalDialog(
        title = "Running Neighbors + UMAP",
        div(
          "Step 1: Computing neighbors...", br(),
          "Step 2: Computing UMAP...", br(),
          tags$small(paste("Neighbors:", max_dims_neighbors, "dims | UMAP:", max_dims_umap, "dims"))
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      
      # Step 1: Find neighbors with reproducible seed
      message(paste("Step 1: Finding neighbors with", max_dims_neighbors, "dimensions..."))
      obj <- findNeighbors_reproducible(
        object = obj,
        dims = 1:max_dims_neighbors,
        seed = 42,
        verbose = FALSE
      )
      message("FindNeighbors computation completed")
      
      # Step 2: Run UMAP with reproducible seed
      message(paste("Step 2: Computing UMAP with", max_dims_umap, "dimensions..."))
      obj <- runUMAP_reproducible(
        object = obj,
        dims = 1:max_dims_umap,
        seed = 42,
        verbose = FALSE
      )
      message("RunUMAP computation completed")
      
      
      # Update object
      spatial_obj(obj)
      message("Spatial object updated with neighbors and UMAP")
      
      # Generate plots
      output$spatial_umap_clusters <- renderPlot({
        DimPlot(obj, reduction = "umap", raster = FALSE) + 
          ggtitle(paste("UMAP (", max_dims_umap, "dims) - Ready for clustering")) +
          theme(plot.title = element_text(hjust = 0.5))
      })
      
      output$spatial_tissue_clusters <- renderPlot({
        if (length(obj@images) > 0) {
          tryCatch({
            SpatialDimPlot(obj) + 
              ggtitle("Spatial View - Ready for clustering") +
              theme(plot.title = element_text(hjust = 0.5))
          }, error = function(e) {
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, label = "Spatial coordinates available"), size = 6) +
              theme_void()
          })
        } else {
          ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5, label = "No spatial coordinates"), size = 6) +
            theme_void()
        }
      })
      
      success_message <- paste0(
        "Neighbors + UMAP completed!\n",
        "Neighbors: ✓ (", max_dims_neighbors, " dims)\n",
        "UMAP: ✓ (", max_dims_umap, " dims)\n",
        "Ready for clustering"
      )
      
      showNotification(success_message, type = "message", duration = 8)
      removeModal()
      
      message("=== NEIGHBORS + UMAP PIPELINE COMPLETED ===")
      
    }, error = function(e) {
      message(paste("=== NEIGHBORS + UMAP ERROR ===", e$message))
      removeModal()
      showNotification(paste("Neighbors + UMAP failed:", e$message), type = "error", duration = 10)
    })
  })
  
  # Find clusters function
  observeEvent(input$spatial_find_clusters, {
    message("=== FIND CLUSTERS BUTTON CLICKED ===")
    
    if (is.null(spatial_obj())) {
      message("ERROR: spatial_obj() is NULL")
      showNotification("No spatial data loaded!", type = "error")
      return()
    }
    
    obj <- spatial_obj()
    
    # Check if neighbors have been computed
    if (is.null(obj@graphs) || length(obj@graphs) == 0) {
      message("ERROR: No neighbor graph found")
      showNotification("Please run 'Run Neighbors + UMAP' first!", type = "warning")
      return()
    }
    
    message("spatial_obj() exists and neighbors computed, proceeding...")
    
    tryCatch({
      # Get and validate parameters
      resolution <- as.numeric(input$spatial_resolution)
      algorithm <- as.integer(input$spatial_cluster_algorithm)
      
      # Validate inputs
      if (is.na(resolution) || resolution <= 0) {
        resolution <- 0.5
        message("Invalid resolution, using default 0.5")
      }
      if (is.na(algorithm) || !algorithm %in% 1:4) {
        algorithm <- 1
        message("Invalid algorithm, using default Louvain")
      }
      
      message(paste("Using resolution:", resolution, "and algorithm:", algorithm))
      
      showModal(modalDialog(
        title = "Finding Clusters",
        div(
          paste("Clustering spots with resolution", resolution, "..."),
          br(),
          tags$small("This may take a few minutes for large datasets")
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      
      message("=== STARTING FIND CLUSTERS ===")
      
      # Get the correct graph name
      graph_name <- names(obj@graphs)[1]
      message(paste("Using graph:", graph_name))
      
      # Find clusters with explicit parameters
      message("Starting FindClusters computation...")
      obj <- FindClusters(obj, 
                          resolution = resolution,
                          algorithm = algorithm,
                          graph.name = graph_name,
                          verbose = FALSE)
      message("FindClusters computation completed")
      
      # Ensure cluster identities are properly formatted
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        obj@meta.data$seurat_clusters <- as.factor(obj@meta.data$seurat_clusters)
        message("Cluster identities converted to factors")
      }
      
      # Update our reactive object and state
      update_spatial_object(obj, list(clustered = TRUE))
      message("Spatial object updated with clusters")
      
      n_clusters <- length(unique(obj$seurat_clusters))
      message(paste("Found", n_clusters, "clusters"))
      
      success_message <- paste0(
        "Clustering completed successfully!\n",
        "Found ", n_clusters, " clusters\n",
        "Resolution: ", resolution, "\n",
        "Algorithm: ", c("Louvain", "Louvain (multilevel)", "SLM")[algorithm]
      )
      
      showNotification(success_message, type = "message", duration = 10)
      removeModal()
      
      message("=== FIND CLUSTERS COMPLETED ===")
      
    }, error = function(e) {
      message(paste("=== FIND CLUSTERS ERROR ===", e$message))
      removeModal()
      
      showNotification(paste("Clustering failed:", e$message), 
                       type = "error", duration = 12)
    })
  })
  
  # Input validation observer
  observe({
    # Validate spatial resolution input
    if (!is.null(input$spatial_resolution)) {
      if (is.na(input$spatial_resolution) || input$spatial_resolution <= 0) {
        updateNumericInput(session, "spatial_resolution", value = 0.5)
        showNotification("Invalid resolution value, reset to 0.5", type = "warning", duration = 3)
      }
    }
    
    # Validate clustering algorithm input
    if (!is.null(input$spatial_cluster_algorithm)) {
      if (is.na(input$spatial_cluster_algorithm) || !input$spatial_cluster_algorithm %in% 1:4) {
        updateSelectInput(session, "spatial_cluster_algorithm", selected = 1)
        showNotification("Invalid algorithm selection, reset to Louvain", type = "warning", duration = 3)
      }
    }
  })
  
  # Generate clustering plots after successful clustering
  observe({
    req(spatial_obj())
    
    if (processing_states$clustered) {
      obj <- spatial_obj()
      
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        n_clusters <- length(unique(obj$seurat_clusters))
        
        # Generate UMAP clustering plot
        output$spatial_umap_clusters <- renderPlot({
          if ("umap" %in% names(obj@reductions)) {
            p <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", 
                         label = input$spatial_plot_labels, pt.size = 0.8, 
                         label.size = 4, raster = FALSE)
            if (!input$spatial_plot_labels) p <- p + NoLegend()
            p + ggtitle(paste("UMAP Clusters (n =", n_clusters, ")")) +
              theme(plot.title = element_text(hjust = 0.5, size = 14)) +
              guides(color = guide_legend(override.aes = list(size = 3)))
          } else {
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, label = "UMAP not available\nRun UMAP first"), 
                        size = 6, color = "red") +
              theme_void()
          }
        })
        
        # Generate spatial clustering plot
        output$spatial_tissue_clusters <- renderPlot({
          if (length(obj@images) > 0) {
            tryCatch({
              SpatialDimPlot(obj, group.by = "seurat_clusters", 
                             label = input$spatial_plot_labels,
                             pt.size.factor = 1.2, alpha = c(0.3, 1),
                             label.size = 3) +
                ggtitle(paste("Spatial Clusters (n =", n_clusters, ")")) +
                theme(plot.title = element_text(hjust = 0.5, size = 14),
                      legend.text = element_text(size = 10)) +
                guides(fill = guide_legend(override.aes = list(size = 3)))
            }, error = function(e) {
              message(paste("Spatial plot error:", e$message))
              ggplot() + 
                geom_text(aes(x = 0.5, y = 0.5, 
                              label = paste("Spatial plot unavailable\n", n_clusters, "clusters found")), 
                          size = 6) +
                theme_void()
            })
          } else {
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, 
                            label = paste("No spatial coordinates\n", n_clusters, "clusters found")), 
                        size = 6) +
              theme_void()
          }
        })
        
        # Generate cluster statistics table
        output$spatial_cluster_stats <- renderDT({
          tryCatch({
            cluster_stats <- obj@meta.data %>%
              group_by(seurat_clusters) %>%
              summarise(
                n_spots = n(),
                mean_genes = round(mean(nFeature_Spatial, na.rm = TRUE), 1),
                mean_umis = round(mean(nCount_Spatial, na.rm = TRUE), 1),
                mean_mt = round(mean(percent.mt, na.rm = TRUE), 2),
                .groups = 'drop'
              )
            
            datatable(cluster_stats, 
                      options = list(pageLength = 15, scrollX = TRUE, dom = 'tip'),
                      rownames = FALSE) %>%
              formatStyle(columns = 1:ncol(cluster_stats), fontSize = '12px')
          }, error = function(e) {
            message(paste("Cluster stats error:", e$message))
            datatable(data.frame(Message = "Statistics unavailable"))
          })
        })
      }
    } else {
      # Show placeholder plots when not clustered
      output$spatial_umap_clusters <- renderPlot({
        if (!is.null(spatial_obj())) {
          obj <- spatial_obj()
          if ("umap" %in% names(obj@reductions)) {
            DimPlot(obj, reduction = "umap", raster = FALSE) + 
              ggtitle("UMAP (Ready for clustering)") +
              theme(plot.title = element_text(hjust = 0.5))
          } else {
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, label = "Run neighbors & UMAP first"), size = 6) +
              theme_void()
          }
        } else {
          ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5, label = "No data loaded"), size = 6) +
            theme_void()
        }
      })
      
      output$spatial_tissue_clusters <- renderPlot({
        if (!is.null(spatial_obj())) {
          obj <- spatial_obj()
          if (length(obj@images) > 0) {
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, label = "Ready for clustering"), size = 6) +
              theme_void()
          } else {
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, label = "No spatial coordinates"), size = 6) +
              theme_void()
          }
        } else {
          ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5, label = "No data loaded"), size = 6) +
            theme_void()
        }
      })
    }
  })
  
  # Status displays for processing state
  output$spatial_processing_status <- renderInfoBox({
    status <- get_processing_status()
    
    if (grepl("Clustered", status)) {
      color <- "green"
      icon_name <- "check-circle"
    } else if (grepl("PCA", status)) {
      color <- "blue"
      icon_name <- "chart-line"
    } else if (grepl("Normalized", status)) {
      color <- "yellow"
      icon_name <- "balance-scale"
    } else {
      color <- "red"
      icon_name <- "exclamation-triangle"
    }
    
    infoBox(
      title = "Processing Status",
      value = status,
      icon = icon(icon_name),
      color = color,
      width = 12
    )
  })
  
  # Individual status outputs
  output$spatial_neighbors_status <- renderText({
    if (!is.null(spatial_obj())) {
      obj <- spatial_obj()
      if (!is.null(obj@graphs) && length(obj@graphs) > 0) {
        "✓ Neighbors computed"
      } else {
        "⚪ Neighbors not computed"
      }
    } else {
      "⚪ No data loaded"
    }
  })
  
  output$spatial_clusters_status <- renderText({
    if (!is.null(spatial_obj())) {
      obj <- spatial_obj()
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        n_clusters <- length(unique(obj$seurat_clusters))
        paste("✓ Clusters found:", n_clusters)
      } else {
        "⚪ Clusters not computed"
      }
    } else {
      "⚪ No data loaded"
    }
  })
  
  output$spatial_pca_status <- renderText({
    if (processing_states$pca_computed) {
      "✓ PCA completed"
    } else {
      "⚪ PCA not computed"
    }
  })
  
  output$spatial_umap_status <- renderText({
    if (!is.null(spatial_obj())) {
      obj <- spatial_obj()
      if ("umap" %in% names(obj@reductions)) {
        "✓ UMAP completed"
      } else {
        "⚪ UMAP not computed"
      }
    } else {
      "⚪ No data loaded"
    }
  })
  
  
  # Add this helper function to validate spatial data integrity
  validate_spatial_clustering_data <- function(obj) {
    issues <- list()
    
    # Check cluster assignments
    if ("seurat_clusters" %in% colnames(obj@meta.data)) {
      cluster_data <- obj@meta.data$seurat_clusters
      na_count <- sum(is.na(cluster_data))
      inf_count <- sum(is.infinite(cluster_data))
      
      if (na_count > 0) issues <- append(issues, paste(na_count, "NA cluster assignments"))
      if (inf_count > 0) issues <- append(issues, paste(inf_count, "infinite cluster assignments"))
    }
    
    # Check spatial coordinates
    if (length(obj@images) > 0) {
      tryCatch({
        coords <- GetTissueCoordinates(obj@images[[1]])
        if (!is.null(coords)) {
          coord_na <- sum(is.na(coords))
          coord_inf <- sum(is.infinite(as.matrix(coords)))
          
          if (coord_na > 0) issues <- append(issues, paste(coord_na, "NA coordinates"))
          if (coord_inf > 0) issues <- append(issues, paste(coord_inf, "infinite coordinates"))
        }
      }, error = function(e) {
        issues <- append(issues, "Coordinate extraction failed")
      })
    }
    
    return(list(valid = length(issues) == 0, issues = issues))
  }
  
  # Helper function to check if spatial coordinates are available
  check_spatial_coordinates <- function(obj) {
    if (length(obj@images) == 0) {
      return(list(available = FALSE, message = "No spatial images found"))
    }
    
    tryCatch({
      coords <- GetTissueCoordinates(obj@images[[1]])
      if (is.null(coords) || nrow(coords) == 0) {
        return(list(available = FALSE, message = "No tissue coordinates found"))
      }
      
      # Check if coordinates match the cells
      cell_names <- colnames(obj)
      coord_names <- rownames(coords)
      
      overlap <- intersect(cell_names, coord_names)
      if (length(overlap) == 0) {
        return(list(available = FALSE, message = "No matching coordinates for cells"))
      }
      
      return(list(available = TRUE, 
                  message = paste("Coordinates available for", length(overlap), "spots"),
                  coords = coords))
      
    }, error = function(e) {
      return(list(available = FALSE, message = paste("Error accessing coordinates:", e$message)))
    })
  }
  
  # Update the spatial clustering plot to handle missing coordinates better
  observe({
    req(spatial_obj())
    
    if (processing_states$clustered) {
      obj <- spatial_obj()
      
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        n_clusters <- length(unique(obj$seurat_clusters))
        
        # Generate spatial clustering plot with better error handling
        output$spatial_tissue_clusters <- renderPlot({
          coord_check <- check_spatial_coordinates(obj)
          
          if (!coord_check$available) {
            # Show informative message
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, 
                            label = paste("Spatial visualization unavailable\n", 
                                          coord_check$message, "\n",
                                          n_clusters, "clusters found")), 
                        size = 5, color = "orange") +
              theme_void() +
              labs(title = "Use UMAP view for cluster visualization")
          } else {
            tryCatch({
              SpatialDimPlot(obj, group.by = "seurat_clusters", 
                             label = input$spatial_plot_labels,
                             pt.size.factor = 1.2, alpha = c(0.3, 1),
                             label.size = 3) +
                ggtitle(paste("Spatial Clusters (n =", n_clusters, ")")) +
                theme(plot.title = element_text(hjust = 0.5, size = 14),
                      legend.text = element_text(size = 10)) +
                guides(fill = guide_legend(override.aes = list(size = 3)))
            }, error = function(e) {
              message(paste("Spatial plot error:", e$message))
              ggplot() + 
                geom_text(aes(x = 0.5, y = 0.5, 
                              label = paste("Spatial plot error\n", 
                                            "Try using UMAP view\n",
                                            n_clusters, "clusters found")), 
                          size = 5, color = "red") +
                theme_void()
            })
          }
        })
      }
    }
  })
  
  ##########################Interactive Visualization##############################
  
  # Update cluster choices
  observe({
    req(spatial_obj())
    obj <- spatial_obj()
    
    if (processing_states$clustered) {
      cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) {
        "annotated_clusters"
      } else {
        "seurat_clusters"
      }
      
      cluster_choices <- c("None", as.character(unique(obj@meta.data[[cluster_col]])))
      
      updateSelectInput(session, "interactive_cluster_highlight", 
                        choices = cluster_choices,
                        selected = "None")
    }
  })
  
  # Reactive function to generate colors based on palette selection
  generate_palette_colors <- reactive({
    req(spatial_obj())
    obj <- spatial_obj()
    
    cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) {
      "annotated_clusters"
    } else {
      "seurat_clusters"
    }
    
    clusters <- sort(unique(obj@meta.data[[cluster_col]]))
    n_clusters <- length(clusters)
    palette_choice <- input$interactive_color_palette %||% "default"
    
    colors <- switch(palette_choice,
                     "default" = {
                       # Use stored colors or Seurat default
                       if (!is.null(obj@misc$cluster_colors)) {
                         obj@misc$cluster_colors
                       } else {
                         scales::hue_pal()(n_clusters)
                       }
                     },
                     "polychrome" = {
                       # High contrast palette
                       if (n_clusters <= 36) {
                         Polychrome::palette36.colors(n_clusters)
                       } else {
                         Polychrome::createPalette(n_clusters, seedcolors = c("#ff0000", "#00ff00", "#0000ff"))
                       }
                     },
                     "set3" = {
                       colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_clusters)
                     },
                     "paired" = {
                       colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_clusters)
                     },
                     "dark2" = {
                       colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_clusters)
                     },
                     scales::hue_pal()(n_clusters)  # Fallback
    )
    
    names(colors) <- clusters
    return(colors)
  })
  
  # Dynamic UI for plot with variable height
  output$interactive_spatial_plot_ui <- renderUI({
    plot_height <- input$interactive_plot_height %||% 700
    plotOutput("interactive_spatial_plot", height = paste0(plot_height, "px"))
  })
  
  # Main spatial plot with highlight functionality
  output$interactive_spatial_plot <- renderPlot({
    req(spatial_obj())
    obj <- spatial_obj()
    
    cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) {
      "annotated_clusters" 
    } else {
      "seurat_clusters"
    }
    
    if (length(obj@images) == 0) {
      return(ggplot() + geom_text(aes(0.5, 0.5, label = "No spatial image"), size = 6) + theme_void())
    }
    
    # Get alpha adjustment
    plot_height <- input$interactive_plot_height %||% 700
    base_alpha <- input$interactive_image_alpha %||% 0.7
    height_factor <- min(plot_height / 700, 2)
    adjusted_alpha <- min(base_alpha * height_factor, 1)
    
    # Get colors from palette
    cluster_colors <- generate_palette_colors()
    
    # Check if highlighting a specific cluster
    if (!is.null(input$interactive_cluster_highlight) && 
        input$interactive_cluster_highlight != "None") {
      
      # Create highlight version
      obj_highlight <- obj
      obj_highlight@meta.data$highlight_group <- ifelse(
        obj_highlight@meta.data[[cluster_col]] == input$interactive_cluster_highlight,
        as.character(obj_highlight@meta.data[[cluster_col]]),
        "Other"
      )
      
      # Colors: highlighted cluster keeps its color, others are gray
      selected_cluster <- input$interactive_cluster_highlight
      highlight_colors <- c("Other" = "#CCCCCC")
      highlight_colors[selected_cluster] <- cluster_colors[selected_cluster]
      
      p <- Seurat::SpatialDimPlot(
        obj_highlight,
        group.by       = "highlight_group",
        label          = input$interactive_show_labels,
        label.size     = 4,
        pt.size.factor = 1.6,
        alpha          = c(adjusted_alpha, 0.9),
        cols           = highlight_colors
      ) +
        labs(title = NULL, fill = NULL)
      
    } else {
      # Normal plot with all clusters
      p <- Seurat::SpatialDimPlot(
        obj,
        group.by       = cluster_col,
        label          = input$interactive_show_labels,
        label.size     = 4,
        pt.size.factor = 1.6,
        alpha          = c(adjusted_alpha, 1),
        cols           = cluster_colors
      ) +
        labs(title = NULL, fill = NULL)
    }
    
    if(isTRUE(input$interactive_remove_legend)) {
      p <- p + NoLegend()
    }
    
    p
  })
  
  interactive_spatial_viz_plot <- reactiveVal(NULL)
  
  output$download_interactive_spatial <- downloadHandler(
    filename = function() {
      paste0("spatial_tissue_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(spatial_obj())
      obj <- spatial_obj()
      
      cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) {
        "annotated_clusters" 
      } else {
        "seurat_clusters"
      }
      
      image_alpha <- input$interactive_image_alpha %||% 0.7
      cluster_colors <- generate_palette_colors()
      
      p <- Seurat::SpatialDimPlot(
        obj,
        group.by       = cluster_col,
        label          = input$interactive_show_labels,
        label.size     = 4,
        pt.size.factor = 1.6,
        alpha          = c(1, 1),  # Full opacity for download
        cols           = cluster_colors
      ) +
        labs(title = NULL, fill = NULL)
      
      if(isTRUE(input$interactive_remove_legend)) {
        p <- p + NoLegend()
      }
      
      ggsave(file, plot = p, width = 12, height = 10, dpi = 300, bg = "white")
    }
  )
  ############## Gene selection plot #####################
  # Update gene choices when spatial data is loaded
  observe({
    req(spatial_obj())
    
    if (processing_states$normalized || !is.null(spatial_obj())) {
      obj <- spatial_obj()
      gene_choices <- rownames(obj)
      
      # Update gene selection picker
      updatePickerInput(session, "spatial_genes_select",
                        choices = gene_choices,
                        selected = NULL)
      
      # Update cluster choices if clusters exist
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        cluster_choices <- unique(obj@meta.data$seurat_clusters)
        cluster_choices <- sort(as.character(cluster_choices))
        
        updateSelectInput(session, "spatial_cluster_order",
                          choices = cluster_choices,
                          selected = cluster_choices)
      }
    }
  })
  
  observe({
    req(spatial_obj())
    
    obj <- spatial_obj()
    
    # Update gene choices for all assays
    all_genes <- character(0)
    current_assay <- input$spatial_viz_assay
    
    if (!is.null(current_assay) && current_assay %in% names(obj@assays)) {
      all_genes <- rownames(obj[[current_assay]])
    } else {
      # Use default assay genes
      all_genes <- rownames(obj[[DefaultAssay(obj)]])
    }
    
    message(paste("Updating gene choices. Available genes:", length(all_genes)))
    
    # Update gene selection picker
    updatePickerInput(session, "spatial_genes_select",
                      choices = all_genes,
                      selected = NULL)
    
    # Update assay choices to include all available assays
    available_assays_list <- names(obj@assays)
    current_default <- DefaultAssay(obj)
    
    updateSelectInput(session, "spatial_viz_assay", 
                      choices = available_assays_list, 
                      selected = current_default)
    
    message(paste("Available assays:", paste(available_assays_list, collapse = ", ")))
    message(paste("Default assay:", current_default))
  })
  
  
  
  # Sync between text inputs and picker selection - SIMPLIFIED AND SAFE
  observeEvent(input$spatial_genes_select, {
    # Only update if the picker has actual selections
    if (!is.null(input$spatial_genes_select) && length(input$spatial_genes_select) > 0) {
      
      # Simple replacement - less complex, more stable
      gene_text <- paste(input$spatial_genes_select, collapse = ", ")
      
      # Update text inputs safely
      updateTextInput(session, "spatial_genes_text", value = gene_text)
      updateTextInput(session, "spatial_violin_genes_text", value = gene_text)
      updateTextInput(session, "spatial_dotplot_genes_text", value = gene_text)
    }
  }, ignoreInit = TRUE)
  
  # Reactive values to store generated plots
  spatial_feature_plot_reactive <- reactiveVal(NULL)
  spatial_violin_plot_reactive <- reactiveVal(NULL)
  spatial_dot_plot_reactive <- reactiveVal(NULL)
  
  # Generate spatial feature plots - AVEC JUSTE LUMINOSITE AJOUTEE
  observeEvent(input$show_spatial_features, {
    genes_from_text <- character(0)
    if (!is.null(input$spatial_genes_text) && nchar(trimws(input$spatial_genes_text)) > 0) {
      genes_from_text <- trimws(strsplit(input$spatial_genes_text, ",")[[1]])
      genes_from_text <- genes_from_text[genes_from_text != ""]
    }
    
    if (length(genes_from_text) == 0) {
      showNotification("Please enter at least one gene", type = "warning")
      return()
    }
    
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      selected_genes <- genes_from_text
      current_assay <- input$spatial_viz_assay %||% "Spatial"
      
      message(paste("=== SPATIAL FEATURE PLOT ==="))
      message(paste("Using assay:", current_assay))
      message(paste("Selected genes:", paste(selected_genes, collapse = ", ")))
      
      # Set the assay temporarily
      original_assay <- DefaultAssay(obj)
      DefaultAssay(obj) <- current_assay
      
      # Validate genes exist in current assay
      available_genes <- rownames(obj[[current_assay]])
      valid_genes <- selected_genes[selected_genes %in% available_genes]
      missing_genes <- selected_genes[!selected_genes %in% available_genes]
      
      if (length(missing_genes) > 0) {
        showNotification(paste("Missing genes:", paste(missing_genes, collapse = ", ")), 
                         type = "warning", duration = 8)
      }
      
      if (length(valid_genes) == 0) {
        showNotification("None of the selected genes found in current assay", type = "error")
        DefaultAssay(obj) <- original_assay
        return()
      }
      
      showModal(modalDialog(
        title = "Generating Spatial Feature Plots",
        paste("Creating plots for", length(valid_genes), "genes using", current_assay, "assay..."),
        easyClose = FALSE, footer = NULL
      ))
      
      # Generate the spatial feature plot with dynamic sizing
      output$spatial_feature_plots <- renderPlot({
        if (length(obj@images) == 0) {
          p <- ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5, label = "No spatial coordinates available"), size = 8) + 
            theme_void()
          
          # Stocker même le plot vide
          spatial_feature_plot_reactive(p)
          return(p)
          
        } else {
          tryCatch({
            n_genes <- length(valid_genes)
            
            # Get image size factor from slider
            image_size_factor <- input$spatial_image_size %||% 1
            
            # Get brightness factor
            brightness_factor <- input$spatial_brightness %||% 1
            
            # Determine number of columns
            if (input$spatial_combine_plots && n_genes > 1) {
              ncol_val <- input$spatial_plot_ncol %||% min(3, n_genes)
            } else {
              ncol_val <- 1
            }
            
            message(paste("Creating spatial plot with", n_genes, "genes in", ncol_val, "columns"))
            message(paste("Image size factor:", image_size_factor))
            message(paste("Brightness factor:", brightness_factor))
            
            # Create plots based on number of genes
            if (n_genes == 1) {
              # Single gene plot
              p <- SpatialFeaturePlot(
                object = obj, 
                features = valid_genes[1],
                crop = input$spatial_crop %||% TRUE,
                slot = "data",
                pt.size.factor = (input$spatial_pt_size %||% 1.6) * image_size_factor,
                alpha = c(0.1, (input$spatial_alpha %||% 1) * brightness_factor),
                image.alpha = 1,
                stroke = 0.25 * image_size_factor
              ) +
                ggtitle(paste("Spatial Expression:", valid_genes[1], "(", current_assay, ")")) +
                theme(plot.title = element_text(hjust = 0.5, size = 14 * image_size_factor))
              
            } else if (input$spatial_combine_plots) {
              # Multi-gene combined plot
              p <- SpatialFeaturePlot(
                object = obj, 
                features = valid_genes,
                crop = input$spatial_crop %||% TRUE,
                slot = "data",
                ncol = ncol_val,
                combine = TRUE,
                pt.size.factor = (input$spatial_pt_size %||% 1.6) * image_size_factor,
                alpha = c(0.1, (input$spatial_alpha %||% 1) * brightness_factor),
                image.alpha = 0.8,
                stroke = 0.25 * image_size_factor
              )
              
            } else {
              # Show first gene only if not combining
              p <- SpatialFeaturePlot(
                object = obj, 
                features = valid_genes[1],
                crop = input$spatial_crop %||% TRUE,
                slot = "data",
                pt.size.factor = (input$spatial_pt_size %||% 1.6) * image_size_factor,
                alpha = c(0.1, (input$spatial_alpha %||% 1) * brightness_factor),
                image.alpha = 1,
                stroke = 0.25 * image_size_factor
              ) +
                ggtitle(paste("Spatial Expression:", valid_genes[1], "(", current_assay, ")")) +
                theme(plot.title = element_text(hjust = 0.5, size = 14 * image_size_factor))
            }
            
            # IMPORTANT: STOCKER LE PLOT DANS LA REACTIVE VALUE ICI
            spatial_feature_plot_reactive(p)
            message("Feature plot stored in reactive value")
            
            return(p)
            
          }, error = function(e) {
            message(paste("Spatial plot error:", e$message))
            error_plot <- ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, label = paste("Error:", substr(e$message, 1, 50))), 
                        size = 4, color = "red") +
              theme_void()
            
            # Stocker même le plot d'erreur
            spatial_feature_plot_reactive(error_plot)
            
            return(error_plot)
          })
        }
      }, 
      height = function() {
        size_factor <- input$spatial_image_size %||% 1
        n_genes <- length(valid_genes)
        
        # Marge pour titre + légende + padding
        title_legend_margin <- 150
        
        if (n_genes == 1) {
          # Hauteur du plot pur + marge pour titre/légende
          plot_height <- 600 * size_factor
          total_height <- plot_height + title_legend_margin
          return(max(550, total_height))  # Minimum 550px
        } else {
          ncol_val <- input$spatial_plot_ncol %||% min(3, n_genes)
          nrow_val <- ceiling(n_genes / ncol_val)
          plot_height <- 400 * nrow_val * size_factor
          total_height <- plot_height + title_legend_margin
          return(max(650, min(1500, total_height)))
        }
      })
      
      # Store plot parameters for download
      spatial_plots$feature_plots <- list(
        genes = valid_genes, 
        assay = current_assay,
        parameters = list(
          pt_size = (input$spatial_pt_size %||% 1.6) * (input$spatial_image_size %||% 1), 
          alpha = (input$spatial_alpha %||% 1) * (input$spatial_brightness %||% 1),
          crop = input$spatial_crop %||% TRUE, 
          combine = input$spatial_combine_plots,
          image_size = input$spatial_image_size %||% 1
        )
      )
      
      # Restore original assay
      DefaultAssay(obj) <- original_assay
      
      removeModal()
      showNotification(paste("Spatial plots generated for", length(valid_genes), "genes"), type = "message")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error generating spatial plots:", e$message), type = "error")
    })
  })
  
  # Generate spatial dot plot by clusters
  observeEvent(input$show_spatial_dotplot, {
    genes_from_text <- character(0)
    if (!is.null(input$spatial_dotplot_genes_text) && nchar(trimws(input$spatial_dotplot_genes_text)) > 0) {
      genes_from_text <- trimws(strsplit(input$spatial_dotplot_genes_text, ",")[[1]])
      genes_from_text <- genes_from_text[genes_from_text != ""]
    }
    
    if (length(genes_from_text) == 0) {
      showNotification("Please enter at least one gene for dot plot", type = "warning")
      return()
    }
    
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      
      # Use annotated clusters if available
      cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) {
        "annotated_clusters"
      } else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        "seurat_clusters"
      } else {
        showNotification("Please run clustering first", type = "warning")
        return()
      }
      
      selected_genes <- genes_from_text
      selected_assay <- input$spatial_viz_assay %||% "Spatial"
      
      # Validate genes exist in assay
      available_genes <- rownames(obj[[selected_assay]])
      valid_genes <- selected_genes[selected_genes %in% available_genes]
      
      if (length(valid_genes) == 0) {
        showNotification("Selected genes not found in current assay", type = "error")
        return()
      }
      
      # Set assay temporarily
      original_assay <- DefaultAssay(obj)
      DefaultAssay(obj) <- selected_assay
      
      # Generate plot
      # Generate plot
      p <- DotPlot(obj, 
                   features = valid_genes, 
                   group.by = cluster_col,
                   scale = input$spatial_dot_scale %||% TRUE, 
                   scale.by = input$spatial_dot_scale_by %||% "radius",
                   dot.min = (input$spatial_dot_min %||% 0)/100, 
                   dot.scale = input$spatial_dot_max %||% 6,
                   cluster.idents = input$spatial_dot_cluster_idents %||% FALSE, 
                   cols = input$spatial_dot_colors %||% "RdYlBu") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
              plot.title = element_text(hjust = 0.5))
      
      # Inverser les axes SI la checkbox est cochée
      if (input$spatial_dot_cluster_idents) {
        p <- p + coord_flip() + labs(x = "Clusters", y = "Genes")
      } else {
        p <- p + labs(x = "Genes", y = "Clusters")
      }
      
      # Store and display
      spatial_dot_plot_reactive(p)
      output$spatial_dotplot <- renderPlot({ 
        p 
      })      
      # Restore original assay
      DefaultAssay(obj) <- original_assay
      
      showNotification("Dot plot generated", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error generating dot plot:", e$message), type = "error")
    })
  })
  
  # Generate spatial violin plot by clusters - PROPERLY FIXED
  observeEvent(input$show_spatial_violin, {
    genes_from_text <- character(0)
    if (!is.null(input$spatial_violin_genes_text) && nchar(trimws(input$spatial_violin_genes_text)) > 0) {
      genes_from_text <- trimws(strsplit(input$spatial_violin_genes_text, ",")[[1]])
      genes_from_text <- genes_from_text[genes_from_text != ""]
    }
    
    if (length(genes_from_text) == 0) {
      showNotification("Please enter at least one gene for violin plot", type = "warning")
      return()
    }
    
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      
      # Use annotated clusters if available
      cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) {
        "annotated_clusters"
      } else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        "seurat_clusters"
      } else {
        showNotification("Please run clustering first", type = "warning")
        return()
      }
      
      selected_genes <- genes_from_text
      selected_assay <- input$spatial_viz_assay %||% "Spatial"
      
      # Validate genes exist in assay
      available_genes <- rownames(obj[[selected_assay]])
      valid_genes <- selected_genes[selected_genes %in% available_genes]
      missing_genes <- selected_genes[!selected_genes %in% available_genes]
      
      if (length(missing_genes) > 0) {
        showNotification(paste("Missing genes:", paste(missing_genes, collapse = ", ")), 
                         type = "warning", duration = 8)
      }
      
      if (length(valid_genes) == 0) {
        showNotification("Selected genes not found in current assay", type = "error")
        return()
      }
      
      # Set assay temporarily
      original_assay <- DefaultAssay(obj)
      DefaultAssay(obj) <- selected_assay
      
      message(paste("Creating violin plot with genes:", paste(valid_genes, collapse = ", ")))
      message(paste("Using cluster column:", cluster_col))
      message(paste("Using assay:", selected_assay))
      
      # Create violin plot with CORRECT parameters
      p <- VlnPlot(
        object = obj,
        features = valid_genes,
        assay = selected_assay,          # Specify assay explicitly
        group.by = cluster_col,          # Group by cluster column
        pt.size = if (input$spatial_violin_points) input$spatial_violin_pt_size else 0,
        log = input$spatial_violin_log,
        ncol = if (length(valid_genes) > 1) min(3, length(valid_genes)) else NULL,
        combine = TRUE,
        fill.by = "ident",               # Color by identity (clusters)
        layer = "data",                  # Use processed data layer
        flip = FALSE,                    # Keep standard orientation
        same.y.lims = FALSE,            # Allow different y-axis for each gene
        raster = FALSE                   # Don't rasterize for better quality
      )
      
      # Apply custom theming and labels
      if (length(valid_genes) == 1) {
        # Single gene - enhance appearance
        p <- p + 
          ggtitle(paste("Expression of", valid_genes[1], "across", cluster_col)) +
          labs(
            x = paste("Clusters (", cluster_col, ")"),
            y = paste0("Expression Level", if(input$spatial_violin_log) " (log)" else ""),
            subtitle = paste("Assay:", selected_assay)
          ) +
          theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            axis.text.y = element_text(size = 12),
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 12),
            legend.position = "right"
          )
        
      } else {
        # Multiple genes - use patchwork for better control
        p <- p & theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 12, face = "bold")  # Gene names
        )
        
        # Add overall annotation
        p <- p + plot_annotation(
          title = paste("Gene Expression across", cluster_col),
          subtitle = paste("Assay:", selected_assay, "| Genes:", paste(valid_genes, collapse = ", ")),
          theme = theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12)
          )
        )
      }
      
      # Apply cluster order if specified
      cluster_order <- input$spatial_cluster_order
      if (!is.null(cluster_order) && length(cluster_order) > 0) {
        # Get current clusters and apply order
        current_clusters <- levels(as.factor(obj@meta.data[[cluster_col]]))
        valid_order <- cluster_order[cluster_order %in% current_clusters]
        
        if (length(valid_order) > 0) {
          # Add missing clusters at the end
          missing_clusters <- setdiff(current_clusters, valid_order)
          final_order <- c(valid_order, missing_clusters)
          
          # Apply the order
          p <- p + scale_x_discrete(limits = final_order)
        }
      }
      
      # Store in reactive value for reuse
      spatial_violin_plot_reactive(p)
      
      output$spatial_violin_plot <- renderPlot({
        spatial_violin_plot_reactive()
      }, height = function() {
        # Compter le nombre de gènes
        if (!is.null(input$spatial_violin_genes_text) && nchar(trimws(input$spatial_violin_genes_text)) > 0) {
          genes <- trimws(strsplit(input$spatial_violin_genes_text, ",")[[1]])
          n_genes <- length(genes[genes != ""])
        } else {
          n_genes <- 1
        }
        
        # Calculer le nombre de lignes (ncol = 3 par défaut pour Seurat)
        ncol_val <- 3
        n_rows <- ceiling(n_genes / ncol_val)
        
        # 350px par ligne de plots
        return(350 * n_rows)
      }) 
      
      # Restore original assay
      DefaultAssay(obj) <- original_assay
      
      
    }, error = function(e) {
      message(paste("Violin plot error:", e$message))
      showNotification(paste("Error generating violin plot:", e$message), type = "error")
    })
  })
  
  # Download handlers using stored plots
  output$download_spatial_features <- downloadHandler(
    filename = function() {
      paste0("spatial_feature_", Sys.Date(), ".", input$spatial_plot_format %||% "png")
    },
    content = function(file) {
      p <- spatial_feature_plot_reactive()
      
      if (!is.null(p)) {
        plot_width <- input$spatial_plot_width %||% 10
        plot_height <- plot_width * 0.8
        
        # Map format to proper device name
        device_name <- switch(tolower(input$spatial_plot_format %||% "png"),
                              "pdf" = "pdf",
                              "png" = "png", 
                              "tiff" = "tiff",
                              "jpeg" = "jpeg",
                              "svg" = "svg",
                              "png")  # default
        
        ggsave(file, plot = p, 
               width = plot_width, 
               height = plot_height,
               dpi = input$spatial_export_dpi %||% 300, 
               device = device_name)
      }
    }
  )
  output$download_spatial_violin <- downloadHandler(
    filename = function() {
      paste0("spatial_violin_", Sys.Date(), ".", input$spatial_plot_format %||% "png")
    },
    content = function(file) {
      p <- spatial_violin_plot_reactive()
      
      if (!is.null(p)) {
        plot_width <- input$spatial_plot_width %||% 10
        plot_height <- plot_width * 0.6
        
        ggsave(file, plot = p, 
               width = plot_width, 
               height = plot_height,
               dpi = input$spatial_export_dpi %||% 300,
               device = input$spatial_plot_format %||% "png")
      }
    }
  )
  
  output$download_spatial_dotplot <- downloadHandler(
    filename = function() {
      paste0("spatial_dotplot_", Sys.Date(), ".", input$spatial_plot_format %||% "png")
    },
    content = function(file) {
      p <- spatial_dot_plot_reactive()
      
      if (!is.null(p)) {
        plot_width <- input$spatial_plot_width %||% 10
        plot_height <- plot_width * 0.8
        
        ggsave(file, plot = p, 
               width = plot_width, 
               height = plot_height,
               dpi = input$spatial_export_dpi %||% 300,
               device = input$spatial_plot_format %||% "png")
      }
    }
  )
  
  # Save spatial object
  output$save_spatial_object <- downloadHandler(
    filename = function() {
      paste0("spatial_object_", Sys.Date(), ".rds")
    },
    content = function(file) {
      showModal(modalDialog(
        title = "Saving Spatial Object",
        "Saving your spatial data object...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      req(spatial_obj())
      saveRDS(spatial_obj(), file)
      
      removeModal()
      showNotification("Spatial object saved successfully!", type = "message")
    }
  )
  
  
  
  ############################## Cluster Annotation ##############################
  
  # Reactive value to store cluster names
  spatial_cluster_names <- reactiveVal(NULL)
  
  # Initialize cluster names when clusters are available
  observe({
    req(spatial_obj())
    
    if (processing_states$clustered && "seurat_clusters" %in% colnames(spatial_obj()@meta.data)) {
      obj <- spatial_obj()
      unique_clusters <- sort(unique(as.character(obj$seurat_clusters)))
      
      # Initialize with default names if not already set
      if (is.null(spatial_cluster_names())) {
        default_names <- setNames(paste0("Cluster_", unique_clusters), unique_clusters)
        spatial_cluster_names(default_names)
      }
    }
  })
  
  # Generate dynamic input fields for cluster renaming
  output$spatial_cluster_rename_inputs <- renderUI({
    req(spatial_obj())
    req(processing_states$clustered)
    
    obj <- spatial_obj()
    if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
      return(p("No clusters found. Please run clustering first.", style = "color: #666;"))
    }
    
    # Get original cluster numbers
    unique_clusters <- sort(unique(as.character(obj$seurat_clusters)))
    
    # Get current names (either annotated or original)
    current_names <- spatial_cluster_names()
    
    # Check if we have annotated clusters to show current state
    has_annotated <- "annotated_clusters" %in% colnames(obj@meta.data)
    
    # Create input fields for each cluster
    input_list <- lapply(unique_clusters, function(cluster) {
      
      # Determine current display name
      if (has_annotated) {
        # Get the annotated name for this cluster number
        cluster_mask <- obj$seurat_clusters == cluster
        if (any(cluster_mask)) {
          current_display_name <- as.character(obj$annotated_clusters[cluster_mask][1])
        } else {
          current_display_name <- paste0("Cluster_", cluster)
        }
      } else if (!is.null(current_names) && cluster %in% names(current_names)) {
        current_display_name <- current_names[[cluster]]
      } else {
        current_display_name <- paste0("Cluster_", cluster)
      }
      
      # Create the input with better labeling
      div(
        class = "form-group",
        style = "margin-bottom: 15px; padding: 10px; border: 1px solid #ddd; border-radius: 5px;",
        
        # Show both original number and current name
        tags$label(
          paste0("Original Cluster ", cluster), 
          style = "font-weight: bold; color: #333; font-size: 12px;"
        ),
        
        # Show current name if different from default
        if (current_display_name != paste0("Cluster_", cluster)) {
          tags$p(
            paste0("Current name: ", current_display_name), 
            style = "margin: 5px 0; color: #007bff; font-style: italic; font-size: 11px;"
          )
        } else {
          NULL
        },
        
        textInput(
          inputId = paste0("spatial_cluster_name_", cluster),
          label = NULL,
          value = current_display_name,
          placeholder = paste0("Enter name for cluster ", cluster)
        )
      )
    })
    
    do.call(tagList, input_list)
  })
  # Apply all cluster names - FIXED VERSION
  observeEvent(input$apply_all_spatial_names, {
    message("=== APPLY CLUSTER NAMES STARTED ===")
    
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      
      # Check if clusters exist
      if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
        showNotification("No clusters found in the object!", type = "error")
        return()
      }
      
      unique_clusters <- sort(unique(as.character(obj$seurat_clusters)))
      message(paste("Found clusters:", paste(unique_clusters, collapse = ", ")))
      
      # Collect all new names
      new_names <- list()
      for (cluster in unique_clusters) {
        input_id <- paste0("spatial_cluster_name_", cluster)
        new_name <- input[[input_id]]
        
        if (!is.null(new_name) && nchar(trimws(new_name)) > 0) {
          new_names[[cluster]] <- trimws(new_name)
          message(paste("Cluster", cluster, "->", new_names[[cluster]]))
        } else {
          new_names[[cluster]] <- paste0("Cluster_", cluster)
          message(paste("Cluster", cluster, "-> (default)", new_names[[cluster]]))
        }
      }
      
      # Update reactive value
      spatial_cluster_names(new_names)
      message("Reactive values updated")
      
      # Create new annotated clusters column
      obj$annotated_clusters <- as.factor(obj$seurat_clusters)
      
      # Get the current levels
      current_levels <- levels(obj$annotated_clusters)
      message(paste("Current levels:", paste(current_levels, collapse = ", ")))
      
      # Create new level names in the same order
      new_level_names <- sapply(current_levels, function(x) {
        if (x %in% names(new_names)) {
          return(new_names[[x]])
        } else {
          return(paste0("Cluster_", x))
        }
      })
      
      # Apply new names
      levels(obj$annotated_clusters) <- new_level_names
      message("New cluster names applied to factor levels")
      
      # Update the spatial object
      spatial_obj(obj)
      message("Spatial object updated")
      
      # FIXED: Remove the problematic showNotification call or fix it
      
      message("=== APPLY CLUSTER NAMES COMPLETED ===")
      
    }, error = function(e) {
      message(paste("=== APPLY CLUSTER NAMES ERROR ===", e$message))
      showNotification(paste("Error applying cluster names:", e$message), type = "error", duration = 10)
    })
  })
  
  # Reset cluster names
  observeEvent(input$reset_spatial_names, {
    req(spatial_obj())
    
    obj <- spatial_obj()
    unique_clusters <- sort(unique(as.character(obj$seurat_clusters)))
    
    # Reset to default names
    default_names <- setNames(paste0("Cluster_", unique_clusters), unique_clusters)
    spatial_cluster_names(default_names)
    
    # Remove annotated_clusters column if it exists
    if ("annotated_clusters" %in% colnames(obj@meta.data)) {
      obj$annotated_clusters <- NULL
      spatial_obj(obj)
    }
    
  })
  # Variable pour stocker le chemin du fichier temporaire
  spatial_tissue_temp_file <- reactiveVal(NULL)
  # Live update plots when cluster names change
  observe({
    req(spatial_obj())
    req(spatial_cluster_names())
    
    obj <- spatial_obj()
    
    # Determine which cluster column to use
    cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) {
      "annotated_clusters"
    } else {
      "seurat_clusters"
    }
    
    
    # --- tiny helper kept as you had ---
    .dark_theme_bg <- function(base_size = 14) {
      theme_void(base_size = base_size) +
        theme(
          plot.background   = element_rect(fill = "black", color = NA),
          panel.background  = element_rect(fill = "black", color = NA),
          legend.background = element_rect(fill = "black", color = NA),
          legend.key        = element_rect(fill = "black", color = NA),
          legend.text       = element_text(color = "white"),
          legend.title      = element_text(color = "white"),
          plot.title        = element_text(color = "white", hjust = 0.5)
        )
    }
    
    # ---------------- UMAP (annotation) ----------------
    output$spatial_annotation_umap <- renderPlot({
      req(spatial_obj())
      obj <- spatial_obj()
      
      cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) "annotated_clusters" else "seurat_clusters"
      
      has_umap <- "umap" %in% Seurat::Reductions(obj)
      shiny::validate(shiny::need(has_umap, "UMAP not found in this object."))
      
      p <- Seurat::DimPlot(
        obj,
        reduction  = "umap",
        group.by   = cluster_col,
        label      = TRUE,
        label.size = 5,
        repel      = TRUE,
        pt.size    = 1,
        raster     = FALSE
      ) +
        ggtitle("UMAP - Annotated Clusters")
      
      # Appliquer les modifications conditionnelles
      if(isTRUE(input$spatial_remove_legend)) {
        p <- p + NoLegend()
      }
      
      if(isTRUE(input$spatial_remove_axes)) {
        p <- p + NoAxes()
      }
      
      if (isTRUE(input$anno_dark_mode)) p <- p + .dark_theme_bg()
      
      # Stocker le plot pour le téléchargement
      spatial_umap_plot(p)
      p
    })
    # ---------------- Spatial (annotation) ----------------
    
    
    
    # REMPLACER output$spatial_annotation_tissue par renderImage
    output$spatial_annotation_tissue <- renderPlot({
      req(spatial_obj())
      obj <- spatial_obj()
      
      cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) "annotated_clusters" else "seurat_clusters"
      
      if (length(obj@images) == 0) {
        p <- ggplot() +
          geom_text(aes(0.5, 0.5, label = "No spatial image"), size = 6) +
          theme_void()
        if (isTRUE(input$anno_dark_mode)) {
          p <- p +
            theme(
              plot.background  = element_rect(fill = "black", color = NA),
              panel.background = element_rect(fill = "black", color = NA),
              text             = element_text(color = "white")
            )
        }
        spatial_tissue_plot(p)
        return(p)
      }
      
      if (!is.null(input$spatial_cluster_highlight) && input$spatial_cluster_highlight != "None") {
        obj_highlight <- obj
        obj_highlight@meta.data$highlight_group <- ifelse(
          obj_highlight@meta.data[[cluster_col]] == input$spatial_cluster_highlight,
          as.character(obj_highlight@meta.data[[cluster_col]]),
          "Other"
        )
        
        selected_cluster <- input$spatial_cluster_highlight
        colors <- c("Other" = "lightgray")
        colors[selected_cluster] <- "red"
        
        p <- Seurat::SpatialDimPlot(
          obj_highlight,
          group.by       = "highlight_group",
          label          = input$spatial_anno_show_labels,
          label.size     = 4,
          pt.size.factor = 1.6,
          alpha          = c(0.3, 1),
          cols           = colors
        ) +
          ggtitle(paste("Spatial View - Cluster", selected_cluster, "Highlighted"))
        
      } else {
        cluster_colors <- NULL
        if (!is.null(obj@misc$cluster_colors)) {
          cluster_colors <- obj@misc$cluster_colors
        }
        
        # Plot spatial SIMPLE sans les connexions pour Ã©viter les erreurs
        p <- Seurat::SpatialDimPlot(
          obj,
          group.by       = cluster_col,
          label          = input$spatial_anno_show_labels,
          label.size     = 4,
          pt.size.factor = 1.6,
          alpha          = c(0.3, 1),
          cols           = cluster_colors
        ) +
          ggtitle("Spatial View - Annotated Clusters")
      }
      
      if(isTRUE(input$spatial_remove_legend)) {
        p <- p + NoLegend()
      }
      
      if(isTRUE(input$spatial_remove_axes)) {
        p <- p + NoAxes()
      }
      
      if (isTRUE(input$anno_dark_mode)) p <- p + .dark_theme_bg()
      
      spatial_tissue_plot(p)
      p
    })
  })
  
  # Ajouter ces reactiveValues
  spatial_umap_plot <- reactiveVal(NULL)
  spatial_tissue_plot <- reactiveVal(NULL)
  
  # Téléchargement UMAP
  output$download_spatial_umap <- createDownloadHandler(
    reactive_data = spatial_umap_plot,
    object_name_reactive = reactive({ "spatial_umap_annotation" }),
    data_name = "umap_plot",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$spatial_plot_format }),
      width = 12,
      height = 10,
      dpi = 300
    ),
    show_modal = FALSE
  )
  
  output$download_spatial_tissue <- downloadHandler(
    filename = function() {
      paste0("spatial_tissue.", input$spatial_plot_format)
    },
    
    content = function(file) {
      req(spatial_tissue_temp_file())
      
      temp_file <- spatial_tissue_temp_file()
      
      if(input$spatial_plot_format == "png") {
        # Copier directement le fichier PNG affiché
        file.copy(temp_file, file, overwrite = TRUE)
      } else {
        # Pour autres formats, recharger et convertir
        req(spatial_obj())
        obj <- spatial_obj()
        cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) "annotated_clusters" else "seurat_clusters"
        cluster_colors <- if (!is.null(obj@misc$cluster_colors)) obj@misc$cluster_colors else NULL
        
        p <- Seurat::SpatialDimPlot(obj, group.by = cluster_col, label = input$spatial_anno_show_labels, label.size = 6, pt.size.factor = 2.5, alpha = c(0.1, 0.9), cols = cluster_colors) + ggtitle("Spatial View - Annotated Clusters")
        
        if(isTRUE(input$spatial_remove_legend)) p <- p + NoLegend()
        if(isTRUE(input$spatial_remove_axes)) p <- p + NoAxes()
        
        ggsave(file, plot = p, width = 16, height = 14, dpi = 300, bg = "white", device = input$spatial_plot_format)
      }
    }
  )
  
  # Observer pour mettre à jour les choix de clusters avec les noms annotés
  observe({
    req(spatial_obj())
    obj <- spatial_obj()
    
    cluster_col <- get_cluster_column(obj)
    
    if (cluster_col %in% colnames(obj@meta.data)) {
      # Use annotated names if available, otherwise original cluster numbers
      cluster_choices <- get_current_cluster_names(obj)
      
      # Update all cluster selection inputs with annotated names
      updateSelectInput(session, "spatial_cluster_highlight", 
                        choices = c("None", cluster_choices),
                        selected = "None")
      
      updateSelectInput(session, "spatial_plotly_cluster_highlight", 
                        choices = c("None", cluster_choices),
                        selected = "None")
    }
  })
  
  # Function to get current cluster names (either annotated or original)
  get_current_cluster_names <- function(obj) {
    if ("annotated_clusters" %in% colnames(obj@meta.data)) {
      return(levels(obj$annotated_clusters))
    } else {
      return(levels(obj$seurat_clusters))
    }
  }
  
  # Function to get current cluster column name
  get_cluster_column <- function(obj) {
    if ("annotated_clusters" %in% colnames(obj@meta.data)) {
      return("annotated_clusters")
    } else {
      return("seurat_clusters")
    }
  }
  # Observer pour gérer l'affichage des plots
  observeEvent(input$annotation_selected_plot, {
    # Cacher tous les panels
    shinyjs::hide("no_plot_panel")
    shinyjs::hide("violin_plot_panel")
    shinyjs::hide("dot_plot_panel")
    shinyjs::hide("feature_plot_panel")
    
    # Montrer seulement le panel sélectionné
    if (input$annotation_selected_plot == "none") {
      shinyjs::show("no_plot_panel")
    } else if (input$annotation_selected_plot == "violin") {
      shinyjs::show("violin_plot_panel")
    } else if (input$annotation_selected_plot == "dot") {
      shinyjs::show("dot_plot_panel")
    } else if (input$annotation_selected_plot == "feature") {
      shinyjs::show("feature_plot_panel")
    }
  })
  
  # Display selected plot in annotation tab
  # 3. Debug dans l'affichage des plots d'annotation
  output$annotation_selected_plot_display <- renderUI({
    selected_plot <- input$annotation_selected_plot
    
    # Force la réactivité sur le changement de sélection
    req(selected_plot)
    
    if (selected_plot == "none") {
      return(
        div(style = "text-align: center; padding: 50px; color: #666;",
            h4("No plot selected"),
            p("Choose a plot type above to display gene expression data"))
      )
    }
    
    # Force un invalidate quand le bouton refresh est cliqué
    input$refresh_annotation_plot
    
    # Force la mise à jour quand les plots changent (sans isolate!)
    violin_exists <- !is.null(spatial_violin_plot_reactive())
    dot_exists <- !is.null(spatial_dot_plot_reactive())
    feature_exists <- !is.null(spatial_feature_plot_reactive())
    
    plot_height <- "350px"
    
    if (selected_plot == "violin") {
      if (violin_exists) {
        plotOutput("annotation_violin_display", height = plot_height)
      } else {
        div(style = "text-align: center; padding: 30px; color: #FF9800;",
            icon("exclamation-triangle", style = "font-size: 2em;"),
            h5("No violin plot available"),
            p("Generate a violin plot in the 'Gene Expression Visualization' tab first"))
      }
    } else if (selected_plot == "dot") {
      if (dot_exists) {
        plotOutput("annotation_dot_display", height = plot_height)
      } else {
        div(style = "text-align: center; padding: 30px; color: #FF9800;",
            icon("exclamation-triangle", style = "font-size: 2em;"),
            h5("No dot plot available"),
            p("Generate a dot plot in the 'Gene Expression Visualization' tab first"))
      }
    } else if (selected_plot == "feature") {
      if (feature_exists) {
        plotOutput("annotation_feature_display", height = plot_height)
      } else {
        div(style = "text-align: center; padding: 30px; color: #FF9800;",
            icon("exclamation-triangle", style = "font-size: 2em;"),
            h5("No feature plot available"),
            p("Generate a feature plot in the 'Gene Expression Visualization' tab first"))
      }
    } else {
      div(style = "text-align: center; padding: 30px; color: #F44336;",
          icon("times-circle", style = "font-size: 2em;"),
          h5("Unknown plot type selected"))
    }
  })
  
  
  output$annotation_feature_display <- renderPlot({
    plot_obj <- spatial_feature_plot_reactive()
    if (!is.null(plot_obj)) {
      plot_obj
    } else {
      # Plot vide avec message
      ggplot() + 
        theme_void() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No feature plot available\nGenerate one in Gene Expression Visualization tab",
                 size = 5, color = "gray50")
    }
  })
  
  
  
  output$annotation_violin_display <- renderPlot({
    plot_obj <- spatial_violin_plot_reactive()
    if (!is.null(plot_obj)) {
      plot_obj
    } else {
      ggplot() + 
        theme_void() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No violin plot available\nGenerate one in Gene Expression Visualization tab",
                 size = 5, color = "gray50")
    }
  })
  
  output$annotation_dot_display <- renderPlot({
    plot_obj <- spatial_dot_plot_reactive()
    if (!is.null(plot_obj)) {
      plot_obj
    } else {
      ggplot() + 
        theme_void() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = "No dot plot available\nGenerate one in Gene Expression Visualization tab",
                 size = 5, color = "gray50")
    }
  })
  
  # 5. Observer pour surveiller les changements des reactive values
  observe({
    feature_val <- spatial_feature_plot_reactive()
    if (!is.null(feature_val)) {
      print(">>> Feature plot reactive updated! Not NULL anymore")
    }
  })
  
  observe({
    violin_val <- spatial_violin_plot_reactive()
    if (!is.null(violin_val)) {
      print(">>> Violin plot reactive updated! Not NULL anymore")
    }
  })
  
  # Save annotated spatial object - VERSION CORRIGÉE avec createDownloadHandler
  output$save_annotated_spatial_object <- createDownloadHandler(
    reactive_data = spatial_obj,  # Ton reactive value spatial_obj
    object_name_reactive = reactive({ 
      # Générer un nom approprié
      if (!is.null(spatial_obj())) {
        paste0("spatial_annotated_", Sys.Date())
      } else {
        "spatial_object"
      }
    }),
    data_name = "spatial_annotated",
    download_type = "seurat",
    show_modal = TRUE  
  )
  
  
  
  
  ############################## Spatial Tissue Viewer Server ##############################
  
  # Shared window for sync (full-res pixel space), define early
  sync_view <- reactiveVal(NULL)
  
  
  
  ## If you upload large files, consider raising the upload limit at app start:
  ## options(shiny.maxRequestSize = 1024^2 * 1024)  # ~1 GB
  
  # --- Reactive state for the H&E viewer ---
  hne_state <- reactiveValues(
    # Source and backend
    src_path    = NULL,     # original file path for TIFF (used by VIPS)
    is_tiff     = FALSE,    # TRUE when using VIPS path for .tif/.tiff
    
    # For PNG/JPG backend
    img_magick  = NULL,     # FULL RES magick image object (PNG/JPG)
    
    # Current display image (cropped/scaled for the viewport)
    img_display = NULL,
    
    # Level-0 (full-res) dimensions of the whole image
    w = NULL,               # width  in pixels
    h = NULL,               # height in pixels
    
    # Current view window in pixel space: [xmin, xmax, ymin, ymax]
    view = NULL,
    
    # Target display size for the longest edge (keeps rendering snappy)
    display_size = 1200
  )
  
  # ============================== Helpers ==============================
  
  # Keep the viewing window inside image bounds and enforce a minimum size
  clamp_view <- function(win, w, h, min_frac = 0.02) {
    if (is.null(win)) return(c(0, w, 0, h))
    xmin <- max(0, min(win[1], win[2])); xmax <- min(w, max(win[1], win[2]))
    ymin <- max(0, min(win[3], win[4])); ymax <- min(h, max(win[3], win[4]))
    
    min_w <- max(16, w * min_frac); min_h <- max(16, h * min_frac)
    if ((xmax - xmin) < min_w) { cx <- (xmin + xmax)/2; xmin <- cx - min_w/2; xmax <- cx + min_w/2 }
    if ((ymax - ymin) < min_h) { cy <- (ymin + ymax)/2; ymin <- cy - min_h/2; ymax <- cy + min_h/2 }
    
    xmin <- max(0, xmin); ymin <- max(0, ymin); xmax <- min(w, xmax); ymax <- min(h, ymax)
    c(xmin, xmax, ymin, ymax)
  }
  
  # Compute a zoomed window centered on (x,y) by a given factor
  zoom_center <- function(win, w, h, x, y, factor) {
    if (is.null(win)) win <- c(0, w, 0, h)
    cur_w <- (win[2] - win[1]) / factor
    cur_h <- (win[4] - win[3]) / factor
    c(x - cur_w/2, x + cur_w/2, y - cur_h/2, y + cur_h/2)
  }
  
  # Prepare cropped & scaled display image from a FULL magick image (PNG/JPG backend)
  prepare_display_image <- function(img_full, view_window, display_size) {
    crop_geo <- sprintf("%dx%d+%d+%d",
                        as.integer(view_window[2] - view_window[1]),
                        as.integer(view_window[4] - view_window[3]),
                        as.integer(view_window[1]),
                        as.integer(view_window[3]))
    img_crop <- magick::image_crop(img_full, geometry = crop_geo)
    info <- magick::image_info(img_crop)
    
    if (max(info$width, info$height) > display_size) {
      scale <- display_size / max(info$width, info$height)
      img_crop <- magick::image_scale(
        img_crop,
        geometry = sprintf("%dx%d", as.integer(info$width * scale), as.integer(info$height * scale))
      )
    }
    img_crop
  }
  
  # Check that vips is in PATH
  ensure_vips <- function() {
    out <- tryCatch(system2("vips", "--version", stdout = TRUE, stderr = TRUE), error = function(e) character())
    nzchar(paste(out, collapse = ""))
  }
  
  # Probe TIFF size quickly via vips header
  vips_tiff_size <- function(path) {
    # Try `vipsheader -a`
    out <- tryCatch(system2("vipsheader", c("-a", path), stdout = TRUE), error = function(e) character())
    if (length(out)) {
      w <- suppressWarnings(as.numeric(sub(".*width:\\s*",  "", out[grep("width:",  out, ignore.case = TRUE)[1]])))
      h <- suppressWarnings(as.numeric(sub(".*height:\\s*", "", out[grep("height:", out, ignore.case = TRUE)[1]])))
      if (is.finite(w) && is.finite(h) && w > 0 && h > 0) return(c(w, h))
    }
    
    # Try plain `vipsheader`
    out2 <- tryCatch(system2("vipsheader", path, stdout = TRUE), error = function(e) character())
    if (length(out2)) {
      # Typical line: "in: 100000x80000 uchar, 3 bands, ..."
      m <- regmatches(out2, regexpr("\\b(\\d+)x(\\d+)\\b", out2))
      if (length(m) && nzchar(m[1])) {
        nums <- as.numeric(strsplit(m[1], "x")[[1]])
        if (length(nums) == 2 && all(is.finite(nums)) && all(nums > 0)) return(nums)
      }
    }
    
    # Fallback (can be slower on huge files): use magick to read meta
    img <- tryCatch(magick::image_read(path), error = function(e) NULL)
    if (!is.null(img)) {
      info <- magick::image_info(img)
      w <- as.numeric(info$width[1]); h <- as.numeric(info$height[1])
      if (is.finite(w) && is.finite(h) && w > 0 && h > 0) return(c(w, h))
    }
    
    c(NA_real_, NA_real_)
  }
  
  # Read the current window from pyramidal TIFF using VIPS, with robust clamping
  read_tiff_region_via_vips <- function(path, view, display_size, full_w, full_h) {
    # ---- clamp to image bounds (integers) ----
    xmin <- max(0L, as.integer(floor(min(view[1], view[2]))))
    xmax <- min(as.integer(ceiling(max(view[1], view[2]))), as.integer(full_w))
    ymin <- max(0L, as.integer(floor(min(view[3], view[4]))))
    ymax <- min(as.integer(ceiling(max(view[3], view[4]))), as.integer(full_h))
    
    # ensure positive size
    w <- max(1L, xmax - xmin)
    h <- max(1L, ymax - ymin)
    
    # ---- enforce a minimum crop size (avoid 1x1) ----
    min_side <- 64L  # you can tweak (32/64/128)
    if (w < min_side || h < min_side) {
      cx <- (xmin + xmax) %/% 2
      cy <- (ymin + ymax) %/% 2
      w <- max(min_side, w); h <- max(min_side, h)
      xmin <- max(0L, cx - w %/% 2); xmax <- min(as.integer(full_w), xmin + w)
      ymin <- max(0L, cy - h %/% 2); ymax <- min(as.integer(full_h), ymin + h)
      # recompute in case we hit borders
      w <- max(1L, xmax - xmin); h <- max(1L, ymax - ymin)
    }
    
    # ---- scale so the larger side maps to display_size ----
    scale <- display_size / max(w, h)
    
    tmp_png  <- tempfile(fileext = ".png")
    tmp_crop <- tempfile(fileext = ".tif")  # intermediate cropped tiff
    
    # 1) Crop (use 'vips crop' which exists on older VIPS too)
    cmd_crop <- c("crop", path, tmp_crop, xmin, ymin, w, h)
    status1 <- tryCatch(system2("vips", cmd_crop, stdout = TRUE, stderr = TRUE),
                        error = function(e) paste("err:", e$message))
    if (!file.exists(tmp_crop)) stop("VIPS crop failed: ", paste(status1, collapse = " "))
    
    # 2) Resize to display size
    cmd_resz <- c("resize", tmp_crop, tmp_png, scale)
    status2 <- tryCatch(system2("vips", cmd_resz, stdout = TRUE, stderr = TRUE),
                        error = function(e) paste("err:", e$message))
    if (!file.exists(tmp_png)) {
      unlink(tmp_crop, force = TRUE)
      stop("VIPS resize failed: ", paste(status2, collapse = " "))
    }
    
    on.exit({ unlink(tmp_crop, force = TRUE); unlink(tmp_png, force = TRUE) }, add = TRUE)
    magick::image_read(tmp_png)
  }
  
  
  # ============================== File upload ==============================
  
  observeEvent(input$hne_file, {
    req(input$hne_file)
    nm  <- input$hne_file$name
    ext <- tolower(tools::file_ext(nm))
    fp  <- input$hne_file$datapath
    message("[HNE] Incoming file: ", nm)
    
    # Reset state
    hne_state$src_path    <- NULL
    hne_state$is_tiff     <- FALSE
    hne_state$img_magick  <- NULL
    hne_state$img_display <- NULL
    hne_state$w <- NULL; hne_state$h <- NULL
    hne_state$view <- NULL
    
    # ---- PNG/JPG -> magick backend ----
    if (ext %in% c("png", "jpg", "jpeg")) {
      tryCatch({
        img  <- magick::image_read(fp)
        info <- magick::image_info(img)
        hne_state$img_magick <- img
        hne_state$w <- as.numeric(info$width[1]); hne_state$h <- as.numeric(info$height[1])
        hne_state$view <- c(0, hne_state$w, 0, hne_state$h)
        sync_view(hne_state$view)
        hne_state$img_display <- prepare_display_image(img, hne_state$view, hne_state$display_size)
        showNotification(sprintf("Loaded %s %dx%d px", toupper(ext), hne_state$w, hne_state$h), type = "message")
      }, error = function(e) {
        showNotification(paste("Image load error:", e$message), type = "error")
      })
      return(invisible())
    }
    
    # ---- BigTIFF (.tif/.tiff) -> VIPS region backend ----
    if (ext %in% c("tif", "tiff")) {
      if (!ensure_vips()) {
        showNotification("VIPS not found in PATH. Install libvips to handle BigTIFF pyramidal.", type = "error", duration = 10)
        return(invisible())
      }
      dims <- vips_tiff_size(fp)
      if (any(is.na(dims)) || any(dims == 0)) {
        showNotification("Could not read TIFF dimensions via VIPS.", type = "error")
        return(invisible())
      }
      hne_state$src_path <- fp
      hne_state$is_tiff  <- TRUE
      hne_state$w <- dims[1]; hne_state$h <- dims[2]
      hne_state$view <- c(0, hne_state$w, 0, hne_state$h)
      sync_view(hne_state$view)
      
      # Prepare first display image from base view
      hne_state$img_display <- tryCatch(
        read_tiff_region_via_vips(
          path         = hne_state$src_path,
          view         = hne_state$view,
          display_size = hne_state$display_size,
          full_w       = hne_state$w,
          full_h       = hne_state$h
        ),
        error = function(e) { showNotification(paste("VIPS region read failed:", e$message), type = "error"); NULL }
      )
      
      if (!is.null(hne_state$img_display)) {
        showNotification(sprintf("BigTIFF detected %dx%d px (pyramidal). Using VIPS for region rendering.", hne_state$w, hne_state$h),
                         type = "message", duration = 8)
      }
      return(invisible())
    }
    
    showNotification("Unsupported format. Use PNG/JPG/TIFF.", type = "error")
  })
  
  # ============================== Interactions ==============================
  
  # petit helper : plancher adaptatif (en px) pour éviter les fenêtres ~1 px
  min_window_px <- function(w, h) {
    # au moins 32 px, et ~0.1% de la dimension (utile pour lames géantes)
    max(32L, as.integer(0.001 * max(w, h)))
  }
  
  ############################## Spatial Tissue Viewer Server (interactions + renderer) ##############################
  
  # NOTE: keep your existing hne_state, clamp_view, zoom_center, prepare_display_image,
  # ensure_vips, vips_tiff_size, read_tiff_region_via_vips, min_window_px helpers as you had.
  
  # ---- File upload handling (unchanged logic) ----
  # (Tu peux garder exactement le bloc de "File upload" qui marchait chez toi.)
  
  # ------------------ Interactions ------------------
  
  observeEvent(input$hne_view_brush, {
    req(hne_state$img_display, hne_state$view, hne_state$w, hne_state$h)
    b <- input$hne_view_brush
    if (is.null(b)) return()
    
    # taille du raster affiché
    dinfo  <- magick::image_info(hne_state$img_display)
    disp_w <- as.numeric(dinfo$width)
    disp_h <- as.numeric(dinfo$height)
    
    # clamp + ordre en pixels d'affichage
    bxmin <- max(0, min(disp_w, b$xmin)); bxmax <- max(0, min(disp_w, b$xmax))
    bymin <- max(0, min(disp_h, b$ymin)); bymax <- max(0, min(disp_h, b$ymax))
    if (bxmin > bxmax) { tmp <- bxmin; bxmin <- bxmax; bxmax <- tmp }
    if (bymin > bymax) { tmp <- bymin; bymin <- bymax; bymax <- tmp }
    
    # ignorer les micro-sélections
    if ((bxmax - bxmin) < 20 || (bymax - bymin) < 20) return()
    
    # fenêtre full-res courante
    cur <- hne_state$view
    vw  <- cur[2] - cur[1]
    vh  <- cur[4] - cur[3]
    
    # mapping linéaire affichage(px) -> full-res(px) dans la vue courante
    xmin_new <- cur[1] + (bxmin / disp_w) * vw
    xmax_new <- cur[1] + (bxmax / disp_w) * vw
    ymin_new <- cur[3] + (bymin / disp_h) * vh
    ymax_new <- cur[3] + (bymax / disp_h) * vh
    
    new_view <- clamp_view(c(xmin_new, xmax_new, ymin_new, ymax_new), hne_state$w, hne_state$h, min_frac = 0)
    
    # plancher absolu (évite le "mega-zoom blanc")
    minpx <- max(32L, as.integer(0.001 * max(hne_state$w, hne_state$h)))
    ww <- new_view[2] - new_view[1]; hh <- new_view[4] - new_view[3]
    if (ww < minpx || hh < minpx) {
      cx <- (new_view[1] + new_view[2]) / 2
      cy <- (new_view[3] + new_view[4]) / 2
      new_view <- clamp_view(c(cx - minpx/2, cx + minpx/2, cy - minpx/2, cy + minpx/2),
                             hne_state$w, hne_state$h, min_frac = 0)
    }
    
    hne_state$view <- new_view
    sync_view(new_view)
    
    session$resetBrush("hne_view_brush")
  })
  
  
  
  # Buttons zoom in/out/reset
  observeEvent(input$hne_zoom_in, {
    req(hne_state$view, hne_state$w, hne_state$h)
    f <- input$hne_zoom_factor %||% 1.5
    win <- hne_state$view; cx <- (win[1] + win[2]) / 2; cy <- (win[3] + win[4]) / 2
    new_view <- zoom_center(win, hne_state$w, hne_state$h, cx, cy, f)
    hne_state$view <- clamp_view(new_view, hne_state$w, hne_state$h, min_frac = 0.001)
    sync_view(new_view)
    
  })
  observeEvent(input$hne_zoom_out, {
    req(hne_state$view, hne_state$w, hne_state$h)
    f <- input$hne_zoom_factor %||% 1.5
    win <- hne_state$view; cx <- (win[1] + win[2]) / 2; cy <- (win[3] + win[4]) / 2
    new_view <- zoom_center(win, hne_state$w, hne_state$h, cx, cy, 1/f)
    hne_state$view <- clamp_view(new_view, hne_state$w, hne_state$h, min_frac = 0)
    sync_view(new_view)
    
  })
  
  
  observeEvent(input$hne_reset, {
    req(hne_state$w, hne_state$h)
    hne_state$view <- c(0, hne_state$w, 0, hne_state$h)
    sync_view(hne_state$view)   # <-- was sync_view(new_view) -> FIX
  })
  
  
  # Pan buttons (move by a fraction of the current window size)
  pan_by <- function(dx_frac = 0, dy_frac = 0) {
    req(hne_state$view, hne_state$w, hne_state$h)
    win <- hne_state$view
    ww  <- win[2] - win[1]
    hh  <- win[4] - win[3]
    dx  <- dx_frac * ww
    dy  <- dy_frac * hh
    new_view <- c(win[1] + dx, win[2] + dx, win[3] + dy, win[4] + dy)
    hne_state$view <- clamp_view(new_view, hne_state$w, hne_state$h, min_frac = 0)
    sync_view(hne_state$view)    
    
  }
  
  observeEvent(input$hne_pan_left,  { pan_by(dx_frac = -0.40, dy_frac =  0.00) })
  observeEvent(input$hne_pan_right, { pan_by(dx_frac =  0.40, dy_frac =  0.00) })
  observeEvent(input$hne_pan_up,    { pan_by(dx_frac =  0.00, dy_frac = -0.40) }) # up = y-
  observeEvent(input$hne_pan_down,  { pan_by(dx_frac =  0.00, dy_frac =  0.40) }) # down = y+
  
  
  # ------------------ Recompute display image on view change ------------------
  
  observe({
    req(hne_state$view, hne_state$w, hne_state$h)
    if (isTRUE(hne_state$is_tiff)) {
      req(hne_state$src_path)
      hne_state$img_display <- tryCatch(
        read_tiff_region_via_vips(
          hne_state$src_path,
          hne_state$view,
          hne_state$display_size,
          full_w = hne_state$w,
          full_h = hne_state$h
        ),
        error = function(e) { showNotification(paste("VIPS region read failed:", e$message), type = "error"); NULL }
      )
    } else if (!is.null(hne_state$img_magick)) {
      hne_state$img_display <- prepare_display_image(hne_state$img_magick, hne_state$view, hne_state$display_size)
    }
  })
  
  # ------------------ Renderer (normalized canvas 0..1) ------------------
  
  output$hne_view <- renderPlot({
    if (is.null(hne_state$img_display)) {
      plot.new(); text(0.5, 0.5, "No image loaded\nPlease select an H&E image file", cex = 1.4, col = "gray50")
      return(invisible())
    }
    info   <- magick::image_info(hne_state$img_display)
    disp_w <- as.numeric(info$width)
    disp_h <- as.numeric(info$height)
    
    par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
    plot.new()
    plot.window(xlim = c(0, disp_w), ylim = c(0, disp_h), asp = 1)  # 1:1 => pas d’écrasement aux bords
    
    rasterImage(as.raster(hne_state$img_display), 0, 0, disp_w, disp_h, interpolate = FALSE)
    
    # (optionnel) overlay
    if (!is.null(hne_state$view) && !is.null(hne_state$w)) {
      vw <- hne_state$view[2] - hne_state$view[1]
      text(disp_w*0.05, disp_h*0.95, paste0("Zoom: ", round(hne_state$w / max(1, vw), 1), "x"),
           col = "red", cex = 1.0, font = 2, pos = 4)
    }
  }, height = 650)
  
  
  
  
  ## =========================== Shared window for sync (full-res pixel space) ===========================
  # Whenever H&E view changes, we copy it into this reactiveVal so the right panel follows.
  sync_view <- reactiveVal(NULL)
  
  # each time you set/update hne_state$view, also do:
  # sync_view(hne_state$view)
  # (On intial load too, after you set hne_state$view <- c(0, w, 0, h))
  
  ## =========================== Extract spot coordinates in image pixel space ===========================
  # Returns a data.frame with columns: barcode, x_img, y_img (+ the selected feature)
  spot_coords <- reactive({
    req(spatial_obj())
    obj <- spatial_obj()
    
    imgs <- Seurat::Images(obj)
    if (is.null(imgs) || length(imgs) < 1) return(NULL)
    imgname <- imgs[1]
    
    coords <- tryCatch(Seurat::GetTissueCoordinates(object = obj, image = imgname),
                       error = function(e) NULL)
    if (is.null(coords)) {
      coords <- tryCatch(Seurat::GetTissueCoordinates(object = obj[[imgname]]),
                         error = function(e) NULL)
    }
    if (is.null(coords)) return(NULL)
    
    coords <- as.data.frame(coords)
    if (nrow(coords) == 0) return(NULL)
    
    cn <- tolower(colnames(coords))
    pick <- function(opts) {
      hit <- match(tolower(opts), cn)
      if (all(is.na(hit))) return(NULL)
      colnames(coords)[na.omit(hit)[1]]
    }
    x_col <- pick(c("imagecol","pxl_col_in_fullres","pxl_col_in_hires","x","col"))
    y_col <- pick(c("imagerow","pxl_row_in_fullres","pxl_row_in_hires","y","row"))
    if (is.null(x_col) || is.null(y_col)) return(NULL)
    
    out <- data.frame(
      barcode = rownames(coords),
      x_img = as.numeric(coords[[x_col]]),
      y_img = as.numeric(coords[[y_col]]),
      stringsAsFactors = FALSE
    )
    
    md <- obj@meta.data
    md$barcode <- rownames(md)
    keep <- intersect(c("barcode","seurat_clusters","nFeature_Spatial","nCount_Spatial"), colnames(md))
    if (length(keep)) {
      out <- merge(out, md[, keep, drop = FALSE], by = "barcode", all.x = TRUE, sort = FALSE)
    }
    out
  })
  
  
  output$trans_view <- renderPlot({
    req(spatial_obj())
    if (is.null(sync_view())) {
      if (!is.null(hne_state$w) && !is.null(hne_state$h)) {
        sync_view(c(0, hne_state$w, 0, hne_state$h))
      } else {
        plot.new(); text(0.5, 0.5, "Waiting for H&E view...", col="gray40"); return(invisible())
      }
    }
    sc <- spot_coords()
    if (is.null(sc) || nrow(sc) == 0) {
      plot.new(); text(0.5, 0.5, "No pixel-space coordinates (Seurat image coords missing)", col="gray40")
      return(invisible())
    }
    win <- sync_view()
    xmin <- win[1]; xmax <- win[2]; ymin <- win[3]; ymax <- win[4]
    inside <- sc$x_img >= xmin & sc$x_img <= xmax & sc$y_img >= ymin & sc$y_img <= ymax
    sc_in  <- sc[inside, , drop = FALSE]
    feat <- input$trans_feature %||% "seurat_clusters"
    suppressPackageStartupMessages(library(ggplot2))
    p <- ggplot(sc_in, aes(x = x_img, y = y_img)) +
      {
        if (feat == "seurat_clusters" && "seurat_clusters" %in% colnames(sc_in)) {
          geom_point(aes(color = factor(seurat_clusters)), size = input$trans_point_size %||% 0.6)
        } else if (feat == "nFeature_Spatial" && "nFeature_Spatial" %in% colnames(sc_in)) {
          geom_point(aes(color = nFeature_Spatial), size = input$trans_point_size %||% 0.6)
        } else if (feat == "nCount_Spatial" && "nCount_Spatial" %in% colnames(sc_in)) {
          geom_point(aes(color = nCount_Spatial), size = input$trans_point_size %||% 0.6)
        } else {
          geom_point(size = input$trans_point_size %||% 0.6, color = "steelblue")
        }
      } +
      coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
      scale_y_reverse() +     # single flip (don’t also invert the limits)
      theme_void(base_size = 11) +
      theme(legend.position = "right", plot.margin = margin(5,5,5,5))
    if (feat == "seurat_clusters" && "seurat_clusters" %in% colnames(sc_in)) {
      p <- p + guides(color = guide_legend(title = "Cluster"))
    } else if (feat %in% c("nFeature_Spatial","nCount_Spatial")) {
      p <- p + guides(color = guide_colorbar(title = feat))
    }
    p
  }, height = 650)
  
  
  ## =========================== Render transcriptomics plot ===========================
  
  ## Pixel-space spot coordinates (robust to Seurat v4/v5 + VisiumV2)
  spot_coords <- reactive({
    req(spatial_obj())
    obj <- spatial_obj()
    
    # 1) pick first spatial image
    imgs <- Seurat::Images(obj)
    if (is.null(imgs) || length(imgs) < 1) return(NULL)
    imgname <- imgs[1]
    
    # 2) get coordinates via official API (works with VisiumV2)
    coords <- tryCatch(
      Seurat::GetTissueCoordinates(object = obj, image = imgname),
      error = function(e) NULL
    )
    # some older builds also accept passing the image object itself:
    if (is.null(coords)) {
      coords <- tryCatch(
        Seurat::GetTissueCoordinates(object = obj[[imgname]]),
        error = function(e) NULL
      )
    }
    if (is.null(coords)) return(NULL)
    
    coords <- as.data.frame(coords)
    if (nrow(coords) == 0) return(NULL)
    
    # 3) choose the pixel-space columns robustly
    cn <- tolower(colnames(coords))
    get_col <- function(opts) {
      # return first present among candidate names
      hit <- match(tolower(opts), cn)
      if (all(is.na(hit))) return(NULL)
      colnames(coords)[na.omit(hit)[1]]
    }
    
    x_col <- get_col(c("imagecol", "pxl_col_in_fullres", "pxl_col_in_hires", "x", "col"))
    y_col <- get_col(c("imagerow", "pxl_row_in_fullres", "pxl_row_in_hires", "y", "row"))
    if (is.null(x_col) || is.null(y_col)) return(NULL)
    
    out <- data.frame(
      barcode = rownames(coords),
      x_img   = as.numeric(coords[[x_col]]),
      y_img   = as.numeric(coords[[y_col]]),
      stringsAsFactors = FALSE
    )
    
    # 4) attach basic meta (if present)
    md <- obj@meta.data
    md$barcode <- rownames(md)
    keep <- intersect(c("barcode", "seurat_clusters", "nFeature_Spatial", "nCount_Spatial"), colnames(md))
    if (length(keep)) {
      out <- merge(out, md[, keep, drop = FALSE], by = "barcode", all.x = TRUE, sort = FALSE)
    }
    
    out
  })
  
  
  output$trans_view <- renderPlot({
    req(spatial_obj(), sync_view())
    sc <- spot_coords()
    if (is.null(sc) || nrow(sc) == 0) {
      plot.new(); text(0.5, 0.5, "No pixel-space spot coordinates.\n(Seurat image coords missing)", cex = 1.1, col = "gray40")
      return(invisible())
    }
    
    win <- sync_view()
    xmin <- win[1]; xmax <- win[2]; ymin <- win[3]; ymax <- win[4]
    inside <- sc$x_img >= xmin & sc$x_img <= xmax & sc$y_img >= ymin & sc$y_img <= ymax
    sc_in  <- sc[inside, , drop = FALSE]
    
    feat <- input$trans_feature %||% "seurat_clusters"
    suppressPackageStartupMessages(library(ggplot2))
    p <- ggplot(sc_in, aes(x = x_img, y = y_img)) +
      {
        if (feat == "seurat_clusters" && "seurat_clusters" %in% colnames(sc_in)) {
          geom_point(aes(color = factor(seurat_clusters)), size = input$trans_point_size %||% 0.6)
        } else if (feat == "nFeature_Spatial" && "nFeature_Spatial" %in% colnames(sc_in)) {
          geom_point(aes(color = nFeature_Spatial), size = input$trans_point_size %||% 0.6)
        } else if (feat == "nCount_Spatial" && "nCount_Spatial" %in% colnames(sc_in)) {
          geom_point(aes(color = nCount_Spatial), size = input$trans_point_size %||% 0.6)
        } else {
          geom_point(size = input$trans_point_size %||% 0.6, color = "steelblue")
        }
      } +
      # Match H&E window orientation
      coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymax, ymin), expand = FALSE) +
      scale_y_reverse() +
      theme_void(base_size = 11) +
      theme(legend.position = "right", plot.margin = margin(5,5,5,5))
    
    if (feat == "seurat_clusters" && "seurat_clusters" %in% colnames(sc_in)) {
      p <- p + guides(color = guide_legend(title = "Cluster"))
    } else if (feat %in% c("nFeature_Spatial","nCount_Spatial")) {
      p <- p + guides(color = guide_colorbar(title = feat))
    }
    
    p
  }, height = 650)
  
  
  
  ############################### SPATIAL CO-LOCALIZATION — SERVER ##############################
  
  # Pré-requis
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(DT)
    library(Seurat)
  })
  
  coloc_results <- reactiveVal(NULL)
  
  # Conversion µm -> bins (résolution exprimée en µm par bin)
  convert_um_to_bins <- function(distance_um, bin_size_um) distance_um / bin_size_um
  
  # Create visual bar for distance metrics
  create_distance_bar <- function(distance_um, max_distance_um) {
    if (is.na(distance_um) || is.na(max_distance_um) || max_distance_um == 0) {
      return("<span style='color: grey;'>N/A</span>")
    }
    pct <- round(100 * distance_um / max_distance_um, 1)
    pct <- min(max(pct, 0), 100)
    if (pct < 25) {
      color <- "#28a745"  
    } else if (pct < 50) {
      color <- "#85c88a"  
    } else if (pct < 75) {
      color <- "#ffc107"  
    } else {
      color <- "#fd7e14"  
    }
    
    bar_html <- paste0(
      '<div style="width: 100%; height: 20px; border: 1px solid #ddd; ',
      'border-radius: 4px; overflow: hidden; background-color: #f8f9fa;">',
      '<div style="width: ', pct, '%; height: 100%; ',
      'background: linear-gradient(90deg, ', color, ' 0%, ', 
      adjustcolor(color, alpha.f = 0.7), ' 100%); ',
      'display: flex; align-items: center; justify-content: center; ',
      'transition: width 0.3s ease;" ',
      'title="', distance_um, ' µm (', pct, '% of max distance)">',
      '<span style="font-size: 10px; color: white; font-weight: bold; ',
      'text-shadow: 1px 1px 2px rgba(0,0,0,0.5);">',
      distance_um, ' µm',
      '</span>',
      '</div></div>'
    )
    
    return(bar_html)
  }
  
  
  # Récup clusters (colonne utilisée)
  get_cluster_info <- function(spatial_obj) {
    cluster_col <- if ("annotated_clusters" %in% colnames(spatial_obj@meta.data)) {
      "annotated_clusters"
    } else if ("seurat_clusters" %in% colnames(spatial_obj@meta.data)) {
      "seurat_clusters"
    } else {
      NULL
    }
    if (is.null(cluster_col)) {
      spots_clusters <- setNames(rep("No_clusters", ncol(spatial_obj)), colnames(spatial_obj))
      cluster_names  <- "No_clusters"
    } else {
      spots_clusters <- setNames(as.character(spatial_obj@meta.data[[cluster_col]]), colnames(spatial_obj))
      cluster_names  <- unique(spots_clusters)
    }
    list(spots = spots_clusters, names = cluster_names, column_used = cluster_col)
  }
  
  # Stat inter-cluster
  # Stat inter-cluster - VERSION AMÉLIORÉE avec pourcentages bidirectionnels
  calculate_inter_cluster_stats <- function(cluster_pairs_df, g1_spots, g2_spots, cluster_info) {
    if (is.null(cluster_pairs_df) || nrow(cluster_pairs_df) == 0) return(data.frame())
    
    # Group stats comme avant
    stats <- cluster_pairs_df %>%
      group_by(ref_cluster, target_cluster) %>%
      summarise(
        n_pairs = n(),
        mean_distance_um   = round(mean(distance_um,   na.rm = TRUE), 1),
        median_distance_um = round(median(distance_um, na.rm = TRUE), 1),
        min_distance_um    = round(min(distance_um,    na.rm = TRUE), 1),
        max_distance_um    = round(max(distance_um,    na.rm = TRUE), 1),
        same_cluster       = (ref_cluster[1] == target_cluster[1]),
        .groups = "drop"
      )
    
    # ADD: Calculate percentages for each direction
    # Get total spots per cluster
    g1_by_cluster <- table(cluster_info$spots[g1_spots])
    g2_by_cluster <- table(cluster_info$spots[g2_spots])
    
    # For each row, calculate what % of gene1 in ref_cluster are close to gene2 in target_cluster
    stats$gene1_pct <- sapply(1:nrow(stats), function(i) {
      ref_cl <- stats$ref_cluster[i]
      # Count unique gene1 spots from ref_cluster that have neighbors in target_cluster
      unique_gene1_spots <- length(unique(cluster_pairs_df$ref_spot[
        cluster_pairs_df$ref_cluster == ref_cl & 
          cluster_pairs_df$target_cluster == stats$target_cluster[i]
      ]))
      total_gene1 <- as.numeric(g1_by_cluster[ref_cl])
      if (is.na(total_gene1) || total_gene1 == 0) return(0)
      round(100 * unique_gene1_spots / total_gene1, 1)
    })
    
    stats <- stats %>%
      arrange(desc(n_pairs)) %>%
      select(ref_cluster, target_cluster, n_pairs, gene1_pct, 
             mean_distance_um, median_distance_um, min_distance_um, max_distance_um, same_cluster)
    
    return(stats)
  }
  
  # Résumé par cluster
  calculate_cluster_summary <- function(g1_spots, g2_spots, g1_to_g2, g2_to_g1, cluster_info) {
    g1_by_cl <- table(cluster_info$spots[g1_spots])
    g2_by_cl <- table(cluster_info$spots[g2_spots])
    all_cl   <- cluster_info$names
    df <- data.frame(
      cluster    = all_cl,
      gene1_spots= as.numeric(g1_by_cl[all_cl]),
      gene2_spots= as.numeric(g2_by_cl[all_cl]),
      stringsAsFactors = FALSE
    )
    df[is.na(df)] <- 0
    add_unique_spots <- function(df, pairs_df, ref_tag) {
      if (!is.null(pairs_df) && nrow(pairs_df) > 0) {
        for (cl in df$cluster) {
          sub <- pairs_df[pairs_df$ref_cluster == cl, ]
          if (nrow(sub) > 0) {
            unique_spots_same <- length(unique(sub$ref_spot[sub$same_cluster]))
            unique_spots_diff <- length(unique(sub$ref_spot[!sub$same_cluster]))
            
            df[df$cluster == cl, paste0(ref_tag, "_coloc_same_cluster")] <- unique_spots_same
            df[df$cluster == cl, paste0(ref_tag, "_coloc_different_cluster")] <- unique_spots_diff
          } else {
            df[df$cluster == cl, paste0(ref_tag, "_coloc_same_cluster")] <- 0
            df[df$cluster == cl, paste0(ref_tag, "_coloc_different_cluster")] <- 0
          }
        }
      } else {
        df[, paste0(ref_tag, "_coloc_same_cluster")] <- 0
        df[, paste0(ref_tag, "_coloc_different_cluster")] <- 0
      }
      df
    }
    df <- add_unique_spots(df, g1_to_g2$cluster_pairs, "gene1")
    df <- add_unique_spots(df, if (is.null(g2_to_g1)) NULL else g2_to_g1$cluster_pairs, "gene2")
    df$gene1_total_coloc <- df$gene1_coloc_same_cluster + df$gene1_coloc_different_cluster
    df$gene2_total_coloc <- df$gene2_coloc_same_cluster + df$gene2_coloc_different_cluster
    safe_pct <- function(num, den) ifelse(den > 0, round(100 * num / den, 1), 0)
    df$gene1_coloc_rate_same <- safe_pct(df$gene1_coloc_same_cluster, df$gene1_spots)
    df$gene1_coloc_rate_diff <- safe_pct(df$gene1_coloc_different_cluster, df$gene1_spots)
    df$gene1_total_coloc_rate<- safe_pct(df$gene1_total_coloc, df$gene1_spots)
    df$gene2_coloc_rate_same <- safe_pct(df$gene2_coloc_same_cluster, df$gene2_spots)
    df$gene2_coloc_rate_diff <- safe_pct(df$gene2_coloc_different_cluster, df$gene2_spots)
    df$gene2_total_coloc_rate<- safe_pct(df$gene2_total_coloc, df$gene2_spots)
    
    df
  }
  
  analyze_gene_colocalization <- function(sp_obj, gene_list,
                                          max_distance_um = 50,
                                          bin_size_um = 8,
                                          assay = "Spatial",
                                          bidirectional = TRUE) {
    max_distance_um <- as.numeric(max_distance_um)
    bin_size_um <- as.numeric(bin_size_um)
    n_genes <- length(gene_list)
    if(n_genes == 2) {
      gene1 <- gene_list[1]
      gene2 <- gene_list[2]
      message(sprintf("=== Co-localization: %s vs %s | %.0f µm, %.0f µm/bin ===",
                      gene1, gene2, max_distance_um, bin_size_um))
      if (length(sp_obj@images) == 0) stop("No spatial image / coordinates.")
      if (!assay %in% names(sp_obj@assays)) stop("Assay not found.")
      if (!gene1 %in% rownames(sp_obj[[assay]])) stop(sprintf("Gene %s not found", gene1))
      if (!gene2 %in% rownames(sp_obj[[assay]])) stop(sprintf("Gene %s not found", gene2))
      coords <- Seurat::GetTissueCoordinates(sp_obj@images[[1]])
      if (is.null(coords) || nrow(coords) == 0) stop("Empty tissue coordinates.")
      coords <- as.data.frame(coords)
      expr <- Seurat::GetAssayData(sp_obj, assay = assay, slot = "data")
      g1 <- expr[gene1, ]; g2 <- expr[gene2, ]
      g1_spots <- names(g1)[g1 > 0]
      g2_spots <- names(g2)[g2 > 0]
      g1_spots <- intersect(g1_spots, rownames(coords))
      g2_spots <- intersect(g2_spots, rownames(coords))
      max_dist_bins <- convert_um_to_bins(max_distance_um, bin_size_um)
      cl_info <- get_cluster_info(sp_obj)
      calc_pairs <- function(ref_spots, tgt_spots, ref_gene, tgt_gene) {
        if (length(ref_spots) == 0 || length(tgt_spots) == 0)
          return(list(count = 0, pairs = list(), cluster_pairs = data.frame()))
        pairs_list <- list(); out_rows <- list(); ref_with_any <- 0
        for (rs in ref_spots) {
          rx <- coords[rs, 1]; ry <- coords[rs, 2]; rcl <- cl_info$spots[rs]
          tx <- coords[tgt_spots, 1]; ty <- coords[tgt_spots, 2]
          d  <- sqrt((rx - tx)^2 + (ry - ty)^2)
          keep <- which(d <= max_dist_bins & tgt_spots != rs)
          if (length(keep)) {
            ref_with_any <- ref_with_any + 1
            pairs_list[[rs]] <- tgt_spots[keep]
            for (k in keep) {
              tspot <- tgt_spots[k]; tcl <- cl_info$spots[tspot]
              out_rows[[length(out_rows)+1]] <- data.frame(
                ref_spot = rs,
                ref_cluster = rcl,
                ref_gene = ref_gene,
                target_spot = tspot,
                target_cluster = tcl,
                target_gene = tgt_gene,
                same_cluster = rcl == tcl,
                distance_bins = round(d[k], 2),
                distance_um   = round(d[k] * bin_size_um, 1),
                stringsAsFactors = FALSE
              )
            }
          }
        }
        cp <- if (length(out_rows)) do.call(rbind, out_rows) else data.frame()
        list(count = ref_with_any, pairs = pairs_list, cluster_pairs = cp)
      }
      g1_to_g2 <- calc_pairs(g1_spots, g2_spots, gene1, gene2)
      g2_to_g1 <- if (isTRUE(bidirectional)) calc_pairs(g2_spots, g1_spots, gene2, gene1) else NULL
      cl_sum <- calculate_cluster_summary(g1_spots, g2_spots, g1_to_g2, g2_to_g1, cl_info)
      inter12 <- calculate_inter_cluster_stats(g1_to_g2$cluster_pairs, g1_spots, g2_spots, cl_info)
      inter21 <- if (!is.null(g2_to_g1)) calculate_inter_cluster_stats(g2_to_g1$cluster_pairs, g2_spots, g1_spots, cl_info) else data.frame()
      spatial_df <- data.frame(
        spot_id = rownames(coords),
        x = coords[, 1], y = coords[, 2],
        gene1_expr = g1[rownames(coords)],
        gene2_expr = g2[rownames(coords)],
        gene1_positive = rownames(coords) %in% g1_spots,
        gene2_positive = rownames(coords) %in% g2_spots,
        gene1_has_neighbor = rownames(coords) %in% names(g1_to_g2$pairs),
        gene2_has_neighbor = if (!is.null(g2_to_g1)) rownames(coords) %in% names(g2_to_g1$pairs) else FALSE,
        cluster = cl_info$spots[rownames(coords)],
        stringsAsFactors = FALSE
      )
      safe_pct <- function(num, den) ifelse(den > 0, round(100 * num / den, 1), NA)
      return(list(
        analysis_type = "bidirectional",
        parameters = list(
          gene1 = gene1, gene2 = gene2,
          max_distance_um = max_distance_um,
          bin_size_um = bin_size_um,
          assay = assay, bidirectional = bidirectional
        ),
        gene1_total_spots = length(g1_spots),
        gene2_total_spots = length(g2_spots),
        gene1_near_gene2  = g1_to_g2$count,
        gene1_coloc_rate  = safe_pct(g1_to_g2$count, length(g1_spots)),
        gene2_near_gene1  = if (!is.null(g2_to_g1)) g2_to_g1$count else NA,
        gene2_coloc_rate  = if (!is.null(g2_to_g1)) safe_pct(g2_to_g1$count, length(g2_spots)) else NA,
        cluster_pairs_gene1_to_gene2 = g1_to_g2$cluster_pairs,
        cluster_pairs_gene2_to_gene1 = if (!is.null(g2_to_g1)) g2_to_g1$cluster_pairs else data.frame(),
        cluster_summary = cl_sum,
        inter_cluster_gene1_to_gene2 = inter12,
        inter_cluster_gene2_to_gene1 = inter21,
        spatial_data = spatial_df
      ))
    } else {
      message(sprintf("=== Multi-gene co-localization: %s | %.0f µm, %.0f µm/bin ===",
                      paste(gene_list, collapse = " + "), max_distance_um, bin_size_um))
      if (length(sp_obj@images) == 0) stop("No spatial image / coordinates.")
      if (!assay %in% names(sp_obj@assays)) stop("Assay not found.")
      for(g in gene_list) {
        if (!g %in% rownames(sp_obj[[assay]])) stop(sprintf("Gene %s not found", g))
      }
      coords <- Seurat::GetTissueCoordinates(sp_obj@images[[1]])
      if (is.null(coords) || nrow(coords) == 0) stop("Empty tissue coordinates.")
      coords <- as.data.frame(coords)
      expr <- Seurat::GetAssayData(sp_obj, assay = assay, slot = "data")
      max_dist_bins <- convert_um_to_bins(max_distance_um, bin_size_um)
      gene_spots <- lapply(gene_list, function(g) {
        spots <- names(expr[g, ])[expr[g, ] > 0]
        intersect(spots, rownames(coords))
      })
      names(gene_spots) <- gene_list
      colocalized_counts <- sapply(gene_list, function(gene) {
        other_genes <- setdiff(gene_list, gene)
        count <- 0
        
        for(spot in gene_spots[[gene]]) {
          all_present <- TRUE
          for(other_gene in other_genes) {
            if(length(gene_spots[[other_gene]]) == 0) {
              all_present <- FALSE
              break
            }
            distances <- sqrt((coords[spot, 1] - coords[gene_spots[[other_gene]], 1])^2 + 
                                (coords[spot, 2] - coords[gene_spots[[other_gene]], 2])^2)
            if(!any(distances <= max_dist_bins)) {
              all_present <- FALSE
              break
            }
          }
          if(all_present) count <- count + 1
        }
        return(count)
      })
      
      # CORRECTION : Calculer les pourcentages correctement
      total_spots_per_gene <- sapply(gene_list, function(g) length(gene_spots[[g]]))
      colocalization_rates <- sapply(gene_list, function(g) {
        total <- total_spots_per_gene[g]
        coloc <- colocalized_counts[g]
        if(total > 0 && !is.na(total) && !is.na(coloc)) {
          round(100 * coloc / total, 1)
        } else {
          0
        }
      })
      
      # Format de retour pour multi-gènes - CORRECTION
      return(list(
        analysis_type = "simultaneous",
        parameters = list(
          gene_list = gene_list,
          max_distance_um = max_distance_um,
          bin_size_um = bin_size_um,
          assay = assay
        ),
        total_spots_per_gene = total_spots_per_gene,
        colocalized_spots_per_gene = colocalized_counts,
        colocalization_rates = colocalization_rates  # CORRECTION : maintenant calculé correctement
      ))
    }
  }
  
  # Lancement analyse
  observeEvent(input$run_coloc_analysis, {
    req(spatial_obj(), input$coloc_genes)
    
    # Parser les gènes séparés par virgules
    gene_list <- strsplit(input$coloc_genes, ",")[[1]]
    gene_list <- trimws(gene_list)
    gene_list <- gene_list[gene_list != ""]
    
    if(length(gene_list) < 2) {
      showNotification("Please enter at least 2 genes", type = "warning")
      return()
    }
    
    # Vérifier existence
    available_genes <- rownames(spatial_obj()[[input$coloc_assay]])
    missing_genes <- gene_list[!gene_list %in% available_genes]
    
    if(length(missing_genes) > 0) {
      showNotification(paste("Genes not found:", paste(missing_genes, collapse = ", ")), type = "error")
      return()
    }
    
    tryCatch({
      showModal(modalDialog(title = "Running Co-localization Analysis", p("Analyzing..."), easyClose = FALSE, footer = NULL))
      
      results <- analyze_gene_colocalization(
        sp_obj = spatial_obj(),  # CHANGÉ : sp_obj au lieu de spatial_obj
        gene_list = gene_list,
        max_distance_um = input$coloc_distance_um,
        bin_size_um = input$spatial_resolution,
        assay = input$coloc_assay
      )
      
      coloc_results(results)
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Flag résultats prêts
  output$coloc_results_ready <- reactive({ !is.null(coloc_results()) })
  outputOptions(output, "coloc_results_ready", suspendWhenHidden = FALSE)
  
  
  # Quick stats summary - CLEAR text with full explanations
  output$coloc_quick_stats <- renderText({
    req(coloc_results())
    r <- coloc_results()
    if (r$analysis_type == "simultaneous") {
      genes <- r$parameters$gene_list
      threshold <- r$parameters$max_distance_um
      bin_size <- r$parameters$bin_size_um
      header <- sprintf("=== Multi-gene Co-localization Analysis ===\nGenes: %s\nDistance threshold: %d µm | Bin size: %d µm\n\n",
                        paste(genes, collapse = " + "), threshold, bin_size)
      gene_lines <- sapply(genes, function(g) {
        total <- r$total_spots_per_gene[g]
        coloc <- r$colocalized_spots_per_gene[g]
        rate <- if (total > 0) round(100 * coloc / total, 1) else 0
        sprintf("%s: %d total spots, %d co-localized (%0.1f%% of %s spots are near all other genes within %d µm)",
                g, total, coloc, rate, g, threshold)
      })
      paste(c(header, gene_lines), collapse = "\n")
    } else {
      params <- r$parameters
      gene1 <- params$gene1
      gene2 <- params$gene2
      threshold <- params$max_distance_um
      bin_size <- params$bin_size_um
      header <- sprintf("=== Co-localization Analysis ===\nGenes: %s vs %s\nDistance threshold: %d µm | Bin size: %d µm\n",
                        gene1, gene2, threshold, bin_size)
      counts <- sprintf("\n%s: %s spots\n%s: %s spots\n",
                        gene1, format(r$gene1_total_spots, big.mark = ","),
                        gene2, format(r$gene2_total_spots, big.mark = ","))
      gene1_rate <- if (!is.na(r$gene1_coloc_rate)) r$gene1_coloc_rate else 0
      gene1_text <- sprintf("\n%0.1f%% of %s spots have at least one %s neighbor within %d µm\n(%d out of %s %s spots)",
                            gene1_rate, gene1, gene2, threshold,
                            r$gene1_near_gene2, format(r$gene1_total_spots, big.mark = ","), gene1)
      gene2_text <- ""
      if (isTRUE(params$bidirectional) && !is.na(r$gene2_coloc_rate)) {
        gene2_rate <- r$gene2_coloc_rate
        gene2_text <- sprintf("\n\n%0.1f%% of %s spots have at least one %s neighbor within %d µm\n(%d out of %s %s spots)",
                              gene2_rate, gene2, gene1, threshold,
                              r$gene2_near_gene1, format(r$gene2_total_spots, big.mark = ","), gene2)
      }
      
      paste0(header, counts, gene1_text, gene2_text)
    }
  })
  
  # Cluster-level co-localization table with visual bars and clear column names
  output$coloc_cluster_table <- renderDT({
    req(coloc_results())
    r <- coloc_results()
    if(r$analysis_type == "simultaneous") {
      genes <- r$parameters$gene_list
      percentages <- sapply(genes, function(g) {
        total <- r$total_spots_per_gene[g]
        coloc <- r$colocalized_spots_per_gene[g] 
        if(total > 0 && !is.na(total) && !is.na(coloc)) {
          round(100 * coloc / total, 1)
        } else {
          0
        }
      })
      summary_df <- data.frame(
        Gene = genes,
        Total_spots = r$total_spots_per_gene,
        Colocalized_spots = r$colocalized_spots_per_gene,
        Coloc_bar = sapply(percentages, create_expression_percentage_bar),
        stringsAsFactors = FALSE
      )
      datatable(
        summary_df, 
        options = list(
          pageLength = 10, 
          scrollX = TRUE,
          columnDefs = list(
            list(className = 'dt-center', targets = '_all'),
            list(width = '150px', targets = 3)
          )
        ),
        rownames = FALSE,
        caption = paste("Multi-gene co-localization:", paste(genes, collapse = " + ")),
        escape = FALSE,
        colnames = c('Gene', 'Total Spots', 'Co-localized Spots', 'Co-localization Rate')
      )
    } else {
      params <- r$parameters
      cs <- r$cluster_summary
      if (is.null(cs) || nrow(cs) == 0) {
        return(datatable(data.frame(Message = "No cluster information available")))
      }
      gene1 <- params$gene1
      gene2 <- params$gene2
      max_dist_threshold <- params$max_distance_um
      inter <- r$inter_cluster_gene1_to_gene2
      intra_distances <- NULL
      if (!is.null(inter) && nrow(inter) > 0 && "ref_cluster" %in% names(inter)) {
        intra_distances <- inter %>%
          filter(same_cluster == TRUE) %>%
          select(
            cluster = ref_cluster,
            intra_mean_dist = mean_distance_um,
            intra_median_dist = median_distance_um
          )
      }
      out <- cs %>%
        transmute(
          Cluster = cluster,
          Gene1_total = gene1_spots,
          Gene1_same = gene1_coloc_same_cluster,
          Gene1_diff = gene1_coloc_different_cluster,
          Gene1_pct_same = gene1_coloc_rate_same,
          Gene1_pct_diff = gene1_coloc_rate_diff,
          Gene1_bar_same = sapply(gene1_coloc_rate_same, create_expression_percentage_bar),
          Gene1_bar_diff = sapply(gene1_coloc_rate_diff, create_expression_percentage_bar),
          Gene2_total = gene2_spots,
          Gene2_same = gene2_coloc_same_cluster,
          Gene2_diff = gene2_coloc_different_cluster,
          Gene2_pct_same = gene2_coloc_rate_same,
          Gene2_pct_diff = gene2_coloc_rate_diff,
          Gene2_bar_same = sapply(gene2_coloc_rate_same, create_expression_percentage_bar),
          Gene2_bar_diff = sapply(gene2_coloc_rate_diff, create_expression_percentage_bar)
        )
      if (!is.null(intra_distances) && nrow(intra_distances) > 0) {
        out <- out %>%
          left_join(intra_distances, by = c("Cluster" = "cluster"))
        out <- out %>%
          mutate(
            intra_mean_bar = sapply(intra_mean_dist, function(d) {
              create_distance_bar(d, max_dist_threshold)  # No need for round() or ifelse
            }),
            intra_median_bar = sapply(intra_median_dist, function(d) {
              create_distance_bar(d, max_dist_threshold)  # Function handles NA internally
            })
          )
      } else {
        out <- out %>%
          mutate(
            intra_mean_dist = NA,
            intra_median_dist = NA,
            intra_mean_bar = "<span style='color: grey;'>N/A</span>",
            intra_median_bar = "<span style='color: grey;'>N/A</span>"
          )
      }
      colnames(out) <- gsub("Gene1", gene1, colnames(out))
      colnames(out) <- gsub("Gene2", gene2, colnames(out))
      hidden_cols <- c(
        paste0(gene1, "_pct_same"), paste0(gene1, "_pct_diff"),
        paste0(gene2, "_pct_same"), paste0(gene2, "_pct_diff"),
        "intra_mean_dist", "intra_median_dist"
      )
      hidden_indices <- which(colnames(out) %in% hidden_cols) - 1
      col_display <- colnames(out)
      col_display[col_display == "Cluster"] <- "Cluster"
      col_display[col_display == paste0(gene1, "_total")] <- paste(gene1, "Total")
      col_display[col_display == paste0(gene1, "_same")] <- paste(gene1, "Same")
      col_display[col_display == paste0(gene1, "_diff")] <- paste(gene1, "Diff")
      col_display[col_display == paste0(gene1, "_bar_same")] <- paste(gene1, "Within %")
      col_display[col_display == paste0(gene1, "_bar_diff")] <- paste(gene1, "Across %")
      col_display[col_display == paste0(gene2, "_total")] <- paste(gene2, "Total")
      col_display[col_display == paste0(gene2, "_same")] <- paste(gene2, "Same")
      col_display[col_display == paste0(gene2, "_diff")] <- paste(gene2, "Diff")
      col_display[col_display == paste0(gene2, "_bar_same")] <- paste(gene2, "Within %")
      col_display[col_display == paste0(gene2, "_bar_diff")] <- paste(gene2, "Across %")
      col_display[col_display == "intra_mean_bar"] <- "Mean Dist"
      col_display[col_display == "intra_median_bar"] <- "Median Dist"
      bar_cols <- grep("_bar_|_bar$", colnames(out))
      number_cols <- which(colnames(out) %in% c(
        "Cluster",
        paste0(gene1, "_total"), paste0(gene1, "_same"), paste0(gene1, "_diff"),
        paste0(gene2, "_total"), paste0(gene2, "_same"), paste0(gene2, "_diff")
      ))
      datatable(
        out, 
        options = list(
          pageLength = 15, 
          scrollX = TRUE,
          columnDefs = list(
            list(visible = FALSE, targets = hidden_indices),
            list(className = 'dt-center', targets = '_all'),
            # Bar columns = WIDE (150px)
            list(width = '150px', targets = bar_cols - 1),
            # Number columns = NARROW (70px)
            list(width = '70px', targets = number_cols - 1)
          ),
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ), 
        rownames = FALSE,
        caption = paste0(
          "Cluster co-localization: ", gene1, " ↔ ", gene2,
          " (≤ ", max_dist_threshold, " µm)"
        ),
        escape = FALSE,
        colnames = col_display
      ) %>%
        formatRound(
          columns = c('intra_mean_dist', 'intra_median_dist'), 
          digits = 1
        )
    }
  })
  
  # Update cluster pair choices
  observe({
    req(coloc_results())
    inter <- coloc_results()$inter_cluster_gene1_to_gene2
    
    if (!is.null(inter) && nrow(inter) > 0 && "ref_cluster" %in% names(inter)) {
      pair_labels <- paste(inter$ref_cluster, "→", inter$target_cluster)
      names(pair_labels) <- pair_labels
      
      updateSelectInput(session, "inter_cluster_pairs_select", 
                        choices = pair_labels,
                        selected = NULL)
    }
  })
  
  # Inter-cluster table - SIMPLIFIED
  output$coloc_inter_cluster_table <- renderDT({
    req(coloc_results())
    r <- coloc_results()
    inter_g1_to_g2 <- r$inter_cluster_gene1_to_gene2
    inter_g2_to_g1 <- r$inter_cluster_gene2_to_gene1
    
    if (is.null(inter_g1_to_g2) || nrow(inter_g1_to_g2) == 0) {
      return(datatable(data.frame(Message = "No inter-cluster interactions detected")))
    }
    
    if (!"gene1_pct" %in% colnames(inter_g1_to_g2)) {
      return(datatable(data.frame(Message = "Error: Percentage column not found. Please re-run colocalization analysis.")))
    }
    
    # Filter out intra-cluster
    g1g2 <- inter_g1_to_g2 %>% filter(ref_cluster != target_cluster)
    g2g1 <- if (!is.null(inter_g2_to_g1) && "gene1_pct" %in% colnames(inter_g2_to_g1)) {
      inter_g2_to_g1 %>% filter(ref_cluster != target_cluster)
    } else {
      data.frame()
    }
    
    if (nrow(g1g2) == 0) {
      return(datatable(data.frame(Message = "No inter-cluster interactions")))
    }
    
    # Create unified table
    unique_pairs <- unique(g1g2[, c("ref_cluster", "target_cluster")])
    
    combined <- lapply(1:nrow(unique_pairs), function(i) {
      clA <- unique_pairs$ref_cluster[i]
      clB <- unique_pairs$target_cluster[i]
      
      row_AB_g1g2 <- g1g2[g1g2$ref_cluster == clA & g1g2$target_cluster == clB, ]
      row_BA_g1g2 <- g1g2[g1g2$ref_cluster == clB & g1g2$target_cluster == clA, ]
      row_AB_g2g1 <- if (nrow(g2g1) > 0) g2g1[g2g1$ref_cluster == clA & g2g1$target_cluster == clB, ] else data.frame()
      row_BA_g2g1 <- if (nrow(g2g1) > 0) g2g1[g2g1$ref_cluster == clB & g2g1$target_cluster == clA, ] else data.frame()
      
      data.frame(
        clusterA = clA,
        clusterB = clB,
        n_pairs_AB = if(nrow(row_AB_g1g2) > 0) row_AB_g1g2$n_pairs[1] else 0,
        g1A_to_g2B_pct = if(nrow(row_AB_g1g2) > 0) row_AB_g1g2$gene1_pct[1] else 0,
        g2B_to_g1A_pct = if(nrow(row_BA_g2g1) > 0) row_BA_g2g1$gene1_pct[1] else 0,
        n_pairs_BA = if(nrow(row_BA_g1g2) > 0) row_BA_g1g2$n_pairs[1] else 0,
        g1B_to_g2A_pct = if(nrow(row_BA_g1g2) > 0) row_BA_g1g2$gene1_pct[1] else 0,
        g2A_to_g1B_pct = if(nrow(row_AB_g2g1) > 0) row_AB_g2g1$gene1_pct[1] else 0,
        mean_dist = if(nrow(row_AB_g1g2) > 0) row_AB_g1g2$mean_distance_um[1] else NA,
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()
    
    if (!isTRUE(input$inter_cluster_show_all) && 
        !is.null(input$inter_cluster_pairs_select) && 
        length(input$inter_cluster_pairs_select) > 0) {
      selected_pairs <- input$inter_cluster_pairs_select
      source_targets <- strsplit(selected_pairs, " → ")
      sources <- sapply(source_targets, function(x) x[1])
      targets <- sapply(source_targets, function(x) x[2])
      keep_rows <- rep(FALSE, nrow(combined))
      for (i in seq_along(sources)) {
        keep_rows <- keep_rows | (combined$clusterA == sources[i] & combined$clusterB == targets[i])
      }
      combined <- combined[keep_rows, ]
    }
    
    if (nrow(combined) == 0) {
      return(datatable(data.frame(Message = "No interactions match filters")))
    }
    
    gene1_name <- r$parameters$gene1
    gene2_name <- r$parameters$gene2
    
    # Use EXISTING create_expression_percentage_bar function
    display_df <- combined %>%
      mutate(
        Pair = paste(clusterA, "↔", clusterB),
        Bar1 = sapply(g1A_to_g2B_pct, create_expression_percentage_bar),
        Bar2 = sapply(g2B_to_g1A_pct, create_expression_percentage_bar),
        Bar3 = sapply(g1B_to_g2A_pct, create_expression_percentage_bar),
        Bar4 = sapply(g2A_to_g1B_pct, create_expression_percentage_bar)
      ) %>%
      select(Pair, n_pairs_AB, n_pairs_BA, Bar1, Bar2, Bar3, Bar4, mean_dist)
    
    colnames(display_df) <- c(
      "Cluster Pair", "n A→B", "n B→A",
      paste0(gene1_name, "(A)→", gene2_name, "(B)"),
      paste0(gene2_name, "(B)→", gene1_name, "(A)"),
      paste0(gene1_name, "(B)→", gene2_name, "(A)"),
      paste0(gene2_name, "(A)→", gene1_name, "(B)"),
      "Mean Dist (µm)"
    )
    
    datatable(
      display_df,
      options = list(
        pageLength = 20, 
        scrollX = TRUE,
        columnDefs = list(list(width = '150px', targets = c(3, 4, 5, 6))),
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      rownames = FALSE,
      caption = paste0("Inter-cluster interactions (≤ ", r$parameters$max_distance_um, " µm)"),
      escape = FALSE
    ) %>%
      formatRound(columns = ncol(display_df), digits = 1)
  })
  
  
  # Update cluster pair choices - FIXED: exclude same→same and ensure all pairs shown
  observe({
    req(coloc_results())
    inter <- coloc_results()$inter_cluster_gene1_to_gene2
    
    if (!is.null(inter) && nrow(inter) > 0 && "ref_cluster" %in% names(inter)) {
      # Filter out intra-cluster pairs (same → same)
      inter_only <- inter %>%
        filter(ref_cluster != target_cluster)
      
      if (nrow(inter_only) > 0) {
        # Create unique pair labels
        pair_df <- inter_only %>%
          select(ref_cluster, target_cluster) %>%
          distinct() %>%
          arrange(ref_cluster, target_cluster)
        
        # Create labels
        pair_labels <- paste(pair_df$ref_cluster, "→", pair_df$target_cluster)
        names(pair_labels) <- pair_labels
        
        updateSelectInput(session, "inter_cluster_pairs_select", 
                          choices = pair_labels,
                          selected = NULL)
      } else {
        # No inter-cluster pairs, only intra-cluster
        updateSelectInput(session, "inter_cluster_pairs_select", 
                          choices = c("No inter-cluster pairs found" = ""),
                          selected = NULL)
      }
    }
  })
  
  # Add reset button handler
  observeEvent(input$reset_cluster_filter, {
    updateCheckboxInput(session, "inter_cluster_show_all", value = TRUE)
    updateSelectInput(session, "inter_cluster_pairs_select", selected = NULL)
  })
  # Cartes spatiales
  # Plot 1: Gene expression map - IMPROVED with better colors
  output$coloc_spatial_plot <- renderPlot({
    req(coloc_results())
    r <- coloc_results()
    
    if(r$analysis_type == "simultaneous") {
      text_plot <- ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                 label = paste("Multi-gene visualization\n\nGenes analyzed:",
                               paste(r$parameters$gene_list, collapse = ", ")),
                 size = 6, hjust = 0.5) +
        theme_void() +
        labs(title = "Multi-gene Co-localization")
      return(text_plot)
    } else {
      d <- r$spatial_data
      p <- r$parameters
      
      # Categorize spots with clear colors
      d$category <- "No expression"
      d$category[d$gene1_positive & !d$gene2_positive] <- "Gene1 only"
      d$category[!d$gene1_positive & d$gene2_positive] <- "Gene2 only"
      d$category[d$gene1_positive & d$gene2_positive] <- "Both genes"
      
      ggplot(d, aes(x = x, y = y)) +
        geom_point(aes(color = category, alpha = category, size = category)) +
        scale_color_manual(
          values = c(
            "No expression" = "grey90",
            "Gene1 only" = "dodgerblue2",
            "Gene2 only" = "firebrick2",
            "Both genes" = "purple3"
          ),
          labels = c("No expression", 
                     paste(p$gene1, "only"), 
                     paste(p$gene2, "only"), 
                     "Both genes")
        ) +
        scale_alpha_manual(
          values = c("No expression" = 0.2, "Gene1 only" = 0.8, 
                     "Gene2 only" = 0.8, "Both genes" = 1),
          guide = "none"
        ) +
        scale_size_manual(
          values = c("No expression" = 0.3, "Gene1 only" = 1.2, 
                     "Gene2 only" = 1.2, "Both genes" = 1.5),
          guide = "none"
        ) +
        scale_y_reverse() + 
        coord_fixed() +
        labs(
          title = paste("Gene Expression:", p$gene1, "vs", p$gene2),
          x = "X coordinate", 
          y = "Y coordinate", 
          color = "Expression Status"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "bottom"
        )
    }
  }, height = 500)
  
  # Plot 2: Co-localization overlay - IMPROVED colors
  output$coloc_overlay_plot <- renderPlot({
    req(coloc_results())
    r <- coloc_results()
    
    if(r$analysis_type == "simultaneous") {
      genes <- r$parameters$gene_list
      
      plot_data <- data.frame(
        Gene = rep(genes, 2),
        Type = rep(c("Total spots", "Co-localized spots"), each = length(genes)),
        Count = c(r$total_spots_per_gene, r$colocalized_spots_per_gene)
      )
      
      ggplot(plot_data, aes(x = Gene, y = Count, fill = Type)) +
        geom_col(position = "dodge", width = 0.7, alpha = 0.8) +
        scale_fill_manual(values = c("Total spots" = "steelblue", 
                                     "Co-localized spots" = "forestgreen")) +
        labs(
          title = "Spot Counts by Gene",
          subtitle = paste("Distance ≤", r$parameters$max_distance_um, "µm"),
          x = "Gene", y = "Number of spots", fill = "Spot type"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        scale_y_continuous(labels = scales::comma)
      
    } else {
      d <- r$spatial_data
      p <- r$parameters
      
      # Show ALL spots as light background, then highlight co-localized
      ggplot(d, aes(x = x, y = y)) +
        # Background: all spots in light color
        geom_point(color = "grey85", size = 0.5, alpha = 0.4) +
        # Highlight: Gene1 with nearby Gene2
        {if (any(d$gene1_has_neighbor))
          geom_point(data = d[d$gene1_has_neighbor, ], 
                     aes(x = x, y = y), 
                     color = "dodgerblue3", size = 2, alpha = 0.9)} +
        # Highlight: Gene2 with nearby Gene1 (if bidirectional)
        {if (isTRUE(p$bidirectional) && any(d$gene2_has_neighbor))
          geom_point(data = d[d$gene2_has_neighbor, ], 
                     aes(x = x, y = y), 
                     color = "firebrick3", size = 1.6, alpha = 0.9)} +
        scale_y_reverse() + 
        coord_fixed() +
        labs(
          title = paste("Co-localized Spots (≤", p$max_distance_um, "µm)"),
          subtitle = paste0(
            "Blue: ", sum(d$gene1_has_neighbor), " ", p$gene1, " spots with nearby ", p$gene2,
            if (isTRUE(p$bidirectional)) {
              paste0(" | Red: ", sum(d$gene2_has_neighbor), " ", p$gene2, " spots with nearby ", p$gene1)
            } else ""
          ),
          x = "X coordinate", 
          y = "Y coordinate"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11)
        )
    }
  }, height = 500)
  
  # Plot 3: H&E with co-localized spots - IMPROVED visibility
  # Plot 3: Spatial map with NAMED clusters - FIXED colors
  output$coloc_histology_overlay <- renderPlot({
    req(coloc_results(), spatial_obj())
    r <- coloc_results()
    sp <- spatial_obj()
    
    tryCatch({
      # Get spatial data
      d <- r$spatial_data
      
      # Check coordinates
      if (!all(c("x", "y") %in% names(d))) {
        stop("No spatial coordinates in data")
      }
      
      # Get cluster information - PRIORITIZE NAMED CLUSTERS
      if ("annotated_clusters" %in% colnames(sp@meta.data)) {
        clusters <- sp@meta.data$annotated_clusters[match(d$spot_id, colnames(sp))]
        d$cluster <- as.character(clusters)
      } else if ("seurat_clusters" %in% colnames(sp@meta.data)) {
        clusters <- sp@meta.data$seurat_clusters[match(d$spot_id, colnames(sp))]
        d$cluster <- as.character(clusters)
      } else if ("cluster" %in% names(d)) {
        d$cluster <- as.character(d$cluster)
      } else {
        d$cluster <- "No cluster"
      }
      
      # Identify co-localized spots
      if (r$analysis_type == "simultaneous") {
        d$is_colocalized <- FALSE
      } else {
        d$is_colocalized <- d$gene1_has_neighbor | d$gene2_has_neighbor
      }
      
      # Get spot size from slider
      spot_size <- if (!is.null(input$coloc_spot_size)) {
        input$coloc_spot_size
      } else {
        1.5
      }
      
      # Create base plot with clusters - NORMAL ggplot2 colors
      p <- ggplot(d, aes(x = x, y = y)) +
        geom_point(aes(color = cluster), size = 1, alpha = 0.5)
      
      # Add co-localized spots only if checkbox is checked
      if (isTRUE(input$show_coloc_spots) && any(d$is_colocalized)) {
        coloc_spots <- d[d$is_colocalized, ]
        
        # Glow effect
        p <- p + 
          geom_point(
            data = coloc_spots,
            aes(x = x, y = y),
            color = "yellow",
            size = spot_size + 0.5,
            alpha = 0.3
          ) +
          # Center point
          geom_point(
            data = coloc_spots,
            aes(x = x, y = y),
            color = "red",
            size = spot_size,
            alpha = 0.9
          )
      }
      
      # Formatting - NO scale_color_manual to keep default ggplot colors
      p <- p +
        scale_y_reverse() +
        coord_fixed() +
        labs(
          title = if (r$analysis_type == "simultaneous") {
            paste("Spatial Map:", paste(r$parameters$gene_list, collapse = " + "))
          } else {
            paste("Spatial Map:", r$parameters$gene1, "↔", r$parameters$gene2)
          },
          subtitle = if (isTRUE(input$show_coloc_spots)) {
            paste0(sum(d$is_colocalized), " co-localized spots highlighted in red")
          } else {
            "Co-localized spots hidden (use checkbox to show)"
          },
          x = "X coordinate",
          y = "Y coordinate",
          color = "Cluster"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 11),
          legend.position = "right",
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 9)
        )
      
      return(p)
      
    }, error = function(e) {
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                 label = paste("Cannot create overlay plot\n", e$message),
                 size = 5, color = "red") +
        theme_void()
    })
  }, height = 600)
  
  ############################## Marker Analysis ##############################
  
  # Update cluster choices when clusters are available
  observe({
    req(spatial_obj())
    
    if (processing_states$clustered && "seurat_clusters" %in% colnames(spatial_obj()@meta.data)) {
      obj <- spatial_obj()
      
      # Get cluster choices - use annotated names if available
      if ("annotated_clusters" %in% colnames(obj@meta.data)) {
        cluster_choices <- levels(obj$annotated_clusters)
      } else {
        cluster_choices <- levels(obj$seurat_clusters)
      }
      
      updateSelectInput(session, "spatial_marker_cluster", 
                        choices = cluster_choices,
                        selected = cluster_choices[1])
      
      updateSelectInput(session, "spatial_comparison_clusters", 
                        choices = cluster_choices,
                        selected = NULL)
    }
  })
  
  # Update available clusters for comparison based on selected target
  observeEvent(input$spatial_marker_cluster, {
    req(spatial_obj(), input$spatial_marker_cluster)
    
    obj <- spatial_obj()
    
    # Get all clusters
    if ("annotated_clusters" %in% colnames(obj@meta.data)) {
      all_clusters <- levels(obj$annotated_clusters)
    } else {
      all_clusters <- levels(obj$seurat_clusters)
    }
    
    # Remove the selected cluster from comparison options
    available_clusters <- setdiff(all_clusters, input$spatial_marker_cluster)
    
    updateSelectInput(session, "spatial_comparison_clusters", 
                      choices = available_clusters,
                      selected = available_clusters[1])
  })
  
  # Store marker results in reactive value
  spatial_marker_results_reactive <- reactiveVal(NULL)
  
  # Dynamic title for results
  output$marker_results_title <- renderText({
    if (input$spatial_comparison_type == "one_vs_all") {
      paste("Markers for", input$spatial_marker_cluster, "vs All Other Clusters")
    } else {
      comparison_clusters <- input$spatial_comparison_clusters
      if (length(comparison_clusters) > 0) {
        paste("Markers for", input$spatial_marker_cluster, "vs", 
              paste(comparison_clusters, collapse = ", "))
      } else {
        "Select clusters to compare"
      }
    }
  })
  
  # Find markers for selected cluster
  observeEvent(input$run_spatial_markers, {
    message("=== FIND SPATIAL MARKERS STARTED ===")
    
    req(spatial_obj())
    req(input$spatial_marker_cluster)
    
    # Check if comparison clusters are selected when needed
    if (input$spatial_comparison_type == "one_vs_selected" && 
        (is.null(input$spatial_comparison_clusters) || length(input$spatial_comparison_clusters) == 0)) {
      showNotification("Please select at least one cluster to compare against", type = "warning")
      return()
    }
    
    tryCatch({
      obj <- spatial_obj()
      
      # Check if clusters exist
      cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) {
        "annotated_clusters"
      } else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        "seurat_clusters"
      } else {
        showNotification("No clusters found! Please run clustering first.", type = "error")
        return()
      }
      
      selected_cluster <- input$spatial_marker_cluster
      comparison_type <- input$spatial_comparison_type
      
      message(paste("Finding markers for cluster:", selected_cluster))
      message(paste("Comparison type:", comparison_type))
      
      # Prepare modal message
      modal_message <- if (comparison_type == "one_vs_all") {
        paste("Analyzing cluster", selected_cluster, "vs all other clusters...")
      } else {
        paste("Analyzing cluster", selected_cluster, "vs clusters:", 
              paste(input$spatial_comparison_clusters, collapse = ", "))
      }
      
      showModal(modalDialog(
        title = "Finding Marker Genes",
        div(
          modal_message,
          br(),
          tags$small("This may take a few minutes")
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      
      # Set identity to the cluster column
      Idents(obj) <- cluster_col
      
      # Find markers based on comparison type
      if (comparison_type == "one_vs_all") {
        # One vs all comparison
        markers <- FindMarkers(
          obj,
          ident.1 = selected_cluster,
          logfc.threshold = input$spatial_marker_logfc,
          min.pct = input$spatial_marker_min_pct,
          only.pos = input$spatial_marker_only_pos,
          verbose = FALSE
        )
      } else {
        # One vs selected comparison
        markers <- FindMarkers(
          obj,
          ident.1 = selected_cluster,
          ident.2 = input$spatial_comparison_clusters,
          logfc.threshold = input$spatial_marker_logfc,
          min.pct = input$spatial_marker_min_pct,
          only.pos = input$spatial_marker_only_pos,
          verbose = FALSE
        )
      }
      
      # Add gene names as a column
      markers$gene <- rownames(markers)
      markers <- markers[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
      
      # Sort by fold change
      markers <- markers[order(-markers$avg_log2FC), ]
      
      # Round numeric columns
      markers$avg_log2FC <- round(markers$avg_log2FC, 3)
      markers$pct.1 <- round(markers$pct.1, 3)
      markers$pct.2 <- round(markers$pct.2, 3)
      markers$p_val <- formatC(markers$p_val, format = "e", digits = 2)
      markers$p_val_adj <- formatC(markers$p_val_adj, format = "e", digits = 2)
      
      # Store results
      spatial_marker_results_reactive(markers)
      
      message(paste("Found", nrow(markers), "marker genes"))
      
      removeModal()
      
      
      message("=== FIND SPATIAL MARKERS COMPLETED ===")
      
    }, error = function(e) {
      message(paste("=== FIND SPATIAL MARKERS ERROR ===", e$message))
      removeModal()
      showNotification(paste("Error finding markers:", e$message), type = "error")
    })
  })
  
  
  # Display marker results
  output$spatial_marker_results <- renderDT({
    markers <- spatial_marker_results_reactive()
    
    if (is.null(markers)) {
      return(datatable(data.frame(Message = "Click 'Find Markers' to analyze a cluster"),
                       options = list(dom = 't'),
                       rownames = FALSE))
    }
    
    datatable(
      markers,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel'),
        columnDefs = list(
          list(className = 'dt-center', targets = 1:5)
        )
      ),
      rownames = FALSE,
      caption = paste("Marker genes for cluster", input$spatial_marker_cluster)
    ) %>%
      formatStyle('avg_log2FC',
                  backgroundColor = styleInterval(c(0.5, 1), c('white', 'lightblue', 'darkblue')),
                  color = styleInterval(c(1), c('black', 'white')))
  })
  
  # Download marker results - CORRECTED VERSION
  output$download_spatial_markers <- downloadHandler(
    filename = function() {
      cluster_name <- gsub(" ", "_", input$spatial_marker_cluster)
      paste0("markers_cluster_", cluster_name, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      markers <- spatial_marker_results_reactive()
      
      # Simply write the file, no notifications here
      if (!is.null(markers)) {
        write.csv(markers, file, row.names = FALSE)
      } else {
        # Write empty file if no data
        write.csv(data.frame(Message = "No results available"), file, row.names = FALSE)
      }
    }
  )
}