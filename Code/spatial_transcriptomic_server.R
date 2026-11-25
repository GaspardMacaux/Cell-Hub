# ==============================================================================
# spatial_transcriptomic_server.R
# Server logic for the spatial transcriptomics module (10x Visium / Visium HD).
# All reactive code lives here. Pure helpers live in functions/spatial_functions.R
# (ZIP staging and loading, RDS loading, coordinates, QC suggestions, tissue viewer, co-localization).
# Other helpers used: dimensional_reduction_functions.R (seed-42 wrappers),
# download_functions.R (createDownloadHandler), cells_genes_expressions.R (create_expression_percentage_bar).
# ==============================================================================

spatial_transcriptomic_server <- function(input, output, session) {
  
  ############################## Reactive state ##############################
  
  # Main spatial object and processing flags
  spatial_obj <- reactiveVal(NULL)
  sketching_applied <- reactiveVal(FALSE)
  available_assays <- reactiveVal(c("Spatial"))
  # Assay used for scaling / PCA (set through "Apply assay selection")
  selected_assay <- reactiveVal("Spatial")
  processing_states <- reactiveValues(qc_applied = FALSE, normalized = FALSE, pca_computed = FALSE, clustered = FALSE)
  # Individual state flags (kept for compatibility)
  spatial_normalized <- reactiveVal(FALSE)
  spatial_qc_applied <- reactiveVal(FALSE)
  spatial_pca_computed <- reactiveVal(FALSE)
  spatial_clustered <- reactiveVal(FALSE)
  # Plot parameters kept for downloads
  spatial_plots <- reactiveValues()
  # Plots generated in the gene expression tab, reused by downloads and the annotation tab
  spatial_feature_plot_reactive <- reactiveVal(NULL)
  spatial_violin_plot_reactive <- reactiveVal(NULL)
  spatial_dot_plot_reactive <- reactiveVal(NULL)
  # Cluster annotation state and annotation plots (stored for downloads)
  spatial_cluster_names <- reactiveVal(NULL)
  # Analysis results
  coloc_results <- reactiveVal(NULL)
  spatial_marker_results_reactive <- reactiveVal(NULL)
  # Interactive tissue viewer: window in full-resolution pixels c(xmin, xmax, ymin, ymax) (NULL = whole image) and clicked spot
  viewer_view <- reactiveVal(NULL)
  viewer_selected_spot <- reactiveVal(NULL)
  # Decoded tissue image cache (plain environment: the hires PNG is decoded once per image, not at every object update)
  viewer_cache <- new.env()
  # Last gene set pushed to the viewer gene menu (plain variable: the menu refreshes only when the gene set changes)
  last_viewer_genes <- character(0)
  # Last assay set pushed to the menus (plain variable: menus refresh only when the set changes)
  last_assay_choices <- character(0)
  # Spatially variable genes (Moran's I) of the current object
  svg_results <- reactiveVal(NULL)
  # Spots selected on the tissue section map with the mouse (barcodes), removed with "Remove selected spots"
  section_selection <- reactiveVal(character(0))
  
  ############################## Download buttons state ##############################
  # Each download button stays disabled until the result it exports exists (cheap checks: the plots themselves are not drawn here)
  bind_download_state("download_viewer_plot", function() req(spatial_obj(), length(spatial_obj()@images) > 0))
  bind_download_state("download_spatial_features", function() spatial_feature_plot_reactive())
  bind_download_state("download_spatial_violin", function() spatial_violin_plot_reactive())
  bind_download_state("download_spatial_dotplot", function() spatial_dot_plot_reactive())
  bind_download_state("save_spatial_object", function() req(spatial_obj()))
  bind_download_state("download_spatial_umap", function() req(spatial_obj(), "umap" %in% names(spatial_obj()@reductions)))
  bind_download_state("download_spatial_tissue", function() req(spatial_obj(), length(spatial_obj()@images) > 0))
  bind_download_state("save_annotated_spatial_object", function() req(spatial_obj()))
  bind_download_state("download_coloc_map", function() { r <- coloc_results(); req(identical(r$analysis_type, "bidirectional")) })
  bind_download_state("download_coloc_bar", function() req(coloc_results()))
  bind_download_state("download_coloc_results", function() req(coloc_results()))
  bind_download_state("download_svg_table", function() req(svg_results()))
  bind_download_state("download_svg_plot", function() req(spatial_obj(), input$svg_plot_gene))
  bind_download_state("download_spatial_markers", function() req(spatial_marker_results_reactive()))
  
  ############################## Internal helpers ##############################
  
  # Update the object and the processing flags in one call
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
  
  # Human-readable processing status (QC → Normalized → Sketched → PCA → Clustered)
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
  
  # Checks run on every newly loaded object, before it is stored:
  # 1. spot coordinates transposed relative to the tissue image (x / y convention lost between Seurat versions or after saving
  #    and reloading): the stored centroids are swapped (transpose_spot_coordinates), so Seurat plots and Cell-Hub views agree;
  #    if that fails, @misc$spot_coords_transposed corrects the Cell-Hub views only (read by get_spot_coordinates_spatial());
  # 2. infer_size = TRUE (.rds files): spot / bin size, standard Visium from the spot grid, Visium HD from a stored value or the
  #    barcode / image / assay names. An undetermined Visium HD bin size leaves @misc$spot_size_um empty (see show_bin_size_modal).
  check_loaded_spatial_object <- function(obj, infer_size = TRUE) {
    if (length(obj@images) == 0) return(obj)
    orientation <- tryCatch(detect_spot_transposition(obj), error = function(e) { message("Orientation check skipped: ", conditionMessage(e)); NULL })
    if (!is.null(orientation)) {
      message(sprintf("Spots on tissue: %.0f%% as read, %.0f%% with x / y swapped", 100 * orientation$on_tissue_as_read, 100 * orientation$on_tissue_swapped))
      obj@misc$spot_coords_transposed <- NULL
      if (isTRUE(orientation$transposed)) {
        before <- get_spot_coordinates_spatial(obj)
        fixed <- tryCatch(transpose_spot_coordinates(obj), error = function(e) { message("Stored coordinates not corrected: ", conditionMessage(e)); NULL })
        after <- if (!is.null(fixed)) tryCatch(get_spot_coordinates_spatial(fixed), error = function(e) NULL) else NULL
        swapped <- !is.null(after) && identical(after$barcode, before$barcode) && isTRUE(all.equal(after$imagerow, before$imagecol)) && isTRUE(all.equal(after$imagecol, before$imagerow))
        if (swapped) {
          obj <- fixed
          message("Stored spot coordinates transposed (Seurat plots and Cell-Hub views)")
        } else {
          obj@misc$spot_coords_transposed <- TRUE
          message("Stored coordinates unchanged: Cell-Hub views corrected, Seurat plots may stay transposed")
        }
        showNotification(sprintf("Spot coordinates were transposed relative to the tissue image (%.0f%% of spots on tissue, %.0f%% after correction): corrected%s.",
                                 100 * orientation$on_tissue_as_read, 100 * orientation$on_tissue_swapped, if (swapped) "" else " in the Cell-Hub views only"),
                         type = "warning", duration = 12)
      }
    }
    if (infer_size) {
      size <- tryCatch(infer_spot_size_um(obj), error = function(e) list(size_um = NA_real_, technology = "unknown", source = conditionMessage(e)))
      message(sprintf("Spot size: %s (%s, %s)", size$size_um, size$technology, size$source))
      obj@misc$spot_size_um <- if (is.finite(size$size_um)) size$size_um else NULL
      if (is.finite(size$size_um) && size$technology == "Visium HD") showNotification(sprintf("Visium HD object: bin size %g µm (%s)", size$size_um, size$source), type = "message", duration = 8)
    }
    return(obj)
  }
  
  # Co-localization scale for a spot / bin size: spot size field and default distance, the nearest neighbours
  # (100 um centre-to-centre for standard Visium spots, 2 bins for Visium HD; below the grid pitch only the same spot counts)
  update_coloc_scale <- function(spot_size_um) {
    if (is.null(spot_size_um) || !is.finite(spot_size_um)) return(invisible(NULL))
    updateNumericInput(session, "coloc_spot_size_um", value = spot_size_um)
    updateNumericInput(session, "coloc_distance_um", value = if (spot_size_um >= 50) 100 else 2 * spot_size_um)
  }
  
  # Ask for the bin size of a Visium HD object whose size could not be determined (answer handled by input$confirm_rds_bin_size)
  show_bin_size_modal <- function(obj) {
    showModal(modalDialog(
      title = "Visium HD bin size",
      tags$p(sprintf("This object has %s contiguous bins: it is a Visium HD object, but its bin size is not stored in it.", format(ncol(obj), big.mark = ","))),
      tags$p("The bin size sets the micrometre scale used by co-localization distances and spot sizes."),
      radioButtons("rds_bin_size_um", "Bin size:", choices = c("2 µm" = 2, "8 µm" = 8, "16 µm" = 16), selected = 8, inline = TRUE),
      footer = actionButton("confirm_rds_bin_size", "Apply", class = "btn-primary"), easyClose = FALSE
    ))
  }
  
  ############################## Assay management ##############################
  
  # Refresh assay menus only when the set of assays changes, so user selections survive QC / PCA / clustering updates
  observe(run_safely("Interface update", persistent = FALSE, {
    req(spatial_obj())
    obj <- spatial_obj()
    assays <- names(obj@assays)
    available_assays(assays)
    if (identical(assays, last_assay_choices)) return()
    last_assay_choices <<- assays
    # Analysis assay of the PCA card: the default assay (SCT after SCTransform, Spatial otherwise)
    updateSelectInput(session, "spatial_assay_select", choices = assays, selected = DefaultAssay(obj))
    updateSelectInput(session, "spatial_viz_assay", choices = assays, selected = DefaultAssay(obj))
    # The co-localization assay menu had static choices and broke after SCT normalization
    updateSelectInput(session, "coloc_assay", choices = assays, selected = DefaultAssay(obj))
  }))
  
  # Current assay information
  output$spatial_assay_info <- renderText({
    req(spatial_obj())
    obj <- spatial_obj()
    assay_info <- paste("Available assays:", paste(names(obj@assays), collapse = ", "))
    if (sketching_applied()) {
      # Actual number of cells in the sketch assay
      sketch_ncells <- ncol(obj[["sketch"]])
      sketch_info <- paste("Sketched data available with", format(sketch_ncells, big.mark = ","), "cells")
      assay_info <- paste(assay_info, "\n", sketch_info)
    }
    selected_assay_name <- input$spatial_viz_assay
    if (!is.null(selected_assay_name) && selected_assay_name %in% names(obj@assays)) {
      # Number of features and cells of the selected assay
      n_features <- nrow(obj[[selected_assay_name]])
      n_cells <- ncol(obj[[selected_assay_name]])
      current_info <- paste("Current assay:", selected_assay_name, "with",
                            format(n_features, big.mark = ","), "features and",
                            format(n_cells, big.mark = ","), "cells")
      assay_info <- paste(assay_info, "\n", current_info)
    }
    return(assay_info)
  })
  
  # Set the default assay of the Seurat object
  observeEvent(input$set_default_assay, run_safely("Default assay", {
    req(spatial_obj())
    tryCatch({
      obj <- spatial_obj()
      selected_assay <- input$spatial_default_assay
      if (selected_assay %in% names(obj@assays)) {
        # Set as default assay in the Seurat object
        DefaultAssay(obj) <- selected_assay
        # Update reactive values
        spatial_obj(obj)
        selected_assay(selected_assay)
        showNotification(paste("Default assay set to:", selected_assay), type = "message", duration = 5)
        # Update the assay selection menus to reflect the change
        updateSelectInput(session, "spatial_viz_assay", selected = selected_assay)
        updateSelectInput(session, "spatial_assay_select", selected = selected_assay)
      } else {
        showNotification("Selected assay not found", type = "error")
      }
    }, error = function(e) {
      showNotification(paste("Error setting assay:", conditionMessage(e)[1]), type = "error")
    })
  }))
  
  # Keep the default assay menu in sync with the object
  observe(run_safely("Interface update", persistent = FALSE, {
    req(spatial_obj())
    obj <- spatial_obj()
    updateSelectInput(session, "spatial_default_assay", choices = names(obj@assays), selected = DefaultAssay(obj))
  }))
  
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
  
  ############################## Data loading ##############################
  
  observeEvent(input$spatial_file, run_safely("spatial file", {
    req(input$spatial_file)
    message("=== SPATIAL FILE UPLOAD STARTED ===")
    showModal(modalDialog(title = "Loading Spatial Data", "Processing your spatial transcriptomics data...", easyClose = FALSE, footer = NULL))
    # Defined outside tryCatch so the finally clause can remove the staged files (tempfile: unique even for two uploads in the same second)
    temp_dir <- tempfile(pattern = "spatial_")
    tryCatch({
      staged <- prepare_visium_directory(input$spatial_file$datapath, temp_dir, bin_size_um = as.numeric(input$spatial_hd_bin_size %||% 8))
      spatial_data <- load_visium_object(staged, species = input$spatial_species_choice %||% "mouse")
      # Physical spot diameter used for micrometre conversions (55 um Visium spot, bin size for Visium HD)
      spatial_data@misc$spot_size_um <- staged$spot_size_um
      # Orientation of the spots on the tissue image (the bin size is known from the ZIP)
      spatial_data <- check_loaded_spatial_object(spatial_data, infer_size = FALSE)
      # A new dataset invalidates every state derived from the previous one (Visium barcodes repeat across slides,
      # so a stale spot selection would point to spots of the new slide)
      sketching_applied(FALSE); selected_assay("Spatial"); spatial_cluster_names(NULL); coloc_results(NULL)
      section_selection(character(0)); spatial_marker_results_reactive(NULL); viewer_selected_spot(NULL)
      update_spatial_object(spatial_data, list(qc_applied = FALSE, normalized = FALSE, pca_computed = FALSE, clustered = FALSE))
      update_coloc_scale(staged$spot_size_um)
      showNotification(paste0("Spatial data loaded successfully! Spots: ", format(ncol(spatial_data), big.mark = ","), ", Genes: ", format(nrow(spatial_data), big.mark = ",")), type = "message", duration = 8)
      message("=== SPATIAL LOADING COMPLETED ===")
    }, error = function(e) {
      message(paste("=== SPATIAL LOADING ERROR ===", conditionMessage(e)))
      showNotification(paste("Error loading spatial data:", conditionMessage(e)), type = "error", duration = 15)
    }, finally = {
      removeModal()
      # Matrix and image are read eagerly, so the staged files are no longer needed
      unlink(temp_dir, recursive = TRUE, force = TRUE)
    })
  }))
  
  # Load a processed Seurat object from an RDS file.
  # loadSpatialSeuratObject() ALWAYS runs UpdateSeuratObject() (spatial images need it even with Assay5).
  observeEvent(input$load_spatial_seurat_file, run_safely("load spatial seurat file", {
    req(input$load_spatial_seurat_file)
    message("=== LOADING SPATIAL SEURAT OBJECT FROM RDS ===")
    tryCatch({
      showModal(modalDialog(
        title = "Loading Spatial Object",
        "Reading your saved spatial Seurat object...",
        easyClose = FALSE, footer = NULL
      ))
      file_path <- input$load_spatial_seurat_file$datapath
      message(paste("Loading file from:", file_path))
      # Spatial-specific loader: always UpdateSeuratObject, v4->v5, JoinLayers
      loaded_obj <- loadSpatialSeuratObject(
        rds_path = file_path,
        species = input$spatial_species_choice %||% "mouse",
        log_function = message
      )
      # Orientation of the spots on the tissue image and spot / bin size
      loaded_obj <- check_loaded_spatial_object(loaded_obj)
      # A new dataset invalidates every state derived from the previous one: results, spot selection (Visium barcodes
      # repeat across slides), marker table, cluster names (NULL lets the initializer restore the object's own annotation)
      coloc_results(NULL); section_selection(character(0)); spatial_marker_results_reactive(NULL)
      spatial_cluster_names(NULL); viewer_selected_spot(NULL); selected_assay(DefaultAssay(loaded_obj))
      # Detect processing state from object contents
      has_spatial <- "Spatial" %in% names(loaded_obj@assays)
      has_images <- length(loaded_obj@images) > 0
      has_pca <- "pca" %in% names(loaded_obj@reductions)
      has_umap <- "umap" %in% names(loaded_obj@reductions)
      has_clusters <- "seurat_clusters" %in% colnames(loaded_obj@meta.data)
      # SCTransform leaves the Spatial assay with counts only: an SCT assay means the object is normalized
      has_sct <- "SCT" %in% names(loaded_obj@assays)
      # Safe normalization detection: only inspect the "data" layer if it exists
      # (raw objects only carry "counts" and requesting "data" would crash)
      has_normalized <- FALSE
      if (has_spatial) {
        available_layers <- tryCatch(names(loaded_obj[["Spatial"]]@layers),
                                     error = function(e) character(0))
        if ("data" %in% available_layers && !"counts" %in% available_layers) {
          has_normalized <- TRUE
        } else if ("data" %in% available_layers && "counts" %in% available_layers) {
          has_normalized <- tryCatch({
            assay_data <- LayerData(loaded_obj, assay = "Spatial", layer = "data")
            n_vals <- length(assay_data@x)
            if (n_vals == 0) FALSE else max(assay_data@x[seq_len(min(1000, n_vals))]) < 50
          }, error = function(e) FALSE)
        }
      }
      has_normalized <- isTRUE(has_normalized) || has_sct
      # Update reactive state (isTRUE guards against NA)
      update_spatial_object(loaded_obj, list(
        qc_applied   = has_normalized,
        normalized   = has_normalized,
        pca_computed = isTRUE(has_pca),
        clustered    = isTRUE(has_clusters)
      ))
      sketching_applied("sketch" %in% names(loaded_obj@assays))
      # Update UI selectors
      assay_names <- names(loaded_obj@assays)
      available_assays(assay_names)
      updateSelectInput(session, "spatial_assay_select",
                        choices = assay_names, selected = DefaultAssay(loaded_obj))
      updateSelectInput(session, "spatial_viz_assay",
                        choices = assay_names, selected = DefaultAssay(loaded_obj))
      # Build success message
      info_parts <- c()
      if (has_spatial) info_parts <- c(info_parts, "Spatial data")
      if (has_images) info_parts <- c(info_parts, "Images")
      if (has_normalized) info_parts <- c(info_parts, if (has_sct) "Normalized (SCT)" else "Normalized")
      if (has_pca) info_parts <- c(info_parts, "PCA")
      if (has_umap) info_parts <- c(info_parts, "UMAP")
      if (has_clusters) info_parts <- c(info_parts, "Clusters")
      showNotification(
        paste0("Loaded: ", format(ncol(loaded_obj), big.mark = ","), " spots, ",
               format(nrow(loaded_obj), big.mark = ","), " genes",
               if (length(info_parts) > 0) paste(" |", paste(info_parts, collapse = ", ")) else ""),
        type = "message", duration = 8)
      removeModal()
      # Bin size of the co-localization tab, or a question when a Visium HD bin size is unknown (after the loading modal is closed)
      if (is.null(loaded_obj@misc$spot_size_um)) show_bin_size_modal(loaded_obj) else update_coloc_scale(loaded_obj@misc$spot_size_um)
      message("=== SPATIAL OBJECT LOADING COMPLETED ===")
    }, error = function(e) {
      message(paste("=== SPATIAL OBJECT LOADING ERROR ===", conditionMessage(e)[1]))
      removeModal()
      showNotification(paste("Error loading spatial object:", conditionMessage(e)[1]),
                       type = "error", duration = 10)
    })
  }))
  
  # Bin size chosen in the modal opened after loading a Visium HD .rds without a known bin size
  observeEvent(input$confirm_rds_bin_size, run_safely("Set bin size", {
    obj <- spatial_obj(); req(obj)
    bin <- as.numeric(input$rds_bin_size_um)
    obj@misc$spot_size_um <- bin
    spatial_obj(obj)
    update_coloc_scale(bin)
    removeModal()
    showNotification(sprintf("Bin size set to %g µm", bin), type = "message")
  }))
  
  # A new dataset invalidates the spatially variable genes of the previous one
  observeEvent(list(input$spatial_file, input$load_spatial_seurat_file), {
    run_safely("Spatially variable genes reset", persistent = FALSE, {
      svg_results(NULL)
    })
  }, ignoreInit = TRUE)
  
  # Replace the tissue image with the full-resolution image used by Space Ranger (downsampled, aligned via scale factors)
  observeEvent(input$upload_high_res_image, run_safely("High-resolution image upload", {
    req(input$upload_high_res_image, spatial_obj())
    tryCatch({
      message("=== LOADING HIGH-RES IMAGE ===")
      current_obj <- spatial_obj()
      if (length(current_obj@images) == 0) {
        stop("No spatial images found")
      }
      image_name <- names(current_obj@images)[1]
      message("Current image class: ", class(current_obj@images[[image_name]]))
      message("Current image dimensions: ", paste(dim(current_obj@images[[image_name]]@image), collapse = " x "))
      # Downsample the full-resolution image and align it through the lowres/hires scale factors
      current_obj <- replace_spatial_image(current_obj, input$upload_high_res_image$datapath, max_side = 2000, log_function = message)
      new_dims <- dim(current_obj@images[[image_name]]@image)[1:2]
      message("Image loaded: ", paste(dim(current_obj@images[[image_name]]@image), collapse = " x "))
      # Check the tissue coordinate range (full-resolution pixels)
      tissue_coords <- get_spot_coordinates_spatial(current_obj, image_name)
      message("Tissue coordinates range: x[", min(tissue_coords$imagecol), ",", max(tissue_coords$imagecol), "] y[", min(tissue_coords$imagerow), ",", max(tissue_coords$imagerow), "]")
      spatial_obj(current_obj)
      message("=== IMAGE REPLACED SUCCESSFULLY ===")
      showNotification(paste("High-res image loaded:", paste(new_dims, collapse = "x")), type = "message")
    }, error = function(e) {
      message("ERROR: ", e$message)
      message("Traceback: ")
      print(traceback())
      showNotification(paste("Error:", conditionMessage(e)[1]), type = "error")
    })
  }))
  
  # Current tissue image information
  output$current_image_info <- renderText({
    req(spatial_obj())
    if (length(spatial_obj()@images) == 0) {
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
  
  # Keep the lowres / hires selector in sync with the image actually displayed (ZIP load, RDS load, custom image upload)
  observe(run_safely("Interface update", persistent = FALSE, {
    req(spatial_obj())
    active <- spatial_obj()@misc$active_tissue_image
    current <- isolate(input$spatial_image_resolution)
    if (is.null(active) || !active %in% c("lowres", "hires")) {
      # Custom image or object without stored images: no option is active, so any click triggers a switch
      if (!is.null(current)) updateRadioButtons(session, "spatial_image_resolution", selected = character(0))
    } else if (!identical(active, current)) {
      updateRadioButtons(session, "spatial_image_resolution", selected = active)
    }
  }))
  
  # Switch the displayed tissue image between the stored lowres and hires versions (spots stay aligned through the scale factors)
  observeEvent(input$spatial_image_resolution, run_safely("spatial image resolution", {
    req(spatial_obj())
    obj <- spatial_obj()
    if (identical(input$spatial_image_resolution, obj@misc$active_tissue_image)) return()
    showModal(modalDialog(title = "Switching tissue image", "Decoding the selected image...", easyClose = FALSE, footer = NULL))
    tryCatch({
      obj <- set_spatial_image_resolution(obj, input$spatial_image_resolution, log_function = message)
      spatial_obj(obj)
      showNotification(paste("Tissue image switched to", input$spatial_image_resolution), type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", conditionMessage(e)[1]), type = "error", duration = 10)
      # Put the selector back on the image still displayed
      active <- obj@misc$active_tissue_image
      updateRadioButtons(session, "spatial_image_resolution", selected = if (!is.null(active) && active %in% c("lowres", "hires")) active else character(0))
    }, finally = removeModal())
  }), ignoreInit = TRUE)
  
  ############################## Dataset information display ##############################
  
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
  
  # Technology inferred from the physical spot size stored at load time (55 um Visium spot, bin size for Visium HD)
  output$spatial_technology_info <- renderInfoBox({
    spot_um <- if (!is.null(spatial_obj())) spatial_obj()@misc$spot_size_um else NULL
    tech <- if (is.null(spot_um)) "Unknown" else if (spot_um >= 50) "Visium (55 µm spots)" else paste0("Visium HD (", spot_um, " µm bins)")
    infoBox(title = "Technology", value = tech, icon = icon("microscope"), color = "purple", fill = TRUE)
  })
  
  # Dataset summary table (metrics missing from the metadata are reported as NA instead of failing)
  output$spatial_dataset_summary <- renderDT({
    req(spatial_obj())
    obj <- spatial_obj()
    md <- obj@meta.data
    # Median or mean of a metadata column, NA when the column is absent
    meta_stat <- function(col, fun) if (col %in% colnames(md)) fun(md[[col]], na.rm = TRUE) else NA_real_
    med_genes <- meta_stat("nFeature_Spatial", median)
    med_umis <- meta_stat("nCount_Spatial", median)
    mean_mt <- meta_stat("percent.mt", mean)
    spot_um <- obj@misc$spot_size_um
    summary_data <- data.frame(
      Metric = c("Total Spots/Cells", "Total Genes", "Spot / bin size", "Available Assays", "Default Assay", "Sketched Data",
                 "Median Genes per Spot", "Median UMI per Spot", "Mean Mitochondrial %", "Has Tissue Image", "Processing Status"),
      Value = c(format(ncol(obj), big.mark = ","), format(nrow(obj), big.mark = ","),
                if (is.null(spot_um)) "Unknown" else paste(spot_um, "µm"),
                paste(names(obj@assays), collapse = ", "), DefaultAssay(obj),
                if (sketching_applied()) "Yes" else "No",
                if (is.na(med_genes)) "NA" else format(round(med_genes), big.mark = ","),
                if (is.na(med_umis)) "NA" else format(round(med_umis), big.mark = ","),
                if (is.na(mean_mt)) "NA" else paste0(round(mean_mt, 2), "%"),
                if (length(obj@images) > 0) "Yes" else "No",
                get_processing_status()),
      stringsAsFactors = FALSE
    )
    datatable(summary_data, options = list(pageLength = 15, searching = FALSE, lengthChange = FALSE, info = FALSE, paging = FALSE), rownames = FALSE)
  })
  
  ############################## Tissue sections ##############################
  # Helpers (detect_tissue_sections, add_tissue_sections_to_object, assign_tissue_sections, summarize_tissue_sections)
  # live in spatial_functions.R. Sections, samples and conditions are stored in the spot metadata
  # (section_detected, section, sample, condition) and in @misc$tissue_sections, so they survive saving and reloading.
  
  observeEvent(input$detect_sections, {
    run_safely("Tissue section detection", {
      req(spatial_obj())
      showModal(modalDialog(title = "Tissue sections", "Detecting separate tissue sections...", easyClose = FALSE, footer = NULL))
      # finally closes the modal; the error itself is reported by run_safely
      tryCatch({
        sections <- detect_tissue_sections(spatial_obj(), max_gap_pitch = input$section_max_gap %||% 2, min_spots = input$section_min_spots %||% 50)
        spatial_obj(add_tissue_sections_to_object(spatial_obj(), sections))
        n_sections <- sum(sections$summary$section != "Debris")
        n_debris <- sum(sections$summary$n_spots[sections$summary$section == "Debris"])
        showNotification(sprintf("%d section(s) detected%s. Name them and assign sample and condition below.", n_sections,
                                 if (n_debris > 0) sprintf(", %d spots in small fragments (Debris)", n_debris) else ""), type = "message", duration = 8)
      }, finally = removeModal())
    })
  })
  
  # One row per detected section: remove box, name, sample (animal) and condition, pre-filled with the last assignment
  output$section_assignment_ui <- renderUI({
    req(spatial_obj())
    info <- spatial_obj()@misc$tissue_sections
    if (is.null(info) || !"section_detected" %in% colnames(spatial_obj()@meta.data)) return(tags$p("Click 'Detect sections' first.", style = "color: #6c757d;"))
    previous <- info$assignment
    value_of <- function(id, column, default) {
      if (is.null(previous) || !id %in% previous$section) return(default)
      v <- previous[[column]][previous$section == id][1]
      if (is.na(v)) "" else as.character(v)
    }
    # Sections still present in the object (removed ones disappear from the table)
    present <- levels(droplevels(spatial_obj()$section_detected))
    summary <- info$summary[info$summary$section %in% present, , drop = FALSE]
    header <- fluidRow(column(2, tags$b("Detected")), column(1, tags$b("Remove")), column(3, tags$b("Name")), column(3, tags$b("Sample (animal)")), column(3, tags$b("Condition")))
    rows <- lapply(seq_len(nrow(summary)), function(i) {
      id <- summary$section[i]
      fluidRow(
        column(2, tags$p(tags$b(id), sprintf(" (%s spots)", format(summary$n_spots[i], big.mark = ",")), style = "margin-top: 8px;")),
        column(1, checkboxInput(paste0("section_remove_", id), NULL, value = FALSE)),
        column(3, textInput(paste0("section_name_", id), NULL, value = value_of(id, "name", id))),
        column(3, textInput(paste0("section_sample_", id), NULL, value = value_of(id, "sample", if (id == "Debris") "" else id))),
        column(3, textInput(paste0("section_condition_", id), NULL, value = value_of(id, "condition", "")))
      )
    })
    tagList(header, rows,
            tags$p(icon("info-circle"), " Same name for two sections = merged (e.g. a section split by a tear). Same sample for several sections = sections of one animal (not biological replicates); sample and condition are stored with the spots for sample-level analyses outside Cell-Hub. Tick 'Remove' and click 'Remove ticked sections' to delete sections.",                   style = "font-size: 12px; color: #6c757d;"))
  })
  
  # Remove every spot of the ticked sections (every downstream result is recomputed on the remaining spots)
  observeEvent(input$remove_checked_sections, {
    run_safely("Section removal", {
      req(spatial_obj())
      obj <- spatial_obj()
      if (!"section_detected" %in% colnames(obj@meta.data)) {
        showNotification("Detect the tissue sections first", type = "warning")
        return()
      }
      ids <- levels(droplevels(obj$section_detected))
      ticked <- ids[vapply(ids, function(id) isTRUE(input[[paste0("section_remove_", id)]]), logical(1))]
      if (length(ticked) == 0) {
        showNotification("Tick 'Remove' for at least one section", type = "warning")
        return()
      }
      keep <- colnames(obj)[is.na(obj$section_detected) | !as.character(obj$section_detected) %in% ticked]
      if (length(keep) == 0) {
        showNotification("Removing these sections would leave no spot", type = "error")
        return()
      }
      n_removed <- ncol(obj) - length(keep)
      obj <- subset(obj, cells = keep)
      for (col in intersect(c("section", "section_detected"), colnames(obj@meta.data))) obj@meta.data[[col]] <- droplevels(obj@meta.data[[col]])
      spatial_obj(obj)
      section_selection(character(0)); svg_results(NULL); coloc_results(NULL)
      showNotification(sprintf("Section(s) %s removed (%s spots)", paste(ticked, collapse = ", "), format(n_removed, big.mark = ",")), type = "message", duration = 6)
    })
  })
  
  # Remove the spots labelled Debris from the object (every downstream result is recomputed on the remaining spots)
  observeEvent(input$remove_debris_spots, {
    run_safely("Debris removal", {
      req(spatial_obj())
      obj <- spatial_obj()
      if (!"section" %in% colnames(obj@meta.data) || !any(obj$section == "Debris", na.rm = TRUE)) {
        showNotification("No spot is labelled Debris", type = "warning")
        return()
      }
      keep <- colnames(obj)[is.na(obj$section) | obj$section != "Debris"]
      n_removed <- ncol(obj) - length(keep)
      obj <- subset(obj, cells = keep)
      obj$section <- droplevels(obj$section); obj$section_detected <- droplevels(obj$section_detected)
      spatial_obj(obj)
      svg_results(NULL); coloc_results(NULL)
      showNotification(sprintf("%s debris spots removed", format(n_removed, big.mark = ",")), type = "message", duration = 6)
    })
  })
  
  output$section_summary_table <- renderDT({
    req(spatial_obj())
    tab <- summarize_tissue_sections(spatial_obj())
    if (is.null(tab)) return(datatable(data.frame(Message = "Click 'Detect sections'"), options = list(dom = "t"), rownames = FALSE))
    datatable(tab, rownames = FALSE, selection = "none", options = list(pageLength = 20, dom = "t"))
  })
  
  # Section map drawn like the tissue viewer (tissue image + spots in full-resolution pixels), so a mouse rectangle selects spots.
  # Plot y is the negated image row (same convention as build_tissue_viewer_plot). Section labels are drawn here, white boxes
  # at the median position of each section, so they stay readable over the tissue.
  section_plot_reactive <- reactive({
    req(spatial_obj())
    obj <- spatial_obj()
    group <- input$section_plot_by %||% "section"
    shiny::validate(shiny::need(group %in% colnames(obj@meta.data), "Click 'Detect sections' first"),
                    shiny::need(length(obj@images) > 0, "No tissue image in this object"))
    coords <- get_spot_coordinates_spatial(obj)
    values <- as.character(obj@meta.data[coords$barcode, group])
    values[is.na(values) | values == ""] <- "Not assigned"
    spots <- data.frame(barcode = coords$barcode, x = coords$imagecol, y = coords$imagerow, value = factor(values), stringsAsFactors = FALSE)
    tissue <- viewer_tissue()
    view <- c(0, tissue$width_full, 0, tissue$height_full)
    width_px <- session$clientData$output_section_plot_width %||% 1000
    colors <- stats::setNames(scales::hue_pal()(nlevels(spots$value)), levels(spots$value))
    colors[names(colors) %in% c("Debris", "Not assigned")] <- "#9E9E9E"
    p <- build_tissue_viewer_plot(spots = spots, view = view, crop = crop_tissue_image(tissue, view, display_size = 1600), value_label = group,
                                  cluster_colors = colors, spot_diameter_px = obj@images[[1]]@scale.factors$spot %||% 10,
                                  spot_shape = if ((obj@misc$spot_size_um %||% 55) < 50) 15 else 16, size_factor = input$section_pt_size %||% 1.6,
                                  spot_alpha = 0.8, image_alpha = 1, show_labels = FALSE, show_legend = TRUE, panel_px = c(width_px, 650))
    # One readable label per group (section, sample or condition), Debris excluded
    labelled <- spots[!spots$value %in% c("Debris", "Not assigned"), , drop = FALSE]
    if (nrow(labelled) > 0) {
      centers <- stats::aggregate(cbind(x, y) ~ value, data = labelled, FUN = stats::median)
      p <- p + geom_label(data = centers, aes(x = x, y = -y, label = value), inherit.aes = FALSE, fill = "white", colour = "black",
                          fontface = "bold", size = 5, alpha = 0.9, label.padding = grid::unit(0.25, "lines"))
    }
    selected <- spots[spots$barcode %in% section_selection(), , drop = FALSE]
    if (nrow(selected) > 0) {
      p <- p + geom_point(data = selected, aes(x = x, y = -y), inherit.aes = FALSE, colour = "black", size = 0.6, shape = 16) +
        labs(subtitle = sprintf("%s spots selected (black)", format(nrow(selected), big.mark = ",")))
    }
    p
  })
  
  output$section_plot <- renderPlot({ section_plot_reactive() }, height = 650)
  
  # Each mouse rectangle adds its spots to the selection
  observeEvent(input$section_brush, {
    run_safely("Section map selection", persistent = FALSE, {
      b <- input$section_brush
      coords <- get_spot_coordinates_spatial(spatial_obj())
      inside <- coords$imagecol >= b$xmin & coords$imagecol <= b$xmax & -coords$imagerow >= b$ymin & -coords$imagerow <= b$ymax
      section_selection(union(section_selection(), coords$barcode[inside]))
      session$resetBrush("section_brush")
    })
  })
  
  observeEvent(input$clear_section_selection, {
    run_safely("Section map selection", persistent = FALSE, {
      section_selection(character(0))
    })
  })
  
  # Remove the selected spots from the object (every downstream result is recomputed on the remaining spots)
  observeEvent(input$remove_selected_spots, {
    run_safely("Spot removal", {
      req(spatial_obj())
      selected <- intersect(section_selection(), colnames(spatial_obj()))
      if (length(selected) == 0) {
        showNotification("Draw a rectangle on the map to select spots first", type = "warning")
        return()
      }
      obj <- spatial_obj()
      obj <- subset(obj, cells = setdiff(colnames(obj), selected))
      for (col in intersect(c("section", "section_detected"), colnames(obj@meta.data))) obj@meta.data[[col]] <- droplevels(obj@meta.data[[col]])
      spatial_obj(obj)
      section_selection(character(0))
      svg_results(NULL); coloc_results(NULL)
      showNotification(sprintf("%s spots removed (%s remaining)", format(length(selected), big.mark = ","), format(ncol(obj), big.mark = ",")), type = "message", duration = 6)
    })
  })
  
  # Apply the names, samples and conditions typed in the assignment table (same name = merged sections)
  observeEvent(input$apply_sections, run_safely("Tissue section assignment", {
    req(spatial_obj())
    info <- spatial_obj()@misc$tissue_sections
    if (is.null(info) || !"section_detected" %in% colnames(spatial_obj()@meta.data)) {
      showNotification("Detect the tissue sections first", type = "warning")
      return()
    }
    ids <- info$summary$section
    read_inputs <- function(prefix) vapply(ids, function(id) as.character(input[[paste0(prefix, id)]] %||% ""), character(1))
    assignment <- data.frame(section = ids, name = read_inputs("section_name_"), sample = read_inputs("section_sample_"),
                             condition = read_inputs("section_condition_"), stringsAsFactors = FALSE)
    spatial_obj(assign_tissue_sections(spatial_obj(), assignment))
    tab <- summarize_tissue_sections(spatial_obj())
    kept <- tab[tab$section != "Debris", , drop = FALSE]
    showNotification(sprintf("%d section(s), %d sample(s), %d condition(s) assigned", nrow(kept), length(unique(kept$sample)),
                             length(unique(kept$condition[kept$condition != ""]))), type = "message", duration = 6)
  }))
  
  ############################## Quality control ##############################
  
  # Generate QC plots (QC metrics are computed on the counts layer)
  observeEvent(input$spatial_qc_plots, run_safely("QC plots", {
    req(spatial_obj())
    tryCatch({
      obj <- spatial_obj()
      options(warn = -1)  # Suppress warnings for plotting
      # Check available layers in the Spatial assay
      spatial_assay <- obj[["Spatial"]]
      available_layers <- names(spatial_assay@layers)
      message(paste("Available layers in Spatial assay:", paste(available_layers, collapse = ", ")))
      # Use the counts layer explicitly for QC metrics
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
            # Fallback plot if the spatial plot fails
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
      })
      # Standard QC violin plots
      output$spatial_vlnplot <- renderPlot({
        VlnPlot(obj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),
                ncol = 3, pt.size = 0) &
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      })
      output$spatial_qc_violin <- renderPlot({
        # Side by side: three stacked violins are unreadable in a 350 px box
        VlnPlot(obj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),
                ncol = 3, pt.size = 0) &
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
      })
      # Feature scatter plots
      output$spatial_feature_scatter <- renderPlot({
        plot1 <- FeatureScatter(obj, feature1 = "nCount_Spatial", feature2 = "percent.mt",
                                pt.size = 0.5) + NoLegend()
        plot2 <- FeatureScatter(obj, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial",
                                pt.size = 0.5) + NoLegend()
        plot1 + plot2
      })
      output$spatial_scatter1 <- renderPlot({
        FeatureScatter(obj, feature1 = "nCount_Spatial", feature2 = "percent.mt",
                       pt.size = 0.5)
      })
      output$spatial_scatter2 <- renderPlot({
        FeatureScatter(obj, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial",
                       pt.size = 0.5)
      })
      options(warn = 0)
      showNotification("QC plots generated successfully!", type = "message")
    }, error = function(e) {
      options(warn = 0)
      showNotification(paste("Error generating QC plots:", conditionMessage(e)[1]), type = "error")
    })
  }))
  
  # Apply exactly the user's thresholds and record them; data-driven bounds are only suggested on failure
  observeEvent(input$apply_spatial_qc, run_safely("QC filtering", {
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
      # Apply exactly the user's thresholds (never silently replaced by data-driven ones)
      filters <- list(
        min_features = input$spatial_nFeature_range[1],
        max_features = input$spatial_nFeature_range[2],
        min_counts   = input$spatial_nCount_range[1],
        max_counts   = input$spatial_nCount_range[2],
        max_mt       = input$spatial_mt_max
      )
      keep_cells <- obj$nFeature_Spatial >= filters$min_features &
        obj$nFeature_Spatial <= filters$max_features &
        obj$nCount_Spatial   >= filters$min_counts &
        obj$nCount_Spatial   <= filters$max_counts &
        obj$percent.mt       <= filters$max_mt
      keep_cells[is.na(keep_cells)] <- FALSE
      # Nothing passes: stop and suggest data-driven lower bounds instead of applying other thresholds
      if (sum(keep_cells) == 0) {
        suggested <- analyze_data_distribution(obj)$suggested_filters
        stop(sprintf("No spot passes these thresholds. Data-driven lower bounds: genes >= %g, UMIs >= %g, MT <= %g%%",
                     suggested$min_features, suggested$min_counts, suggested$max_mt))
      }
      obj <- obj[, keep_cells]
      # Record the thresholds actually applied, for the report and reproducibility
      obj@misc$qc_params <- c(filters, list(tissue_filter = isTRUE(input$filter_tissue_spots), spots_before = original_spots, spots_after = ncol(obj)))
      update_spatial_object(obj, list(qc_applied = TRUE))
      filtered_spots <- ncol(obj)
      removal_rate <- round((1 - filtered_spots / original_spots) * 100, 1)
      success_message <- paste0(
        "QC completed!\n",
        "Remaining: ", format(filtered_spots, big.mark = ","), " spots (", 100 - removal_rate, "%)"
      )
      showNotification(success_message, type = "message", duration = 10)
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste("QC failed:", conditionMessage(e)[1]), type = "error")
    })
  }))
  
  ############################## Normalization ##############################
  
  # Normalize spatial data (LogNormalize + VST, or SCTransform)
  observeEvent(input$normalize_spatial_data, run_safely("normalize spatial data", {
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
                      format(n_spots, big.mark = ","), method, nvar_in))
      showModal(modalDialog(
        title = "Normalizing Data",
        div(
          paste("Normalizing", format(n_spots, big.mark = ","), "spots using", method, "..."),
          br(), tags$small("This may take several minutes")
        ),
        easyClose = FALSE, footer = NULL
      ))
      # Free memory and silence warnings during normalization
      gc()
      old_warn <- getOption("warn"); options(warn = -1)
      normalization_success <- FALSE
      method_used <- "none"
      if (identical(method, "SCT")) {
        message("=== PATH: SCTransform ===")
        # SCTransform on the Spatial assay creates the 'SCT' assay
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
        # Classic path; cap the number of variable features on very large datasets
        nfeatures_use <- if (n_spots > 10000) {
          message(sprintf("Large dataset detected (>10k). Using min(1000, %d) variable features.", nvar_in))
          min(1000L, as.integer(nvar_in))
        } else {
          as.integer(nvar_in)
        }
        # Normalize + variable features on Spatial, then ScaleData
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
      # Variable features plot
      message("Preparing Variable Features plot...")
      output$spatial_variable_features <- renderPlot({
        tryCatch({
          if (identical(method, "SCT")) {
            # SCT stores variable features in the SCT assay
            if ("SCT" %in% Seurat::Assays(obj) && length(Seurat::VariableFeatures(obj, assay = "SCT")) > 0) {
              Seurat::VariableFeaturePlot(obj, assay = "SCT") +
                ggtitle(paste("Variable Features —", method_used))
            } else {
              ggplot() + geom_text(aes(x = 0.5, y = 0.5, label = "No variable features (SCT)"), size = 6) + theme_void()
            }
          } else {
            # Classic path: variable features in the Spatial assay
            if (length(Seurat::VariableFeatures(obj, assay = "Spatial")) > 0) {
              Seurat::VariableFeaturePlot(obj, assay = "Spatial") +
                ggtitle(paste("Variable Features —", method_used))
            } else {
              ggplot() + geom_text(aes(x = 0.5, y = 0.5, label = "No variable features (Spatial)"), size = 6) + theme_void()
            }
          }
        }, error = function(e) {
          ggplot() + geom_text(aes(x = 0.5, y = 0.5, label = "Normalization completed"), size = 6) + theme_void()
        })
      })
      # Update object and notify
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
      showNotification(paste("Normalization failed:", conditionMessage(e)[1],
                             "\nTry applying stricter QC filters to reduce dataset size"),
                       type = "error")
    })
  }))
  
  ############################## PCA ##############################
  
  # Run PCA on the selected analysis assay
  observeEvent(input$spatial_scale_pca, run_safely("spatial scale pca", {
    req(spatial_obj())
    if (!processing_states$normalized) {
      showNotification("Please normalize the data first!", type = "warning")
      return()
    }
    tryCatch({
      showModal(modalDialog(title = "Running PCA", "Computing principal components...", easyClose = FALSE, footer = NULL))
      obj <- spatial_obj()
      # Assay chosen in the PCA card (kept as the analysis assay shown in the status box)
      assay_for_pca <- input$spatial_assay_select %||% selected_assay()
      if (!assay_for_pca %in% names(obj@assays)) stop(paste("Assay not found:", assay_for_pca))
      selected_assay(assay_for_pca)
      # Determine the scaling strategy
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
      # Number of PCs chosen in the PCA card, bounded by the numbers of spots and features
      n_spots <- ncol(obj)
      n_features <- nrow(obj[[assay_for_pca]])
      requested_pcs <- suppressWarnings(as.integer(input$spatial_pca_dims %||% 30))
      if (is.na(requested_pcs) || requested_pcs < 2) requested_pcs <- 30L
      max_pcs <- min(requested_pcs, n_spots - 1, n_features - 1)
      message(paste("Running PCA with", max_pcs, "components on", assay_for_pca, "assay"))
      # Seed-42 wrapper (dimensional_reduction_functions.R), parameters stored in @misc$pca_params
      obj <- runPCA_reproducible(obj, npcs = max_pcs, seed = 42, assay = assay_for_pca, verbose = FALSE)
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
      showNotification(paste("PCA failed:", conditionMessage(e)[1]), type = "error")
    })
  }))
  
  ############################## Neighbors, UMAP and clustering ##############################
  
  # Combined Neighbors + UMAP (seed-42 wrappers)
  observeEvent(input$spatial_run_neighbors_umap, run_safely("spatial run neighbors umap", {
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
      # Parameters from UI inputs
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
      # Step 1: find neighbors with reproducible seed
      message(paste("Step 1: Finding neighbors with", max_dims_neighbors, "dimensions..."))
      obj <- findNeighbors_reproducible(
        object = obj,
        dims = 1:max_dims_neighbors,
        seed = 42,
        verbose = FALSE
      )
      message("FindNeighbors computation completed")
      # Step 2: run UMAP with reproducible seed
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
      showNotification(paste("Neighbors + UMAP failed:", conditionMessage(e)[1]), type = "error", duration = 10)
    })
  }))
  
  # Find clusters on the SNN graph of the PCA assay (seed-42 wrapper)
  observeEvent(input$spatial_find_clusters, run_safely("spatial find clusters", {
    message("=== FIND CLUSTERS BUTTON CLICKED ===")
    if (is.null(spatial_obj())) {
      message("ERROR: spatial_obj() is NULL")
      showNotification("No spatial data loaded!", type = "error")
      return()
    }
    obj <- spatial_obj()
    # FindNeighbors stores <assay>_nn before <assay>_snn: clustering must run on the SNN graph
    snn_graphs <- grep("_snn$", names(obj@graphs), value = TRUE)
    if (length(snn_graphs) == 0) {
      message("ERROR: No neighbor graph found")
      showNotification("Please run 'Run Neighbors + UMAP' first!", type = "warning")
      return()
    }
    message("spatial_obj() exists and neighbors computed, proceeding...")
    tryCatch({
      # Get and validate parameters
      resolution <- as.numeric(input$spatial_resolution)
      algorithm <- as.integer(input$spatial_cluster_algorithm)
      if (is.na(resolution) || resolution <= 0) {
        resolution <- 0.5
        message("Invalid resolution, using default 0.5")
      }
      # Restrict to 1:3 — Leiden (4) requires the leidenalg Python package, not shipped
      if (is.na(algorithm) || !algorithm %in% 1:3) {
        algorithm <- 1L
        message("Invalid algorithm, using default Louvain")
      }
      algo_labels <- c("Louvain", "Louvain (multilevel)", "SLM")
      message(paste("Using resolution:", resolution, "and algorithm:", algo_labels[algorithm]))
      showModal(modalDialog(
        title = "Finding Clusters",
        div(
          paste("Clustering spots with resolution", resolution,
                "using", algo_labels[algorithm], "..."),
          br(),
          tags$small("This may take a few minutes for large datasets")
        ),
        easyClose = FALSE, footer = NULL
      ))
      message("=== STARTING FIND CLUSTERS ===")
      # Use the SNN graph built from the assay the PCA was computed on
      pca_assay <- tryCatch(DefaultAssay(obj[["pca"]]), error = function(e) DefaultAssay(obj))
      graph_name <- if (paste0(pca_assay, "_snn") %in% snn_graphs) paste0(pca_assay, "_snn") else snn_graphs[1]
      message(paste("Using graph:", graph_name))
      # Find clusters through the seed-42 wrapper (stores parameters in @misc$clustering_params)
      message("Starting FindClusters computation...")
      obj <- findClusters_reproducible(obj,
                                       resolution = resolution,
                                       algorithm = algorithm,
                                       seed = 42,
                                       graph.name = graph_name)
      message("FindClusters computation completed")
      # Ensure cluster identities are properly formatted
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        obj@meta.data$seurat_clusters <- as.factor(obj@meta.data$seurat_clusters)
        message("Cluster identities converted to factors")
      }
      # A new partition invalidates annotations and co-localization results built on the previous one
      if ("annotated_clusters" %in% colnames(obj@meta.data)) {
        obj$annotated_clusters <- NULL
        message("Previous cluster annotations removed (new partition)")
      }
      spatial_cluster_names(NULL)
      coloc_results(NULL)
      # Update reactive object and state
      update_spatial_object(obj, list(clustered = TRUE))
      message("Spatial object updated with clusters")
      n_clusters <- length(unique(obj$seurat_clusters))
      message(paste("Found", n_clusters, "clusters"))
      success_message <- paste0(
        "Clustering completed successfully!\n",
        "Found ", n_clusters, " clusters\n",
        "Resolution: ", resolution, "\n",
        "Algorithm: ", algo_labels[algorithm]
      )
      showNotification(success_message, type = "message", duration = 10)
      removeModal()
      message("=== FIND CLUSTERS COMPLETED ===")
    }, error = function(e) {
      message(paste("=== FIND CLUSTERS ERROR ===", e$message))
      removeModal()
      showNotification(paste("Clustering failed:", conditionMessage(e)[1]),
                       type = "error", duration = 12)
    })
  }))
  
  # Input validation for clustering parameters
  observe(run_safely("Interface update", persistent = FALSE, {
    # Validate clustering resolution
    if (!is.null(input$spatial_resolution)) {
      if (is.na(input$spatial_resolution) || input$spatial_resolution <= 0) {
        updateNumericInput(session, "spatial_resolution", value = 0.5)
        showNotification("Invalid resolution value, reset to 0.5", type = "warning", duration = 3)
      }
    }
    # Validate clustering algorithm (1:3 only — Leiden not supported)
    if (!is.null(input$spatial_cluster_algorithm)) {
      if (is.na(input$spatial_cluster_algorithm) || !input$spatial_cluster_algorithm %in% 1:3) {
        updateSelectInput(session, "spatial_cluster_algorithm", selected = 1)
        showNotification("Invalid algorithm selection, reset to Louvain", type = "warning", duration = 3)
      }
    }
  }))
  
  # Clustering plots after clustering (placeholders before)
  observe(run_safely("Interface update", persistent = FALSE, {
    req(spatial_obj())
    if (processing_states$clustered) {
      obj <- spatial_obj()
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        n_clusters <- length(unique(obj$seurat_clusters))
        # UMAP clustering plot
        output$spatial_umap_clusters <- renderPlot({
          if ("umap" %in% names(obj@reductions)) {
            p <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
                         label = isTRUE(input$spatial_plot_labels), pt.size = 0.8,
                         label.size = 4, raster = FALSE)
            if (!isTRUE(input$spatial_plot_labels)) p <- p + NoLegend()
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
        # Spatial clustering plot
        output$spatial_tissue_clusters <- renderPlot({
          if (length(obj@images) > 0) {
            tryCatch({
              SpatialDimPlot(obj, group.by = "seurat_clusters",
                             label = isTRUE(input$spatial_plot_labels),
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
      }
    } else {
      # Placeholder plots when not clustered
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
  }))
  # Spatial clustering plot with explicit handling of missing coordinates
  observe(run_safely("Interface update", persistent = FALSE, {
    req(spatial_obj())
    if (processing_states$clustered) {
      obj <- spatial_obj()
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        n_clusters <- length(unique(obj$seurat_clusters))
        output$spatial_tissue_clusters <- renderPlot({
          coord_check <- check_spatial_coordinates(obj)
          if (!coord_check$available) {
            # Informative message when coordinates are missing
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
                             label = isTRUE(input$spatial_plot_labels),
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
  }))
  
  # Processing status info box
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
  
  # Assays available in the object (status box of the clustering tab)
  output$spatial_assay_status <- renderInfoBox({
    assays <- if (!is.null(spatial_obj())) names(spatial_obj()@assays) else character(0)
    infoBox(title = "Available Assays", value = if (length(assays)) paste(assays, collapse = ", ") else "None", icon = icon("layer-group"), color = "teal", width = 12)
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
  
  ############################## Interactive tissue viewer ##############################
  # Single viewer replacing the former "Interactive Spatial Visualization" and "Spatial Tissue Viewer" tabs.
  # Server-side zoom: the tissue image stored in the object is cropped to the current window and redrawn with the spots it contains.
  # Helpers (get_tissue_image_spatial, crop_tissue_image, build_tissue_viewer_plot, clamp_view, zoom_center, pan_view)
  # live in spatial_functions.R. The full dataset is always drawn (no sketching for spatial views).
  
  # Printed once per session: confirms the viewer server code is the one running
  message("[viewer] tissue viewer server code loaded")
  
  # Best tissue image of the object, decoded once per image (QC, clustering or annotation update the object, not its image).
  # A decoding error is reported in the console and in the outputs instead of crashing the session through the reset observer.
  viewer_tissue <- reactive({
    req(spatial_obj())
    obj <- spatial_obj()
    shiny::validate(shiny::need(length(obj@images) > 0, "No tissue image in this object"))
    stored <- obj@misc$tissue_images
    key <- if (!is.null(stored$hires)) {
      paste("hires", length(stored$hires$png), stored$hires$scalef)
    } else if (!is.null(stored$lowres)) {
      paste("lowres", length(stored$lowres$png), stored$lowres$scalef)
    } else {
      paste("object", paste(dim(obj@images[[1]]@image), collapse = "x"), obj@images[[1]]@scale.factors$lowres)
    }
    if (!identical(key, viewer_cache$key)) {
      message("[viewer] decoding tissue image: ", key)
      t0 <- Sys.time()
      tissue <- tryCatch(get_tissue_image_spatial(obj), error = function(e) { message("[viewer] tissue image error: ", conditionMessage(e)); NULL })
      shiny::validate(shiny::need(!is.null(tissue), "The tissue image could not be decoded (see the R console)"))
      viewer_cache$tissue <- tissue
      viewer_cache$key <- key
      message(sprintf("Tissue viewer image: %s (%d x %d px, %.1f s)", tissue$source, as.integer(tissue$width_px), as.integer(tissue$height_px),
                      as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
    viewer_cache$tissue
  })
  
  # A new tissue image (new dataset) resets the window to the whole image
  observe(run_safely("Interface update", persistent = FALSE, {
    viewer_tissue()
    if (!identical(viewer_cache$view_key, viewer_cache$key)) {
      viewer_cache$view_key <- viewer_cache$key
      viewer_view(NULL)
      viewer_selected_spot(NULL)
    }
  }))
  
  # Current window in full-resolution pixels (NULL stored value = whole image)
  viewer_window <- reactive({
    tissue <- viewer_tissue()
    viewer_view() %||% c(0, tissue$width_full, 0, tissue$height_full)
  })
  
  # Spot positions (full-resolution pixels) and metadata used for coloring and the spot details panel
  viewer_spots <- reactive({
    req(spatial_obj())
    obj <- spatial_obj()
    coords <- tryCatch(get_spot_coordinates_spatial(obj), error = function(e) NULL)
    shiny::validate(shiny::need(!is.null(coords) && nrow(coords) > 0, "No spot coordinates in this object"))
    md <- obj@meta.data[coords$barcode, , drop = FALSE]
    cluster_col <- get_cluster_column_spatial(obj)
    clusters <- if (cluster_col %in% colnames(md)) factor(as.character(md[[cluster_col]]), levels = get_cluster_names_spatial(obj)) else factor(rep("No clusters", nrow(md)))
    data.frame(barcode = coords$barcode, x = coords$imagecol, y = coords$imagerow, cluster = clusters,
               nCount_Spatial = md$nCount_Spatial, nFeature_Spatial = md$nFeature_Spatial,
               percent.mt = if ("percent.mt" %in% colnames(md)) md$percent.mt else NA_real_, stringsAsFactors = FALSE)
  })
  
  # Values used to color the spots: clusters (factor), gene expression or a QC metric (numeric), or the categories of the
  # last two-gene co-localization analysis (factor with fixed colors)
  viewer_values <- reactive({
    spots <- viewer_spots()
    mode <- input$viewer_color_by %||% "cluster"
    if (mode == "cluster") return(list(value = spots$cluster, label = "Cluster"))
    if (mode == "gene") {
      gene <- input$viewer_gene
      shiny::validate(shiny::need(!is.null(gene) && nzchar(gene), "Select a gene to color the spots"))
      obj <- spatial_obj()
      assay <- DefaultAssay(obj)
      expr <- tryCatch(LayerData(obj, assay = assay, layer = "data")[gene, spots$barcode], error = function(e) NULL)
      shiny::validate(shiny::need(!is.null(expr), paste0("No normalized 'data' layer for ", gene, " in assay ", assay, ": normalize the data first")))
      return(list(value = as.numeric(expr), label = gene))
    }
    if (mode == "coloc") {
      r <- coloc_results()
      shiny::validate(shiny::need(!is.null(r), "Run a co-localization analysis first (Co-localization tab)"))
      shiny::validate(shiny::need(identical(r$analysis_type, "bidirectional"), "The co-localization view shows two-gene analyses only"))
      categories <- colocalization_categories(r$spatial_data, r$parameters, partners_only = isTRUE(input$viewer_coloc_partners_only))
      # Spots absent from the analysis (object changed since) are shown as Neither
      value <- unname(categories$category[spots$barcode])
      value[is.na(value)] <- categories$levels[1]
      return(list(value = factor(value, levels = categories$levels), colors = categories$colors,
                  label = sprintf("%s / %s (%g µm)", r$parameters$gene1, r$parameters$gene2, r$parameters$max_distance_um)))
    }
    shiny::validate(shiny::need(mode %in% colnames(spots), paste(mode, "is not available in this object")))
    return(list(value = spots[[mode]], label = mode))
  })
  
  # Spots inside the current window, with the plotted y (negated image row) used by the click handler.
  # Co-localization view: spots expressing neither gene are not drawn (the tissue stays visible), and the others are
  # sorted so that spots with a partner are drawn on top.
  viewer_visible <- reactive({
    spots <- viewer_spots()
    spots$value <- viewer_values()$value
    win <- viewer_window()
    spots <- spots[spots$x >= win[1] & spots$x <= win[2] & spots$y >= win[3] & spots$y <= win[4], , drop = FALSE]
    if (identical(input$viewer_color_by, "coloc")) {
      spots <- spots[as.integer(spots$value) > 1, , drop = FALSE]
      spots <- spots[order(spots$value), , drop = FALSE]
    }
    spots$y_plot <- -spots$y
    spots
  })
  
  # Cluster colors for the selected palette (stored cluster colors when they cover every cluster, Seurat hues otherwise)
  viewer_cluster_colors <- reactive({
    req(spatial_obj())
    obj <- spatial_obj()
    clusters <- get_cluster_names_spatial(obj)
    n_clusters <- length(clusters)
    if (n_clusters == 0) return(NULL)
    stored <- obj@misc$cluster_colors
    colors <- switch(input$viewer_palette %||% "default",
                     "default" = {
                       if (!is.null(names(stored)) && all(clusters %in% names(stored))) stored[clusters]
                       else if (length(stored) == n_clusters) stored
                       else scales::hue_pal()(n_clusters)
                     },
                     "polychrome" = {
                       if (n_clusters <= 36) {
                         Polychrome::palette36.colors(n_clusters)
                       } else {
                         # createPalette samples colors at random: seed 42 keeps the palette reproducible
                         set.seed(42)
                         Polychrome::createPalette(n_clusters, seedcolors = c("#ff0000", "#00ff00", "#0000ff"))
                       }
                     },
                     "set3" = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_clusters),
                     "paired" = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_clusters),
                     "dark2" = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_clusters),
                     scales::hue_pal()(n_clusters))
    return(setNames(unname(colors), clusters))
  })
  
  # Slider settings debounced: dragging a slider triggers one redraw when it stops, not one per step
  viewer_settings <- debounce(reactive({
    list(spot_size = input$viewer_spot_size %||% 1, spot_alpha = input$viewer_spot_alpha %||% 0.8,
         image_alpha = input$viewer_image_alpha %||% 1, plot_height = input$viewer_plot_height %||% 750)
  }), millis = 500)
  
  # Viewer plot for an output size given in points (screen: pixels at 72 dpi; download: inches x 72), shared by display and export.
  # With "Show spots" unticked only the tissue image is drawn: spots and coloring values are not even computed.
  viewer_plot_for_size <- function(width_pt, height_pt, display_size) {
    obj <- spatial_obj()
    settings <- viewer_settings()
    win <- viewer_window()
    show_spots <- !isFALSE(input$viewer_show_spots)
    spots <- if (show_spots) viewer_visible() else NULL
    vals <- if (show_spots) viewer_values() else list(label = "")
    crop <- crop_tissue_image(viewer_tissue(), win, display_size)
    spot_um <- obj@misc$spot_size_um %||% 55
    color_by_cluster <- identical(input$viewer_color_by %||% "cluster", "cluster")
    build_tissue_viewer_plot(
      spots = spots, view = win, crop = crop, value_label = vals$label,
      # Cluster palette, fixed co-localization colors, or the default scales
      cluster_colors = if (color_by_cluster) viewer_cluster_colors() else vals$colors,
      highlight = if (color_by_cluster && !identical(input$viewer_cluster_highlight, "None")) input$viewer_cluster_highlight else NULL,
      spot_diameter_px = obj@images[[1]]@scale.factors$spot %||% 10,
      # Visium HD bins are squares, standard Visium spots are discs
      spot_shape = if (spot_um < 50) 15 else 16,
      size_factor = settings$spot_size, spot_alpha = settings$spot_alpha, image_alpha = settings$image_alpha,
      show_labels = color_by_cluster && isTRUE(input$viewer_show_labels), show_legend = !isTRUE(input$viewer_remove_legend),
      panel_px = c(width_pt, height_pt)
    )
  }
  
  # Viewer plot; errors are printed in the console and shown in the plot area (req / validate messages pass through unchanged)
  output$viewer_plot <- renderPlot({
    req(spatial_obj())
    width_px <- session$clientData$output_viewer_plot_width %||% 1000
    message("[viewer] drawing plot: ", width_px, " x ", viewer_settings()$plot_height, " px")
    p <- tryCatch(viewer_plot_for_size(width_px, viewer_settings()$plot_height, display_size = 1600),
                  error = function(e) {
                    if (inherits(e, c("validation", "shiny.silent.error"))) stop(e)
                    message("[viewer] plot error: ", conditionMessage(e))
                    NULL
                  })
    shiny::validate(shiny::need(!is.null(p), "The viewer plot failed (see the R console)"))
    p
  }, height = function() viewer_settings()$plot_height)
  
  # Gene menu follows the genes of the default assay (server-side selectize: thousands of genes)
  observe(run_safely("Interface update", persistent = FALSE, {
    req(spatial_obj())
    genes <- rownames(spatial_obj()[[DefaultAssay(spatial_obj())]])
    if (identical(genes, last_viewer_genes)) return()
    last_viewer_genes <<- genes
    current <- isolate(input$viewer_gene)
    updateSelectizeInput(session, "viewer_gene", choices = genes, selected = if (!is.null(current) && current %in% genes) current else character(0), server = TRUE)
  }))
  
  # Highlight menu follows the active cluster names, keeping the current choice while it still exists
  observe(run_safely("Interface update", persistent = FALSE, {
    req(spatial_obj())
    clusters <- get_cluster_names_spatial(spatial_obj())
    message("[viewer] highlight menu: column '", get_cluster_column_spatial(spatial_obj()), "', ", length(clusters), " clusters")
    current <- isolate(input$viewer_cluster_highlight)
    updateSelectInput(session, "viewer_cluster_highlight", choices = c("None", clusters),
                      selected = if (!is.null(current) && current %in% clusters) current else "None")
  }))
  
  # Brush: zoom on the drawn rectangle (plot y is the negated image row)
  observeEvent(input$viewer_brush, run_safely("viewer brush", {
    b <- input$viewer_brush
    tissue <- viewer_tissue()
    viewer_view(clamp_view(c(b$xmin, b$xmax, -b$ymax, -b$ymin), tissue$width_full, tissue$height_full, min_frac = 0.002))
    session$resetBrush("viewer_brush")
  }))
  
  # Zoom by a factor around the window centre; zooming out past the whole image returns to the full view
  viewer_zoom <- function(factor) {
    tissue <- viewer_tissue()
    win <- viewer_window()
    new_view <- zoom_center(win, tissue$width_full, tissue$height_full, mean(win[1:2]), mean(win[3:4]), factor)
    new_view <- clamp_view(pan_view(new_view, tissue$width_full, tissue$height_full), tissue$width_full, tissue$height_full, min_frac = 0.002)
    full_view <- (new_view[2] - new_view[1]) >= tissue$width_full && (new_view[4] - new_view[3]) >= tissue$height_full
    viewer_view(if (full_view) NULL else new_view)
  }
  
  # Pan by a fraction of the current window size, staying inside the image
  viewer_pan <- function(dx_frac, dy_frac) {
    tissue <- viewer_tissue()
    win <- viewer_window()
    viewer_view(pan_view(win, tissue$width_full, tissue$height_full, dx_frac * (win[2] - win[1]), dy_frac * (win[4] - win[3])))
  }
  
  observeEvent(input$viewer_zoom_in, run_safely("viewer zoom in", viewer_zoom(2)))
  observeEvent(input$viewer_zoom_out, run_safely("viewer zoom out", viewer_zoom(0.5)))
  observeEvent(input$viewer_reset, run_safely("viewer reset", viewer_view(NULL)))
  observeEvent(input$viewer_dblclick, run_safely("viewer dblclick", viewer_view(NULL)))
  observeEvent(input$viewer_pan_left, run_safely("viewer pan left", viewer_pan(-0.4, 0)))
  observeEvent(input$viewer_pan_right, run_safely("viewer pan right", viewer_pan(0.4, 0)))
  observeEvent(input$viewer_pan_up, run_safely("viewer pan up", viewer_pan(0, -0.4)))
  observeEvent(input$viewer_pan_down, run_safely("viewer pan down", viewer_pan(0, 0.4)))
  
  # Click: details of the nearest spot (within 10 px); ignored in the tissue-only view
  observeEvent(input$viewer_click, run_safely("viewer click", {
    req(!isFALSE(input$viewer_show_spots))
    hit <- nearPoints(viewer_visible(), input$viewer_click, xvar = "x", yvar = "y_plot", threshold = 10, maxpoints = 1)
    viewer_selected_spot(if (nrow(hit) > 0) hit else NULL)
  }))
  
  output$viewer_spot_info <- renderText({
    s <- viewer_selected_spot()
    if (is.null(s)) return("Click a spot to see its details")
    obj <- spatial_obj()
    um_per_px <- tryCatch(get_um_per_pixel_spatial(obj, obj@misc$spot_size_um %||% 55), error = function(e) NA_real_)
    lines <- c(paste("Barcode:", s$barcode), paste("Cluster:", as.character(s$cluster)),
               paste("UMIs:", format(s$nCount_Spatial, big.mark = ",")), paste("Genes detected:", format(s$nFeature_Spatial, big.mark = ",")),
               sprintf("Mitochondrial: %.2f %%", s$percent.mt))
    if (!is.factor(s$value)) lines <- c(lines, sprintf("%s: %.3f", viewer_values()$label, s$value))
    if (identical(input$viewer_color_by, "coloc")) lines <- c(lines, paste("Co-localization:", as.character(s$value)))
    if (is.finite(um_per_px)) lines <- c(lines, sprintf("Position: %.0f, %.0f µm", s$x * um_per_px, s$y * um_per_px))
    paste(lines, collapse = "\n")
  })
  
  # Zoom level, window size in mm and number of spots in view.
  # Counted from the spot positions only, so the line stays visible while no gene is selected
  # (nrow(viewer_visible()) turned that validation message into an S4 dispatch error).
  output$viewer_view_info <- renderText({
    tissue <- viewer_tissue()
    win <- viewer_window()
    spots <- viewer_spots()
    n_spots <- sum(spots$x >= win[1] & spots$x <= win[2] & spots$y >= win[3] & spots$y <= win[4])
    obj <- spatial_obj()
    um_per_px <- tryCatch(get_um_per_pixel_spatial(obj, obj@misc$spot_size_um %||% 55), error = function(e) NA_real_)
    info <- sprintf("Zoom: x%.1f | %s spots in view", tissue$width_full / (win[2] - win[1]), format(n_spots, big.mark = ","))
    if (is.finite(um_per_px)) info <- paste0(info, sprintf(" | Window: %.2f x %.2f mm", (win[2] - win[1]) * um_per_px / 1000, (win[4] - win[3]) * um_per_px / 1000))
    return(info)
  })
  
  # Export of the current view with the current settings (12 in wide, height following the window shape).
  # req() / validate() failures (no gene selected, data not normalized) are silent in a downloadHandler: they are reported here.
  output$download_viewer_plot <- downloadHandler(
    filename = function() paste0("tissue_viewer_", Sys.Date(), ".", input$viewer_download_format %||% "png"),
    content = function(file) {
      win <- viewer_window()
      width_in <- 12
      height_in <- min(max(width_in * (win[4] - win[3]) / (win[2] - win[1]), 4), 16)
      p <- tryCatch(viewer_plot_for_size(width_in * 72, height_in * 72, display_size = 3000), error = function(e) {
        msg <- if (nzchar(conditionMessage(e))) conditionMessage(e) else "nothing to export yet"
        message("[viewer] export error: ", msg)
        showNotification(paste("Viewer export failed:", msg), type = "warning", duration = 8)
        stop(e)
      })
      ggsave(file, plot = p, width = width_in, height = height_in, dpi = 300, bg = "white", device = input$viewer_download_format %||% "png")
    }
  )
  
  ############################## Gene expression visualization ##############################
  
  # Cluster order choices follow the active cluster column (gene choices are handled by the observer below)
  observe(run_safely("Interface update", persistent = FALSE, {
    req(spatial_obj())
    obj <- spatial_obj()
    # Active cluster column (annotated names take precedence), the one the violin plot groups by
    cluster_choices <- get_cluster_names_spatial(obj)
    if (length(cluster_choices) > 0) {
      updateSelectInput(session, "spatial_cluster_order",
                        choices = cluster_choices,
                        selected = cluster_choices)
    }
  }))
  
  # Gene choices follow the visualization assay
  observe(run_safely("Interface update", persistent = FALSE, {
    req(spatial_obj())
    obj <- spatial_obj()
    all_genes <- character(0)
    current_assay <- input$spatial_viz_assay
    if (!is.null(current_assay) && current_assay %in% names(obj@assays)) {
      all_genes <- rownames(obj[[current_assay]])
    } else {
      # Fall back to the default assay genes
      all_genes <- rownames(obj[[DefaultAssay(obj)]])
    }
    message(paste("Updating gene choices. Available genes:", length(all_genes)))
    # Keep still-valid selected genes instead of clearing the picker on every update
    updatePickerInput(session, "spatial_genes_select",
                      choices = all_genes,
                      selected = intersect(isolate(input$spatial_genes_select), all_genes))
    # spatial_viz_assay is not updated here: doing so snapped it back to the default assay after every user change
    message(paste("Available assays:", paste(names(obj@assays), collapse = ", ")))
    message(paste("Default assay:", DefaultAssay(obj)))
  }))
  
  # Copy the picker selection into the three gene text inputs
  observeEvent(input$spatial_genes_select, run_safely("spatial genes select", {
    # Only update if the picker has actual selections
    if (!is.null(input$spatial_genes_select) && length(input$spatial_genes_select) > 0) {
      gene_text <- paste(input$spatial_genes_select, collapse = ", ")
      updateTextInput(session, "spatial_genes_text", value = gene_text)
      updateTextInput(session, "spatial_violin_genes_text", value = gene_text)
      updateTextInput(session, "spatial_dotplot_genes_text", value = gene_text)
    }
  }), ignoreInit = TRUE)
  
  # Display settings of the spatial feature plot, debounced: a burst of clicks or slider moves triggers a single redraw
  feature_display_params <- debounce(reactive({
    list(image_size = input$spatial_image_size %||% 1, pt_size = input$spatial_pt_size %||% 1.6, brightness = input$spatial_brightness %||% 1,
         alpha = input$spatial_alpha %||% 1, crop = input$spatial_crop %||% TRUE, combine = isTRUE(input$spatial_combine_plots),
         ncol = input$spatial_plot_ncol)
  }), millis = 1000)
  
  # Spatial feature plots (display settings read from feature_display_params, so they redraw without clicking again)
  observeEvent(input$show_spatial_features, run_safely("show spatial features", {
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
      selected_genes <- genes_from_text
      current_assay <- input$spatial_viz_assay %||% "Spatial"
      message(paste("=== SPATIAL FEATURE PLOT ==="))
      message(paste("Using assay:", current_assay))
      message(paste("Selected genes:", paste(selected_genes, collapse = ", ")))
      # Local copy whose default assay is the visualization assay: renderPlot reads it lazily, so it is never restored
      plot_obj <- spatial_obj()
      if (!current_assay %in% names(plot_obj@assays)) current_assay <- DefaultAssay(plot_obj)
      DefaultAssay(plot_obj) <- current_assay
      # Validate genes exist in the current assay
      available_genes <- rownames(plot_obj[[current_assay]])
      valid_genes <- selected_genes[selected_genes %in% available_genes]
      missing_genes <- selected_genes[!selected_genes %in% available_genes]
      if (length(missing_genes) > 0) {
        showNotification(paste("Missing genes:", paste(missing_genes, collapse = ", ")),
                         type = "warning", duration = 8)
      }
      if (length(valid_genes) == 0) {
        showNotification("None of the selected genes found in current assay", type = "error")
        return()
      }
      showModal(modalDialog(
        title = "Generating Spatial Feature Plots",
        paste("Creating plots for", length(valid_genes), "genes using", current_assay, "assay..."),
        easyClose = FALSE, footer = NULL
      ))
      # Spatial feature plot with dynamic sizing
      output$spatial_feature_plots <- renderPlot({
        # Debounced display settings (image size, spot size, brightness, combine)
        prm <- feature_display_params()
        if (length(plot_obj@images) == 0) {
          p <- ggplot() +
            geom_text(aes(x = 0.5, y = 0.5, label = "No spatial coordinates available"), size = 8) +
            theme_void()
          # Store the placeholder plot too
          spatial_feature_plot_reactive(p)
          return(p)
        } else {
          tryCatch({
            n_genes <- length(valid_genes)
            image_size_factor <- prm$image_size
            # Upper alpha bound capped at 1 (higher values have no visible effect)
            max_alpha <- min(1, prm$alpha * prm$brightness)
            ncol_val <- if (prm$combine && n_genes > 1) prm$ncol %||% min(3, n_genes) else 1
            message(paste("Creating spatial plot with", n_genes, "genes in", ncol_val, "columns"))
            if (n_genes > 1 && prm$combine) {
              # Multi-gene combined plot
              p <- SpatialFeaturePlot(object = plot_obj, features = valid_genes, crop = prm$crop, slot = "data", ncol = ncol_val, combine = TRUE,
                                      pt.size.factor = prm$pt_size * image_size_factor, alpha = c(0.1, max_alpha), image.alpha = 0.8,
                                      stroke = 0.25 * image_size_factor)
            } else {
              # One gene, or the first gene only when plots are not combined
              p <- SpatialFeaturePlot(object = plot_obj, features = valid_genes[1], crop = prm$crop, slot = "data",
                                      pt.size.factor = prm$pt_size * image_size_factor, alpha = c(0.1, max_alpha), image.alpha = 1,
                                      stroke = 0.25 * image_size_factor) +
                ggtitle(paste("Spatial Expression:", valid_genes[1], "(", current_assay, ")")) +
                theme(plot.title = element_text(hjust = 0.5, size = 14 * image_size_factor))
            }
            # Store the plot for downloads and the annotation tab
            spatial_feature_plot_reactive(p)
            return(p)
          }, error = function(e) {
            message(paste("Spatial plot error:", e$message))
            error_plot <- ggplot() +
              geom_text(aes(x = 0.5, y = 0.5, label = paste("Error:", substr(conditionMessage(e)[1], 1, 50))),
                        size = 4, color = "red") +
              theme_void()
            # Store the error plot too
            spatial_feature_plot_reactive(error_plot)
            return(error_plot)
          })
        }
      },
      height = function() {
        prm <- feature_display_params()
        size_factor <- prm$image_size
        n_genes <- length(valid_genes)
        # Margin for title + legend + padding
        title_legend_margin <- 150
        if (n_genes == 1 || !prm$combine) {
          # Single panel (one gene, or first gene only when plots are not combined)
          return(max(550, 600 * size_factor + title_legend_margin))
        } else {
          ncol_val <- prm$ncol %||% min(3, n_genes)
          nrow_val <- ceiling(n_genes / ncol_val)
          return(max(650, min(1500, 400 * nrow_val * size_factor + title_legend_margin)))
        }
      })
      # Store plot parameters for download
      spatial_plots$feature_plots <- list(
        genes = valid_genes,
        assay = current_assay,
        parameters = list(
          pt_size = (input$spatial_pt_size %||% 1.6) * (input$spatial_image_size %||% 1),
          alpha = min(1, (input$spatial_alpha %||% 1) * (input$spatial_brightness %||% 1)),
          crop = input$spatial_crop %||% TRUE,
          combine = input$spatial_combine_plots,
          image_size = input$spatial_image_size %||% 1
        )
      )
      removeModal()
      showNotification(paste("Spatial plots generated for", length(valid_genes), "genes"), type = "message")
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error generating spatial plots:", conditionMessage(e)[1]), type = "error")
    })
  }))
  
  # Dot plot by clusters
  observeEvent(input$show_spatial_dotplot, run_safely("show spatial dotplot", {
    if (is.null(spatial_obj())) { showNotification("Load a spatial dataset first.", type = "warning"); return() }
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
      check <- validate_gene_input(input$spatial_dotplot_genes_text, obj, input$spatial_viz_assay %||% "Spatial")
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      selected_assay <- check$assay
      valid_genes <- check$genes
      # Set the assay temporarily
      original_assay <- DefaultAssay(obj)
      DefaultAssay(obj) <- selected_assay
      # Generate plot
      p <- DotPlot(obj,
                   features = valid_genes,
                   group.by = cluster_col,
                   scale = input$spatial_dot_scale %||% TRUE,
                   scale.by = input$spatial_dot_scale_by %||% "radius",
                   dot.min = (input$spatial_dot_min %||% 0)/100,
                   dot.scale = input$spatial_dot_max %||% 6,
                   cluster.idents = FALSE,
                   cols = input$spatial_dot_colors %||% "RdYlBu") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5))
      # Flip axes when the checkbox is ticked
      if (isTRUE(input$spatial_dot_cluster_idents)) {
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
      showNotification(paste("Error generating dot plot:", conditionMessage(e)[1]), type = "error")
    })
  }))
  
  # Violin plot by clusters
  observeEvent(input$show_spatial_violin, run_safely("show spatial violin", {
    if (is.null(spatial_obj())) { showNotification("Load a spatial dataset first.", type = "warning"); return() }
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
      check <- validate_gene_input(input$spatial_violin_genes_text, obj, input$spatial_viz_assay %||% "Spatial")
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      selected_assay <- check$assay
      valid_genes <- check$genes
      # Set the assay temporarily
      original_assay <- DefaultAssay(obj)
      DefaultAssay(obj) <- selected_assay
      message(paste("Creating violin plot with genes:", paste(valid_genes, collapse = ", ")))
      message(paste("Using cluster column:", cluster_col))
      message(paste("Using assay:", selected_assay))
      p <- VlnPlot(
        object = obj,
        features = valid_genes,
        assay = selected_assay,          # Explicit assay
        group.by = cluster_col,          # Group by cluster column
        pt.size = if (input$spatial_violin_points) input$spatial_violin_pt_size else 0,
        log = input$spatial_violin_log,
        ncol = if (length(valid_genes) > 1) min(3, length(valid_genes)) else NULL,
        combine = TRUE,
        fill.by = "ident",               # Color by identity (clusters)
        layer = "data",                  # Normalized data layer
        flip = FALSE,                    # Standard orientation
        same.y.lims = FALSE,             # Independent y-axis per gene
        raster = FALSE                   # No rasterization for better quality
      )
      if (length(valid_genes) == 1) {
        # Single gene: enhanced appearance
        p <- p +
          ggtitle(paste("Expression of", valid_genes[1], "across", cluster_col)) +
          labs(
            x = paste("Clusters (", cluster_col, ")"),
            y = paste0("Expression Level", if (input$spatial_violin_log) " (log)" else ""),
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
        # Multiple genes: theme applied to every patchwork panel
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
        # Overall annotation
        p <- p + plot_annotation(
          title = paste("Gene Expression across", cluster_col),
          subtitle = paste("Assay:", selected_assay, "| Genes:", paste(valid_genes, collapse = ", ")),
          theme = theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12)
          )
        )
      }
      # Apply the cluster order if specified
      cluster_order <- input$spatial_cluster_order
      if (!is.null(cluster_order) && length(cluster_order) > 0) {
        current_clusters <- levels(as.factor(obj@meta.data[[cluster_col]]))
        valid_order <- cluster_order[cluster_order %in% current_clusters]
        if (length(valid_order) > 0) {
          # Clusters missing from the order are appended at the end
          missing_clusters <- setdiff(current_clusters, valid_order)
          final_order <- c(valid_order, missing_clusters)
          # & applies the scale to every panel of a multi-gene patchwork (+ would only change the last panel)
          p <- p & scale_x_discrete(limits = final_order)        }
      }
      # Store in reactive value for reuse
      spatial_violin_plot_reactive(p)
      output$spatial_violin_plot <- renderPlot({
        spatial_violin_plot_reactive()
      }, height = function() {
        # Height follows the genes actually plotted (3 columns, 350 px per row), not the raw text being edited
        350 * ceiling(max(1, length(valid_genes)) / 3)
      })
      # Restore original assay
      DefaultAssay(obj) <- original_assay
    }, error = function(e) {
      message(paste("Violin plot error:", e$message))
      showNotification(paste("Error generating violin plot:", conditionMessage(e)[1]), type = "error")
    })
  }))
  
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
        # Map format to a ggsave device name
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
  
  # Save the spatial object as RDS (the modal is closed on every exit path, errors included)
  output$save_spatial_object <- downloadHandler(
    filename = function() {
      paste0("spatial_object_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(spatial_obj())
      showModal(modalDialog(title = "Saving Spatial Object", "Saving your spatial data object...", easyClose = FALSE, footer = NULL))
      on.exit(removeModal(), add = TRUE)
      saveRDS(spatial_obj(), file)
      showNotification("Spatial object saved successfully!", type = "message")
    }
  )
  
  ############################## Cluster annotation ##############################
  
  # Initialize cluster names when clusters are available.
  # If the loaded object already has an "annotated_clusters" column (prior Cell-Hub session),
  # restore those names instead of resetting to Cluster_0, Cluster_1, ...
  observe(run_safely("Interface update", persistent = FALSE, {
    req(spatial_obj())
    if (processing_states$clustered && "seurat_clusters" %in% colnames(spatial_obj()@meta.data)) {
      obj <- spatial_obj()
      unique_clusters <- sort(unique(as.character(obj$seurat_clusters)))
      # Only initialize once (NULL check)
      if (is.null(spatial_cluster_names())) {
        if ("annotated_clusters" %in% colnames(obj@meta.data)) {
          # Map each seurat_clusters level to its annotated_clusters label
          restored <- setNames(
            sapply(unique_clusters, function(cl) {
              mask <- obj$seurat_clusters == cl
              if (any(mask)) as.character(obj$annotated_clusters[mask][1])
              else paste0("Cluster_", cl)
            }),
            unique_clusters
          )
          spatial_cluster_names(restored)
          message("Restored annotated cluster names from loaded object")
        } else {
          default_names <- setNames(paste0("Cluster_", unique_clusters), unique_clusters)
          spatial_cluster_names(default_names)
        }
      }
    }
  }))
  
  # Dynamic input fields for cluster renaming
  output$spatial_cluster_rename_inputs <- renderUI({
    req(spatial_obj())
    req(processing_states$clustered)
    obj <- spatial_obj()
    if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
      return(p("No clusters found. Please run clustering first.", style = "color: #666;"))
    }
    # Original cluster numbers
    unique_clusters <- sort(unique(as.character(obj$seurat_clusters)))
    # Current names (annotated or original)
    current_names <- spatial_cluster_names()
    has_annotated <- "annotated_clusters" %in% colnames(obj@meta.data)
    input_list <- lapply(unique_clusters, function(cluster) {
      # Current display name for this cluster
      if (has_annotated) {
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
      div(
        class = "form-group",
        style = "margin-bottom: 15px; padding: 10px; border: 1px solid #ddd; border-radius: 5px;",
        # Original cluster number
        tags$label(
          paste0("Original Cluster ", cluster),
          style = "font-weight: bold; color: #333; font-size: 12px;"
        ),
        # Current name when different from the default
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
  
  # Apply all cluster names (builds the annotated_clusters factor)
  observeEvent(input$apply_all_spatial_names, run_safely("apply all spatial names", {
    message("=== APPLY CLUSTER NAMES STARTED ===")
    req(spatial_obj())
    tryCatch({
      obj <- spatial_obj()
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
      spatial_cluster_names(new_names)
      message("Reactive values updated")
      # New annotated clusters column, relabelled level by level (identical names merge clusters)
      obj$annotated_clusters <- as.factor(obj$seurat_clusters)
      current_levels <- levels(obj$annotated_clusters)
      message(paste("Current levels:", paste(current_levels, collapse = ", ")))
      new_level_names <- sapply(current_levels, function(x) {
        if (x %in% names(new_names)) {
          return(new_names[[x]])
        } else {
          return(paste0("Cluster_", x))
        }
      })
      levels(obj$annotated_clusters) <- new_level_names
      message("New cluster names applied to factor levels")
      spatial_obj(obj)
      message("Spatial object updated")
      message("=== APPLY CLUSTER NAMES COMPLETED ===")
    }, error = function(e) {
      message(paste("=== APPLY CLUSTER NAMES ERROR ===", e$message))
      showNotification(paste("Error applying cluster names:", conditionMessage(e)[1]), type = "error", duration = 10)
    })
  }))
  
  # Reset cluster names to Cluster_<n> and drop the annotation column
  observeEvent(input$reset_spatial_names, run_safely("reset spatial names", {
    req(spatial_obj())
    obj <- spatial_obj()
    unique_clusters <- sort(unique(as.character(obj$seurat_clusters)))
    default_names <- setNames(paste0("Cluster_", unique_clusters), unique_clusters)
    spatial_cluster_names(default_names)
    if ("annotated_clusters" %in% colnames(obj@meta.data)) {
      obj$annotated_clusters <- NULL
      spatial_obj(obj)
    }
  }))
  
  # UMAP annotation plot, shared by the display and the download (the download no longer depends on the tab having been shown)
  spatial_umap_plot <- reactive({
    req(spatial_obj())
    obj <- spatial_obj()
    # Force reactivity on cluster names so the plot updates on rename
    spatial_cluster_names()
    cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) "annotated_clusters" else "seurat_clusters"
    has_umap <- "umap" %in% Seurat::Reductions(obj)
    shiny::validate(shiny::need(has_umap, "UMAP not found in this object."))
    p <- Seurat::DimPlot(obj, reduction = "umap", group.by = cluster_col, label = TRUE, label.size = 5, repel = TRUE, pt.size = 1, raster = FALSE) +
      ggtitle("UMAP - Annotated Clusters")
    if (isTRUE(input$spatial_remove_legend)) p <- p + NoLegend()
    if (isTRUE(input$spatial_remove_axes)) p <- p + NoAxes()
    if (isTRUE(input$anno_dark_mode)) p <- p + .dark_theme_bg()
    p
  })
  
  output$spatial_annotation_umap <- renderPlot({ spatial_umap_plot() })
  
  # Spatial tissue annotation plot, shared by the display and the download
  spatial_tissue_plot <- reactive({
    req(spatial_obj())
    obj <- spatial_obj()
    # Force reactivity on cluster names so the plot updates on rename
    spatial_cluster_names()
    cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) "annotated_clusters" else "seurat_clusters"
    if (length(obj@images) == 0) {
      p <- ggplot() + geom_text(aes(0.5, 0.5, label = "No spatial image"), size = 6) + theme_void()
      if (isTRUE(input$anno_dark_mode)) {
        p <- p + theme(plot.background = element_rect(fill = "black", color = NA), panel.background = element_rect(fill = "black", color = NA),
                       text = element_text(color = "white"))
      }
      return(p)
    }
    # Guard all input$ refs — annotation tab may not have been visited yet
    show_labels <- isTRUE(input$spatial_anno_show_labels)
    if (!is.null(input$spatial_cluster_highlight) && input$spatial_cluster_highlight != "None") {
      obj_highlight <- obj
      obj_highlight@meta.data$highlight_group <- ifelse(obj_highlight@meta.data[[cluster_col]] == input$spatial_cluster_highlight,
                                                        as.character(obj_highlight@meta.data[[cluster_col]]), "Other")
      selected_cluster <- input$spatial_cluster_highlight
      colors <- c("Other" = "lightgray")
      colors[selected_cluster] <- "red"
      p <- Seurat::SpatialDimPlot(obj_highlight, group.by = "highlight_group", label = show_labels, label.size = 4, pt.size.factor = 1.6,
                                  alpha = c(0.3, 1), cols = colors) +
        ggtitle(paste("Spatial View - Cluster", selected_cluster, "Highlighted"))
    } else {
      cluster_colors <- if (!is.null(obj@misc$cluster_colors)) obj@misc$cluster_colors else NULL
      p <- Seurat::SpatialDimPlot(obj, group.by = cluster_col, label = show_labels, label.size = 4, pt.size.factor = 1.6, alpha = c(0.3, 1), cols = cluster_colors) +
        ggtitle("Spatial View - Annotated Clusters")
    }
    if (isTRUE(input$spatial_remove_legend)) p <- p + NoLegend()
    if (isTRUE(input$spatial_remove_axes)) p <- p + NoAxes()
    if (isTRUE(input$anno_dark_mode)) p <- p + .dark_theme_bg()
    p
  })
  
  output$spatial_annotation_tissue <- renderPlot({ spatial_tissue_plot() })
  
  # UMAP annotation download
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
  
  # Export the annotation tissue plot exactly as displayed (stored by its renderPlot)
  output$download_spatial_tissue <- downloadHandler(
    filename = function() paste0("spatial_tissue.", input$spatial_plot_format %||% "png"),
    content = function(file) {
      p <- spatial_tissue_plot(); req(p)
      ggsave(file, plot = p, width = 16, height = 14, dpi = 300, bg = if (isTRUE(input$anno_dark_mode)) "black" else "white", device = input$spatial_plot_format %||% "png")
    }
  )
  
  # Cluster highlight menus follow the annotated names
  observe(run_safely("Interface update", persistent = FALSE, {
    req(spatial_obj())
    obj <- spatial_obj()
    cluster_col <- get_cluster_column_spatial(obj)
    if (cluster_col %in% colnames(obj@meta.data)) {
      # Annotated names if available, otherwise original cluster numbers
      cluster_choices <- get_cluster_names_spatial(obj)
      updateSelectInput(session, "spatial_cluster_highlight",
                        choices = c("None", cluster_choices),
                        selected = "None")
      updateSelectInput(session, "spatial_plotly_cluster_highlight",
                        choices = c("None", cluster_choices),
                        selected = "None")
    }
  }))
  
  # Show only the panel of the selected alternative plot
  observeEvent(input$annotation_selected_plot, run_safely("annotation selected plot", {
    # Hide all panels
    shinyjs::hide("no_plot_panel")
    shinyjs::hide("violin_plot_panel")
    shinyjs::hide("dot_plot_panel")
    shinyjs::hide("feature_plot_panel")
    # Show the selected panel only
    if (input$annotation_selected_plot == "none") {
      shinyjs::show("no_plot_panel")
    } else if (input$annotation_selected_plot == "violin") {
      shinyjs::show("violin_plot_panel")
    } else if (input$annotation_selected_plot == "dot") {
      shinyjs::show("dot_plot_panel")
    } else if (input$annotation_selected_plot == "feature") {
      shinyjs::show("feature_plot_panel")
    }
  }))
  
  
  output$annotation_feature_display <- renderPlot({
    plot_obj <- spatial_feature_plot_reactive()
    if (!is.null(plot_obj)) {
      plot_obj
    } else {
      # Empty plot with message
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
  
  # Save the annotated spatial object
  output$save_annotated_spatial_object <- createDownloadHandler(
    reactive_data = spatial_obj,
    object_name_reactive = reactive({
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
  
  ############################## Gene co-localization ##############################
  # Analysis helpers (analyze_gene_colocalization, calculate_cluster_summary, calculate_inter_cluster_stats,
  # get_spot_positions_um, create_distance_bar) live in spatial_functions.R.
  
  # Run analysis
  observeEvent(input$run_coloc_analysis, {
    run_safely("Co-localization", {
      req(spatial_obj(), input$coloc_genes)
      # Parse comma-separated genes
      gene_list <- strsplit(input$coloc_genes, ",")[[1]]
      gene_list <- trimws(gene_list)
      gene_list <- gene_list[gene_list != ""]
      if (length(gene_list) < 2) {
        showNotification("Please enter at least 2 genes", type = "warning")
        return()
      }
      # Check existence
      available_genes <- rownames(spatial_obj()[[input$coloc_assay]])
      missing_genes <- gene_list[!gene_list %in% available_genes]
      if (length(missing_genes) > 0) {
        showNotification(paste("Genes not found:", paste(missing_genes, collapse = ", ")), type = "error")
        return()
      }
      # Physical spot / bin diameter: user value, else the one stored at load time, else the Visium spot (55 um)
      spot_size_um <- input$coloc_spot_size_um
      if (is.null(spot_size_um) || is.na(spot_size_um) || spot_size_um <= 0) spot_size_um <- spatial_obj()@misc$spot_size_um %||% 55
      # Permutations only apply to two genes; invalid values switch the test off
      n_perm <- suppressWarnings(as.integer(input$coloc_n_permutations %||% 0))
      if (is.na(n_perm) || n_perm < 0) n_perm <- 0L
      if (length(gene_list) > 2 && n_perm > 0) showNotification("The permutation test is available for two genes only: it is skipped for 3+ genes", type = "warning", duration = 6)
      showModal(modalDialog(title = "Running Co-localization Analysis",
                            p(if (n_perm > 0 && length(gene_list) == 2) sprintf("Analyzing, with %d permutations...", n_perm) else "Analyzing..."),
                            easyClose = FALSE, footer = NULL))
      # finally closes the modal; the error itself is reported by run_safely
      tryCatch({
        results <- analyze_gene_colocalization(
          sp_obj = spatial_obj(),
          gene_list = gene_list,
          max_distance_um = input$coloc_distance_um,
          spot_size_um = spot_size_um,
          assay = input$coloc_assay,
          bidirectional = isTRUE(input$coloc_bidirectional),
          n_permutations = if (length(gene_list) == 2) n_perm else 0L,
          seed = 42
        )
        coloc_results(results)
      }, finally = removeModal())
    })
  })
  # Results-ready flag (drives the conditionalPanel in the UI)
  output$coloc_results_ready <- reactive({ !is.null(coloc_results()) })
  outputOptions(output, "coloc_results_ready", suspendWhenHidden = FALSE)
  
  # Calibration check shown next to the distance parameters: the grid pitch must match the technology
  output$detected_resolution_info <- renderText({
    req(spatial_obj())
    spot_um <- input$coloc_spot_size_um
    if (is.null(spot_um) || is.na(spot_um) || spot_um <= 0) return("Enter the spot / bin size in µm")
    tryCatch({
      pos <- get_spot_positions_um(spatial_obj(), spot_um)
      sprintf("Scale: %.3f µm/pixel\nGrid pitch: %.1f µm (expected ~100 µm for Visium, ~bin size for Visium HD)", pos$um_per_px, pos$pitch_um)
    }, error = function(e) paste("Calibration unavailable:", conditionMessage(e)))
  })
  
  # Quick stats summary with full explanations (and the permutation test when it was run)
  output$coloc_quick_stats <- renderText({
    req(coloc_results())
    r <- coloc_results()
    if (r$analysis_type == "simultaneous") {
      genes <- r$parameters$gene_list
      threshold <- r$parameters$max_distance_um
      spot_size <- r$parameters$spot_size_um
      header <- sprintf("=== Multi-gene Co-localization Analysis ===\nGenes: %s\nDistance threshold: %g µm | Spot/bin size: %g µm\n\n",
                        paste(genes, collapse = " + "), threshold, spot_size)
      gene_lines <- sapply(genes, function(g) {
        total <- r$total_spots_per_gene[g]
        coloc <- r$colocalized_spots_per_gene[g]
        rate <- if (total > 0) round(100 * coloc / total, 1) else 0
        sprintf("%s: %d total spots, %d co-localized (%0.1f%% of %s spots are near all other genes within %g µm)",
                g, total, coloc, rate, g, threshold)
      })
      paste(c(header, gene_lines, "\nDescriptive measure: no permutation test for 3+ genes."), collapse = "\n")
    } else {
      params <- r$parameters
      gene1 <- params$gene1
      gene2 <- params$gene2
      threshold <- params$max_distance_um
      spot_size <- params$spot_size_um
      header <- sprintf("=== Co-localization Analysis ===\nGenes: %s vs %s\nDistance threshold: %g µm | Spot/bin size: %g µm | Grid pitch: %g µm\n",
                        gene1, gene2, threshold, spot_size, params$pitch_um)
      # Below the grid pitch no neighbouring spot is reachable: only co-expression within the same spot is counted
      if (is.finite(params$pitch_um) && threshold < params$pitch_um) {
        header <- paste0(header, sprintf("WARNING: the distance (%g µm) is smaller than the distance between neighbouring spots (%g µm): only spots expressing both genes are counted. Use at least %g µm to include neighbours.\n",
                                         threshold, params$pitch_um, ceiling(params$pitch_um)))
      }
      counts <- sprintf("\n%s: %s spots\n%s: %s spots\n",
                        gene1, format(r$gene1_total_spots, big.mark = ","),
                        gene2, format(r$gene2_total_spots, big.mark = ","))
      gene1_rate <- if (!is.na(r$gene1_coloc_rate)) r$gene1_coloc_rate else 0
      gene1_text <- sprintf("\n%0.1f%% of %s spots have at least one %s neighbor within %g µm\n(%d out of %s %s spots)",
                            gene1_rate, gene1, gene2, threshold,
                            r$gene1_near_gene2, format(r$gene1_total_spots, big.mark = ","), gene1)
      gene2_text <- ""
      if (isTRUE(params$bidirectional) && !is.na(r$gene2_coloc_rate)) {
        gene2_rate <- r$gene2_coloc_rate
        gene2_text <- sprintf("\n\n%0.1f%% of %s spots have at least one %s neighbor within %g µm\n(%d out of %s %s spots)",
                              gene2_rate, gene2, gene1, threshold,
                              r$gene2_near_gene1, format(r$gene2_total_spots, big.mark = ","), gene2)
      }
      # Permutation test: observed vs expected rate when the target gene is placed at random among all spots
      perm <- r$permutation
      perm_text <- if (is.null(perm)) {
        "\n\nNo permutation test (set the number of permutations above 0 to test the association)."
      } else {
        format_perm <- function(x, ref, tgt) sprintf("%s -> %s: observed %.1f%%, expected %.1f%%, enrichment %s, empirical p = %.3g",
                                                     ref, tgt, x$observed_pct, x$expected_pct, ifelse(is.na(x$enrichment), "NA", format(x$enrichment)), x$p_value)
        lines <- format_perm(perm$gene1_to_gene2, gene1, gene2)
        if (!is.null(perm$gene2_to_gene1)) lines <- c(lines, format_perm(perm$gene2_to_gene1, gene2, gene1))
        paste0("\n\n=== Spatial association test (", perm$n_permutations, " permutations, seed ", perm$seed, ") ===\n",
               paste(lines, collapse = "\n"),
               "\nThe percentages above measure proximity (is the second gene available near the first?). This test asks a different question:",
               "\nare the two genes found together MORE than expected from their abundance, if the second gene were placed at random?",
               "\nA high p-value does not mean the genes are not close: a widely expressed gene is close to everything, so its proximity is not specific.",
               "\nMinimum attainable p = ", signif(1 / (perm$n_permutations + 1), 3), ". Spot-level test: spots are not independent cells.")
      }
      paste0(header, counts, gene1_text, gene2_text, perm_text)
    }
  })
  
  # Cluster-level co-localization table with visual bars
  output$coloc_cluster_table <- renderDT({
    req(coloc_results())
    r <- coloc_results()
    if (r$analysis_type == "simultaneous") {
      genes <- r$parameters$gene_list
      percentages <- sapply(genes, function(g) {
        total <- r$total_spots_per_gene[g]
        coloc <- r$colocalized_spots_per_gene[g]
        if (total > 0 && !is.na(total) && !is.na(coloc)) {
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
        # create_distance_bar handles NA internally
        out <- out %>%
          mutate(
            intra_mean_bar = sapply(intra_mean_dist, function(d) create_distance_bar(d, max_dist_threshold)),
            intra_median_bar = sapply(intra_median_dist, function(d) create_distance_bar(d, max_dist_threshold))
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
            # Bar columns: wide (150px)
            list(width = '150px', targets = bar_cols - 1),
            # Number columns: narrow (70px)
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
  
  # Inter-cluster pair choices (intra-cluster pairs excluded)
  observe(run_safely("Interface update", persistent = FALSE, {
    req(coloc_results())
    inter <- coloc_results()$inter_cluster_gene1_to_gene2
    if (!is.null(inter) && nrow(inter) > 0 && "ref_cluster" %in% names(inter)) {
      # Filter out intra-cluster pairs (same → same)
      inter_only <- inter %>%
        filter(ref_cluster != target_cluster)
      if (nrow(inter_only) > 0) {
        # Unique pair labels
        pair_df <- inter_only %>%
          select(ref_cluster, target_cluster) %>%
          distinct() %>%
          arrange(ref_cluster, target_cluster)
        pair_labels <- paste(pair_df$ref_cluster, "→", pair_df$target_cluster)
        names(pair_labels) <- pair_labels
        updateSelectInput(session, "inter_cluster_pairs_select",
                          choices = pair_labels,
                          selected = NULL)
      } else {
        # No inter-cluster pairs, only intra-cluster ones
        updateSelectInput(session, "inter_cluster_pairs_select",
                          choices = c("No inter-cluster pairs found" = ""),
                          selected = NULL)
      }
    }
  }))
  
  # Inter-cluster table (both directions for each cluster pair)
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
    # Filter out intra-cluster rows
    g1g2 <- inter_g1_to_g2 %>% filter(ref_cluster != target_cluster)
    g2g1 <- if (!is.null(inter_g2_to_g1) && "gene1_pct" %in% colnames(inter_g2_to_g1)) {
      inter_g2_to_g1 %>% filter(ref_cluster != target_cluster)
    } else {
      data.frame()
    }
    if (nrow(g1g2) == 0) {
      return(datatable(data.frame(Message = "No inter-cluster interactions")))
    }
    # Unified table: one row per (A, B) pair with both genes in both directions
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
        n_pairs_AB = if (nrow(row_AB_g1g2) > 0) row_AB_g1g2$n_pairs[1] else 0,
        g1A_to_g2B_pct = if (nrow(row_AB_g1g2) > 0) row_AB_g1g2$gene1_pct[1] else 0,
        g2B_to_g1A_pct = if (nrow(row_BA_g2g1) > 0) row_BA_g2g1$gene1_pct[1] else 0,
        n_pairs_BA = if (nrow(row_BA_g1g2) > 0) row_BA_g1g2$n_pairs[1] else 0,
        g1B_to_g2A_pct = if (nrow(row_BA_g1g2) > 0) row_BA_g1g2$gene1_pct[1] else 0,
        g2A_to_g1B_pct = if (nrow(row_AB_g2g1) > 0) row_AB_g2g1$gene1_pct[1] else 0,
        mean_dist = if (nrow(row_AB_g1g2) > 0) row_AB_g1g2$mean_distance_um[1] else NA,
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
  
  # Reset the inter-cluster pair filter
  observeEvent(input$reset_cluster_filter, run_safely("reset cluster filter", {
    updateCheckboxInput(session, "inter_cluster_show_all", value = TRUE)
    updateSelectInput(session, "inter_cluster_pairs_select", selected = NULL)
  }))
  
  # Co-localization map: tissue image with the spots of the last two-gene analysis coloured by category.
  # Uses the tissue image and spot positions of the interactive viewer, so both views stay aligned.
  coloc_map_plot <- reactive({
    r <- coloc_results(); req(r)
    shiny::validate(shiny::need(identical(r$analysis_type, "bidirectional"), "The map shows two-gene analyses: see the bar chart for 3+ genes"))
    spots <- viewer_spots()[, c("barcode", "x", "y")]
    tissue <- tryCatch(viewer_tissue(), error = function(e) NULL)
    view <- if (!is.null(tissue)) c(0, tissue$width_full, 0, tissue$height_full) else c(range(spots$x), range(spots$y))
    crop <- if (!is.null(tissue)) crop_tissue_image(tissue, view, display_size = 2000) else NULL
    categories <- colocalization_categories(r$spatial_data, r$parameters, partners_only = isTRUE(input$coloc_map_partners_only))
    p <- r$parameters
    # One p-value per direction that was tested (gene2 -> gene1 is absent when the analysis is not bidirectional)
    perm_p <- if (!is.null(r$permutation)) vapply(Filter(Negate(is.null), list(r$permutation$gene1_to_gene2, r$permutation$gene2_to_gene1)), function(x) sprintf("%.3g", x$p_value), character(1)) else character(0)
    subtitle <- sprintf("Distance ≤ %g µm (bin / spot %g µm)%s", p$max_distance_um, p$spot_size_um,
                        if (length(perm_p) > 0) paste0(" | permutation p = ", paste(perm_p, collapse = " / ")) else "")
    build_colocalization_map(spots, categories, crop = crop, view = view, point_size = input$coloc_map_point_size %||% 0.8,
                             title = paste(p$gene1, "/", p$gene2), subtitle = subtitle)
  })
  
  output$coloc_map <- renderPlot({
    coloc_map_plot()
  }, res = 96)
  
  # Where the pairs are: spots with a partner per cluster and direction (two genes), or positive / co-localized spots per gene (3+ genes)
  coloc_cluster_bar_plot <- reactive({
    r <- coloc_results(); req(r)
    if (r$analysis_type == "simultaneous") {
      genes <- r$parameters$gene_list
      d <- data.frame(gene = factor(rep(genes, 2), levels = genes), type = rep(c("Positive spots", "Co-localized spots"), each = length(genes)),
                      n = c(unname(r$total_spots_per_gene), unname(r$colocalized_spots_per_gene)))
      return(ggplot(d, aes(x = gene, y = n, fill = type)) + geom_col(position = "dodge", width = 0.7) +
               scale_fill_manual(values = c("Positive spots" = "#9ECAE1", "Co-localized spots" = "#08519C"), name = NULL) +
               scale_y_continuous(labels = scales::comma) + labs(x = NULL, y = "Number of spots", title = sprintf("Distance ≤ %g µm", r$parameters$max_distance_um)) +
               theme_minimal(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5), plot.background = element_rect(fill = "white", colour = NA)))
    }
    obj <- spatial_obj()
    cluster_col <- get_cluster_column_spatial(obj)
    shiny::validate(shiny::need(cluster_col %in% colnames(obj@meta.data), "No clusters in this object"))
    clusters <- obj@meta.data[[cluster_col]]
    if (!is.factor(clusters)) clusters <- factor(as.character(clusters), levels = get_cluster_names_spatial(obj))
    names(clusters) <- colnames(obj)
    d <- colocalization_by_cluster(r$spatial_data, r$parameters, clusters)
    shiny::validate(shiny::need(nrow(d) > 0, "No spot has a partner within the distance"))
    p <- r$parameters
    colors <- stats::setNames(c("#08519C", "#A50F15"), c(paste(p$gene1, "next to", p$gene2), paste(p$gene2, "next to", p$gene1)))
    ggplot(d, aes(x = n_spots, y = cluster, fill = direction)) + geom_col(position = "dodge", width = 0.75) +
      scale_fill_manual(values = colors, name = NULL) + scale_y_discrete(limits = rev) +
      labs(x = "Spots with a partner within the distance", y = NULL, title = "Where the pairs are") +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom", legend.direction = "vertical", plot.title = element_text(face = "bold", hjust = 0.5),
                                            panel.grid.major.y = element_blank(), plot.background = element_rect(fill = "white", colour = NA))
  })
  
  output$coloc_cluster_bar <- renderPlot({
    coloc_cluster_bar_plot()
  }, res = 96)
  
  output$download_coloc_map <- downloadHandler(
    filename = function() { r <- coloc_results(); paste0("colocalization_map_", r$parameters$gene1, "_", r$parameters$gene2, "_", Sys.Date(), ".png") },
    content = function(file) { ggsave(file, coloc_map_plot(), width = 12, height = 9, dpi = 300, bg = "white") }
  )
  
  output$download_coloc_bar <- downloadHandler(
    filename = function() paste0("colocalization_by_cluster_", Sys.Date(), ".png"),
    content = function(file) { ggsave(file, coloc_cluster_bar_plot(), width = 6, height = 8, dpi = 300, bg = "white") }
  )
  
  # Multi-sheet Excel export of the co-localization results (pair tables only when they fit the row limit)
  output$download_coloc_results <- downloadHandler(
    filename = function() paste0("colocalization_", Sys.Date(), ".xlsx"),
    content = function(file) {
      r <- coloc_results(); req(r)
      sheets <- list(parameters = data.frame(parameter = names(r$parameters), value = vapply(r$parameters, function(v) paste(v, collapse = ", "), character(1)), stringsAsFactors = FALSE))
      if (r$analysis_type == "simultaneous") {
        sheets$summary <- data.frame(gene = names(r$total_spots_per_gene), positive_spots = unname(r$total_spots_per_gene),
                                     colocalized_spots = unname(r$colocalized_spots_per_gene), colocalization_rate_pct = unname(r$colocalization_rates))
      } else {
        if (!is.null(r$permutation)) {
          directions <- Filter(Negate(is.null), list(gene1_to_gene2 = r$permutation$gene1_to_gene2, gene2_to_gene1 = r$permutation$gene2_to_gene1))
          sheets$permutation_test <- data.frame(direction = names(directions),
                                                observed_pct = vapply(directions, function(x) x$observed_pct, numeric(1)),
                                                expected_pct = vapply(directions, function(x) x$expected_pct, numeric(1)),
                                                enrichment = vapply(directions, function(x) as.numeric(x$enrichment), numeric(1)),
                                                empirical_p = vapply(directions, function(x) x$p_value, numeric(1)),
                                                n_permutations = r$permutation$n_permutations, seed = r$permutation$seed, stringsAsFactors = FALSE)
        }
        sheets$cluster_summary <- r$cluster_summary
        sheets$inter_cluster_g1_to_g2 <- r$inter_cluster_gene1_to_gene2
        sheets$inter_cluster_g2_to_g1 <- r$inter_cluster_gene2_to_gene1
        sheets$pairs_g1_to_g2 <- r$cluster_pairs_gene1_to_gene2
        sheets$pairs_g2_to_g1 <- r$cluster_pairs_gene2_to_gene1
      }
      sheets <- Filter(function(s) is.data.frame(s) && ncol(s) > 0, sheets)
      openxlsx::write.xlsx(sheets, file)
    }
  )
  
  ############################## Spatially variable genes (Moran's I) ##############################
  # Helpers (find_spatially_variable_genes, build_spatial_knn_weights, compute_morans_i) live in spatial_functions.R.
  # Always computed on the full dataset (no sketching for spatial analyses).
  
  # Assay menu follows the object, keeping the current choice while it exists
  observe({
    run_safely("Spatially variable genes assay menu", persistent = FALSE, {
      req(spatial_obj())
      assays <- names(spatial_obj()@assays)
      current <- isolate(input$svg_assay)
      updateSelectInput(session, "svg_assay", choices = assays, selected = if (!is.null(current) && current %in% assays) current else DefaultAssay(spatial_obj()))
    })
  })
  
  observeEvent(input$run_svg, {
    run_safely("Spatially variable genes", {
      req(spatial_obj())
      if (!processing_states$normalized) {
        showNotification("Please normalize the data first!", type = "warning")
        return()
      }
      showModal(modalDialog(title = "Spatially variable genes", "Computing Moran's I on the full dataset...", easyClose = FALSE, footer = NULL))
      # finally closes the modal; the error itself is reported by run_safely
      tryCatch({
        res <- find_spatially_variable_genes(spatial_obj(), assay = input$svg_assay %||% DefaultAssay(spatial_obj()),
                                             feature_source = input$svg_feature_source %||% "variable",
                                             n_features = input$svg_n_features %||% 2000, k = input$svg_k %||% 6,
                                             min_pct = (input$svg_min_pct %||% 1) / 100, log_function = message)
        svg_results(res)
        # Parameters kept in the object for the report and reproducibility
        obj <- spatial_obj()
        obj@misc$spatially_variable_params <- attr(res, "params")
        spatial_obj(obj)
        updateSelectizeInput(session, "svg_plot_gene", choices = res$gene, selected = res$gene[1], server = TRUE)
        showNotification(sprintf("Moran's I computed for %d genes: %d with adjusted p < 0.05", nrow(res), sum(res$p_val_adj < 0.05)),
                         type = "message", duration = 8)
      }, finally = removeModal())
    })
  })
  
  # Results table: p-values stored raw, formatted only for display (p-values below the smallest representable double are shown as "< 1e-300")
  output$svg_table <- renderDT({
    res <- svg_results()
    if (is.null(res)) return(datatable(data.frame(Message = "Click 'Find spatially variable genes'"), options = list(dom = "t"), rownames = FALSE))
    show <- res
    show$morans_I <- round(show$morans_I, 4); show$z_score <- round(show$z_score, 2); show$mean_expr <- round(show$mean_expr, 3)
    format_p <- function(p) ifelse(p == 0, "< 1e-300", formatC(p, format = "e", digits = 2))
    show$p_val <- format_p(show$p_val); show$p_val_adj <- format_p(show$p_val_adj)
    p <- attr(res, "params")
    datatable(show, filter = "top", selection = "single", rownames = FALSE,
              caption = sprintf("%d genes tested on %s spots (%s assay, %d-neighbour graph). Rank by Moran's I: with many spots, p-values are tiny for any patterned gene.",
                                p$n_tested, format(p$n_spots, big.mark = ","), p$assay, p$k),
              options = list(pageLength = 20, scrollX = TRUE))
  })
  
  # Clicking a row shows that gene on the tissue
  observeEvent(input$svg_table_rows_selected, {
    run_safely("Spatially variable genes row selection", persistent = FALSE, {
      res <- svg_results(); req(res)
      gene <- res$gene[input$svg_table_rows_selected[1]]
      updateSelectizeInput(session, "svg_plot_gene", choices = res$gene, selected = gene, server = TRUE)
    })
  })
  
  svg_plot_reactive <- reactive({
    req(spatial_obj(), input$svg_plot_gene)
    obj <- spatial_obj()
    assay <- input$svg_assay %||% DefaultAssay(obj)
    shiny::validate(shiny::need(length(obj@images) > 0, "No tissue image in this object"),
                    shiny::need(assay %in% names(obj@assays) && input$svg_plot_gene %in% rownames(obj[[assay]]), "Gene not found in the selected assay"))
    DefaultAssay(obj) <- assay
    res <- svg_results()
    moran <- if (!is.null(res) && input$svg_plot_gene %in% res$gene) sprintf(" (Moran's I = %.3f)", res$morans_I[res$gene == input$svg_plot_gene][1]) else ""
    SpatialFeaturePlot(obj, features = input$svg_plot_gene, pt.size.factor = input$svg_pt_size %||% 1.6) +
      ggtitle(paste0(input$svg_plot_gene, moran)) + theme(plot.title = element_text(hjust = 0.5, size = 14))
  })
  
  output$svg_plot <- renderPlot({ svg_plot_reactive() }, height = 600)
  
  output$download_svg_table <- downloadHandler(
    filename = function() paste0("spatially_variable_genes_", Sys.Date(), ".csv"),
    content = function(file) {
      res <- svg_results()
      write.csv(if (is.null(res)) data.frame(Message = "No results available") else res, file, row.names = FALSE)
    }
  )
  
  output$download_svg_plot <- downloadHandler(
    filename = function() paste0("spatial_", input$svg_plot_gene %||% "gene", "_", Sys.Date(), ".", input$svg_plot_format %||% "png"),
    content = function(file) {
      p <- svg_plot_reactive()
      ggsave(file, plot = p, width = 8, height = 7, dpi = 300, bg = "white", device = input$svg_plot_format %||% "png")
    }
  )
  
  ############################## Marker analysis ##############################
  
  # Cluster choices for marker analysis (annotated names if available)
  observe(run_safely("Interface update", persistent = FALSE, {
    req(spatial_obj())
    if (processing_states$clustered && "seurat_clusters" %in% colnames(spatial_obj()@meta.data)) {
      obj <- spatial_obj()
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
  }))
  
  # Comparison choices exclude the target cluster
  observeEvent(input$spatial_marker_cluster, run_safely("spatial marker cluster", {
    req(spatial_obj(), input$spatial_marker_cluster)
    obj <- spatial_obj()
    if ("annotated_clusters" %in% colnames(obj@meta.data)) {
      all_clusters <- levels(obj$annotated_clusters)
    } else {
      all_clusters <- levels(obj$seurat_clusters)
    }
    available_clusters <- setdiff(all_clusters, input$spatial_marker_cluster)
    updateSelectInput(session, "spatial_comparison_clusters",
                      choices = available_clusters,
                      selected = available_clusters[1])
  }))
  
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
  
  # Find markers for the selected cluster — p-values stored as raw numerics, formatted at render time
  observeEvent(input$run_spatial_markers, run_safely("run spatial markers", {
    message("=== FIND SPATIAL MARKERS STARTED ===")
    req(spatial_obj())
    req(input$spatial_marker_cluster)
    # Comparison clusters are required for one-vs-selected
    if (input$spatial_comparison_type == "one_vs_selected" &&
        (is.null(input$spatial_comparison_clusters) || length(input$spatial_comparison_clusters) == 0)) {
      showNotification("Please select at least one cluster to compare against", type = "warning")
      return()
    }
    tryCatch({
      obj <- spatial_obj()
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
      modal_message <- if (comparison_type == "one_vs_all") {
        paste("Analyzing cluster", selected_cluster, "vs all other clusters...")
      } else {
        paste("Analyzing cluster", selected_cluster, "vs clusters:",
              paste(input$spatial_comparison_clusters, collapse = ", "))
      }
      showModal(modalDialog(
        title = "Finding Marker Genes",
        div(modal_message, br(), tags$small("This may take a few minutes")),
        easyClose = FALSE,
        footer = NULL
      ))
      # Identities = active cluster column
      Idents(obj) <- cluster_col
      if (comparison_type == "one_vs_all") {
        markers <- FindMarkers(obj, ident.1 = selected_cluster, logfc.threshold = input$spatial_marker_logfc,
                               min.pct = input$spatial_marker_min_pct, only.pos = input$spatial_marker_only_pos, verbose = FALSE)
      } else {
        markers <- FindMarkers(obj, ident.1 = selected_cluster, ident.2 = input$spatial_comparison_clusters,
                               logfc.threshold = input$spatial_marker_logfc, min.pct = input$spatial_marker_min_pct,
                               only.pos = input$spatial_marker_only_pos, verbose = FALSE)
      }
      markers$gene <- rownames(markers)
      markers <- markers[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
      # Sort by fold change
      markers <- markers[order(-markers$avg_log2FC), ]
      # Round effect-size columns; p-values stay numeric (slider filter + correct numeric sort)
      markers$avg_log2FC <- round(markers$avg_log2FC, 3)
      markers$pct.1 <- round(markers$pct.1, 3)
      markers$pct.2 <- round(markers$pct.2, 3)
      spatial_marker_results_reactive(markers)
      message(paste("Found", nrow(markers), "marker genes"))
      removeModal()
      message("=== FIND SPATIAL MARKERS COMPLETED ===")
    }, error = function(e) {
      message(paste("=== FIND SPATIAL MARKERS ERROR ===", e$message))
      removeModal()
      showNotification(paste("Error finding markers:", conditionMessage(e)[1]), type = "error")
    })
  }))
  
  # Marker results — p-values kept numeric so the column filters render range sliders
  output$spatial_marker_results <- renderDT({
    markers <- spatial_marker_results_reactive()
    if (is.null(markers)) {
      return(datatable(data.frame(Message = "Click 'Find Markers' to analyze a cluster"),
                       options = list(dom = 't'),
                       rownames = FALSE))
    }
    datatable(
      markers,
      filter = "top",
      options = list(pageLength = 25, scrollX = TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel'),
                     columnDefs = list(list(className = 'dt-center', targets = 1:5))),
      rownames = FALSE,
      caption = paste("Marker genes for cluster", input$spatial_marker_cluster)
    ) %>%
      formatStyle('avg_log2FC',
                  backgroundColor = styleInterval(c(0.5, 1), c('white', 'lightblue', 'darkblue')),
                  color = styleInterval(c(1), c('black', 'white'))) %>%
      formatSignif(c('p_val', 'p_val_adj'), digits = 3)
  })
  
  # Marker results download (raw numeric values)
  output$download_spatial_markers <- downloadHandler(
    filename = function() {
      cluster_name <- gsub(" ", "_", input$spatial_marker_cluster)
      paste0("markers_cluster_", cluster_name, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      markers <- spatial_marker_results_reactive()
      if (!is.null(markers)) {
        write.csv(markers, file, row.names = FALSE)
      } else {
        # Empty file if no data
        write.csv(data.frame(Message = "No results available"), file, row.names = FALSE)
      }
    }
  )
  
}
