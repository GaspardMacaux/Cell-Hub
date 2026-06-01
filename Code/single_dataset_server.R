############################## Script Single Dataset Analysis ##############################

# single_dataset_server.R

single_dataset_server <- function(input, output, session) {
  ############################## Loading Data ##############################
  # Tab 1 : Loading Data
  
  shinyjs::useShinyjs() # button deactivation
  
  # Reactives objects for the single dataset part
  single_dataset_object <- reactiveVal(NULL)
  gene_list <- reactiveValues(features = NULL)
  # Analysis parameter tracking for CSV/PDF export
  qc_stats_single <- reactiveVal(NULL)
  
  # ✅ Create log system for SINGLE module
  log_system <- createLogSystem()
  
  # ✅ Render logs output
  output$single_loading_logs <- renderPrint({
    cat(log_system$log_value())
  })
  
  # File loading with logs (raw data: ZIP, H5, H5AD)
  observeEvent(input$file, {
    log_system$clear_logs()
    log_system$add_log("=== DATA LOADING STARTED ===", "info")
    showNotification("Uploading and processing data...", type = "message")
    tryCatch({
      file_extension <- logFileInfo(
        file_name = input$file$name,
        file_path = input$file$datapath,
        log_function = log_system$add_log
      )
      # Validate file type (ZIP, H5/HDF5, H5AD)
      if (!file_extension %in% c("zip", "h5", "hdf5", "h5ad")) {
        log_system$add_log("Invalid file type! Please upload .zip, .h5/.hdf5, or .h5ad", "error")
        stop("Please upload a .zip, .h5/.hdf5, or .h5ad file.")
      }
      # Inspect H5 / H5AD before loading
      if (file_extension %in% c("h5", "hdf5")) {
        inspectH5File(input$file$datapath, log_system$add_log)
      } else if (file_extension == "h5ad") {
        inspectH5ADFile(input$file$datapath, log_system$add_log)
      }
      # Clean workspace
      log_system$add_log("Cleaning workspace...", "info")
      cleanWorkspace("single")
      # Load data
      log_system$add_log(paste("Loading", input$dataset_type, "data for", input$species_choice), "info")
      log_system$add_log("Reading file... (this may take a while for large files)", "info")
      seuratObj <- loadRawData(
        file_path = input$file$datapath,
        dataset_type = input$dataset_type,
        species = input$species_choice,
        log_function = log_system$add_log
      )
      # Get dataset name
      dataset_name <- getDatasetFileName(list(input$file), "SingleDataset")
      seuratObj@project.name <- dataset_name
      logSeuratObjectStats(seuratObj, log_system$add_log)
      # Memory optimization if requested
      if (isTRUE(input$optimize_memory_single == "slim")) {
        log_system$add_log("Applying memory optimization...", "info")
        seuratObj <- slimSeuratObject(seuratObj, log_function = log_system$add_log)
      }
      logMemoryUsage(log_system$add_log, seuratObj, "Loaded object")
      single_dataset_object(seuratObj)
      showNotification(
        paste0("Data loaded successfully! Cells: ", ncol(seuratObj),
               ", Genes: ", nrow(seuratObj)),
        type = "message"
      )
      updateUIElements()
    }, error = function(e) {
      log_system$add_log(paste("ERROR:", conditionMessage(e)[1]), "error")
      log_system$add_log("Loading failed. Check file format and try again.", "error")
      showNotification(paste("Error processing files:", conditionMessage(e)[1]), type = "error")
      message("Error during file processing: ", conditionMessage(e)[1])
    })
  })
  
  
  # Load processed Seurat object (RDS or H5AD)
  # Load processed Seurat object (RDS or H5AD)
  observeEvent(input$load_seurat_file, {
    log_system$clear_logs()
    log_system$add_log("=== LOADING PROCESSED OBJECT ===", "info")
    tryCatch({
      file_extension <- logFileInfo(
        file_name = input$load_seurat_file$name,
        file_path = input$load_seurat_file$datapath,
        log_function = log_system$add_log
      )
      if (!file_extension %in% c("rds", "h5ad")) {
        log_system$add_log("Invalid file type! Please upload .rds or .h5ad", "error")
        stop("Please upload a .rds or .h5ad file.")
      }
      loaded_seurat <- loadSeuratObject(
        rds_path = input$load_seurat_file$datapath,
        clean_before = TRUE,
        module_type = "single",
        file_format = file_extension,
        log_function = log_system$add_log
      )
      # Memory optimization if requested
      if (isTRUE(input$optimize_memory_single == "slim")) {
        log_system$add_log("Applying memory optimization...", "info")
        loaded_seurat <- slimSeuratObject(loaded_seurat, log_function = log_system$add_log)
      }
      logMemoryUsage(log_system$add_log, loaded_seurat, "Loaded object")
      single_dataset_object(loaded_seurat)
      
      if ("umap" %in% names(loaded_seurat@reductions)) {
        log_system$add_log("Generating initial UMAP plot...", "info")
        plot_result <- generateInitialPlot(
          seurat_obj = loaded_seurat,
          remove_axes = input$remove_axes %||% FALSE,
          remove_legend = input$remove_legend %||% FALSE
        )
        if (!is.null(plot_result)) {
          single_dataset_object(plot_result$seurat_obj)
          clustering_plot(plot_result$plot)
          log_system$add_log("UMAP plot generated successfully", "success")
        }
      } else {
        log_system$add_log("No UMAP reduction found in object", "warning")
      }
      log_system$add_log("=== OBJECT LOADED SUCCESSFULLY ===", "success")
      showNotification("Object loaded successfully!", type = "message")
    }, error = function(e) {
      log_system$add_log(paste("ERROR:", conditionMessage(e)[1]), "error")
      log_system$add_log("Loading failed. Check file format and try again.", "error")
      showNotification(paste("Error loading object:", conditionMessage(e)[1]), type = "error")
    })
  })
  
  # Fonction réactive pour récupérer les clusters depuis l'objet Seurat
  get_clusters <- reactive({
    req(single_dataset_object())
    return(getClusters(single_dataset_object()))
  })



 

  # Function to update UI elements after data loading or dataset changes
  updateUIElements <- function() {
    req(single_dataset_object())  # Ensure the object is not null

    # Update cluster-related inputs
    cluster_choices <- unique(Idents(single_dataset_object()))
    updateSelectInput(session, "ident_1", choices = cluster_choices)
    updateCheckboxGroupInput(session, "ident_2", choices = cluster_choices, selected = "")
    updateSelectInput(session, "cluster_order_vln", choices = cluster_choices)
    updateSelectInput(session, "cluster_order_dot", choices = cluster_choices)


    # Update gene-related inputs
    gene_list <- sort(rownames(LayerData(single_dataset_object(), assay = "RNA", layer = 'counts')))
    updatePickerInput(session, "gene_select", choices = gene_list)
    updatePickerInput(session, "gene_select_genes_analysis", choices = gene_list)
  }

  # Observer for Seurat object changes
  observeEvent(single_dataset_object(), {
    updateUIElements()
  })




  ############################## QC metrics and normalization ##############################
  
  # Helper function for gene detection with different thresholds
  get_gene_count_with_threshold <- function(seurat_obj, min_cells = 3, min_counts = 0) {
    counts_matrix <- LayerData(seurat_obj, assay = "RNA", layer = "counts")
    genes_detected <- rowSums(counts_matrix >= min_counts) >= min_cells
    gene_count <- sum(genes_detected)
    
    message(sprintf("Genes detected (≥%d cells, ≥%d counts): %d", min_cells, min_counts, gene_count))
    return(gene_count)
  }
  
  # Count nuclei before and after QC (UNCHANGED)
  nuclei_count <- reactive({
    req(single_dataset_object())
    seuratRNA_subset <- subset(single_dataset_object(),
                               subset = nFeature_RNA > input$nFeature_range[1] &
                                 nFeature_RNA < input$nFeature_range[2] &
                                 percent.mt < input$percent.mt_max
    )
    if (exists("seuratRNA_subset")) {
      # Log QC parameters for verification
      message(sprintf("QC Parameters used:"))
      message(sprintf("- nFeature range: %d to %d", input$nFeature_range[1], input$nFeature_range[2]))
      message(sprintf("- Max MT%%: %f", input$percent.mt_max))
      message(sprintf("Original nuclei: %d", dim(single_dataset_object())[2]))
      message(sprintf("Filtered nuclei: %d", dim(seuratRNA_subset)[2]))
      
      # Verify filters were actually applied
      if(any(seuratRNA_subset$nFeature_RNA < input$nFeature_range[1]) ||
         any(seuratRNA_subset$nFeature_RNA > input$nFeature_range[2]) ||
         any(seuratRNA_subset$percent.mt > input$percent.mt_max)) {
        warning("Some nuclei outside QC parameters remain!")
      }
      
      return(dim(seuratRNA_subset)[2])
    } else {
      warning("Subsetting failed!")
      return(0)
    }
  })
  
  # DEBUGGING: Test all possible gene counting methods
  debug_gene_counts <- function(seurat_obj) {
    counts_matrix <- LayerData(seurat_obj, assay = "RNA", layer = "counts")
    
    message("=== DEBUGGING GENE COUNTS ===")
    
    # Method 1: All genes in matrix (no filtering)
    all_genes <- nrow(counts_matrix)
    message(sprintf("1. ALL genes in matrix: %d", all_genes))
    
    # Method 2: Genes detected in ≥1 cell
    genes_1_cell <- sum(rowSums(counts_matrix > 0) >= 1)
    message(sprintf("2. Genes in ≥1 cell: %d", genes_1_cell))
    
    # Method 3: Genes detected in ≥3 cells (your current)
    genes_3_cells <- sum(rowSums(counts_matrix > 0) >= 3)
    message(sprintf("3. Genes in ≥3 cells: %d", genes_3_cells))
    
    # Method 4: Genes with any expression > 0 anywhere
    genes_any_expression <- sum(rowSums(counts_matrix) > 0)
    message(sprintf("4. Genes with any UMI > 0: %d", genes_any_expression))
    
    # Method 5: Check for duplicated gene names
    unique_gene_names <- length(unique(rownames(counts_matrix)))
    message(sprintf("5. Unique gene names: %d", unique_gene_names))
    
    # Method 6: Check for empty gene names
    non_empty_names <- sum(rownames(counts_matrix) != "" & !is.na(rownames(counts_matrix)))
    message(sprintf("6. Non-empty gene names: %d", non_empty_names))
    
    target_count <- 16188
    message(sprintf("\nTARGET: %d", target_count))
    message(sprintf("Closest method: %s", 
                    names(which.min(abs(c(all_genes, genes_1_cell, genes_3_cells, genes_any_expression) - target_count)))))
    
    return(list(
      all = all_genes,
      one_cell = genes_1_cell, 
      three_cells = genes_3_cells,
      any_expression = genes_any_expression,
      unique_names = unique_gene_names,
      non_empty_names = non_empty_names
    ))
  }
  
  # Use this function to test:
  detected_genes_count <- reactive({
    req(single_dataset_object())
    
    # Run debugging first
    debug_results <- debug_gene_counts(single_dataset_object())
    
    # Based on debugging, choose the method closest to 16,188
    # CHANGE THIS LINE based on your debug results:
    final_count <- debug_results$one_cell  # Try this first
    
    return(final_count)
  })
  
  # GENES PER NUCLEUS (median/mean - AFTER QC filtering)
  nucleus_gene_metrics <- reactive({
    req(single_dataset_object())
    
    # Apply QC filters for per-nucleus metrics
    seuratRNA_subset <- subset(single_dataset_object(),
                               subset = nFeature_RNA > input$nFeature_range[1] &
                                 nFeature_RNA < input$nFeature_range[2] &
                                 percent.mt < input$percent.mt_max
    )
    
    if (exists("seuratRNA_subset")) {
      # Calculate both mean AND median 
      mean_genes <- mean(seuratRNA_subset$nFeature_RNA, na.rm = TRUE)
      median_genes <- median(seuratRNA_subset$nFeature_RNA, na.rm = TRUE)
      
      # Also get original values for comparison
      original_mean <- mean(single_dataset_object()$nFeature_RNA, na.rm = TRUE)
      original_median <- median(single_dataset_object()$nFeature_RNA, na.rm = TRUE)
      
      message(sprintf("=== GENES PER NUCLEUS ==="))
      message(sprintf("Original median: %.1f, filtered median: %.1f", original_median, median_genes))
      message(sprintf("Original mean: %.1f, filtered mean: %.1f", original_mean, mean_genes))
      
      return(list(
        mean = round(mean_genes, 0),
        median = round(median_genes, 0),
        original_median = round(original_median, 0)
      ))
    } else {
      return(list(mean = 0, median = 0, original_median = 0))
    }
  })
  
  # Render custom stat cards
  output$stats_cards_output <- renderUI({
    req(single_dataset_object())
    
    tryCatch({
      # Get values from the same reactive sources
      count_value <- nuclei_count()
      total_nuclei <- dim(single_dataset_object())[2]
      retention_pct <- round((count_value/total_nuclei)*100, 2)
      
      gene_count <- detected_genes_count()
      
      gene_metrics <- nucleus_gene_metrics()
      median_genes <- gene_metrics$median
      original_median <- gene_metrics$original_median
      
      # Create cards
      tagList(
        # Card 1: Number of Nuclei
        tags$div(
          class = "stat-card stat-card-violet",
          tags$div(class = "stat-icon", icon("dna", class = "fa-3x")),
          tags$div(
            class = "stat-content",
            tags$div(class = "stat-title", "Number of Nuclei"),
            tags$div(class = "stat-value", format(count_value, big.mark = ",")),
            tags$div(class = "stat-subtitle", 
                     paste0(retention_pct, "% retained"))
          )
        ),
        
        # Card 2: Unique Genes
        tags$div(
          class = "stat-card stat-card-orange",
          tags$div(class = "stat-icon", icon("dna", class = "fa-3x")),
          tags$div(
            class = "stat-content",
            tags$div(class = "stat-title", "Unique Genes"),
            tags$div(class = "stat-value", format(gene_count, big.mark = ",")),
            tags$div(class = "stat-subtitle", "detected in dataset")
          )
        ),
        
        # Card 3: Median Genes per Nucleus
        tags$div(
          class = "stat-card stat-card-green",
          tags$div(class = "stat-icon", icon("chart-line", class = "fa-3x")),
          tags$div(
            class = "stat-content",
            tags$div(class = "stat-title", "Median Genes/Nucleus"),
            tags$div(class = "stat-value", format(median_genes, big.mark = ",")),
            tags$div(class = "stat-subtitle", 
                     paste0("was ", format(original_median, big.mark = ","), " before QC"))
          )
        )
      )
      
    }, error = function(e) {
      tags$div(
        style = "padding: 20px; background: #ffebee; border-radius: 8px; border-left: 4px solid #f44336;",
        icon("exclamation-triangle"), " Error loading statistics: ", e$message
      )
    })
  })
  
  # Display QC metrics - Manual plot with cell points (complete Seurat style)
  output$vlnplot <- renderPlot({
    req(input$QCmetrics, single_dataset_object())
    
    tryCatch({
      if (input$QCmetrics) {
        message("=== Creating QC VlnPlot (Seurat-style with points) ===")
        
        seurat_obj <- isolate(single_dataset_object())
        
        # Verify metrics exist
        required_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
        missing_metrics <- setdiff(required_metrics, colnames(seurat_obj@meta.data))
        
        if(length(missing_metrics) > 0) {
          stop(paste("Missing metrics:", paste(missing_metrics, collapse = ", ")))
        }
        
        # Extract metadata
        meta_data <- as.data.frame(seurat_obj@meta.data)
        
        # Add cluster info
        meta_data$cluster <- tryCatch({
          as.character(Idents(seurat_obj))
        }, error = function(e) {
          rep("All", nrow(meta_data))
        })
        
        # Order clusters
        cluster_levels <- unique(meta_data$cluster)
        numeric_clusters <- cluster_levels[!is.na(suppressWarnings(as.numeric(cluster_levels)))]
        text_clusters <- cluster_levels[is.na(suppressWarnings(as.numeric(cluster_levels)))]
        
        if (length(numeric_clusters) > 0) {
          numeric_clusters <- as.character(sort(as.numeric(numeric_clusters)))
        }
        if (length(text_clusters) > 0) {
          text_clusters <- sort(text_clusters)
        }
        
        meta_data$cluster <- factor(meta_data$cluster, levels = c(numeric_clusters, text_clusters))
        
        # Reshape for plotting
        plot_data <- tidyr::pivot_longer(
          meta_data,
          cols = all_of(required_metrics),
          names_to = "metric",
          values_to = "value"
        )
        
        # Keep original metric names
        plot_data$metric <- factor(plot_data$metric, levels = required_metrics)
        
        # Seurat default colors
        n_clusters <- length(unique(meta_data$cluster))
        colors <- scales::hue_pal()(n_clusters)
        
        # Create Seurat-style violin plot WITH POINTS
        plot <- ggplot(plot_data, aes(x = cluster, y = value, fill = cluster)) +
          geom_violin(scale = "width", trim = TRUE) +
          geom_jitter(height = 0, width = 0.3, size = 0.1, alpha = 0.4) +  # ← ADD POINTS
          facet_wrap(~metric, scales = "free_y", ncol = 3) +
          theme_classic() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none",
            strip.background = element_blank(),
            strip.text = element_text(size = 11, face = "bold")
          ) +
          labs(x = NULL, y = NULL) +
          scale_fill_manual(values = colors)
        
        message("✓ VlnPlot created successfully")
        return(plot)
      }
    }, error = function(e) {
      message("ERROR: ", e$message)
      showNotification(paste0("Error: ", e$message), type = "error")
      
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, 
                 label = paste("Error:\n", e$message), 
                 color = "red", size = 5) +
        theme_void()
    })
  })
  
  # Display scatter plots with verification
  output$scatter_plot1 <- renderPlot({
    req(input$show_plots, single_dataset_object())
    tryCatch({
      if (input$show_plots) {
        # Create subset with verification
        seuratRNA_subset <- subset(single_dataset_object(),
                                   subset = nFeature_RNA > input$nFeature_range[1] &
                                     nFeature_RNA < input$nFeature_range[2] &
                                     percent.mt < input$percent.mt_max
        )

        # Verify subsetting worked
        if(ncol(seuratRNA_subset) == 0) {
          stop("No cells pass the current filters!")
        }

        FeatureScatter(seuratRNA_subset,
                       feature1 = "nCount_RNA",
                       feature2 = "percent.mt")
      }
    }, error = function(e) {
      showNotification(paste0("Error in scatter plot 1: ", e$message), type = "error")
    })
  })

  # Display scatter plot 2 with same verifications
  output$scatter_plot2 <- renderPlot({
    req(input$show_plots, single_dataset_object())
    tryCatch({
      if (input$show_plots) {
        seuratRNA_subset <- subset(single_dataset_object(),
                                   subset = nFeature_RNA > input$nFeature_range[1] &
                                     nFeature_RNA < input$nFeature_range[2] &
                                     percent.mt < input$percent.mt_max
        )

        if(ncol(seuratRNA_subset) == 0) {
          stop("No cells pass the current filters!")
        }

        FeatureScatter(seuratRNA_subset,
                       feature1 = "nCount_RNA",
                       feature2 = "nFeature_RNA")
      }
    }, error = function(e) {
      showNotification(paste0("Error in scatter plot 2: ", e$message), type = "error")
    })
  })

  # Apply QC filters with thorough verification
  # Apply QC filters with thorough verification
  observeEvent(input$apply_qc, {
    req(single_dataset_object())
    
    tryCatch({
      showModal(modalDialog(
        title = "Please Wait",
        "Applying QC filters...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      seurat_obj <- single_dataset_object()
      
      # Capture pre-QC stats BEFORE any filtering
      pre_n      <- ncol(seurat_obj)
      pre_genes  <- round(median(seurat_obj$nFeature_RNA, na.rm = TRUE))
      pre_counts <- round(median(seurat_obj$nCount_RNA,   na.rm = TRUE))
      pre_mt     <- round(median(seurat_obj$percent.mt,   na.rm = TRUE), 2)
      
      # Record mean metrics for notification (unchanged)
      orig_mt_mean      <- mean(seurat_obj$percent.mt)
      orig_feature_mean <- mean(seurat_obj$nFeature_RNA)
      
      # Apply filters with verification
      seurat_obj <- subset(seurat_obj,
                           subset = nFeature_RNA > input$nFeature_range[1] &
                             nFeature_RNA < input$nFeature_range[2] &
                             percent.mt < input$percent.mt_max)
      
      # Verify filtering
      if (any(seurat_obj$nFeature_RNA <= input$nFeature_range[1]) ||
          any(seurat_obj$nFeature_RNA >= input$nFeature_range[2]) ||
          any(seurat_obj$percent.mt   >= input$percent.mt_max)) {
        stop("QC filtering failed - some nuclei outside parameters remain")
      }
      
      single_dataset_object(seurat_obj)
      
      # Store pre/post QC stats for export
      qc_stats_single(list(
        n_before             = pre_n,
        n_after              = ncol(seurat_obj),
        median_genes_before  = pre_genes,
        median_genes_after   = round(median(seurat_obj$nFeature_RNA, na.rm = TRUE)),
        median_counts_before = pre_counts,
        median_counts_after  = round(median(seurat_obj$nCount_RNA,   na.rm = TRUE)),
        median_mt_before     = pre_mt,
        median_mt_after      = round(median(seurat_obj$percent.mt,   na.rm = TRUE), 2)
      ))
      
      num_nuclei_after_qc <- ncol(seurat_obj)
      msg <- paste0(
        "QC filters applied successfully:\n",
        sprintf("- Nuclei: %d → %d (%d%% retained)\n",
                pre_n, num_nuclei_after_qc,
                round(num_nuclei_after_qc / pre_n * 100)),
        sprintf("- Mean MT%%: %.2f%% → %.2f%%\n",
                orig_mt_mean, mean(seurat_obj$percent.mt)),
        sprintf("- Mean Features: %.1f → %.1f",
                orig_feature_mean, mean(seurat_obj$nFeature_RNA))
      )
      
      showNotification(msg, type = "message", duration = 10)
      removeModal()
      updateUIElements()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Error applying QC: ", e$message), type = "error")
    })
  })

  observeEvent(input$normalize_data, {
    req(input$num_var_features, input$normalization_method_single)
    tryCatch({
      # Show modal dialog
      showModal(modalDialog(
        title = "Please Wait",
        paste0("Normalizing data using ", input$normalization_method_single, " method..."),
        easyClose = FALSE,
        footer = NULL
      ))
      
      if (input$normalization_method_single == "SCTransform") {
        # SCTransform replaces NormalizeData, ScaleData, and FindVariableFeatures
        normalized_seurat <- SCTransform(
          single_dataset_object(), 
          variable.features.n = input$num_var_features,
          verbose = FALSE
        )
        
        # IMPORTANT: Set default assay to SCT for downstream analysis
        DefaultAssay(normalized_seurat) <- "SCT"
        
        print("SCTransform normalization completed")
        print(paste0("Default assay set to: ", DefaultAssay(normalized_seurat)))
        print(paste0("Variable features found: ", length(VariableFeatures(normalized_seurat))))
        
      } else {
        # Standard normalization workflow
        normalized_seurat <- NormalizeData(
          single_dataset_object(), 
          normalization.method = input$normalization_method_single, 
          scale.factor = input$scale_factor
        )
        print(paste0("Data normalization completed using method: ", input$normalization_method_single))
        
        # Identification of variable features
        normalized_seurat <- FindVariableFeatures(
          normalized_seurat, 
          selection.method = "vst", 
          nfeatures = input$num_var_features
        )
        print("Variable features identification completed")
        
        # Keep default assay as RNA
        DefaultAssay(normalized_seurat) <- "RNA"
      }
      
      single_dataset_object(normalized_seurat)
      
      # Close modal dialog
      removeModal()
      
      showNotification(
        paste0("Normalization completed using ", input$normalization_method_single, " method!"),
        type = "message"
      )
      
    }, error = function(e) {
      removeModal()
      showNotification(
        paste0("Error during data normalization: ", e$message), 
        type = "error",
        duration = 10
      )
      print(paste("Normalization error:", e$message))
    })
  })
  
  # Plot rendering of variable features
  output$variable_feature_plot <- renderPlot({
    req(input$show_plots, single_dataset_object())
    
    # Vérifications de sécurité
    obj <- single_dataset_object()
    if (is.null(obj)) return(NULL)
    
    # Vérifier que les variable features existent
    var_features <- VariableFeatures(obj)
    if (is.null(var_features) || length(var_features) == 0) {
      return(NULL)
    }
    
    tryCatch({
      plot1 <- VariableFeaturePlot(obj)
      top10 <- head(var_features, 10)
      
      if (length(top10) > 0) {
        plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
        return(plot2)
      } else {
        return(plot1)
      }
    }, error = function(e) {
      message("Error creating variable feature plot: ", e$message)
      return(NULL)
    })
  })

  ############################## Scaling, PCA and elbow plot ##############################
  #  Tab 3: Scaling and PCA and elbow plot

  # Scaling data
  observeEvent(input$scale_button, {
    showNotification("Scaling data...", type = "message")
    tryCatch({
      showModal(modalDialog(
        title = "Please Wait",
        "Running PCA on the dataset...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      seurat_obj <- single_dataset_object()
      current_assay <- DefaultAssay(seurat_obj)
      
      # Check if SCTransform was used
      if (current_assay == "SCT" || "SCT" %in% names(seurat_obj@assays)) {
        print("SCTransform detected - data already scaled, skipping ScaleData")
        # SCTransform already scaled the data, just run PCA
        # Make sure we're using the SCT assay
        DefaultAssay(seurat_obj) <- "SCT"
        
      } else {
        # Standard workflow: scale RNA data
        print("Scaling data for all genes in RNA assay")
        seurat_obj <- ScaleData(
          seurat_obj,
          assay = "RNA", 
          features = rownames(seurat_obj)
        )
        print("Data scaling completed")
      }
      print(paste0("Running PCA on assay: ", DefaultAssay(seurat_obj)))
      
      pca_result <- RunPCA(
        seurat_obj, 
        features = VariableFeatures(seurat_obj),
        verbose = FALSE
      )
      
      # Memory optimization: drop scale.data after PCA (no longer needed, recomputed for heatmaps)
      if (isTRUE(input$optimize_memory_single == "slim")) {
        log_system$add_log("Applying post-PCA memory optimization...", "info")
        pca_result <- slimSeuratObject(pca_result, keep_counts = TRUE,
                                       keep_scale_data = FALSE,
                                       log_function = log_system$add_log)
      }
      single_dataset_object(pca_result)
      logMemoryUsage(log_system$add_log, pca_result, "Post-PCA object")
      print("PCA completed successfully")
      
      shinyjs::enable("pca_button")
      removeModal()
      showNotification("PCA completed successfully.", type = "message")
      
    }, error = function(e) {
      removeModal()
      print(paste("Error during scaling/PCA:", e$message))
      showNotification(paste("Error during scaling/PCA:", e$message), type = "error", duration = 10)
    })
  })


  #Printing PCA results
  output$pca_results <- renderPrint({
    req(!is.null(single_dataset_object()))
    if ("pca" %in% names(single_dataset_object())) {
      print(single_dataset_object()$pca, dims = 1:5, nfeatures = 5)
    }
  })

  #Vizdimloading plot
  output$loading_plot <- renderPlot({
    req(!is.null(single_dataset_object()))
    tryCatch(VizDimLoadings(single_dataset_object(), dims = 1:2, reduction = "pca"), error = function(e){})
  })

  # Dimplot
  output$dim_plot <- renderPlot({
    req(!is.null(single_dataset_object()))
    tryCatch(DimPlot(single_dataset_object(), reduction = "pca"), error = function(e){})
  })

  # Show the ElbowPlot
  observeEvent(input$run_elbow, {
    tryCatch({
      output$elbow_plot <- renderPlot({
        req(single_dataset_object())
        ElbowPlot(single_dataset_object(), ndims = 50)
      })}, error = function(e) {
        showNotification(paste("Error running elbow plot:", e$message), type = "error")
      })
  })

  ############################## Neighbors calculation and clustering ##############################

  clustering_plot <- reactiveVal()

  # Finding neighbors
  observeEvent(input$run_neighbors, {
    tryCatch({
      req(!is.null(single_dataset_object()))
      showModal(modalDialog(
        title = "Please Wait",
        "Looking for neighbors and calculating UMAP...",
        easyClose = FALSE,
        footer = NULL
      ))
      seurat_tmp <- single_dataset_object()
      seurat_tmp <- findNeighbors_reproducible(
        object = seurat_tmp,
        dims = 1:input$dimension_1,
        seed = 42,
        verbose = FALSE
      )
      seurat_tmp <- runUMAP_reproducible(
        object = seurat_tmp,
        dims = 1:input$dimension_1,
        seed = 42,
        verbose = FALSE
      )
      single_dataset_object(seurat_tmp)
      print("Neighbors found and UMAP calculated.")
      plot <- DimPlot(single_dataset_object(), group.by = "ident") + ggtitle(NULL)
      if(input$remove_axes) { plot <- plot + NoAxes() }
      if(input$remove_legend) { plot <- plot + NoLegend() }
      clustering_plot(plot)
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Neighbor search error: ", e$message), type = "error")
    })
  })

# Clustering
  observeEvent(input$run_clustering, {
    tryCatch({
      req(!is.null(single_dataset_object()))
      showModal(modalDialog(
        title = "Please Wait",
        "Clustering process started...",
        easyClose = FALSE,
        footer = NULL
      ))
      seurat_tmp <- findClusters_reproducible(
        object     = single_dataset_object(),
        resolution = input$resolution_step1,
        algorithm  = as.integer(input$algorithm_select),
        seed       = 42,
        verbose    = TRUE
      )
      single_dataset_object(seurat_tmp)
      
      plot <- DimPlot(single_dataset_object(), group.by = "ident", label = FALSE) + ggtitle("")
      if (input$remove_axes)   { plot <- plot + NoAxes() }
      if (input$remove_legend) { plot <- plot + NoLegend() }
      clustering_plot(plot)
      
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Clustering error: ", e$message), type = "error")
    })
  })



  output$clustering_plot <- renderPlot({
    req(clustering_plot())
    clustering_plot()
  })

  #SAVING UMAP PLOT
  output$downloadUMAP <- createDownloadHandler(
    reactive_data = clustering_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "UMAP_plot")
    }),
    data_name = "UMAP",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$umap_export_format }),
      width = reactive({ ifelse(input$umap_export_format == "pdf", 11, 10) }),
      height = reactive({ ifelse(input$umap_export_format == "pdf", 8, 6) }),
      dpi = reactive({ input$dpi_umap })
    )
  )
  
######################Doublet Finder#########################
  
  # Reactive to store doublet results
  doublet_results_data <- reactiveVal(NULL)
  
  observeEvent(input$run_doubletfinder, {
    req(single_dataset_object())
    showModal(modalDialog(
      title = "Running DoubletFinder",
      p("Processing your data..."),
      footer = NULL,
      easyClose = FALSE
    ))
    tryCatch({
      seu <- single_dataset_object()
      meta <- seu@meta.data
      seu@meta.data <- data.frame(
        lapply(meta, function(x) {
          if (is.data.frame(x)) x[[1]] else as.vector(x)
        }),
        row.names = rownames(meta),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      sweep.res.list <- paramSweep(seu, PCs = 1:input$pc_use, sct = FALSE)
      sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
      bcmvn <- find.pK(sweep.stats)
      nExp_poi <- round(input$doublet_rate / 100 * ncol(seu))
      pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
      message("[DoubletFinder] nExp_poi: ", nExp_poi, " | optimal pK: ", pK)
      # reuse.pANN must be NULL (not FALSE) to trigger fresh pANN computation —
      # doubletFinder checks !is.null(reuse.pANN) so FALSE enters the wrong branch
      # and crashes xtfrm() trying to order seu@meta.data[, FALSE]
      seu <- doubletFinder(seu,
                           PCs = 1:input$pc_use,
                           pN = input$pN_value,
                           pK = pK,
                           nExp = nExp_poi,
                           reuse.pANN = NULL,
                           sct = FALSE)
      df_col <- grep("DF.classifications", colnames(seu@meta.data), value = TRUE)[1]
      message("[DoubletFinder] df_col: ", df_col)
      message("[DoubletFinder] classifications: ", paste(unique(seu@meta.data[[df_col]]), collapse = ", "))
      output$doublet_umap <- renderPlot({
        DimPlot(seu,
                reduction = "umap",
                group.by = df_col,
                cols = c("Singlet" = "grey", "Doublet" = "red")) +
          ggtitle("Doublets Identification")
      })
      stats <- table(seu@meta.data[[df_col]])
      stats_df <- data.frame(
        Category = names(stats),
        Count = as.numeric(stats),
        Percentage = round(as.numeric(stats) / sum(stats) * 100, 2)
      )
      output$doublet_stats <- renderDT({ stats_df })
      doublet_results_data(stats_df)
      single_dataset_object(seu)
      removeModal()
      showNotification("DoubletFinder analysis completed!", type = "message")
    }, error = function(e) {
      message("[DoubletFinder] ERROR: ", e$message)
      removeModal()
      showNotification(paste("Error in DoubletFinder:", e$message), type = "error")
    })
  })
  
  # Download handler for doublet results
  output$download_doublet_results <- downloadHandler(
    filename = function() {
      obj_name <- getObjectNameForDownload(single_dataset_object(), default_name = "Dataset")
      paste0(obj_name, "_doublet_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(doublet_results_data())
      write.csv(doublet_results_data(), file, row.names = FALSE)
    }
  )
  
  # Remove doublets
  observeEvent(input$remove_doublets, {
    req(single_dataset_object())
    tryCatch({
      seu <- single_dataset_object()
      df_col <- grep("DF.classifications", colnames(seu@meta.data), value = TRUE)[1]
      if (is.null(df_col)) {
        showNotification("Please run doublet detection first", type = "warning")
        return()
      }
      # Use colnames(seu) — guaranteed alignment with Seurat v5 subset()
      cells_keep <- colnames(seu)[seu@meta.data[[df_col]] == "Singlet"]
      message("[RemoveDoublets] singlets to keep: ", length(cells_keep))
      if (length(cells_keep) == 0) {
        showNotification("No singlets found — check DoubletFinder classification column", type = "error")
        return()
      }
      seu_filtered <- subset(seu, cells = cells_keep)
      single_dataset_object(seu_filtered)
      output$doublet_umap <- renderPlot({
        DimPlot(seu_filtered, reduction = "umap", group.by = "seurat_clusters") +
          ggtitle("UMAP (Doublets Removed)")
      })
      showNotification(
        sprintf("Removed %d doublets from dataset", ncol(seu) - ncol(seu_filtered)),
        type = "message"
      )
    }, error = function(e) {
      message("[RemoveDoublets] ERROR: ", e$message)
      showNotification(paste("Error removing doublets:", e$message), type = "error")
    })
  })
  
  ############################## Visualization of expressed genes ##############################


  # Observer pour mettre à jour les choix d'assays
  observe({
    req(single_dataset_object())
    updateAssayChoices(session, single_dataset_object())
  })
  
  # Observer pour mettre à jour les choix de gènes
  observeEvent(c(single_dataset_object(), input$viz_assay), {
    req(single_dataset_object(), input$viz_assay)
    updateGeneChoices(session, single_dataset_object(), input$viz_assay)
  })
  
  # Observer pour mettre à jour les sélecteurs de cluster
  observe({
    req(single_dataset_object())
    updateClusterChoices(session, single_dataset_object())
  })
  
  
  # Observer to update text inputs when genes are selected in pickerInput
  observeEvent(input$gene_select, {
    selected_genes <- input$gene_select
    
    if (!is.null(selected_genes) && length(selected_genes) > 0) {
      # Convert selected genes to comma-separated string
      genes_text <- paste(selected_genes, collapse = ", ")
      
      # Update all gene text inputs with selected genes
      updateTextInput(session, "gene_list_feature", value = genes_text)
      updateTextInput(session, "gene_list_dotplot", value = genes_text)
      updateTextInput(session, "gene_list_ridge_plot", value = genes_text)
      updateTextInput(session, "gene_list_vln", value = genes_text)  # If you have this one too
    }
  })
  
  ############################## FeaturePlot - SINGLE DATASET with 3D ##############################
  
  # FeaturePlot avec options 2D/3D
  feature_plot <- reactiveVal()
  
  observeEvent(input$show_feature, {
    tryCatch({
      req(input$gene_list_feature, single_dataset_object())
      genes <- unique(trimws(unlist(strsplit(input$gene_list_feature, ","))))
      seurat_object <- single_dataset_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay
      present_genes <- genes[genes %in% rownames(LayerData(seurat_object, assay = "RNA", layer = "counts"))]
      
      if (length(present_genes) > 0) {
        print("Creating plot")
        
        # Get cutoffs - CORRECTED VERSION
        # Cutoffs can be quantiles (e.g., "q10") or numeric values
        min_cut <- if (!is.null(input$min_cutoff) && !is.na(input$min_cutoff) && input$min_cutoff != "") {
          # If it starts with 'q', keep it as a string (quantile)
          # Otherwise convert to numeric
          if (grepl("^q", input$min_cutoff)) {
            input$min_cutoff  # Keep as string for Seurat to interpret as quantile
          } else {
            as.numeric(input$min_cutoff)
          }
        } else {
          "q0"  # Default: no minimum cutoff
        }
        
        max_cut <- if (!is.null(input$max_cutoff) && !is.na(input$max_cutoff) && input$max_cutoff != "") {
          # If it starts with 'q', keep it as a string (quantile)
          # Otherwise convert to numeric
          if (grepl("^q", input$max_cutoff)) {
            input$max_cutoff  # Keep as string for Seurat to interpret as quantile
          } else {
            as.numeric(input$max_cutoff)
          }
        } else {
          NA  # Default: no maximum cutoff (show full range)
        }
        
        message(paste("FeaturePlot cutoffs: min =", min_cut, ", max =", max_cut))
        
        # Check if 3D mode is enabled
        if(!is.null(input$enable_3d_feature) && input$enable_3d_feature) {
          # 3D MODE
          print("Creating 3D FeaturePlot")
          
          # Check if 3D reduction exists (using existing function name: umap3d)
          if(!"umap3d" %in% names(seurat_object@reductions)) {
            showNotification(
              "3D UMAP not found. Please run 'Run UMAP 3D' first in the dimensional reduction section.",
              type = "warning",
              duration = 5
            )
            return()
          }
          
          # Limit to 3 genes for 3D
          if(length(present_genes) > 3) {
            showNotification(
              "Maximum 3 genes allowed for 3D plot. Using first 3 genes.",
              type = "warning",
              duration = 3
            )
            present_genes <- present_genes[1:3]
          }
          
          # Get display options
          hide_grid <- !is.null(input$hide_grid_feature) && input$hide_grid_feature
          hide_axes <- !is.null(input$add_noaxes_feature) && input$add_noaxes_feature
          dark_mode <- !is.null(input$dark_mode_feature) && input$dark_mode_feature
          
          # Create 3D plot with display options
          plot <- create3DFeaturePlot(
            seurat_object = seurat_object,
            genes = present_genes,
            reduction = "umap3d",  # Using existing reduction name
            pt_size = 3,
            min_cutoff = min_cut,
            max_cutoff = max_cut,
            order = TRUE,
            hide_grid = hide_grid,
            hide_axes = hide_axes,
            dark_mode = dark_mode
          )
          
          feature_plot(plot)
          print("3D plot created successfully")
          
        } else {
          # 2D MODE (original)
          print("Creating 2D FeaturePlot")
          
          # Configuration de base du thème
          base_theme <- theme(
            axis.text = element_text(size = input$axis_text_size),
            axis.title = element_text(size = input$title_text_size),
            plot.title = element_text(size = input$title_text_size),
            legend.text = element_text(size = input$axis_text_size),
            axis.line = element_line(linewidth = input$axis_line_width),
            axis.ticks = element_line(linewidth = input$axis_line_width)
          )
          
          # Créer le plot avec les paramètres appropriés
          plot <- FeaturePlot(
            seurat_object,
            features = present_genes,
            blend = input$show_coexpression && length(present_genes) > 1,
            blend.threshold = 1,
            order = TRUE,
            min.cutoff = min_cut,
            max.cutoff = max_cut
          ) + base_theme
          
          # Ajout des modifications conditionnelles
          if (input$add_noaxes_feature) {
            print("Adding NoAxes")
            plot <- plot + NoAxes()
          }
          if (input$add_nolegend_feature) {
            print("Adding NoLegend")
            plot <- plot + NoLegend()
          }
          
          # Dark mode for 2D (if checkbox exists)
          if(!is.null(input$dark_mode_feature) && input$dark_mode_feature) {
            plot <- plot + 
              theme(
                panel.background = element_rect(fill = "#1a1a1a"),
                plot.background = element_rect(fill = "#1a1a1a"),
                panel.grid = element_blank(),
                text = element_text(color = "white"),
                axis.text = element_text(color = "white"),
                axis.title = element_text(color = "white"),
                legend.background = element_rect(fill = "#1a1a1a"),
                legend.text = element_text(color = "white"),
                legend.title = element_text(color = "white")
              )
          }
          
          feature_plot(plot)
          print("2D plot created successfully")
        }
      }
    }, error = function(e) {
      showNotification(paste("Error in FeaturePlot:", e$message), type = "error")
      print(paste("Error details:", e$message))
    })
  })
  
  # Render output - handle both 2D (ggplot) and 3D (plotly)
  output$feature_plot <- renderPlot({
    req(feature_plot())
    plot_obj <- feature_plot()
    
    # Check if it's a ggplot object (2D)
    if(inherits(plot_obj, "gg") || inherits(plot_obj, "patchwork")) {
      plot_obj
    }
  })
  
  # Add plotly output for 3D
  output$feature_plot_3d <- renderPlotly({
    req(feature_plot())
    plot_obj <- feature_plot()
    
    # Check if it's a plotly object (3D)
    if(inherits(plot_obj, "plotly")) {
      plot_obj
    }
  })
  
  # DotPlot reactive value
  dot_plot <- reactiveVal()
  # DotPlot observer - using modular function
  observeEvent(input$show_dot, {
    tryCatch({
      req(single_dataset_object())
      genes <- unique(trimws(strsplit(input$gene_list_dotplot, ",")[[1]]))
      seurat_object <- single_dataset_object()
      
      # Filter to selected clusters if specified
      seurat_subset <- seurat_object
      if (!is.null(input$cluster_order_dot) && length(input$cluster_order_dot) > 0) {
        Idents(seurat_subset) <- factor(
          Idents(seurat_subset),
          levels = input$cluster_order_dot
        )
        seurat_subset <- subset(seurat_subset, idents = input$cluster_order_dot)
        message(paste("Filtered to clusters:", paste(input$cluster_order_dot, collapse=", ")))
      }
      
      # Get parameters
      color_palette <- input$color_palette_dotplot_single %||% "default"
      dot_scale <- input$dot_scale_single %||% 1
      
      # Create DotPlot using modular function
      plot <- create_dotplot(
        seurat_object = seurat_subset,
        genes = genes,
        group_by = NULL,  # Use Idents
        split_by = NULL,
        color_palette = color_palette,
        dot_scale = dot_scale,
        assay = input$viz_assay,
        rotate_axis = TRUE,
        axis_text_size = input$axis_text_size,
        title_text_size = input$title_text_size,
        axis_line_width = input$axis_line_width
      )
      
      # Apply optional modifications
      if (input$add_noaxes_dot) plot <- plot + NoAxes()
      if (input$add_nolegend_dot) plot <- plot + NoLegend()
      if (input$invert_axes) plot <- plot + coord_flip()
      
      dot_plot(plot)
      
      message(paste("DotPlot created - palette:", color_palette, "| dot.scale:", dot_scale))
      
    }, error = function(e) {
      showNotification(paste("Error in DotPlot:", e$message), type = "error")
      message(paste("DotPlot error details:", e$message))
    })
  })
  
  # Display the DotPlot
  output$dot_plot <- renderPlot({
    req(dot_plot())
    dot_plot()
  })

  # SINGLE DATASET - VlnPlot and RidgePlot REPLACEMENTS
  # Replace the existing observeEvent blocks in single_dataset_server_newarch.R
  
  ############################## VlnPlot - SINGLE DATASET ##############################
  
  # VlnPlot
  vln_plot <- reactiveVal()
  observeEvent(input$show_vln, {
    tryCatch({
      req(input$gene_list_vln)
      print("Starting VlnPlot generation")
      
      genes <- unique(trimws(unlist(strsplit(input$gene_list_vln, ","))))
      print(paste("Processing genes:", paste(genes, collapse = ", ")))
      
      seurat_object <- single_dataset_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay
      
      if (!is.null(input$cluster_order_vln) && length(input$cluster_order_vln) > 0) {
        print(paste("Selected clusters:", paste(input$cluster_order_vln, collapse=", ")))
        
        seurat_subset <- seurat_object
        Idents(seurat_subset) <- factor(
          Idents(seurat_subset),
          levels = input$cluster_order_vln
        )
        seurat_subset <- subset(seurat_subset, idents = input$cluster_order_vln)
        
        if(is.null(seurat_subset) || ncol(seurat_subset) == 0) {
          showNotification("No cells found with selected clusters.", type = "error")
          return()
        }
        
        # Use manual VlnPlot function
        plot <- createManualVlnPlot(
          seurat_object = seurat_subset,
          genes = genes,
          group_by = "ident",
          pt_size = ifelse(input$hide_vln_points, 0, 1),
          axis_text_size = ifelse(is.null(input$axis_text_size), 12, input$axis_text_size),
          title_text_size = ifelse(is.null(input$title_text_size), 14, input$title_text_size),
          axis_line_width = ifelse(is.null(input$axis_line_width), 1, input$axis_line_width),
          add_noaxes = !is.null(input$add_noaxes_vln) && input$add_noaxes_vln,
          add_nolegend = !is.null(input$add_nolegend_vln) && input$add_nolegend_vln
        )
        
      } else {
        # No cluster selection
        plot <- createManualVlnPlot(
          seurat_object = seurat_object,
          genes = genes,
          group_by = "ident",
          pt_size = ifelse(input$hide_vln_points, 0, 1),
          axis_text_size = ifelse(is.null(input$axis_text_size), 12, input$axis_text_size),
          title_text_size = ifelse(is.null(input$title_text_size), 14, input$title_text_size),
          axis_line_width = ifelse(is.null(input$axis_line_width), 1, input$axis_line_width),
          add_noaxes = !is.null(input$add_noaxes_vln) && input$add_noaxes_vln,
          add_nolegend = !is.null(input$add_nolegend_vln) && input$add_nolegend_vln
        )
      }
      
      vln_plot(plot)
      print("Plot successfully stored")
      
    }, error = function(e) {
      showNotification(paste("Error in VlnPlot:", e$message), type = "error")
      print(paste("Error in VlnPlot:", e$message))
      print("Error details:")
      print(e)
    })
  })
  
  output$vln_plot <- renderPlot({
    req(vln_plot())
    vln_plot()
  })
  
  
  ############################## RidgePlot - SINGLE DATASET ##############################
  
  # RidgePlot
  ridge_plot <- reactiveVal()
  
  observeEvent(input$show_ridge, {
    tryCatch({
      req(input$gene_list_ridge_plot, single_dataset_object())
      
      genes <- unique(trimws(unlist(strsplit(input$gene_list_ridge_plot, ","))))
      print(paste("Genes for RidgePlot:", paste(genes, collapse = ", ")))
      
      seurat_object <- single_dataset_object()
      req(seurat_object)
      
      if(is.null(seurat_object)) {
        showNotification("Seurat object is NULL. Please reload your data.", type = "error")
        print("ERROR: seurat_object is NULL")
        return()
      }
      
      print(paste("Object cells:", ncol(seurat_object)))
      print(paste("Available assays:", paste(names(seurat_object@assays), collapse = ", ")))
      print(paste("Requested assay:", input$viz_assay))
      
      DefaultAssay(seurat_object) <- input$viz_assay
      print("DefaultAssay set successfully")
      
      # Use manual RidgePlot function
      plot <- createManualRidgePlot(
        seurat_object = seurat_object,
        genes = genes,
        group_by = "ident",
        axis_text_size = ifelse(is.null(input$axis_text_size), 12, input$axis_text_size),
        title_text_size = ifelse(is.null(input$title_text_size), 14, input$title_text_size),
        axis_line_width = ifelse(is.null(input$axis_line_width), 1, input$axis_line_width),
        add_noaxes = !is.null(input$add_noaxes_ridge) && input$add_noaxes_ridge,
        add_nolegend = !is.null(input$add_nolegend_ridge) && input$add_nolegend_ridge
      )
      
      ridge_plot(plot)
      print("RidgePlot stored successfully")
      
    }, error = function(e) {
      showNotification(paste("Error in RidgePlot:", e$message), type = "error")
      print(paste("RidgePlot error:", e$message))
      print(traceback())
    })
  })
  
  output$ridge_plot <- renderPlot({
    req(ridge_plot())
    ridge_plot()
  })

  # Download handler for RidgePlot
  output$download_ridge_plot <- createDownloadHandler(
    reactive_data = ridge_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "RidgePlot")
    }),
    data_name = "RidgePlot",
    download_type = "plot",
    plot_params = list(
      file_type = reactive(input$plot_format),  
      width = 10,
      height = 8,
      dpi = reactive(input$dpi_plot)  
    )
  )
  
  # Download handler for DotPlot
  output$downloadDotPlot <- createDownloadHandler(
    reactive_data = dot_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "DotPlot")
    }),
    data_name = "DotPlot",
    download_type = "plot",
    plot_params = list(
      file_type = reactive(input$plot_format),  
      width = 10,
      height = 8,
      dpi = reactive(input$dpi_plot)  
    )
  )
  
  # Download handler for ViolinPlot
  output$downloadVlnPlot <- createDownloadHandler(
    reactive_data = vln_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "VlnPlot")
    }),
    data_name = "VlnPlot",
    download_type = "plot",
    plot_params = list(
      file_type = reactive(input$plot_format), 
      width = 10,
      height = 8,
      dpi = reactive(input$dpi_plot) 
    )
  )
  
  # Download handler for FeaturePlot
  output$downloadFeaturePlot <- createDownloadHandler(
    reactive_data = feature_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "FeaturePlot")
    }),
    data_name = "FeaturePlot",
    download_type = "plot",
    plot_params = list(
      file_type = reactive(input$plot_format), 
      width = 10,
      height = 8,
      dpi = reactive(input$dpi_plot) 
    )
  )
  
  # Analyse de l'expression des gènes
  number_of_nuclei <- reactiveVal(NULL)

  observeEvent(input$analyze_btn, {
    req(input$gene_select_genes_analysis, single_dataset_object())
    
    tryCatch({
      result <- analyze_gene_expression(
        seurat_obj = single_dataset_object(),
        selected_genes = input$gene_select_genes_analysis,
        assay_name = input$viz_assay,
        expression_threshold = input$expression_threshold %||% 0.1,
        is_integrated = FALSE
      )
      
      number_of_nuclei(result$data)
      
      output$expression_summary <- renderDT({
        render_expression_table(result, "expression_summary")
      })
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  
  # Rendre le plot actuel dans l'UI
  output$selected_plot_display <- renderPlot({
    req(current_plot())
    current_plot()
  })

  # Réactive pour stocker et mettre à jour le plot actuel basé sur la sélection de l'utilisateur
  current_plot <- reactiveVal()
  
  
  
  
  observe({
    # Mise à jour du plot actuel en fonction de la sélection
    plot_type <- input$plot_type_select
    if (plot_type == "FeaturePlot") {
      current_plot(feature_plot())
    } else if (plot_type == "VlnPlot") {
      current_plot(vln_plot())
    } else if (plot_type == "DotPlot") {
      current_plot(dot_plot())
    } else if (plot_type == "RidgePlot") {
      current_plot(ridge_plot())
    }
  })
  
  output$download_genes_number_expresion <- createDownloadHandler(
    reactive_data = number_of_nuclei,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "GeneExpression")
    }),
    data_name = "number_of_gene_expression",
    download_type = "csv"
  )
  
  # Download handler for Seurat object with modal
  output$save_seurat_object_2 <- createDownloadHandler(
    reactive_data = single_dataset_object,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset")
    }),
    data_name = "seurat_object",
    download_type = "seurat",
    show_modal = TRUE  
  )

  ############################## Heatmap single dataset ##############################
  
  heatmap_plot <- reactiveVal()
  heatmap_n_genes <- reactiveVal(0)  # Track number of genes for dynamic height
  
  # Update selectors when data changes
  observeEvent(single_dataset_object(), {
    req(single_dataset_object())
    seurat_obj <- single_dataset_object()
    
    # Update clusters
    clusters <- levels(Idents(seurat_obj))
    updateSelectizeInput(session, "clusters_heatmap_single", 
                         choices = clusters, 
                         selected = clusters,
                         server = TRUE)
    
    # Update split_by with filtered metadata
    excluded_patterns <- c(
      "percent.mt", "nCount_", "nFeature_", "^RNA_snn_", "^RNA_nn_",              
      "^Spatial_snn_", "^Spatial_nn_", "original_clusters", "^pANN", "^DF",
      "dataset_origin", "original_seurat_clusters", "original_ClusterIdents",
      "cluster_name_only", "source_format", "^integrated", "^integrated_snn",
      "^merge_method$", "^integration_method$", "^original_clusterIdents$", "^cluster_names$" 
    )
    pattern <- paste(excluded_patterns, collapse = "|")
    meta_cols <- colnames(seurat_obj@meta.data)
    filtered_cols <- meta_cols[!grepl(pattern, meta_cols)]
    
    updateSelectInput(session, "split_by_heatmap_single", 
                      choices = c("None" = "None", filtered_cols),
                      selected = "None")
  })
  
  # Store top markers table for CSV export (NULL when using manual gene input)
  heatmap_markers_data_single <- reactiveVal(NULL)
  
  # Generate heatmap - SIMPLE average expression only
  # Store top markers table for CSV export (NULL when using manual gene input)
  heatmap_markers_data_single <- reactiveVal(NULL)
  
  # Generate heatmap
  observeEvent(input$generateCustomHeatmap, {
    showModal(modalDialog(
      title = "Processing",
      div(h4("Generating Heatmap...", style = "text-align: center;")),
      footer = NULL,
      easyClose = FALSE
    ))
    req(single_dataset_object())
    seurat_object <- single_dataset_object()
    tryCatch({
      selected_assay <- input$assay_select_heatmap_single
      if (is.null(selected_assay)) selected_assay <- "RNA"
      if (!(selected_assay %in% names(seurat_object@assays))) selected_assay <- "RNA"
      
      # Drop NA idents at server level before anything else
      idents_vec <- as.character(Idents(seurat_object))
      valid_cells <- colnames(seurat_object)[!is.na(idents_vec)]
      if (length(valid_cells) < ncol(seurat_object)) {
        message("[Heatmap] Dropping ", ncol(seurat_object) - length(valid_cells), " cells with NA idents")
        seurat_object <- subset(seurat_object, cells = valid_cells)
      }
      Idents(seurat_object) <- droplevels(Idents(seurat_object))
      
      message("[DEBUG] Idents levels: ", paste(levels(Idents(seurat_object)), collapse = ", "))
      message("[DEBUG] Idents NA count: ", sum(is.na(as.character(Idents(seurat_object)))))
      message("[DEBUG] counts class: ", class(LayerData(seurat_object, assay = "RNA", layer = "counts")))
      message("[DEBUG] assays: ", paste(names(seurat_object@assays), collapse = ", "))
      
      DefaultAssay(seurat_object) <- selected_assay
      clusters_selected <- input$clusters_heatmap_single
      if (is.null(clusters_selected) || length(clusters_selected) == 0) {
        clusters_selected <- levels(Idents(seurat_object))
      }
      if (input$use_top10_genes) {
        n_genes <- ifelse(is.null(input$n_top_genes_single), 10, input$n_top_genes_single)
        markers_result <- get_top_markers(seurat_object, n_genes = n_genes, assay = selected_assay)
        genes <- markers_result$genes
        heatmap_markers_data_single(markers_result$markers_table)
        message(paste("Using top", n_genes, "genes per cluster"))
      } else {
        gene_text <- trimws(input$gene_select_heatmap)
        if (nchar(gene_text) == 0) {
          showNotification("Please enter genes.", type = "error")
          removeModal()
          return()
        }
        gene_result <- parse_and_validate_genes(gene_text, seurat_object, selected_assay)
        genes <- gene_result$valid_genes
        heatmap_markers_data_single(NULL)
        if (length(gene_result$missing_genes) > 0) {
          showNotification(paste(length(gene_result$missing_genes), "gene(s) not found"),
                           type = "warning", duration = 6)
        }
      }
      heatmap_n_genes(length(genes))
      color_palette <- input$color_palette_heatmap_single
      if (is.null(color_palette)) color_palette <- "viridis"
      scale_rows <- TRUE
      split_by <- input$split_by_heatmap_single
      if (is.null(split_by) || split_by == "None") split_by <- NULL
      plot <- generate_split_heatmap(
        seurat_object = seurat_object,
        genes = genes,
        clusters = clusters_selected,
        assay = selected_assay,
        split_by = split_by,
        scale_rows = scale_rows,
        color_palette = color_palette
      )
      heatmap_plot(plot)
      showNotification(paste0("Heatmap with ", length(genes), " genes"), type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", paste(e$message, collapse = " ")), type = "error")
      message(paste("Full error:", e$message))
    })
    removeModal()
  })
  # Download top markers table as CSV (only available in top markers mode)
  output$download_heatmap_markers_single <- downloadHandler(
    filename = function() {
      obj_name <- getObjectNameForDownload(single_dataset_object(), default_name = "object")
      paste0(obj_name, "_top_markers_per_cluster_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      markers_df <- heatmap_markers_data_single()
      req(markers_df)
      write.csv(markers_df, file, row.names = FALSE)
    }
  )
  
  # Render heatmap with dynamic height
  output$heatmap_single <- renderPlot({
    req(heatmap_plot())
    heatmap_plot()
  }, height = function() {
    # Dynamic height calculation: minimum 400px, then 20px per gene
    n_genes <- heatmap_n_genes()
    if (n_genes == 0) return(600)  # Default height
    max(400, min(20 * n_genes, 2000))  # Min 400px, max 2000px, 20px per gene
  })
  
  # Download handler with DYNAMIC HEIGHT
  output$download_heatmap_single <- createDownloadHandler(
    reactive_data = heatmap_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "Heatmap")
    }),
    data_name = "heatmap",
    download_type = "plot",
    plot_params = list(
      file_type = reactive(input$format_heatmap_single), 
      width = 10,
      height = reactive({
        # Dynamic height based on number of genes
        n_genes <- heatmap_n_genes()
        if (n_genes == 0) return(8)  # Default height
        # Calculate: minimum 6, then 0.15 inch per gene, max 50
        max(6, min(0.15 * n_genes, 50))
      }),
      dpi = reactive(input$dpi_heatmap_single)  
    )
  )
  
  
  # Mettez à jour les choix de selectInput pour les gènes/features et les clusters
  observe({
    req(single_dataset_object())
    updateGeneChoices(session, single_dataset_object(), "RNA", c("feature1_select", "feature2_select"))
    updateClusterTextInputs(session, single_dataset_object(), c("scatter_text_clusters", "text_clusters"))
  })

  scatter_plot_single <- reactiveVal()
  
  # Update color_by and split_by choices when data changes
  observeEvent(single_dataset_object(), {
    req(single_dataset_object())
    seurat_obj <- single_dataset_object()
    
    # Filter metadata columns
    excluded_patterns <- c(
      "percent.mt", "nCount_", "nFeature_", "^RNA_snn_", "^RNA_nn_",              
      "^Spatial_snn_", "^Spatial_nn_", "original_clusters", "^pANN", "^DF",
      "dataset_origin", "original_seurat_clusters", "original_ClusterIdents",
      "cluster_name_only", "source_format", "^integrated", "^integrated_snn"
    )
    pattern <- paste(excluded_patterns, collapse = "|")
    meta_cols <- colnames(seurat_obj@meta.data)
    filtered_cols <- meta_cols[!grepl(pattern, meta_cols)]
    
    # Update color_by
    updateSelectInput(session, "color_by_scatter_single", 
                      choices = c("Cluster" = "cluster", filtered_cols),
                      selected = "cluster")
    
    # Update split_by
    updateSelectInput(session, "split_by_scatter_single", 
                      choices = c("None" = "None", filtered_cols),
                      selected = "None")
  })
  
  # Generate scatter plot
  observeEvent(input$generateScatter_single, {
    req(single_dataset_object())
    
    showModal(modalDialog(
      title = "Processing",
      div(h4("Generating Scatter Plot...", style = "text-align: center;")),
      footer = NULL,
      easyClose = FALSE
    ))
    
    seurat_obj <- single_dataset_object()
    
    tryCatch({
      # Get genes
      gene1 <- trimws(input$feature1_select)
      gene2 <- trimws(input$feature2_select)
      
      if (gene1 == "" || gene2 == "") {
        showNotification("Please enter both gene names.", type = "error")
        removeModal()
        return()
      }
      
      selected_assay <- "RNA"
      
      # Handle cluster filtering
      cluster_text <- trimws(input$scatter_text_clusters)
      if (cluster_text == "") {
        clusters_to_use <- NULL
      } else {
        clusters_to_use <- trimws(unlist(strsplit(cluster_text, ",")))
        clusters_to_use <- clusters_to_use[nchar(clusters_to_use) > 0]
        if (length(clusters_to_use) == 0) clusters_to_use <- NULL
      }
      
      # Get display options
      color_by <- input$color_by_scatter_single
      if (is.null(color_by)) color_by <- "cluster"
      
      split_by <- input$split_by_scatter_single
      if (is.null(split_by) || split_by == "None") split_by <- NULL
      
      # Get thresholds
      threshold1 <- input$threshold_gene1_scatter_single
      threshold2 <- input$threshold_gene2_scatter_single
      if (is.null(threshold1)) threshold1 <- 0
      if (is.null(threshold2)) threshold2 <- 0
      
      # Generate plot
      plot <- generate_feature_scatter(
        seurat_object = seurat_obj,
        gene1 = gene1,
        gene2 = gene2,
        clusters = clusters_to_use,
        assay = selected_assay,
        color_by = color_by,
        split_by = split_by,
        threshold_gene1 = threshold1,
        threshold_gene2 = threshold2
      )
      
      scatter_plot_single(plot)
      showNotification("Scatter plot generated successfully.", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      message(paste("Full error:", e$message))
    })
    
    removeModal()
  })
  
  
  
  output$scatter_plot_single <- renderPlot({
    req(scatter_plot_single())
    scatter_plot_single()
  })
  
  output$download_scatter_single <- createDownloadHandler(
    reactive_data = scatter_plot_single,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "FeatureScatter")
    }),
    data_name = "scatter",
    download_type = "plot",
    plot_params = list(
      file_type = reactive(input$format_scatter_single), 
      width = 10,
      height = 8,
      dpi = reactive(input$dpi_scatter_single)  
    )
  )


  ############################## Final UMAP ##############################
  # Tab 8: Final UMAP

  # Save initial cluster state for undo functionality
  observe({
    req(single_dataset_object())
    if(length(cluster_history_single$states) == 0) {
      save_cluster_state_single(single_dataset_object(), "Initial state")
    }
  })
  
 ##################### Cluster History System for Undo ####################
  
  # Simple history system for undo functionality
  cluster_history_single <- reactiveValues(
    states = list(),
    max_history = 10
  )
  
  # Function to save current cluster state
  save_cluster_state_single <- function(seurat_obj, action_description = "Cluster modification") {
    current_clusters <- as.character(Idents(seurat_obj))
    names(current_clusters) <- names(Idents(seurat_obj))
    
    new_state <- list(
      clusters = current_clusters,
      action = action_description
    )
    
    cluster_history_single$states <- append(cluster_history_single$states, list(new_state))
    
    # Limit history size
    if(length(cluster_history_single$states) > cluster_history_single$max_history) {
      cluster_history_single$states <- cluster_history_single$states[-1]
    }
  }
  
  # Function to undo last action
  undo_last_action_single <- function() {
    if(length(cluster_history_single$states) <= 1) return(NULL)
    
    # Remove last state
    cluster_history_single$states <- cluster_history_single$states[-length(cluster_history_single$states)]
    
    # Return previous state
    if(length(cluster_history_single$states) > 0) {
      return(cluster_history_single$states[[length(cluster_history_single$states)]]$clusters)
    }
    return(NULL)
  }
  
  
  
  
  
  #Reactive variable for that tab
  cluster_colours <- reactiveVal()  # Initialize a list to store cluster colors

  observe({
    req(single_dataset_object())
    updateClusterChoices(session, single_dataset_object(), list(select = c("select_cluster", "cluster_select")))
  })


  # Function for renaming each cluster
  observeEvent(input$rename_single_cluster_button, {
    req(input$select_cluster, input$rename_single_cluster, single_dataset_object())
    
    # SAVE STATE BEFORE MODIFICATION
    save_cluster_state_single(single_dataset_object())
    
    updated_seurat <- single_dataset_object()
    
    if (input$rename_single_cluster %in% unique(Idents(updated_seurat))) {
      cells_to_merge <- which(Idents(updated_seurat) %in% c(input$select_cluster, input$rename_single_cluster))
      Idents(updated_seurat, cells = cells_to_merge) <- input$rename_single_cluster
      showNotification(paste("Clusters merged under the name:", input$rename_single_cluster), type = "message")
    } else {
      Idents(updated_seurat, cells = which(Idents(updated_seurat) == input$select_cluster)) <- input$rename_single_cluster
      showNotification(paste("Cluster renamed to:", input$rename_single_cluster), type = "message")
    }
    
    updated_seurat$cluster_names <- as.character(Idents(updated_seurat))
    
    if(!all(sort(unique(as.character(Idents(updated_seurat)))) == 
            sort(unique(as.character(updated_seurat$cluster_names))))) {
      showNotification("Warning: Cluster names mismatch between Idents and metadata", type = "warning")
    } else {
      showNotification("Cluster names synchronized in metadata", type = "message")
    }
    
    single_dataset_object(updated_seurat)
    updateSelectInput(session, "select_cluster", choices = unique(Idents(single_dataset_object())))
    updateSelectInput(session, "cluster_select", choices = unique(Idents(single_dataset_object())))
  })
  
  # Observer for Undo button
  observeEvent(input$undo_cluster_single, {
    req(single_dataset_object())
    
    restored_clusters <- undo_last_action_single()
    
    if(!is.null(restored_clusters)) {
      updated_seurat <- single_dataset_object()
      Idents(updated_seurat) <- restored_clusters
      updated_seurat$cluster_names <- as.character(Idents(updated_seurat))
      single_dataset_object(updated_seurat)
      
      # Update UI
      updateSelectInput(session, "select_cluster", choices = unique(Idents(updated_seurat)))
      updateSelectInput(session, "cluster_select", choices = unique(Idents(updated_seurat)))
      
      showNotification("Undid last action", type = "message")
    } else {
      showNotification("Nothing to undo", type = "warning")
    }
  })
  
  get_cluster_colors <- function(seurat_object) {
    if (is.null(seurat_object)) return(NULL)
    colors <- find_stored_colors(seurat_object, cluster_names = levels(Idents(seurat_object)))
    if (!is.null(colors)) {
      message("Using stored cluster colors")
      return(colors)
    }
    message("Using Seurat default colors")
    return(NULL)
  }
  
  # Update color picker to reflect the stored color of the selected cluster
  observeEvent(input$cluster_select, {
    req(single_dataset_object(), input$cluster_select)
    colors <- get_cluster_colors(single_dataset_object())
    if (!is.null(colors) && input$cluster_select %in% names(colors)) {
      updateColourInput(session, "cluster_colour", value = colors[[input$cluster_select]])
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$update_colour, {
    message("Update the color of the selected cluster")
    updated_seurat <- single_dataset_object()
    cluster_colors <- get_cluster_colors(updated_seurat)
    
    if (is.null(cluster_colors)) {
      message("No colors stored, creating default colors")
      clusters      <- levels(Idents(updated_seurat))
      n_clusters    <- length(clusters)
      default_colors <- scales::hue_pal()(n_clusters)
      names(default_colors) <- clusters
      cluster_colors <- default_colors
      message("Created default colors for ", n_clusters, " clusters")
    }
    
    if (!is.null(input$cluster_select) && input$cluster_select %in% names(cluster_colors)) {
      cluster_colors[input$cluster_select] <- input$cluster_colour
      message(paste("Updating color for cluster", input$cluster_select, "to", input$cluster_colour))
    } else {
      message("Selected cluster: ", input$cluster_select)
      message("Available clusters: ", paste(names(cluster_colors), collapse = ", "))
      showNotification("Selected cluster is not valid.", type = "error")
      return()
    }
    
    updated_seurat@misc$cluster_colors <- cluster_colors
    single_dataset_object(updated_seurat)
    message("Current cluster colors:")
    print(updated_seurat@misc$cluster_colors)
    showNotification("Cluster color updated!", type = "message")
  })
  
  # Display final UMAP with 2D/3D support
  output$umap_finale <- renderPlotly({
    req(single_dataset_object())
    message("Generating the final UMAP")
    
    seurat_obj <- single_dataset_object()
    
    # Check if 3D view is requested
    use_3d <- isTRUE(input$umap_3d_toggle)
    
    # Get cluster colors
    cluster_colors <- get_cluster_colors(seurat_obj)
    if (is.null(cluster_colors)) {
      n_clusters <- length(unique(Idents(seurat_obj)))
      cluster_colors <- scales::hue_pal()(n_clusters)
      names(cluster_colors) <- levels(Idents(seurat_obj))
      message("Using Seurat default colors")
    } else {
      message("Using custom cluster colors")
    }
    
    # Generate plot based on dimensionality
    if (use_3d && hasValid3DUMAP(seurat_obj)) {
      # Extract 3D UMAP coordinates
      umap_coords <- as.data.frame(Embeddings(seurat_obj, reduction = "umap3d"))
      colnames(umap_coords) <- c("UMAP_1", "UMAP_2", "UMAP_3")
      umap_coords$cluster <- Idents(seurat_obj)
      
      # Calculate cluster centers for labels
      centers <- umap_coords %>%
        group_by(cluster) %>%
        summarise(
          UMAP_1 = median(UMAP_1),
          UMAP_2 = median(UMAP_2),
          UMAP_3 = median(UMAP_3),
          .groups = 'drop'
        )
      
      # Get display options
      hide_grid <- isTRUE(input$umap_3d_hide_grid)
      hide_axes <- isTRUE(input$umap_3d_hide_axes)
      dark_mode <- isTRUE(input$umap_3d_dark_mode)
      
      # Set colors based on mode
      bg_color <- if(dark_mode) "#1a1a1a" else "white"
      grid_color <- if(dark_mode) "#404040" else "#e0e0e0"
      axis_color <- if(dark_mode) "#ffffff" else "#000000"
      text_color <- if(dark_mode) "#ffffff" else "#000000"
      
      # Configure axis settings
      axis_config <- list(
        title = "",
        showgrid = !hide_grid,
        showline = !hide_axes,
        zeroline = !hide_axes,
        showticklabels = !hide_axes,
        gridcolor = grid_color,
        gridwidth = 1,
        linecolor = axis_color,
        tickcolor = axis_color,
        tickfont = list(color = text_color)
      )
      
      # Add axis titles if axes are visible
      if (!hide_axes) {
        axis_config$title <- list(text = "", font = list(color = text_color))
      }
      
      # Create axis configs for each dimension
      xaxis_config <- axis_config
      yaxis_config <- axis_config
      zaxis_config <- axis_config
      
      if (!hide_axes) {
        xaxis_config$title <- list(text = "UMAP_1", font = list(color = text_color))
        yaxis_config$title <- list(text = "UMAP_2", font = list(color = text_color))
        zaxis_config$title <- list(text = "UMAP_3", font = list(color = text_color))
      }
      
      # Create 3D plotly
      p <- plot_ly() %>%
        add_trace(
          data = umap_coords,
          x = ~UMAP_1, 
          y = ~UMAP_2, 
          z = ~UMAP_3,
          color = ~cluster,
          colors = cluster_colors,
          type = "scatter3d",
          mode = "markers",
          marker = list(size = input$pt_size, opacity = 0.7),
          hoverinfo = "text",
          text = ~paste("Cluster:", cluster)
        ) %>%
        add_trace(
          data = centers,
          x = ~UMAP_1,
          y = ~UMAP_2,
          z = ~UMAP_3,
          type = "scatter3d",
          mode = "text",
          text = ~cluster,
          textfont = list(size = input$label_font_size, color = text_color),
          hoverinfo = "none",
          showlegend = FALSE
        ) %>%
        layout(
          scene = list(
            xaxis = xaxis_config,
            yaxis = yaxis_config,
            zaxis = zaxis_config,
            bgcolor = bg_color,
            camera = list(
              eye = list(x = 1.5, y = 1.5, z = 1.5)
            )
          ),
          paper_bgcolor = bg_color,
          plot_bgcolor = bg_color,
          showlegend = FALSE,
          title = list(
            text = input$plot_title, 
            font = list(size = 24, color = text_color)
          )
        )
      
      message("3D UMAP generated with custom display options")
      return(p)
      
    } else if (use_3d && !hasValid3DUMAP(seurat_obj)) {
      # 3D requested but not computed yet
      empty_plot <- plotly_empty() %>%
        layout(
          title = list(
            text = "Please compute 3D UMAP first using the button above",
            font = list(size = 16)
          )
        )
      return(empty_plot)
      
    } else {
      # 2D UMAP (existing code - unchanged)
      plot_data <- DimPlot(
        seurat_obj,
        group.by = "ident",
        pt.size = input$pt_size,
        label = TRUE,
        label.size = input$label_font_size
      ) +
        NoLegend() +
        theme(axis.line = element_line(size = 0.5))
      
      if (!is.null(cluster_colors)) {
        plot_data <- plot_data + scale_color_manual(values = cluster_colors)
      }
      
      # Convert to interactive plotly
      interactive_plot <- ggplotly(plot_data, tooltip = "text") %>%
        layout(
          title = list(text = input$plot_title, font = list(size = 24)),
          hovermode = "closest"
        )
      
      message("2D UMAP generated")
      return(interactive_plot)
    }
  })
  
  output$save_seurat_object_1 <- createDownloadHandler(
    reactive_data = single_dataset_object,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset") }),
    data_name = "seurat_object",
    download_type = "seurat",
    show_modal = TRUE
  )
  
  ############################## Find markers for a specific cluster ##############################

  # Variables réactives pour cet onglet
  markers_global_single <- reactiveVal()
  umap_plot <- reactiveVal()
  
  # Function to calculate all differential markers once only
  clean_gene_names_for_html <- function(gene_names) {
    tryCatch({
      original_names <- rownames(single_dataset_object())
      cleaned_names <- sapply(gene_names, function(gene) {
        base_name <- gsub("\\.\\d+$", "", gene) 
        if (base_name %in% original_names) {
          return(base_name)
        }
        return(gene)  
      })
      
      return(cleaned_names)
    }, error = function(e) {
      showNotification(paste("Error cleaning gene names:", e$message), type = "error")
      return(gene_names) 
    })
  }
  
  # UMAP pour un cluster spécifique avec couleurs mises à jour
  output$umap_plot <- renderPlot({
    req(single_dataset_object())
    message("UMAP generation for a specific cluster")
    label_option <- input$show_labels  
    plot_data <- DimPlot(
      single_dataset_object(),
      group.by = "ident",
      pt.size = input$pt_size,
      label = label_option, 
      label.size = input$label_font_size
    ) +
      theme(axis.line = element_line(size = 0.5)) +
      theme_void() +
      theme(legend.position = "none") 
    
    cluster_colors <- get_cluster_colors(single_dataset_object())
    if (!is.null(cluster_colors)) {
      plot_data <- plot_data + scale_color_manual(values = cluster_colors)
      message("Using custom cluster colors")
    } else {
      message("Using Seurat default colors")
    }
    
    plot_data <- plot_data + NoLegend() + ggtitle(NULL)
    umap_plot(plot_data)
    message("UMAP for a specific cluster generated")
    return(plot_data)
  })
  
  # Update cluster choices for global comparison
  observe({
    req(single_dataset_object())
    updateSelectInput(session, "target_cluster_global_single",
                      choices = levels(single_dataset_object()))
  })
  
  # Update DE assay choices when object is loaded
  observe({
    req(single_dataset_object())
    available_assays <- names(single_dataset_object()@assays)
    updateSelectInput(session, "assay_de_single",
                      choices = available_assays,
                      selected = if ("RNA" %in% available_assays) "RNA" else available_assays[1])
  })
  
  
  # Download handler for UMAP cluster plot
  output$downloadUMAPCluster <- createDownloadHandler(
    reactive_data = umap_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "UMAP_plot") 
    }),
    data_name = "UMAP",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$umap_cluster_format }),  
      width = reactive({ ifelse(input$umap_cluster_format == "pdf", 11, 10) }),  
      height = reactive({ ifelse(input$umap_cluster_format == "pdf", 8, 6) }),  
      dpi = reactive({ input$dpi_umap_cluster })  
    )
  )
  
  # Download handler for Seurat object
  output$save_seurat_object_3 <- createDownloadHandler(
    reactive_data = single_dataset_object,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset") 
    }),
    data_name = "seurat_object",
    download_type = "seurat",
    show_modal = TRUE
  )
  
  calculate_markers_for_cluster <- function() {
    req(single_dataset_object(),
        input$target_cluster_global_single,
        input$min_pct_single,
        input$logfc_threshold_single)
    
    tryCatch({
      seurat_obj_global <- JoinLayers(single_dataset_object(), assay = input$assay_de_single)
      markers <- FindMarkers(
        seurat_obj_global,
        ident.1 = input$target_cluster_global_single,
        min.pct = input$min_pct_single,
        logfc.threshold = input$logfc_threshold_single,
        assay = input$assay_de_single,
        return.thresh = 0.05
      )
      
      # Sort by absolute log2 fold change (most differential genes first)
      markers <- markers[order(abs(markers$avg_log2FC), decreasing = TRUE), ]
      
      # Store in global list for Venn diagrams BEFORE any modifications (rule #7)
      markers_copy_venn <- markers
      markers_copy_venn$gene <- rownames(markers_copy_venn)
      
      table_name <- paste0("SingleDataset_Cluster_",
                           input$target_cluster_global_single,
                           "_vs_All_",
                           format(Sys.time(), "%H%M%S"))
      description <- paste0("Cluster ", input$target_cluster_global_single, " vs All Others")
      parameters <- list(
        min_pct = input$min_pct_single,
        logfc_threshold = input$logfc_threshold_single
      )
      
      result <- storeDETable(single_gene_table_storage(),
                             markers_copy_venn,
                             table_name,
                             description,
                             "cluster_group",
                             parameters)
      if (result$success) {
        single_gene_table_storage(result$storage)
      }
      
      # Store raw numeric markers — formatting and filtering happen at render time
      markers$gene <- rownames(markers)
      markers_global_single(markers)
      
      return(TRUE)
    }, error = function(e) {
      showNotification(paste0("Error finding markers: ", e$message), type = "error")
      return(FALSE)
    })
  }
  
  
  # Observer for the find markers button
  observeEvent(input$find_markers_global_single, {
    tryCatch({
      showModal(modalDialog(
        title = "Finding Markers",
        "Comparing cluster against all others...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      calculate_markers_for_cluster()
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Error finding markers: ", e$message), type = "error")
    })
  })
  
  # Compute 3D UMAP when button is clicked
  observeEvent(input$compute_3d_umap, {
    req(single_dataset_object())
    
    showModal(modalDialog(
      title = "Computing 3D UMAP",
      "Calculating 3D UMAP coordinates... This may take a moment.",
      easyClose = FALSE,
      footer = NULL
    ))
    
    tryCatch({
      seurat_obj <- single_dataset_object()
      
      # Get dims from existing PCA
      available_pcs <- ncol(seurat_obj@reductions$pca@cell.embeddings)
      umap_dims <- min(30, available_pcs)
      
      # Compute 3D UMAP
      seurat_obj <- runUMAP3D_reproducible(
        object = seurat_obj,
        dims = 1:umap_dims,
        seed = 42,
        verbose = TRUE
      )
      
      # Update reactive value
      single_dataset_object(seurat_obj)
      
      removeModal()
      showNotification("3D UMAP computed successfully!", type = "message")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error computing 3D UMAP:", e$message), type = "error")
      message("Error in 3D UMAP computation: ", e$message)
    })
  })
  
  # Display results table
  output$table_global_single <- renderDataTable({
    req(markers_global_single())
    df <- markers_global_single()
    threshold <- input$pval_adj_filter_global_single
    if (!is.null(threshold) && !is.na(threshold)) {
      df <- df[df$p_val_adj <= threshold, , drop = FALSE]
    }
    cleaned <- clean_gene_names_for_html(df$gene)
    df$gene_link <- paste0('<a href="#" class="gene-name" data-gene="', cleaned, '">', cleaned, '</a>')
    df$p_val     <- format_pvalue_robust(df$p_val)
    df$p_val_adj <- format_pvalue_robust(df$p_val_adj)
    df <- df[, c("gene_link", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
    colnames(df) <- c("Gene", "Log2FC", "Pct.1", "Pct.2", "P-value", "Adj. P-value")
    datatable(
      df,
      escape = FALSE,
      options = list(
        pageLength = 10,
        lengthMenu = c(10, 25, 50, 100, 200, 500),
        order = list(list(1, "desc")),
        dom = "Blfrtip",
        scrollX = TRUE,
        columnDefs = list(
          list(className = "dt-center", targets = 1:5),
          list(width = "150px", targets = 0)
        )
      ),
      class = "cell-border stripe",
      rownames = FALSE
    )
  })
  
  # Download handler for global comparison (UN SEUL!)
  output$download_markers_global_single <- createDownloadHandler(
    reactive_data = reactive({
      df <- markers_global_single()
      req(df)
      threshold <- input$pval_adj_filter_global_single
      if (!is.null(threshold) && !is.na(threshold)) {
        df <- df[df$p_val_adj <= threshold, , drop = FALSE]
      }
      df$p_val     <- format_pvalue_robust(df$p_val)
      df$p_val_adj <- format_pvalue_robust(df$p_val_adj)
      colnames(df)[colnames(df) == "gene"]        <- "Gene"
      colnames(df)[colnames(df) == "avg_log2FC"]  <- "Log2FC"
      colnames(df)[colnames(df) == "pct.1"]       <- "Pct.1"
      colnames(df)[colnames(df) == "pct.2"]       <- "Pct.2"
      colnames(df)[colnames(df) == "p_val"]       <- "P-value"
      colnames(df)[colnames(df) == "p_val_adj"]   <- "Adj. P-value"
      return(df)
    }),
    object_name_reactive = reactive({
      getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset")
    }),
    data_name = reactive({
      paste0("Cluster_", input$target_cluster_global_single, "_vs_All_markers")
    }),
    download_type = "csv"
  )
  
  ############Pairwise Cluster Comparison###########
  
  # Initialize reactive variable
  markers_pairwise_single <- reactiveVal()
  
  # Update cluster choices
  observe({
    req(single_dataset_object())
    updateSelectInput(session, "cluster1_pairwise_single", 
                      choices = levels(single_dataset_object()))
  })
  
  observeEvent(input$cluster1_pairwise_single, {
    updateSelectInput(session, "cluster2_pairwise_single",
                      choices = setdiff(levels(single_dataset_object()), 
                                        input$cluster1_pairwise_single))
  })
  
  # Function to calculate markers for pairwise comparison
  calculate_markers_for_comparison <- function() {
    req(single_dataset_object(),
        input$cluster1_pairwise_single,
        input$cluster2_pairwise_single,
        input$min_pct_pairwise_single,
        input$logfc_threshold_pairwise_single)
    
    tryCatch({
      # Check for overlap between groups
      if (any(input$cluster1_pairwise_single %in% input$cluster2_pairwise_single)) {
        showNotification("Cluster groups must be distinct", type = "error")
        return(FALSE)
      }
      
      seurat_obj <- single_dataset_object()
      
      # Create new identity for comparison groups
      new_idents <- as.character(Idents(seurat_obj))
      new_idents[new_idents %in% input$cluster1_pairwise_single] <- "group1"
      new_idents[new_idents %in% input$cluster2_pairwise_single] <- "group2"
      Idents(seurat_obj) <- new_idents
      
      seurat_obj <- JoinLayers(seurat_obj, assay = input$assay_de_single)
      markers <- FindMarkers(seurat_obj,
                             ident.1 = "group1",
                             ident.2 = "group2",
                             min.pct = input$min_pct_pairwise_single,
                             logfc.threshold = input$logfc_threshold_pairwise_single,
                             assay = input$assay_de_single,
                             return.thresh = 0.05)
      
      # Sort by absolute log2 fold change
      markers <- markers[order(abs(markers$avg_log2FC), decreasing = TRUE), ]
      
      # Store in global list for Venn diagrams BEFORE any modifications (rule #7)
      markers_copy_venn <- markers
      markers_copy_venn$gene <- rownames(markers_copy_venn)
      
      group1_text <- paste(input$cluster1_pairwise_single, collapse = "_")
      group2_text <- paste(input$cluster2_pairwise_single, collapse = "_")
      table_name <- paste0("SingleDataset_Clusters_", group1_text, "_vs_",
                           group2_text, "_", format(Sys.time(), "%H%M%S"))
      description <- paste0("Clusters [", group1_text, "] vs [", group2_text, "]")
      parameters <- list(
        min_pct = input$min_pct_pairwise_single,
        logfc_threshold = input$logfc_threshold_pairwise_single,
        group1 = input$cluster1_pairwise_single,
        group2 = input$cluster2_pairwise_single
      )
      
      result <- storeDETable(single_gene_table_storage(),
                             markers_copy_venn,
                             table_name,
                             description,
                             "cluster_group",
                             parameters)
      if (result$success) {
        single_gene_table_storage(result$storage)
      }
      
      # Store raw numeric markers with comparison label — formatting and filtering happen at render time
      markers$gene <- rownames(markers)
      markers$comparison <- sprintf("Group1(%s) vs Group2(%s)",
                                    paste(input$cluster1_pairwise_single, collapse = ","),
                                    paste(input$cluster2_pairwise_single, collapse = ","))
      markers_pairwise_single(markers)
      
      return(TRUE)
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      return(FALSE)
    })
  }
  
  # Observer for the compare markers button
  observeEvent(input$compare_markers_pairwise_single, {
    tryCatch({
      showModal(modalDialog(
        title = "Comparing Clusters",
        "Finding differential markers between cluster groups...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      calculate_markers_for_comparison()
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Error comparing markers: ", e$message), type = "error")
    })
  })
  
  # Display results table
  output$table_pairwise_single <- renderDataTable({
    req(markers_pairwise_single())
    df <- markers_pairwise_single()
    threshold <- input$pval_adj_filter_pairwise_single
    if (!is.null(threshold) && !is.na(threshold)) {
      df <- df[df$p_val_adj <= threshold, , drop = FALSE]
    }
    cleaned <- clean_gene_names_for_html(df$gene)
    df$gene_link <- paste0('<a href="#" class="gene-name" data-gene="', cleaned, '">', cleaned, '</a>')
    df$p_val     <- format_pvalue_robust(df$p_val)
    df$p_val_adj <- format_pvalue_robust(df$p_val_adj)
    df <- df[, c("gene_link", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj", "comparison")]
    colnames(df) <- c("Gene", "Log2FC", "Pct.1", "Pct.2", "P-value", "Adj. P-value", "Comparison")
    datatable(
      df,
      escape = FALSE,
      options = list(
        pageLength = 10,
        lengthMenu = c(10, 25, 50, 100, 200, 500),
        order = list(list(1, "desc")),
        dom = "Blfrtip",
        scrollX = TRUE,
        columnDefs = list(
          list(className = "dt-center", targets = 1:6),
          list(width = "150px", targets = 0)
        )
      ),
      class = "cell-border stripe",
      rownames = FALSE
    )
  })
  
  # Download handler for pairwise comparison
  output$download_markers_pairwise_single <- createDownloadHandler(
    reactive_data = reactive({
      df <- markers_pairwise_single()
      req(df)
      threshold <- input$pval_adj_filter_pairwise_single
      if (!is.null(threshold) && !is.na(threshold)) {
        df <- df[df$p_val_adj <= threshold, , drop = FALSE]
      }
      df$p_val     <- format_pvalue_robust(df$p_val)
      df$p_val_adj <- format_pvalue_robust(df$p_val_adj)
      colnames(df)[colnames(df) == "gene"]        <- "Gene"
      colnames(df)[colnames(df) == "avg_log2FC"]  <- "Log2FC"
      colnames(df)[colnames(df) == "pct.1"]       <- "Pct.1"
      colnames(df)[colnames(df) == "pct.2"]       <- "Pct.2"
      colnames(df)[colnames(df) == "p_val"]       <- "P-value"
      colnames(df)[colnames(df) == "p_val_adj"]   <- "Adj. P-value"
      colnames(df)[colnames(df) == "comparison"]  <- "Comparison"
      return(df)
    }),
    object_name_reactive = reactive({
      getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset")
    }),
    data_name = reactive({
      paste0("Clusters_",
             paste(input$cluster1_pairwise_single, collapse = "_"),
             "_vs_",
             paste(input$cluster2_pairwise_single, collapse = "_"),
             "_markers")
    }),
    download_type = "csv"
  )
  
  
  
  
  #####################Exclusive Biomarkers Discovery###############################
  
  # Initialize reactive variable for exclusive markers
  exclusive_markers_single <- reactiveVal()
  
  # Update target cluster choices for exclusive biomarkers tab
  observe({
    req(single_dataset_object())
    updateSelectInput(session, "exclusive_target_cluster_single",
                      choices = levels(single_dataset_object()))
  })
  
  # Find exclusive biomarkers when button is clicked
  observeEvent(input$find_exclusive_markers_single, {
    req(single_dataset_object(), input$exclusive_target_cluster_single)
    tryCatch({
      showModal(modalDialog(
        title = "Finding Exclusive Biomarkers",
        "Analyzing gene expression patterns across clusters...",
        easyClose = FALSE,
        footer = NULL
      ))
      markers <- find_exclusive_biomarkers(
        seurat_obj = single_dataset_object(),
        target_cluster = input$exclusive_target_cluster_single,
        min_pct_target = input$exclusive_min_pct_target_single,
        max_pct_other = input$exclusive_max_pct_other_single,
        min_log2fc = input$exclusive_min_log2fc_single,
        detection_threshold = input$exclusive_detection_threshold_single,
        min_mean_expr_target = input$exclusive_min_mean_expr_single,
        statistical_test = input$exclusive_statistical_test_single,
        max_pvalue = input$exclusive_max_pvalue_single,
        pvalue_adjustment = "BH",
        assay_name = DefaultAssay(single_dataset_object()),
        top_n = NULL,
        verbose = TRUE
      )
      removeModal()
      if (nrow(markers) == 0) {
        showNotification("No exclusive biomarkers found with current parameters. Try relaxing the filters.",
                         type = "error",
                         duration = 10)
        exclusive_markers_single(NULL)
        return()
      }
      # Store results for Venn diagrams BEFORE formatting p-values
      markers_copy_venn <- markers
      markers_copy_venn$p_val <- ifelse("pvalue" %in% colnames(markers_copy_venn),
                                        markers_copy_venn$pvalue, NA)
      markers_copy_venn$p_val_adj <- ifelse("pvalue_adjusted" %in% colnames(markers_copy_venn),
                                            markers_copy_venn$pvalue_adjusted, NA)
      markers_copy_venn$avg_log2FC <- markers_copy_venn$log2_fold_change
      cluster_text <- paste(input$exclusive_target_cluster_single, collapse = "_")
      table_name <- paste0("SingleDataset_ExclusiveMarkers_Cluster",
                           cluster_text, "_", format(Sys.time(), "%H%M%S"))
      description <- paste0("Exclusive markers for cluster(s): ",
                            paste(input$exclusive_target_cluster_single, collapse = ", "))
      parameters <- list(
        min_pct_target = input$exclusive_min_pct_target_single,
        max_pct_other = input$exclusive_max_pct_other_single,
        min_log2fc = input$exclusive_min_log2fc_single,
        detection_threshold = input$exclusive_detection_threshold_single,
        statistical_test = input$exclusive_statistical_test_single
      )
      result <- storeDETable(single_gene_table_storage(),
                             markers_copy_venn,
                             table_name,
                             description,
                             "exclusive_markers",
                             parameters)
      if (result$success) {
        single_gene_table_storage(result$storage)
        showNotification(result$message, type = "message")
      } else {
        showNotification(result$message, type = "warning")
      }
      if ("pvalue" %in% colnames(markers)) {
        markers$pvalue <- format_pvalue_robust(markers$pvalue)
      }
      if ("pvalue_adjusted" %in% colnames(markers)) {
        markers$pvalue_adjusted <- format_pvalue_robust(markers$pvalue_adjusted)
      }
      cleaned_gene_names <- clean_gene_names_for_html(markers$gene)
      markers$gene_link <- paste0('<a href="#" class="gene-name" data-gene="',
                                  cleaned_gene_names, '">',
                                  cleaned_gene_names, '</a>')
      if (input$exclusive_statistical_test_single != "none") {
        display_cols <- c("gene_link", "pct_detected_target", "pct_detected_other",
                          "detection_difference", "mean_expr_target", "mean_expr_other",
                          "log2_fold_change", "fold_change", "specificity_score",
                          "pvalue", "pvalue_adjusted")
        col_names <- c("Gene", "% Target", "% Others", "% Diff",
                       "Mean Target", "Mean Others", "Log2FC", "FC",
                       "Specificity", "P-value", "Adj. P-value")
      } else {
        display_cols <- c("gene_link", "pct_detected_target", "pct_detected_other",
                          "detection_difference", "mean_expr_target", "mean_expr_other",
                          "log2_fold_change", "fold_change", "specificity_score")
        col_names <- c("Gene", "% Target", "% Others", "% Diff",
                       "Mean Target", "Mean Others", "Log2FC", "FC", "Specificity")
      }
      markers_display <- markers[, display_cols]
      colnames(markers_display) <- col_names
      exclusive_markers_single(markers_display)
      showNotification(paste("Found", nrow(markers), "exclusive biomarkers!"), type = "message")
    }, error = function(e) {
      removeModal()
      showNotification(
        paste("Error finding exclusive biomarkers:", paste(e$message, collapse = " ")),
        type = "error",
        duration = 10
      )
      message("Error finding exclusive biomarkers: ", paste(e$message, collapse = " "))
    })
  })
  
  # Display results table
  output$table_exclusive_markers_single <- renderDataTable({
    req(exclusive_markers_single())
    df <- exclusive_markers_single()
    df <- df[!is.na(df$Gene) & df$Gene != "", ]
    
    spec_vals  <- as.numeric(df$Specificity)
    log2_vals  <- as.numeric(df$Log2FC)
    spec_range <- range(spec_vals, na.rm = TRUE, finite = TRUE)
    log2_range <- range(log2_vals, na.rm = TRUE, finite = TRUE)
    if (any(!is.finite(spec_range))) spec_range <- c(0, 1)
    if (any(!is.finite(log2_range))) log2_range <- c(0, 1)
    
    datatable(
      df,
      escape = FALSE,
      options = list(
        pageLength = 15,
        lengthMenu = c(15, 25, 50, 100, 200),
        order = list(list(8, 'desc')),
        dom = 'Blfrtip',
        scrollX = TRUE,
        columnDefs = list(
          list(className = 'dt-center', targets = 1:(ncol(df) - 1)),
          list(width = '120px', targets = 0)
        )
      ),
      class = 'cell-border stripe',
      rownames = FALSE
    ) %>%
      formatStyle(
        'Specificity',
        background = styleColorBar(spec_range, '#90EE90'),
        backgroundSize = '80% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      ) %>%
      formatStyle(
        'Log2FC',
        background = styleColorBar(log2_range, '#ADD8E6'),
        backgroundSize = '80% 90%',
        backgroundRepeat = 'no-repeat',
        backgroundPosition = 'center'
      )
  })
  
  # Download handler for exclusive markers
  output$download_exclusive_markers_single <- createDownloadHandler(
    reactive_data = reactive({
      markers <- exclusive_markers_single()
      req(markers)
      
      # Remove HTML links for download
      markers$Gene <- gsub('.*data-gene="([^"]+)".*', '\\1', markers$Gene)
      return(markers)
    }),
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset") 
    }),
    data_name = reactive({ 
      paste0("ExclusiveMarkers_Cluster", paste(input$exclusive_target_cluster_single, collapse = "_"))
    }),
    download_type = "csv"
  )
  
  #####################Exclusive Biomarkers Visualization###############################
  
  # Reactive value to store the plot
  exclusive_plot_single <- reactiveVal()
  
  # Generate plot when button is clicked
  observeEvent(input$generate_exclusive_plot_single, {
    req(single_dataset_object(), input$exclusive_genes_to_plot_single)
    
    tryCatch({
      # Parse genes from text input
      genes_input <- trimws(strsplit(input$exclusive_genes_to_plot_single, ",")[[1]])
      genes_input <- genes_input[genes_input != ""]
      
      if (length(genes_input) == 0) {
        showNotification("Please enter at least one gene name", type = "warning")
        return()
      }
      
      if (length(genes_input) > 10) {
        showNotification("Please enter maximum 10 genes for visualization", type = "warning")
        genes_input <- genes_input[1:10]
      }
      
      # Validate genes exist in dataset
      seurat_obj <- single_dataset_object()
      available_genes <- rownames(seurat_obj[[DefaultAssay(seurat_obj)]])
      missing_genes <- setdiff(genes_input, available_genes)
      
      if (length(missing_genes) > 0) {
        showNotification(
          paste("Genes not found in dataset:", paste(missing_genes, collapse = ", ")),
          type = "error",
          duration = 10
        )
        genes_input <- intersect(genes_input, available_genes)
        if (length(genes_input) == 0) {
          return()
        }
      }
      
      showModal(modalDialog(
        title = "Generating Plot",
        paste("Creating", input$exclusive_plot_type_single, "for", length(genes_input), "gene(s)..."),
        easyClose = FALSE,
        footer = NULL
      ))
      
      # Generate plot based on selected type
      plot_obj <- NULL
      
      if (input$exclusive_plot_type_single == "dotplot") {
        # Dot plot showing expression across clusters
        plot_obj <- DotPlot(
          seurat_obj,
          features = genes_input,
          cols = c("lightgrey", "blue"),
          dot.scale = 8
        ) +
          RotatedAxis() +
          labs(
            title = paste("Expression of Exclusive Biomarkers across Clusters"),
            subtitle = paste("Target cluster:", paste(input$exclusive_target_cluster_single, collapse = ", ")),
            x = "Genes",
            y = "Cluster"
          ) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            legend.position = "right"
          )
        
      } else if (input$exclusive_plot_type_single == "violin") {
        # Violin plot
        plot_obj <- VlnPlot(
          seurat_obj,
          features = genes_input,
          ncol = min(3, length(genes_input)),
          pt.size = 0.1
        ) +
          labs(
            title = paste("Expression Distribution of Exclusive Biomarkers"),
            subtitle = paste("Target cluster:", paste(input$exclusive_target_cluster_single, collapse = ", "))
          ) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30")
          )
        
      } else if (input$exclusive_plot_type_single == "feature") {
        # Feature plot on UMAP
        plot_obj <- FeaturePlot(
          seurat_obj,
          features = genes_input,
          ncol = min(3, length(genes_input)),
          reduction = "umap",
          pt.size = 0.5
        ) +
          labs(
            title = paste("Spatial Expression of Exclusive Biomarkers"),
            subtitle = paste("Target cluster:", paste(input$exclusive_target_cluster_single, collapse = ", "))
          ) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30")
          )
      }
      
      removeModal()
      
      if (!is.null(plot_obj)) {
        exclusive_plot_single(plot_obj)
        showNotification("Plot generated successfully!", type = "message")
      }
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error generating plot:", e$message), type = "error", duration = 10)
      print(paste("Detailed error:", e$message))
    })
  })
  
  # Render the plot
  output$plot_exclusive_markers_single <- renderPlot({
    req(exclusive_plot_single())
    print(exclusive_plot_single())
  }, height = 600)
  
  # Download handler for the plot
  output$download_exclusive_plot_single <- downloadHandler(
    filename = function() {
      object_name <- getObjectNameForDownload(single_dataset_object(), default_name = "ExclusiveMarkers")
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      genes_text <- gsub(",", "_", input$exclusive_genes_to_plot_single)
      genes_text <- gsub(" ", "", genes_text)
      genes_text <- substr(genes_text, 1, 50)  # Limit filename length
      
      paste0(object_name, "_ExclusiveMarkers_", genes_text, "_", 
             input$exclusive_plot_type_single, "_", timestamp, ".", 
             input$exclusive_plot_format_single)
    },
    content = function(file) {
      req(exclusive_plot_single())
      
      tryCatch({
        plot_obj <- exclusive_plot_single()
        
        # Determine plot dimensions based on number of genes
        genes_input <- trimws(strsplit(input$exclusive_genes_to_plot_single, ",")[[1]])
        genes_input <- genes_input[genes_input != ""]
        n_genes <- length(genes_input)
        
        # Adjust dimensions based on plot type
        if (input$exclusive_plot_type_single == "dotplot") {
          width <- max(8, n_genes * 1.2)
          height <- 6
        } else if (input$exclusive_plot_type_single == "violin") {
          ncol <- min(3, n_genes)
          nrow <- ceiling(n_genes / ncol)
          width <- ncol * 4
          height <- nrow * 4
        } else {  # feature
          ncol <- min(3, n_genes)
          nrow <- ceiling(n_genes / ncol)
          width <- ncol * 5
          height <- nrow * 4
        }
        
        # Save based on format
        if (input$exclusive_plot_format_single == "pdf") {
          pdf(file, width = width, height = height)
          print(plot_obj)
          dev.off()
        } else if (input$exclusive_plot_format_single == "svg") {
          svg(file, width = width, height = height)
          print(plot_obj)
          dev.off()
        } else {
          # For raster formats (tiff, png, jpeg)
          ggsave(
            filename = file,
            plot = plot_obj,
            width = width,
            height = height,
            dpi = input$exclusive_plot_dpi_single,
            device = input$exclusive_plot_format_single
          )
        }
        
        showNotification("Plot downloaded successfully!", type = "message")
        
      }, error = function(e) {
        showNotification(paste("Error downloading plot:", e$message), type = "error")
        print(paste("Download error:", e$message))
      })
    }
  )
  
  
  
######################################### Venn Diagramm ################################
  
  # Initialize reactive storage
  single_gene_table_storage <- reactiveVal(list())
  single_current_gene_lists <- reactiveVal(NULL)
  single_venn_plot_rendered <- reactiveVal(NULL)
  
  # Update select inputs
  observe({
    updateVennSelectInputs(session, single_gene_table_storage(), c("venn_table_1_single", "venn_table_2_single", "venn_table_3_single"))
  })
  
  # Enable/disable generate button
  observe({
    tables <- single_gene_table_storage()
    if (length(tables) < 2) {
      shinyjs::disable("generate_venn_btn_single")
    } else {
      shinyjs::enable("generate_venn_btn_single")
    }
  })
  
  # Generate Venn diagram
  observeEvent(input$generate_venn_btn_single, {
    showModal(modalDialog(title = "Generating Venn Diagram", "Processing...", easyClose = FALSE, footer = NULL))
    result <- processVennGeneration(
      table_storage = single_gene_table_storage(),
      selected_tables = c(input$venn_table_1_single, input$venn_table_2_single, input$venn_table_3_single),
      filter_params = list(
        significant_only = c(input$significant_only_venn_1_single, input$significant_only_venn_2_single, input$significant_only_venn_3_single),
        log_fc_threshold = c(input$log_fc_threshold_venn_1_single, input$log_fc_threshold_venn_2_single, input$log_fc_threshold_venn_3_single),
        p_val_threshold = input$p_val_threshold_venn_single,
        use_adjusted_p = input$use_adjusted_p_venn_single,
        direction = input$venn_direction_single
      ),
      colors = c(input$venn_color_1_single, input$venn_color_2_single, input$venn_color_3_single)
    )
    
    removeModal()
    
    if (result$success) {
      single_venn_plot_rendered(result$venn_plot)
      single_current_gene_lists(result$overlaps)
      updateSelectInput(session, "selected_gene_set_single", choices = names(result$overlaps))
      shinyjs::enable("download_venn_diagram_single")
    } else {
      showNotification(result$message, type = "error")
    }
  })
  
  # Render Venn diagram
  output$venn_plot_single <- renderPlot({
    req(single_venn_plot_rendered())
    grid.draw(single_venn_plot_rendered())
  })
  
  # Display gene table
  output$venn_gene_table_single <- renderDT({
    req(single_current_gene_lists(), input$selected_gene_set_single)
    overlaps <- single_current_gene_lists()
    selected_genes <- overlaps[[input$selected_gene_set_single]]
    if (length(selected_genes) == 0) {
      return(data.frame(Gene = character(0)))
    }
    gene_df <- data.frame(Gene = selected_genes)
    datatable(gene_df, options = list(pageLength = 15, scrollX = TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')), rownames = FALSE)
  })
  
  
  
  output$download_venn_diagram_single <- createDownloadHandler(
    reactive_data = single_venn_plot_rendered,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "VennDiagram") }),
    data_name = "venn_comparison",
    download_type = "plot",
    plot_params = list(
      file_type = input$venn_diagram_format_single,
      width = 8,
      height = 6,
      dpi = input$venn_diagram_dpi_single
    )
  )
  

  output$download_venn_gene_lists_single <- createDownloadHandler(
    reactive_data = single_current_gene_lists,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "VennDiagram") }),
    data_name = "gene_lists",
    download_type = "ods"
  )
  
##########################CLuster Composition#########################
  
  # Add this reactive value with other reactive values in single dataset server
  cluster_composition_single <- reactiveVal(NULL)
  
  # Observer for generating cluster composition table for single dataset
  observeEvent(input$generate_cluster_composition_single, {
    req(single_dataset_object())
    
    tryCatch({
      cluster_composition <- create_cluster_composition_table(
        single_dataset_object(),
        is_integrated = FALSE
      )
      cluster_composition_single(cluster_composition)
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  output$cluster_composition_single <- renderDT({
    req(cluster_composition_single())
    render_cluster_composition_table(cluster_composition_single(), is_integrated = FALSE)
  })
  
  # Download handler for single dataset cluster composition
  output$download_cluster_composition_single <- createDownloadHandler(
    reactive_data = reactive({ 
      data <- cluster_composition_single()
      data$Size_Bar <- NULL  
      return(data)
    }),
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "ClusterComposition") }),
    data_name = "cluster_composition",
    download_type = "csv"
  )
  
  
  #####################Co-expression###############################
  # Reactive values for storing results
  gene_coexpression_data_single <- reactiveVal(NULL)
  coexpression_plot_single <- reactiveVal(NULL)
  
  # Analyze co-expression when button is clicked
  observeEvent(input$analyze_coexpression_single, {
    req(single_dataset_object(), input$gene_text_coexpression_single)
    
    tryCatch({
      # Parse and validate genes
      genes_input <- trimws(strsplit(input$gene_text_coexpression_single, ",")[[1]])
      genes_input <- genes_input[genes_input != ""]
      
      if (length(genes_input) != 2) {
        showNotification("Please enter exactly 2 gene names for co-expression analysis", type = "error")
        return()
      }
      
      # Create thresholds vector
      thresholds <- c(input$gene1_threshold_single, input$gene2_threshold_single)
      names(thresholds) <- genes_input
      
      # Validate genes exist
      available_genes <- rownames(single_dataset_object()[[DefaultAssay(single_dataset_object())]])
      missing_genes <- setdiff(genes_input, available_genes)
      if (length(missing_genes) > 0) {
        showNotification(paste("Genes not found:", paste(missing_genes, collapse = ", ")), type = "error")
        return()
      }
      
      # Use modular function from cells_genes_expressions_newarch.R
      coexpr_results <- analyze_gene_coexpression(
        seurat_obj = single_dataset_object(), 
        genes = genes_input, 
        assay_name = DefaultAssay(single_dataset_object()), 
        expression_thresholds = thresholds, 
        is_integrated = FALSE
      )
      
      # Store results and create plot using modular functions
      gene_coexpression_data_single(coexpr_results)
      coexpression_plot_single(create_coexpression_plot(coexpr_results$data, coexpr_results$genes_analyzed))
      
      showNotification("Co-expression analysis completed!", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Render results using modular functions
  output$gene_coexpression_table_single <- renderDT({
    req(gene_coexpression_data_single())
    render_coexpression_table(gene_coexpression_data_single(), "gene_coexpression_table_single")
  })
  
  output$gene_coexpression_plot_single <- renderPlot({
    req(coexpression_plot_single())
    coexpression_plot_single()
  }, height = 600)
  
  # Download handlers for coexpression analysis
  output$download_coexpression_table_single <- createDownloadHandler(
    reactive_data = gene_coexpression_data_single,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "Coexpression") 
    }),
    data_name = "coexpression_analysis",
    download_type = "csv"
  )
  
  output$download_coexpression_plot_single <- createDownloadHandler(
    reactive_data = coexpression_plot_single,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "Coexpression") 
    }),
    data_name = "coexpression_plot",
    download_type = "plot",
    plot_params = list(
      file_type = reactive(input$format_coexpr_single),
      width = 10,
      height = 8,
      dpi = reactive(input$dpi_coexpr_single)
    )
  )
  
  
  ########################## Volcano Plot ##########################
  
  # Reactive value for storing volcano plot
  volcano_plot_single <- reactiveVal(NULL)
  
  # Update volcano plot source selector when DEG results are available
  observe({
    tables <- single_gene_table_storage()
    if (length(tables) > 0) {
      table_choices <- c("None" = "", names(tables))
      updateSelectInput(session, "volcano_deg_source_single", choices = table_choices)
    } else {
      updateSelectInput(session, "volcano_deg_source_single", choices = c("None" = ""))
    }
  })
  
  # Generate volcano plot
  observeEvent(input$generate_volcano_btn_single, {
    req(input$volcano_deg_source_single, input$volcano_deg_source_single != "")
    
    tryCatch({
      showModal(modalDialog(
        title = "Generating Volcano Plot", 
        "Processing...", 
        easyClose = FALSE, 
        footer = NULL
      ))
      
      # Get selected DEG results from storage
      tables <- single_gene_table_storage()
      table_entry <- tables[[input$volcano_deg_source_single]]
      
      # Extract the actual data from the table entry
      deg_data <- table_entry$data
      
      if (is.null(deg_data) || nrow(deg_data) == 0) {
        removeModal()
        showNotification("Selected DEG table is empty or not found", type = "error")
        return()
      }
      
      # Generate volcano plot using function from volcano_plot_functions_newarch.R
      volcano <- generateVolcanoPlot(
        deg_results = deg_data,
        log2fc_threshold = input$volcano_log2fc_threshold_single,
        pval_threshold = input$volcano_pval_threshold_single,
        color_up = input$volcano_color_up_single,
        color_down = input$volcano_color_down_single,
        color_ns = input$volcano_color_ns_single,
        label_top_genes = input$volcano_label_genes_single,
        point_size = input$volcano_point_size_single,
        point_alpha = input$volcano_point_alpha_single
      )
      
      volcano_plot_single(volcano)
      removeModal()
      showNotification("Volcano plot generated successfully!", type = "message")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error generating volcano plot:", e$message), type = "error")
      message("DEBUG volcano error: ", e$message)
      print(traceback())
    })
  })
  
  # Render volcano plot
  output$volcano_plot_single <- renderPlot({
    req(volcano_plot_single())
    volcano_plot_single()
  })
  
  # Download volcano plot
  output$download_volcano_plot_single <- createDownloadHandler(
    reactive_data = volcano_plot_single,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "VolcanoPlot") 
    }),
    data_name = "volcano_plot",
    download_type = "plot",
    plot_params = list(
      file_type = input$volcano_plot_format_single,
      width = 10,
      height = 8,
      dpi = input$volcano_plot_dpi_single
    )
  )
  
  
  
  
  ############################## Subseting ##############################

  # Variables wich stock the subset of single_dataset_object()
  subset_seurat <- reactiveVal()

  observe({
    { subset_seurat(single_dataset_object())
    }
  })

  observe({
    if (!is.null(single_dataset_object())) {
      updateSelectInput(session, "select_ident_subset", choices = unique(Idents(single_dataset_object())))
    }
  })

  #Run subset
  observeEvent(input$apply_subset, {
    req(single_dataset_object())
    tryCatch({
      subset_seurat_temp <- single_dataset_object()

      if (length(input$select_ident_subset) > 0) {
        subset_seurat_temp <- subset(x = subset_seurat_temp, idents = input$select_ident_subset)

        if (nrow(subset_seurat_temp@meta.data) == 0) {
          showNotification("No cells found with the selected identities.", type = "error")
          return()
        }
      }

      subset_seurat(subset_seurat_temp)



    }, error = function(e) {
      showNotification(paste("Error applying subset: ", e$message), type = "error")
    })
  })


  observeEvent(input$apply_gene_subset, {
    tryCatch({
      req(single_dataset_object())
      
      gene_list <- trimws(unlist(strsplit(input$gene_list_subset, ",")))
      gene_list <- gene_list[gene_list != ""]
      
      if (length(gene_list) == 0) {
        showNotification("Please enter valid gene names.", type = "error")
        return()
      }
      
      message("=== GENE SUBSET ANALYSIS ===")
      message("Genes requested: ", paste(gene_list, collapse = ", "))
      message("Active assay: ", DefaultAssay(single_dataset_object()))
      message("Expression threshold: ", input$expression_threshold)
      message("Negative selection: ", input$negative_gene_subset)
      
      available_genes <- rownames(single_dataset_object())
      missing_genes   <- setdiff(gene_list, available_genes)
      
      if (length(missing_genes) > 0) {
        showNotification(
          paste("Genes not found:", paste(missing_genes, collapse = ", ")),
          type = "error"
        )
        return()
      }
      
      # Fetch normalized expression data using Seurat v5 API
      expression_data <- FetchData(
        single_dataset_object(),
        vars  = gene_list,
        layer = "data"
      )
      
      for (gene in gene_list) {
        gene_expr         <- expression_data[[gene]]
        n_cells_total     <- length(gene_expr)
        n_cells_expressed <- sum(gene_expr > 0)
        n_cells_above     <- sum(gene_expr >= input$expression_threshold)
        max_expr          <- max(gene_expr)
        mean_expr         <- mean(gene_expr[gene_expr > 0])
        
        message(sprintf(
          "Gene %s: %d/%d cells express it (%.1f%%), %d above threshold (%.1f%%), max=%.2f, mean(expr>0)=%.2f",
          gene,
          n_cells_expressed, n_cells_total, 100 * n_cells_expressed / n_cells_total,
          n_cells_above, 100 * n_cells_above / n_cells_total,
          max_expr, mean_expr
        ))
        
        if (n_cells_expressed == 0) {
          showNotification(paste("WARNING:", gene, "is not expressed in any cells!"), type = "warning", duration = 5)
        } else if (n_cells_above == 0) {
          showNotification(paste("WARNING:", gene, "is expressed but no cells above threshold", input$expression_threshold), type = "warning", duration = 5)
        }
      }
      
      expression_matrix <- sapply(gene_list, function(gene) {
        expression_data[[gene]] >= input$expression_threshold
      })
      
      n_genes_per_cell <- rowSums(expression_matrix)
      
      if (isTRUE(input$negative_gene_subset)) {
        cells_to_keep <- colnames(single_dataset_object())[
          n_genes_per_cell < input$num_genes_to_express
        ]
        message(sprintf("Negative selection: keeping %d cells NOT expressing required genes", length(cells_to_keep)))
      } else {
        cells_to_keep <- colnames(single_dataset_object())[
          n_genes_per_cell >= input$num_genes_to_express
        ]
        message(sprintf("Positive selection: keeping %d cells expressing required genes", length(cells_to_keep)))
      }
      
      message(sprintf(
        "Cells kept: %d/%d (%.1f%%)",
        length(cells_to_keep),
        ncol(single_dataset_object()),
        100 * length(cells_to_keep) / ncol(single_dataset_object())
      ))
      
      if (length(cells_to_keep) == 0) {
        showNotification(
          "No cells found matching criteria. Try adjusting the expression threshold or number of genes required.",
          type = "error",
          duration = 10
        )
        return()
      }
      
      subset_seurat_temp <- subset(single_dataset_object(), cells = cells_to_keep)
      subset_seurat(subset_seurat_temp)
      
      showNotification(
        paste("Subset created:", length(cells_to_keep), "cells retained"),
        type = "message"
      )
      
    }, error = function(e) {
      message("ERROR in gene subset: ", conditionMessage(e))
      showNotification(paste("Error:", conditionMessage(e)[1]), type = "error")
    })
  })
  

  # Downloading Seurat subset object
  output$download_subset_seurat <- createDownloadHandler(
    reactive_data = subset_seurat,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "SeuratSubset") }),
    data_name = "seurat_subset",
    download_type = "seurat",
    show_modal = TRUE
  )
  
  # Umap plot with all clusters
  reactivePlotAll <- reactive({
    req(single_dataset_object())
    plot_data_all <- DimPlot(single_dataset_object(), group.by = "ident", label = TRUE) +
      theme(axis.line = element_line(size = 0.5)) +
      NoLegend()+ggtitle(NULL)
    return(plot_data_all)
  })

  # Umap plot with filtered data
  reactivePlotSubset <- reactive({
    req(subset_seurat())
    plot_data_subset <- DimPlot(subset_seurat(), group.by = "ident", label = TRUE) +
      theme(axis.line = element_line(size = 0.5)) +
      NoLegend()+ggtitle(NULL)
    return(plot_data_subset)
  })

  # Rendering Umap plot with all clusters
  output$global_umap <- renderPlot({
    plot_data_all <- reactivePlotAll()
    print(plot_data_all)
  })

  # Rendering Umap plot with filtered data
  output$subset_umap <- renderPlot({
    plot_data_subset <- reactivePlotSubset()
    print(plot_data_subset)
  })

  # Download analysis parameters as CSV
  output$download_params_single <- downloadHandler(
    filename = function() {
      obj_name <- getObjectNameForDownload(single_dataset_object(), default_name = "object")
      paste0(obj_name, "_analysis_parameters_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      params_df <- collect_analysis_params(
        input      = input,
        module     = "single",
        qc_stats   = qc_stats_single(),
        seurat_obj = single_dataset_object()
      )
      write.csv(params_df, file, row.names = FALSE)
    }
  )
  
  output$download_report_single <- downloadHandler(
    filename = function() {
      obj_name <- getObjectNameForDownload(single_dataset_object(), default_name = "object")
      paste0(obj_name, "_report_", format(Sys.Date(), "%Y%m%d"), ".pdf")
    },
    content = function(file) {
      # Resolve script path locally — no external variable dependency
      script_candidates <- c(file.path(getwd(), "generate_report.py"), "/app/generate_report.py")
      script_path <- script_candidates[file.exists(script_candidates)][1]
      if (is.na(script_path)) stop("generate_report.py not found in project root or /app/")
      
      csv_tmp <- tempfile(fileext = ".csv")
      on.exit(unlink(csv_tmp), add = TRUE)
      params_df <- collect_analysis_params(
        input      = input,
        module     = "single",
        qc_stats   = qc_stats_single(),
        seurat_obj = single_dataset_object()
      )
      write.csv(params_df, csv_tmp, row.names = FALSE)
      result <- system2(
        "python3",
        args   = c(shQuote(script_path),
                   "--csv",    shQuote(csv_tmp),
                   "--output", shQuote(file),
                   "--module", "single"),
        stdout = TRUE,
        stderr = TRUE
      )
      if (!file.exists(file)) {
        stop("PDF generation failed:\n", paste(result, collapse = "\n"))
      }
    }
  )
  
  
  
}

