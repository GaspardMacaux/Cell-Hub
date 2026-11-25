############################## Script Single Dataset Analysis ##############################

# single_dataset_server.R

single_dataset_server <- function(input, output, session) {

  ############################## Shared reactive state ##############################
  # Objects used by several sections below

  # Reactives objects for the single dataset part
  single_dataset_object <- reactiveVal(NULL)

  # Register in the session-wide registry so other modules (CellChat, trajectory) can import it
  session$userData$single_dataset_object <- single_dataset_object

  gene_list <- reactiveValues(features = NULL)

  # Analysis parameter tracking for CSV/PDF export
  qc_stats_single <- reactiveVal(NULL)

  # ✅ Create log system for SINGLE module
  log_system <- createLogSystem()

  clustering_plot <- reactiveVal()

  # FeaturePlot avec options 2D/3D
  feature_plot <- reactiveVal()

  # DotPlot: same input checks; create_dotplot() also joins layers (rule 1)
  dot_plot <- reactiveVal()

  # VlnPlot: empty field, unknown genes and empty cluster selections are reported instead of failing silently
  vln_plot <- reactiveVal()

  forced_col_single <- reactiveVal(NULL)

  # Simple history system for undo functionality
  cluster_history_single <- reactiveValues(
    states = list(),
    max_history = 10
  )

  # Initialize reactive storage
  single_gene_table_storage <- reactiveVal(list())

  ############################## Download buttons state ##############################
  # Each download button stays disabled until the result it exports exists
  bind_download_state("downloadUMAP", function() req(clustering_plot()))
  bind_download_state("download_doublet_results", function() req(doublet_results_data()))
  bind_download_state("download_ridge_plot", function() ridge_plot())
  bind_download_state("downloadDotPlot", function() dot_plot())
  bind_download_state("downloadVlnPlot", function() vln_plot())
  bind_download_state("downloadFeaturePlot", function() feature_plot())
  bind_download_state("download_genes_number_expresion", function() number_of_nuclei())
  bind_download_state("download_heatmap_markers_single", function() { markers_df <- heatmap_markers_data_single(); req(markers_df) })
  bind_download_state("download_heatmap_single", function() heatmap_plot())
  bind_download_state("download_scatter_single", function() scatter_plot_single())
  bind_download_state("download_transfer_clusters_single", function() req(transfer_result_single()))
  bind_download_state("download_transfer_cells_single", function() req(transfer_result_single()))
  bind_download_state("download_singler_clusters_single", function() req(singler_result_single()))
  bind_download_state("download_singler_cells_single", function() req(singler_result_single()))
  bind_download_state("downloadUMAPCluster", function() umap_plot())
  bind_download_state("download_markers_global_single", function() { df <- markers_global_single(); req(df) })
  bind_download_state("download_markers_pairwise_single", function() { df <- markers_pairwise_single(); req(df) })
  bind_download_state("download_enrichment_ora_plot_single", function() { r <- enrichment_results_single(); req(!is.null(r$ora)) })
  bind_download_state("download_enrichment_gsea_plot_single", function() { r <- enrichment_results_single(); req(!is.null(r$gsea)) })
  bind_download_state("download_enrichment_ora_table_single", function() { r <- enrichment_results_single(); req(r$ora) })
  bind_download_state("download_enrichment_gsea_table_single", function() { r <- enrichment_results_single(); req(r$gsea) })
  bind_download_state("download_exclusive_markers_single", function() { markers <- exclusive_markers_single(); req(markers) })
  bind_download_state("download_exclusive_plot_single", function() req(exclusive_plot_single()))
  bind_download_state("download_venn_diagram_single", function() single_venn_plot_rendered())
  bind_download_state("download_venn_gene_lists_single", function() single_current_gene_lists())
  bind_download_state("download_cluster_composition_single", function() cluster_composition_single())
  bind_download_state("download_coexpression_table_single", function() { res <- gene_coexpression_data_single(); req(res) })
  bind_download_state("download_coexpression_plot_single", function() { res <- gene_coexpression_data_single(); req(res, length(res$genes_analyzed) >= 2) })
  bind_download_state("download_coexpr_de_single", function() req(coexpr_de_results_single()))
  bind_download_state("download_volcano_plot_single", function() volcano_plot_single())
  bind_download_state("download_subset_seurat", function() subset_seurat())
  bind_download_state("download_session_single", function() req(single_dataset_object()))
  
  ############################## Loading Data ##############################

  # ---- Helpers ----

  # Function to update UI elements after data loading or dataset changes
  updateUIElements <- function() {
    req(single_dataset_object())
    obj <- single_dataset_object()
    cluster_col <- if ("ClusterIdents" %in% colnames(obj@meta.data)) "ClusterIdents" else if ("seurat_clusters" %in% colnames(obj@meta.data)) "seurat_clusters" else NULL
    cluster_choices <- if (!is.null(cluster_col)) unique(as.character(obj@meta.data[[cluster_col]])) else as.character(unique(Idents(obj)))
    updateSelectInput(session, "cluster_order_vln", choices = cluster_choices)
    updateSelectInput(session, "cluster_order_dot", choices = cluster_choices)
    # Gene names read from the assay itself: LayerData(layer = "counts") fails on SCT-only objects and on split counts layers
    gene_assay <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
    gene_list <- sort(rownames(obj[[gene_assay]]))
    updatePickerInput(session, "gene_select", choices = gene_list)
    updatePickerInput(session, "gene_select_genes_analysis", choices = gene_list)
    return(invisible(NULL))
  }

  # ---- Observers and computations ----

  # Tab 1 : Loading Data

  shinyjs::useShinyjs() # button deactivation

  # File loading with logs (raw data: ZIP, H5, H5AD)
  observeEvent(input$file, run_safely("file", {
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
      if (isTRUE(input$optimize_memory_single %in% c("slim", "aggressive"))) {
        log_system$add_log("Applying memory optimization...", "info")
        seuratObj <- slimSeuratObject(seuratObj, level = input$optimize_memory_single, log_function = log_system$add_log)
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
  }))

  # Load processed Seurat object (RDS or H5AD)
  observeEvent(input$load_seurat_file, run_safely("load seurat file", {
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
      if (isTRUE(input$optimize_memory_single %in% c("slim", "aggressive"))) {
        log_system$add_log("Applying memory optimization...", "info")
        loaded_seurat <- slimSeuratObject(loaded_seurat, level = input$optimize_memory_single, log_function = log_system$add_log)
      }
      logMemoryUsage(log_system$add_log, loaded_seurat, "Loaded object")
      # Restore the saved session (stored tables + QC statistics); without one, clear the previous
      # dataset's tables so Venn diagrams never mix results from two different objects
      session_state <- unpack_session_state(loaded_seurat)
      single_gene_table_storage(session_state$gene_tables)
      qc_stats_single(session_state$qc_stats)
      if (session_state$found) {
        log_system$add_log(sprintf("Session restored: %d stored table(s)%s",
                                   length(session_state$gene_tables),
                                   if (!is.null(session_state$qc_stats)) ", QC statistics" else ""), "success")
      }
      if (!is.null(get_umap_reduction(loaded_seurat))) {
        log_system$add_log("Generating initial UMAP plot...", "info")
        plot_result <- generateInitialPlot(
          seurat_obj    = loaded_seurat,
          remove_axes   = FALSE,
          remove_legend = FALSE,
          cluster_col   = if ("ClusterIdents" %in% colnames(loaded_seurat@meta.data)) "ClusterIdents"
          else if ("seurat_clusters" %in% colnames(loaded_seurat@meta.data)) "seurat_clusters"
          else NULL
        )
        if (!is.null(plot_result)) {
          loaded_seurat <- plot_result$seurat_obj
          clustering_plot(plot_result$plot)
          log_system$add_log("UMAP plot generated successfully", "success")
        }
      } else {
        log_system$add_log("No UMAP reduction found in object", "warning")
      }
      # Assign once after all modifications — triggers reactive chain only once
      single_dataset_object(loaded_seurat)
      log_system$add_log("=== OBJECT LOADED SUCCESSFULLY ===", "success")
      showNotification(
        if (session_state$found) sprintf("Object loaded successfully! Session restored (%d stored table(s)).", length(session_state$gene_tables))
        else "Object loaded successfully!",
        type = "message"
      )
    }, error = function(e) {
      log_system$add_log(paste("ERROR:", conditionMessage(e)[1]), "error")
      log_system$add_log("Loading failed. Check file format and try again.", "error")
      showNotification(paste("Error loading object:", conditionMessage(e)[1]), type = "error")
    })
  }))

  # Fonction réactive pour récupérer les clusters depuis l'objet Seurat
  get_clusters <- reactive({
    req(single_dataset_object())
    return(getClusters(single_dataset_object()))
  })

  # ---- UI and text outputs ----

  # ✅ Render logs output
  output$single_loading_logs <- renderPrint({
    cat(log_system$log_value())
  })

  ############################## QC metrics and normalization ##############################

  # ---- Helpers ----

  # Helper function for gene detection with different thresholds
  get_gene_count_with_threshold <- function(seurat_obj, min_cells = 3, min_counts = 0) {
    counts_matrix <- LayerData(seurat_obj, assay = "RNA", layer = "counts")
    genes_detected <- rowSums(counts_matrix >= min_counts) >= min_cells
    gene_count <- sum(genes_detected)

    message(sprintf("Genes detected (≥%d cells, ≥%d counts): %d", min_cells, min_counts, gene_count))
    return(gene_count)
  }

  # ---- Observers and computations ----

  # Count nuclei before and after QC (UNCHANGED)
  nuclei_count <- reactive({
    req(single_dataset_object())
    meta <- single_dataset_object()@meta.data
    sum(qc_pass_cells(meta, input$nFeature_range[1], input$nFeature_range[2], input$percent.mt_max))
  })

  # Genes detected in at least one cell (union over split counts layers); computed once per object change
  detected_genes_count <- reactive({
    req(single_dataset_object())
    obj <- single_dataset_object()
    assay <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)
    layers <- SeuratObject::Layers(obj[[assay]], search = "counts")
    if (length(layers) == 0) return(nrow(obj[[assay]]))
    detected <- unique(unlist(lapply(layers, function(l) {
      m <- LayerData(obj, assay = assay, layer = l)
      rownames(m)[Matrix::rowSums(m > 0) > 0]
    })))
    return(length(detected))
  })

  # GENES PER NUCLEUS (median/mean - AFTER QC filtering)
  nucleus_gene_metrics <- reactive({
    req(single_dataset_object())

    # Apply QC filters for per-nucleus metrics
    seuratRNA_subset <- subset(single_dataset_object(), cells = colnames(single_dataset_object())[
      qc_pass_cells(single_dataset_object()@meta.data, input$nFeature_range[1], input$nFeature_range[2], input$percent.mt_max)])

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

      # Apply filters (qc_pass_cells: genes and mitochondrial %, strict inequalities)
      thresholds <- list(min_features = input$nFeature_range[1], max_features = input$nFeature_range[2], max_mt_percent = input$percent.mt_max)
      keep <- qc_pass_cells(seurat_obj@meta.data, thresholds$min_features, thresholds$max_features, thresholds$max_mt_percent)
      if (!any(keep)) stop("No nuclei pass the current QC thresholds")
      seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[keep])

      single_dataset_object(seurat_obj)

      # Store pre/post QC stats and the thresholds applied for export (saved with the session)
      qc_stats_single(list(
        n_before             = pre_n,
        n_after              = ncol(seurat_obj),
        median_genes_before  = pre_genes,
        median_genes_after   = round(median(seurat_obj$nFeature_RNA, na.rm = TRUE)),
        median_counts_before = pre_counts,
        median_counts_after  = round(median(seurat_obj$nCount_RNA,   na.rm = TRUE)),
        median_mt_before     = pre_mt,
        median_mt_after      = round(median(seurat_obj$percent.mt,   na.rm = TRUE), 2),
        thresholds           = thresholds
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
        # Standard normalization workflow, always on the RNA assay: after an earlier SCTransform the default assay is SCT,
        # and NormalizeData / FindVariableFeatures would otherwise normalize the SCT assay
        normalized_seurat <- single_dataset_object()
        DefaultAssay(normalized_seurat) <- "RNA"
        normalized_seurat <- NormalizeData(
          normalized_seurat,
          assay = "RNA",
          normalization.method = input$normalization_method_single,
          scale.factor = input$scale_factor
        )
        print(paste0("Data normalization completed using method: ", input$normalization_method_single))

        # Identification of variable features
        normalized_seurat <- FindVariableFeatures(
          normalized_seurat,
          assay = "RNA",
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

  # ---- UI and text outputs ----

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

  # ---- Plots ----

  # QC violin plot rendered automatically; violins use every cell, optional jitter points are capped for drawing speed
  output$vlnplot <- renderPlot({
    req(single_dataset_object())
    tryCatch({
      seurat_obj <- single_dataset_object()
      required_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
      missing_metrics <- setdiff(required_metrics, colnames(seurat_obj@meta.data))
      if (length(missing_metrics) > 0) stop(paste("Missing metrics:", paste(missing_metrics, collapse = ", ")))
      meta_data <- as.data.frame(seurat_obj@meta.data[, required_metrics, drop = FALSE])
      meta_data$cell <- rownames(meta_data)
      meta_data$cluster <- tryCatch(as.character(Idents(seurat_obj)), error = function(e) rep("All", nrow(meta_data)))
      # Numeric clusters in numeric order, then named clusters alphabetically
      cluster_levels <- unique(meta_data$cluster)
      is_num <- !is.na(suppressWarnings(as.numeric(cluster_levels)))
      ordered_levels <- c(as.character(sort(as.numeric(cluster_levels[is_num]))), sort(cluster_levels[!is_num]))
      meta_data$cluster <- factor(meta_data$cluster, levels = ordered_levels)
      plot_data <- tidyr::pivot_longer(meta_data, cols = all_of(required_metrics), names_to = "metric", values_to = "value")
      plot_data$metric <- factor(plot_data$metric, levels = required_metrics)
      colors <- scales::hue_pal()(length(ordered_levels))
      p <- ggplot(plot_data, aes(x = cluster, y = value, fill = cluster)) +
        geom_violin(scale = "width", trim = TRUE)
      caption_text <- NULL
      # Points shown by default (checkbox absent = TRUE); at most 20,000 cells drawn (random subset, seed 42)
      if (!isFALSE(input$qc_show_points_single)) {
        max_points <- 20000
        n_cells <- nrow(meta_data)
        subsampled <- n_cells > max_points
        point_cells <- if (subsampled) { set.seed(42); sample(meta_data$cell, max_points) } else meta_data$cell
        point_data <- plot_data[plot_data$cell %in% point_cells, , drop = FALSE]
        p <- p + geom_jitter(data = point_data, height = 0, width = 0.3, size = 0.1, alpha = 0.4)
        if (subsampled) {
          caption_text <- sprintf("Points: random subset of %s of %s cells (seed 42). Violins use all cells.",
                                  format(max_points, big.mark = ","), format(n_cells, big.mark = ","))
        }
      }
      p +
        facet_wrap(~metric, scales = "free_y", ncol = 3) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
              strip.background = element_blank(), strip.text = element_text(size = 11, face = "bold"),
              plot.caption = element_text(size = 10, color = "grey40", hjust = 0)) +
        labs(x = NULL, y = NULL, caption = caption_text) +
        scale_fill_manual(values = colors)
    }, error = function(e) {
      message("QC violin plot error: ", conditionMessage(e)[1])
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:\n", conditionMessage(e)[1]), color = "red", size = 5) + theme_void()
    })
  })

  # Display scatter plots with verification; on error the message is drawn in the plot area (renderPlot must return a plot)
  output$scatter_plot1 <- renderPlot({
    req(input$show_plots, single_dataset_object())
    tryCatch({
      seuratRNA_subset <- subset(single_dataset_object(), cells = colnames(single_dataset_object())[
        qc_pass_cells(single_dataset_object()@meta.data, input$nFeature_range[1], input$nFeature_range[2], input$percent.mt_max)])
      if (ncol(seuratRNA_subset) == 0) stop("No cells pass the current filters!")
      FeatureScatter(seuratRNA_subset, feature1 = "nCount_RNA", feature2 = "percent.mt")
    }, error = function(e) {
      showNotification(paste0("Error in scatter plot 1: ", e$message), type = "error")
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:\n", conditionMessage(e)[1]), color = "red", size = 5) + theme_void()
    })
  })

  # Display scatter plot 2 with same verifications
  output$scatter_plot2 <- renderPlot({
    req(input$show_plots, single_dataset_object())
    tryCatch({
      seuratRNA_subset <- subset(single_dataset_object(), cells = colnames(single_dataset_object())[
        qc_pass_cells(single_dataset_object()@meta.data, input$nFeature_range[1], input$nFeature_range[2], input$percent.mt_max)])
      if (ncol(seuratRNA_subset) == 0) stop("No cells pass the current filters!")
      FeatureScatter(seuratRNA_subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    }, error = function(e) {
      showNotification(paste0("Error in scatter plot 2: ", e$message), type = "error")
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:\n", conditionMessage(e)[1]), color = "red", size = 5) + theme_void()
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

  # ---- Observers and computations ----

  #  Tab 3: Scaling and PCA and elbow plot

  # Scaling data
  observeEvent(input$scale_button, run_safely("scale button", {
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
      if (current_assay == "SCT") {
        print("SCTransform detected - data already scaled, skipping ScaleData")
        # SCTransform already scaled the data, just run PCA
        # Make sure we're using the SCT assay
        DefaultAssay(seurat_obj) <- "SCT"

      } else {
        # Only the variable features enter the PCA: scaling every gene would build a dense genes x cells matrix (~24 GB for 100k cells)
        print(paste("Scaling", length(VariableFeatures(seurat_obj)), "variable features in RNA assay"))
        seurat_obj <- ScaleData(seurat_obj, assay = "RNA", features = VariableFeatures(seurat_obj))
        print("Data scaling completed")
      }
      print(paste0("Running PCA on assay: ", DefaultAssay(seurat_obj)))

      # Seed-42 wrapper (dimensional_reduction_functions.R), same 50 PCs as the RunPCA default
      pca_result <- runPCA_reproducible(seurat_obj, npcs = 50, seed = 42, features = VariableFeatures(seurat_obj), verbose = FALSE)

      # Memory optimization: drop scale.data after PCA (no longer needed, recomputed for heatmaps)
      if (isTRUE(input$optimize_memory_single %in% c("slim", "aggressive"))) {
        log_system$add_log("Applying post-PCA memory optimization...", "info")
        pca_result <- slimSeuratObject(pca_result, level = input$optimize_memory_single, keep_counts = TRUE,
                                       keep_scale_data = FALSE,
                                       log_function = log_system$add_log)
      }
      single_dataset_object(pca_result)
      logMemoryUsage(log_system$add_log, pca_result, "Post-PCA object")
      print("PCA completed successfully")

      removeModal()
      showNotification("PCA completed successfully.", type = "message")

    }, error = function(e) {
      removeModal()
      print(paste("Error during scaling/PCA:", e$message))
      showNotification(paste("Error during scaling/PCA:", e$message), type = "error", duration = 10)
    })
  }))

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

  # ---- UI and text outputs ----

  #Printing PCA results
  output$pca_results <- renderPrint({
    req(!is.null(single_dataset_object()))
    if ("pca" %in% names(single_dataset_object())) {
      print(single_dataset_object()$pca, dims = 1:5, nfeatures = 5)
    }
  })

  # ---- Plots ----

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

  ############################## Neighbors calculation and clustering ##############################

  # ---- Observers and computations ----

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
      # Raw plot stored without labels or display options: styling is applied at render time
      clustering_plot(DimPlot(single_dataset_object(), group.by = "ident", label = FALSE) + ggtitle(NULL))
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
      showModal(modalDialog(title = "Please Wait", "Clustering process started...", easyClose = FALSE, footer = NULL))
      seurat_tmp <- findClusters_reproducible(object = single_dataset_object(), resolution = input$resolution_step1,
                                              algorithm = as.integer(input$algorithm_select_single), seed = 42, verbose = TRUE)
      # Cluster selectors (markers, heatmap, renaming, final UMAP) switch to the new clusters; otherwise they keep the column
      # selected before clustering (orig.ident on raw data, or an outdated ClusterIdents after re-clustering)
      forced_col_single("seurat_clusters")
      single_dataset_object(seurat_tmp)
      # Raw plot stored without labels or display options: styling is applied at render time
      clustering_plot(DimPlot(single_dataset_object(), group.by = "ident", label = FALSE) + ggtitle(NULL))
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Clustering error: ", e$message), type = "error")
    })
  })

  # Display options applied on top of the stored raw plot, shared by the screen and the download
  clustering_plot_styled <- reactive({
    req(clustering_plot())
    tryCatch(
      style_cluster_umap(clustering_plot(), label_mode = input$cluster_label_mode_single %||% "legend",
                         label_size = input$cluster_label_size_single %||% 5, bold_labels = isTRUE(input$cluster_bold_labels_single),
                         dark_mode = isTRUE(input$cluster_dark_mode_single), remove_axes = isTRUE(input$remove_axes)),
      error = function(e) {
        message("Clustering plot styling error: ", conditionMessage(e)[1])
        showNotification(paste("Display options could not be applied:", conditionMessage(e)[1]), type = "error", duration = 8)
        clustering_plot()
      })
  })

  # ---- Plots ----

  output$clustering_plot <- renderPlot({
    clustering_plot_styled()
  })

  # ---- Downloads ----

  #SAVING UMAP PLOT
  output$downloadUMAP <- createDownloadHandler(
    reactive_data = clustering_plot_styled,
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

  # ---- Reactive state ----

  # Reactive to store doublet results
  doublet_results_data <- reactiveVal(NULL)

  # ---- Observers and computations ----

  observeEvent(input$run_doubletfinder, run_safely("run doubletfinder", {
    req(single_dataset_object())
    showModal(modalDialog(title = "Running DoubletFinder", p("Processing your data..."), footer = NULL, easyClose = FALSE))
    tryCatch({
      seu <- single_dataset_object()
      # Validate prerequisites up front so failures surface as explicit messages, not internal DoubletFinder errors
      if (!"pca" %in% names(seu@reductions)) stop("PCA not found: run normalization, scaling and PCA before DoubletFinder.")
      n_pcs_available <- ncol(Embeddings(seu, reduction = "pca"))
      if (is.null(input$pc_use) || is.na(input$pc_use) || input$pc_use < 1) stop("Number of PCs must be at least 1.")
      if (input$pc_use > n_pcs_available) stop(sprintf("Number of PCs (%d) exceeds the %d PCs computed.", input$pc_use, n_pcs_available))
      if (is.null(input$pN_value) || is.na(input$pN_value) || input$pN_value <= 0 || input$pN_value >= 1) stop("pN must be strictly between 0 and 1.")
      if (is.null(input$doublet_rate) || is.na(input$doublet_rate) || input$doublet_rate <= 0) stop("Expected doublet rate must be greater than 0.")
      pk_mode <- input$doublet_pK_mode %||% "auto"
      if (identical(pk_mode, "manual") && (is.null(input$pK_value) || is.na(input$pK_value) || input$pK_value <= 0 || input$pK_value >= 1)) {
        stop("Manual pK must be strictly between 0 and 1.")
      }
      # Seurat v5: DoubletFinder reads a single counts matrix, so split counts layers must be joined first
      if (inherits(seu[["RNA"]], "Assay5") && length(Layers(seu[["RNA"]], search = "counts")) > 1) seu <- JoinLayers(seu, assay = "RNA")
      # Drop pANN/classification columns from previous runs so the lookup below only sees the current run
      old_df_cols <- grep("^pANN_|^DF\\.classifications_", colnames(seu@meta.data), value = TRUE)
      if (length(old_df_cols) > 0) seu@meta.data <- seu@meta.data[, setdiff(colnames(seu@meta.data), old_df_cols), drop = FALSE]
      # Flatten nested data.frame columns that break DoubletFinder metadata handling; factors keep their levels
      meta <- seu@meta.data
      seu@meta.data <- data.frame(lapply(meta, function(x) if (is.data.frame(x)) x[[1]] else if (is.factor(x)) x else as.vector(x)),
                                  row.names = rownames(meta), check.names = FALSE, stringsAsFactors = FALSE)
      pcs <- seq_len(input$pc_use)
      nExp_poi <- round(input$doublet_rate / 100 * ncol(seu))
      if (nExp_poi < 1) stop("Expected number of doublets is 0: increase the expected doublet rate.")
      if (identical(pk_mode, "manual")) {
        # User-supplied pK: the pN-pK sweep is skipped entirely
        pK <- as.numeric(input$pK_value)
        pk_source <- "manual"
      } else {
        # Artificial doublets are sampled randomly in the sweep: fixed seed for reproducibility
        set.seed(42)
        sweep.res.list <- paramSweep(seu, PCs = pcs, sct = FALSE)
        sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
        if (length(pK) != 1 || !is.finite(pK)) stop("pK optimization failed (no valid BCmetric): check PCA and number of cells.")
        pk_source <- "BCmvn maximum (pN-pK sweep)"
      }
      message("[DoubletFinder] nExp_poi: ", nExp_poi, " | optimal pK: ", pK)
      # reuse.pANN is left at the package default on purpose: older DoubletFinder releases test
      # `reuse.pANN != FALSE` (passing NULL raises "argument is of length zero"), newer releases test
      # `!is.null(reuse.pANN)` — the installed version's own default is valid in both cases
      set.seed(42)
      seu <- doubletFinder(seu, PCs = pcs, pN = input$pN_value, pK = pK, nExp = nExp_poi, sct = FALSE)
      df_col <- tail(grep("^DF\\.classifications_", colnames(seu@meta.data), value = TRUE), 1)
      if (length(df_col) == 0) stop("DoubletFinder returned no classification column.")
      message("[DoubletFinder] df_col: ", df_col)
      message("[DoubletFinder] classifications: ", paste(unique(seu@meta.data[[df_col]]), collapse = ", "))
      # Record the parameters actually applied (read by the CSV/PDF export; travels with the saved object)
      seu@misc$doubletfinder_params <- list(doublet_rate = input$doublet_rate, n_pcs = input$pc_use, pN = input$pN_value, pK = pK,
                                            pK_source = pk_source, nExp = nExp_poi, n_doublets = sum(seu@meta.data[[df_col]] == "Doublet"),
                                            n_cells = ncol(seu), seed = 42, package_version = as.character(utils::packageVersion("DoubletFinder")))
      # Fall back to PCA when no UMAP has been computed yet
      plot_reduction <- get_umap_reduction(seu) %||% "pca"
      output$doublet_umap <- renderPlot({
        DimPlot(seu, reduction = plot_reduction, group.by = df_col, cols = c("Singlet" = "grey", "Doublet" = "red")) + ggtitle("Doublets Identification")
      })
      stats <- table(seu@meta.data[[df_col]])
      stats_df <- data.frame(Category = names(stats), Count = as.numeric(stats), Percentage = round(as.numeric(stats) / sum(stats) * 100, 2))
      output$doublet_stats <- renderDT({ stats_df })
      doublet_results_data(stats_df)
      single_dataset_object(seu)
      removeModal()
      showNotification("DoubletFinder analysis completed!", type = "message")
    }, error = function(e) {
      message("[DoubletFinder] ERROR: ", conditionMessage(e)[1])
      removeModal()
      # Persistent notification so the error cannot be missed (default 5 s duration is easy to overlook)
      showNotification(paste("Error in DoubletFinder:", conditionMessage(e)[1]), type = "error", duration = NULL)
    })
  }))

  # Remove doublets
  observeEvent(input$remove_doublets, {
    req(single_dataset_object())
    tryCatch({
      seu <- single_dataset_object()
      # Column of the last DoubletFinder run (grep()[1] returned NA, never NULL, so the missing-run check never fired)
      df_col <- tail(grep("^DF\\.classifications_", colnames(seu@meta.data), value = TRUE), 1)
      if (length(df_col) == 0) {
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
      # Same reduction rule as the DoubletFinder plot: UMAP when computed, PCA otherwise; colored by the active identities
      plot_reduction <- get_umap_reduction(seu_filtered) %||% "pca"
      output$doublet_umap <- renderPlot({
        DimPlot(seu_filtered, reduction = plot_reduction) +
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

  # ---- Downloads ----

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

  ############################## Visualization of expressed genes ##############################

  # Observer pour mettre à jour les choix d'assays
  observe(run_safely("Interface update", persistent = FALSE, {
    req(single_dataset_object())
    updateAssayChoices(session, single_dataset_object())
  }))

  # Observer pour mettre à jour les choix de gènes
  observeEvent(c(single_dataset_object(), input$viz_assay), run_safely("Action", {
    req(single_dataset_object(), input$viz_assay)
    updateGeneChoices(session, single_dataset_object(), input$viz_assay)
  }))

  # Observer pour mettre à jour les sélecteurs de cluster
  observe(run_safely("Interface update", persistent = FALSE, {
    req(single_dataset_object())
    updateClusterChoices(session, single_dataset_object())
  }))

  # Observer to update text inputs when genes are selected in pickerInput
  observeEvent(input$gene_select, run_safely("gene select", {
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
  }))

  ############################## FeaturePlot - SINGLE DATASET with 3D ##############################

  # ---- Observers and computations ----

  observeEvent(input$show_feature, {
    tryCatch({
      seurat_object <- single_dataset_object()
      if (is.null(seurat_object)) { showNotification("Load a dataset first.", type = "warning"); return() }
      DefaultAssay(seurat_object) <- input$viz_assay
      check <- validate_gene_input(input$gene_list_feature, seurat_object, input$viz_assay)
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      present_genes <- check$genes
      if (length(present_genes) == 0) return()
      # Cutoffs are quantiles ("q10", kept as text for Seurat) or numeric values; empty = no cutoff
      parse_cutoff <- function(value, default) {
        if (is.null(value) || is.na(value) || !nzchar(value)) return(default)
        if (grepl("^q", value)) value else as.numeric(value)
      }
      min_cut <- parse_cutoff(input$min_cutoff, "q0")
      max_cut <- parse_cutoff(input$max_cutoff, NA)
      message(paste("FeaturePlot cutoffs: min =", min_cut, ", max =", max_cut))
      if (isTRUE(input$enable_3d_feature)) {
        # 3D mode: needs the umap3d reduction, at most 3 genes
        if (!"umap3d" %in% names(seurat_object@reductions)) {
          showNotification("3D UMAP not found. Please run 'Run UMAP 3D' first in the dimensional reduction section.", type = "warning", duration = 5)
          return()
        }
        if (length(present_genes) > 3) {
          showNotification("Maximum 3 genes allowed for 3D plot. Using first 3 genes.", type = "warning", duration = 3)
          present_genes <- present_genes[1:3]
        }
        plot <- create3DFeaturePlot(seurat_object = seurat_object, genes = present_genes, reduction = "umap3d", pt_size = 3,
                                    min_cutoff = min_cut, max_cutoff = max_cut, order = TRUE, hide_grid = isTRUE(input$hide_grid_feature),
                                    hide_axes = isTRUE(input$add_noaxes_feature), dark_mode = isTRUE(input$dark_mode_feature))
        feature_plot(plot)
      } else {
        # Seurat blends exactly two features: with any other number the co-expression option is ignored with a warning
        blend <- isTRUE(input$show_coexpression) && length(present_genes) == 2
        if (isTRUE(input$show_coexpression) && length(present_genes) != 2) {
          showNotification("Co-expression (blend) needs exactly 2 genes: genes shown separately.", type = "warning", duration = 6)
        }
        base_theme <- theme(axis.text = element_text(size = input$axis_text_size), axis.title = element_text(size = input$title_text_size),
                            plot.title = element_text(size = input$title_text_size), legend.text = element_text(size = input$axis_text_size),
                            axis.line = element_line(linewidth = input$axis_line_width), axis.ticks = element_line(linewidth = input$axis_line_width))
        plot <- FeaturePlot(seurat_object, features = present_genes, blend = blend, blend.threshold = 1, order = TRUE,
                            min.cutoff = min_cut, max.cutoff = max_cut) + base_theme
        if (isTRUE(input$add_noaxes_feature)) plot <- plot + NoAxes()
        if (isTRUE(input$add_nolegend_feature)) plot <- plot + NoLegend()
        if (isTRUE(input$dark_mode_feature)) {
          plot <- plot + theme(panel.background = element_rect(fill = "#1a1a1a"), plot.background = element_rect(fill = "#1a1a1a"),
                               panel.grid = element_blank(), text = element_text(color = "white"), axis.text = element_text(color = "white"),
                               axis.title = element_text(color = "white"), legend.background = element_rect(fill = "#1a1a1a"),
                               legend.text = element_text(color = "white"), legend.title = element_text(color = "white"))
        }
        feature_plot(plot)
      }
    }, error = function(e) {
      showNotification(paste("Error in FeaturePlot:", e$message), type = "error")
      message("FeaturePlot error: ", e$message)
    })
  })

  observeEvent(input$show_dot, {
    tryCatch({
      seurat_object <- single_dataset_object()
      if (is.null(seurat_object)) { showNotification("Load a dataset first.", type = "warning"); return() }
      check <- validate_gene_input(input$gene_list_dotplot, seurat_object, input$viz_assay)
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      if (!is.null(input$cluster_order_dot) && length(input$cluster_order_dot) > 0) {
        valid_clusters <- intersect(input$cluster_order_dot, levels(Idents(seurat_object)))
        if (length(valid_clusters) == 0) {
          showNotification("The selected clusters do not match the current cluster names (renamed?). Reselect them.", type = "error", duration = 8)
          return()
        }
        seurat_object <- subset(seurat_object, idents = valid_clusters)
        Idents(seurat_object) <- factor(Idents(seurat_object), levels = valid_clusters)
      }
      plot <- create_dotplot(
        seurat_object = seurat_object, genes = check$genes, group_by = NULL, split_by = NULL,
        color_palette = input$color_palette_dotplot_single %||% "default", dot_scale = input$dot_scale_single %||% 1,
        assay = check$assay, scale_expression = isTRUE(input$dotplot_scale_single), rotate_axis = TRUE,
        axis_text_size = input$axis_text_size, title_text_size = input$title_text_size, axis_line_width = input$axis_line_width
      )
      if (isTRUE(input$add_noaxes_dot)) plot <- plot + NoAxes()
      if (isTRUE(input$add_nolegend_dot)) plot <- plot + NoLegend()
      if (isTRUE(input$invert_axes_single)) plot <- plot + coord_flip()
      dot_plot(plot)
    }, error = function(e) {
      showNotification(paste("Error in DotPlot:", conditionMessage(e)[1]), type = "error", duration = 10)
      message("DotPlot error: ", conditionMessage(e)[1])
    })
  })

  # ---- Plots ----

  output$feature_plot <- renderPlot({
    req(feature_plot())
    plot_obj <- feature_plot()
    if(inherits(plot_obj, "gg") || inherits(plot_obj, "patchwork")) {
      txt_face <- if (isTRUE(input$bold_axis_text)) "bold" else "plain"
      plot_obj + theme(
        axis.text    = element_text(size = input$axis_text_size, face = txt_face),
        axis.title   = element_text(size = input$title_text_size, face = txt_face),
        plot.title   = element_text(size = input$title_text_size, face = txt_face),
        strip.text   = element_text(size = input$title_text_size, face = txt_face),
        legend.text  = element_text(size = input$axis_text_size),
        legend.title = element_text(size = input$axis_text_size, face = txt_face),
        axis.line    = element_line(linewidth = input$axis_line_width),
        axis.ticks   = element_line(linewidth = input$axis_line_width)
      )
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

  # Display the DotPlot — theme applied here so font/line inputs are reactive
  output$dot_plot <- renderPlot({
    req(dot_plot())
    txt_face <- if (isTRUE(input$bold_axis_text)) "bold" else "plain"
    dot_plot() +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = input$axis_text_size, face = txt_face),
        axis.text.y  = element_text(size = input$axis_text_size, face = txt_face),
        axis.title   = element_text(size = input$title_text_size, face = txt_face),
        plot.title   = element_text(size = input$title_text_size, face = txt_face),
        legend.text  = element_text(size = input$axis_text_size),
        legend.title = element_text(size = input$axis_text_size, face = txt_face),
        axis.line    = element_line(linewidth = input$axis_line_width),
        axis.ticks   = element_line(linewidth = input$axis_line_width)
      )
  })

  # SINGLE DATASET - VlnPlot and RidgePlot REPLACEMENTS
  # Replace the existing observeEvent blocks in single_dataset_server_newarch.R

  ############################## VlnPlot - SINGLE DATASET ##############################

  # ---- Observers and computations ----

  observeEvent(input$show_vln, {
    tryCatch({
      seurat_object <- single_dataset_object()
      if (is.null(seurat_object)) { showNotification("Load a dataset first.", type = "warning"); return() }
      DefaultAssay(seurat_object) <- input$viz_assay
      check <- validate_gene_input(input$gene_list_vln, seurat_object, input$viz_assay)
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      if (!is.null(input$cluster_order_vln) && length(input$cluster_order_vln) > 0) {
        valid_clusters <- intersect(input$cluster_order_vln, levels(Idents(seurat_object)))
        if (length(valid_clusters) == 0) {
          showNotification("The selected clusters do not match the current cluster names (renamed?). Reselect them.", type = "error", duration = 8)
          return()
        }
        seurat_object <- subset(seurat_object, idents = valid_clusters)
        Idents(seurat_object) <- factor(Idents(seurat_object), levels = valid_clusters)
      }
      plot <- createManualVlnPlot(
        seurat_object = seurat_object, genes = check$genes, group_by = "ident",
        pt_size = ifelse(isTRUE(input$hide_vln_points), 0, 1),
        axis_text_size = input$axis_text_size %||% 12, title_text_size = input$title_text_size %||% 14,
        axis_line_width = input$axis_line_width %||% 1,
        add_noaxes = isTRUE(input$add_noaxes_vln), add_nolegend = isTRUE(input$add_nolegend_vln)
      )
      vln_plot(plot)
    }, error = function(e) {
      showNotification(paste("Error in VlnPlot:", conditionMessage(e)[1]), type = "error", duration = 10)
      message("VlnPlot error: ", conditionMessage(e)[1])
    })
  })

  # ---- Plots ----

  output$vln_plot <- renderPlot({
    req(vln_plot())
    txt_face <- if (isTRUE(input$bold_axis_text)) "bold" else "plain"
    vln_plot() + theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = input$axis_text_size, face = txt_face),
      axis.text.y  = element_text(size = input$axis_text_size, face = txt_face),
      axis.title   = element_text(size = input$title_text_size, face = txt_face),
      plot.title   = element_text(size = input$title_text_size, face = txt_face),
      strip.text   = element_text(size = input$title_text_size, face = txt_face),
      legend.text  = element_text(size = input$axis_text_size),
      legend.title = element_text(size = input$axis_text_size, face = txt_face),
      axis.line    = element_line(linewidth = input$axis_line_width),
      axis.ticks   = element_line(linewidth = input$axis_line_width)
    )
  })

  ############################## RidgePlot - SINGLE DATASET ##############################

  # ---- Reactive state ----

  # RidgePlot: same input checks as the VlnPlot
  ridge_plot <- reactiveVal()

  # Analyse de l'expression des gènes
  number_of_nuclei <- reactiveVal(NULL)

  # Réactive pour stocker et mettre à jour le plot actuel basé sur la sélection de l'utilisateur
  current_plot <- reactiveVal()

  # ---- Observers and computations ----

  observeEvent(input$show_ridge, {
    tryCatch({
      seurat_object <- single_dataset_object()
      if (is.null(seurat_object)) { showNotification("Load a dataset first.", type = "warning"); return() }
      DefaultAssay(seurat_object) <- input$viz_assay
      check <- validate_gene_input(input$gene_list_ridge_plot, seurat_object, input$viz_assay)
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      plot <- createManualRidgePlot(
        seurat_object = seurat_object, genes = check$genes, group_by = "ident",
        axis_text_size = input$axis_text_size %||% 12, title_text_size = input$title_text_size %||% 14,
        axis_line_width = input$axis_line_width %||% 1,
        add_noaxes = isTRUE(input$add_noaxes_ridge), add_nolegend = isTRUE(input$add_nolegend_ridge)
      )
      ridge_plot(plot)
    }, error = function(e) {
      showNotification(paste("Error in RidgePlot:", conditionMessage(e)[1]), type = "error", duration = 10)
      message("RidgePlot error: ", conditionMessage(e)[1])
    })
  })

  observeEvent(input$analyze_btn, {
    req(input$gene_select_genes_analysis, single_dataset_object())

    tryCatch({
      result <- analyze_gene_expression(
        seurat_obj = single_dataset_object(),
        selected_genes = input$gene_select_genes_analysis,
        assay_name = input$viz_assay,
        expression_threshold = input$expression_threshold_single %||% 0.1,
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

  observe(run_safely("Interface update", persistent = FALSE, {
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
  }))

  # ---- Plots ----

  output$ridge_plot <- renderPlot({
    req(ridge_plot())
    txt_face <- if (isTRUE(input$bold_axis_text)) "bold" else "plain"
    ridge_plot() + theme(
      axis.text.x  = element_text(size = input$axis_text_size, face = txt_face),
      axis.text.y  = element_text(size = input$axis_text_size, face = txt_face),
      axis.title   = element_text(size = input$title_text_size, face = txt_face),
      plot.title   = element_text(size = input$title_text_size, face = txt_face),
      strip.text   = element_text(size = input$title_text_size, face = txt_face),
      legend.text  = element_text(size = input$axis_text_size),
      legend.title = element_text(size = input$axis_text_size, face = txt_face),
      axis.line    = element_line(linewidth = input$axis_line_width),
      axis.ticks   = element_line(linewidth = input$axis_line_width)
    )
  })

  # Rendre le plot actuel dans l'UI
  output$selected_plot_display <- renderPlot({
    req(current_plot())
    current_plot()
  })

  # ---- Downloads ----

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
    ),
    plot_modifier = function(p) {
      txt_face <- if (isTRUE(input$bold_axis_text)) "bold" else "plain"
      p + theme(
        axis.text.x  = element_text(size = input$axis_text_size, face = txt_face),
        axis.text.y  = element_text(size = input$axis_text_size, face = txt_face),
        axis.title   = element_text(size = input$title_text_size, face = txt_face),
        plot.title   = element_text(size = input$title_text_size, face = txt_face),
        strip.text   = element_text(size = input$title_text_size, face = txt_face),
        legend.text  = element_text(size = input$axis_text_size),
        legend.title = element_text(size = input$axis_text_size, face = txt_face),
        axis.line    = element_line(linewidth = input$axis_line_width),
        axis.ticks   = element_line(linewidth = input$axis_line_width)
      )
    }
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
    ),
    plot_modifier = function(p) {
      txt_face  <- if (isTRUE(input$bold_axis_text)) "bold" else "plain"
      gene_face <- if (isTRUE(input$bold_axis_text)) "bold.italic" else "italic"
      p + theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = input$axis_text_size, face = txt_face),
        axis.text.y  = element_text(size = input$axis_text_size, face = gene_face),
        axis.title   = element_text(size = input$title_text_size, face = txt_face),
        plot.title   = element_text(size = input$title_text_size, face = txt_face),
        legend.text  = element_text(size = input$axis_text_size),
        legend.title = element_text(size = input$axis_text_size, face = txt_face),
        axis.line    = element_line(linewidth = input$axis_line_width),
        axis.ticks   = element_line(linewidth = input$axis_line_width)
      )
    }
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
    ),
    plot_modifier = function(p) {
      txt_face  <- if (isTRUE(input$bold_axis_text)) "bold" else "plain"
      p + theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = input$axis_text_size, face = txt_face),
        axis.text.y  = element_text(size = input$axis_text_size, face = txt_face),
        axis.title   = element_text(size = input$title_text_size, face = txt_face),
        plot.title   = element_text(size = input$title_text_size, face = txt_face),
        strip.text   = element_text(size = input$title_text_size, face = txt_face),
        legend.text  = element_text(size = input$axis_text_size),
        legend.title = element_text(size = input$axis_text_size, face = txt_face),
        axis.line    = element_line(linewidth = input$axis_line_width),
        axis.ticks   = element_line(linewidth = input$axis_line_width)
      )
    }
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
    ),
    plot_modifier = function(p) {
      txt_face  <- if (isTRUE(input$bold_axis_text)) "bold" else "plain"
      p + theme(
        axis.text    = element_text(size = input$axis_text_size, face = txt_face),
        axis.title   = element_text(size = input$title_text_size, face = txt_face),
        plot.title   = element_text(size = input$title_text_size, face = txt_face),
        legend.text  = element_text(size = input$axis_text_size),
        legend.title = element_text(size = input$axis_text_size, face = txt_face),
        axis.line    = element_line(linewidth = input$axis_line_width),
        axis.ticks   = element_line(linewidth = input$axis_line_width)
      )
    }
  )

  output$download_genes_number_expresion <- createDownloadHandler(
    reactive_data = number_of_nuclei,
    object_name_reactive = reactive({
      getObjectNameForDownload(single_dataset_object(), default_name = "GeneExpression")
    }),
    data_name = "number_of_gene_expression",
    download_type = "csv"
  )

  ############################## Heatmap single dataset ##############################

  # ---- Reactive state ----
  heatmap_plot <- reactiveVal()
  heatmap_n_genes <- reactiveVal(0)  # Track number of genes for dynamic height

  # Store top genes table for CSV export (NULL when using manual gene input)
  heatmap_markers_data_single <- reactiveVal(NULL)

  scatter_plot_single <- reactiveVal()

  # ---- Observers and computations ----

  # Update selectors when data changes
  observeEvent(single_dataset_object(), run_safely("Interface update", persistent = FALSE, {
    req(single_dataset_object())
    seurat_obj <- single_dataset_object()

    # Update clusters
    clusters <- isolate(active_clusters_single())
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
  }))

  # Update the Group 1 / Group 2 choices when the split column changes
  observeEvent(input$split_by_heatmap_single, run_safely("split by heatmap single", persistent = FALSE, {
    req(single_dataset_object())
    split_by <- input$split_by_heatmap_single
    values <- character(0)
    if (!is.null(split_by) && !split_by %in% c("", "None") && split_by %in% colnames(single_dataset_object()@meta.data)) {
      values <- sort(unique(trimws(as.character(single_dataset_object()@meta.data[[split_by]]))))
      values <- values[!is.na(values) & values != ""]
    }
    updateSelectizeInput(session, "heatmap_group1_single", choices = values, selected = character(0), server = TRUE)
    updateSelectizeInput(session, "heatmap_group2_single", choices = values, selected = character(0), server = TRUE)
  }))

  # Generate heatmap: clusters, or cluster x split value with the second axis drawn as panels (build_grouped_expression_heatmap)
  observeEvent(input$generateCustomHeatmap, run_safely("Heatmap (single dataset)", {
    if (is.null(single_dataset_object())) stop("Load a dataset first")
    showModal(modalDialog(title = "Generating heatmap...", "Please wait while generating average expression heatmap", easyClose = FALSE, footer = NULL))
    obj <- single_dataset_object()
    split_by <- input$split_by_heatmap_single
    if (is.null(split_by) || split_by %in% c("", "None")) split_by <- NULL
    # Clusters read from the column that fills the cluster selector, so renamed clusters match the selection
    cluster_col <- active_cluster_column_single()
    cluster_values <- if (!is.null(cluster_col)) as.character(obj@meta.data[[cluster_col]]) else as.character(Idents(obj))
    res <- build_grouped_expression_heatmap(obj, cluster_values = cluster_values, clusters = input$clusters_heatmap_single, split_by = split_by,
                                            gene_source = input$gene_source_heatmap_single %||% "manual", gene_text = input$gene_select_heatmap,
                                            n_top = input$n_top_genes_single %||% 10, group1 = input$heatmap_group1_single, group2 = input$heatmap_group2_single,
                                            direction = input$heatmap_de_direction_single %||% "both", order_by = input$heatmap_order_single %||% "cluster",
                                            scale_mode = input$heatmap_scale_single %||% "zscore", color_palette = input$color_palette_heatmap_single %||% "viridis",
                                            assay = if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj))
    heatmap_markers_data_single(res$markers_table)
    heatmap_n_genes(length(res$genes))
    heatmap_plot(res$plot)
    removeModal()
    if (!is.null(res$gene_warning)) showNotification(res$gene_warning, type = "warning", duration = 8)
    if (length(res$skipped) > 0) showNotification(paste("Clusters without result:", paste(res$skipped, collapse = "; ")), type = "warning", duration = 15)
    if (identical(input$heatmap_scale_single, "zscore") && res$n_columns == 2) showNotification("Row z-score over 2 columns only shows the direction (every value is +/-0.71). Choose 'Average expression' to compare magnitudes.", type = "warning", duration = 12)
    showNotification(sprintf("Heatmap with %d genes and %d columns", length(res$genes), res$n_columns), type = "message")
  }))

  # Mettez à jour les choix de selectInput pour les gènes/features et les clusters
  observe(run_safely("Interface update", persistent = FALSE, {
    req(single_dataset_object())
    updateGeneChoices(session, single_dataset_object(), "RNA", c("feature1_select", "feature2_select"))
    updateClusterTextInputs(session, single_dataset_object(), "scatter_text_clusters")
    }))

  # Update color_by and split_by choices when data changes
  observeEvent(single_dataset_object(), run_safely("Interface update", persistent = FALSE, {
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
  }))

  # Generate scatter plot
  observeEvent(input$generateScatter_single, run_safely("generateScatter single", {
    req(single_dataset_object())

    showModal(modalDialog(
      title = "Processing",
      div(h4("Generating Scatter Plot...", style = "text-align: center;")),
      footer = NULL,
      easyClose = FALSE
    ))

    seurat_obj <- single_dataset_object()

    tryCatch({
      # One gene per field, case-insensitive (the first name is used if several are typed)
      selected_assay <- "RNA"
      check1 <- validate_gene_input(input$feature1_select, seurat_obj, selected_assay)
      check2 <- validate_gene_input(input$feature2_select, seurat_obj, selected_assay)
      if (!check1$ok || !check2$ok) {
        showNotification(paste(c(if (!check1$ok) paste("First gene:", check1$message), if (!check2$ok) paste("Second gene:", check2$message)), collapse = " | "),
                         type = "error", duration = 8)
        removeModal()
        return()
      }
      selected_assay <- check1$assay
      gene1 <- check1$genes[1]
      gene2 <- check2$genes[1]

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
  }))

  # ---- Plots ----

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

  output$scatter_plot_single <- renderPlot({
    req(scatter_plot_single())
    scatter_plot_single()
  })

  # ---- Downloads ----

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

  # Reactive column name — resolves which metadata column backs the cluster selector.
  active_cluster_column_single <- reactive({
    req(single_dataset_object())
    obj <- single_dataset_object()
    col_input <- if (!is.null(input$cluster_col_single) && nzchar(input$cluster_col_single)) {
      input$cluster_col_single
    } else if (!is.null(input$umap_col_single) && nzchar(input$umap_col_single)) {
      input$umap_col_single
    } else {
      NULL
    }
    if (!is.null(col_input) && col_input %in% colnames(obj@meta.data)) {
      col_input
    } else if ("ClusterIdents" %in% colnames(obj@meta.data)) {
      "ClusterIdents"
    } else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
      "seurat_clusters"
    } else {
      NULL
    }
  })

  active_clusters_single <- reactive({
    req(single_dataset_object())
    obj <- single_dataset_object()
    col <- active_cluster_column_single()
    if (!is.null(col)) unique(as.character(obj@meta.data[[col]]))
    else as.character(unique(Idents(obj)))
  })

  observe(run_safely("Interface update", persistent = FALSE, {
    req(single_dataset_object())
    forced <- isolate(forced_col_single())
    col <- if (!is.null(forced)) {
      forced_col_single(NULL)
      forced
    } else {
      isolate(input$cluster_col_single)
    }
    meta_cols <- colnames(single_dataset_object()@meta.data)
    selected <- if (!is.null(col) && nzchar(col) && col %in% meta_cols) col else NULL
    update_umap_selectors(session, single_dataset_object(), module = "single",
                          selected_col = selected)
  }))

  # Keep the two "Color by" selectors in sync across tabs (cluster_col_single ↔ umap_col_single)
  # and the two "UMAP reduction" selectors (cluster_reduction_single ↔ umap_reduction_single).
  observeEvent(input$cluster_col_single, run_safely("cluster col single", {
    req(nzchar(input$cluster_col_single))
    if (!identical(input$umap_col_single, input$cluster_col_single))
      updateSelectInput(session, "umap_col_single", selected = input$cluster_col_single)
  }), ignoreInit = TRUE)

  observeEvent(input$umap_col_single, run_safely("umap col single", {
    req(nzchar(input$umap_col_single))
    if (!identical(input$cluster_col_single, input$umap_col_single))
      updateSelectInput(session, "cluster_col_single", selected = input$umap_col_single)
  }), ignoreInit = TRUE)

  observeEvent(input$cluster_reduction_single, run_safely("cluster reduction single", {
    req(nzchar(input$cluster_reduction_single))
    if (!identical(input$umap_reduction_single, input$cluster_reduction_single))
      updateSelectInput(session, "umap_reduction_single", selected = input$cluster_reduction_single)
  }), ignoreInit = TRUE)

  observeEvent(input$umap_reduction_single, run_safely("umap reduction single", {
    req(nzchar(input$umap_reduction_single))
    if (!identical(input$cluster_reduction_single, input$umap_reduction_single))
      updateSelectInput(session, "cluster_reduction_single", selected = input$umap_reduction_single)
  }), ignoreInit = TRUE)

  # Save initial cluster state once on load
  observe(run_safely("Interface update", persistent = FALSE, {
    req(single_dataset_object())
    if (length(cluster_history_single$states) == 0) {
      save_cluster_state_single(single_dataset_object(), "Initial state")
    }
  }))

  # Central cluster observer — no implicit dependency on umap_col_single
  observe(run_safely("Interface update", persistent = FALSE, {
    req(single_dataset_object(), active_clusters_single())
    clusters <- active_clusters_single()
    updateSelectInput(session, "select_cluster",                  choices = clusters)
    updateSelectInput(session, "cluster_select",                  choices = clusters)
    updateSelectInput(session, "target_cluster_global_single",    choices = clusters)
    updateSelectInput(session, "cluster1_pairwise_single",        choices = clusters)
    updateSelectInput(session, "exclusive_target_cluster_single", choices = clusters)
    updateSelectizeInput(session, "clusters_heatmap_single", choices = clusters, selected = clusters, server = TRUE)
    # VlnPlot and DotPlot subset cells with subset(idents = ...): their order selectors list the active identities,
    # refreshed after clustering, renaming and annotation; still-valid selections are kept
    ident_levels <- levels(droplevels(as.factor(Idents(single_dataset_object()))))
    updateSelectInput(session, "cluster_order_vln", choices = ident_levels, selected = intersect(isolate(input$cluster_order_vln), ident_levels))
    updateSelectizeInput(session, "cluster_order_dot", choices = ident_levels, selected = intersect(isolate(input$cluster_order_dot), ident_levels))
    message("Central cluster observer — ", length(clusters), " clusters — col: ", isolate(input$cluster_col_single) %||% "auto-detected")
  }))

 ##################### Cluster History System for Undo ####################

  # ---- Reactive state ----

  #Reactive variable for that tab
  cluster_colours <- reactiveVal()  # Initialize a list to store cluster colors

  # ---- Helpers ----

  # Snapshot of the cluster labels taken before a change, named by cell barcode (labels default to the active identities)
  save_cluster_state_single <- function(seurat_obj, action_description = "Cluster modification", labels = NULL) {
    if (is.null(labels)) labels <- as.character(Idents(seurat_obj))
    new_state <- list(clusters = stats::setNames(as.character(labels), colnames(seurat_obj)), action = action_description)
    cluster_history_single$states <- append(cluster_history_single$states, list(new_state))
    # Limit history size
    if (length(cluster_history_single$states) > cluster_history_single$max_history) cluster_history_single$states <- cluster_history_single$states[-1]
  }

  # Undo: restore the last snapshot (taken before the last change) and drop it; the first entry is the initial state and is kept
  undo_last_action_single <- function() {
    n_states <- length(cluster_history_single$states)
    if (n_states <= 1) return(NULL)
    restored <- cluster_history_single$states[[n_states]]$clusters
    cluster_history_single$states <- cluster_history_single$states[-n_states]
    restored
  }

  # ---- Observers and computations ----

  # Function for renaming each cluster: labels come from the column shown in the cluster selector; only the identities,
  # ClusterIdents and cluster_names are written (other metadata columns such as genotype or sample are never modified)
  observeEvent(input$rename_single_cluster_button, run_safely("rename single cluster button", {
    req(input$select_cluster, single_dataset_object())
    new_name <- trimws(input$rename_single_cluster %||% "")
    if (!nzchar(new_name)) stop("Enter a new name for the cluster")
    updated_seurat <- single_dataset_object()
    cluster_col <- active_cluster_column_single()
    labels <- if (!is.null(cluster_col)) as.character(updated_seurat@meta.data[[cluster_col]]) else as.character(Idents(updated_seurat))
    target <- !is.na(labels) & labels == input$select_cluster
    if (!any(target)) stop(sprintf("Cluster '%s' not found in column %s: reselect it", input$select_cluster, cluster_col %||% "Idents()"))
    save_cluster_state_single(updated_seurat, paste("Rename", input$select_cluster), labels = labels)
    # Renaming to an existing name merges the two clusters
    merged <- new_name %in% labels[!target]
    labels[target] <- new_name
    Idents(updated_seurat) <- factor(stats::setNames(labels, colnames(updated_seurat)))
    updated_seurat$cluster_names <- labels
    updated_seurat$ClusterIdents <- labels
    updateClusterSelectors(session, updated_seurat, cluster_col = "ClusterIdents", rename_id = "select_cluster", color_id = "cluster_select")
    forced_col_single("ClusterIdents")
    single_dataset_object(updated_seurat)
    showNotification(if (merged) paste("Clusters merged under the name:", new_name) else paste("Cluster renamed to:", new_name), type = "message")
  }))

  # Undo: labels restored by cell barcode (Idents<- would recycle them by position once cells have been removed);
  # cells absent from the snapshot keep their current label; other metadata columns are never modified
  observeEvent(input$undo_cluster_single, run_safely("undo cluster single", {
    req(single_dataset_object())
    restored_clusters <- undo_last_action_single()
    if (is.null(restored_clusters)) {
      showNotification("Nothing to undo", type = "warning")
      return()
    }
    updated_seurat <- single_dataset_object()
    labels <- unname(restored_clusters[colnames(updated_seurat)])
    if (all(is.na(labels))) stop("The saved state belongs to another object (dataset reloaded?): nothing restored")
    current <- as.character(Idents(updated_seurat))
    labels[is.na(labels)] <- current[is.na(labels)]
    Idents(updated_seurat) <- factor(stats::setNames(labels, colnames(updated_seurat)))
    updated_seurat$cluster_names <- labels
    updated_seurat$ClusterIdents <- labels
    updateClusterSelectors(session, updated_seurat, cluster_col = "ClusterIdents", rename_id = "select_cluster", color_id = "cluster_select")
    forced_col_single("ClusterIdents")
    single_dataset_object(updated_seurat)
    showNotification("Undid last action", type = "message")
  }))

  ############################## Automated annotation (reference label transfer + SingleR) ##############################

  # ---- Reactive state ----
  transfer_result_single <- reactiveVal(NULL)
  singler_result_single <- reactiveVal(NULL)
  annotation_reference_single <- reactiveVal(NULL)

  # celldex references are cached on disk (pre-filled in the Docker image) so they are downloaded only once
  singler_cache_dir <- Sys.getenv("CELLHUB_REF_CACHE", file.path(tools::R_user_dir("CellHub", "cache"), "singler_references"))

  # ---- Helpers ----

  # Cluster vector used for the per-cluster vote (same column as the cluster selectors)
  annotation_clusters_single <- function(obj) {
    cluster_col <- active_cluster_column_single()
    clusters <- if (!is.null(cluster_col)) as.character(obj@meta.data[[cluster_col]]) else as.character(Idents(obj))
    list(column = cluster_col, values = clusters)
  }

  # Cluster summary recomputed with the method's current agreement threshold (key = "transfer" or "singler")
  annotation_summary_single <- function(result, key) {
    min_agreement <- input[[paste0(key, "_min_agreement_single")]] %||% 50
    summarize_annotation_by_cluster(result$cells, min_agreement / 100, result$params$max_unassigned)
  }

  # Rename clusters from an annotation result (undoable, like manual renaming)
  apply_annotation_labels_single <- function(result, method_name, key) {
    obj <- single_dataset_object()
    cluster_col <- result$params$cluster_col
    current <- if (!is.null(cluster_col) && cluster_col %in% colnames(obj@meta.data)) {
      as.character(obj@meta.data[[cluster_col]])
    } else {
      as.character(Idents(obj))
    }
    new_names <- build_annotation_cluster_names(annotation_summary_single(result, key),
                                                keep_low_confidence = isTRUE(input[[paste0(key, "_keep_low_single")]]),
                                                make_unique = isTRUE(input[[paste0(key, "_unique_single")]]))
    labels <- unname(new_names[current])
    labels[is.na(labels)] <- current[is.na(labels)]
    save_cluster_state_single(obj, paste(method_name, "labels"))
    Idents(obj) <- factor(stats::setNames(labels, colnames(obj)))
    obj$cluster_names <- labels
    obj$ClusterIdents <- labels
    updateClusterSelectors(session, obj,
                           cluster_col = "ClusterIdents",
                           rename_id   = "select_cluster",
                           color_id    = "cluster_select")
    forced_col_single("ClusterIdents")
    single_dataset_object(obj)
    showNotification(sprintf("%s labels applied to %d of %d clusters", method_name,
                             sum(new_names != names(new_names)), length(new_names)), type = "message", duration = 6)
  }

  # Cluster summary table with confidence colouring
  render_annotation_summary_single <- function(result, key) {
    df <- annotation_summary_single(result, key)
    colnames(df) <- c("Cluster", "Cells", "Predicted label", "Agreement (%)", "2nd label", "2nd (%)",
                      "Unassigned (%)", "Median score", "Median delta", "Confidence")
    datatable(df, rownames = FALSE, options = list(pageLength = 25, scrollX = TRUE)) %>%
      formatStyle("Confidence", target = "row", backgroundColor = styleEqual(c("High", "Low"), c("#e8f5e9", "#ffebee")))
  }

  # Status line under the Run button: not run yet, or a one-line summary of the last run
  render_annotation_status_single <- function(result, key, extra = NULL) {
    if (is.null(result)) {
      return(tags$p(icon("circle-info"), " Not run yet.", style = "margin-top: 10px; font-size: 12px; color: #6c757d;"))
    }
    s <- annotation_summary_single(result, key)
    n_low <- sum(s$confidence == "Low")
    tags$div(
      style = "margin-top: 10px; padding: 8px 10px; border-radius: 6px; font-size: 12px; background-color: #e8f5e9; color: #2e7d32;",
      icon("circle-check"), sprintf(" %d clusters annotated: %d confident, %d to check.", nrow(s), nrow(s) - n_low, n_low),
      if (!is.null(extra)) tags$div(extra, style = "color: #6c757d; margin-top: 4px;")
    )
  }

  # UMAP coloured by a per-cell annotation column or by its confidence score
  plot_annotation_umap_single <- function(column, title, continuous = FALSE) {
    obj <- single_dataset_object()
    reduc <- get_umap_reduction(obj)
    shiny::validate(shiny::need(!is.null(reduc), "No UMAP found: run the clustering step first."))
    req(column %in% colnames(obj@meta.data))
    if (continuous) {
      FeaturePlot(obj, features = column, reduction = reduc, order = TRUE) + ggtitle(title)
    } else {
      DimPlot(obj, reduction = reduc, group.by = column, label = TRUE, repel = TRUE, label.size = 3) + ggtitle(title)
    }
  }

  # ---- Observers and computations ----

  # A new dataset invalidates previous annotation results and the cluster history (snapshots of another object)
  observeEvent(list(input$file, input$load_seurat_file), run_safely("Action", {
    singler_result_single(NULL)
    transfer_result_single(NULL)
    cluster_history_single$states <- list()
  }), ignoreInit = TRUE)

  ############ Label transfer from a user reference ############

  # ---- Observers and computations ----

  observeEvent(input$annotation_reference_file_single, run_safely("annotation reference file single", {
    showModal(modalDialog(title = "Loading reference", "Reading the annotated reference dataset...", easyClose = FALSE, footer = NULL))
    tryCatch({
      ref <- read_reference_object(input$annotation_reference_file_single$datapath)
      cols <- get_reference_label_columns(ref)
      if (length(cols) == 0) stop("No categorical metadata column with 2-300 values found in the reference")
      preferred <- intersect(c("cell_type", "celltype", "CellType", "cell_types", "annotation", "ClusterIdents", "cluster_names"), cols)
      annotation_reference_single(ref)
      updateSelectInput(session, "annotation_reference_label_single", choices = cols,
                        selected = if (length(preferred) > 0) preferred[1] else cols[1])
      transfer_result_single(NULL)
      removeModal()
      showNotification(sprintf("Reference loaded: %s cells, %s genes", format(ncol(ref), big.mark = ","), format(nrow(ref), big.mark = ",")),
                       type = "message")
    }, error = function(e) {
      removeModal()
      annotation_reference_single(NULL)
      showNotification(paste("Error loading reference:", conditionMessage(e)[1]), type = "error", duration = 10)
    })
  }))

  observeEvent(input$run_transfer_single, run_safely("run transfer single", {
    req(single_dataset_object())
    if (is.null(annotation_reference_single())) {
      showNotification("Load an annotated reference first (step 1).", type = "warning")
      return()
    }
    req(input$annotation_reference_label_single)
    showModal(modalDialog(title = "Label transfer from reference",
                          "Finding anchors between the reference and your dataset... This can take several minutes.",
                          easyClose = FALSE, footer = NULL))
    tryCatch({
      obj <- single_dataset_object()
      cl <- annotation_clusters_single(obj)
      ref <- prepare_reference_object(annotation_reference_single(), input$annotation_reference_label_single,
                                      max_cells_per_label = input$transfer_max_cells_single)
      result <- run_label_transfer(obj, ref, label_column = input$annotation_reference_label_single, clusters = cl$values,
                                   dims = input$transfer_dims_single, min_score = input$transfer_min_score_single,
                                   min_agreement = input$transfer_min_agreement_single / 100)
      result$params$reference <- input$annotation_reference_file_single$name
      result$params$cluster_col <- cl$column
      single_dataset_object(add_annotation_metadata(obj, result, prefix = "RefTransfer"))
      transfer_result_single(result)
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste("Label transfer error:", conditionMessage(e)[1]), type = "error", duration = 10)
      message("Label transfer error: ", conditionMessage(e)[1])
    })
  }))

  observeEvent(input$apply_transfer_labels_single, run_safely("apply transfer labels single", {
    req(single_dataset_object())
    if (is.null(transfer_result_single())) {
      showNotification("Run the label transfer first (step 2).", type = "warning")
      return()
    }
    tryCatch(apply_annotation_labels_single(transfer_result_single(), "Reference", "transfer"),
             error = function(e) showNotification(paste("Error applying labels:", conditionMessage(e)[1]), type = "error", duration = 10))
  }))

  # ---- UI and text outputs ----

  # Reference summary: size and the labels found in the chosen column
  output$annotation_reference_info_single <- renderUI({
    ref <- annotation_reference_single()
    if (is.null(ref)) {
      return(tags$p(icon("circle-info"), " No reference loaded.", style = "font-size: 12px; color: #6c757d;"))
    }
    col <- input$annotation_reference_label_single
    labels <- if (!is.null(col) && col %in% colnames(ref@meta.data)) sort(unique(na.omit(as.character(ref@meta.data[[col]])))) else character(0)
    tags$div(
      style = "font-size: 12px; background-color: #f3f0ff; padding: 8px 10px; border-radius: 6px;",
      tags$p(tags$strong(format(ncol(ref), big.mark = ",")), " cells, ", tags$strong(format(nrow(ref), big.mark = ",")), " genes, ",
             tags$strong(length(labels)), " labels", style = "margin: 0 0 4px 0;"),
      tags$p(paste(head(labels, 12), collapse = ", "), if (length(labels) > 12) sprintf(" ... (+%d)", length(labels) - 12),
             style = "margin: 0; color: #6c757d;")
    )
  })

  output$transfer_status_single <- renderUI({
    res <- transfer_result_single()
    extra <- if (!is.null(res)) sprintf("%s anchors, %s reference cells, %s shared genes", format(res$params$n_anchors, big.mark = ","),
                                        format(res$params$n_reference_cells, big.mark = ","), format(res$params$n_common_genes, big.mark = ",")) else NULL
    render_annotation_status_single(res, "transfer", extra)
  })

  # ---- Plots ----

  output$transfer_umap_labels_single <- renderPlot({
    req(transfer_result_single())
    plot_annotation_umap_single("RefTransfer_pruned_label", "Transferred labels (per cell)")
  })

  output$transfer_umap_score_single <- renderPlot({
    req(transfer_result_single())
    plot_annotation_umap_single("RefTransfer_score", "Prediction score (max)", continuous = TRUE)
  })

  # ---- Tables ----

  output$transfer_summary_single <- renderDT({
    req(transfer_result_single())
    render_annotation_summary_single(transfer_result_single(), "transfer")
  })

  # ---- Downloads ----

  output$download_transfer_clusters_single <- downloadHandler(
    filename = function() paste0(getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset"),
                                 "_RefTransfer_cluster_summary_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      req(transfer_result_single())
      res <- transfer_result_single()
      out <- annotation_summary_single(res, "transfer")
      out$reference <- res$params$reference
      out$label_column <- res$params$label_column
      out$min_score <- res$params$min_score
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$download_transfer_cells_single <- downloadHandler(
    filename = function() paste0(getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset"),
                                 "_RefTransfer_per_cell_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      req(transfer_result_single())
      write.csv(transfer_result_single()$cells, file, row.names = FALSE)
    }
  )

  ############ SingleR (built-in celldex references) ############

  # ---- Helpers ----

  # Stored cluster colors, returned only when they cover every value of the displayed cluster column
  # (keyed by cluster name, they would leave points uncoloured when the UMAP is coloured by another column)
  get_cluster_colors <- function(seurat_object) {
    if (is.null(seurat_object)) return(NULL)
    colors <- find_stored_colors(seurat_object, cluster_names = levels(Idents(seurat_object)))
    if (is.null(colors)) {
      message("Using Seurat default colors")
      return(NULL)
    }
    shown <- tryCatch(active_clusters_single(), error = function(e) NULL)
    if (!is.null(shown) && !all(shown[!is.na(shown)] %in% names(colors))) {
      message("Stored colors do not cover the displayed column: using default colors")
      return(NULL)
    }
    message("Using stored cluster colors")
    return(colors)
  }

  # ---- Observers and computations ----

  # Reference choices follow the selected species
  observeEvent(input$singler_species_single, run_safely("singler species single", {
    refs <- get_singler_references(input$singler_species_single)
    updateSelectInput(session, "singler_reference_single", choices = refs, selected = if (length(refs) > 0) refs[[1]] else NULL)
  }))

  # Default the SingleR species to the one chosen when loading raw data
  observeEvent(input$species_choice, run_safely("species choice", {
    if (input$species_choice %in% c("human", "mouse")) {
      updateSelectInput(session, "singler_species_single", selected = input$species_choice)
    }
  }))

  observeEvent(input$run_singler_single, run_safely("run singler single", {
    req(single_dataset_object(), input$singler_reference_single)
    showModal(modalDialog(title = "Automated annotation (SingleR)",
                          "Loading the reference and annotating cells... This can take several minutes for large datasets.",
                          easyClose = FALSE, footer = NULL))
    tryCatch({
      obj <- single_dataset_object()
      cl <- annotation_clusters_single(obj)
      reference <- load_singler_reference(input$singler_reference_single, singler_cache_dir)
      result <- run_singler_annotation(obj, reference, label_column = input$singler_label_level_single, clusters = cl$values,
                                       mode = input$singler_mode_single, assay = DefaultAssay(obj),
                                       min_agreement = input$singler_min_agreement_single / 100)
      result$params$reference <- input$singler_reference_single
      result$params$cluster_col <- cl$column
      single_dataset_object(add_annotation_metadata(obj, result, prefix = "SingleR"))
      singler_result_single(result)
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste("SingleR error:", conditionMessage(e)[1]), type = "error", duration = 10)
      message("SingleR error: ", conditionMessage(e)[1])
    })
  }))

  observeEvent(input$apply_singler_labels_single, run_safely("apply singler labels single", {
    req(single_dataset_object())
    if (is.null(singler_result_single())) {
      showNotification("Run SingleR first (step 2).", type = "warning")
      return()
    }
    tryCatch(apply_annotation_labels_single(singler_result_single(), "SingleR", "singler"),
             error = function(e) showNotification(paste("Error applying labels:", conditionMessage(e)[1]), type = "error", duration = 10))
  }))

  # Update color picker to reflect the stored color of the selected cluster
  observeEvent(input$cluster_select, run_safely("cluster select", {
    req(single_dataset_object(), input$cluster_select)
    colors <- get_cluster_colors(single_dataset_object())
    if (!is.null(colors) && input$cluster_select %in% names(colors)) {
      updateColourInput(session, "cluster_colour", value = colors[[input$cluster_select]])
    }
  }), ignoreInit = TRUE)

  observeEvent(input$update_colour, run_safely("update colour", {
    message("Update the color of the selected cluster")
    updated_seurat <- single_dataset_object()
    cluster_colors <- get_cluster_colors(updated_seurat)

    if (is.null(cluster_colors)) {
      message("No colors stored, creating default colors")
      clusters       <- active_clusters_single()
      n_clusters     <- length(clusters)
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
  }))

  # ---- UI and text outputs ----

  output$singler_status_single <- renderUI({
    res <- singler_result_single()
    extra <- if (!is.null(res)) sprintf("%s, %s shared genes", res$params$reference, format(res$params$n_common_genes, big.mark = ",")) else NULL
    render_annotation_status_single(res, "singler", extra)
  })

  # ---- Plots ----

  output$singler_umap_labels_single <- renderPlot({
    req(singler_result_single())
    plot_annotation_umap_single("SingleR_pruned_label", "SingleR labels (per cell)")
  })

  output$singler_umap_delta_single <- renderPlot({
    req(singler_result_single())
    plot_annotation_umap_single("SingleR_delta", "SingleR confidence (delta to next-best label)", continuous = TRUE)
  })

  # Display final UMAP with 2D/3D support
  output$umap_finale <- renderPlotly({
    req(single_dataset_object())
    seurat_obj  <- single_dataset_object()
    reduc_input <- if (!is.null(input$cluster_reduction_single) && nzchar(input$cluster_reduction_single)) {
      input$cluster_reduction_single
    } else if (!is.null(input$umap_reduction_single) && nzchar(input$umap_reduction_single)) {
      input$umap_reduction_single
    } else {
      NULL
    }
    umap_reduc <- if (!is.null(reduc_input) && reduc_input %in% names(seurat_obj@reductions)) {
      reduc_input
    } else {
      get_umap_reduction(seurat_obj)
    }
    if (is.null(umap_reduc)) return(NULL)
    use_3d <- isTRUE(input$umap_3d_toggle)
    cluster_col <- if (!is.null(input$umap_col_single) &&
                       nzchar(input$umap_col_single) &&
                       input$umap_col_single %in% colnames(seurat_obj@meta.data)) {
      input$umap_col_single
    } else {
      NULL
    }
    cluster_values <- if (!is.null(cluster_col)) {
      as.character(seurat_obj@meta.data[[cluster_col]])
    } else {
      as.character(Idents(seurat_obj))
    }
    cluster_colors <- get_cluster_colors(seurat_obj)
    if (is.null(cluster_colors)) {
      n <- length(unique(cluster_values))
      cluster_colors <- setNames(scales::hue_pal()(n), unique(cluster_values))
      message("Using Seurat default colors")
    } else {
      message("Using custom cluster colors")
    }
    if (use_3d && hasValid3DUMAP(seurat_obj)) {
      umap_coords <- as.data.frame(Embeddings(seurat_obj, reduction = "umap3d"))
      colnames(umap_coords) <- c("UMAP_1", "UMAP_2", "UMAP_3")
      umap_coords$cluster <- cluster_values
      centers <- umap_coords %>%
        group_by(cluster) %>%
        summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2),
                  UMAP_3 = median(UMAP_3), .groups = "drop")
      hide_grid  <- isTRUE(input$umap_3d_hide_grid)
      hide_axes  <- isTRUE(input$umap_3d_hide_axes)
      dark_mode  <- isTRUE(input$umap_3d_dark_mode)
      bg_color   <- if (dark_mode) "#1a1a1a" else "white"
      grid_color <- if (dark_mode) "#404040" else "#e0e0e0"
      axis_color <- if (dark_mode) "#ffffff" else "#000000"
      text_color <- if (dark_mode) "#ffffff" else "#000000"
      axis_config <- list(
        title = "", showgrid = !hide_grid, showline = !hide_axes,
        zeroline = !hide_axes, showticklabels = !hide_axes,
        gridcolor = grid_color, gridwidth = 1, linecolor = axis_color,
        tickcolor = axis_color, tickfont = list(color = text_color)
      )
      xaxis_config <- axis_config; yaxis_config <- axis_config; zaxis_config <- axis_config
      if (!hide_axes) {
        xaxis_config$title <- list(text = "UMAP_1", font = list(color = text_color))
        yaxis_config$title <- list(text = "UMAP_2", font = list(color = text_color))
        zaxis_config$title <- list(text = "UMAP_3", font = list(color = text_color))
      }
      p <- plot_ly() %>%
        add_trace(data = umap_coords, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
                  color = ~cluster, colors = cluster_colors,
                  type = "scatter3d", mode = "markers",
                  marker = list(size = input$pt_size, opacity = 0.7),
                  hoverinfo = "text", text = ~paste("Cluster:", cluster)) %>%
        add_trace(data = centers, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
                  type = "scatter3d", mode = "text", text = ~cluster,
                  textfont = list(size = input$label_font_size, color = text_color),
                  hoverinfo = "none", showlegend = FALSE) %>%
        layout(
          scene = list(xaxis = xaxis_config, yaxis = yaxis_config, zaxis = zaxis_config,
                       bgcolor = bg_color, camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))),
          paper_bgcolor = bg_color, plot_bgcolor = bg_color, showlegend = FALSE,
          title = list(text = input$plot_title, font = list(size = 24, color = text_color))
        )
      message("3D UMAP generated")
      return(p)
    } else if (use_3d && !hasValid3DUMAP(seurat_obj)) {
      return(plotly_empty() %>%
               layout(title = list(text = "Please compute 3D UMAP first using the button above",
                                   font = list(size = 16))))
    } else {
      plot_data <- DimPlot(seurat_obj, reduction = umap_reduc,
                           group.by = cluster_col,
                           pt.size  = input$pt_size,
                           label    = TRUE,
                           label.size = input$label_font_size) +
        NoLegend() + theme(axis.line = element_line(linewidth = 0.5))
      if (!is.null(cluster_colors)) {
        plot_data <- plot_data + scale_color_manual(values = cluster_colors)
      }
      interactive_plot <- ggplotly(plot_data, tooltip = "text") %>%
        layout(title = list(text = input$plot_title, font = list(size = 24)),
               hovermode = "closest")
      message("2D UMAP generated — reduction: ", umap_reduc, " — column: ",
              if (!is.null(cluster_col)) cluster_col else "Idents()")
      return(interactive_plot)
    }
  })

  # ---- Tables ----

  output$singler_summary_single <- renderDT({
    req(singler_result_single())
    render_annotation_summary_single(singler_result_single(), "singler")
  })

  # ---- Downloads ----

  output$download_singler_clusters_single <- downloadHandler(
    filename = function() paste0(getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset"),
                                 "_SingleR_cluster_summary_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      req(singler_result_single())
      res <- singler_result_single()
      out <- annotation_summary_single(res, "singler")
      out$reference <- res$params$reference
      out$label_column <- res$params$label_column
      out$mode <- res$params$mode
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$download_singler_cells_single <- downloadHandler(
    filename = function() paste0(getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset"),
                                 "_SingleR_per_cell_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      req(singler_result_single())
      write.csv(singler_result_single()$cells, file, row.names = FALSE)
    }
  )

  ############################## Find markers for a specific cluster ##############################

  # ---- Reactive state ----

  # Variables réactives pour cet onglet
  markers_global_single <- reactiveVal()

  umap_plot <- reactiveVal()

  # ---- Helpers ----

  calculate_markers_for_cluster <- function() {
    req(single_dataset_object(), input$target_cluster_global_single, input$min_pct_single, input$logfc_threshold_single)
    tryCatch({
      seurat_obj_global <- single_dataset_object()
      de_assay <- input$assay_de_single
      if (!de_assay %in% names(seurat_obj_global@assays)) stop(sprintf("Assay '%s' not found", de_assay))
      # Split layers exist only in Seurat v5 assays: JoinLayers has no method for SCT or v3 (Seurat v4) assays
      if (inherits(seurat_obj_global[[de_assay]], "Assay5")) seurat_obj_global <- JoinLayers(seurat_obj_global, assay = de_assay)
      # Selectors are filled from the active metadata column, so identities must come from the same column
      cluster_col <- active_cluster_column_single()
      cluster_vec <- if (!is.null(cluster_col)) as.character(seurat_obj_global@meta.data[[cluster_col]]) else as.character(Idents(seurat_obj_global))
      names(cluster_vec) <- colnames(seurat_obj_global)
      Idents(seurat_obj_global) <- factor(cluster_vec)
      if (!input$target_cluster_global_single %in% cluster_vec) {
        showNotification(paste0("Cluster '", input$target_cluster_global_single, "' not found in column ",
                                if (!is.null(cluster_col)) cluster_col else "Idents()"), type = "error")
        return(FALSE)
      }
      markers <- FindMarkers(seurat_obj_global, ident.1 = input$target_cluster_global_single, min.pct = input$min_pct_single,
                             logfc.threshold = input$logfc_threshold_single, assay = de_assay, return.thresh = 0.05)
      # Sort by absolute log2 fold change (most differential genes first)
      markers <- markers[order(abs(markers$avg_log2FC), decreasing = TRUE), ]
      # Store in global list for Venn diagrams BEFORE any modifications (rule #7)
      markers_copy_venn <- markers
      markers_copy_venn$gene <- rownames(markers_copy_venn)
      table_name <- paste0("SingleDataset_Cluster_", input$target_cluster_global_single, "_vs_All_", format(Sys.time(), "%H%M%S"))
      description <- paste0("Cluster ", input$target_cluster_global_single, " vs All Others")
      parameters <- list(min_pct = input$min_pct_single, logfc_threshold = input$logfc_threshold_single,
                         cluster_column = if (!is.null(cluster_col)) cluster_col else "Idents", assay = de_assay)
      result <- storeDETable(single_gene_table_storage(), markers_copy_venn, table_name, description, "cluster_group", parameters)
      if (result$success) single_gene_table_storage(result$storage)
      # Store raw numeric markers — formatting and filtering happen at render time
      markers$gene <- rownames(markers)
      markers_global_single(markers)
      return(TRUE)
    }, error = function(e) {
      showNotification(paste0("Error finding markers: ", e$message), type = "error")
      return(FALSE)
    })
  }

  # ---- Observers and computations ----

  # Update DE assay choices when object is loaded
  observe(run_safely("Interface update", persistent = FALSE, {
    req(single_dataset_object())
    available_assays <- names(single_dataset_object()@assays)
    updateSelectInput(session, "assay_de_single",
                      choices = available_assays,
                      selected = if ("RNA" %in% available_assays) "RNA" else available_assays[1])
  }))

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

  # Compute 3D UMAP from the same PCA dimensions as the 2D UMAP (Clustering step)
  observeEvent(input$compute_3d_umap, run_safely("compute 3d umap", {
    req(single_dataset_object())
    showModal(modalDialog(title = "Computing 3D UMAP", "Calculating 3D UMAP coordinates... This may take a moment.", easyClose = FALSE, footer = NULL))
    tryCatch({
      seurat_obj <- single_dataset_object()
      if (!"pca" %in% names(seurat_obj@reductions)) stop("Run PCA first (Scaling & PCA step)")
      n_dims <- if (!is.null(input$dimension_1) && !is.na(input$dimension_1)) input$dimension_1 else 30
      seurat_obj <- runUMAP3D_reproducible(object = seurat_obj, dims = 1:n_dims, reduction = "pca", seed = 42, verbose = TRUE)
      single_dataset_object(seurat_obj)
      removeModal()
      showNotification(sprintf("3D UMAP computed (PCA, %d dimensions): shown below the clustering plot.", n_dims), type = "message")
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error computing 3D UMAP:", conditionMessage(e)[1]), type = "error")
      message("Error in 3D UMAP computation: ", conditionMessage(e)[1])
    })
  }))

  # ---- UI and text outputs ----

  # Flag used by the UI to show the 3D panel only once a 3D UMAP exists
  output$has_umap3d_single <- reactive({
    obj <- single_dataset_object()
    !is.null(obj) && "umap3d" %in% names(obj@reductions)
  })

  outputOptions(output, "has_umap3d_single", suspendWhenHidden = FALSE)

  # ---- Plots ----

  # UMAP overview on the Cluster Comparison tab — respects both reduction selectors
  output$umap_plot <- renderPlot({
    req(single_dataset_object())
    seurat_obj  <- single_dataset_object()
    reduc_input <- if (!is.null(input$umap_reduction_single) && nzchar(input$umap_reduction_single)) {
      input$umap_reduction_single
    } else if (!is.null(input$cluster_reduction_single) && nzchar(input$cluster_reduction_single)) {
      input$cluster_reduction_single
    } else {
      NULL
    }
    umap_reduc <- if (!is.null(reduc_input) && reduc_input %in% names(seurat_obj@reductions)) {
      reduc_input
    } else {
      get_umap_reduction(seurat_obj)
    }
    if (is.null(umap_reduc)) return(NULL)
    cluster_col <- if (!is.null(input$umap_col_single) &&
                       nzchar(input$umap_col_single) &&
                       input$umap_col_single %in% colnames(seurat_obj@meta.data)) {
      input$umap_col_single
    } else {
      NULL
    }
    plot_data <- DimPlot(seurat_obj, reduction = umap_reduc,
                         group.by = cluster_col,
                         pt.size  = input$pt_size,
                         label    = input$show_labels,
                         label.size = input$label_font_size) +
      theme_void() + NoLegend() + ggtitle(NULL)
    cluster_colors <- get_cluster_colors(seurat_obj)
    if (!is.null(cluster_colors)) {
      plot_data <- plot_data + scale_color_manual(values = cluster_colors)
    }
    umap_plot(plot_data)
    plot_data
  })

  # Interactive 3D UMAP in the Clustering step, coloured by the current clusters
  output$clustering_umap3d_single <- renderPlotly({
    req(single_dataset_object())
    obj <- single_dataset_object()
    req("umap3d" %in% names(obj@reductions))
    build_umap3d_plotly(obj, cluster_values = as.character(Idents(obj)), cluster_colors = get_cluster_colors(obj),
                        pt_size = if (isTRUE(input$umap3d_pt_size_clust_single > 0)) input$umap3d_pt_size_clust_single else 2, label_size = 14,
                        hide_grid = isTRUE(input$umap3d_hide_grid_clust_single),
                        hide_axes = isTRUE(input$umap3d_hide_axes_clust_single),
                        dark_mode = isTRUE(input$umap3d_dark_clust_single))
  })

  # ---- Tables ----

  # Display results table — p-values stay numeric, formatted client-side by build_de_datatable
  output$table_global_single <- renderDataTable({
    req(markers_global_single())
    df <- markers_global_single()
    threshold <- input$pval_adj_filter_global_single
    if (!is.null(threshold) && !is.na(threshold)) {
      df <- df[df$p_val_adj <= threshold, , drop = FALSE]
    }
    cleaned <- clean_gene_names_for_html(df$gene)
    df$gene_link <- paste0('<a href="#" class="gene-name" data-gene="', cleaned, '">', cleaned, '</a>')
    df <- df[, c("gene_link", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
    colnames(df) <- c("Gene", "Log2FC", "Pct.1", "Pct.2", "P-value", "Adj. P-value")
    build_de_datatable(df, order_col = "Log2FC")
  })

  # ---- Downloads ----

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

  # ---- Reactive state ----

  # Initialize reactive variable
  markers_pairwise_single <- reactiveVal()

  # ---- Helpers ----

  # Function to calculate markers for pairwise comparison
  calculate_markers_for_comparison <- function() {
    req(single_dataset_object(), input$cluster1_pairwise_single, input$cluster2_pairwise_single,
        input$min_pct_pairwise_single, input$logfc_threshold_pairwise_single)
    tryCatch({
      # Check for overlap between groups
      if (any(input$cluster1_pairwise_single %in% input$cluster2_pairwise_single)) {
        showNotification("Cluster groups must be distinct", type = "error")
        return(FALSE)
      }
      seurat_obj <- single_dataset_object()
      de_assay <- input$assay_de_single
      if (!de_assay %in% names(seurat_obj@assays)) stop(sprintf("Assay '%s' not found", de_assay))
      # Selectors are filled from the active metadata column, so groups must be built from the same column
      cluster_col <- active_cluster_column_single()
      cluster_vec <- if (!is.null(cluster_col)) as.character(seurat_obj@meta.data[[cluster_col]]) else as.character(Idents(seurat_obj))
      new_idents <- cluster_vec
      new_idents[cluster_vec %in% input$cluster1_pairwise_single] <- "group1"
      new_idents[cluster_vec %in% input$cluster2_pairwise_single] <- "group2"
      n_group1 <- sum(new_idents == "group1", na.rm = TRUE)
      n_group2 <- sum(new_idents == "group2", na.rm = TRUE)
      if (n_group1 == 0 || n_group2 == 0) {
        showNotification(sprintf("No cells found for %s in column %s", if (n_group1 == 0) "group 1" else "group 2",
                                 if (!is.null(cluster_col)) cluster_col else "Idents()"), type = "error")
        return(FALSE)
      }
      names(new_idents) <- colnames(seurat_obj)
      Idents(seurat_obj) <- factor(new_idents)
      # Split layers exist only in Seurat v5 assays: JoinLayers has no method for SCT or v3 (Seurat v4) assays
      if (inherits(seurat_obj[[de_assay]], "Assay5")) seurat_obj <- JoinLayers(seurat_obj, assay = de_assay)
      markers <- FindMarkers(seurat_obj, ident.1 = "group1", ident.2 = "group2", min.pct = input$min_pct_pairwise_single,
                             logfc.threshold = input$logfc_threshold_pairwise_single, assay = de_assay, return.thresh = 0.05)
      # Sort by absolute log2 fold change
      markers <- markers[order(abs(markers$avg_log2FC), decreasing = TRUE), ]
      # Store in global list for Venn diagrams BEFORE any modifications (rule #7)
      markers_copy_venn <- markers
      markers_copy_venn$gene <- rownames(markers_copy_venn)
      group1_text <- paste(input$cluster1_pairwise_single, collapse = "_")
      group2_text <- paste(input$cluster2_pairwise_single, collapse = "_")
      table_name <- paste0("SingleDataset_Clusters_", group1_text, "_vs_", group2_text, "_", format(Sys.time(), "%H%M%S"))
      description <- paste0("Clusters [", group1_text, "] vs [", group2_text, "]")
      parameters <- list(min_pct = input$min_pct_pairwise_single, logfc_threshold = input$logfc_threshold_pairwise_single,
                         group1 = input$cluster1_pairwise_single, group2 = input$cluster2_pairwise_single,
                         cluster_column = if (!is.null(cluster_col)) cluster_col else "Idents", assay = de_assay)
      result <- storeDETable(single_gene_table_storage(), markers_copy_venn, table_name, description, "cluster_group", parameters)
      if (result$success) single_gene_table_storage(result$storage)
      # Store raw numeric markers with comparison label — formatting and filtering happen at render time
      markers$gene <- rownames(markers)
      markers$comparison <- sprintf("Group1(%s) vs Group2(%s)", paste(input$cluster1_pairwise_single, collapse = ","),
                                    paste(input$cluster2_pairwise_single, collapse = ","))
      markers_pairwise_single(markers)
      return(TRUE)
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      return(FALSE)
    })
  }

  # ---- Observers and computations ----

  # Update cluster choices
  observeEvent(input$cluster1_pairwise_single, run_safely("cluster1 pairwise single", {
    req(single_dataset_object())
    updateSelectInput(session, "cluster2_pairwise_single",
                      choices = setdiff(active_clusters_single(), input$cluster1_pairwise_single))
  }))

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

  # ---- Tables ----

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
    df <- df[, c("gene_link", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj", "comparison")]
    colnames(df) <- c("Gene", "Log2FC", "Pct.1", "Pct.2", "P-value", "Adj. P-value", "Comparison")
    build_de_datatable(df, order_col = "Log2FC")
  })

  # ---- Downloads ----

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

  ############################## Pathway enrichment ##############################

  # ---- Reactive state ----

  # ORA / GSEA on a new comparison (FindMarkers, logfc.threshold = 0) or ORA on a stored DE table; helpers in enrichment_functions.R
  enrichment_results_single <- reactiveVal(NULL)

  enrichment_gmt_single <- reactiveVal(NULL)

  # ---- Observers and computations ----

  # A new object invalidates the results computed on the previous one
  observeEvent(single_dataset_object(), run_safely("Enrichment reset", persistent = FALSE, { enrichment_results_single(NULL) }), ignoreInit = TRUE)

  # Grouping columns: metadata columns with 2 to 200 values
  observe(run_safely("Enrichment interface update", persistent = FALSE, {
    obj <- single_dataset_object(); req(obj)
    cols <- enrichment_grouping_columns(obj@meta.data)
    current <- isolate(input$enrichment_group_col_single)
    default <- if ("seurat_clusters" %in% cols) "seurat_clusters" else cols[1]
    updateSelectInput(session, "enrichment_group_col_single", choices = cols, selected = if (!is.null(current) && current %in% cols) current else default)
  }))

  # Group and restriction values follow the grouping column
  observeEvent(list(single_dataset_object(), input$enrichment_group_col_single), run_safely("Enrichment interface update", persistent = FALSE, {
    obj <- single_dataset_object(); req(obj, input$enrichment_group_col_single, input$enrichment_group_col_single %in% colnames(obj@meta.data))
    values <- enrichment_group_values(obj@meta.data[[input$enrichment_group_col_single]])
    updateSelectizeInput(session, "enrichment_group1_single", choices = values, selected = character(0))
    updateSelectizeInput(session, "enrichment_group2_single", choices = values, selected = character(0))
    # Restriction to clusters of the active identities (annotated names when the clusters were renamed)
    updateSelectizeInput(session, "enrichment_restrict_single", choices = levels(droplevels(as.factor(Idents(obj)))), selected = character(0))
  }))

  # Stored DE tables (same storage as the Venn diagrams)
  observe(run_safely("Enrichment interface update", persistent = FALSE, {
    updateSelectInput(session, "enrichment_stored_table_single", choices = stored_table_choices(single_gene_table_storage()))
  }))

  # GMT file read on upload and matched to the dataset genes; the md5 checksum identifies the file version in the parameters
  observeEvent(input$enrichment_gmt_file_single, run_safely("Gene set file", {
    obj <- single_dataset_object(); req(obj)
    sets <- match_gene_sets_to_dataset(read_gmt_gene_sets(input$enrichment_gmt_file_single$datapath), rownames(obj))
    if (sets$n_sets == 0) stop("No gene of the file matches the dataset genes")
    enrichment_gmt_single(c(sets, list(file = input$enrichment_gmt_file_single$name, md5 = unname(tools::md5sum(input$enrichment_gmt_file_single$datapath)))))
  }))

  observeEvent(input$enrichment_run_single, run_safely("Pathway enrichment", {
    obj <- single_dataset_object()
    if (is.null(obj)) stop("Load a dataset first")
    settings <- list(source = input$enrichment_source_single, group_col = input$enrichment_group_col_single, group1 = input$enrichment_group1_single,
                     group2 = input$enrichment_group2_single, restrict_values = input$enrichment_restrict_single, min_pct = input$enrichment_min_pct_single,
                     stored_table = input$enrichment_stored_table_single, gene_sets = input$enrichment_gene_sets_single, ontology = input$enrichment_ontology_single,
                     min_size = input$enrichment_min_size_single, max_size = input$enrichment_max_size_single, methods = input$enrichment_methods_single,
                     direction = input$enrichment_direction_single, ora_max_p_adj = input$enrichment_ora_max_p_adj_single, ora_min_log2fc = input$enrichment_ora_min_log2fc_single)
    showModal(modalDialog(title = "Pathway enrichment", "Testing gene sets...", easyClose = FALSE, footer = NULL))
    res <- tryCatch(run_enrichment_analysis(obj, settings, stored_tables = single_gene_table_storage(), gmt = enrichment_gmt_single()), finally = removeModal())
    enrichment_results_single(res)
    # Every run is stored with the DE tables, like the other analyses: saved with the session (Save & Session, where
    # unwanted runs can be deleted) and summarized with its parameters in the CSV / PDF export
    tables <- single_gene_table_storage()
    tables[[paste0("Enrichment_", format(Sys.time(), "%Y%m%d_%H%M%S"))]] <- enrichment_storage_entry(res, input$enrichment_show_max_p_adj_single %||% 0.05)
    single_gene_table_storage(tables)
    if (res$gsea_skipped) showNotification("GSEA needs all tested genes ranked: it runs only on a new comparison. ORA only for stored tables.", type = "warning", duration = 8)
    showNotification(paste0(enrichment_done_message(res, input$enrichment_show_max_p_adj_single %||% 0.05), ". Results stored (Save & Session)."), type = "message", duration = 8)
  }))

  enrichment_ora_plot_obj_single <- reactive({
    r <- enrichment_results_single(); shiny::validate(shiny::need(!is.null(r$ora), "Run ORA to see results"))
    enrichment_plot_or_message(r$ora$table, "ORA", input$enrichment_show_max_p_adj_single %||% 0.05, input$enrichment_ora_n_terms_single %||% 15, paste("ORA —", r$label), r$gene_set_label)
  })

  enrichment_gsea_plot_obj_single <- reactive({
    r <- enrichment_results_single(); shiny::validate(shiny::need(!is.null(r$gsea), "Run GSEA (new comparison) to see results"))
    enrichment_plot_or_message(r$gsea$table, "GSEA", input$enrichment_show_max_p_adj_single %||% 0.05, input$enrichment_gsea_n_terms_single %||% 15, paste("GSEA —", r$label), r$gene_set_label)
  })

  # ---- UI and text outputs ----

  output$enrichment_gmt_info_single <- renderText({
    s <- enrichment_gmt_single(); req(s)
    sprintf("%s: %d sets, %.0f%% of set genes found in the dataset", s$file, s$n_sets, 100 * s$fraction_matched)
  })

  output$enrichment_summary_single <- renderText({ r <- enrichment_results_single(); req(r); enrichment_summary_text(r) })

  # ---- Plots ----
  output$enrichment_ora_plot_single <- renderPlot({ enrichment_ora_plot_obj_single() }, res = 96)
  output$enrichment_gsea_plot_single <- renderPlot({ enrichment_gsea_plot_obj_single() }, res = 96)

  # ---- Tables ----
  output$enrichment_ora_table_single <- DT::renderDT({ r <- enrichment_results_single(); req(r$ora, nrow(r$ora$table) > 0); enrichment_datatable(r$ora$table, "ORA", input$enrichment_show_max_p_adj_single %||% 0.05, isTRUE(input$enrichment_show_all_single)) })
  output$enrichment_gsea_table_single <- DT::renderDT({ r <- enrichment_results_single(); req(r$gsea, nrow(r$gsea$table) > 0); enrichment_datatable(r$gsea$table, "GSEA", input$enrichment_show_max_p_adj_single %||% 0.05, isTRUE(input$enrichment_show_all_single)) })
  output$enrichment_params_table_single <- renderTable({ r <- enrichment_results_single(); req(r); r$params })

  # ---- Downloads ----

  output$download_enrichment_ora_plot_single <- downloadHandler(filename = function() paste0(enrichment_file_stub(enrichment_results_single()$label), "_ORA.png"),
                                                                content = function(file) ggplot2::ggsave(file, enrichment_ora_plot_obj_single(), width = 10, height = 7, dpi = 300, bg = "white"))

  output$download_enrichment_gsea_plot_single <- downloadHandler(filename = function() paste0(enrichment_file_stub(enrichment_results_single()$label), "_GSEA.png"),
                                                                 content = function(file) ggplot2::ggsave(file, enrichment_gsea_plot_obj_single(), width = 10, height = 7, dpi = 300, bg = "white"))

  # Result tables as Excel files: all tested terms sorted by signed significance, shaded red (up) / blue (down), parameters in a second sheet
  output$download_enrichment_ora_table_single <- downloadHandler(filename = function() paste0(enrichment_file_stub(enrichment_results_single()$label), "_ORA.xlsx"),
                                                                 content = function(file) { r <- enrichment_results_single(); req(r$ora); write_enrichment_xlsx(r, "ORA", file, input$enrichment_show_max_p_adj_single %||% 0.05) })

  output$download_enrichment_gsea_table_single <- downloadHandler(filename = function() paste0(enrichment_file_stub(enrichment_results_single()$label), "_GSEA.xlsx"),
                                                                  content = function(file) { r <- enrichment_results_single(); req(r$gsea); write_enrichment_xlsx(r, "GSEA", file, input$enrichment_show_max_p_adj_single %||% 0.05) })

  output$download_enrichment_params_single <- downloadHandler(filename = function() paste0(enrichment_file_stub(enrichment_results_single()$label), "_parameters.csv"),
                                                              content = function(file) { r <- enrichment_results_single(); req(r); utils::write.csv(r$params, file, row.names = FALSE) })

  #####################Exclusive Biomarkers Discovery###############################

  # ---- Reactive state ----

  # Initialize reactive variable for exclusive markers
  exclusive_markers_single <- reactiveVal()

  # ---- Observers and computations ----

  # Populate group-by column selector: any categorical metadata column
  observe(run_safely("Interface update", persistent = FALSE, {
    req(single_dataset_object())
    meta <- single_dataset_object()@meta.data
    categorical_cols <- colnames(meta)[sapply(meta, function(x) is.character(x) || is.factor(x))]
    cluster_col <- active_cluster_column_single()
    selected <- if (!is.null(cluster_col) && cluster_col %in% categorical_cols) cluster_col else categorical_cols[1]
    updateSelectInput(session, "exclusive_groupby_col_single", choices = categorical_cols, selected = selected)
  }))

  # Populate Target and Compare Against when group-by column changes
  observeEvent(input$exclusive_groupby_col_single, run_safely("exclusive groupby col single", {
    req(single_dataset_object(), nzchar(input$exclusive_groupby_col_single))
    col <- input$exclusive_groupby_col_single
    req(col %in% colnames(single_dataset_object()@meta.data))
    values <- sort(unique(as.character(single_dataset_object()@meta.data[[col]])))
    updateSelectInput(session, "exclusive_target_cluster_single", choices = values, selected = NULL)
    updateSelectInput(session, "exclusive_compare_clusters_single", choices = values, selected = NULL)
  }))

  # Keep Compare Against choices distinct from Target selection
  observeEvent(input$exclusive_target_cluster_single, run_safely("exclusive target cluster single", {
    req(single_dataset_object(), nzchar(input$exclusive_groupby_col_single))
    col <- input$exclusive_groupby_col_single
    req(col %in% colnames(single_dataset_object()@meta.data))
    values <- sort(unique(as.character(single_dataset_object()@meta.data[[col]])))
    updateSelectInput(session, "exclusive_compare_clusters_single",
                      choices = setdiff(values, input$exclusive_target_cluster_single))
  }), ignoreInit = TRUE)

  # Find exclusive biomarkers when button is clicked
  observeEvent(input$find_exclusive_markers_single, {
    req(single_dataset_object(), input$exclusive_groupby_col_single, input$exclusive_target_cluster_single)
    tryCatch({
      showModal(modalDialog(
        title = "Finding Exclusive Biomarkers",
        "Analyzing gene expression patterns across clusters...",
        easyClose = FALSE,
        footer = NULL
      ))
      seurat_obj  <- single_dataset_object()
      groupby_col <- input$exclusive_groupby_col_single
      target_groups <- input$exclusive_target_cluster_single
      compare_groups <- if (!is.null(input$exclusive_compare_clusters_single) && length(input$exclusive_compare_clusters_single) > 0) {
        input$exclusive_compare_clusters_single
      } else {
        setdiff(unique(as.character(seurat_obj@meta.data[[groupby_col]])), target_groups)
      }
      if (!groupby_col %in% colnames(seurat_obj@meta.data)) {
        showNotification(paste("Column", groupby_col, "not found in metadata!"), type = "error")
        removeModal()
        return()
      }
      if (any(target_groups %in% compare_groups)) {
        showNotification("Target and comparison groups must not overlap!", type = "error")
        removeModal()
        return()
      }
      both_groups   <- c(target_groups, compare_groups)
      cells_to_keep <- colnames(seurat_obj)[as.character(seurat_obj@meta.data[[groupby_col]]) %in% both_groups]
      if (length(cells_to_keep) == 0) {
        showNotification("No cells found for the selected groups!", type = "error")
        removeModal()
        return()
      }
      seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
      Idents(seurat_obj) <- factor(as.character(seurat_obj@meta.data[[groupby_col]]))
      message(paste0("Comparing ", groupby_col, " [", paste(target_groups, collapse = ", "),
                     "] vs [", paste(compare_groups, collapse = ", "), "]"))
      markers <- find_exclusive_biomarkers(
        seurat_obj = seurat_obj,
        target_cluster = target_groups,
        min_pct_target = input$exclusive_min_pct_target_single,
        max_pct_other = input$exclusive_max_pct_other_single,
        min_log2fc = input$exclusive_min_log2fc_single,
        detection_threshold = input$exclusive_detection_threshold_single,
        min_mean_expr_target = input$exclusive_min_mean_expr_single,
        statistical_test = input$exclusive_statistical_test_single,
        max_pvalue = input$exclusive_max_pvalue_single,
        pvalue_adjustment = "BH",
        assay_name = DefaultAssay(seurat_obj),
        top_n = NULL,
        verbose = TRUE
      )
      removeModal()
      if (nrow(markers) == 0) {
        showNotification("No exclusive biomarkers found with current parameters. Try relaxing the filters.",
                         type = "error", duration = 10)
        exclusive_markers_single(NULL)
        return()
      }
      # Copy whole p-value vectors (ifelse() with a scalar test returned only the first element, recycled to every gene)
      markers_copy_venn <- markers
      markers_copy_venn$p_val <- if ("pvalue" %in% colnames(markers_copy_venn)) markers_copy_venn$pvalue else NA_real_
      markers_copy_venn$p_val_adj <- if ("pvalue_adjusted" %in% colnames(markers_copy_venn)) markers_copy_venn$pvalue_adjusted else NA_real_
      markers_copy_venn$avg_log2FC <- markers_copy_venn$log2_fold_change
      table_name <- paste0("SingleDataset_ExclusiveMarkers_", groupby_col, "_",
                           paste(target_groups, collapse = "_"), "_vs_",
                           paste(compare_groups, collapse = "_"), "_", format(Sys.time(), "%H%M%S"))
      description <- paste0("Exclusive markers: ", groupby_col,
                            " [", paste(target_groups, collapse = ", "), "] vs [",
                            paste(compare_groups, collapse = ", "), "]")
      parameters <- list(
        groupby_col    = groupby_col,
        target_groups  = target_groups,
        compare_groups = compare_groups,
        min_pct_target = input$exclusive_min_pct_target_single,
        max_pct_other  = input$exclusive_max_pct_other_single,
        min_log2fc     = input$exclusive_min_log2fc_single,
        detection_threshold = input$exclusive_detection_threshold_single,
        statistical_test = input$exclusive_statistical_test_single
      )
      result <- storeDETable(single_gene_table_storage(), markers_copy_venn, table_name,
                             description, "exclusive_markers", parameters)
      if (result$success) {
        single_gene_table_storage(result$storage)
        showNotification(result$message, type = "message")
      } else {
        showNotification(result$message, type = "warning")
      }
      # p-values stay numeric in the display table: formatted client-side so column filters render sliders
      cleaned_gene_names <- clean_gene_names_for_html(markers$gene)
      markers$gene_link <- paste0('<a href="#" class="gene-name" data-gene="',
                                  cleaned_gene_names, '">', cleaned_gene_names, '</a>')
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
        type = "error", duration = 10
      )
      message("Error finding exclusive biomarkers: ", paste(e$message, collapse = " "))
    })
  })

  # ---- Tables ----

  # Display results table
  # Display results table — all numeric columns (incl. p-values) get range sliders
  output$table_exclusive_markers_single <- renderDataTable({
    req(exclusive_markers_single())
    df <- exclusive_markers_single()
    df <- df[!is.na(df$Gene) & df$Gene != "", ]
    spec_range <- range(as.numeric(df$Specificity), na.rm = TRUE, finite = TRUE)
    log2_range <- range(as.numeric(df$Log2FC), na.rm = TRUE, finite = TRUE)
    if (any(!is.finite(spec_range))) spec_range <- c(0, 1)
    if (any(!is.finite(log2_range))) log2_range <- c(0, 1)
    build_de_datatable(df, order_col = "Specificity", page_length = 15,
                       length_menu = c(15, 25, 50, 100, 200), gene_col_width = "120px") %>%
      formatStyle('Specificity', background = styleColorBar(spec_range, '#90EE90'), backgroundSize = '80% 90%',
                  backgroundRepeat = 'no-repeat', backgroundPosition = 'center') %>%
      formatStyle('Log2FC', background = styleColorBar(log2_range, '#ADD8E6'), backgroundSize = '80% 90%',
                  backgroundRepeat = 'no-repeat', backgroundPosition = 'center')
  })

  # ---- Downloads ----

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

  # ---- Reactive state ----

  # Reactive value to store the plot
  exclusive_plot_single <- reactiveVal()

  # ---- Observers and computations ----

  # Generate plot when button is clicked
  observeEvent(input$generate_exclusive_plot_single, {
    req(single_dataset_object())
    tryCatch({
      # Genes present in the default assay (case-insensitive), at most 10
      seurat_obj <- single_dataset_object()
      check <- validate_gene_input(input$exclusive_genes_to_plot_single, seurat_obj, DefaultAssay(seurat_obj))
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      genes_input <- check$genes
      if (length(genes_input) > 10) {
        showNotification("Maximum 10 genes for visualization: only the first 10 are plotted.", type = "warning")
        genes_input <- genes_input[1:10]
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
          reduction = get_umap_reduction(seurat_obj) %||% "pca",
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

  # ---- Plots ----

  # Render the plot
  output$plot_exclusive_markers_single <- renderPlot({
    req(exclusive_plot_single())
    print(exclusive_plot_single())
  }, height = 600)

  # ---- Downloads ----

  # Download handler for the plot
  output$download_exclusive_plot_single <- downloadHandler(
    filename = function() {
      object_name <- getObjectNameForDownload(single_dataset_object(), default_name = "ExclusiveMarkers")
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      genes_text <- substr(paste(parse_gene_tokens(input$exclusive_genes_to_plot_single), collapse = "_"), 1, 50)
      paste0(object_name, "_ExclusiveMarkers_", genes_text, "_",
             input$exclusive_plot_type_single, "_", timestamp, ".",
             input$exclusive_plot_format_single)
    },
    content = function(file) {
      req(exclusive_plot_single())

      tryCatch({
        plot_obj <- exclusive_plot_single()

        # Determine plot dimensions based on number of genes
        n_genes <- max(1, min(length(parse_gene_tokens(input$exclusive_genes_to_plot_single)), 10))

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

  # ---- Reactive state ----
  single_current_gene_lists <- reactiveVal(NULL)
  single_venn_plot_rendered <- reactiveVal(NULL)

  # ---- Observers and computations ----

  # Update select inputs (enrichment results left out)
  observe(run_safely("Interface update", persistent = FALSE, {
    updateVennSelectInputs(session, de_tables_only(single_gene_table_storage()), c("venn_table_1_single", "venn_table_2_single", "venn_table_3_single"))
  }))

  # Enable/disable generate button
  observe(run_safely("Interface update", persistent = FALSE, {
    tables <- single_gene_table_storage()
    if (length(tables) < 2) {
      shinyjs::disable("generate_venn_btn_single")
    } else {
      shinyjs::enable("generate_venn_btn_single")
    }
  }))

  # Generate Venn diagram
  observeEvent(input$generate_venn_btn_single, run_safely("generate venn btn single", {
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
  }))

  # ---- Plots ----

  # Render Venn diagram
  output$venn_plot_single <- renderPlot({
    req(single_venn_plot_rendered())
    grid.draw(single_venn_plot_rendered())
  })

  # ---- Tables ----

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

  # ---- Downloads ----

  output$download_venn_diagram_single <- createDownloadHandler(
    reactive_data = single_venn_plot_rendered,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "VennDiagram") }),
    data_name = "venn_comparison",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$venn_diagram_format_single }),
      width = 8,
      height = 6,
      dpi = reactive({ input$venn_diagram_dpi_single })
    )
  )

  output$download_venn_gene_lists_single <- createDownloadHandler(
    reactive_data = single_current_gene_lists,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "VennDiagram") }),
    data_name = "gene_lists",
    download_type = "ods"
  )

##########################CLuster Composition#########################

  # ---- Reactive state ----

  # Add this reactive value with other reactive values in single dataset server
  cluster_composition_single <- reactiveVal(NULL)

  # ---- Observers and computations ----

  # Observer for generating cluster composition table for single dataset
  observeEvent(input$generate_cluster_composition_single, {
    req(single_dataset_object())
    tryCatch({
      cluster_composition <- create_cluster_composition_table(
        single_dataset_object(),
        is_integrated = FALSE,
        cluster_col   = if (!is.null(input$umap_col_single) && nzchar(input$umap_col_single))
          input$umap_col_single else NULL
      )
      cluster_composition_single(cluster_composition)
    }, error = function(e) {
      showNotification(paste("Error:", conditionMessage(e)[1]), type = "error")
    })
  })

  # ---- Tables ----

  output$cluster_composition_single <- renderDT({
    req(cluster_composition_single())
    render_cluster_composition_table(cluster_composition_single(), is_integrated = FALSE)
  })

  # ---- Downloads ----

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

  # ---- Reactive state ----

  # Reactive values for storing results
  gene_coexpression_data_single <- reactiveVal(NULL)

  # Estimate each gene's threshold from its level in background cells and fill the threshold fields
  coexpr_background_estimate_single <- reactiveVal(NULL)

  # ---- Helpers ----

  # Thresholds of the validated genes, read from the per-token fields of the selected unit (any detection when a field is empty)
  coexpr_thresholds_single <- function(check) {
    unit <- input$coexpr_threshold_unit_single %||% "counts"
    default_value <- if (unit == "counts") 1 else 0
    values <- vapply(seq_along(check$tokens), function(i) {
      val <- input[[paste0("coexpr_thresh_single_", unit, "_", i)]]
      if (is.null(val) || is.na(val)) default_value else val
    }, numeric(1))
    align_gene_values(check, values)
  }

  # ---- Observers and computations ----

  # Plot rebuilt from the stored results whenever the visualization selector changes (no new analysis needed)
  coexpression_plot_single <- reactive({
    res <- gene_coexpression_data_single()
    req(res, length(res$genes_analyzed) >= 2)
    create_coexpression_plot(res$data, res$genes_analyzed, plot_type = input$coexpr_plot_type %||% "bar")
  })

  # Background clusters for the ambient-RNA threshold estimate (same identities as the DE cluster selector)
  observe(run_safely("Interface update", persistent = FALSE, {
    req(single_dataset_object())
    updateSelectInput(session, "coexpr_background_clusters_single", choices = sort(as.character(unique(Idents(single_dataset_object())))))
  }))

  observeEvent(input$estimate_coexpr_thresholds_single, run_safely("Co-expression background thresholds (single dataset)", {
    req(single_dataset_object())
    seurat_obj <- single_dataset_object()
    background <- input$coexpr_background_clusters_single
    if (length(background) == 0) stop("Select the cluster(s) that should not express these genes")
    assay <- DefaultAssay(seurat_obj)
    check <- validate_gene_input(input$gene_text_coexpression_single, seurat_obj, assay)
    if (!check$ok) stop(check$message)
    unit <- input$coexpr_threshold_unit_single %||% "counts"
    cells <- colnames(seurat_obj)[as.character(Idents(seurat_obj)) %in% background]
    est <- estimate_coexpression_thresholds(seurat_obj, check$genes, cells, assay_name = assay, threshold_layer = unit, quantile = 0.95)
    for (i in seq_along(check$tokens)) {
      g <- check$matched[i]
      if (!is.na(g)) updateNumericInput(session, paste0("coexpr_thresh_single_", unit, "_", i), value = est$threshold[est$gene == g])
    }
    coexpr_background_estimate_single(list(table = est, clusters = background, unit = unit, n_cells = length(cells)))
  }))

  # Populate group-by column and filter groups for co-expression
  observe(run_safely("Interface update", persistent = FALSE, {
    req(single_dataset_object())
    meta <- single_dataset_object()@meta.data
    categorical_cols <- colnames(meta)[sapply(meta, function(x) is.character(x) || is.factor(x))]
    cluster_col <- active_cluster_column_single()
    selected <- if (!is.null(cluster_col) && cluster_col %in% categorical_cols) cluster_col else categorical_cols[1]
    updateSelectInput(session, "coexpr_groupby_col_single", choices = categorical_cols, selected = selected)
  }))

  observeEvent(input$coexpr_groupby_col_single, run_safely("coexpr groupby col single", {
    req(single_dataset_object(), nzchar(input$coexpr_groupby_col_single))
    col <- input$coexpr_groupby_col_single
    req(col %in% colnames(single_dataset_object()@meta.data))
    values <- sort(unique(as.character(single_dataset_object()@meta.data[[col]])))
    updateSelectInput(session, "coexpr_filter_groups_single", choices = values, selected = NULL)
  }))

  observeEvent(input$analyze_coexpression_single, run_safely("Co-expression (single dataset)", {
    req(single_dataset_object())
    seurat_obj <- single_dataset_object()
    check <- validate_gene_input(input$gene_text_coexpression_single, seurat_obj, DefaultAssay(seurat_obj))
    if (!check$ok) stop(check$message)
    if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
    thresholds <- coexpr_thresholds_single(check)
    genes_input <- names(thresholds)
    if (length(genes_input) < 2) stop("Co-expression needs at least 2 different genes found in the dataset")
    # Restrict to the selected groups of the chosen metadata column, if any
    groupby_col <- input$coexpr_groupby_col_single
    filter_groups <- input$coexpr_filter_groups_single
    if (length(filter_groups) > 0 && !is.null(groupby_col) && nzchar(groupby_col) && groupby_col %in% colnames(seurat_obj@meta.data)) {
      cells_to_keep <- colnames(seurat_obj)[as.character(seurat_obj@meta.data[[groupby_col]]) %in% filter_groups]
      if (length(cells_to_keep) == 0) stop("No cells found for the selected groups")
      seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
      message(sprintf("Co-expression restricted to %s [%s]: %d cells", groupby_col, paste(filter_groups, collapse = ", "), ncol(seurat_obj)))
    }
    coexpr_results <- analyze_gene_coexpression(seurat_obj = seurat_obj, genes = genes_input, assay_name = DefaultAssay(seurat_obj),
                                                expression_thresholds = thresholds, is_integrated = FALSE,
                                                threshold_layer = input$coexpr_threshold_unit_single %||% "counts")
    gene_coexpression_data_single(coexpr_results)
    showNotification("Co-expression analysis completed", type = "message")
  }))

  # ---- UI and text outputs ----

  # Per-gene threshold fields: one numericInput per token (field i <-> token i); separate IDs per unit so a UMI count is
  # never reused as a log-normalized value; typed values survive re-renders
  output$coexpr_thresholds_ui_single <- renderUI({
    genes <- parse_gene_tokens(input$gene_text_coexpression_single)
    if (length(genes) == 0) return(helpText("Enter gene names to set individual thresholds"))
    unit <- input$coexpr_threshold_unit_single %||% "counts"
    tagList(lapply(seq_along(genes), function(i) {
      input_id <- paste0("coexpr_thresh_single_", unit, "_", i)
      if (unit == "counts") {
        numericInput(input_id, label = paste0(genes[i], ": min. UMIs per cell"), value = isolate(input[[input_id]]) %||% 1, min = 1, max = 100, step = 1)
      } else {
        numericInput(input_id, label = paste0(genes[i], ": normalized expression >"), value = isolate(input[[input_id]]) %||% 0, min = 0, max = 10, step = 0.1)
      }
    }))
  })

  output$coexpr_background_info_single <- renderUI({
    r <- coexpr_background_estimate_single()
    if (is.null(r)) return(NULL)
    lines <- sprintf("%s: background 95th percentile = %s -> %s (%.1f%% of background cells still positive)", r$table$gene,
                     r$table$background_quantile_value, format_threshold_label(r$table$gene, r$table$threshold, r$unit), r$table$pct_background_positive)
    tags$div(style = "font-size: 11px; color: #555; margin-top: 6px;",
             tags$p(sprintf("Background: %s (%d cells)", paste(r$clusters, collapse = ", "), r$n_cells), style = "margin: 0 0 2px 0; font-weight: bold;"),
             lapply(lines, function(l) tags$p(l, style = "margin: 0;")))
  })

  # ---- Plots ----

  output$gene_coexpression_plot_single <- renderPlot({
    req(coexpression_plot_single())
    coexpression_plot_single()
  }, height = 600)

  # ---- Tables ----

  # Render results using modular functions
  output$gene_coexpression_table_single <- renderDT({
    req(gene_coexpression_data_single())
    render_coexpression_table(gene_coexpression_data_single(), "gene_coexpression_table_single")
  })

  # ---- Downloads ----

  # Download handlers for coexpression analysis
  output$download_coexpression_table_single <- createDownloadHandler(
    reactive_data = reactive({
      res <- gene_coexpression_data_single()
      req(res)
      res$data[, !grepl("_Bar$|^Coexpression_Visual$", colnames(res$data)), drop = FALSE]
    }),
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "Coexpression") }),
    data_name = "coexpression_analysis",
    download_type = "csv"
  )

  # Download of the co-expression plot (the tab has no format / resolution selector: PNG at 300 dpi)
  output$download_coexpression_plot_single <- createDownloadHandler(
    reactive_data = coexpression_plot_single,
    object_name_reactive = reactive({
      getObjectNameForDownload(single_dataset_object(), default_name = "Coexpression")
    }),
    data_name = "coexpression_plot",
    download_type = "plot",
    plot_params = list(
      file_type = reactive(input$format_coexpr_single %||% "png"),
      width = 10,
      height = 8,
      dpi = reactive(input$dpi_coexpr_single %||% 300)
    )
  )

  #####################DE: Expressors vs Rest - Single Dataset###############################

  # ---- Reactive state ----
  coexpr_de_results_single <- reactiveVal(NULL)

  # ---- Observers and computations ----

  # Populate cluster selector when object is available
  observe({
    req(single_dataset_object())
    tryCatch({
      clusters <- sort(as.character(unique(Idents(single_dataset_object()))))
      updateSelectInput(session, "coexpr_de_cluster_single", choices = clusters)
    }, error = function(e) NULL)
  })

  observeEvent(input$run_coexpr_de_single, {
    req(single_dataset_object())
    tryCatch({
      check <- validate_gene_input(input$gene_text_coexpression_single, single_dataset_object(), DefaultAssay(single_dataset_object()))
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      selected_clusters <- input$coexpr_de_cluster_single
      if (is.null(selected_clusters) || length(selected_clusters) == 0) {
        showNotification("Select at least one cluster", type = "warning")
        return()
      }
      de_thresholds <- coexpr_thresholds_single(check)
      genes_input <- names(de_thresholds)
      showNotification("Running within-cluster DE...", type = "message", duration = 3)
      result <- run_coexpression_de(
        seurat_obj            = single_dataset_object(),
        genes                 = genes_input,
        clusters              = selected_clusters,
        assay_name            = DefaultAssay(single_dataset_object()),
        expression_thresholds = de_thresholds,
        threshold_layer       = input$coexpr_threshold_unit_single %||% "counts"
      )
      coexpr_de_results_single(result)
      label <- if (length(genes_input) == 1) "expressors" else "co-expressors"
      showNotification(
        sprintf("DE complete — %d %s vs %d others, %d DE genes (%s)",
                result$n_coexpressors, label, result$n_non_coexpressors,
                nrow(result$markers), result$cluster_label),
        type = "message", duration = 5
      )
    }, error = function(e) {
      showNotification(paste("DE error:", conditionMessage(e)[1]), type = "error", duration = 8)
      message("coexpr_de single error: ", conditionMessage(e)[1])
    })
  })

  # ---- UI and text outputs ----

  output$coexpr_de_info_single <- renderUI({
    r <- coexpr_de_results_single()
    if (is.null(r)) return(NULL)
    bd <- r$breakdown
    n_genes <- length(r$genes_used)
    if (n_genes == 1) {
      breakdown_txt <- ""
      pos_label <- paste0(r$genes_used, "+")
      neg_label <- paste0(r$genes_used, "-")
    } else if (!is.null(bd$only_gene1)) {
      breakdown_txt <- sprintf(" (only %s: %d | only %s: %d | neither: %d)",
                               r$genes_used[1], bd$only_gene1,
                               r$genes_used[2], bd$only_gene2, bd$neither)
      pos_label <- "Co-expressors (all)"
      neg_label <- "Non-co-expressors"
    } else {
      breakdown_txt <- sprintf(" (partial: %d | none: %d)", bd$partial, bd$none)
      pos_label <- "Co-expressors (all)"
      neg_label <- "Non-co-expressors"
    }
    # Thresholds in the unit used for the run: "Pax7 >= 2 UMI" or "Pax7 > 0.5"
    thresh_txt <- paste(format_threshold_label(names(r$thresholds_used), r$thresholds_used, r$threshold_layer %||% "data"), collapse = ", ")
    tags$div(
      tags$p(sprintf("Cluster(s): %s | %s | %d DE genes", r$cluster_label, thresh_txt, nrow(r$markers)),
             style = "color: #aaa; font-size: 12px; margin-bottom: 2px;"),
      tags$p(sprintf("✓ %s: %d cells  |  ✗ %s: %d cells%s",
                     pos_label, r$n_coexpressors, neg_label, r$n_non_coexpressors, breakdown_txt),
             style = "color: #aaa; font-size: 12px; margin: 0;")
    )
  })

  # ---- Tables ----

  output$coexpr_de_table_single <- renderDT({
    req(coexpr_de_results_single())
    m <- coexpr_de_results_single()$markers
    m$avg_log2FC <- round(m$avg_log2FC, 3)
    m$pct.1      <- round(m$pct.1, 3)
    m$pct.2      <- round(m$pct.2, 3)
    m$p_val_adj  <- formatC(m$p_val_adj, format = "e", digits = 2)
    m$p_val      <- formatC(m$p_val,     format = "e", digits = 2)
    DT::datatable(
      m[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "p_val")],
      colnames  = c("Gene", "log2FC", "pct.expr", "pct.non", "p_adj", "p_val"),
      rownames  = FALSE, filter = "top", extensions = "Buttons",
      options   = list(pageLength = 25, dom = "Bfrtip",
                       buttons = list(
                         list(extend = "csv",   filename = paste0("ExprDE_", format(Sys.time(), "%Y%m%d")), text = "CSV"),
                         list(extend = "excel", filename = paste0("ExprDE_", format(Sys.time(), "%Y%m%d")), text = "Excel")
                       ))
    ) |> DT::formatStyle("avg_log2FC", backgroundColor = DT::styleInterval(0, c("#c0392b22", "#1a572622")))
  })

  # ---- Downloads ----

  output$download_coexpr_de_single <- downloadHandler(
    filename = function() paste0("ExprDE_", paste(input$coexpr_de_cluster_single, collapse = "-"),
                                 "_", format(Sys.time(), "%Y%m%d"), ".csv"),
    content  = function(file) {
      req(coexpr_de_results_single())
      write.csv(coexpr_de_results_single()$markers, file, row.names = FALSE)
    }
  )

  ########################## Volcano Plot ##########################

  # ---- Reactive state ----

  # Reactive value for storing volcano plot
  volcano_plot_single <- reactiveVal(NULL)

  # ---- Observers and computations ----

  # Update volcano plot source selector when DEG results are available (enrichment results left out)
  observe(run_safely("Interface update", persistent = FALSE, {
    tables <- de_tables_only(single_gene_table_storage())
    if (length(tables) > 0) {
      table_choices <- c("None" = "", names(tables))
      updateSelectInput(session, "volcano_deg_source_single", choices = table_choices)
    } else {
      updateSelectInput(session, "volcano_deg_source_single", choices = c("None" = ""))
    }
  }))

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

  # ---- Plots ----

  # Render volcano plot
  output$volcano_plot_single <- renderPlot({
    req(volcano_plot_single())
    volcano_plot_single()
  })

  # ---- Downloads ----

  # Download volcano plot
  output$download_volcano_plot_single <- createDownloadHandler(
    reactive_data = volcano_plot_single,
    object_name_reactive = reactive({
      getObjectNameForDownload(single_dataset_object(), default_name = "VolcanoPlot")
    }),
    data_name = "volcano_plot",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$volcano_plot_format_single }),
      width = 10,
      height = 8,
      dpi = reactive({ input$volcano_plot_dpi_single })
    )
  )

  ############################## Subseting ##############################

  # ---- Reactive state ----

  # Variables wich stock the subset of single_dataset_object()
  subset_seurat <- reactiveVal()

  # ---- Observers and computations ----

  observe(run_safely("Interface update", persistent = FALSE, {
    { subset_seurat(single_dataset_object())
    }
  }))

  observe(run_safely("Interface update", persistent = FALSE, {
    if (!is.null(single_dataset_object())) {
      updateSelectInput(session, "select_ident_subset", choices = unique(Idents(single_dataset_object())))
    }
  }))

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
      # Every gene must exist: the "number of genes to express" rule depends on the full list
      check <- validate_gene_input(input$gene_list_subset, single_dataset_object(), DefaultAssay(single_dataset_object()))
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (length(check$missing) > 0) {
        showNotification(sprintf("Genes not found in assay '%s': %s", check$assay, paste(check$missing, collapse = ", ")), type = "error", duration = 8)
        return()
      }
      if (!is.null(check$message)) showNotification(check$message, type = "message", duration = 6)
      gene_list <- check$genes

      message("=== GENE SUBSET ANALYSIS ===")
      message("Genes requested: ", paste(gene_list, collapse = ", "))
      message("Active assay: ", DefaultAssay(single_dataset_object()))
      message("Expression threshold: ", input$expression_threshold)
      message("Negative selection: ", input$negative_gene_subset)

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

  # ---- Plots ----

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

  # ---- Downloads ----

  # Downloading Seurat subset object
  output$download_subset_seurat <- createDownloadHandler(
    reactive_data = subset_seurat,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "SeuratSubset") }),
    data_name = "seurat_subset",
    download_type = "seurat",
    show_modal = TRUE
  )

  # Download analysis parameters as CSV
  output$download_params_single <- downloadHandler(
    filename = function() {
      obj_name <- getObjectNameForDownload(single_dataset_object(), default_name = "object")
      paste0(obj_name, "_analysis_parameters_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      # Session parameters + one row per stored analysis (DE tables, exclusive markers, enrichment runs)
      params_df <- rbind(
        collect_analysis_params(
          module     = "single",
          qc_stats   = qc_stats_single(),
          seurat_obj = single_dataset_object()
        ),
        collect_stored_analyses(single_gene_table_storage())
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
      # Session parameters + one row per stored analysis (DE tables, exclusive markers, enrichment runs)
      params_df <- rbind(
        collect_analysis_params(
          module     = "single",
          qc_stats   = qc_stats_single(),
          seurat_obj = single_dataset_object()
        ),
        collect_stored_analyses(single_gene_table_storage())
      )
      result <- generate_pdf_report(params_df, file, module = "single")
      if (!result$success) {
        message("[Report] PDF generation failed:\n", result$message)
        # Surface the Python error in the interface, not only in the console
        showNotification(paste("PDF report generation failed:", result$message), type = "error", duration = NULL)
        stop("PDF generation failed")
      }
    }
  )

  ############################## Session: save, contents and stored tables ##############################

  # ---- Reactive state ----

  # IDs captured when the confirmation opens, so a selection change in between cannot delete other tables
  pending_table_deletion_single <- reactiveVal(character(0))

  # ---- Observers and computations ----

  # Seurat object + session state: every "Download Seurat Object" button saves the full analysis
  session_object_single <- reactive({
    req(single_dataset_object())
    pack_session_state(single_dataset_object(), gene_tables = single_gene_table_storage(),
                       qc_stats = qc_stats_single(), module = "single")
  })

  # Stored result tables, listed by their human-readable description
  session_tables_view_single <- reactive({
    summarize_gene_tables(single_gene_table_storage())
  })

  observeEvent(input$delete_session_tables_single, run_safely("delete session tables single", {
    selected <- input$session_tables_single_rows_selected
    if (length(selected) == 0) {
      showNotification("Select at least one table to delete.", type = "warning")
      return()
    }
    view <- session_tables_view_single()
    pending_table_deletion_single(view$ID[selected])
    showModal(modalDialog(
      title = "Delete stored tables",
      tags$p(sprintf("Delete %d table(s)? Venn diagrams and volcano plots will no longer be able to use them.", length(selected))),
      tags$ul(lapply(view$Table[selected], tags$li)),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_delete_session_tables_single", "Delete", class = "btn-danger")
      )
    ))
  }))

  observeEvent(input$confirm_delete_session_tables_single, run_safely("confirm delete session tables single", {
    ids <- pending_table_deletion_single()
    tables <- single_gene_table_storage()
    tables[intersect(ids, names(tables))] <- NULL
    single_gene_table_storage(tables)
    pending_table_deletion_single(character(0))
    removeModal()
    showNotification(sprintf("%d table(s) deleted", length(ids)), type = "message")
  }))

  # ---- UI and text outputs ----

  output$session_size_info_single <- renderText({
    tables <- single_gene_table_storage()
    sprintf("%d stored table(s), %s in total. Plots are never stored.",
            length(tables), format(object.size(tables), units = "auto"))
  })

  # ---- Tables ----

  # What the object contains
  output$session_object_summary_single <- renderTable({
    req(single_dataset_object())
    summarize_seurat_object(single_dataset_object())
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$session_tables_single <- DT::renderDT({
    DT::datatable(
      session_tables_view_single(),
      rownames  = FALSE,
      selection = "multiple",
      options   = list(dom = "tip", pageLength = 10,
                       # ID column is hidden: used only to delete the right entries
                       columnDefs = list(list(visible = FALSE, targets = 0)))
    )
  })

  # ---- Downloads ----

  output$download_session_single <- createDownloadHandler(
    reactive_data = session_object_single,
    object_name_reactive = reactive({
      getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset")
    }),
    data_name = "CellHub_session",
    download_type = "seurat",
    show_modal = TRUE
  )
}

