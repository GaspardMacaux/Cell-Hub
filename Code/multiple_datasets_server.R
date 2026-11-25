
############################## Script Multiple Datasets Analysis ##############################

# multiple_datasets_server.R

multiple_datasets_server <- function(input, output, session) {  

  ############################## Shared reactive state ##############################
  # Objects used by several sections below

  # Variables for the merge dataset part
  multiple_datasets_object <- reactiveVal()

  # Register in the session-wide registry so other modules (CellChat, trajectory) can import it
  session$userData$multiple_datasets_object <- multiple_datasets_object

  seurat_objects <- reactiveValues()
  data_loaded <- reactiveValues()
  rv_metadata <- reactiveValues(num_fields = 1)
  registered_observers <- new.env(parent = emptyenv())  # tracks which slot observers exist

  # Analysis parameter tracking for CSV export
  qc_stats_multi          <- reactiveVal(NULL)

  integration_method_used <- reactiveVal(NULL)  # NULL = method read from the object at export (extract_integration_info) 
  harmony_vars_used       <- reactiveVal(NULL)

  # og system pour le panneau permanent
  log_system <- createLogSystem()

  feature_plot_merge <- reactiveVal()
  dot_plot_merge <- reactiveVal()
  vln_plot_merge <- reactiveVal()
  ridge_plot_merge <- reactiveVal()

  # Historique très simple
  cluster_history_merge <- reactiveValues(
    states = list(),
    max_history = 10
  )

  forced_col_merge <- reactiveVal(NULL)

  # Reactive variables for this tab
  diff_genes_compare <- reactiveVal(NULL)          # Stores DE genes for a cluster vs all others

  diff_genes_compare_cluster <- reactiveVal(NULL)  # Stores DE genes for group1 vs group2
  diff_genes_compare_datasets <- reactiveVal(NULL) # Stores DE genes for dataset1 vs dataset2
  gene_table_storage <- reactiveVal(list())        # Stores all DE gene tables for Venn diagrams
  current_gene_lists <- reactiveVal(NULL)          # Stores current lists for Venn diagrams

  # Reactive value to store cluster composition data
  cluster_composition_multiple <- reactiveVal(NULL)

  
  
  ############################## Download buttons state ##############################
  # Each download button stays disabled until the result it exports exists
  bind_download_state("downloadUMAP_merge", function() req(clustering_plot_merge()))
  bind_download_state("downloadVlnPlotMerge", function() vln_plot_merge())
  bind_download_state("downloadFeaturePlotMerge", function() feature_plot_merge())
  bind_download_state("downloadDotPlotMerge", function() dot_plot_merge())
  bind_download_state("downloadRidgePlotMerge", function() ridge_plot_merge())
  bind_download_state("download_genes_number_expression_merge", function() number_of_nuclei_merge())
  bind_download_state("download_heatmap_multi", function() heatmap_plot_multidataset())
  bind_download_state("download_heatmap_markers_multi", function() { markers_df <- heatmap_markers_data_multi(); req(markers_df) })
  bind_download_state("download_scatter_multi", function() scatter_plot_multi())
  bind_download_state("download_transfer_clusters_merge", function() req(transfer_result_merge()))
  bind_download_state("download_transfer_cells_merge", function() req(transfer_result_merge()))
  bind_download_state("download_singler_clusters_merge", function() req(singler_result_merge()))
  bind_download_state("download_singler_cells_merge", function() req(singler_result_merge()))
  bind_download_state("download_filtered_umap_plot", function() filtered_umap_plot())
  bind_download_state("download_markers_single_cluster_merge", function() { df <- diff_genes_compare(); req(df) })
  bind_download_state("download_markers_multiple_clusters_merge", function() { df <- diff_genes_compare_cluster(); req(df) })
  bind_download_state("download_diff_dataset_cluster", function() diff_genes_compare_datasets())
  bind_download_state("download_cluster_composition_multiple", function() cluster_composition_multiple())
  bind_download_state("download_coexpression_table_multiple", function() gene_coexpression_data_multiple())
  bind_download_state("download_coexpression_plot_multiple", function() { res <- gene_coexpression_data_multiple(); req(res, length(res$genes_analyzed) >= 2) })
  bind_download_state("download_coexpr_de_multiple", function() req(coexpr_de_results_multiple()))
  bind_download_state("download_exclusive_markers_merge", function() { markers <- exclusive_markers_merge(); req(markers) })
  bind_download_state("download_enrichment_ora_plot_merge", function() { r <- enrichment_results_merge(); req(!is.null(r$ora)) })
  bind_download_state("download_enrichment_gsea_plot_merge", function() { r <- enrichment_results_merge(); req(!is.null(r$gsea)) })
  bind_download_state("download_enrichment_ora_table_merge", function() { r <- enrichment_results_merge(); req(r$ora) })
  bind_download_state("download_enrichment_gsea_table_merge", function() { r <- enrichment_results_merge(); req(r$gsea) })
  bind_download_state("download_exclusive_plot_merge", function() req(exclusive_plot_merge()))
  bind_download_state("download_venn_diagram", function() venn_plot_rendered())
  bind_download_state("download_venn_gene_lists", function() current_gene_lists())
  bind_download_state("download_pie_chart_merge", function() cluster_pie_plot_merge())
  bind_download_state("download_volcano_plot_merge", function() volcano_plot_merge())
  bind_download_state("download_subset_merge", function() subset_seurat_merge())
  bind_download_state("download_session_merge", function() req(multiple_datasets_object()))
  
  ############################## Loading Data ##############################

  # ---- Reactive state ----
  merged_gene_tables <- reactiveValues()

  # log system pour le modal d'intégration
  modal_log_system <- createLogSystem()

  # ---- Observers and computations ----
  shinyjs::disable("add_field")
  shinyjs::disable("add_metadata")

  # Initialize metadata split column with default
  observe(run_safely("Interface update", persistent = FALSE, {
    if (!is.null(multiple_datasets_object())) {
      if (is.null(input$metadata_split_column)) {
        # Set initial value
        if ("integration_group" %in% colnames(multiple_datasets_object()@meta.data)) {
          updateSelectInput(session, "metadata_split_column", selected = "integration_group")
        } else if ("dataset" %in% colnames(multiple_datasets_object()@meta.data)) {
          updateSelectInput(session, "metadata_split_column", selected = "dataset")
        }
      }
    }
  }))

  # ---- UI and text outputs ----

  output$datasets_loaded <- reactive({
    !is.null(multiple_datasets_object())
  })

  output$merge_loading_logs <- renderPrint({
    cat(log_system$log_value())
  })

  output$merge_preintegrated_logs <- renderPrint({
    cat(log_system$log_value())
  })

  ############################## Force Memory Cleanup ##############################

  # ---- Helpers ----

  # Processing individual datasets using modular functions
  observe_file_input <- function(index) {
    key <- paste0("obs_", index)
    if (isTRUE(registered_observers[[key]])) return()
    registered_observers[[key]] <- TRUE
    observeEvent(input[[paste0("merge", index)]], run_safely("Action", {
      req(input[[paste0("merge", index)]])
      # Reset this slot so a re-upload always triggers a fresh load
      seurat_objects[[paste0("seurat_object", index)]] <- NULL
      data_loaded[[paste0("loaded", index)]] <- FALSE
      log_system$add_log("", "info")
      log_system$add_log(paste("=== LOADING DATASET", index, "==="), "info")
      withProgress(message = paste("Processing dataset", index), value = 0, {
        tryCatch({
          file_input <- input[[paste0("merge", index)]]
          file_extension <- logFileInfo(
            file_name    = file_input$name,
            file_path    = file_input$datapath,
            log_function = log_system$add_log
          )
          qc_params    <- extractQCParams(input)
          dataset_type <- input[[paste0("dataset_type_merge", index)]]
          dataset_name <- paste0("Dataset_", index)
          log_system$add_log(paste("Dataset type:", dataset_type), "info")
          log_system$add_log(paste("Placeholder name:", dataset_name, "(final name applied at integration)"), "info")
          # Validate file type: ZIP, H5 (Cell Ranger) or H5AD for raw types, RDS or H5AD for processed
          if (dataset_type %in% c("snRNA_merge", "multiome_merge")) {
            if (!file_extension %in% c("zip", "h5", "hdf5", "h5ad")) {
              log_system$add_log("Error: Expected ZIP, H5 or H5AD file for raw data", "error")
              showNotification(paste("Dataset", index, ": Please upload a ZIP (10X folder), H5 (Cell Ranger) or H5AD file"), type = "error")
              return(NULL)
            }
          } else if (dataset_type == "seurat_object_merge") {
            if (!file_extension %in% c("rds", "h5ad")) {
              log_system$add_log("Error: Expected RDS or H5AD file for processed object", "error")
              showNotification(paste("Dataset", index, ": Please upload an RDS or H5AD file"), type = "error")
              return(NULL)
            }
          }
          seurat_object <- preprocessRawDataset(
            file_path            = file_input$datapath,
            dataset_type         = dataset_type,
            species              = input$species_choice_merge,
            dataset_name         = dataset_name,
            qc_params            = qc_params,
            normalization_method = input$normalization_method_merge,
            log_function         = log_system$add_log
          )
          # Memory optimization if requested
          if (isTRUE(input$optimize_memory_merge %in% c("slim", "aggressive"))) {
            log_system$add_log(paste("Optimizing memory for dataset", index, "..."), "info")
            seurat_object <- slimSeuratObject(seurat_object, level = input$optimize_memory_merge, log_function = log_system$add_log)
          }
          logMemoryUsage(log_system$add_log, seurat_object, paste("Dataset", index))
          seurat_objects[[paste0("seurat_object", index)]] <- seurat_object
          data_loaded[[paste0("loaded", index)]] <- TRUE
          log_system$add_log(paste("=== DATASET", index, "READY FOR INTEGRATION ==="), "success")
          message(paste("Dataset", index, "processed and stored successfully"))
          showNotification(
            paste0("Dataset ", index, " processed: ",
                   format(ncol(seurat_object), big.mark = ","), " cells"),
            type = "message"
          )
          shinyjs::enable("add_field")
          shinyjs::enable("add_metadata")
        }, error = function(e) {
          log_system$add_log(paste("ERROR:", conditionMessage(e)[1]), "error")
          message(paste("Error in dataset", index, ":", conditionMessage(e)[1]))
          showNotification(paste("Error processing dataset", index, ":", conditionMessage(e)[1]), type = "error")
        })
      })
    }), ignoreInit = TRUE)
  }

  # Gene picker of the visualization tab, refilled after every new object (group_by_select is filled by OBSERVER 1)
  updateUIElements <- function() {
    if (is.null(multiple_datasets_object())) return(invisible(NULL))
    gene_list <- reactive_gene_list_merge()
    updatePickerInput(session, "geneInput_merge", choices = gene_list)
    message(paste("UI elements updated successfully - loaded", length(gene_list) - 1, "genes"))
    invisible(NULL)
  }

  # ---- Observers and computations ----

  observeEvent(input$force_cleanup_merge, {
    tryCatch({
      showModal(modalDialog(
        title = "Cleaning Memory",
        tags$div(
          style = "text-align: center;",
          tags$h4(icon("broom"), " Freeing memory..."),
          tags$p("Please wait a few seconds"),
          tags$br(),
          tags$div(class = "progress progress-striped active",
                   tags$div(class = "progress-bar progress-bar-warning", 
                            role = "progressbar",
                            style = "width: 100%"))
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      
      message("=== FORCING AGGRESSIVE MEMORY CLEANUP ===")
      
      # Get RAM before cleanup
      if (.Platform$OS.type == "windows") {
        ram_before <- memory.size() / 1024
        message(paste("RAM before cleanup:", round(ram_before, 1), "GB"))
      }
      
      # Force multiple garbage collections
      gc()
      Sys.sleep(0.5)  # Give system time to cleanup
      gc()
      Sys.sleep(0.5)
      gc()
      
      # Clean temporary files
      temp_files <- list.files(tempdir(), full.names = TRUE, pattern = "\\.rds$")
      if (length(temp_files) > 0) {
        file.remove(temp_files)
        message(paste("Removed", length(temp_files), "temporary .rds files"))
      }
      
      # Get RAM after cleanup
      if (.Platform$OS.type == "windows") {
        ram_after <- memory.size() / 1024
        freed <- ram_before - ram_after
        message(paste("RAM after cleanup:", round(ram_after, 1), "GB"))
        message(paste("Memory freed:", round(freed, 1), "GB"))
        
        Sys.sleep(1)  # Let user see the modal
        removeModal()
        
        showNotification(
          HTML(paste0(
            "<strong>✓ Memory cleaned successfully!</strong><br>",
            "Freed: ", round(freed, 1), " GB<br>",
            "Current usage: ", round(ram_after, 1), " GB / ", 
            round(memory.limit() / 1024, 1), " GB"
          )),
          type = "message",
          duration = 8
        )
      } else {
        # For Mac/Linux (can't easily check memory)
        Sys.sleep(1)
        removeModal()
        
        showNotification(
          "✓ Memory cleanup completed!",
          type = "message",
          duration = 5
        )
      }
      
      message("=== CLEANUP COMPLETED ===")
      
    }, error = function(e) {
      removeModal()
      showNotification(
        paste("Cleanup error:", e$message),
        type = "error",
        duration = 5
      )
    })
  })

  # Load pre-integrated Seurat object
  observeEvent(input$load_seurat_file_merge, {
    tryCatch({
      file_extension <- tolower(tools::file_ext(input$load_seurat_file_merge$name))
      if (!file_extension %in% c("rds", "h5ad")) {
        showNotification("Please upload a .rds or .h5ad file", type = "error")
        return(NULL)
      }
      loaded_seurat <- loadSeuratObject(
        rds_path = input$load_seurat_file_merge$datapath,
        add_dataset_column = TRUE,
        dataset_name = tools::file_path_sans_ext(basename(input$load_seurat_file_merge$name)),
        clean_before = TRUE,
        module_type = "multiple",
        file_format = file_extension,
        log_function = log_system$add_log
      )
      # Cluster names held only in Idents are otherwise invisible to metadata-based selectors
      loaded_seurat <- materializeIdentsAsColumn(loaded_seurat, column_name = "ClusterIdents", log_function = log_system$add_log)
      # Memory optimization if requested
      if (isTRUE(input$optimize_memory_preintegrated %in% c("slim", "aggressive"))) {
        log_system$add_log("Applying memory optimization...", "info")
        loaded_seurat <- slimSeuratObject(loaded_seurat, level = input$optimize_memory_preintegrated, log_function = log_system$add_log)
      }
      logMemoryUsage(log_system$add_log, loaded_seurat, "Pre-integrated object")
      # Restore the saved session (stored tables + QC statistics); without one, clear the previous
      # object's tables so Venn diagrams never mix results from two different objects
      session_state <- unpack_session_state(loaded_seurat)
      gene_table_storage(session_state$gene_tables)
      qc_stats_multi(session_state$qc_stats)
      # Integration method of a previous run in this session does not describe this object: read it from the object at export
      integration_method_used(NULL)
      harmony_vars_used(NULL)
      if (session_state$found) {
        log_system$add_log(sprintf("Session restored: %d stored table(s)%s",
                                   length(session_state$gene_tables),
                                   if (!is.null(session_state$qc_stats)) ", QC statistics" else ""), "success")
      }
      multiple_datasets_object(loaded_seurat)
      init_clustering_plot_merge(loaded_seurat)
      updateUIElements()
      shinyjs::enable("add_field")
      shinyjs::enable("add_metadata")
      showNotification(
        if (session_state$found) sprintf("Pre-integrated object loaded successfully! Session restored (%d stored table(s)).", length(session_state$gene_tables))
        else "Pre-integrated object loaded successfully!",
        type = "message"
      )
    }, error = function(e) {
      showNotification(paste("Error loading object:", conditionMessage(e)[1]), type = "error")
    })
  })

  observeEvent(input$simple_merge, {
    tryCatch({
      showModal(modalDialog(
        title = "Processing Datasets",
        "Merging datasets without integration...",
        easyClose = FALSE,
        footer = NULL
      ))
      # Apply dataset names from UI at integration time (user may have changed them after upload)
      seurat_objects <- applyDatasetNamesFromInputs(seurat_objects, input$num_datasets, input)
      # QC thresholds, statistics and normalization of each dataset, kept before merge() drops them (saved with the session)
      qc_stats_multi(build_merge_qc_stats(lapply(seq_len(input$num_datasets), function(i) seurat_objects[[paste0("seurat_object", i)]])))
      # Collect all processed Seurat objects
      seurat_list <- list()
      for (i in 1:input$num_datasets) {
        seurat_obj <- seurat_objects[[paste0("seurat_object", i)]]
        if (is.null(seurat_obj)) {
          stop(paste("Dataset", i, "is not processed yet"))
        }
        seurat_list[[i]] <- seurat_obj
      }
      # Use modular integration function
      merged_object <- performDataIntegration(
        seurat_list = seurat_list,
        integration_method = "simple"
      )
      # Store integrated object
      multiple_datasets_object(merged_object)
      init_clustering_plot_merge(merged_object)
      # Method read from the object at export (merge_method column); clears the value of a previous run
      integration_method_used(NULL)
      harmony_vars_used(NULL)
      removeModal()
      updateUIElements()
      shinyjs::enable("add_field")
      shinyjs::enable("add_metadata")
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error during merge:", e$message), type = "error")
    })
  })

  # fastMNN Integration Handler
  observeEvent(input$integrate_fastmnn, {
    tryCatch({
      showModal(modalDialog(
        title = "fastMNN Integration",
        tags$div(
          style = "text-align: center;",
          tags$h4("Running fastMNN integration..."),
          tags$p("Method: batchelor::fastMNN on log-normalized counts"),
          tags$p("Skipping ScaleData/PCA — fastMNN runs its own multiBatchPCA"),
          tags$br(),
          tags$div(class = "progress progress-striped active",
                   tags$div(class = "progress-bar progress-bar-info", 
                            role = "progressbar",
                            style = "width: 100%"))
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      seurat_objects <- applyDatasetNamesFromInputs(seurat_objects, input$num_datasets, input)
      # QC thresholds, statistics and normalization of each dataset, kept before merge() drops them (saved with the session)
      qc_stats_multi(build_merge_qc_stats(lapply(seq_len(input$num_datasets), function(i) seurat_objects[[paste0("seurat_object", i)]])))
      seurat_list <- list()
      for (i in 1:input$num_datasets) {
        seurat_obj <- seurat_objects[[paste0("seurat_object", i)]]
        if (is.null(seurat_obj)) {
          stop(paste("Dataset", i, "is not processed yet"))
        }
        seurat_list[[i]] <- seurat_obj
      }
      integrated_object <- performDataIntegration(
        seurat_list = seurat_list,
        integration_method = "fastmnn",
        mnn_dims = 50,
        mnn_k = 20
      )
      # Memory optimization after integration
      if (isTRUE(input$optimize_memory_merge == "slim")) {
        log_system$add_log("Applying post-fastMNN memory optimization...", "info")
        integrated_object <- slimSeuratObject(integrated_object, keep_counts = TRUE,
                                              keep_scale_data = FALSE,
                                              log_function = log_system$add_log)
      }
      multiple_datasets_object(integrated_object)
      init_clustering_plot_merge(integrated_object)
      # Same label as extract_integration_info(), used by the parameter export and the PDF report
      integration_method_used("fastMNN")
      harmony_vars_used(NULL)
      logMemoryUsage(log_system$add_log, integrated_object, "FastMNN integrated object")
      shinyjs::enable("add_field")
      shinyjs::enable("add_metadata")
      removeModal()
      # Gene picker of the visualization tab (otherwise empty after fastMNN)
      updateUIElements()
      cleanupIntegrationMemory()
      showNotification(
        paste("fastMNN integration completed successfully!",
              ncol(integrated_object), "cells integrated"),
        type = "message",
        duration = 10
      )
    }, error = function(e) {
      removeModal()
      showNotification(
        paste("fastMNN integration failed:", conditionMessage(e)[1], 
              "\nTry Simple Merge or check your data."), 
        type = "error",
        duration = 15
      )
    })
  })
  
  # RPCA Integration Handler (Seurat v5 reciprocal PCA, split by the 'dataset' column)
  observeEvent(input$integrate_rpca, {
    tryCatch({
      showModal(modalDialog(
        title = "RPCA Integration",
        tags$div(
          style = "text-align: center;",
          tags$h4("Running RPCA integration..."),
          tags$p("Method: Seurat IntegrateLayers + RPCAIntegration (30 dimensions)"),
          tags$br(),
          tags$div(class = "progress progress-striped active",
                   tags$div(class = "progress-bar progress-bar-danger", role = "progressbar", style = "width: 100%"))
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      # Apply dataset names from UI at integration time (user may have changed them after upload)
      seurat_objects <- applyDatasetNamesFromInputs(seurat_objects, input$num_datasets, input)
      # QC thresholds, statistics and normalization of each dataset, kept before merge() drops them (saved with the session)
      qc_stats_multi(build_merge_qc_stats(lapply(seq_len(input$num_datasets), function(i) seurat_objects[[paste0("seurat_object", i)]])))
      seurat_list <- list()
      for (i in 1:input$num_datasets) {
        seurat_obj <- seurat_objects[[paste0("seurat_object", i)]]
        if (is.null(seurat_obj)) stop(paste("Dataset", i, "is not processed yet"))
        seurat_list[[i]] <- seurat_obj
      }
      integrated_object <- performDataIntegration(seurat_list = seurat_list, integration_method = "rpca", rpca_dims = 30, rpca_k_anchor = 5)
      # Memory optimization after integration
      if (isTRUE(input$optimize_memory_merge == "slim")) {
        log_system$add_log("Applying post-RPCA memory optimization...", "info")
        integrated_object <- slimSeuratObject(integrated_object, keep_counts = TRUE, keep_scale_data = FALSE, log_function = log_system$add_log)
      }
      multiple_datasets_object(integrated_object)
      init_clustering_plot_merge(integrated_object)
      # Same label as extract_integration_info(), used by the parameter export and the PDF report
      integration_method_used("RPCA")
      harmony_vars_used(NULL)
      logMemoryUsage(log_system$add_log, integrated_object, "RPCA integrated object")
      shinyjs::enable("add_field")
      shinyjs::enable("add_metadata")
      removeModal()
      # Gene picker of the visualization tab
      updateUIElements()
      cleanupIntegrationMemory()
      showNotification(paste("RPCA integration completed successfully!", ncol(integrated_object), "cells integrated"), type = "message", duration = 10)
    }, error = function(e) {
      removeModal()
      showNotification(paste("RPCA integration failed:", conditionMessage(e)[1], "\nTry Harmony or Simple Merge, or check your data."), type = "error", duration = 15)
    })
  })

  # Modal for dataset loading
  observeEvent(input$open_file_input_modal, run_safely("open file input modal", {
    # Reset all slot states so a new load cycle starts clean
    for (i in seq_len(20)) {
      data_loaded[[paste0("loaded", i)]]           <- NULL
      seurat_objects[[paste0("seurat_object", i)]] <- NULL
    }
    
    showModal(modalDialog(
      title = NULL,
      size = "l",
      easyClose = TRUE,
      footer = NULL,
      
      # Header with close button
      tags$div(
        style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
             padding: 15px 20px; margin: -15px -15px 20px -15px; 
             border-radius: 5px 5px 0 0; color: white;",
        tags$div(
          style = "display: flex; justify-content: space-between; align-items: center;",
          tags$h3("Load Multiple Datasets", style = "margin: 0; font-weight: bold;"),
          tags$button(
            type = "button",
            class = "close",
            "data-dismiss" = "modal",
            style = "color: white; opacity: 1; font-size: 28px;",
            "×"
          )
        )
      ),
      
      tagList(
        # Dataset configuration box
        tags$div(
          style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; 
               border: 2px solid #dee2e6; margin-bottom: 15px;",
          tags$h5("Dataset Configuration",
                  style = "color: #495057; margin-top: 0; margin-bottom: 10px; font-weight: bold;"),
          fluidRow(
            column(4, numericInput('num_datasets', 'Number of datasets:', value = 2, min = 1)),
            column(4, selectInput("species_choice_merge", "Species:",
                                  choices = c("Mouse" = "mouse", "Human" = "human", "Rat" = "rat"),
                                  selected = "mouse")),
            column(4, selectInput("normalization_method_merge", "Normalization Method:",
                                  choices = c("LogNormalize (standard)" = "LogNormalize",
                                              "SCTransform (regularized negative binomial)" = "SCTransform",
                                              "CLR (centered log-ratio)" = "CLR",
                                              "RC (relative counts)" = "RC"),
                                  selected = "LogNormalize"))
          ),
          fluidRow(
            column(12,
                   tags$div(
                     style = "margin-top: 10px; background-color: white; padding: 10px; border-radius: 6px; border-left: 4px solid #11998e;",
                     radioButtons("optimize_memory_merge", NULL,
                                  choices = list(
                                    "Safe (recommended)" = "slim",
                                    "Aggressive" = "aggressive",
                                    "Keep full object" = "full"
                                  ),
                                  selected = "slim", inline = TRUE),
                     tags$p(icon("info-circle"),
                            " Safe: removes RNA scale.data only (recomputed before every PCA). Aggressive: also removes SCT and integrated assays, whose corrected values are then lost.",
                            style = "font-size: 11px; color: #6c757d; margin: 0;")
                   )
            )
          )
        ),
        
        # QC Parameters box
        tags$div(
          style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; 
               border: 2px solid #dee2e6; margin-bottom: 15px;",
          tags$h5("Quality Control Parameters",
                  style = "color: #495057; margin-top: 0; margin-bottom: 10px; font-weight: bold;"),
          fluidRow(
            column(4, numericInput("min_features_merge", "Min features/cell:", value = 200, min = 0)),
            column(4, numericInput("max_features_merge", "Max features/cell:", value = 3500, min = 0)),
            column(4, numericInput("max_mt_percent_merge", "Max MT %:", value = 5, min = 0, max = 100))
          ),
          tags$p("These parameters filter cells during processing.",
                 style = "font-size: 11px; color: #6c757d; margin-bottom: 0;")
        ),
        
        # Dynamic file inputs box
        tags$div(
          style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; 
               border: 2px solid #dee2e6; margin-bottom: 15px;",
          tags$h5("Dataset Files",
                  style = "color: #495057; margin-top: 0; margin-bottom: 10px; font-weight: bold;"),
          uiOutput("fileInputs")
        ),
        
        # Integration Method
        tags$div(
          style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; 
               border: 2px solid #dee2e6; margin-bottom: 15px;",
          tags$h5("Integration Method",
                  style = "color: #495057; margin-top: 0; margin-bottom: 10px; font-weight: bold;"),
          
          div(
            style = "display: grid; grid-template-columns: repeat(5, minmax(0, 1fr)); gap: 10px;",
            # Standard Integration
            div(
              style = "flex: 1 1 160px; min-width: 0; background-color: #e3f2fd; padding: 12px; border-radius: 6px; 
                   border: 2px solid #2196F3; min-height: 200px; display: flex; 
                   flex-direction: column;",
              tags$div(
                style = "flex: 1;",
                tags$div(
                  style = "text-align: center; margin-bottom: 8px;",
                  tags$h6(icon("project-diagram"), strong(" Standard"),
                          style = "color: #1976D2; margin: 0; font-size: 14px;")
                ),
                tags$p(style = "font-size: 11px; margin: 0 0 8px 0; text-align: center; color: #555;",
                       "Best for 2-3 datasets"),
                tags$div(
                  style = "font-size: 10px; margin-bottom: 10px; line-height: 1.4;",
                  "✓ Highest quality", tags$br(),
                  "✓ Strong batch correction", tags$br(),
                  "⚠ Slow (5-15 min)", tags$br(),
                  tags$span("⚠ May crash >4 datasets", style = "color: #d32f2f;")
                )
              ),
              tags$div(
                style = "margin-top: auto;",
                actionButton('integrate', 'Run Standard Integration',
                             class = 'btn-primary', icon = icon("play"),
                             disabled = TRUE,
                             style = "width: 100%; font-size: 11px; font-weight: bold; white-space: normal;")
              )
            ),
            
            # Harmony Integration
            div(
              style = "flex: 1 1 160px; min-width: 0; background-color: #e8f5e9; padding: 12px; border-radius: 6px; 
                   border: 2px solid #4CAF50; min-height: 200px; display: flex; 
                   flex-direction: column;",
              tags$div(
                style = "flex: 1;",
                tags$div(
                  style = "text-align: center; margin-bottom: 8px;",
                  tags$h6(icon("bolt"), strong(" Harmony"),
                          style = "color: #2E7D32; margin: 0; font-size: 14px;")
                ),
                tags$p(style = "font-size: 11px; margin: 0 0 8px 0; text-align: center; color: #555;",
                       "Best for 4+ datasets"),
                tags$div(
                  style = "font-size: 10px; margin-bottom: 10px; line-height: 1.4;",
                  "✓ Fast (2-5 min)", tags$br(),
                  "✓ Scalable 50+ datasets", tags$br(),
                  "✓ Stable & reliable", tags$br(),
                  tags$span(icon("thumbs-up"), " Recommended",
                            style = "color: #2E7D32; font-weight: bold;")
                )
              ),
              tags$div(
                style = "margin-top: auto;",
                # Batch variable corrected by Harmony (choices filled by the observer when the first dataset is loaded)
                selectInput("harmony_vars_integration", "Batch variable:", choices = NULL, width = "100%"),
                actionButton('integrate_harmony', 'Run Harmony Integration',
                             class = 'btn-success', icon = icon("play"),
                             disabled = TRUE,
                             style = "width: 100%; font-size: 11px; font-weight: bold; white-space: normal;")
              )
            ),
            
            # fastMNN Integration
            div(
              style = "flex: 1 1 160px; min-width: 0; background-color: #ede7f6; padding: 12px; border-radius: 6px; 
       border: 2px solid #673AB7; min-height: 200px; display: flex; 
       flex-direction: column;",
              tags$div(
                style = "flex: 1;",
                tags$div(
                  style = "text-align: center; margin-bottom: 8px;",
                  tags$h6(
                    icon("link"), 
                    strong(" fastMNN"),
                    style = "color: #4527A0; margin: 0; font-size: 14px;"
                  )
                ),
                tags$p(
                  style = "font-size: 11px; margin: 0 0 8px 0; text-align: center; color: #555;",
                  "Mutual Nearest Neighbors"
                ),
                tags$div(
                  style = "font-size: 10px; margin-bottom: 10px; line-height: 1.4;",
                  "✓ Fast (3-8 min)", tags$br(),
                  "✓ Skips ScaleData/PCA", tags$br(),
                  "✓ Robust to composition shifts", tags$br(),
                  tags$span("Good for technical batches", style = "color: #4527A0;")
                )
              ),
              tags$div(
                style = "margin-top: auto;",
                actionButton(
                  'integrate_fastmnn',
                  'Run fastMNN Integration',
                  class = 'btn-info',
                  icon = icon("play"),
                  disabled = TRUE,
                  style = "width: 100%; font-size: 11px; font-weight: bold; white-space: normal;"
                )
              )
            ),
            # RPCA Integration
            div(
              style = "flex: 1 1 160px; min-width: 0; background-color: #fce4ec; padding: 12px; border-radius: 6px; 
                   border: 2px solid #E91E63; min-height: 200px; display: flex; 
                   flex-direction: column;",
              tags$div(
                style = "flex: 1;",
                tags$div(
                  style = "text-align: center; margin-bottom: 8px;",
                  tags$h6(icon("exchange-alt"), strong(" RPCA"),
                          style = "color: #AD1457; margin: 0; font-size: 14px;")
                ),
                tags$p(style = "font-size: 11px; margin: 0 0 8px 0; text-align: center; color: #555;",
                       "Reciprocal PCA (Seurat)"),
                tags$div(
                  style = "font-size: 10px; margin-bottom: 10px; line-height: 1.4;",
                  "✓ Faster than Standard", tags$br(),
                  "✓ Less memory than Standard", tags$br(),
                  "✓ Conservative correction", tags$br(),
                  tags$span("Good when cell types differ", style = "color: #AD1457;")
                )
              ),
              tags$div(
                style = "margin-top: auto;",
                actionButton('integrate_rpca', 'Run RPCA Integration',
                             class = 'btn-danger', icon = icon("play"),
                             disabled = TRUE,
                             style = "width: 100%; font-size: 11px; font-weight: bold; white-space: normal;")
              )
            ),
            # Simple Merge
            div(
              style = "flex: 1 1 160px; min-width: 0; background-color: #fff3e0; padding: 12px; border-radius: 6px; 
                   border: 2px solid #FF9800; min-height: 200px; display: flex; 
                   flex-direction: column;",
              tags$div(
                style = "flex: 1;",
                tags$div(
                  style = "text-align: center; margin-bottom: 8px;",
                  tags$h6(icon("layer-group"), strong(" Simple Merge"),
                          style = "color: #E65100; margin: 0; font-size: 14px;")
                ),
                tags$p(style = "font-size: 11px; margin: 0 0 8px 0; text-align: center; color: #555;",
                       "No batch correction"),
                tags$div(
                  style = "font-size: 10px; margin-bottom: 10px; line-height: 1.4;",
                  "✓ Instant (<1 min)", tags$br(),
                  "✓ Preserves all data", tags$br(),
                  "⚠ No batch correction", tags$br(),
                  tags$span("Only for similar experiments", style = "color: #E65100;")
                )
              ),
              tags$div(
                style = "margin-top: auto;",
                actionButton('simple_merge', 'Merge Without Integration',
                             class = 'btn-warning', icon = icon("play"),
                             disabled = TRUE,
                             style = "width: 100%; font-size: 11px; font-weight: bold; white-space: normal;")
              )
            )
          )
        ),
        
        # Dataset Loading Logs
        tags$div(
          style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; 
               border: 2px solid #dee2e6; margin-bottom: 15px;",
          tags$h5(icon("terminal"), " Dataset Loading Logs",
                  style = "color: #495057; margin-top: 0; margin-bottom: 10px; font-weight: bold;"),
          tags$pre(
            style = "background-color: #1e272e; color: #00ff00; padding: 15px; 
                 border-radius: 6px; font-family: 'Courier New', monospace; 
                 font-size: 10px; max-height: 300px; overflow-y: auto; 
                 margin: 0; white-space: pre-wrap;",
            verbatimTextOutput("merge_loading_logs")
          ),
          tags$p(icon("info-circle"), " Monitor dataset loading progress here. Each dataset will be logged as it's uploaded.",
                 style = "font-size: 11px; color: #6c757d; margin: 10px 0 0 0;")
        )
      )
    ))
    
    # Dynamic file inputs generation
    observeEvent(input$num_datasets, run_safely("num datasets", {
      req(input$num_datasets > 0)
      output$fileInputs <- renderUI({
        lapply(1:input$num_datasets, function(i) {
          fluidRow(
            column(4, textInput(paste0('dataset_name', i), paste0('Dataset ', i, ' name:'),
                                value = paste0("Dataset_", i))),
            column(4, selectInput(paste0("dataset_type_merge", i), "Type:",
                                  choices = list("snRNA-seq" = "snRNA_merge",
                                                 "Multiome" = "multiome_merge",
                                                 "Seurat Object" = "seurat_object_merge"))),
            column(4, fileInput(paste0('merge', i), 'File:',
                                accept = c('.rds', '.zip', '.h5', '.hdf5', '.h5ad', 'application/x-hdf5')))
          )
        })
      })
    }))
  }))

  # Setup file input observers
  observeEvent(input$num_datasets, run_safely("num datasets", {
    lapply(1:input$num_datasets, function(i) {
      observe_file_input(i)
    })
  }))

  # Enable/disable all three buttons based on file status
  observe(run_safely("Interface update", persistent = FALSE, {
    req(input$num_datasets)
    
    # Check file processing status
    files_status <- sapply(1:input$num_datasets, function(i) {
      dataset_type <- input[[paste0("dataset_type_merge", i)]]
      file_input <- input[[paste0("merge", i)]]
      is_processed <- !is.null(data_loaded[[paste0("loaded", i)]]) && 
        data_loaded[[paste0("loaded", i)]]
      
      if (!is.null(dataset_type)) {
        if (dataset_type == "seurat_object_merge") {
          return(!is.null(file_input))
        } else {
          return(is_processed)
        }
      }
      return(FALSE)
    })
    
    all_loaded <- all(files_status)
    enough_files <- length(files_status) >= 2
    should_enable <- all_loaded && enough_files
    
    # Enable every integration button
    shinyjs::toggleState("integrate", condition = should_enable)  
    shinyjs::toggleState("integrate_harmony", condition = should_enable)  
    shinyjs::toggleState("simple_merge", condition = should_enable)  
    shinyjs::toggleState("integrate_fastmnn", condition = should_enable)
    shinyjs::toggleState("integrate_rpca", condition = should_enable)
  }))

  # Standard integration using modular functions
  observeEvent(input$integrate, {
    tryCatch({
      showModal(modalDialog(
        title = "Please Wait",
        "Integrating datasets...",
        easyClose = FALSE,
        footer = NULL
      ))
      # Apply dataset names from UI at integration time (user may have changed them after upload)
      seurat_objects <- applyDatasetNamesFromInputs(seurat_objects, input$num_datasets, input)
      # QC thresholds, statistics and normalization of each dataset, kept before merge() drops them (saved with the session)
      qc_stats_multi(build_merge_qc_stats(lapply(seq_len(input$num_datasets), function(i) seurat_objects[[paste0("seurat_object", i)]])))
      # Collect all processed Seurat objects
      seurat_list <- list()
      for (i in 1:input$num_datasets) {
        seurat_obj <- seurat_objects[[paste0("seurat_object", i)]]
        if (is.null(seurat_obj)) {
          stop(paste("Dataset", i, "is not processed yet"))
        }
        seurat_list[[i]] <- seurat_obj
      }
      # Use modular integration function
      integrated_object <- performDataIntegration(
        seurat_list = seurat_list,
        integration_method = "standard"
      )
      # Memory optimization after integration
      if (isTRUE(input$optimize_memory_merge == "slim")) {
        log_system$add_log("Applying post-integration memory optimization...", "info")
        integrated_object <- slimSeuratObject(integrated_object, keep_counts = TRUE,
                                              keep_scale_data = FALSE,
                                              log_function = log_system$add_log)
      }
      # Store integrated object
      multiple_datasets_object(integrated_object)
      init_clustering_plot_merge(integrated_object)
      # Method read from the object at export (anchor integration); clears the value of a previous run
      integration_method_used(NULL)
      harmony_vars_used(NULL)
      logMemoryUsage(log_system$add_log, integrated_object, "Integrated object")
      # Enable UI controls
      shinyjs::enable("add_field")
      shinyjs::enable("add_metadata")
      removeModal()
      updateUIElements()
      # Clean up individual datasets + temp dirs
      cleanupIntegrationMemory()
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error during integration:", e$message), type = "error")
    })
  })

  # Harmony Integration Handler 
  observeEvent(input$integrate_harmony, {
    tryCatch({
      showModal(modalDialog(
        title = "Harmony Integration",
        tags$div(
          style = "text-align: center;",
          tags$h4("Running Harmony integration..."),
          tags$p("Correcting for: ", input$harmony_vars_integration %||% "dataset"),
          tags$p("Using 30 dimensions"),
          tags$br(),
          tags$div(class = "progress progress-striped active",
                   tags$div(class = "progress-bar progress-bar-success", 
                            role = "progressbar",
                            style = "width: 100%"))
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      # Apply dataset names from UI at integration time (user may have changed them after upload)
      seurat_objects <- applyDatasetNamesFromInputs(seurat_objects, input$num_datasets, input)
      # QC thresholds, statistics and normalization of each dataset, kept before merge() drops them (saved with the session)
      qc_stats_multi(build_merge_qc_stats(lapply(seq_len(input$num_datasets), function(i) seurat_objects[[paste0("seurat_object", i)]])))
      # Collect all processed Seurat objects
      seurat_list <- list()
      for (i in 1:input$num_datasets) {
        seurat_obj <- seurat_objects[[paste0("seurat_object", i)]]
        if (is.null(seurat_obj)) {
          stop(paste("Dataset", i, "is not processed yet"))
        }
        seurat_list[[i]] <- seurat_obj
      }
      harmony_vars_to_use <- if (!is.null(input$harmony_vars_integration) && 
                                 input$harmony_vars_integration != "") {
        input$harmony_vars_integration
      } else {
        if ("integration_group" %in% colnames(seurat_list[[1]]@meta.data)) {
          "integration_group"
        } else {
          "dataset"
        }
      }
      message(paste("✓ Using batch correction variable:", harmony_vars_to_use))
      integrated_object <- performDataIntegration(
        seurat_list = seurat_list,
        integration_method = "harmony",
        harmony_vars = harmony_vars_to_use,
        harmony_dims = 30
      )
      # Memory optimization after integration
      if (isTRUE(input$optimize_memory_merge == "slim")) {
        log_system$add_log("Applying post-Harmony memory optimization...", "info")
        integrated_object <- slimSeuratObject(integrated_object, keep_counts = TRUE,
                                              keep_scale_data = FALSE,
                                              log_function = log_system$add_log)
      }
      # Store integrated object
      multiple_datasets_object(integrated_object)
      init_clustering_plot_merge(integrated_object)
      # Method and batch variable recorded for the parameter export and the PDF report
      integration_method_used("Harmony")
      harmony_vars_used(harmony_vars_to_use)
      logMemoryUsage(log_system$add_log, integrated_object, "Harmony integrated object")
      # Enable UI controls
      shinyjs::enable("add_field")
      shinyjs::enable("add_metadata")
      removeModal()
      # Gene picker of the visualization tab (otherwise empty after Harmony)
      updateUIElements()
      cleanupIntegrationMemory()
      showNotification(
        paste("Harmony integration completed successfully!",
              ncol(integrated_object), "cells integrated"),
        type = "message",
        duration = 10
      )
    }, error = function(e) {
      removeModal()
      showNotification(
        paste("Harmony integration failed:", e$message, 
              "\nTry Simple Merge or check your data."), 
        type = "error",
        duration = 15
      )
    })
  })

  observe(run_safely("Interface update", persistent = FALSE, {
    req(input$num_datasets)
    
    first_loaded <- NULL
    for (i in 1:input$num_datasets) {
      seurat_obj <- seurat_objects[[paste0("seurat_object", i)]]
      if (!is.null(seurat_obj)) {
        first_loaded <- seurat_obj
        break
      }
    }
    
    if (!is.null(first_loaded)) {
      meta_cols <- colnames(first_loaded@meta.data)
      
      is_secondary_integration <- "integration_group" %in% meta_cols
      
      choices_list <- list()
      
      if (is_secondary_integration) {
        choices_list[["Integration groups (current level)"]] <- "integration_group"
        choices_list[["Original datasets (all samples)"]] <- "dataset"
        default_selection <- "integration_group"
        
        message("✓ Secondary integration detected - offering both integration levels")
      } else {
        choices_list[["Datasets"]] <- "dataset"
        default_selection <- "dataset"
        
        message("✓ First integration - using dataset column")
      }
      
      useful_vars <- meta_cols[!meta_cols %in% c(
        "nCount_RNA", "nFeature_RNA", "percent.mt", 
        "orig.ident", "seurat_clusters", "dataset", "integration_group",
        "dataset_origin", "integration_level"
      )]
      
      for (var in useful_vars) {
        n_levels <- length(unique(first_loaded@meta.data[[var]]))
        if (n_levels >= 2 && n_levels <= 20) {
          choices_list[[paste0(var, " (", n_levels, " levels)")]] <- var
        }
      }
      
      updateSelectInput(
        session, 
        "harmony_vars_integration",
        choices = choices_list,
        selected = default_selection
      )
    }
  }))

  # Metadata management using modular functions
  observeEvent(input$add_field, run_safely("add field", {
    rv_metadata$num_fields <- rv_metadata$num_fields + 1
  }))

  reactive_metadata_fields <- reactive({
    req(multiple_datasets_object())
    all_metadata_fields <- colnames(multiple_datasets_object()@meta.data)
    
    # Fields to exclude
    exclude_fields <- c("orig.ident", "nCount_RNA", "nCount_ATAC", "nFeature_RNA",
                        "nFeature_ATAC", "percent.mt")
    
    # Pattern-based exclusions
    exclude_pattern <- "^snn_res|^pANN|^PC_|^RNA_snn|^ATAC_snn|^integrated_snn"
    exclude_fields <- c(exclude_fields, all_metadata_fields[grepl(exclude_pattern, all_metadata_fields)])
    
    include_fields <- setdiff(all_metadata_fields, exclude_fields)
    return(include_fields)
  })

  # Process metadata using modular function
  observeEvent(input$add_metadata, {
    tryCatch({
      req(multiple_datasets_object())
      req(input$metadata_split_column)  # NEW: require split column
      
      # Validate that at least one field has a name
      has_valid_field <- FALSE
      for (j in 1:rv_metadata$num_fields) {
        field_name <- input[[paste0("metadata_name_", j)]]
        if (!is.null(field_name) && nchar(trimws(field_name)) > 0) {
          has_valid_field <- TRUE
          break
        }
      }
      
      if (!has_valid_field) {
        showNotification(
          "Please enter at least one metadata field name before applying.",
          type = "warning",
          duration = 5
        )
        return()
      }
      
      # Use modular metadata processing function with split column
      updated_seurat <- processMetadataFromUI(
        seurat_object = multiple_datasets_object(),
        input = input,
        num_fields = rv_metadata$num_fields,
        split_column = input$metadata_split_column  # NEW: pass split column
      )
      
      # Update the reactive object
      multiple_datasets_object(updated_seurat)
      
      # Show success notification
      split_values <- unique(updated_seurat@meta.data[[input$metadata_split_column]])
      fields_added <- 0
      for (j in 1:rv_metadata$num_fields) {
        field_name <- input[[paste0("metadata_name_", j)]]
        if (!is.null(field_name) && nchar(trimws(field_name)) > 0) {
          fields_added <- fields_added + 1
        }
      }
      
      showNotification(
        HTML(paste0(
          "<strong>✓ Metadata applied successfully!</strong><br>",
          "Added ", fields_added, " metadata column(s) split by '", input$metadata_split_column, "' across ", 
          length(split_values), " groups."
        )),
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      showNotification(
        paste("Error adding metadata:", e$message), 
        type = "error",
        duration = 10
      )
    })
  })

  # Reactive gene list for multiple datasets
  reactive_gene_list_merge <- reactive({
    req(multiple_datasets_object())
    
    tryCatch({
      obj <- multiple_datasets_object()
      
      # Detect available assays
      available_assays <- names(obj@assays)
      
      # Try to find the right assay
      if ("RNA" %in% available_assays) {
        assay_to_use <- "RNA"
      } else if ("Spatial" %in% available_assays) {
        assay_to_use <- "Spatial"
      } else {
        # Use the default assay if neither RNA nor Spatial
        assay_to_use <- DefaultAssay(obj)
      }
      
      message(paste("Using assay:", assay_to_use, "for gene list"))
      
      # Get genes from the appropriate assay
      # Features of the assay: LayerData(layer = "counts") reads only the first layer (counts.1) of a merged object
      unique_genes <- rownames(obj[[assay_to_use]])
      return(c("", unique_genes))
      
    }, error = function(e) {
      showNotification(
        paste("Error loading gene list:", e$message, 
              "\nPlease ensure your Seurat object has either 'RNA' or 'Spatial' assay."),
        type = "error",
        duration = 10
      )
      return(c(""))  # Return empty list on error
    })
  })

  # ---- UI and text outputs ----

  # Display loaded datasets in metadata management tab (add before output$metadata_inputs)
  output$loaded_datasets_display_metadata <- renderUI({
    req(multiple_datasets_object())
    
    datasets <- unique(multiple_datasets_object()@meta.data$dataset)
    n_cells_per_dataset <- sapply(datasets, function(ds) {
      sum(multiple_datasets_object()@meta.data$dataset == ds)
    })
    
    tags$div(
      style = "display: grid; grid-template-columns: repeat(auto-fill, minmax(220px, 1fr)); gap: 12px;",
      lapply(seq_along(datasets), function(i) {
        tags$div(
          style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
             padding: 15px; border-radius: 8px; color: white;
             box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
          tags$div(
            style = "font-weight: bold; font-size: 16px; margin-bottom: 5px;",
            icon("layer-group"), " ", datasets[i]
          ),
          tags$div(
            style = "font-size: 13px; opacity: 0.9;",
            format(n_cells_per_dataset[i], big.mark = ","), " cells"
          )
        )
      })
    )
  })

  # Display existing metadata columns
  output$existing_metadata_display <- renderUI({
    req(multiple_datasets_object())
    
    # Get all metadata columns
    all_cols <- colnames(multiple_datasets_object()@meta.data)
    
    # Exclude standard columns
    exclude_cols <- c("orig.ident", "nCount_RNA", "nCount_ATAC", "nFeature_RNA",
                      "nFeature_ATAC", "percent.mt", "dataset")
    
    # Pattern-based exclusions
    exclude_pattern <- "^snn_res|^pANN|^PC_|^RNA_snn|^ATAC_snn|^integrated_snn|^seurat_clusters"
    
    # Filter columns
    metadata_cols <- all_cols[!all_cols %in% exclude_cols]
    metadata_cols <- metadata_cols[!grepl(exclude_pattern, metadata_cols)]
    
    if (length(metadata_cols) == 0) {
      return(
        tags$div(
          style = "text-align: center; padding: 20px; color: #6c757d;",
          icon("inbox", style = "font-size: 32px; margin-bottom: 10px;"),
          tags$p("No custom metadata columns yet.", style = "margin: 0; font-size: 14px;")
        )
      )
    }
    
    # Display metadata columns with their unique values
    tags$div(
      style = "display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); gap: 12px;",
      lapply(metadata_cols, function(col) {
        unique_vals <- unique(multiple_datasets_object()@meta.data[[col]])
        unique_vals <- unique_vals[!is.na(unique_vals)]
        
        tags$div(
          style = "background-color: #f8f9fa; padding: 12px; border-radius: 6px; 
             border-left: 4px solid #f093fb;",
          tags$div(
            style = "font-weight: bold; color: #495057; margin-bottom: 8px; font-size: 14px;",
            icon("tag"), " ", col
          ),
          tags$div(
            style = "font-size: 12px; color: #6c757d;",
            if (length(unique_vals) > 0) {
              paste(length(unique_vals), "unique values:", paste(head(unique_vals, 3), collapse = ", "),
                    if (length(unique_vals) > 3) "..." else "")
            } else {
              "No values"
            }
          )
        )
      })
    )
  })

  # Metadata split column selector
  output$metadata_split_selector <- renderUI({
    req(multiple_datasets_object())
    
    # Get available columns from metadata (excluding technical columns)
    all_cols <- colnames(multiple_datasets_object()@meta.data)
    
    # Exclude technical/standard columns that shouldn't be used for splitting
    exclude_cols <- c("orig.ident", "nCount_RNA", "nCount_ATAC", "nFeature_RNA",
                      "nFeature_ATAC", "percent.mt")
    exclude_pattern <- "^snn_res|^pANN|^PC_|^RNA_snn|^ATAC_snn|^integrated_snn"
    
    available_cols <- all_cols[!all_cols %in% exclude_cols]
    available_cols <- available_cols[!grepl(exclude_pattern, available_cols)]
    
    # Default to "dataset" if it exists
    default_col <- if ("dataset" %in% available_cols) "dataset" else available_cols[1]
    
    tags$div(
      style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
             padding: 20px; border-radius: 8px; margin-bottom: 20px; 
             box-shadow: 0 4px 6px rgba(0,0,0,0.1);",
      tags$div(
        style = "display: flex; align-items: center; margin-bottom: 15px;",
        icon("layer-group", style = "font-size: 24px; color: white; margin-right: 12px;"),
        tags$h5(
          "Select Split Column",
          style = "color: white; margin: 0; font-weight: bold;"
        )
      ),
      tags$p(
        "Choose which column to use for grouping your datasets. 
       Each unique value in this column will get its own metadata fields.",
        style = "color: rgba(255,255,255,0.9); margin-bottom: 15px; font-size: 14px;"
      ),
      selectInput(
        "metadata_split_column",
        NULL,
        choices = available_cols,
        selected = default_col,
        width = "100%"
      ),
      tags$div(
        style = "margin-top: 10px; padding: 10px; background-color: rgba(255,255,255,0.1); 
               border-radius: 5px; border-left: 3px solid white;",
        icon("info-circle", style = "color: white;"),
        tags$span(
          paste0(" Currently splitting by: ", default_col),
          style = "color: white; font-size: 13px; margin-left: 5px;"
        )
      )
    )
  })

  output$metadata_inputs <- renderUI({
    req(multiple_datasets_object())
    req(input$metadata_split_column)  # NEW: require split column selection
    
    # Get the column to split by
    split_column <- input$metadata_split_column
    
    # Get unique values from the selected column
    split_values <- unique(multiple_datasets_object()@meta.data[[split_column]])
    split_values <- sort(split_values[!is.na(split_values)])
    
    # If no fields added yet
    if (rv_metadata$num_fields == 0) {
      return(
        tags$div(
          style = "text-align: center; padding: 40px; background-color: #f8f9fa; 
             border-radius: 8px; border: 2px dashed #dee2e6;",
          icon("plus-circle", style = "font-size: 48px; color: #6c757d; margin-bottom: 15px;"),
          tags$h5("No metadata fields added yet", style = "color: #495057; margin-bottom: 10px;"),
          tags$p("Click 'Add New Metadata Field' button below to create your first metadata column.", 
                 style = "color: #6c757d; margin: 0;")
        )
      )
    }
    
    # Generate metadata input fields
    lapply(1:rv_metadata$num_fields, function(j) {
      tags$div(
        style = "background: linear-gradient(to right, #f8f9fa 0%, #ffffff 100%); 
           padding: 20px; border-radius: 8px; margin-bottom: 20px; 
           border: 2px solid #dee2e6;",
        
        # Field header with number
        tags$div(
          style = "display: flex; justify-content: space-between; align-items: center; 
             margin-bottom: 15px; padding-bottom: 10px; border-bottom: 2px solid #dee2e6;",
          tags$h5(
            icon("edit"), " Metadata Field ", j,
            style = "color: #495057; margin: 0; font-weight: bold;"
          ),
          tags$span(
            icon("info-circle"), " This will create a new column in your dataset",
            style = "color: #6c757d; font-size: 12px;"
          )
        ),
        
        # Column name input
        tags$div(
          style = "margin-bottom: 20px;",
          tags$label(
            "Column Name:",
            style = "font-weight: 600; color: #495057; margin-bottom: 8px; 
               display: block; font-size: 14px;"
          ),
          textInput(
            paste0("metadata_name_", j), 
            NULL,
            value = "",
            placeholder = "e.g., condition, timepoint, treatment, replicate, batch..."
          ),
          tags$small(
            icon("lightbulb"), " Use descriptive names like 'condition', 'timepoint', 'replicate'",
            style = "color: #6c757d; font-size: 12px;"
          )
        ),
        
        # Values for each split value
        tags$div(
          tags$label(
            paste0("Assign value for each ", split_column, ":"),
            style = "font-weight: 600; color: #495057; margin-bottom: 12px; 
               display: block; font-size: 14px;"
          ),
          tags$div(
            style = "display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); 
               gap: 12px;",
            lapply(split_values, function(split_value) {
              # Get number of cells for this value
              n_cells <- sum(multiple_datasets_object()@meta.data[[split_column]] == split_value)
              
              tags$div(
                style = "background-color: white; padding: 15px; border-radius: 6px; 
                   border-left: 4px solid #f093fb; box-shadow: 0 2px 4px rgba(0,0,0,0.05);",
                tags$div(
                  style = "display: flex; justify-content: space-between; align-items: center; 
                     margin-bottom: 8px;",
                  tags$div(
                    style = "font-weight: bold; color: #495057; font-size: 13px;",
                    icon("tag"), " ", split_value
                  ),
                  tags$span(
                    format(n_cells, big.mark = ","), " cells",
                    style = "font-size: 11px; color: #6c757d;"
                  )
                ),
                textInput(
                  paste0("metadata_value_", split_value, "_", j), 
                  NULL,
                  value = "",
                  placeholder = paste0("Value for ", split_value, "...")
                )
              )
            })
          )
        )
      )
    })
  })

  outputOptions(output, 'datasets_loaded', suspendWhenHidden = FALSE)

  ############################## Scaling and PCA reduction ##############################

  # ---- Reactive state ----

  # Observer forfind Clusters and display UMAP
  clustering_plot_merge <- reactiveVal()

  # Metadata column the stored raw plot is grouped by (needed to place cluster names on the UMAP)
  clustering_group_col_merge <- reactiveVal(NULL)

  # ---- Helpers ----

  # Reduction of the Cluster Cells plot: shared UMAP selector when it exists in this object, otherwise the default UMAP
  selected_reduction_merge <- function(seurat_obj) {
    choice <- isolate(input$cluster_cells_reduction_merge)
    if (!is.null(choice) && nzchar(choice) && choice %in% names(seurat_obj@reductions)) return(choice)
    return(get_umap_reduction(seurat_obj))
  }

  # Raw plot of Cluster Cells in the selected reduction; group_col NULL = active identities when clusters exist, samples otherwise
  init_clustering_plot_merge <- function(seurat_obj, group_col = NULL) {
    reduction <- selected_reduction_merge(seurat_obj)
    if (is.null(reduction)) {
      # No UMAP yet: clear the plot of the previous object
      clustering_plot_merge(NULL)
      clustering_group_col_merge(NULL)
      return(invisible(NULL))
    }
    meta_cols <- colnames(seurat_obj@meta.data)
    if (is.null(group_col) || !group_col %in% c(meta_cols, "ident")) {
      # Active identities follow FindClusters, renaming, undo and applied annotation labels
      group_col <- if (any(c("ClusterIdents", "seurat_clusters") %in% meta_cols)) "ident"
      else if ("integration_group" %in% meta_cols) "integration_group"
      else if ("dataset" %in% meta_cols) "dataset"
      else "ident"
    }
    tryCatch({
      plot <- DimPlot(seurat_obj, reduction = reduction, group.by = group_col, label = FALSE) + ggtitle(NULL)
      # Stored cluster colours reused only when they cover every group of the plotted column
      cluster_colors <- get_cluster_colors_merge(seurat_obj)
      group_values <- unique(as.character(plot$data[[group_col]]))
      if (!is.null(cluster_colors) && all(group_values %in% names(cluster_colors))) plot <- plot + scale_color_manual(values = cluster_colors)
      clustering_group_col_merge(group_col)
      clustering_plot_merge(plot)
      message("Cluster Cells plot built — reduction: ", reduction, " — column: ", group_col)
    }, error = function(e) {
      message("Cluster Cells plot error: ", conditionMessage(e)[1])
      showNotification(paste("The plot could not be displayed in Cluster Cells:", conditionMessage(e)[1]), type = "warning", duration = 8)
    })
    return(invisible(NULL))
  }

  # ---- Observers and computations ----

  # Scaling, PCA and Elbowplot
  observeEvent(input$runScalePCA, {
    req(multiple_datasets_object())
    tryCatch({
      showModal(modalDialog(
        title = "Please Wait",
        "Processing...",
        easyClose = FALSE,
        footer = NULL
      ))
      seurat_object_temp <- multiple_datasets_object()
      current_assay <- DefaultAssay(seurat_object_temp)
      if ("harmony" %in% names(seurat_object_temp@reductions)) {
        message("Harmony reduction detected - skipping scaling and PCA")
        output$elbow_plot2 <- renderPlot({
          ElbowPlot(seurat_object_temp, reduction = "harmony", ndims = 50) +
            ggtitle("Harmony Elbow Plot") +
            theme_minimal()
        })
        showNotification(
          "Harmony integration already completed. Showing Harmony elbow plot. 
        Scaling and PCA not needed - proceed to Find Neighbors.",
          type = "message",
          duration = 8
        )
      } else if (current_assay == "SCT" || "SCT" %in% names(seurat_object_temp@assays)) {
        message("SCTransform detected - data already scaled, running PCA on SCT assay")
        showNotification("SCTransform detected - running PCA...", type = "message")
        DefaultAssay(seurat_object_temp) <- "SCT"
        # Seed-42 wrapper (dimensional_reduction_functions.R)
        seurat_object_temp <- runPCA_reproducible(seurat_object_temp, npcs = 50, seed = 42, features = VariableFeatures(seurat_object_temp), verbose = FALSE)
        multiple_datasets_object(seurat_object_temp)
        req(seurat_object_temp[["pca"]])
        output$elbow_plot2 <- renderPlot({
          ElbowPlot(seurat_object_temp, reduction = "pca", ndims = 50) +
            ggtitle("PCA Elbow Plot (SCT assay)") +
            theme_minimal()
        })
        showNotification("PCA completed on SCTransform data.", type = "message")
      } else {
        message("Running scaling and PCA for integrated data (RNA assay)")
        showNotification("Scaling and PCA have begun...", type = "message")
        seurat_object_temp <- FindVariableFeatures(seurat_object_temp, selection.method = "vst", nfeatures = 3000)
        # Only the variable features enter the PCA: scaling every gene of the merged object would build a dense genes x cells matrix
        seurat_object_temp <- ScaleData(seurat_object_temp, features = VariableFeatures(seurat_object_temp))
        seurat_object_temp <- runPCA_reproducible(seurat_object_temp, npcs = 50, seed = 42, features = VariableFeatures(seurat_object_temp), verbose = FALSE)
        multiple_datasets_object(seurat_object_temp)
        req(seurat_object_temp[["pca"]])
        output$elbow_plot2 <- renderPlot({
          ElbowPlot(seurat_object_temp, reduction = "pca", ndims = 50) +
            ggtitle("PCA Elbow Plot") +
            theme_minimal()
        })
        showNotification("Scaling and PCA are complete.", type = "message")
      }
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Scaling and PCA errors:", e$message), type = "error")
    })
  })

  # Any change of the object (renaming, undo, applied labels, colours) refreshes the Cluster Cells plot;
  # sample colouring (after Find Neighbors) is kept, any cluster colouring switches to the active identities
  observeEvent(multiple_datasets_object(), run_safely("Cluster Cells refresh", persistent = FALSE, {
    obj <- multiple_datasets_object()
    req(obj, !is.null(isolate(clustering_plot_merge())))
    current_col <- isolate(clustering_group_col_merge())
    init_clustering_plot_merge(obj, group_col = if (isTRUE(current_col %in% c("integration_group", "dataset"))) current_col else "ident")
  }), ignoreInit = TRUE)

  # Changing the reduction rebuilds the plot with the same colouring column and propagates the choice to the other UMAP selectors
  observeEvent(input$cluster_cells_reduction_merge, run_safely("cluster cells reduction merge", {
    obj <- multiple_datasets_object()
    req(obj, nzchar(input$cluster_cells_reduction_merge))
    init_clustering_plot_merge(obj, group_col = clustering_group_col_merge())
    if (!identical(input$umap_reduction_merge, input$cluster_cells_reduction_merge))
      updateSelectInput(session, "umap_reduction_merge", selected = input$cluster_cells_reduction_merge)
  }), ignoreInit = TRUE)

  observeEvent(input$runFindNeighbors, {
    tryCatch({
      showModal(modalDialog(
        title = "Please Wait",
        "Finding neighbors and running UMAP...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      req(multiple_datasets_object())
      seurat_object_temp <- multiple_datasets_object()
      
      # NOUVEAU: Join layers if needed (Seurat v5) BEFORE any analysis
      if (packageVersion("Seurat") >= "5.0.0") {
        message("Checking for Seurat v5 layers...")
        
        tryCatch({
          for (assay_name in names(seurat_object_temp@assays)) {
            assay_obj <- seurat_object_temp[[assay_name]]
            
            # Check if layers exist and are multiple
            if ("layers" %in% slotNames(assay_obj)) {
              layer_names <- names(assay_obj@layers)
              if (length(layer_names) > 1) {
                message(paste("Joining", length(layer_names), "layers in", assay_name))
                seurat_object_temp[[assay_name]] <- JoinLayers(assay_obj)
              }
            }
          }
        }, error = function(e) {
          message(paste("Note: Could not join layers:", e$message))
        })
      }
      
      # Determine correct reduction
      available_reductions <- names(seurat_object_temp@reductions)
      message(paste("Available reductions:", paste(available_reductions, collapse = ", ")))
      
      reduction_to_use <- if ("integrated.cca" %in% available_reductions) {
        "integrated.cca"
      } else if ("integrated.rpca" %in% available_reductions) {
        "integrated.rpca" 
      } else if ("integrated.mnn" %in% available_reductions) {
        "integrated.mnn"
      } else if ("harmony" %in% available_reductions) {
        "harmony"
      } else if ("pca" %in% available_reductions) {
        "pca"
      } else {
        stop("No suitable reduction found")
      }
      
      umap_dims <- input$dimension_2
      message(paste("Using reduction:", reduction_to_use, "with", umap_dims, "dimensions"))
      
      # FindNeighbors
      message("Running FindNeighbors...")
      seurat_object_temp <- findNeighbors_reproducible(
        object = seurat_object_temp,
        reduction = reduction_to_use,
        dims = 1:umap_dims,
        seed = 42,
        verbose = FALSE
      )
      
      # UMAP
      message("Running UMAP...")
      seurat_object_temp <- runUMAP_reproducible(
        object = seurat_object_temp,
        reduction = reduction_to_use,
        dims = 1:umap_dims,
        seed = 42,
        verbose = FALSE
      )
      
      multiple_datasets_object(seurat_object_temp)
      
      # Raw plot coloured by sample in the selected reduction; display options are applied at render time
      init_clustering_plot_merge(seurat_object_temp, group_col = if ("integration_group" %in% colnames(seurat_object_temp@meta.data)) "integration_group" else "dataset")
      
      showNotification("Finding neighbors and UMAP completed.", type = "message")
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", e$message), type = "error")
    })
  })

  observeEvent(input$runFindClusters, {
    tryCatch({
      showModal(modalDialog(
        title = "Please Wait",
        "Finding clusters...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      req(multiple_datasets_object())
      seurat_integrated_temp <- multiple_datasets_object()
      
      if (is.null(seurat_integrated_temp@graphs) || length(seurat_integrated_temp@graphs) == 0) {
        showNotification("No neighbor graph found! Run 'Find Neighbors' first.", type = "error")
        removeModal()
        return()
      }
      
      available_graphs <- names(seurat_integrated_temp@graphs)
      message(paste("Available graphs:", paste(available_graphs, collapse = ", ")))
      # SNN graph of the last FindNeighbors call (FindNeighbors names its graphs <assay>_nn and <assay>_snn);
      # the kNN graph (<assay>_nn) is never used for clustering
      graph_name <- seurat_integrated_temp@misc$neighbors_params$snn_graph
      if (is.null(graph_name) || is.na(graph_name) || !graph_name %in% available_graphs) {
        snn_graphs <- grep("snn$", available_graphs, value = TRUE)
        if (length(snn_graphs) == 0) {
          showNotification("No shared nearest-neighbour (SNN) graph found! Run 'Find Neighbors' first.", type = "error")
          removeModal()
          return()
        }
        graph_name <- tail(snn_graphs, 1)
      }
      message(paste("Using graph:", graph_name))
      
      metadata_backup <- seurat_integrated_temp@meta.data
      
      seurat_integrated_temp <- findClusters_reproducible(
        object     = seurat_integrated_temp,
        resolution = input$resolution_step2,
        algorithm  = as.integer(input$algorithm_select_merge),
        seed       = 42,
        verbose    = TRUE,
        graph.name = graph_name
      )
      
      if (!"seurat_clusters" %in% colnames(seurat_integrated_temp@meta.data)) {
        warning("seurat_clusters not found, trying to restore...")
        seurat_integrated_temp@meta.data <- metadata_backup
        stop("Clustering failed to create seurat_clusters column")
      }
      
      n_clusters <- length(unique(seurat_integrated_temp$seurat_clusters))
      message(paste("Created", n_clusters, "clusters"))
      
      multiple_datasets_object(seurat_integrated_temp)
      
      # Raw plot coloured by the new clusters (active identities) in the selected reduction; display options are applied at render time
      init_clustering_plot_merge(seurat_integrated_temp, group_col = "ident")
      
      showNotification(paste0("Clustering completed! Found ", n_clusters, " clusters."), type = "message")
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Error: ", e$message), type = "error")
      message(paste("Full error:", e))
    })
  })

  # Display options applied on top of the stored raw plot, shared by the screen and the download
  clustering_plot_merge_styled <- reactive({
    req(clustering_plot_merge())
    tryCatch(
      style_cluster_umap(clustering_plot_merge(), label_mode = input$cluster_label_mode_merge %||% "both",
                         label_size = input$cluster_label_size_merge %||% 5, bold_labels = isTRUE(input$cluster_bold_labels_merge),
                         dark_mode = isTRUE(input$dark_mode_umap_cluster_merge), remove_axes = isTRUE(input$remove_axes_umap_merge),
                         id = clustering_group_col_merge() %||% "seurat_clusters"),
      error = function(e) {
        message("Clustering plot styling error (merge): ", conditionMessage(e)[1])
        showNotification(paste("Display options could not be applied:", conditionMessage(e)[1]), type = "error", duration = 8)
        clustering_plot_merge()
      })
  })

  # ---- Plots ----

  output$UMAPPlot_cluster_merge <- renderPlot({
    clustering_plot_merge_styled()
  })

  # ---- Downloads ----

  output$downloadUMAP_merge <- createDownloadHandler(
    reactive_data = clustering_plot_merge_styled,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedUMAP") 
    }),
    data_name = "UMAP_plot",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$umap_merge_format }),
      width = reactive({ ifelse(input$umap_merge_format == "pdf", 11, 10) }),
      height = reactive({ ifelse(input$umap_merge_format == "pdf", 8, 6) }),
      dpi = reactive({ input$dpi_umap_merge })
    )
  )

  ############################## Visualize genes expressions ##############################
  
  ############################## OBSERVERS FOR METADATA SELECTION ##############################
  ############################## Gene Selection Synchronization ##############################

  # ---- Observers and computations ----

  # Observer: When genes are selected in pickerInput, update all text inputs
  observeEvent(input$geneInput_merge, run_safely("geneInput merge", {
    selected_genes <- input$geneInput_merge
    
    if (!is.null(selected_genes) && length(selected_genes) > 0) {
      genes_text <- paste(selected_genes, collapse = ", ")
      
      # Update all gene textInput fields
      updateTextInput(session, "gene_list_vln_merge", value = genes_text)
      updateTextInput(session, "gene_list_feature_merge", value = genes_text)
      updateTextInput(session, "gene_list_dot_merge", value = genes_text)
      updateTextInput(session, "gene_list_ridge_merge", value = genes_text)
      updateTextInput(session, "gene_list_genes_expression_merge", value = genes_text)
      
      message(paste("Updated all gene inputs with:", genes_text))
    }
  }))

  # OBSERVER 1: Update group_by_select choices (main axis of plots - X axis)
  # OBSERVER 1: Update group_by_select choices (main axis of plots - X axis)
  observeEvent(multiple_datasets_object(), run_safely("Interface update", persistent = FALSE, {
    seurat_object <- multiple_datasets_object()
    req(seurat_object)
    # Exclude technical columns; cluster_name_only is kept because it is the only column holding cluster names when the source objects had no ClusterIdents
    excluded_patterns <- c(
      "percent.mt", "nCount_ATAC", "nFeature_ATAC", "nFeature_RNA", "nCount_RNA",
      "nFeature_Spatial", "nCount_Spatial", "^RNA_snn_", "^RNA_nn_",
      "^Spatial_snn_", "^Spatial_nn_", "original_clusters", "^pANN", "^DF",
      "dataset_origin", "original_seurat_clusters", "original_ClusterIdents",
      "source_format", "^integrated", "^integrated_snn",
      "^merge_method$", "^integration_method$", "^original_clusterIdents$"
    )
    pattern <- paste(excluded_patterns, collapse = "|")
    metadata_fields <- colnames(seurat_object@meta.data)
    display_fields <- metadata_fields[!grepl(pattern, metadata_fields)]
    if (length(display_fields) == 0) {
      message("No valid metadata fields found for group_by")
      return()
    }
    # Default priority: current selection > ClusterIdents > cluster_name_only > seurat_clusters
    current_selection <- isolate(input$group_by_select)
    default_choice <- if (!is.null(current_selection) && current_selection != "" && current_selection %in% display_fields) {
      current_selection
    } else if ("ClusterIdents" %in% display_fields) {
      "ClusterIdents"
    } else if ("cluster_name_only" %in% display_fields) {
      "cluster_name_only"
    } else if ("seurat_clusters" %in% display_fields) {
      "seurat_clusters"
    } else {
      ""
    }
    choices_with_empty <- c("Select grouping..." = "", display_fields)
    updateSelectizeInput(session, "group_by_select", choices = choices_with_empty, selected = default_choice, server = TRUE)
    message(paste("Updated group_by_select. Available:", paste(display_fields, collapse = ", "), "| default:", default_choice))
  }), priority = 10, ignoreNULL = FALSE)

  # OBSERVER 2: Update metadata_to_compare choices (for split.by visual separation)
  observeEvent(list(multiple_datasets_object(), input$group_by_select), run_safely("Action", {
    seurat_object <- multiple_datasets_object()
    req(seurat_object)
    
    # Same exclusions as group_by
    excluded_patterns <- c(
      "percent.mt", "nCount_ATAC", "nFeature_ATAC", "nFeature_RNA", "nCount_RNA", 
      "nFeature_Spatial", "nCount_Spatial", "^RNA_snn_", "^RNA_nn_", 
      "^Spatial_snn_", "^Spatial_nn_", "original_clusters", "^pANN", "^DF", 
      "dataset_origin", "original_seurat_clusters", "original_ClusterIdents",
      "^SCT_snn", "^integrated_snn"
    )
    
    metadata_fields <- colnames(seurat_object@meta.data)
    pattern <- paste(excluded_patterns, collapse = "|")
    display_fields <- metadata_fields[!grepl(pattern, metadata_fields)]
    
    # Remove the current group_by column from split choices
    group_by_val <- input$group_by_select
    if (!is.null(group_by_val) && group_by_val != "") {
      split_fields <- display_fields[display_fields != group_by_val]
    } else {
      split_fields <- display_fields
    }
    
    # Use isolate to get current selection without triggering
    current_selection <- isolate(input$metadata_to_compare)
    
    # NEW LOGIC: Default to "None" unless user has already made a selection
    if (!is.null(current_selection) && current_selection != "" && current_selection %in% split_fields) {
      # Keep user's current selection if it's still valid
      default_split <- current_selection
    } else {
      # Always default to "None" for new selections
      default_split <- ""
    }
    
    updateSelectizeInput(session, "metadata_to_compare",
                         choices = c("None (no split)" = "", split_fields),
                         selected = default_split,
                         server = TRUE)
    
    message(paste("Updated metadata_to_compare (split.by options). Available:", paste(head(split_fields, 3), collapse=", ")))
    
  }), priority = 9, ignoreNULL = FALSE)

  # OBSERVER 3: Update cluster_order choices based on group_by_select values
  observeEvent(input$group_by_select, run_safely("group by select", {
    req(multiple_datasets_object())
    
    # If no group_by selected, clear the cluster_order selectors
    if (is.null(input$group_by_select) || input$group_by_select == "") {
      updateSelectizeInput(session, "cluster_order_dotplot_merge",
                           choices = character(0),
                           selected = character(0),
                           server = TRUE)
      updateSelectizeInput(session, "cluster_order_vln_merge",
                           choices = character(0),
                           selected = character(0),
                           server = TRUE)
      message("No group_by selected - cluster_order cleared")
      return()
    }
    
    seurat_object <- multiple_datasets_object()
    group_by <- input$group_by_select
    
    # Get unique values from the group_by column
    clusters <- tryCatch({
      if (group_by == "seurat_clusters") {
        # Use Idents for seurat_clusters
        cluster_idents <- Idents(object = seurat_object)
        as.character(unique(cluster_idents))
      } else if (group_by %in% colnames(seurat_object@meta.data)) {
        # Use metadata column
        as.character(unique(seurat_object@meta.data[[group_by]]))
      } else {
        character(0)
      }
    }, error = function(e) {
      message(paste("Error getting cluster values:", e$message))
      character(0)
    })
    
    # Clean NAs and empty strings
    clusters <- clusters[!is.na(clusters) & clusters != ""]
    
    if (length(clusters) > 0) {
      # Update all cluster order selectors
      updateSelectizeInput(session, "cluster_order_dotplot_merge",
                           choices = clusters,
                           selected = clusters,
                           server = TRUE)
      updateSelectizeInput(session, "cluster_order_vln_merge",
                           choices = clusters,
                           selected = clusters,
                           server = TRUE)
      
      message(paste("Updated cluster_order from", group_by, ":", paste(head(clusters, 3), collapse=", ")))
    } else {
      message(paste("No valid clusters found for", group_by))
    }
    
  }), priority = 8, ignoreNULL = FALSE)

  # ---- UI and text outputs ----

  # OBSERVER 4: Dynamic UI to filter split values
  output$split_values_filter_ui_gene_viz <- renderUI({
    req(multiple_datasets_object(), input$metadata_to_compare)
    
    # Only show if a split column is selected
    if (is.null(input$metadata_to_compare) || input$metadata_to_compare == "") {
      return(NULL)
    }
    
    split_col <- input$metadata_to_compare
    seurat_obj <- multiple_datasets_object()
    
    # Check column exists
    if (!split_col %in% colnames(seurat_obj@meta.data)) {
      return(NULL)
    }
    
    # Get unique values from the split column
    available_values <- unique(seurat_obj@meta.data[[split_col]])
    available_values <- sort(available_values[!is.na(available_values)])
    
    if (length(available_values) == 0) {
      return(NULL)
    }
    
    selectizeInput("split_values_filter_gene_viz",
                   paste("Filter", split_col, "values (optional):"),
                   choices = available_values,
                   selected = NULL,
                   multiple = TRUE,
                   options = list(
                     plugins = list('remove_button'),
                     placeholder = 'Select values (empty = show all)'
                   ))
  })

  ############################## FEATURE PLOT ##############################

  # ---- Observers and computations ----

  # Generate feature plot
  observeEvent(input$show_feature_plot_merge, {
    tryCatch({
      seurat_object <- multiple_datasets_object()
      if (is.null(seurat_object)) { showNotification("Load or integrate datasets first.", type = "warning"); return() }
      DefaultAssay(seurat_object) <- input$viz_assay_merge
      # Empty field and unknown genes are reported; genes are checked in the selected assay
      check <- validate_gene_input(input$gene_list_feature_merge, seurat_object, input$viz_assay_merge)
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      present_genes <- check$genes
      
      message(paste("Plotting genes:", paste(present_genes, collapse = ", ")))
      
      # Apply filters if specified
      if(!is.null(input$metadata_to_compare) && input$metadata_to_compare != "") {
        split_col <- input$metadata_to_compare
        
        if(split_col %in% colnames(seurat_object@meta.data)) {
          # Get selected values
          selected_vals <- input$split_values_filter_gene_viz          
          if(!is.null(selected_vals) && length(selected_vals) > 0) {
            cells_to_keep <- seurat_object@meta.data[[split_col]] %in% selected_vals
            if(sum(cells_to_keep) == 0) {
              showNotification("No cells match the selected filter values", type = "warning")
              return()
            }
            seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
            message(paste("Filtered to", ncol(seurat_object), "cells"))
          }
        } else {
          showNotification("Selected metadata column not found", type = "warning")
          return()
        }
      } else {
        # Check for cells_to_keep from group_by filtering
        if(!is.null(input$group_by_select) && input$group_by_select != "") {
          group_by_col <- input$group_by_select
          selected_groups <- input$cluster_order_feature_merge
          
          if(!is.null(selected_groups) && length(selected_groups) > 0 && 
             group_by_col %in% colnames(seurat_object@meta.data)) {
            cells_to_keep <- seurat_object@meta.data[[group_by_col]] %in% selected_groups
            if(sum(cells_to_keep) > 0) {
              seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
              message(paste("Filtered to", ncol(seurat_object), "cells by group_by"))
            }
          }
        }
        
        if(ncol(seurat_object) == 0) {
          showNotification("No cells match the selected filter values", type = "warning")
          return()
        }
      }
      
      # Cutoffs handling
      min_cut <- if (!is.null(input$min_cutoff_feature_merge) && 
                     !is.na(input$min_cutoff_feature_merge) && 
                     input$min_cutoff_feature_merge != "") {
        if (grepl("^q", input$min_cutoff_feature_merge)) {
          input$min_cutoff_feature_merge
        } else {
          as.numeric(input$min_cutoff_feature_merge)
        }
      } else {
        "q0"
      }
      
      max_cut <- if (!is.null(input$max_cutoff_feature_merge) && 
                     !is.na(input$max_cutoff_feature_merge) && 
                     input$max_cutoff_feature_merge != "") {
        if (grepl("^q", input$max_cutoff_feature_merge)) {
          input$max_cutoff_feature_merge
        } else {
          as.numeric(input$max_cutoff_feature_merge)
        }
      } else {
        NA
      }
      
      message(paste("FeaturePlot cutoffs: min =", min_cut, ", max =", max_cut))
      
      # Check if 3D mode
      if(!is.null(input$enable_3d_feature_merge) && input$enable_3d_feature_merge) {
        # 3D MODE
        if(!"umap3d" %in% names(seurat_object@reductions)) {
          showNotification("3D UMAP not found. Please run 'Run UMAP 3D' first.", type = "warning", duration = 5)
          return()
        }
        
        if(length(present_genes) > 3) {
          showNotification("Maximum 3 genes for 3D. Using first 3.", type = "warning", duration = 3)
          present_genes <- present_genes[1:3]
        }
        
        hide_grid <- !is.null(input$hide_grid_feature_merge) && input$hide_grid_feature_merge
        hide_axes <- !is.null(input$add_noaxes_feature_merge) && input$add_noaxes_feature_merge
        dark_mode <- !is.null(input$dark_mode_feature_plot_merge) && input$dark_mode_feature_plot_merge
        
        combined_plot <- create3DFeaturePlot(
          seurat_object = seurat_object,
          genes = present_genes,
          reduction = "umap3d",
          pt_size = 3,
          min_cutoff = min_cut,
          max_cutoff = max_cut,
          order = TRUE,
          hide_grid = hide_grid,
          hide_axes = hide_axes,
          dark_mode = dark_mode
        )
        
        attr(combined_plot, "n_rows") <- 1
        feature_plot_merge(combined_plot)
        
      } else {
        # 2D MODE
        base_theme <- theme(
          plot.title = element_text(size = input$title_text_size_merge %||% 14, face = "bold"),
          axis.text.x = element_text(size = input$axis_text_size_merge %||% 12),
          axis.text.y = element_text(size = input$axis_text_size_merge %||% 12),
          axis.title = element_text(size = input$title_text_size_merge %||% 14),
          legend.text = element_text(size = input$axis_text_size_merge %||% 12),
          axis.line = element_line(linewidth = input$axis_line_width_merge %||% 0.5),
          axis.ticks = element_line(linewidth = input$axis_line_width_merge %||% 0.5)
        )
        
        # Determine split option
        split_option <- input$metadata_to_compare
        use_split <- !is.null(split_option) && split_option != ""
        
        # Blend mode requires exactly 2 genes and is incompatible with split.by
        use_blend <- isTRUE(input$show_coexpression_merge) && length(present_genes) == 2
        
        # Warn and disable split if blend is active — Seurat does not support both
        if (use_blend && use_split) {
          showNotification(
            "Blend mode is not compatible with split.by — split disabled for dual expression.",
            type = "warning", duration = 5
          )
          use_split <- FALSE
        }
        
        # Number of rows/columns — not applicable in blend mode
        n_genes <- length(present_genes)
        n_cols <- if(n_genes == 1) 1 else if(n_genes == 2) 2 else 3
        n_rows <- ceiling(n_genes / n_cols)
        
        if(use_split) {
          combined_plot <- FeaturePlot(
            seurat_object,
            features = present_genes,
            split.by = split_option,
            order = TRUE,
            min.cutoff = min_cut,
            max.cutoff = max_cut,
            ncol = n_cols
          ) + base_theme
        } else if (use_blend) {
          # Blend mode: ncol is handled internally by Seurat (returns 3 panels)
          combined_plot <- FeaturePlot(
            seurat_object,
            features = present_genes,
            blend = TRUE,
            blend.threshold = 1,
            order = TRUE,
            min.cutoff = min_cut,
            max.cutoff = max_cut
          ) + base_theme
          # Blend always produces 3 panels side by side
          n_rows <- 1
        } else {
          combined_plot <- FeaturePlot(
            seurat_object,
            features = present_genes,
            order = TRUE,
            min.cutoff = min_cut,
            max.cutoff = max_cut,
            ncol = n_cols
          ) + base_theme
        }
        
        # Add conditional modifications
        if(!is.null(input$add_noaxes_feature_merge) && input$add_noaxes_feature_merge) {
          combined_plot <- combined_plot + NoAxes()
        }
        if(!is.null(input$add_nolegend_feature_merge) && input$add_nolegend_feature_merge) {
          combined_plot <- combined_plot + NoLegend()
        }
        
        # Dark mode
        if(!is.null(input$dark_mode_feature_plot_merge) && input$dark_mode_feature_plot_merge) {
          combined_plot <- combined_plot +
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
        
        attr(combined_plot, "n_rows") <- n_rows
        feature_plot_merge(combined_plot)
      }
      
    }, error = function(e) {
      showNotification(paste("Error in FeaturePlot:", e$message), type = "error")
      message(paste("FeaturePlot error:", e$message))
    })
  })

  # ---- Plots ----

  output$FeaturePlot2 <- renderPlot({
    req(feature_plot_merge())
    plot_obj <- feature_plot_merge()
    if(inherits(plot_obj, "gg") || inherits(plot_obj, "patchwork")) {
      txt_face <- if (isTRUE(input$bold_axis_text_merge)) "bold" else "plain"
      styled <- plot_obj +
        theme(
          axis.text    = element_text(size = input$axis_text_size_merge, face = txt_face),
          axis.title   = element_text(size = input$title_text_size_merge, face = txt_face),
          plot.title   = element_text(size = input$title_text_size_merge, face = txt_face),
          strip.text   = element_text(size = input$title_text_size_merge, face = txt_face),
          legend.text  = element_text(size = input$axis_text_size_merge),
          legend.title = element_text(size = input$axis_text_size_merge, face = txt_face),
          axis.line    = element_line(linewidth = input$axis_line_width_merge),
          axis.ticks   = element_line(linewidth = input$axis_line_width_merge)
        )
      if(!is.null(input$dark_mode_feature_plot_merge) && isTRUE(input$dark_mode_feature_plot_merge)) {
        styled <- styled & 
          theme(
            panel.background = element_rect(fill = "black"),
            plot.background = element_rect(fill = "black"),
            panel.grid = element_blank(),
            text = element_text(color = "white"),
            axis.text = element_text(color = "white"),
            legend.background = element_rect(fill = "black"),
            legend.text = element_text(color = "white")
          )
      }
      styled
    }
  }, height = function() {
    plot_obj <- feature_plot_merge()
    if(!is.null(plot_obj)) {
      n_rows <- attr(plot_obj, "n_rows") %||% 1
      max(400, n_rows * 350)
    } else {
      400
    }
  })

  output$FeaturePlot2_3d <- renderPlotly({
    req(feature_plot_merge())
    plot_obj <- feature_plot_merge()
    if(inherits(plot_obj, "plotly")) {
      plot_obj
    }
  })

  ############################## DOT PLOT ##############################

  # ---- Observers and computations ----

  # DotPlot observer - using modular function
  observeEvent(input$runDotPlot, {
    tryCatch({
      seurat_object <- multiple_datasets_object()
      if (is.null(seurat_object)) { showNotification("Load or integrate datasets first.", type = "warning"); return() }
      # Check that group_by is selected
      if (is.null(input$group_by_select) || input$group_by_select == "") {
        showNotification("Please select a grouping column (Group By) first", type = "warning", duration = 5)
        return()
      }
      DefaultAssay(seurat_object) <- input$viz_assay_merge
      # Empty field and unknown genes are reported; genes are checked in the selected assay
      check <- validate_gene_input(input$gene_list_dot_merge, seurat_object, input$viz_assay_merge)
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      genes <- check$genes
      message("=== DotPlot ===")
      message("group_by: ", input$group_by_select)
      message("split.by: ", input$metadata_to_compare %||% "none")
      message("Initial cells: ", ncol(seurat_object))
      # STEP 1: Filter by split column values if specified
      split_by_column <- NULL
      if (!is.null(input$metadata_to_compare) && input$metadata_to_compare != "") {
        split_by_column <- input$metadata_to_compare
        if (!is.null(input$split_values_filter_gene_viz) && length(input$split_values_filter_gene_viz) > 0) {
          cells_to_keep <- which(seurat_object@meta.data[[split_by_column]] %in% input$split_values_filter_gene_viz)
          if (length(cells_to_keep) > 0) {
            seurat_object <- subset(seurat_object, cells = cells_to_keep)
            message(paste("Filtered to", length(cells_to_keep), "cells by split values"))
          } else {
            showNotification("No cells match the selected split filter values", type = "warning")
            return()
          }
        }
      }
      # STEP 2: Determine group_by column
      group_by <- input$group_by_select
      # STEP 3: Filter and reorder by cluster_order
      if (!is.null(input$cluster_order_dotplot_merge) && length(input$cluster_order_dotplot_merge) > 0) {
        if (group_by == "seurat_clusters") {
          cells_to_keep <- as.character(Idents(seurat_object)) %in% input$cluster_order_dotplot_merge
        } else {
          cells_to_keep <- seurat_object@meta.data[[group_by]] %in% input$cluster_order_dotplot_merge
        }
        # Selection made before a renaming no longer matches any cell: report it instead of a subset() error
        if (sum(cells_to_keep, na.rm = TRUE) == 0) {
          showNotification("The selected groups do not match the current values of the Group By column (renamed?). Reselect them.", type = "error", duration = 8)
          return()
        }
        seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[which(cells_to_keep)])
        message(paste("Filtered to", ncol(seurat_object), "cells by group_by values"))
        # Force factor levels after subset
        if (group_by == "seurat_clusters") {
          Idents(seurat_object) <- factor(as.character(Idents(seurat_object)), levels = input$cluster_order_dotplot_merge)
        } else {
          seurat_object@meta.data[[group_by]] <- factor(as.character(seurat_object@meta.data[[group_by]]), levels = input$cluster_order_dotplot_merge)
        }
        message("cluster_order_dotplot_merge input: ", paste(input$cluster_order_dotplot_merge, collapse = ", "))
        message("Levels after reorder: ", paste(
          if (group_by == "seurat_clusters") levels(Idents(seurat_object)) else levels(seurat_object@meta.data[[group_by]]),
          collapse = ", "))
      }
      # Check we still have cells
      if (ncol(seurat_object) == 0) {
        showNotification("No cells found with selected parameters", type = "error")
        return()
      }
      message(paste("Final cells:", ncol(seurat_object)))
      # STEP 4: Create DotPlot using modular function (JoinLayers handled inside the helper)
      palette_name <- input$color_palette_dotplot %||% "default"
      dot_scale <- input$dot_scale_dotplot %||% 1
      plot <- create_dotplot(
        seurat_object = seurat_object,
        genes = genes,
        group_by = if (group_by == "seurat_clusters") NULL else group_by,
        split_by = split_by_column,
        color_palette = palette_name,
        dot_scale = dot_scale,
        assay = input$viz_assay_merge,
        scale_expression = isTRUE(input$dotplot_scale_merge),
        rotate_axis = FALSE
      )
      # STEP 5: Apply legend customization
      if (input$add_nolegend_dot_merge) {
        plot <- plot + guides(size = "none")
        message("Removed size legend (kept color gradient legend)")
      }
      # STEP 6: Apply axes customization
      if (input$add_noaxes_dot_merge) {
        plot <- plot + NoAxes()
        message("Removed axes")
      }
      # STEP 7: Invert axes if requested
      if (input$invert_axes_merge) {
        plot <- plot + coord_flip()
        message("Inverted axes")
      }
      dot_plot_merge(plot)
      message("DotPlot successfully created")
    }, error = function(e) {
      showNotification(paste("Error in DotPlot:", e$message), type = "error")
      message(paste("DotPlot error:", e$message))
    })
  })

  # ---- Plots ----

  # Render DotPlot — theme applied here so font/line inputs are reactive
  output$DotPlot2 <- renderPlot({
    req(dot_plot_merge())
    txt_face <- if (isTRUE(input$bold_axis_text_merge)) "bold" else "plain"
    dot_plot_merge() +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = input$axis_text_size_merge, face = txt_face),
        axis.text.y  = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        plot.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        legend.text  = element_text(size = input$axis_text_size_merge),
        legend.title = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.line    = element_line(linewidth = input$axis_line_width_merge),
        axis.ticks   = element_line(linewidth = input$axis_line_width_merge)
      )
  })

  ############################## VIOLIN PLOT ##############################

  # ---- Observers and computations ----

  observeEvent(input$runVlnPlot, {
    tryCatch({
      seurat_object <- multiple_datasets_object()
      if (is.null(seurat_object)) { showNotification("Load or integrate datasets first.", type = "warning"); return() }
      # Check that group_by is selected
      if (is.null(input$group_by_select) || input$group_by_select == "") {
        showNotification("Please select a grouping column (Group By) first", type = "warning", duration = 5)
        return()
      }
      DefaultAssay(seurat_object) <- input$viz_assay_merge
      # Empty field and unknown genes are reported; genes are checked in the selected assay
      check <- validate_gene_input(input$gene_list_vln_merge, seurat_object, input$viz_assay_merge)
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      genes <- check$genes
      
      message("=== VlnPlot ===")
      message("group_by: ", input$group_by_select)
      message("split.by: ", input$metadata_to_compare %||% "none")
      message("Initial cells: ", ncol(seurat_object))
      
      # STEP 1: Filter by split column values if specified
      split_by_column <- NULL
      if (!is.null(input$metadata_to_compare) && input$metadata_to_compare != "") {
        split_by_column <- input$metadata_to_compare
        
        if (!is.null(input$split_values_filter_gene_viz) && length(input$split_values_filter_gene_viz) > 0) {
          cells_to_keep <- which(seurat_object@meta.data[[split_by_column]] %in% input$split_values_filter_gene_viz)
          
          if (length(cells_to_keep) > 0) {
            seurat_object <- subset(seurat_object, cells = cells_to_keep)
            message(paste("Filtered to", length(cells_to_keep), "cells by split values"))
          } else {
            showNotification("No cells match the selected split filter values", type = "warning")
            return()
          }
        }
      }
      
      # STEP 2: Determine group_by column
      group_by <- input$group_by_select
      
      # STEP 3: Filter by group_by values (cluster_order) if specified
      if (!is.null(input$cluster_order_vln_merge) && length(input$cluster_order_vln_merge) > 0) {
        
        # Determine which cells to keep based on group_by column
        if (group_by == "seurat_clusters") {
          cells_to_keep <- as.character(Idents(seurat_object)) %in% input$cluster_order_vln_merge
        } else {
          cells_to_keep <- seurat_object@meta.data[[group_by]] %in% input$cluster_order_vln_merge
        }
        
        # Selection made before a renaming no longer matches any cell: report it instead of a subset() error
        if (sum(cells_to_keep, na.rm = TRUE) == 0) {
          showNotification("The selected groups do not match the current values of the Group By column (renamed?). Reselect them.", type = "error", duration = 8)
          return()
        }
        seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[which(cells_to_keep)])
        message(paste("Filtered to", ncol(seurat_object), "cells by group_by values"))
        
        # Set factor levels for ordering - CRITICAL: Don't touch Idents for non-seurat_clusters
        if (group_by == "seurat_clusters") {
          Idents(seurat_object) <- factor(Idents(seurat_object), levels = input$cluster_order_vln_merge)
        } else {
          seurat_object@meta.data[[group_by]] <- factor(seurat_object@meta.data[[group_by]], levels = input$cluster_order_vln_merge)
        }
      }
      
      # Check we still have cells
      if (ncol(seurat_object) == 0) {
        showNotification("No cells found with selected parameters", type = "error")
        return()
      }
      
      message(paste("Final cells:", ncol(seurat_object)))
      message(paste("group.by will be:", group_by))
      
      # Determine point size
      pt_size <- if (input$hide_vln_points_merge) 0 else 0.1
      
      # STEP 4: Create VlnPlot - Handle split.by
      if (!is.null(split_by_column)) {
        
        if (length(genes) > 1) {
          plot_list <- lapply(genes, function(gene) {
            if (group_by == "seurat_clusters") {
              p <- VlnPlot(seurat_object, features = gene, split.by = split_by_column, pt.size = pt_size)
            } else {
              p <- VlnPlot(seurat_object, features = gene, group.by = group_by, split.by = split_by_column, pt.size = pt_size)
            }
            if (input$add_noaxes_vln_merge) p <- p + NoAxes()
            return(p)
          })
          
          plot <- wrap_plots(plot_list, ncol = 2)
          attr(plot, "n_rows") <- ceiling(length(genes) / 2)
          
        } else {
          if (group_by == "seurat_clusters") {
            plot <- VlnPlot(seurat_object, features = genes, split.by = split_by_column, pt.size = pt_size)
          } else {
            plot <- VlnPlot(seurat_object, features = genes, group.by = group_by, split.by = split_by_column, pt.size = pt_size)
          }
          
          if (input$add_noaxes_vln_merge) plot <- plot + NoAxes()
          attr(plot, "n_rows") <- 1
        }
        
        message("VlnPlot created with split.by")
        
      } else {
        
        if (group_by == "seurat_clusters") {
          seurat_object$plot_group <- as.character(Idents(seurat_object))
          plot_group_by <- "plot_group"
        } else {
          plot_group_by <- group_by
        }
        
        plot <- createManualVlnPlotMultiple(seurat_object = seurat_object, genes = genes, group_by = plot_group_by, pt_size = pt_size, add_noaxes = input$add_noaxes_vln_merge, add_nolegend = input$add_nolegend_vln_merge, ncol = 2)
        
        message("VlnPlot created without split (manual function)")
      }
      
      vln_plot_merge(plot)
      message("VlnPlot successfully created")
      
    }, error = function(e) {
      showNotification(paste("Error in VlnPlot:", e$message), type = "error")
      message(paste("VlnPlot error:", e$message))
    })
  })

  # ---- Plots ----

  output$VlnPlot2 <- renderPlot({
    req(vln_plot_merge())
    txt_face <- if (isTRUE(input$bold_axis_text_merge)) "bold" else "plain"
    vln_plot_merge() +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = input$axis_text_size_merge, face = txt_face),
        axis.text.y  = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        plot.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        strip.text   = element_text(size = input$title_text_size_merge, face = txt_face),
        legend.text  = element_text(size = input$axis_text_size_merge),
        legend.title = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.line    = element_line(linewidth = input$axis_line_width_merge),
        axis.ticks   = element_line(linewidth = input$axis_line_width_merge)
      )
  }, height = function() {
    plot_obj <- vln_plot_merge()
    if (!is.null(plot_obj)) {
      n_rows <- attr(plot_obj, "n_rows") %||% 1
      max(400, n_rows * 400)
    } else {
      400
    }
  })

  ############################## RIDGE PLOT ##############################

  # ---- Observers and computations ----

  observeEvent(input$runRidgePlot, {
    tryCatch({
      seurat_object <- multiple_datasets_object()
      if (is.null(seurat_object)) { showNotification("Load or integrate datasets first.", type = "warning"); return() }
      # Check that group_by is selected
      if (is.null(input$group_by_select) || input$group_by_select == "") {
        showNotification("Please select a grouping column (Group By) first", type = "warning", duration = 5)
        return()
      }
      DefaultAssay(seurat_object) <- input$viz_assay_merge
      # Empty field and unknown genes are reported; genes are checked in the selected assay
      check <- validate_gene_input(input$gene_list_ridge_merge, seurat_object, input$viz_assay_merge)
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      genes <- check$genes
      
      # Filter by split column values if specified
      if (!is.null(input$metadata_to_compare) && 
          input$metadata_to_compare != "" &&
          !is.null(input$split_values_filter_gene_viz) && 
          length(input$split_values_filter_gene_viz) > 0) {
        
        cells_to_keep <- which(seurat_object@meta.data[[input$metadata_to_compare]] %in% input$split_values_filter_gene_viz)
        
        if (length(cells_to_keep) > 0) {
          seurat_object <- subset(seurat_object, cells = cells_to_keep)
          message(paste("Filtered to", length(cells_to_keep), "cells for RidgePlot"))
        } else {
          showNotification("No cells match the selected filter values", type = "warning")
          return()
        }
      }
      
      if(input$group_by_select != "ident" && 
         !input$group_by_select %in% colnames(seurat_object@meta.data)) {
        showNotification(paste("Column", input$group_by_select, "not found in metadata."), type = "error")
        return()
      }
      
      if(ncol(seurat_object) == 0) {
        showNotification("No cells found with selected parameters", type = "error")
        return()
      }
      
      plot <- createManualRidgePlotMultiple(
        seurat_object = seurat_object,
        genes = genes,
        group_by = input$group_by_select,
        add_noaxes = !is.null(input$add_noaxes_ridge_merge) && input$add_noaxes_ridge_merge,
        add_nolegend = !is.null(input$add_nolegend_ridge_merge) && input$add_nolegend_ridge_merge
      )
      
      ridge_plot_merge(plot)
      message("RidgePlot generated successfully")
      
    }, error = function(e) {
      showNotification(paste("Error in RidgePlot:", e$message), type = "error")
      message(paste("RidgePlot error:", e$message))
    })
  })

  # ---- Plots ----

  output$Ridge_plot_merge <- renderPlot({
    req(ridge_plot_merge())
    txt_face <- if (isTRUE(input$bold_axis_text_merge)) "bold" else "plain"
    ridge_plot_merge() +
      theme(
        axis.text.x  = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.text.y  = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        plot.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        strip.text   = element_text(size = input$title_text_size_merge, face = txt_face),
        legend.text  = element_text(size = input$axis_text_size_merge),
        legend.title = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.line    = element_line(linewidth = input$axis_line_width_merge),
        axis.ticks   = element_line(linewidth = input$axis_line_width_merge)
      )
  })

  # ---- Downloads ----

  # Download handler for VlnPlot Merge
  output$downloadVlnPlotMerge <- createDownloadHandler(
    reactive_data = vln_plot_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets") 
    }),
    data_name = reactive({ 
      paste0("VlnPlot_", paste(parse_gene_tokens(input$gene_list_vln_merge), collapse="_")) 
    }),
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$plot_format_merge }),
      width = 15, 
      height = reactive({ min(49, (attr(vln_plot_merge(), "n_rows") %||% 1) * 8) }),
      dpi = reactive({ input$dpi_input_merge })
    ),
    plot_modifier = function(p) {
      txt_face <- if (isTRUE(input$bold_axis_text_merge)) "bold" else "plain"
      p + theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = input$axis_text_size_merge, face = txt_face),
        axis.text.y  = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        plot.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        strip.text   = element_text(size = input$title_text_size_merge, face = txt_face),
        legend.text  = element_text(size = input$axis_text_size_merge),
        legend.title = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.line    = element_line(linewidth = input$axis_line_width_merge),
        axis.ticks   = element_line(linewidth = input$axis_line_width_merge)
      )
    }
  )

  # Download handler for FeaturePlot Merge
  output$downloadFeaturePlotMerge <- createDownloadHandler(
    reactive_data = feature_plot_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets") 
    }),
    data_name = reactive({ 
      paste0("FeaturePlot_", paste(parse_gene_tokens(input$gene_list_feature_merge), collapse="_")) 
    }),
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$plot_format_merge }),
      width = reactive({ min(12 * length(trimws(strsplit(input$gene_list_feature_merge, ",")[[1]])), 24) }), 
      height = reactive({ min(49, (attr(feature_plot_merge(), "n_rows") %||% 1) * 8) }),
      dpi = reactive({ input$dpi_input_merge })
    ),
    plot_modifier = function(p) {
      txt_face <- if (isTRUE(input$bold_axis_text_merge)) "bold" else "plain"
      p + theme(
        axis.text    = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        plot.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        legend.text  = element_text(size = input$axis_text_size_merge),
        legend.title = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.line    = element_line(linewidth = input$axis_line_width_merge),
        axis.ticks   = element_line(linewidth = input$axis_line_width_merge)
      )
    }
  )

  # Download handler for DotPlot Merge
  output$downloadDotPlotMerge <- createDownloadHandler(
    reactive_data = dot_plot_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets") 
    }),
    data_name = reactive({ 
      paste0("DotPlot_", paste(parse_gene_tokens(input$gene_list_dot_merge), collapse="_")) 
    }),
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$plot_format_merge }),
      width = reactive({ min(12 + length(trimws(strsplit(input$gene_list_dot_merge, ",")[[1]])) * 0.5, 24) }), 
      height = 8, 
      dpi = reactive({ input$dpi_input_merge })
    ),
    plot_modifier = function(p) {
      txt_face  <- if (isTRUE(input$bold_axis_text_merge)) "bold" else "plain"
      gene_face <- if (isTRUE(input$bold_axis_text_merge)) "bold.italic" else "italic"
      p + theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = input$axis_text_size_merge, face = txt_face),
        axis.text.y  = element_text(size = input$axis_text_size_merge, face = gene_face),
        axis.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        plot.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        legend.text  = element_text(size = input$axis_text_size_merge),
        legend.title = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.line    = element_line(linewidth = input$axis_line_width_merge),
        axis.ticks   = element_line(linewidth = input$axis_line_width_merge)
      )
    }
  )

  # Download handler for RidgePlot Merge
  output$downloadRidgePlotMerge <- createDownloadHandler(
    reactive_data = ridge_plot_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets") 
    }),
    data_name = reactive({ 
      paste0("RidgePlot_", paste(parse_gene_tokens(input$gene_list_ridge_merge), collapse="_")) 
    }),
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$plot_format_merge }),
      width = reactive({ min(10 + length(trimws(strsplit(input$gene_list_ridge_merge, ",")[[1]])) * 0.5, 20) }), 
      height = 8, 
      dpi = reactive({ input$dpi_input_merge })
    ),
    plot_modifier = function(p) {
      txt_face <- if (isTRUE(input$bold_axis_text_merge)) "bold" else "plain"
      p + theme(
        axis.text.x  = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.text.y  = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        plot.title   = element_text(size = input$title_text_size_merge, face = txt_face),
        strip.text   = element_text(size = input$title_text_size_merge, face = txt_face),
        legend.text  = element_text(size = input$axis_text_size_merge),
        legend.title = element_text(size = input$axis_text_size_merge, face = txt_face),
        axis.line    = element_line(linewidth = input$axis_line_width_merge),
        axis.ticks   = element_line(linewidth = input$axis_line_width_merge)
      )
    }
  )

  ################# Number of nuclei expressing a gene ###################

  # ---- Reactive state ----
  number_of_nuclei_merge <- reactiveVal(NULL)

  # ---- Observers and computations ----

  # Observer pour l'analyse d'expression dans les datasets intégrés
  observeEvent(input$analyze_btn_genes_expression_merge, {
    tryCatch({
      if (is.null(multiple_datasets_object())) { showNotification("Load or integrate datasets first.", type = "warning"); return() }
      check <- validate_gene_input(input$gene_list_genes_expression_merge, multiple_datasets_object(), input$viz_assay_merge)
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      genes_vector <- check$genes
      
      message(paste("Analyzing genes:", paste(genes_vector, collapse = ", ")))
      
      # Now pass the vector (not the string) to the function
      result <- analyze_gene_expression(
        seurat_obj = multiple_datasets_object(),
        selected_genes = genes_vector,  # ✅ Pass vector, not string
        assay_name = input$viz_assay_merge,
        expression_threshold = input$logfc_threshold_genes_expression_merge %||% 0.1,
        is_integrated = TRUE
      )
      
      # Store results
      number_of_nuclei_merge(result$data)
      
      # Display table
      output$expression_summary_merge <- renderDT({
        render_expression_table(result, "expression_summary_merge")
      })
      
      showNotification(
        paste("Analysis completed for", length(genes_vector), "genes"),
        type = "message"
      )
      
    }, error = function(e) {
      showNotification(paste("Error processing expression data:", e$message), type = "error")
      message("Gene expression error: ", e$message)
    })
  })

  # ---- Downloads ----

  # Pour le téléchargement des données
  output$download_genes_number_expression_merge <- createDownloadHandler(
    reactive_data = number_of_nuclei_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "IntegratedExpression") 
    }),
    data_name = "integrated_expression_analysis",
    download_type = "csv"
  )

  ############################## Heatmap and dual expression multi dataset ##############################

  # ---- Reactive state ----
  heatmap_plot_multidataset <- reactiveVal()
  heatmap_n_genes_multi <- reactiveVal(0)  # Track number of genes for dynamic height

  # Store top genes table for CSV export (NULL when using manual gene input)
  heatmap_markers_data_multi <- reactiveVal(NULL)

  # ---- Observers and computations ----

  # Update selectors when data changes; selections still valid in the new object are kept
  observeEvent(multiple_datasets_object(), run_safely("Interface update", persistent = FALSE, {
    req(multiple_datasets_object())
    seurat_obj <- multiple_datasets_object()
    # Update clusters
    clusters <- levels(droplevels(Idents(seurat_obj)))
    current_clusters <- intersect(isolate(input$clusters_heatmap_multi), clusters)
    updateSelectizeInput(session, "clusters_heatmap_multi", 
                         choices = clusters, 
                         selected = if (length(current_clusters) > 0) current_clusters else clusters,
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
    current_split <- isolate(input$split_by_heatmap_multi)
    updateSelectInput(session, "split_by_heatmap_multi", 
                      choices = c("None" = "None", filtered_cols),
                      selected = if (!is.null(current_split) && current_split %in% filtered_cols) current_split else "None")
  }))

  # Update the Group 1 / Group 2 choices when the split column changes
  observeEvent(input$split_by_heatmap_multi, run_safely("split by heatmap multi", persistent = FALSE, {
    req(multiple_datasets_object())
    split_by <- input$split_by_heatmap_multi
    values <- character(0)
    if (!is.null(split_by) && !split_by %in% c("", "None") && split_by %in% colnames(multiple_datasets_object()@meta.data)) {
      values <- sort(unique(trimws(as.character(multiple_datasets_object()@meta.data[[split_by]]))))
      values <- values[!is.na(values) & values != ""]
    }
    updateSelectizeInput(session, "heatmap_group1_multi", choices = values, selected = character(0), server = TRUE)
    updateSelectizeInput(session, "heatmap_group2_multi", choices = values, selected = character(0), server = TRUE)
  }))

  # Generate heatmap: columns = clusters, or cluster x split value with the second axis drawn as panels
  observeEvent(input$generateCustomHeatmap_multi, run_safely("Heatmap (multiple datasets)", {
    # Configuration checked before the modal opens; run_safely closes the modal on any later error
    if (is.null(multiple_datasets_object())) stop("Load or integrate datasets first")
    gene_source <- input$gene_source_heatmap_multi %||% "manual"
    split_by <- input$split_by_heatmap_multi
    if (is.null(split_by) || split_by %in% c("", "None")) split_by <- NULL
    group1 <- input$heatmap_group1_multi
    group2 <- input$heatmap_group2_multi
    if (gene_source == "condition_de") {
      if (is.null(split_by)) stop("Choose in 'Split by' the column holding the conditions (e.g. genotype)")
      if (length(group1) == 0 || length(group2) == 0) stop("Select at least one value for Group 1 and for Group 2")
      if (length(intersect(group1, group2)) > 0) stop("A value cannot be in both groups: ", paste(intersect(group1, group2), collapse = ", "))
    }
    showModal(modalDialog(title = "Generating heatmap...", "Please wait while generating average expression heatmap", easyClose = FALSE, footer = NULL))
    seurat_object <- multiple_datasets_object()
    # Average expression and DE read the RNA assay (log-normalized counts); the interface has no assay selector
    selected_assay <- if ("RNA" %in% names(seurat_object@assays)) "RNA" else DefaultAssay(seurat_object)
    all_clusters <- levels(droplevels(Idents(seurat_object)))
    clusters_selected <- input$clusters_heatmap_multi
    if (length(clusters_selected) == 0) clusters_selected <- all_clusters
    clusters_selected <- clusters_selected[clusters_selected %in% all_clusters]
    if (length(clusters_selected) == 0) stop("None of the selected clusters exists in the current object")
    n_top <- max(1L, as.integer(input$n_top_genes_multi %||% 10))
    # Cells displayed: selected clusters and, with a split column, the chosen values (both groups in condition mode)
    keep <- as.character(Idents(seurat_object)) %in% clusters_selected
    split_levels <- NULL
    if (!is.null(split_by)) {
      split_vals <- trimws(as.character(seurat_object@meta.data[[split_by]]))
      split_levels <- if (gene_source == "condition_de") c(group1, group2) else sort(unique(split_vals[keep & !is.na(split_vals) & split_vals != ""]))
      keep <- keep & split_vals %in% split_levels
    }
    if (!any(keep)) stop("No cells left for the selected clusters and groups")
    obj <- prepare_rna_assay_for_de(subset(seurat_object, cells = colnames(seurat_object)[keep]), selected_assay)$obj
    # Gene selection
    if (gene_source == "cluster_markers") {
      # Markers of each cluster against all other clusters of the full object
      markers_result <- get_top_markers(seurat_object, n_genes = n_top, assay = selected_assay, clusters_to_include = clusters_selected)
      genes <- markers_result$genes
      heatmap_markers_data_multi(markers_result$markers_table)
      subtitle <- sprintf("Top %d markers per cluster (vs all other clusters)", n_top)
    } else if (gene_source == "condition_de") {
      de_result <- get_top_condition_de_genes(obj, split_by, group1, group2, clusters = clusters_selected, n_genes = n_top,
                                              direction = input$heatmap_de_direction_multi %||% "both", assay = selected_assay)
      genes <- de_result$genes
      heatmap_markers_data_multi(de_result$markers_table)
      if (length(de_result$skipped) > 0) showNotification(paste("Clusters without result:", paste(de_result$skipped, collapse = "; ")), type = "warning", duration = 15)
      subtitle <- sprintf("Top %d DE genes per cluster, %s vs %s (cell-level Wilcoxon, adj. p < 0.05)", n_top, paste(group1, collapse = "+"), paste(group2, collapse = "+"))
    } else {
      gene_result <- validate_gene_input(input$gene_select_heatmap_multi, obj, selected_assay)
      if (!gene_result$ok) stop(gene_result$message)
      if (!is.null(gene_result$message)) showNotification(gene_result$message, type = "warning", duration = 8)
      genes <- gene_result$genes
      heatmap_markers_data_multi(NULL)
      subtitle <- NULL
    }
    if (length(genes) == 0) stop("No gene to display")
    # Column layout: one column per cluster, or per cluster x split value with panels on the primary axis
    cluster_vals <- as.character(Idents(obj))
    order_by <- input$heatmap_order_multi %||% "cluster"
    if (is.null(split_by)) {
      col_id <- cluster_vals
      ids <- clusters_selected[clusters_selected %in% col_id]
      col_info <- data.frame(id = ids, x = ids, facet = NA_character_, stringsAsFactors = FALSE)
      x_title <- "Cluster"
    } else {
      col_id <- paste(cluster_vals, trimws(as.character(obj@meta.data[[split_by]])), sep = "___")
      # expand.grid varies its first column fastest: split values inside each cluster, or clusters inside each split value
      grid <- if (order_by == "cluster") expand.grid(split = split_levels, cluster = clusters_selected, stringsAsFactors = FALSE) else expand.grid(cluster = clusters_selected, split = split_levels, stringsAsFactors = FALSE)
      grid$id <- paste(grid$cluster, grid$split, sep = "___")
      grid <- grid[grid$id %in% col_id, , drop = FALSE]
      col_info <- if (order_by == "cluster") data.frame(id = grid$id, x = grid$split, facet = grid$cluster, stringsAsFactors = FALSE) else data.frame(id = grid$id, x = grid$cluster, facet = grid$split, stringsAsFactors = FALSE)
      x_title <- if (order_by == "cluster") split_by else "Cluster"
    }
    avg <- compute_group_average_expression(obj, genes, factor(col_id, levels = col_info$id), assay = selected_assay)
    scale_mode <- input$heatmap_scale_multi %||% "zscore"
    if (scale_mode == "zscore" && nrow(col_info) == 2) showNotification("Row z-score over 2 columns only shows the direction (every value is +/-0.71). Choose 'Average expression' to compare magnitudes.", type = "warning", duration = 12)
    plot <- plot_grouped_expression_heatmap(avg, col_info, scale_mode = scale_mode, color_palette = input$color_palette_heatmap_multi %||% "viridis",
                                            title = paste("Heatmap -", selected_assay, "assay"), subtitle = subtitle, x_title = x_title)
    heatmap_n_genes_multi(length(genes))
    heatmap_plot_multidataset(plot)
    removeModal()
    showNotification(sprintf("Heatmap with %d genes and %d columns", length(genes), nrow(col_info)), type = "message")
  }))

  # ---- Plots ----

  # Render heatmap with dynamic height
  output$heatmap_plot_multi <- renderPlot({
    req(heatmap_plot_multidataset())
    heatmap_plot_multidataset()
  }, height = function() {
    n_genes <- heatmap_n_genes_multi()
    if (n_genes == 0) return(600)
    max(400, min(20 * n_genes, 2000))
  })

  # ---- Downloads ----

  # Download plot handler with dynamic height
  output$download_heatmap_multi <- createDownloadHandler(
    reactive_data = heatmap_plot_multidataset,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MultiDatasets") 
    }),
    data_name = "heatmap",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$format_heatmap_multi %||% "png" }), 
      width = 10, 
      height = reactive({
        n_genes <- heatmap_n_genes_multi()
        if (n_genes == 0) return(8)
        max(6, min(0.15 * n_genes, 50))
      }), 
      dpi = reactive({ input$dpi_heatmap_multi %||% 300 })
    )
  )

  # Download top markers table as CSV (only populated in top markers mode)
  output$download_heatmap_markers_multi <- downloadHandler(
    filename = function() {
      paste0("top_markers_per_cluster_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      markers_df <- heatmap_markers_data_multi()
      req(markers_df)
      write.csv(markers_df, file, row.names = FALSE)
    }
  )

  ############################## Scatter Plot Multi Dataset ##############################

  # ---- Reactive state ----
  scatter_plot_multi <- reactiveVal()

  # Variable réactive pour stocker le plot actuel
  current_plot_merge <- reactiveVal()

  # ---- Observers and computations ----

  # Update color_by and split_by choices when data changes
  observeEvent(multiple_datasets_object(), run_safely("Interface update", persistent = FALSE, {
    req(multiple_datasets_object())
    seurat_obj <- multiple_datasets_object()
    
    # Filter metadata columns
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
    
    # Update color_by
    updateSelectInput(session, "color_by_scatter_multi", 
                      choices = c("Cluster" = "cluster", filtered_cols),
                      selected = "cluster")
    
    # Update split_by
    updateSelectInput(session, "split_by_scatter_multi", 
                      choices = c("None" = "None", filtered_cols),
                      selected = "None")
  }))

  # Generate scatter plot
  observeEvent(input$generateScatter_multi, run_safely("generateScatter multi", {
    req(multiple_datasets_object())
    
    showModal(modalDialog(
      title = "Processing",
      div(h4("Generating Scatter Plot...", style = "text-align: center;")),
      footer = NULL,
      easyClose = FALSE
    ))
    
    seurat_obj <- multiple_datasets_object()
    
    tryCatch({
      # One gene per field, case-insensitive (the first name is used if several are typed)
      selected_assay <- input$assay_select_scatter_multi %||% "RNA"
      check1 <- validate_gene_input(input$feature1_select_multi, seurat_obj, selected_assay)
      check2 <- validate_gene_input(input$feature2_select_multi, seurat_obj, selected_assay)
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
      cluster_text <- trimws(input$scatter_text_clusters_multi)
      if (cluster_text == "") {
        clusters_to_use <- NULL
      } else {
        clusters_to_use <- trimws(unlist(strsplit(cluster_text, ",")))
        clusters_to_use <- clusters_to_use[nchar(clusters_to_use) > 0]
        if (length(clusters_to_use) == 0) clusters_to_use <- NULL
      }
      
      # Get display options
      color_by <- input$color_by_scatter_multi
      if (is.null(color_by)) color_by <- "cluster"
      
      split_by <- input$split_by_scatter_multi
      if (is.null(split_by) || split_by == "None") split_by <- NULL
      
      # Get thresholds
      threshold1 <- input$threshold_gene1_multi
      threshold2 <- input$threshold_gene2_multi
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
      
      scatter_plot_multi(plot)
      showNotification("Scatter plot generated successfully.", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      message(paste("Full error:", e$message))
    })
    
    removeModal()
  }))

  # Observer les changements dans la sélection du type de plot et mettre à jour le plot
  observe(run_safely("Interface update", persistent = FALSE, {
    plot_type <- input$plot_type_select_merge
    if (plot_type == "FeaturePlot") {
      current_plot_merge(feature_plot_merge())
    } else if (plot_type == "VlnPlot") {
      current_plot_merge(vln_plot_merge())
    } else if (plot_type == "DotPlot") {
      current_plot_merge(dot_plot_merge())
    } else if (plot_type == "RidgePlot") {
      current_plot_merge(ridge_plot_merge())
    }
  }))

  # ---- Plots ----

  # Render scatter plot
  output$scatter_plot_multi <- renderPlot({
    req(scatter_plot_multi())
    scatter_plot_multi()
  })

  # Afficher le plot sélectionné
  output$selected_plot_display_merge <- renderPlot({
    req(current_plot_merge())
    current_plot_merge()
  })

  # ---- Downloads ----

  # Download handler
  output$download_scatter_multi <- createDownloadHandler(
    reactive_data = scatter_plot_multi,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MultiDatasets_FeatureScatter")
    }),
    data_name = "scatter",
    download_type = "plot",
    plot_params = list(
      file_type = reactive(input$format_scatter_multi), 
      width = 10,
      height = 8,
      dpi = reactive(input$dpi_scatter_multi)  
    )
  )

  ############################## Assigning cell identity merge ##############################

  # ---- Reactive state ----

  # Reactive variable for that tab
  cluster_colours_merge <- reactiveVal()

  # ---- Helpers ----

  # Sauvegarder l'état
  save_cluster_state_merge <- function(seurat_obj, action_description = "Cluster modification") {
    current_clusters <- as.character(Idents(seurat_obj))
    names(current_clusters) <- names(Idents(seurat_obj))
    
    new_state <- list(
      clusters = current_clusters,
      action = action_description
    )
    
    cluster_history_merge$states <- append(cluster_history_merge$states, list(new_state))
    
    if(length(cluster_history_merge$states) > cluster_history_merge$max_history) {
      cluster_history_merge$states <- cluster_history_merge$states[-1]
    }
  }

  # Undo: restore the last snapshot (taken before the last change) and drop it; the first entry is the initial state and is kept
  undo_last_action <- function() {
    n_states <- length(cluster_history_merge$states)
    if (n_states <= 1) return(NULL)
    restored <- cluster_history_merge$states[[n_states]]$clusters
    cluster_history_merge$states <- cluster_history_merge$states[-n_states]
    restored
  }

  # Function to apply updated colors to Seurat object and update reactive value
  apply_cluster_colours <- function(seurat_object) {
    cluster_colors <- cluster_colours_merge()
    if (!is.null(cluster_colors)) {
      seurat_object@misc$cluster_colors <- cluster_colors
    }
    return(seurat_object)
  }

  # ---- Observers and computations ----

  # Central cluster observer — updates rename/color selectors when object or column changes
  observe(run_safely("Interface update", persistent = FALSE, {
    req(multiple_datasets_object())
    col <- if (!is.null(input$filtered_umap_col_merge) &&
               nzchar(input$filtered_umap_col_merge) &&
               input$filtered_umap_col_merge %in% colnames(multiple_datasets_object()@meta.data)) {
      input$filtered_umap_col_merge
    } else {
      get_cluster_column(multiple_datasets_object())
    }
    updateClusterSelectors(session, multiple_datasets_object(),
                           cluster_col = col,
                           rename_id   = "select_cluster_merge",
                           color_id    = "select_color_merge")
    get_cluster_colors_merge(multiple_datasets_object())
  }))

  # Update color picker to reflect the stored color of the selected cluster
  observeEvent(input$select_color_merge, run_safely("select color merge", {
    req(multiple_datasets_object(), input$select_color_merge)
    colors <- get_cluster_colors_merge(multiple_datasets_object())
    if (!is.null(colors) && input$select_color_merge %in% names(colors)) {
      updateColourInput(session, "select_cluster_merge_color", value = colors[[input$select_color_merge]])
    }
  }), ignoreInit = TRUE)

  observe(run_safely("Interface update", persistent = FALSE, {
    req(multiple_datasets_object())
    forced <- isolate(forced_col_merge())
    col <- if (!is.null(forced)) {
      forced_col_merge(NULL)
      forced
    } else {
      isolate(input$filtered_umap_col_merge)
    }
    meta_cols <- colnames(multiple_datasets_object()@meta.data)
    selected <- if (!is.null(col) && nzchar(col) && col %in% meta_cols) col else NULL
    update_umap_selectors(session, multiple_datasets_object(), module = "merge",
                          selected_col = selected)
  }))

  observeEvent(input$rename_single_cluster_merge_button, run_safely("rename single cluster merge button", {
    req(input$select_cluster_merge, input$rename_single_cluster_merge, multiple_datasets_object())
    updated_seurat <- multiple_datasets_object()
    new_name <- trimws(input$rename_single_cluster_merge)
    req(nzchar(new_name))
    old_idents <- as.character(Idents(updated_seurat))
    # The rename acts on the active identities: a value of another column picked through 'Color by' matches no cell
    if (!input$select_cluster_merge %in% old_idents) {
      showNotification(sprintf("'%s' is not a current cluster identity. Set 'Color by' to the cluster column (e.g. ClusterIdents) and select the cluster again.", input$select_cluster_merge), type = "error", duration = 10)
      return()
    }
    save_cluster_state_merge(updated_seurat)
    merged <- new_name %in% old_idents
    new_idents <- old_idents
    new_idents[old_idents == input$select_cluster_merge] <- new_name
    # Cluster order kept: the renamed cluster takes the place of the old one (or joins the existing one when merged)
    old_levels <- levels(Idents(updated_seurat))
    new_levels <- unique(ifelse(old_levels == input$select_cluster_merge, new_name, old_levels))
    Idents(updated_seurat) <- factor(stats::setNames(new_idents, colnames(updated_seurat)), levels = new_levels[new_levels %in% new_idents])
    # Renamed labels always go to ClusterIdents, the stable column read by the cluster selectors
    updated_seurat$ClusterIdents <- new_idents
    # The 'Color by' column is synchronized only when it held the identities before the rename (a cluster column);
    # any other column (dataset, genotype, QC metric) is never overwritten
    active_col <- input$filtered_umap_col_merge
    if (!is.null(active_col) && nzchar(active_col) && active_col != "ClusterIdents" && active_col %in% colnames(updated_seurat@meta.data) &&
        identical(as.character(updated_seurat@meta.data[[active_col]]), old_idents)) {
      updated_seurat@meta.data[[active_col]] <- new_idents
    }
    updateClusterSelectors(session, updated_seurat,
                           cluster_col = "ClusterIdents",
                           rename_id   = "select_cluster_merge",
                           color_id    = "select_color_merge")
    # Force the observer to select ClusterIdents when it fires after object update
    forced_col_merge("ClusterIdents")
    multiple_datasets_object(updated_seurat)
    showNotification(if (merged) paste("Clusters merged under the name:", new_name) else paste("Cluster renamed to:", new_name), type = "message")
  }))

  observeEvent(input$update_colour_merge_button, run_safely("update colour merge button", {
    message("Update the color of the selected cluster for multiple datasets")
    updated_seurat <- multiple_datasets_object()
    cluster_colors <- get_cluster_colors_merge(updated_seurat)
    
    if (is.null(cluster_colors)) {
      message("Initializing cluster colors with Seurat defaults")
      all_clusters   <- levels(Idents(updated_seurat))
      n_clusters     <- length(all_clusters)
      default_colors <- scales::hue_pal()(n_clusters)
      cluster_colors <- setNames(default_colors, all_clusters)
      message(paste("Initialized", n_clusters, "cluster colors"))
    }
    
    if (is.null(input$select_color_merge) ||
        !input$select_color_merge %in% names(cluster_colors)) {
      showNotification(
        paste("Selected cluster", input$select_color_merge, "is not valid.",
              "Available clusters:", paste(names(cluster_colors), collapse = ", ")),
        type = "error", duration = 10
      )
      return()
    }
    
    cluster_colors[input$select_color_merge] <- input$select_cluster_merge_color
    message(paste("Updating color for cluster", input$select_color_merge,
                  "to", input$select_cluster_merge_color))
    
    updated_seurat@misc$cluster_colors <- cluster_colors
    multiple_datasets_object(updated_seurat)
    message("Current cluster colors:")
    print(updated_seurat@misc$cluster_colors)
    showNotification(
      paste("Cluster", input$select_color_merge, "color updated successfully"),
      type = "message"
    )
  }))

  ############################## Automated annotation (reference label transfer + SingleR) ##############################

  # ---- Reactive state ----
  transfer_result_merge <- reactiveVal(NULL)
  singler_result_merge <- reactiveVal(NULL)
  annotation_reference_merge <- reactiveVal(NULL)

  # celldex references are cached on disk (pre-filled in the Docker image) so they are downloaded only once
  singler_cache_dir_merge <- Sys.getenv("CELLHUB_REF_CACHE", file.path(tools::R_user_dir("CellHub", "cache"), "singler_references"))

  # ---- Helpers ----

  # Cluster vector used for the per-cluster vote (same column as the cluster selectors; "ident" falls back to Idents)
  annotation_clusters_merge <- function(obj) {
    cluster_col <- active_cluster_column_merge()
    if (is.null(cluster_col) || !cluster_col %in% colnames(obj@meta.data)) cluster_col <- NULL
    clusters <- if (!is.null(cluster_col)) as.character(obj@meta.data[[cluster_col]]) else as.character(Idents(obj))
    list(column = cluster_col, values = clusters)
  }

  # Cluster summary recomputed with the method's current agreement threshold (key = "transfer" or "singler")
  annotation_summary_merge <- function(result, key) {
    min_agreement <- input[[paste0(key, "_min_agreement_merge")]] %||% 50
    summarize_annotation_by_cluster(result$cells, min_agreement / 100, result$params$max_unassigned)
  }

  # Rename clusters from an annotation result (undoable, like manual renaming)
  apply_annotation_labels_merge <- function(result, method_name, key) {
    obj <- multiple_datasets_object()
    cluster_col <- result$params$cluster_col
    current <- if (!is.null(cluster_col) && cluster_col %in% colnames(obj@meta.data)) {
      as.character(obj@meta.data[[cluster_col]])
    } else {
      as.character(Idents(obj))
    }
    new_names <- build_annotation_cluster_names(annotation_summary_merge(result, key),
                                                keep_low_confidence = isTRUE(input[[paste0(key, "_keep_low_merge")]]),
                                                make_unique = isTRUE(input[[paste0(key, "_unique_merge")]]))
    labels <- unname(new_names[current])
    labels[is.na(labels)] <- current[is.na(labels)]
    save_cluster_state_merge(obj, paste(method_name, "labels"))
    Idents(obj) <- factor(stats::setNames(labels, colnames(obj)))
    obj$ClusterIdents <- labels
    get_cluster_colors_merge(obj)
    updateClusterSelectors(session, obj,
                           cluster_col = "ClusterIdents",
                           rename_id   = "select_cluster_merge",
                           color_id    = "select_color_merge")
    forced_col_merge("ClusterIdents")
    multiple_datasets_object(obj)
    showNotification(sprintf("%s labels applied to %d of %d clusters", method_name,
                             sum(new_names != names(new_names)), length(new_names)), type = "message", duration = 6)
  }

  # Cluster summary table with confidence colouring
  render_annotation_summary_merge <- function(result, key) {
    df <- annotation_summary_merge(result, key)
    colnames(df) <- c("Cluster", "Cells", "Predicted label", "Agreement (%)", "2nd label", "2nd (%)",
                      "Unassigned (%)", "Median score", "Median delta", "Confidence")
    datatable(df, rownames = FALSE, options = list(pageLength = 25, scrollX = TRUE)) %>%
      formatStyle("Confidence", target = "row", backgroundColor = styleEqual(c("High", "Low"), c("#e8f5e9", "#ffebee")))
  }

  # Status line under the Run button: not run yet, or a one-line summary of the last run
  render_annotation_status_merge <- function(result, key, extra = NULL) {
    if (is.null(result)) {
      return(tags$p(icon("circle-info"), " Not run yet.", style = "margin-top: 10px; font-size: 12px; color: #6c757d;"))
    }
    s <- annotation_summary_merge(result, key)
    n_low <- sum(s$confidence == "Low")
    tags$div(
      style = "margin-top: 10px; padding: 8px 10px; border-radius: 6px; font-size: 12px; background-color: #e8f5e9; color: #2e7d32;",
      icon("circle-check"), sprintf(" %d clusters annotated: %d confident, %d to check.", nrow(s), nrow(s) - n_low, n_low),
      tags$div(sprintf("Cluster column used for the vote: %s", result$params$cluster_col %||% "active identities"), style = "color: #6c757d; margin-top: 4px;"),
      if (!is.null(extra)) tags$div(extra, style = "color: #6c757d; margin-top: 2px;")
    )
  }

  # UMAP coloured by a per-cell annotation column or by its confidence score (same reduction as the final UMAP)
  plot_annotation_umap_merge <- function(column, title, continuous = FALSE) {
    obj <- multiple_datasets_object()
    reduc <- if (!is.null(input$umap_reduction_merge) && nzchar(input$umap_reduction_merge) &&
                 input$umap_reduction_merge %in% names(obj@reductions)) input$umap_reduction_merge else get_umap_reduction(obj)
    shiny::validate(shiny::need(!is.null(reduc), "No UMAP found: run the clustering step first."))
    req(column %in% colnames(obj@meta.data))
    if (continuous) {
      FeaturePlot(obj, features = column, reduction = reduc, order = TRUE) + ggtitle(title)
    } else {
      DimPlot(obj, reduction = reduc, group.by = column, label = TRUE, repel = TRUE, label.size = 3) + ggtitle(title)
    }
  }

  # ---- Observers and computations ----

  # A new merged/integrated object invalidates previous annotation results and the undo history of the previous object
  observeEvent(list(input$load_seurat_file_merge, input$integrate, input$integrate_harmony, input$integrate_fastmnn, input$integrate_rpca, input$simple_merge), run_safely("Action", {
    singler_result_merge(NULL)
    transfer_result_merge(NULL)
    cluster_history_merge$states <- list()
  }), ignoreInit = TRUE)

  ############ Label transfer from a user reference ############

  # ---- Observers and computations ----

  observeEvent(input$annotation_reference_file_merge, run_safely("annotation reference file merge", {
    showModal(modalDialog(title = "Loading reference", "Reading the annotated reference dataset...", easyClose = FALSE, footer = NULL))
    tryCatch({
      ref <- read_reference_object(input$annotation_reference_file_merge$datapath)
      cols <- get_reference_label_columns(ref)
      if (length(cols) == 0) stop("No categorical metadata column with 2-300 values found in the reference")
      preferred <- intersect(c("cell_type", "celltype", "CellType", "cell_types", "annotation", "ClusterIdents", "cluster_names"), cols)
      annotation_reference_merge(ref)
      updateSelectInput(session, "annotation_reference_label_merge", choices = cols,
                        selected = if (length(preferred) > 0) preferred[1] else cols[1])
      transfer_result_merge(NULL)
      removeModal()
      showNotification(sprintf("Reference loaded: %s cells, %s genes", format(ncol(ref), big.mark = ","), format(nrow(ref), big.mark = ",")),
                       type = "message")
    }, error = function(e) {
      removeModal()
      annotation_reference_merge(NULL)
      showNotification(paste("Error loading reference:", conditionMessage(e)[1]), type = "error", duration = 10)
    })
  }))

  observeEvent(input$run_transfer_merge, run_safely("run transfer merge", {
    req(multiple_datasets_object())
    if (is.null(annotation_reference_merge())) {
      showNotification("Load an annotated reference first (step 1).", type = "warning")
      return()
    }
    req(input$annotation_reference_label_merge)
    showModal(modalDialog(title = "Label transfer from reference",
                          "Finding anchors between the reference and your integrated dataset... This can take several minutes.",
                          easyClose = FALSE, footer = NULL))
    tryCatch({
      obj <- multiple_datasets_object()
      cl <- annotation_clusters_merge(obj)
      ref <- prepare_reference_object(annotation_reference_merge(), input$annotation_reference_label_merge,
                                      max_cells_per_label = input$transfer_max_cells_merge)
      result <- run_label_transfer(obj, ref, label_column = input$annotation_reference_label_merge, clusters = cl$values,
                                   dims = input$transfer_dims_merge, min_score = input$transfer_min_score_merge,
                                   min_agreement = input$transfer_min_agreement_merge / 100)
      result$params$reference <- input$annotation_reference_file_merge$name
      result$params$cluster_col <- cl$column
      multiple_datasets_object(add_annotation_metadata(obj, result, prefix = "RefTransfer"))
      transfer_result_merge(result)
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste("Label transfer error:", conditionMessage(e)[1]), type = "error", duration = 10)
      message("Label transfer error (merge): ", conditionMessage(e)[1])
    })
  }))

  observeEvent(input$apply_transfer_labels_merge, run_safely("apply transfer labels merge", {
    req(multiple_datasets_object())
    if (is.null(transfer_result_merge())) {
      showNotification("Run the label transfer first (step 2).", type = "warning")
      return()
    }
    tryCatch(apply_annotation_labels_merge(transfer_result_merge(), "Reference", "transfer"),
             error = function(e) showNotification(paste("Error applying labels:", conditionMessage(e)[1]), type = "error", duration = 10))
  }))

  # ---- UI and text outputs ----

  # Reference summary: size and the labels found in the chosen column
  output$annotation_reference_info_merge <- renderUI({
    ref <- annotation_reference_merge()
    if (is.null(ref)) {
      return(tags$p(icon("circle-info"), " No reference loaded.", style = "font-size: 12px; color: #6c757d;"))
    }
    col <- input$annotation_reference_label_merge
    labels <- if (!is.null(col) && col %in% colnames(ref@meta.data)) sort(unique(na.omit(as.character(ref@meta.data[[col]])))) else character(0)
    tags$div(
      style = "font-size: 12px; background-color: #f3f0ff; padding: 8px 10px; border-radius: 6px;",
      tags$p(tags$strong(format(ncol(ref), big.mark = ",")), " cells, ", tags$strong(format(nrow(ref), big.mark = ",")), " genes, ",
             tags$strong(length(labels)), " labels", style = "margin: 0 0 4px 0;"),
      tags$p(paste(head(labels, 12), collapse = ", "), if (length(labels) > 12) sprintf(" ... (+%d)", length(labels) - 12),
             style = "margin: 0; color: #6c757d;")
    )
  })

  output$transfer_status_merge <- renderUI({
    res <- transfer_result_merge()
    extra <- if (!is.null(res)) sprintf("%s anchors, %s reference cells, %s shared genes", format(res$params$n_anchors, big.mark = ","),
                                        format(res$params$n_reference_cells, big.mark = ","), format(res$params$n_common_genes, big.mark = ",")) else NULL
    render_annotation_status_merge(res, "transfer", extra)
  })

  # ---- Plots ----

  output$transfer_umap_labels_merge <- renderPlot({
    req(transfer_result_merge())
    plot_annotation_umap_merge("RefTransfer_pruned_label", "Transferred labels (per cell)")
  })

  output$transfer_umap_score_merge <- renderPlot({
    req(transfer_result_merge())
    plot_annotation_umap_merge("RefTransfer_score", "Prediction score (max)", continuous = TRUE)
  })

  # ---- Tables ----

  output$transfer_summary_merge <- renderDT({
    req(transfer_result_merge())
    render_annotation_summary_merge(transfer_result_merge(), "transfer")
  })

  # ---- Downloads ----

  output$download_transfer_clusters_merge <- downloadHandler(
    filename = function() paste0(getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets"),
                                 "_RefTransfer_cluster_summary_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      req(transfer_result_merge())
      res <- transfer_result_merge()
      out <- annotation_summary_merge(res, "transfer")
      out$reference <- res$params$reference
      out$label_column <- res$params$label_column
      out$min_score <- res$params$min_score
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$download_transfer_cells_merge <- downloadHandler(
    filename = function() paste0(getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets"),
                                 "_RefTransfer_per_cell_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      req(transfer_result_merge())
      write.csv(transfer_result_merge()$cells, file, row.names = FALSE)
    }
  )

  ############ SingleR (built-in celldex references) ############

  # ---- Observers and computations ----

  # Reference choices follow the selected species
  observeEvent(input$singler_species_merge, run_safely("singler species merge", {
    refs <- get_singler_references(input$singler_species_merge)
    updateSelectInput(session, "singler_reference_merge", choices = refs, selected = if (length(refs) > 0) refs[[1]] else NULL)
  }))

  # Default the SingleR species to the one chosen in the loading modal
  observeEvent(input$species_choice_merge, run_safely("species choice merge", {
    if (input$species_choice_merge %in% c("human", "mouse")) {
      updateSelectInput(session, "singler_species_merge", selected = input$species_choice_merge)
    }
  }))

  observeEvent(input$run_singler_merge, run_safely("run singler merge", {
    req(multiple_datasets_object(), input$singler_reference_merge)
    showModal(modalDialog(title = "Automated annotation (SingleR)",
                          "Loading the reference and annotating cells... This can take several minutes for large integrated datasets.",
                          easyClose = FALSE, footer = NULL))
    tryCatch({
      obj <- multiple_datasets_object()
      cl <- annotation_clusters_merge(obj)
      reference <- load_singler_reference(input$singler_reference_merge, singler_cache_dir_merge)
      result <- run_singler_annotation(obj, reference, label_column = input$singler_label_level_merge, clusters = cl$values,
                                       mode = input$singler_mode_merge, assay = annotation_assay_merge(obj),
                                       min_agreement = input$singler_min_agreement_merge / 100)
      result$params$reference <- input$singler_reference_merge
      result$params$cluster_col <- cl$column
      multiple_datasets_object(add_annotation_metadata(obj, result, prefix = "SingleR"))
      singler_result_merge(result)
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste("SingleR error:", conditionMessage(e)[1]), type = "error", duration = 10)
      message("SingleR error (merge): ", conditionMessage(e)[1])
    })
  }))

  observeEvent(input$apply_singler_labels_merge, run_safely("apply singler labels merge", {
    req(multiple_datasets_object())
    if (is.null(singler_result_merge())) {
      showNotification("Run SingleR first (step 2).", type = "warning")
      return()
    }
    tryCatch(apply_annotation_labels_merge(singler_result_merge(), "SingleR", "singler"),
             error = function(e) showNotification(paste("Error applying labels:", conditionMessage(e)[1]), type = "error", duration = 10))
  }))

  # Compute 3D UMAP from the same (integrated) reduction and dimensions as the 2D UMAP
  observeEvent(input$compute_3d_umap_merge, run_safely("compute 3d umap merge", {
    req(multiple_datasets_object())
    showModal(modalDialog(title = "Computing 3D UMAP", "Calculating 3D UMAP coordinates... This may take a moment.", easyClose = FALSE, footer = NULL))
    tryCatch({
      seurat_obj <- multiple_datasets_object()
      reduction <- get_integration_reduction(seurat_obj)
      if (is.null(reduction)) stop("No PCA or integrated reduction found: run the integration and Scaling & PCA steps first")
      n_dims <- if (!is.null(input$dimension_2) && !is.na(input$dimension_2)) input$dimension_2 else 30
      seurat_obj <- runUMAP3D_reproducible(object = seurat_obj, dims = 1:n_dims, reduction = reduction, seed = 42, verbose = TRUE)
      multiple_datasets_object(seurat_obj)
      removeModal()
      showNotification(sprintf("3D UMAP computed (%s, %d dimensions): shown below the clustering plot.", reduction, n_dims), type = "message")
      }, error = function(e) {
      removeModal()
      showNotification(paste("Error computing 3D UMAP:", conditionMessage(e)[1]), type = "error")
      message("Error in 3D UMAP computation: ", conditionMessage(e)[1])
    })
  }))

  # Observer pour Undo
  observeEvent(input$undo_cluster_merge, run_safely("undo cluster merge", {
    req(multiple_datasets_object())
    restored_clusters <- undo_last_action()
    if (is.null(restored_clusters)) {
      showNotification("Nothing to undo", type = "warning")
      return()
    }
    updated_seurat <- multiple_datasets_object()
    # A snapshot taken on another object (different cells) cannot be restored on this one
    if (!setequal(names(restored_clusters), colnames(updated_seurat))) {
      cluster_history_merge$states <- list()
      showNotification("The undo history belonged to a previous object and was cleared.", type = "warning")
      return()
    }
    old_idents <- as.character(Idents(updated_seurat))
    restored <- unname(restored_clusters[colnames(updated_seurat)])
    Idents(updated_seurat) <- factor(stats::setNames(restored, colnames(updated_seurat)))
    updated_seurat$ClusterIdents <- restored
    # Same rule as the rename: only a 'Color by' column that held the identities is synchronized
    active_col <- input$filtered_umap_col_merge
    if (!is.null(active_col) && nzchar(active_col) && active_col != "ClusterIdents" && active_col %in% colnames(updated_seurat@meta.data) &&
        identical(as.character(updated_seurat@meta.data[[active_col]]), old_idents)) {
      updated_seurat@meta.data[[active_col]] <- restored
    }
    updateClusterSelectors(session, updated_seurat,
                           cluster_col = "ClusterIdents",
                           rename_id   = "select_cluster_merge",
                           color_id    = "select_color_merge")
    forced_col_merge("ClusterIdents")
    multiple_datasets_object(updated_seurat)
    showNotification("Undid last action", type = "message")
  }))

  # ---- UI and text outputs ----

  output$singler_status_merge <- renderUI({
    res <- singler_result_merge()
    extra <- if (!is.null(res)) sprintf("%s, %s shared genes", res$params$reference, format(res$params$n_common_genes, big.mark = ",")) else NULL
    render_annotation_status_merge(res, "singler", extra)
  })

  # Flag used by the UI to show the 3D panel only once a 3D UMAP exists
  output$has_umap3d_merge <- reactive({
    obj <- multiple_datasets_object()
    !is.null(obj) && "umap3d" %in% names(obj@reductions)
  })

  outputOptions(output, "has_umap3d_merge", suspendWhenHidden = FALSE)

  # ---- Plots ----

  output$singler_umap_labels_merge <- renderPlot({
    req(singler_result_merge())
    plot_annotation_umap_merge("SingleR_pruned_label", "SingleR labels (per cell)")
  })

  output$singler_umap_delta_merge <- renderPlot({
    req(singler_result_merge())
    plot_annotation_umap_merge("SingleR_delta", "SingleR confidence (delta to next-best label)", continuous = TRUE)
  })

  # Interactive 3D UMAP in the Clustering step, coloured by the current clusters
  output$clustering_umap3d_merge <- renderPlotly({
    req(multiple_datasets_object())
    obj <- multiple_datasets_object()
    req("umap3d" %in% names(obj@reductions))
    build_umap3d_plotly(obj, cluster_values = as.character(Idents(obj)), cluster_colors = get_cluster_colors_merge(obj),
                        pt_size = if (isTRUE(input$umap3d_pt_size_clust_merge > 0)) input$umap3d_pt_size_clust_merge else 2, label_size = 14,
                        hide_grid = isTRUE(input$umap3d_hide_grid_clust_merge),
                        hide_axes = isTRUE(input$umap3d_hide_axes_clust_merge),
                        dark_mode = isTRUE(input$umap3d_dark_clust_merge))
  })

  # Display final UMAP with 2D/3D support
  output$umap_finale_merge <- renderPlotly({
    req(multiple_datasets_object())
    seurat_obj  <- multiple_datasets_object()
    umap_reduc <- if (!is.null(input$umap_reduction_merge) &&
                      nzchar(input$umap_reduction_merge) &&
                      input$umap_reduction_merge %in% names(seurat_obj@reductions)) {
      input$umap_reduction_merge
    } else {
      get_umap_reduction(seurat_obj)
    }
    if (is.null(umap_reduc)) return(NULL)
    message("Generating the final UMAP for multiple datasets — reduction: ", umap_reduc)
    use_3d         <- isTRUE(input$umap_3d_toggle_merge)
    cluster_colors <- get_cluster_colors_merge(seurat_obj)
    cluster_col <- if (!is.null(input$umap_finale_col_merge) &&
                       nzchar(input$umap_finale_col_merge) &&
                       input$umap_finale_col_merge %in% colnames(seurat_obj@meta.data)) {
      input$umap_finale_col_merge
    } else {
      get_cluster_column(seurat_obj)
    }
    cluster_values <- if (cluster_col == "ident") {
      as.character(Idents(seurat_obj))
    } else {
      as.character(seurat_obj@meta.data[[cluster_col]])
    }
    if (use_3d && hasValid3DUMAP(seurat_obj)) {
      umap_coords <- as.data.frame(Embeddings(seurat_obj, reduction = "umap3d"))
      colnames(umap_coords) <- c("UMAP_1", "UMAP_2", "UMAP_3")
      umap_coords$cluster <- cluster_values
      if (!all(levels(as.factor(umap_coords$cluster)) %in% names(cluster_colors))) {
        missing_clusters <- setdiff(levels(as.factor(umap_coords$cluster)), names(cluster_colors))
        default_colors <- scales::hue_pal()(length(missing_clusters))
        names(default_colors) <- missing_clusters
        cluster_colors <- c(cluster_colors, default_colors)
      }
      centers <- data.frame(
        cluster = levels(as.factor(umap_coords$cluster)),
        UMAP_1 = vapply(levels(as.factor(umap_coords$cluster)), function(cl) {
          median(umap_coords$UMAP_1[umap_coords$cluster == cl], na.rm = TRUE)
        }, numeric(1)),
        UMAP_2 = vapply(levels(as.factor(umap_coords$cluster)), function(cl) {
          median(umap_coords$UMAP_2[umap_coords$cluster == cl], na.rm = TRUE)
        }, numeric(1)),
        UMAP_3 = vapply(levels(as.factor(umap_coords$cluster)), function(cl) {
          median(umap_coords$UMAP_3[umap_coords$cluster == cl], na.rm = TRUE)
        }, numeric(1))
      )
      hide_grid  <- isTRUE(input$umap_3d_hide_grid_merge)
      hide_axes  <- isTRUE(input$umap_3d_hide_axes_merge)
      dark_mode  <- isTRUE(input$umap_3d_dark_mode_merge)
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
                  color = ~cluster, colors = cluster_colors, type = "scatter3d", mode = "markers",
                  marker = list(size = input$pt_size_merge, opacity = 0.7),
                  hoverinfo = "text", text = ~paste("Cluster:", cluster)) %>%
        add_trace(data = centers, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
                  type = "scatter3d", mode = "text", text = ~cluster,
                  textfont = list(size = input$label_font_size_merge * 4, color = text_color),
                  hoverinfo = "none", showlegend = FALSE) %>%
        layout(
          scene = list(xaxis = xaxis_config, yaxis = yaxis_config, zaxis = zaxis_config,
                       bgcolor = bg_color, camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))),
          paper_bgcolor = bg_color, plot_bgcolor = bg_color, showlegend = FALSE,
          title = list(text = input$plot_title_merge, font = list(size = 24, color = text_color))
        )
      message("3D UMAP generated")
      return(p)
    } else if (use_3d && !hasValid3DUMAP(seurat_obj)) {
      return(plotly_empty() %>%
               layout(title = list(text = "Please compute 3D UMAP first using the button above",
                                   font = list(size = 16))))
    } else {
      umap_data <- as.data.frame(Embeddings(seurat_obj, reduction = umap_reduc))
      message("UMAP column names: ", paste(colnames(umap_data), collapse = ", "))
      if (ncol(umap_data) < 2) {
        showNotification("UMAP reduction has incorrect dimensions", type = "error")
        return(NULL)
      }
      colnames(umap_data)[1:2] <- c("UMAP_1", "UMAP_2")
      umap_data$cluster <- cluster_values
      if (!all(levels(as.factor(umap_data$cluster)) %in% names(cluster_colors))) {
        missing_clusters <- setdiff(levels(as.factor(umap_data$cluster)), names(cluster_colors))
        default_colors <- scales::hue_pal()(length(missing_clusters))
        names(default_colors) <- missing_clusters
        cluster_colors <- c(cluster_colors, default_colors)
      }
      centers <- data.frame(
        cluster = levels(as.factor(umap_data$cluster)),
        UMAP_1 = vapply(levels(as.factor(umap_data$cluster)), function(cl) {
          median(umap_data$UMAP_1[umap_data$cluster == cl], na.rm = TRUE)
        }, numeric(1)),
        UMAP_2 = vapply(levels(as.factor(umap_data$cluster)), function(cl) {
          median(umap_data$UMAP_2[umap_data$cluster == cl], na.rm = TRUE)
        }, numeric(1))
      )
      p <- plot_ly() %>%
        add_trace(data = umap_data, x = ~UMAP_1, y = ~UMAP_2,
                  color = ~cluster, colors = cluster_colors,
                  type = "scatter", mode = "markers",
                  marker = list(size = input$pt_size_merge, opacity = 0.7),
                  hoverinfo = "text", text = ~paste("Cluster:", cluster)) %>%
        layout(title = list(text = input$plot_title_merge, font = list(size = 24)),
               xaxis = list(title = "UMAP 1", zeroline = FALSE),
               yaxis = list(title = "UMAP 2", zeroline = FALSE),
               hovermode = "closest", showlegend = FALSE)
      for (i in seq_len(nrow(centers))) {
        p <- p %>% add_annotations(
          x = centers$UMAP_1[i], y = centers$UMAP_2[i],
          text = as.character(centers$cluster[i]),
          showarrow = FALSE,
          font = list(size = input$label_font_size_merge * 4, color = "black")
        )
      }
      message("2D UMAP generated — reduction: ", umap_reduc, " — column: ", cluster_col)
      return(p)
    }
  })

  # ---- Tables ----

  output$singler_summary_merge <- renderDT({
    req(singler_result_merge())
    render_annotation_summary_merge(singler_result_merge(), "singler")
  })

  # ---- Downloads ----

  output$download_singler_clusters_merge <- downloadHandler(
    filename = function() paste0(getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets"),
                                 "_SingleR_cluster_summary_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      req(singler_result_merge())
      res <- singler_result_merge()
      out <- annotation_summary_merge(res, "singler")
      out$reference <- res$params$reference
      out$label_column <- res$params$label_column
      out$mode <- res$params$mode
      write.csv(out, file, row.names = FALSE)
    }
  )

  output$download_singler_cells_merge <- downloadHandler(
    filename = function() paste0(getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets"),
                                 "_SingleR_per_cell_", format(Sys.Date(), "%Y%m%d"), ".csv"),
    content = function(file) {
      req(singler_result_merge())
      write.csv(singler_result_merge()$cells, file, row.names = FALSE)
    }
  )

  ############################# Differential Gene Expression Analysis  #############################

  # ---- Reactive state ----
  filtered_umap_plot <- reactiveVal(NULL)          # Stores the filtered UMAP plot

  # ---- Observers and computations ----

  # Disable download buttons at startup
  shinyjs::disable("download_markers_single_cluster_merge")

  shinyjs::disable("download_markers_multiple_clusters_merge")
  shinyjs::disable("download_diff_dataset_cluster")
  shinyjs::disable("download_venn_diagram")

  # Reactive column name backing the cluster selector — same resolution logic used everywhere
  # this dropdown is populated, so calculate_DE can group.by the exact same column (fixes
  # "Cannot find the following identities" when the color-by selector points to a non-cluster
  # metadata column, e.g. a condition/genotype column).
  active_cluster_column_merge <- reactive({
    req(multiple_datasets_object())
    obj <- multiple_datasets_object()
    col_input <- if (!is.null(input$filtered_umap_col_merge) && nzchar(input$filtered_umap_col_merge)) {
      input$filtered_umap_col_merge
    } else if (!is.null(input$umap_finale_col_merge) && nzchar(input$umap_finale_col_merge)) {
      input$umap_finale_col_merge
    } else {
      NULL
    }
    if (!is.null(col_input) && col_input %in% colnames(obj@meta.data)) {
      col_input
    } else {
      get_cluster_column(obj)
    }
  })

  # Reactive cluster list — reads either color selector (kept in sync via observers below).
  active_clusters_merge <- reactive({
    req(multiple_datasets_object())
    obj <- multiple_datasets_object()
    col <- active_cluster_column_merge()
    if (!is.null(col) && col %in% colnames(obj@meta.data)) {
      unique(as.character(obj@meta.data[[col]]))
    } else {
      as.character(unique(Idents(obj)))
    }
  })

  # Keep the two "Color by" selectors in sync across tabs (umap_finale_col_merge ↔ filtered_umap_col_merge).
  observeEvent(input$umap_finale_col_merge, run_safely("umap finale col merge", {
    req(nzchar(input$umap_finale_col_merge))
    if (!identical(input$filtered_umap_col_merge, input$umap_finale_col_merge))
      updateSelectInput(session, "filtered_umap_col_merge", selected = input$umap_finale_col_merge)
  }), ignoreInit = TRUE)

  observeEvent(input$filtered_umap_col_merge, run_safely("filtered umap col merge", {
    req(nzchar(input$filtered_umap_col_merge))
    if (!identical(input$umap_finale_col_merge, input$filtered_umap_col_merge))
      updateSelectInput(session, "umap_finale_col_merge", selected = input$filtered_umap_col_merge)
  }), ignoreInit = TRUE)

  observeEvent(input$umap_reduction_merge, run_safely("umap reduction merge", {
    req(nzchar(input$umap_reduction_merge))
    if (!identical(input$filtered_umap_reduction_merge, input$umap_reduction_merge))
      updateSelectInput(session, "filtered_umap_reduction_merge", selected = input$umap_reduction_merge)
    if (!identical(input$cluster_cells_reduction_merge, input$umap_reduction_merge))
      updateSelectInput(session, "cluster_cells_reduction_merge", selected = input$umap_reduction_merge)
  }), ignoreInit = TRUE)

  observeEvent(input$filtered_umap_reduction_merge, run_safely("filtered umap reduction merge", {
    req(nzchar(input$filtered_umap_reduction_merge))
    if (!identical(input$umap_reduction_merge, input$filtered_umap_reduction_merge))
      updateSelectInput(session, "umap_reduction_merge", selected = input$filtered_umap_reduction_merge)
  }), ignoreInit = TRUE)

  # ---- UI and text outputs ----

  # Dataset filter for the filtered UMAP (values read by output$filtered_umap_plot and its download handler)
  output$dataset_filter_ui <- renderUI({
    req(multiple_datasets_object())
    obj <- multiple_datasets_object()
    if (!"dataset" %in% colnames(obj@meta.data)) return(NULL)
    datasets <- sort(unique(as.character(obj@meta.data$dataset)))
    current <- isolate(input$dataset_filter)
    selected <- if (!is.null(current) && all(current %in% c("All Datasets", datasets))) current else "All Datasets"
    selectInput("dataset_filter", "Datasets shown:", choices = c("All Datasets", datasets), selected = selected, multiple = TRUE)
  })

  # ---- Plots ----

  output$filtered_umap_plot <- renderPlot({
    req(multiple_datasets_object(), input$filtered_umap_col_merge)
    seurat_obj <- multiple_datasets_object()
    reduc_input <- if (!is.null(input$filtered_umap_reduction_merge) && nzchar(input$filtered_umap_reduction_merge)) {
      input$filtered_umap_reduction_merge
    } else if (!is.null(input$umap_reduction_merge) && nzchar(input$umap_reduction_merge)) {
      input$umap_reduction_merge
    } else {
      NULL
    }
    umap_reduc <- if (!is.null(reduc_input) && reduc_input %in% names(seurat_obj@reductions)) {
      reduc_input
    } else {
      get_umap_reduction(seurat_obj)
    }
    if (is.null(umap_reduc)) return(NULL)
    show_labels <- input$show_labels_merge
    label_size  <- ifelse(input$bold_labels_merge, 8, 5)
    cluster_col <- if (nzchar(input$filtered_umap_col_merge) &&
                       input$filtered_umap_col_merge %in% colnames(seurat_obj@meta.data)) {
      input$filtered_umap_col_merge
    } else {
      get_cluster_column(seurat_obj)
    }
    # Apply dataset filter only if the column exists and a filter is selected
    has_dataset_col <- "dataset" %in% colnames(seurat_obj@meta.data)
    seurat_to_plot <- if (has_dataset_col &&
                          !is.null(input$dataset_filter) &&
                          !"All Datasets" %in% input$dataset_filter) {
      valid_datasets <- input$dataset_filter %in% unique(seurat_obj@meta.data$dataset)
      subset(seurat_obj, subset = dataset %in% input$dataset_filter[valid_datasets])
    } else {
      seurat_obj
    }
    plot <- DimPlot(seurat_to_plot, reduction = umap_reduc,
                    group.by = cluster_col, label = show_labels,
                    label.size = label_size) +
      NoAxes() + NoLegend() + ggtitle(NULL)
    cluster_colors <- get_cluster_colors_merge(seurat_obj)
    if (!is.null(cluster_colors)) {
      plot <- plot + scale_color_manual(values = cluster_colors)
    }
    if (input$dark_mode_filtered_umap) {
      plot <- plot +
        theme(panel.background = element_rect(fill = "black"),
              plot.background  = element_rect(fill = "black"),
              text             = element_text(color = "white"),
              panel.grid       = element_blank())
      if (show_labels) {
        for (i in seq_along(plot$layers)) {
          if (inherits(plot$layers[[i]]$geom, c("GeomText", "GeomTextRepel"))) {
            plot$layers[[i]]$aes_params$colour <- "white"
          }
        }
      }
    }
    filtered_umap_plot(plot)
    plot
  })

  # ---- Downloads ----

  output$download_filtered_umap_plot <- createDownloadHandler(
    reactive_data = filtered_umap_plot,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedUMAP") }),
    data_name = reactive({ paste0("filtered_UMAP", if(!is.null(input$dataset_filter)) paste0("_", paste(input$dataset_filter, collapse = "_")) else "") }),
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$filtered_umap_format }),
      width = reactive({ if(input$filtered_umap_format == "pdf") 11 else 10 }),
      height = reactive({ if(input$filtered_umap_format == "pdf") 8 else 6 }),
      dpi = reactive({ input$filtered_umap_plot_dpi })
    )
  )

  ############################# Cluster vs All Comparison #############################

  # ---- Observers and computations ----

  # Cluster vs All
  observe(run_safely("Interface update", persistent = FALSE, {
    req(multiple_datasets_object())
    cluster_choices <- active_clusters_merge()
    updateSelectInput(session, "selected_cluster", choices = cluster_choices, selected = cluster_choices[1])
  }))

  # Reactive function for markers
  # Reactive function for markers
  observeEvent(input$calculate_DE, {
    tryCatch({
      showModal(modalDialog(
        title = "Please wait",
        "Calculating differentially expressed genes...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      showNotification("Calculating differentially expressed genes...", type = "message")
      
      subset_seurat <- multiple_datasets_object()
      cluster_col   <- active_cluster_column_merge()
      # Split Seurat v5 data layers must be joined before FindMarkers (rule 1)
      de_assay <- input$assay_de_merge %||% DefaultAssay(subset_seurat)
      if (inherits(subset_seurat[[de_assay]], "Assay5") && length(SeuratObject::Layers(subset_seurat[[de_assay]], search = "data")) > 1) subset_seurat <- JoinLayers(subset_seurat, assay = de_assay)
      
      markers <- FindMarkers(subset_seurat,
                             ident.1 = input$selected_cluster,
                             group.by = cluster_col,
                             min.pct = input$min_pct_merge,
                             logfc.threshold = input$logfc_threshold_merge,
                             layer = "data",
                             assay = input$assay_de_merge,
                             return.thresh = 0.05)
      
      if (nrow(markers) == 0) {
        showNotification("No differentially expressed genes found for the selected cluster.", type = "message")
        removeModal()
        return()
      }
      
      # Store raw numeric results for Venn diagrams BEFORE any formatting (rule #7)
      markers$gene <- rownames(markers)
      table_name  <- paste0("Cluster_", input$selected_cluster, "_vs_All_", format(Sys.time(), "%H%M%S"))
      description <- paste0("Cluster ", input$selected_cluster, " vs All Others")
      parameters  <- list(min_pct = input$min_pct_merge, logfc_threshold = input$logfc_threshold_merge, group_by = cluster_col)
      
      store_de_table(
        table_data  = markers,
        table_name  = table_name,
        description = description,
        type        = "single_cluster",
        parameters  = parameters
      )
      
      # Store raw markers — formatting and filtering happen at render time
      diff_genes_compare(markers)
      shinyjs::enable("download_markers_single_cluster_merge")
      
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Error calculating differentially expressed genes: ", e$message), type = "error")
    })
  })

  # ---- Tables ----

  output$DE_genes_table <- renderDataTable({
    req(diff_genes_compare())
    tryCatch({
      df <- diff_genes_compare()
      threshold <- input$pval_adj_filter_single_cluster_merge
      if (!is.null(threshold) && !is.na(threshold)) {
        df <- df[df$p_val_adj <= threshold, , drop = FALSE]
      }
      cleaned <- clean_gene_names_for_html(df$gene)
      df$gene_link <- paste0('<a href="#" class="gene-name" data-gene="', cleaned, '">', cleaned, '</a>')
      # p-values stay numeric, formatted client-side by build_de_datatable
      df <- df[, c("gene_link", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
      colnames(df) <- c("Gene", "Log2FC", "Pct.1", "Pct.2", "P-value", "Adj. P-value")
      build_de_datatable(df, order_col = "Log2FC")
    }, error = function(e) {
      showNotification(paste0("Error displaying gene table: ", e$message), type = "error")
    })
  })

  # ---- Downloads ----

  output$download_markers_single_cluster_merge <- createDownloadHandler(
    reactive_data = reactive({
      df <- diff_genes_compare()
      req(df)
      threshold <- input$pval_adj_filter_single_cluster_merge
      if (!is.null(threshold) && !is.na(threshold)) {
        df <- df[df$p_val_adj <= threshold, , drop = FALSE]
      }
      df <- df[order(df$p_val_adj), ]
      df$p_val     <- format_pvalue_robust(df$p_val)
      df$p_val_adj <- format_pvalue_robust(df$p_val_adj)
      colnames(df)[colnames(df) == "gene"]       <- "Gene"
      colnames(df)[colnames(df) == "avg_log2FC"] <- "Log2FC"
      colnames(df)[colnames(df) == "pct.1"]      <- "Pct.1"
      colnames(df)[colnames(df) == "pct.2"]      <- "Pct.2"
      colnames(df)[colnames(df) == "p_val"]      <- "P-value"
      colnames(df)[colnames(df) == "p_val_adj"]  <- "Adj. P-value"
      return(df)
    }),
    object_name_reactive = reactive({
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets")
    }),
    data_name = reactive({ paste0("cluster_", input$selected_cluster, "_vs_all_others") }),
    download_type = "csv"
  )

  ############################# Group vs Group Comparison #############################

  # ---- Observers and computations ----

  # Observer to update second group choices
  observeEvent(input$cluster1_compare, run_safely("cluster1 compare", {
    req(multiple_datasets_object())
    remaining_clusters <- setdiff(active_clusters_merge(), input$cluster1_compare)
    updateSelectInput(session, "cluster2_compare", choices = remaining_clusters)
  }))

  # Cluster comparison
  observeEvent(input$compare_clusters_button, {
    tryCatch({
      showModal(modalDialog(
        title = "Please wait",
        "Finding differentially expressed genes...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      req(input$cluster1_compare, input$cluster2_compare, multiple_datasets_object())
      
      if (any(input$cluster1_compare %in% input$cluster2_compare)) {
        showNotification("Groups must not have overlapping clusters!", type = "error")
        removeModal()
        return()
      }
      
      seurat_obj <- multiple_datasets_object()
      # Groups built from the column the cluster selectors are filled from (renamed / annotated clusters), not from Idents()
      cluster_col <- active_cluster_column_merge()
      cluster_vec <- if (!is.null(cluster_col) && cluster_col %in% colnames(seurat_obj@meta.data)) as.character(seurat_obj@meta.data[[cluster_col]]) else as.character(Idents(seurat_obj))
      new_idents <- cluster_vec
      new_idents[cluster_vec %in% input$cluster1_compare] <- "group1"
      new_idents[cluster_vec %in% input$cluster2_compare] <- "group2"
      if (sum(new_idents == "group1") == 0 || sum(new_idents == "group2") == 0) {
        showNotification("No cells found for one of the groups in the active cluster column", type = "error")
        removeModal()
        return()
      }
      Idents(seurat_obj) <- factor(stats::setNames(new_idents, colnames(seurat_obj)))
      # Split Seurat v5 data layers must be joined before FindMarkers (rule 1)
      de_assay <- input$assay_de_merge %||% DefaultAssay(seurat_obj)
      if (inherits(seurat_obj[[de_assay]], "Assay5") && length(SeuratObject::Layers(seurat_obj[[de_assay]], search = "data")) > 1) seurat_obj <- JoinLayers(seurat_obj, assay = de_assay)      
      temp_res <- FindMarkers(seurat_obj,
                              ident.1 = "group1",
                              ident.2 = "group2",
                              min.pct = input$min_pct_compare_merge,
                              logfc.threshold = input$logfc_threshold_compare_merge,
                              assay = input$assay_de_merge,
                              return.thresh = 0.05)
      
      # Store raw numeric results for Venn diagrams BEFORE any formatting (rule #7)
      temp_res$gene <- rownames(temp_res)
      temp_res$comparison <- sprintf("Group1(%s) vs Group2(%s)",
                                     paste(input$cluster1_compare, collapse = ","),
                                     paste(input$cluster2_compare, collapse = ","))
      
      group1_text <- paste(input$cluster1_compare, collapse = "_")
      group2_text <- paste(input$cluster2_compare, collapse = "_")
      table_name  <- paste0("Clusters_", group1_text, "_vs_", group2_text, "_", format(Sys.time(), "%H%M%S"))
      description <- paste0("Clusters [", group1_text, "] vs [", group2_text, "]")
      parameters  <- list(
        min_pct = input$min_pct_compare_merge,
        logfc.threshold = input$logfc_threshold_compare_merge,
        group1 = input$cluster1_compare,
        group2 = input$cluster2_compare
      )
      
      store_de_table(
        table_data  = temp_res,
        table_name  = table_name,
        description = description,
        type        = "cluster_group",
        parameters  = parameters
      )
      
      # Store raw markers — formatting and filtering happen at render time
      diff_genes_compare_cluster(temp_res)
      shinyjs::enable("download_markers_multiple_clusters_merge")
      
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste("Cluster comparison error:", e$message), type = "error")
    })
  })

  # ---- UI and text outputs ----

  # Compare Two Clusters
  output$cluster1_compare_ui <- renderUI({
    req(multiple_datasets_object())
    selectInput("cluster1_compare",
                label    = "Select clusters for group 1",
                choices  = active_clusters_merge(),
                selected = NULL,
                multiple = TRUE)
  })

  output$cluster2_compare_ui <- renderUI({
    req(multiple_datasets_object())
    remaining_clusters <- setdiff(active_clusters_merge(), input$cluster1_compare)
    selectInput("cluster2_compare",
                label    = "Select clusters for group 2",
                choices  = remaining_clusters,
                selected = NULL,
                multiple = TRUE)
  })

  # ---- Tables ----

  output$diff_genes_table_compare <- renderDataTable({
    req(diff_genes_compare_cluster())
    tryCatch({
      df <- diff_genes_compare_cluster()
      threshold <- input$pval_adj_filter_multiple_clusters_merge
      if (!is.null(threshold) && !is.na(threshold)) {
        df <- df[df$p_val_adj <= threshold, , drop = FALSE]
      }
      cleaned <- clean_gene_names_for_html(df$gene)
      df$gene_link <- paste0('<a href="#" class="gene-name" data-gene="', cleaned, '">', cleaned, '</a>')
      # p-values stay numeric, formatted client-side by build_de_datatable
      df <- df[, c("gene_link", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj", "comparison")]
      colnames(df) <- c("Gene", "Log2FC", "Pct.1", "Pct.2", "P-value", "Adj. P-value", "Comparison")
      build_de_datatable(df, order_col = "Log2FC")
    }, error = function(e) {
      showNotification(paste0("Error displaying gene table: ", e$message), type = "error")
    })
  })

  # ---- Downloads ----

  output$download_markers_multiple_clusters_merge <- createDownloadHandler(
    reactive_data = reactive({
      df <- diff_genes_compare_cluster()
      req(df)
      threshold <- input$pval_adj_filter_multiple_clusters_merge
      if (!is.null(threshold) && !is.na(threshold)) {
        df <- df[df$p_val_adj <= threshold, , drop = FALSE]
      }
      df <- df[order(df$p_val_adj), ]
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
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets")
    }),
    data_name = reactive({
      paste0("diff_genes_comparison_",
             paste(input$cluster1_compare, collapse = "-"),
             "_VS_",
             paste(input$cluster2_compare, collapse = "-"))
    }),
    download_type = "csv"
  )

  ############################# Group Comparison (Metadata-based) #############################

  # ---- Observers and computations ----

  # Observer to update group 2 choices when group 1 changes
  observeEvent(input$dataset1_compare, run_safely("dataset1 compare", {
    req(multiple_datasets_object(), input$metadata_column_compare)
    
    metadata_values <- unique(multiple_datasets_object()@meta.data[[input$metadata_column_compare]])
    metadata_values <- sort(metadata_values[!is.na(metadata_values)])
    remaining_values <- setdiff(metadata_values, input$dataset1_compare)
    
    updateSelectInput(session, "dataset2_compare", choices = remaining_values)
  }))

  # Samples x groups design of the pending pseudobulk comparison, shown before the analysis is run
  pseudobulk_design_preview <- reactive({
    req(identical(input$de_method_compare, "pseudobulk"), multiple_datasets_object(), input$sample_column_compare,
        input$metadata_column_compare, input$cluster_column_compare, input$dataset1_compare, input$dataset2_compare)
    seurat_obj <- multiple_datasets_object()
    req(all(c(input$sample_column_compare, input$metadata_column_compare) %in% colnames(seurat_obj@meta.data)))
    clusters <- if (isTRUE(input$all_clusters)) NULL else input$cluster_compare
    req(isTRUE(input$all_clusters) || length(clusters) > 0)
    min_cells <- if (is.numeric(input$pb_min_cells_compare) && !is.na(input$pb_min_cells_compare)) input$pb_min_cells_compare else 10
    condition <- compare_cell_condition(seurat_obj, input$metadata_column_compare, input$cluster_column_compare, input$dataset1_compare, input$dataset2_compare, clusters)
    summarize_pseudobulk_design(seurat_obj@meta.data[[input$sample_column_compare]], condition, min_cells)
  })

  observeEvent(input$compare_datasets_button, {
    run_safely("Group comparison", {
      req(input$metadata_column_compare, input$dataset1_compare, input$dataset2_compare, input$cluster_column_compare, multiple_datasets_object())
      if (any(input$dataset1_compare %in% input$dataset2_compare)) {
        showNotification("Groups must not have overlapping values!", type = "error")
        return()
      }
      seurat_obj   <- multiple_datasets_object()
      metadata_col <- input$metadata_column_compare
      cluster_col  <- input$cluster_column_compare
      use_pseudobulk <- identical(input$de_method_compare, "pseudobulk")
      message(paste("Using metadata column:", metadata_col))
      message(paste("Using cluster column:", cluster_col))
      if (!metadata_col %in% colnames(seurat_obj@meta.data)) {
        showNotification(paste("Metadata column", metadata_col, "not found!"), type = "error")
        return()
      }
      if (!cluster_col %in% colnames(seurat_obj@meta.data) && cluster_col != "seurat_clusters") {
        showNotification(paste("Cluster column", cluster_col, "not found!"), type = "error")
        return()
      }
      cluster_values <- if (cluster_col == "seurat_clusters") as.character(Idents(object = seurat_obj)) else as.character(seurat_obj@meta.data[[cluster_col]])
      clusters_to_analyze <- if (input$all_clusters) unique(cluster_values[!is.na(cluster_values) & cluster_values != ""]) else input$cluster_compare
      if (!input$all_clusters && (is.null(clusters_to_analyze) || length(clusters_to_analyze) == 0)) {
        showNotification("Please select at least one cluster or check 'Compare all clusters'", type = "error")
        return()
      }
      if (use_pseudobulk && (is.null(input$sample_column_compare) || !nzchar(input$sample_column_compare))) {
        # Checked before the modal opens: a silent req() inside run_safely would leave the modal on screen
        stop("Select the sample column (one value = one animal / biological replicate) for the pseudobulk test")
      }
      showModal(modalDialog(title = "Please wait", if (use_pseudobulk) "Pseudobulk DESeq2 on samples..." else "Comparing groups and clusters...", easyClose = FALSE, footer = NULL))
      if (use_pseudobulk) {
        # Sample-level test: raw RNA counts summed per sample, DESeq2 dispersion estimated between samples of the same group
        min_cells <- if (is.numeric(input$pb_min_cells_compare) && !is.na(input$pb_min_cells_compare)) input$pb_min_cells_compare else 10
        condition <- compare_cell_condition(seurat_obj, metadata_col, cluster_col, input$dataset1_compare, input$dataset2_compare, clusters_to_analyze)
        pb_res <- run_pseudobulk_de(seurat_obj, condition, input$sample_column_compare, assay = "RNA", min_cells = min_cells, min_samples = 2)
        design <- attr(pb_res, "design")
        n_samples <- c(sum(design$used & design$condition == "Group1"), sum(design$used & design$condition == "Group2"))
        message(sprintf("Pseudobulk DESeq2: %d vs %d samples, %d genes tested", n_samples[1], n_samples[2], nrow(pb_res)))
        # BH correction was done on all tested genes; the user thresholds only filter the output
        temp_res <- filter_pseudobulk_results(pb_res, logfc_threshold = input$logfc_threshold_datasets, min_pct = input$min_pct_compare_dataset_merge, max_p = 0.05)
      } else {
        # Cell-level test on the log-normalized RNA assay, never on the integrated (CCA) or SCT assay
        seurat_obj <- prepare_rna_assay_for_de(seurat_obj)$obj
        seurat_obj$temp_grouping <- paste(seurat_obj@meta.data[[metadata_col]], cluster_values, sep = "_")
        group1_ids <- as.vector(outer(input$dataset1_compare, clusters_to_analyze, paste, sep = "_"))
        group2_ids <- as.vector(outer(input$dataset2_compare, clusters_to_analyze, paste, sep = "_"))
        available_groups <- unique(seurat_obj$temp_grouping)
        missing_group1 <- setdiff(group1_ids, available_groups)
        missing_group2 <- setdiff(group2_ids, available_groups)
        if (length(missing_group1) > 0) showNotification(paste("Group 1 combinations not found:", paste(missing_group1, collapse = ", ")), type = "warning", duration = 10)
        if (length(missing_group2) > 0) showNotification(paste("Group 2 combinations not found:", paste(missing_group2, collapse = ", ")), type = "warning", duration = 10)
        group1_ids <- intersect(group1_ids, available_groups)
        group2_ids <- intersect(group2_ids, available_groups)
        if (length(group1_ids) == 0 || length(group2_ids) == 0) {
          removeModal()
          showNotification("No valid group combinations found for comparison!", type = "error")
          return()
        }
        message(paste("Comparing", length(group1_ids), "vs", length(group2_ids), "group-cluster combinations"))
        temp_res <- FindMarkers(object = seurat_obj, ident.1 = group1_ids, ident.2 = group2_ids, group.by = "temp_grouping", assay = "RNA",
                                min.pct = input$min_pct_compare_dataset_merge, logfc.threshold = input$logfc_threshold_datasets, return.thresh = 0.05)
      }
      if (nrow(temp_res) == 0) {
        removeModal()
        showNotification("No differentially expressed genes found with these parameters!", type = "error")
        return()
      }
      # Store raw numeric results for Venn diagrams BEFORE any formatting (rule #7)
      method_tag   <- if (use_pseudobulk) "PB" else "WX"
      group1_text  <- paste(input$dataset1_compare, collapse = "_")
      group2_text  <- paste(input$dataset2_compare, collapse = "_")
      cluster_text <- if (input$all_clusters) "AllClusters" else paste(clusters_to_analyze, collapse = "_")
      table_name   <- paste0("Metadata_", method_tag, "_", group1_text, "_vs_", group2_text, "_", cluster_text, "_", format(Sys.time(), "%H%M%S"))
      description  <- paste0(metadata_col, ": [", group1_text, "] vs [", group2_text, "] | ", cluster_col, ": ", cluster_text,
                             if (use_pseudobulk) sprintf(" | pseudobulk DESeq2 (%d vs %d samples)", n_samples[1], n_samples[2]) else " | Wilcoxon")
      parameters   <- list(metadata_col = metadata_col, cluster_col = cluster_col, group1 = input$dataset1_compare, group2 = input$dataset2_compare,
                           clusters = clusters_to_analyze, min_pct = input$min_pct_compare_dataset_merge, logfc_threshold = input$logfc_threshold_datasets,
                           method = if (use_pseudobulk) "pseudobulk_DESeq2" else "wilcoxon",
                           sample_col = if (use_pseudobulk) input$sample_column_compare else NA, n_samples = if (use_pseudobulk) n_samples else NA)
      store_de_table(table_data = temp_res, table_name = table_name, description = description, type = "metadata_group", parameters = parameters)
      # Add comparison label then store raw — formatting and filtering happen at render time
      cluster_label <- paste(if (input$all_clusters) "all" else clusters_to_analyze, collapse = ", ")
      temp_res$comparison <- sprintf("Group1(%s: %s | %s: %s) vs Group2(%s: %s | %s: %s) [%s]",
                                     metadata_col, paste(input$dataset1_compare, collapse = ", "), cluster_col, cluster_label,
                                     metadata_col, paste(input$dataset2_compare, collapse = ", "), cluster_col, cluster_label,
                                     if (use_pseudobulk) sprintf("pseudobulk DESeq2, %d vs %d samples", n_samples[1], n_samples[2]) else "Wilcoxon, cell level")
      diff_genes_compare_datasets(temp_res)
      shinyjs::enable("download_diff_dataset_cluster")
      removeModal()
      showNotification(paste0("Found ", nrow(temp_res), " differentially expressed genes"), type = "message", duration = 5)
    })
  })

  # ---- UI and text outputs ----

  # UI to select which metadata column to use for grouping
  output$metadata_column_compare_ui <- renderUI({
    req(multiple_datasets_object())
    
    # Get all metadata columns
    meta_cols <- colnames(multiple_datasets_object()@meta.data)
    
    # Exclude technical/unwanted columns that don't make sense for grouping
    exclude_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", 
                      "S.Score", "G2M.Score", "Phase", "ClusterIdents",
                      "RNA_snn_res.0.5", "RNA_snn_res.0.8", "RNA_snn_res.1",
                      "seurat_clusters", "ident", "orig.ident")
    
    # Keep only categorical/factor columns or character columns
    available_cols <- setdiff(meta_cols, exclude_cols)
    
    # Default to "dataset" if it exists, otherwise first available column
    default_col <- if("dataset" %in% available_cols) "dataset" else available_cols[1]
    
    selectInput("metadata_column_compare",
                label = "Grouping variable (metadata column):",
                choices = available_cols,
                selected = default_col)
  })

  # UI to select values for group 1 based on selected metadata column
  output$dataset1_compare_ui <- renderUI({
    req(multiple_datasets_object(), input$metadata_column_compare)
    
    # Get unique values from selected metadata column
    metadata_values <- unique(multiple_datasets_object()@meta.data[[input$metadata_column_compare]])
    metadata_values <- sort(metadata_values[!is.na(metadata_values)])  # Remove NA and sort
    
    selectInput("dataset1_compare",
                label = "Select values for Group 1:",
                choices = metadata_values,
                selected = NULL,
                multiple = TRUE)
  })

  # UI to select values for group 2
  output$dataset2_compare_ui <- renderUI({
    req(multiple_datasets_object(), input$metadata_column_compare, input$dataset1_compare)
    
    # Get remaining values not in group 1
    metadata_values <- unique(multiple_datasets_object()@meta.data[[input$metadata_column_compare]])
    metadata_values <- sort(metadata_values[!is.na(metadata_values)])
    remaining_values <- setdiff(metadata_values, input$dataset1_compare)
    
    selectInput("dataset2_compare",
                label = "Select values for Group 2:",
                choices = remaining_values,
                selected = NULL,
                multiple = TRUE)
  })

  # NEW: Selector for which clustering column to use
  output$cluster_column_selector_ui <- renderUI({
    req(multiple_datasets_object())
    
    seurat_obj <- multiple_datasets_object()
    meta_cols <- colnames(seurat_obj@meta.data)
    
    # Find potential cluster columns
    # Look for columns that are factors or have reasonable number of unique values
    potential_cluster_cols <- meta_cols[sapply(meta_cols, function(col) {
      col_data <- seurat_obj@meta.data[[col]]
      is.factor(col_data) || (length(unique(col_data)) < 100 && length(unique(col_data)) > 1)
    })]
    
    # Prioritize common cluster column names
    priority_cols <- c("ClusterIdents", "seurat_clusters", "annotated_clusters", 
                       "cell_type", "celltype", "clusters")
    priority_cols <- intersect(priority_cols, potential_cluster_cols)
    other_cols <- setdiff(potential_cluster_cols, priority_cols)
    
    cluster_col_choices <- c(priority_cols, other_cols)
    
    # Default selection
    default_col <- if ("ClusterIdents" %in% cluster_col_choices) {
      "ClusterIdents"
    } else if ("seurat_clusters" %in% cluster_col_choices) {
      "seurat_clusters"
    } else {
      cluster_col_choices[1]
    }
    
    selectInput("cluster_column_compare",
                label = "Select clustering column to use:",
                choices = cluster_col_choices,
                selected = default_col)
  })

  # UI to select clusters to analyze - CORRECTED VERSION
  output$cluster_compare_ui <- renderUI({
    message("=== cluster_compare_ui triggered ===")
    
    req(multiple_datasets_object(), input$cluster_column_compare)
    
    if (input$all_clusters) {
      message("cluster_compare_ui: all_clusters is TRUE, returning NULL")
      return(NULL)
    }
    
    seurat_obj <- multiple_datasets_object()
    cluster_col <- input$cluster_column_compare
    
    message(paste("cluster_compare_ui: Using cluster column:", cluster_col))
    
    # Get cluster values using the SAME logic as gene expression plotting
    cluster_values <- tryCatch({
      if (cluster_col == "seurat_clusters") {
        # For seurat_clusters, use Idents
        cluster_idents <- Idents(object = seurat_obj)
        as.character(cluster_idents)
      } else if (cluster_col %in% colnames(seurat_obj@meta.data)) {
        # For other metadata columns
        col_data <- seurat_obj@meta.data[[cluster_col]]
        # Handle factors properly
        if (is.factor(col_data)) {
          as.character(col_data)
        } else {
          as.character(col_data)
        }
      } else {
        message(paste("ERROR: Cluster column", cluster_col, "not found!"))
        return(tags$p(paste("Column", cluster_col, "not found"), style="color: red;"))
      }
    }, error = function(e) {
      message(paste("ERROR getting cluster values:", e$message))
      return(NULL)
    })
    
    if (is.null(cluster_values)) {
      return(tags$p("Error loading clusters", style="color: red;"))
    }
    
    message(paste("cluster_compare_ui: Raw cluster values (first 10):", 
                  paste(head(cluster_values, 10), collapse=", ")))
    
    # Get unique clusters and clean
    cluster_choices <- unique(cluster_values)
    cluster_choices <- cluster_choices[!is.na(cluster_choices) & cluster_choices != ""]
    
    message(paste("cluster_compare_ui: Unique clusters:", 
                  paste(cluster_choices, collapse=", ")))
    
    # Sort: numeric first, then alphabetic
    numeric_clusters <- cluster_choices[!is.na(suppressWarnings(as.numeric(cluster_choices)))]
    text_clusters <- cluster_choices[is.na(suppressWarnings(as.numeric(cluster_choices)))]
    
    if (length(numeric_clusters) > 0) {
      numeric_clusters <- as.character(sort(as.numeric(numeric_clusters)))
    }
    if (length(text_clusters) > 0) {
      text_clusters <- sort(text_clusters)
    }
    
    cluster_choices <- c(numeric_clusters, text_clusters)
    
    if (length(cluster_choices) == 0) {
      message("ERROR: No clusters found!")
      return(tags$p("No clusters found", style="color: red;"))
    }
    
    selectizeInput("cluster_compare",
                   label = paste0("Select clusters to analyze (from '", cluster_col, "'):"),
                   choices = cluster_choices,
                   selected = NULL,
                   multiple = TRUE,
                   options = list(
                     plugins = list('remove_button'),
                     placeholder = 'Select multiple clusters'
                   ))
  })

  # Sample column for pseudobulk: identifies the biological replicate (one animal = one sample)
  output$sample_column_compare_ui <- renderUI({
    req(multiple_datasets_object())
    meta <- multiple_datasets_object()@meta.data
    categorical <- names(meta)[vapply(meta, function(x) is.character(x) || is.factor(x), logical(1))]
    default_col <- if ("dataset" %in% categorical) "dataset" else if ("orig.ident" %in% categorical) "orig.ident" else categorical[1]
    selectInput("sample_column_compare", label = "Sample column (one value = one animal / biological replicate):", choices = categorical, selected = default_col)
  })

  output$pseudobulk_design_status <- renderUI({
    design <- pseudobulk_design_preview()
    problem <- check_pseudobulk_design(design, min_samples = 2)
    if (!is.null(problem)) return(tags$p(icon("exclamation-triangle"), " ", problem, style = "color: #c62828; margin: 5px 0;"))
    n1 <- sum(design$used & design$condition == "Group1"); n2 <- sum(design$used & design$condition == "Group2")
    note <- if (min(n1, n2) < 3) " At least 3 samples per group are recommended: with 2 the variance estimate relies mostly on the other genes." else ""
    tags$p(icon("check-circle"), sprintf(" Design OK: %d vs %d samples.%s", n1, n2, note), style = "color: #2e7d32; margin: 5px 0;")
  })

  # ---- Tables ----

  output$pseudobulk_design_compare <- renderTable({
    design <- pseudobulk_design_preview()
    req(nrow(design) > 0)
    group_label <- c(Group1 = paste("Group 1:", paste(input$dataset1_compare, collapse = ", ")), Group2 = paste("Group 2:", paste(input$dataset2_compare, collapse = ", ")))
    data.frame(Sample = design$sample, Group = unname(group_label[design$condition]), Cells = design$n_cells, Used = ifelse(design$used, "yes", "no (too few cells)"))
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$diff_dataset_cluster <- renderDataTable({
    req(diff_genes_compare_datasets())
    df <- diff_genes_compare_datasets()
    threshold <- input$pval_adj_filter_datasets_merge
    if (!is.null(threshold) && !is.na(threshold)) {
      df <- df[df$p_val_adj <= threshold, , drop = FALSE]
    }
    cleaned <- clean_gene_names_for_html(rownames(df))
    df$gene_link <- paste0('<a href="#" class="gene-name" data-gene="', cleaned, '">', cleaned, '</a>')
    # p-values stay numeric, formatted client-side by build_de_datatable
    display_cols <- c("gene_link", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj", "comparison")
    df <- df[, display_cols]
    colnames(df) <- c("Gene", "Log2FC", "Pct.1", "Pct.2", "P-value", "Adj. P-value", "Comparison")
    build_de_datatable(df, order_col = "Log2FC")
  })

  # ---- Downloads ----

  output$download_diff_dataset_cluster <- downloadHandler(
    filename = function() {
      dataset1_text <- paste(input$dataset1_compare, collapse = "-")
      dataset2_text <- paste(input$dataset2_compare, collapse = "-")
      cluster_text  <- if (input$all_clusters) "AllClusters" else paste(input$cluster_compare, collapse = "-")
      paste0("diff-datasets-", dataset1_text, "-vs-", dataset2_text, "-", cluster_text, "-", Sys.Date(), ".csv")
    },
    content = function(file) {
      tryCatch({
        df <- diff_genes_compare_datasets()
        req(!is.null(df))
        threshold <- input$pval_adj_filter_datasets_merge
        if (!is.null(threshold) && !is.na(threshold)) {
          df <- df[df$p_val_adj <= threshold, , drop = FALSE]
        }
        df <- df[order(df$p_val_adj), ]
        df$p_val     <- format_pvalue_robust(df$p_val)
        df$p_val_adj <- format_pvalue_robust(df$p_val_adj)
        colnames(df)[colnames(df) == "avg_log2FC"] <- "Log2FC"
        colnames(df)[colnames(df) == "pct.1"]      <- "Pct.1"
        colnames(df)[colnames(df) == "pct.2"]      <- "Pct.2"
        colnames(df)[colnames(df) == "p_val"]      <- "P-value"
        colnames(df)[colnames(df) == "p_val_adj"]  <- "Adj. P-value"
        colnames(df)[colnames(df) == "comparison"] <- "Comparison"
        write.csv(df, file)
      }, error = function(e) {
        showNotification(paste0("Error downloading gene comparison table: ", e$message), type = "error")
      })
    },
    contentType = "text/csv"
  )

  ############################## Cluster Composition Multiple ##############################

  ############################# Metadata Column Selectors for Composition Analysis #############################

  # ---- Observers and computations ----

  # Observer for generating cluster composition table for multiple datasets
  observeEvent(input$generate_cluster_table_multiple, {
    tryCatch({
      req(multiple_datasets_object(), input$metadata_column_composition)
      
      showModal(modalDialog(
        title = "Generating Cluster Composition Table",
        "Analyzing cluster composition across groups...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      seurat_obj <- multiple_datasets_object()
      metadata_col <- input$metadata_column_composition
      
      # Validate metadata column exists
      if (!metadata_col %in% colnames(seurat_obj@meta.data)) {
        showNotification(paste("Metadata column", metadata_col, "not found!"), type = "error")
        removeModal()
        return()
      }
      
      # Use modular function with custom metadata column
      cluster_composition <- create_cluster_composition_table(
        seurat_obj, 
        is_integrated = TRUE,
        metadata_column = metadata_col
      )
      
      # Store in reactive value
      cluster_composition_multiple(cluster_composition)
      
      removeModal()
      showNotification(
        paste("Cluster composition table generated for:", metadata_col), 
        type = "message"
      )
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error generating table:", e$message), type = "error", duration = 10)
    })
  })

  # ---- UI and text outputs ----

  # Metadata column selector for cluster composition table
  output$metadata_column_composition_ui <- renderUI({
    req(multiple_datasets_object())
    
    meta_cols <- colnames(multiple_datasets_object()@meta.data)
    
    # Exclude technical columns
    excluded_patterns <- c(
      "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo",
      "S.Score", "G2M.Score", "Phase", "ClusterIdents",
      "^RNA_snn_res", "^integrated", "seurat_clusters", "ident"
    )
    
    pattern <- paste(excluded_patterns, collapse = "|")
    available_cols <- meta_cols[!grepl(pattern, meta_cols)]
    
    # Default to "dataset" if it exists
    default_col <- if("dataset" %in% available_cols) "dataset" else available_cols[1]
    
    selectInput("metadata_column_composition",
                label = "Group by:",
                choices = available_cols,
                selected = default_col)
  })

  # Metadata column selector for pie chart
  output$metadata_column_pie_ui <- renderUI({
    req(multiple_datasets_object())
    
    meta_cols <- colnames(multiple_datasets_object()@meta.data)
    
    # Exclude technical columns
    excluded_patterns <- c(
      "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo",
      "S.Score", "G2M.Score", "Phase", "ClusterIdents",
      "^RNA_snn_res", "^integrated", "seurat_clusters", "ident"
    )
    
    pattern <- paste(excluded_patterns, collapse = "|")
    available_cols <- meta_cols[!grepl(pattern, meta_cols)]
    
    # Default to "dataset" if it exists
    default_col <- if("dataset" %in% available_cols) "dataset" else available_cols[1]
    
    selectInput("metadata_column_pie",
                label = "Group by:",
                choices = available_cols,
                selected = default_col)
  })

  # ---- Tables ----

  # Render cluster composition table for multiple datasets using existing function
  output$cluster_table_multiple <- renderDT({
    req(cluster_composition_multiple())
    render_cluster_composition_table(cluster_composition_multiple(), is_integrated = TRUE)
  })

  # ---- Downloads ----

  # Download handler for multiple datasets cluster composition using modular functions
  output$download_cluster_composition_multiple <- createDownloadHandler(
    reactive_data = reactive({ data <- cluster_composition_multiple(); data$Size_Bar <- NULL; return(data) }),
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets") }),
    data_name = "cluster_composition",
    download_type = "csv"
  )

  ###############Co-expression calculation#############

  # ---- Reactive state ----
  gene_coexpression_data_multiple <- reactiveVal(NULL)

  # Estimate each gene's threshold from its level in background cells and fill the threshold fields
  coexpr_background_estimate_multiple <- reactiveVal(NULL)

  coexpr_de_results_multiple <- reactiveVal(NULL)

  # ---- Helpers ----

  # Thresholds of the validated genes, read from the per-token fields of the selected unit (any detection when a field is empty)
  coexpr_thresholds_multi <- function(check) {
    unit <- input$coexpr_threshold_unit_multiple %||% "counts"
    default_value <- if (unit == "counts") 1 else 0
    values <- vapply(seq_along(check$tokens), function(i) {
      val <- input[[paste0("coexpr_thresh_multi_", unit, "_", i)]]
      if (is.null(val) || is.na(val)) default_value else val
    }, numeric(1))
    align_gene_values(check, values)
  }

  # ---- Observers and computations ----

  # Plot rebuilt from the stored results whenever the grouping or visualization selectors change (no new analysis needed)
  coexpression_plot_multiple <- reactive({
    res <- gene_coexpression_data_multiple()
    req(res, length(res$genes_analyzed) >= 2)
    group_by <- input$coexpr_group_by_multiple %||% "dataset"
    secondary_split <- if (identical(group_by, "dataset")) input$coexpr_secondary_split_multiple_alt %||% "none" else input$coexpr_secondary_split_multiple %||% "none"
    create_coexpression_plot(res$data, res$genes_analyzed, group_by = group_by, secondary_split = secondary_split, plot_type = input$coexpr_plot_type_multiple %||% "bar", split_label = res$split_label %||% "Dataset")
  })

  # Background clusters for the ambient-RNA threshold estimate (same identities as the DE cluster selector)
  observe(run_safely("Interface update", persistent = FALSE, {
    req(multiple_datasets_object())
    updateSelectInput(session, "coexpr_background_clusters_multiple", choices = sort(as.character(unique(Idents(multiple_datasets_object())))))
  }))

  observeEvent(input$estimate_coexpr_thresholds_multiple, run_safely("Co-expression background thresholds (multiple datasets)", {
    req(multiple_datasets_object())
    seurat_obj <- multiple_datasets_object()
    background <- input$coexpr_background_clusters_multiple
    if (length(background) == 0) stop("Select the cluster(s) that should not express these genes")
    assay <- input$coexpr_assay_multiple %||% "RNA"
    check <- validate_gene_input(input$gene_text_coexpression_multiple, seurat_obj, assay)
    if (!check$ok) stop(check$message)
    unit <- input$coexpr_threshold_unit_multiple %||% "counts"
    cells <- colnames(seurat_obj)[as.character(Idents(seurat_obj)) %in% background]
    est <- estimate_coexpression_thresholds(seurat_obj, check$genes, cells, assay_name = assay, threshold_layer = unit, quantile = 0.95)
    for (i in seq_along(check$tokens)) {
      g <- check$matched[i]
      if (!is.na(g)) updateNumericInput(session, paste0("coexpr_thresh_multi_", unit, "_", i), value = est$threshold[est$gene == g])
    }
    coexpr_background_estimate_multiple(list(table = est, clusters = background, unit = unit, n_cells = length(cells)))
  }))

  # Update assay choices when object is loaded
  observe(run_safely("Interface update", persistent = FALSE, {
    req(multiple_datasets_object())
    updateSelectInput(session, "coexpr_assay_multiple", choices = names(multiple_datasets_object()@assays), selected = DefaultAssay(multiple_datasets_object()))
  }))

  # Populate group-by column and filter groups for co-expression
  observe(run_safely("Interface update", persistent = FALSE, {
    req(multiple_datasets_object())
    meta <- multiple_datasets_object()@meta.data
    categorical_cols <- colnames(meta)[sapply(meta, function(x) is.character(x) || is.factor(x))]
    updateSelectInput(session, "coexpr_groupby_col_merge", choices = categorical_cols,
                      selected = if ("dataset" %in% categorical_cols) "dataset" else categorical_cols[1])
  }))

  observeEvent(input$coexpr_groupby_col_merge, run_safely("coexpr groupby col merge", {
    req(multiple_datasets_object(), nzchar(input$coexpr_groupby_col_merge))
    col <- input$coexpr_groupby_col_merge
    req(col %in% colnames(multiple_datasets_object()@meta.data))
    values <- sort(unique(as.character(multiple_datasets_object()@meta.data[[col]])))
    updateSelectInput(session, "coexpr_filter_groups_merge", choices = values, selected = NULL)
  }))

  observeEvent(input$analyze_coexpression_multiple, {
    req(multiple_datasets_object())
    tryCatch({
      seurat_obj <- multiple_datasets_object()
      check <- validate_gene_input(input$gene_text_coexpression_multiple, seurat_obj, input$coexpr_assay_multiple %||% "RNA")
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      thresholds <- coexpr_thresholds_multi(check)
      genes_input <- names(thresholds)
      if (length(genes_input) < 2) {
        showNotification("Co-expression needs at least 2 different genes found in the dataset.", type = "error", duration = 8)
        return()
      }
      # Split column: second level of the table and plot breakdown, and column used by the "Keep only" filter
      split_col <- input$coexpr_groupby_col_merge
      if (is.null(split_col) || !nzchar(split_col) || !split_col %in% colnames(seurat_obj@meta.data)) {
        showNotification("Choose a metadata column in 'Split by column'.", type = "error", duration = 8)
        return()
      }
      # Optional restriction to some values of the split column
      filter_groups <- input$coexpr_filter_groups_merge
      if (length(filter_groups) > 0) {
        cells_to_keep <- colnames(seurat_obj)[as.character(seurat_obj@meta.data[[split_col]]) %in% filter_groups]
        if (length(cells_to_keep) == 0) {
          showNotification("No cells found for the selected values!", type = "error")
          return()
        }
        seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
        message(paste0("Co-expression restricted to ", split_col, " [", paste(filter_groups, collapse = ", "), "] — ", ncol(seurat_obj), " cells"))
      }
      coexpr_results <- analyze_gene_coexpression(
        seurat_obj            = seurat_obj,
        genes                 = genes_input,
        assay_name            = input$coexpr_assay_multiple %||% "RNA",
        expression_thresholds = thresholds,
        is_integrated         = TRUE,
        threshold_layer       = input$coexpr_threshold_unit_multiple %||% "counts",
        split_col             = split_col
      )
      gene_coexpression_data_multiple(coexpr_results)
      showNotification(paste0("Co-expression analysis completed (split by ", split_col, ")"), type = "message")
    }, error = function(e) {
      showNotification(paste("Co-expression error:", conditionMessage(e)[1]), type = "error")
      message("Co-expression error details: ", conditionMessage(e)[1])
    })
  })

  # Populate cluster selector when object is available
  observe({
    req(multiple_datasets_object())
    tryCatch({
      clusters <- sort(as.character(unique(Idents(multiple_datasets_object()))))
      updateSelectInput(session, "coexpr_de_cluster_multiple", choices = clusters)
    }, error = function(e) NULL)
  })

  observeEvent(input$run_coexpr_de_multiple, {
    req(multiple_datasets_object())
    tryCatch({
      check <- validate_gene_input(input$gene_text_coexpression_multiple, multiple_datasets_object(), input$coexpr_assay_multiple %||% "RNA")
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      selected_clusters <- input$coexpr_de_cluster_multiple
      if (is.null(selected_clusters) || length(selected_clusters) == 0) {
        showNotification("Select at least one cluster", type = "warning")
        return()
      }
      seurat_obj <- multiple_datasets_object()
      # Clusters are read from the active identities (same source as the cluster selector): no seurat_clusters column needed
      if (!any(as.character(Idents(seurat_obj)) %in% selected_clusters)) stop("The selected clusters no longer match the current identities (renamed?): reselect them")
      # Per-gene thresholds from the dynamic inputs (same as the main analysis)
      de_thresholds <- coexpr_thresholds_multi(check)
      genes_input <- names(de_thresholds)
      showNotification("Running within-cluster DE...", type = "message", duration = 3)
      result <- run_coexpression_de(
        seurat_obj            = seurat_obj,
        genes                 = genes_input,
        clusters              = selected_clusters,
        assay_name            = input$coexpr_assay_multiple %||% "RNA",
        expression_thresholds = de_thresholds,
        threshold_layer       = input$coexpr_threshold_unit_multiple %||% "counts"
      )
      coexpr_de_results_multiple(result)
      label <- if (length(genes_input) == 1) "expressors" else "co-expressors"
      showNotification(
        sprintf("DE complete — %d %s vs %d others, %d DE genes (%s)",
                result$n_coexpressors, label, result$n_non_coexpressors,
                nrow(result$markers), result$cluster_label),
        type = "message", duration = 5
      )
    }, error = function(e) {
      showNotification(paste("DE error:", conditionMessage(e)[1]), type = "error", duration = 8)
      message("coexpr_de error: ", conditionMessage(e)[1])
    })
  })

  # ---- UI and text outputs ----

  # Per-gene threshold fields: one numericInput per token, in the order typed (field i <-> token i); typed values survive re-renders
  # never reused as a log-normalized value; typed values survive re-renders
  output$coexpr_thresholds_ui_multiple <- renderUI({
    genes <- parse_gene_tokens(input$gene_text_coexpression_multiple)
    if (length(genes) == 0) return(helpText("Enter gene names to set individual thresholds"))
    unit <- input$coexpr_threshold_unit_multiple %||% "counts"
    tagList(lapply(seq_along(genes), function(i) {
      input_id <- paste0("coexpr_thresh_multi_", unit, "_", i)
      if (unit == "counts") {
        numericInput(input_id, label = paste0(genes[i], ": min. UMIs per cell"), value = isolate(input[[input_id]]) %||% 1, min = 1, max = 100, step = 1)
      } else {
        numericInput(input_id, label = paste0(genes[i], ": normalized expression >"), value = isolate(input[[input_id]]) %||% 0, min = 0, max = 10, step = 0.1)
      }
    }))
  })

  output$coexpr_background_info_multiple <- renderUI({
    r <- coexpr_background_estimate_multiple()
    if (is.null(r)) return(NULL)
    lines <- sprintf("%s: background 95th percentile = %s -> %s (%.1f%% of background cells still positive)", r$table$gene,
                     r$table$background_quantile_value, format_threshold_label(r$table$gene, r$table$threshold, r$unit), r$table$pct_background_positive)
    tags$div(style = "font-size: 11px; color: #555; margin-top: 6px;",
             tags$p(sprintf("Background: %s (%d cells)", paste(r$clusters, collapse = ", "), r$n_cells), style = "margin: 0 0 2px 0; font-weight: bold;"),
             lapply(lines, function(l) tags$p(l, style = "margin: 0;")))
  })

  output$coexpr_de_info_multiple <- renderUI({
    r <- coexpr_de_results_multiple()
    if (is.null(r)) return(NULL)
    
    bd <- r$breakdown
    n_genes <- length(r$genes_used)
    # Breakdown text depends on number of genes
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
      tags$p(
        sprintf("Cluster(s): %s | %s | %d DE genes", r$cluster_label, thresh_txt, nrow(r$markers)),
        style = "color: #aaa; font-size: 12px; margin-bottom: 2px;"
      ),
      tags$p(
        sprintf("\u2713 %s: %d cells  |  \u2717 %s: %d cells%s",
                pos_label, r$n_coexpressors, neg_label, r$n_non_coexpressors, breakdown_txt),
        style = "color: #aaa; font-size: 12px; margin: 0;"
      )
    )
  })

  # ---- Plots ----

  output$gene_coexpression_plot_multiple <- renderPlot({
    req(coexpression_plot_multiple())
    coexpression_plot_multiple()
  }, height = 700)

  # ---- Tables ----

  # Render results using modular functions
  output$gene_coexpression_table_multiple <- renderDT({
    req(gene_coexpression_data_multiple())
    render_coexpression_table(gene_coexpression_data_multiple(), "gene_coexpression_table_multiple")
  })

  # Render summary statistics using modular approach
  output$coexpression_summary_stats_multiple <- renderTable({
    req(gene_coexpression_data_multiple())
    create_coexpression_summary_stats(gene_coexpression_data_multiple()$data, gene_coexpression_data_multiple()$genes_analyzed, split_label = gene_coexpression_data_multiple()$split_label %||% "Dataset")
  }, striped = TRUE, hover = TRUE, spacing = 'l')

  output$coexpr_de_table_multiple <- renderDT({
    req(coexpr_de_results_multiple())
    m <- coexpr_de_results_multiple()$markers
    m$avg_log2FC <- round(m$avg_log2FC, 3)
    m$pct.1      <- round(m$pct.1, 3)
    m$pct.2      <- round(m$pct.2, 3)
    # p-values stay numeric so the column filter renders a range slider; display formatting is client-side
    DT::datatable(
      m[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "p_val")],
      colnames  = c("Gene", "log2FC", "pct.coexpr", "pct.non", "p_adj", "p_val"),
      rownames  = FALSE,
      filter    = "top",
      extensions = "Buttons",
      options   = list(
        pageLength = 25, dom = "Bfrtip",
        buttons = list(
          list(extend = "csv",   filename = paste0("CoexprDE_", format(Sys.time(), "%Y%m%d")), text = "CSV"),
          list(extend = "excel", filename = paste0("CoexprDE_", format(Sys.time(), "%Y%m%d")), text = "Excel")
        )
      )
    ) |>
      DT::formatStyle("avg_log2FC", backgroundColor = DT::styleInterval(0, c("#c0392b22", "#1a572622"))) |>
      DT::formatSignif(c("p_val_adj", "p_val"), digits = 3)
  })

  # ---- Downloads ----

  # Download handlers using modular functions
  output$download_coexpression_table_multiple <- createDownloadHandler(
    reactive_data = reactive({ data <- gene_coexpression_data_multiple()$data; data$Coexpression_Visual <- NULL; return(data) }),
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets") }),
    data_name = "coexpression_analysis",
    download_type = "csv"
  )

  output$download_coexpression_plot_multiple <- createDownloadHandler(
    reactive_data = coexpression_plot_multiple,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets") }),
    data_name = "coexpression_plot",
    download_type = "plot",
    plot_params = list(file_type = "pdf", width = 12, height = 8, dpi = 300)
  )

  output$download_coexpr_de_multiple <- downloadHandler(
    filename = function() paste0("CoexprDE_",
                                 paste(input$coexpr_de_cluster_multiple, collapse = "-"),
                                 "_", format(Sys.time(), "%Y%m%d"), ".csv"),
    content  = function(file) {
      req(coexpr_de_results_multiple())
      write.csv(coexpr_de_results_multiple()$markers, file, row.names = FALSE)
    }
  )

  #####################Exclusive Biomarkers Discovery - Multiple Datasets###############################

  # ---- Reactive state ----

  # Initialize reactive variable for exclusive markers
  exclusive_markers_merge <- reactiveVal()

  # ---- Observers and computations ----

  # Populate group-by column selector: any categorical metadata column
  observe(run_safely("Interface update", persistent = FALSE, {
    req(multiple_datasets_object())
    meta <- multiple_datasets_object()@meta.data
    categorical_cols <- colnames(meta)[sapply(meta, function(x) is.character(x) || is.factor(x))]
    updateSelectInput(session, "exclusive_groupby_col_merge", choices = categorical_cols,
                      selected = if ("dataset" %in% categorical_cols) "dataset" else categorical_cols[1])
  }))

  # Populate Target group(s) and Compare Against group(s) when group-by column changes
  observeEvent(input$exclusive_groupby_col_merge, run_safely("exclusive groupby col merge", {
    req(multiple_datasets_object(), nzchar(input$exclusive_groupby_col_merge))
    obj <- multiple_datasets_object()
    col <- input$exclusive_groupby_col_merge
    req(col %in% colnames(obj@meta.data))
    values <- sort(unique(as.character(obj@meta.data[[col]])))
    updateSelectInput(session, "exclusive_target_groups_merge", choices = values, selected = NULL)
    updateSelectInput(session, "exclusive_compare_groups_merge", choices = values, selected = NULL)
  }))

  # Keep Compare Against choices distinct from Target selection
  observeEvent(input$exclusive_target_groups_merge, run_safely("exclusive target groups merge", {
    req(multiple_datasets_object(), nzchar(input$exclusive_groupby_col_merge))
    obj <- multiple_datasets_object()
    col <- input$exclusive_groupby_col_merge
    req(col %in% colnames(obj@meta.data))
    values <- sort(unique(as.character(obj@meta.data[[col]])))
    updateSelectInput(session, "exclusive_compare_groups_merge",
                      choices = setdiff(values, input$exclusive_target_groups_merge))
  }))

  # Populate cluster selector from current active cluster column
  observe(run_safely("Interface update", persistent = FALSE, {
    req(multiple_datasets_object())
    updateSelectInput(session, "exclusive_target_cluster_merge",
                      choices = active_clusters_merge())
  }))

  # Find exclusive biomarkers when button is clicked
  observeEvent(input$find_exclusive_markers_merge, {
    # Missing settings are reported instead of stopping silently (req() gave no feedback after the click)
    if (is.null(multiple_datasets_object())) { showNotification("Load and integrate datasets first.", type = "warning"); return() }
    if (is.null(input$exclusive_groupby_col_merge) || !nzchar(input$exclusive_groupby_col_merge)) { showNotification("Select a grouping column.", type = "warning"); return() }
    if (length(input$exclusive_target_groups_merge) == 0) { showNotification("Select at least one target group.", type = "warning"); return() }
    tryCatch({
      showModal(modalDialog(
        title = "Finding Exclusive Biomarkers",
        "Analyzing gene expression patterns across integrated datasets...",
        easyClose = FALSE,
        footer = NULL
      ))
      seurat_obj <- multiple_datasets_object()
      groupby_col <- input$exclusive_groupby_col_merge
      target_groups <- input$exclusive_target_groups_merge
      compare_groups <- if (!is.null(input$exclusive_compare_groups_merge) && length(input$exclusive_compare_groups_merge) > 0) {
        input$exclusive_compare_groups_merge
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
      if (!is.null(input$exclusive_target_cluster_merge) && length(input$exclusive_target_cluster_merge) > 0) {
        cluster_col   <- active_cluster_column_merge()
        if (!is.null(cluster_col) && cluster_col %in% colnames(seurat_obj@meta.data)) {
          cluster_cells <- colnames(seurat_obj)[as.character(seurat_obj@meta.data[[cluster_col]]) %in% input$exclusive_target_cluster_merge]
        } else {
          cluster_cells <- colnames(seurat_obj)[as.character(Idents(seurat_obj)) %in% input$exclusive_target_cluster_merge]
        }
        cells_to_keep <- intersect(cells_to_keep, cluster_cells)
      }
      if (length(cells_to_keep) == 0) {
        showNotification("No cells found for the selected group/cluster combination!", type = "error")
        removeModal()
        return()
      }
      seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
      # Detection rates and fold changes on the log-normalized RNA assay (the integrated assay has ~2,000 corrected genes)
      seurat_obj <- prepare_rna_assay_for_de(seurat_obj)$obj
      Idents(seurat_obj) <- factor(as.character(seurat_obj@meta.data[[groupby_col]]))
      cluster_text <- if (!is.null(input$exclusive_target_cluster_merge) && length(input$exclusive_target_cluster_merge) > 0) {
        paste(input$exclusive_target_cluster_merge, collapse = "_")
      } else {
        "AllClusters"
      }
      message(paste0("Comparing ", groupby_col, " [", paste(target_groups, collapse = ", "),
                     "] vs [", paste(compare_groups, collapse = ", "), "]",
                     if (cluster_text != "AllClusters") paste0(" | clusters: ", cluster_text) else ""))
      markers <- find_exclusive_biomarkers(
        seurat_obj = seurat_obj,
        target_cluster = target_groups,
        min_pct_target = input$exclusive_min_pct_target_merge,
        max_pct_other = input$exclusive_max_pct_other_merge,
        min_log2fc = input$exclusive_min_log2fc_merge,
        detection_threshold = input$exclusive_detection_threshold_merge,
        min_mean_expr_target = input$exclusive_min_mean_expr_merge,
        statistical_test = input$exclusive_statistical_test_merge,
        max_pvalue = input$exclusive_max_pvalue_merge,
        pvalue_adjustment = "BH",
        assay_name = "RNA",
        top_n = NULL,
        verbose = TRUE
      )
      removeModal()
      if (nrow(markers) == 0) {
        showNotification("No exclusive biomarkers found with current parameters. Try relaxing the filters.", 
                         type = "error", duration = 10)
        exclusive_markers_merge(NULL)
        return()
      }
      # Copy whole p-value vectors (ifelse() with a scalar test returned only the first element, recycled to every gene)
      markers_copy_venn <- markers
      markers_copy_venn$p_val <- if ("pvalue" %in% colnames(markers_copy_venn)) markers_copy_venn$pvalue else NA_real_
      markers_copy_venn$p_val_adj <- if ("pvalue_adjusted" %in% colnames(markers_copy_venn)) markers_copy_venn$pvalue_adjusted else NA_real_
      markers_copy_venn$avg_log2FC <- markers_copy_venn$log2_fold_change
      table_name <- paste0("ExclusiveMarkers_", groupby_col, "_",
                           paste(target_groups, collapse = "_"), "_vs_",
                           paste(compare_groups, collapse = "_"), "_",
                           cluster_text, "_", format(Sys.time(), "%H%M%S"))
      description <- paste0("Exclusive markers: ", groupby_col,
                            " [", paste(target_groups, collapse = ", "), "] vs [",
                            paste(compare_groups, collapse = ", "), "]",
                            if (cluster_text != "AllClusters") paste0(" | clusters: ", cluster_text) else "")
      parameters <- list(
        groupby_col     = groupby_col,
        target_groups   = target_groups,
        compare_groups  = compare_groups,
        target_clusters = input$exclusive_target_cluster_merge,
        min_pct_target  = input$exclusive_min_pct_target_merge,
        max_pct_other   = input$exclusive_max_pct_other_merge,
        min_log2fc      = input$exclusive_min_log2fc_merge,
        detection_threshold = input$exclusive_detection_threshold_merge,
        statistical_test = input$exclusive_statistical_test_merge
      )
      result <- storeDETable(gene_table_storage(), markers_copy_venn, table_name,
                             description, "exclusive_markers", parameters)
      if (result$success) gene_table_storage(result$storage)
      # p-values stay numeric in the display table: formatted client-side so column filters render sliders
      cleaned_gene_names <- clean_gene_names_for_html(markers$gene)
      markers$gene_link <- paste0('<a href="#" class="gene-name" data-gene="',
                                  cleaned_gene_names, '">', cleaned_gene_names, '</a>')
      if (input$exclusive_statistical_test_merge != "none") {
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
      exclusive_markers_merge(markers_display)
      showNotification(paste("Found", nrow(markers), "exclusive biomarkers!"), type = "message")
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error finding exclusive biomarkers:", conditionMessage(e)[1]),
                       type = "error", duration = 10)
      message("=== EXCLUSIVE BIOMARKERS ERROR ===\n", conditionMessage(e)[1])
    })
  })
  # ---- Tables ----

  # Display results table — all numeric columns (incl. p-values) get range sliders
  output$table_exclusive_markers_merge <- renderDataTable({
    req(exclusive_markers_merge())
    df <- exclusive_markers_merge()
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
  output$download_exclusive_markers_merge <- createDownloadHandler(
    reactive_data = reactive({
      markers <- exclusive_markers_merge()
      req(markers)
      
      # Remove HTML links for download
      markers$Gene <- gsub('.*data-gene="([^"]+)".*', '\\1', markers$Gene)
      return(markers)
    }),
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets") 
    }),
    data_name = reactive({ 
      dataset_text <- if (!is.null(input$exclusive_target_dataset_merge) && length(input$exclusive_target_dataset_merge) > 0) {
        paste0("Dataset", paste(input$exclusive_target_dataset_merge, collapse = "_"), "_")
      } else {
        ""
      }
      paste0("ExclusiveMarkers_", dataset_text, "Cluster", paste(input$exclusive_target_cluster_merge, collapse = "_"))
    }),
    download_type = "csv"
  )

  ############################## Pathway enrichment ##############################

  # ---- Reactive state ----

  # ORA / GSEA on a new comparison (FindMarkers, logfc.threshold = 0) or ORA on a stored DE table; helpers in enrichment_functions.R
  enrichment_results_merge <- reactiveVal(NULL)

  enrichment_gmt_merge <- reactiveVal(NULL)

  # ---- Observers and computations ----

  # A new object invalidates the results computed on the previous one
  observeEvent(multiple_datasets_object(), run_safely("Enrichment reset", persistent = FALSE, { enrichment_results_merge(NULL) }), ignoreInit = TRUE)

  # Grouping columns: metadata columns with 2 to 200 values
  observe(run_safely("Enrichment interface update", persistent = FALSE, {
    obj <- multiple_datasets_object(); req(obj)
    cols <- enrichment_grouping_columns(obj@meta.data)
    current <- isolate(input$enrichment_group_col_merge)
    default <- if ("seurat_clusters" %in% cols) "seurat_clusters" else cols[1]
    updateSelectInput(session, "enrichment_group_col_merge", choices = cols, selected = if (!is.null(current) && current %in% cols) current else default)
  }))

  # Group and restriction values follow the grouping column
  observeEvent(list(multiple_datasets_object(), input$enrichment_group_col_merge), run_safely("Enrichment interface update", persistent = FALSE, {
    obj <- multiple_datasets_object(); req(obj, input$enrichment_group_col_merge, input$enrichment_group_col_merge %in% colnames(obj@meta.data))
    values <- enrichment_group_values(obj@meta.data[[input$enrichment_group_col_merge]])
    updateSelectizeInput(session, "enrichment_group1_merge", choices = values, selected = character(0))
    updateSelectizeInput(session, "enrichment_group2_merge", choices = values, selected = character(0))
    # Restriction to clusters of the active identities (annotated names when the clusters were renamed)
    updateSelectizeInput(session, "enrichment_restrict_merge", choices = levels(droplevels(as.factor(Idents(obj)))), selected = character(0))
  }))

  # Stored DE tables (same storage as the Venn diagrams)
  observe(run_safely("Enrichment interface update", persistent = FALSE, {
    updateSelectInput(session, "enrichment_stored_table_merge", choices = stored_table_choices(gene_table_storage()))
  }))

  # GMT file read on upload and matched to the dataset genes; the md5 checksum identifies the file version in the parameters
  observeEvent(input$enrichment_gmt_file_merge, run_safely("Gene set file", {
    obj <- multiple_datasets_object(); req(obj)
    sets <- match_gene_sets_to_dataset(read_gmt_gene_sets(input$enrichment_gmt_file_merge$datapath), rownames(obj))
    if (sets$n_sets == 0) stop("No gene of the file matches the dataset genes")
    enrichment_gmt_merge(c(sets, list(file = input$enrichment_gmt_file_merge$name, md5 = unname(tools::md5sum(input$enrichment_gmt_file_merge$datapath)))))
  }))

  # Sample column of the pseudobulk option: categorical metadata columns, 'dataset' by default
  observe(run_safely("Enrichment interface update", persistent = FALSE, {
    obj <- multiple_datasets_object(); req(obj)
    cols <- enrichment_grouping_columns(obj@meta.data)
    current <- isolate(input$enrichment_sample_col_merge)
    default <- if ("dataset" %in% cols) "dataset" else if ("orig.ident" %in% cols) "orig.ident" else cols[1]
    updateSelectInput(session, "enrichment_sample_col_merge", choices = cols, selected = if (!is.null(current) && current %in% cols) current else default)
  }))

  observeEvent(input$enrichment_run_merge, run_safely("Pathway enrichment", {
    obj <- multiple_datasets_object()
    if (is.null(obj)) stop("Load a dataset first")
    settings <- list(source = input$enrichment_source_merge, group_col = input$enrichment_group_col_merge, group1 = input$enrichment_group1_merge,
                     group2 = input$enrichment_group2_merge, restrict_values = input$enrichment_restrict_merge, min_pct = input$enrichment_min_pct_merge,
                     stored_table = input$enrichment_stored_table_merge, gene_sets = input$enrichment_gene_sets_merge, ontology = input$enrichment_ontology_merge,
                     min_size = input$enrichment_min_size_merge, max_size = input$enrichment_max_size_merge, methods = input$enrichment_methods_merge,
                     direction = input$enrichment_direction_merge, ora_max_p_adj = input$enrichment_ora_max_p_adj_merge, ora_min_log2fc = input$enrichment_ora_min_log2fc_merge,
                     de_method = input$enrichment_de_method_merge, sample_col = input$enrichment_sample_col_merge, pb_min_cells = input$enrichment_pb_min_cells_merge)
    showModal(modalDialog(title = "Pathway enrichment", "Testing gene sets...", easyClose = FALSE, footer = NULL))
    res <- tryCatch(run_enrichment_analysis(obj, settings, stored_tables = gene_table_storage(), gmt = enrichment_gmt_merge()), finally = removeModal())
    enrichment_results_merge(res)
    # Every run is stored with the DE tables, like the other analyses: saved with the session (Save & Session, where
    # unwanted runs can be deleted) and summarized with its parameters in the CSV / PDF export
    tables <- gene_table_storage()
    tables[[paste0("Enrichment_", format(Sys.time(), "%Y%m%d_%H%M%S"))]] <- enrichment_storage_entry(res, input$enrichment_show_max_p_adj_merge %||% 0.05)
    gene_table_storage(tables)
    if (res$gsea_skipped) showNotification("GSEA needs all tested genes ranked: it runs only on a new comparison. ORA only for stored tables.", type = "warning", duration = 8)
    showNotification(paste0(enrichment_done_message(res, input$enrichment_show_max_p_adj_merge %||% 0.05), ". Results stored (Save & Session)."), type = "message", duration = 8)
  }))

  enrichment_ora_plot_obj_merge <- reactive({
    r <- enrichment_results_merge(); shiny::validate(shiny::need(!is.null(r$ora), "Run ORA to see results"))
    enrichment_plot_or_message(r$ora$table, "ORA", input$enrichment_show_max_p_adj_merge %||% 0.05, input$enrichment_ora_n_terms_merge %||% 15, paste("ORA —", r$label), r$gene_set_label)
  })

  enrichment_gsea_plot_obj_merge <- reactive({
    r <- enrichment_results_merge(); shiny::validate(shiny::need(!is.null(r$gsea), "Run GSEA (new comparison) to see results"))
    enrichment_plot_or_message(r$gsea$table, "GSEA", input$enrichment_show_max_p_adj_merge %||% 0.05, input$enrichment_gsea_n_terms_merge %||% 15, paste("GSEA —", r$label), r$gene_set_label)
  })

  # ---- UI and text outputs ----

  output$enrichment_gmt_info_merge <- renderText({
    s <- enrichment_gmt_merge(); req(s)
    sprintf("%s: %d sets, %.0f%% of set genes found in the dataset", s$file, s$n_sets, 100 * s$fraction_matched)
  })

  output$enrichment_summary_merge <- renderText({ r <- enrichment_results_merge(); req(r); enrichment_summary_text(r) })

  # ---- Plots ----
  output$enrichment_ora_plot_merge <- renderPlot({ enrichment_ora_plot_obj_merge() }, res = 96)
  output$enrichment_gsea_plot_merge <- renderPlot({ enrichment_gsea_plot_obj_merge() }, res = 96)

  # ---- Tables ----
  output$enrichment_ora_table_merge <- DT::renderDT({ r <- enrichment_results_merge(); req(r$ora, nrow(r$ora$table) > 0); enrichment_datatable(r$ora$table, "ORA", input$enrichment_show_max_p_adj_merge %||% 0.05, isTRUE(input$enrichment_show_all_merge)) })
  output$enrichment_gsea_table_merge <- DT::renderDT({ r <- enrichment_results_merge(); req(r$gsea, nrow(r$gsea$table) > 0); enrichment_datatable(r$gsea$table, "GSEA", input$enrichment_show_max_p_adj_merge %||% 0.05, isTRUE(input$enrichment_show_all_merge)) })
  output$enrichment_params_table_merge <- renderTable({ r <- enrichment_results_merge(); req(r); r$params })

  # ---- Downloads ----

  output$download_enrichment_ora_plot_merge <- downloadHandler(filename = function() paste0(enrichment_file_stub(enrichment_results_merge()$label), "_ORA.png"),
                                                               content = function(file) ggplot2::ggsave(file, enrichment_ora_plot_obj_merge(), width = 10, height = 7, dpi = 300, bg = "white"))

  output$download_enrichment_gsea_plot_merge <- downloadHandler(filename = function() paste0(enrichment_file_stub(enrichment_results_merge()$label), "_GSEA.png"),
                                                                content = function(file) ggplot2::ggsave(file, enrichment_gsea_plot_obj_merge(), width = 10, height = 7, dpi = 300, bg = "white"))

  # Result tables as Excel files: all tested terms sorted by signed significance, shaded red (up) / blue (down), parameters in a second sheet
  output$download_enrichment_ora_table_merge <- downloadHandler(filename = function() paste0(enrichment_file_stub(enrichment_results_merge()$label), "_ORA.xlsx"),
                                                                content = function(file) { r <- enrichment_results_merge(); req(r$ora); write_enrichment_xlsx(r, "ORA", file, input$enrichment_show_max_p_adj_merge %||% 0.05) })

  output$download_enrichment_gsea_table_merge <- downloadHandler(filename = function() paste0(enrichment_file_stub(enrichment_results_merge()$label), "_GSEA.xlsx"),
                                                                 content = function(file) { r <- enrichment_results_merge(); req(r$gsea); write_enrichment_xlsx(r, "GSEA", file, input$enrichment_show_max_p_adj_merge %||% 0.05) })

  output$download_enrichment_params_merge <- downloadHandler(filename = function() paste0(enrichment_file_stub(enrichment_results_merge()$label), "_parameters.csv"),
                                                             content = function(file) { r <- enrichment_results_merge(); req(r); utils::write.csv(r$params, file, row.names = FALSE) })

  #####################Exclusive Biomarkers Visualization - Multiple Datasets###############################

  # ---- Reactive state ----

  # Reactive value to store the plot
  exclusive_plot_merge <- reactiveVal()

  # ---- Observers and computations ----

  observeEvent(input$generate_exclusive_plot_merge, {
    req(multiple_datasets_object())
    tryCatch({
      # Genes present in the default assay (case-insensitive), at most 10
      seurat_obj <- multiple_datasets_object()
      check <- validate_gene_input(input$exclusive_genes_to_plot_merge, seurat_obj, DefaultAssay(seurat_obj))
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (!is.null(check$message)) showNotification(check$message, type = "warning", duration = 8)
      genes_input <- check$genes
      if (length(genes_input) > 10) {
        showNotification("Maximum 10 genes for visualization: only the first 10 are plotted.", type = "warning")
        genes_input <- genes_input[1:10]
      }
      plot_type <- input$exclusive_plot_type_merge %||% "dotplot"
      # Optional split column, shared by the three plot types
      split_column <- input$exclusive_split_by_column_merge
      if (is.null(split_column) || !nzchar(split_column) || !split_column %in% colnames(seurat_obj@meta.data)) split_column <- NULL
      # Cell filter on the split values, applied before the modal opens so that an early return never leaves it on screen
      plot_seurat <- seurat_obj
      if (!is.null(split_column) && length(input$exclusive_split_values_filter) > 0) {
        cells_to_keep <- colnames(seurat_obj)[as.character(seurat_obj@meta.data[[split_column]]) %in% input$exclusive_split_values_filter]
        if (length(cells_to_keep) == 0) {
          showNotification("No cells match the selected filter values", type = "warning")
          return()
        }
        plot_seurat <- subset(seurat_obj, cells = cells_to_keep)
        message(paste("Filtered to", length(cells_to_keep), "cells"))
      }
      showModal(modalDialog(
        title = "Generating Plot",
        paste("Creating", plot_type, "for", length(genes_input), "gene(s)..."),
        easyClose = FALSE,
        footer = NULL
      ))
      title_text <- if (!is.null(split_column)) paste("Expression of Exclusive Biomarkers by", split_column) else NULL
      title_theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), legend.position = "right")
      axis_theme <- theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 10))
      plot_obj <- if (plot_type == "dotplot") {
        # With split.by, DotPlot needs one colour per split value: two colours fail from 3 values on
        dot_cols <- if (is.null(split_column)) c("lightgrey", "blue") else scales::hue_pal()(length(unique(as.character(plot_seurat@meta.data[[split_column]]))))
        DotPlot(plot_seurat, features = genes_input, cols = dot_cols, dot.scale = 8, split.by = split_column) + RotatedAxis() +
          labs(title = title_text %||% "Expression of Exclusive Biomarkers across Clusters", x = "Genes", y = "Cluster") + title_theme + axis_theme
      } else if (plot_type == "violin") {
        VlnPlot(plot_seurat, features = genes_input, split.by = split_column, pt.size = 0) +
          labs(title = title_text %||% "Expression of Exclusive Biomarkers across Clusters", x = "Cluster", y = "Expression Level") + title_theme + axis_theme
      } else {
        FeaturePlot(plot_seurat, features = genes_input, split.by = split_column, pt.size = 0.5, ncol = 2) +
          labs(title = title_text %||% "Expression of Exclusive Biomarkers on UMAP") + title_theme
      }
      removeModal()
      exclusive_plot_merge(plot_obj)
      showNotification("Plot generated successfully!", type = "message")
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error generating plot:", e$message), type = "error", duration = 10)
      message("Exclusive biomarker plot error: ", e$message)
    })
  })

  # ---- UI and text outputs ----

  # Dynamic UI to filter values of the split column
  output$exclusive_split_values_filter_ui <- renderUI({
    req(multiple_datasets_object(), input$exclusive_split_by_column_merge)
    
    # Only show if a split column is selected
    if (is.null(input$exclusive_split_by_column_merge) || 
        input$exclusive_split_by_column_merge == "") {
      return(NULL)
    }
    
    split_col <- input$exclusive_split_by_column_merge
    seurat_obj <- multiple_datasets_object()
    
    # Check column exists
    if (!split_col %in% colnames(seurat_obj@meta.data)) {
      return(NULL)
    }
    
    # Get unique values from the split column
    available_values <- unique(seurat_obj@meta.data[[split_col]])
    available_values <- sort(available_values[!is.na(available_values)])
    
    if (length(available_values) == 0) {
      return(NULL)
    }
    
    selectizeInput("exclusive_split_values_filter",
                   paste("Filter", split_col, "values:"),
                   choices = available_values,
                   selected = NULL,  # Empty = show all
                   multiple = TRUE,
                   options = list(
                     plugins = list('remove_button'),
                     placeholder = 'Select values (empty = show all)'
                   ))
  })

  # ---- Plots ----

  # Render the plot
  output$plot_exclusive_markers_merge <- renderPlot({
    req(exclusive_plot_merge())
    print(exclusive_plot_merge())
  }, height = 600)

  # ---- Downloads ----

  # Download handler for the plot
  output$download_exclusive_plot_merge <- downloadHandler(
    filename = function() {
      object_name <- getObjectNameForDownload(multiple_datasets_object(), default_name = "ExclusiveMarkers")
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      genes_text <- substr(paste(parse_gene_tokens(input$exclusive_genes_to_plot_merge), collapse = "_"), 1, 50)  
      
      paste0(object_name, "_ExclusiveMarkers_", genes_text, "_", 
             input$exclusive_plot_type_merge, "_", timestamp, ".", 
             input$exclusive_plot_format_merge)
    },
    content = function(file) {
      req(exclusive_plot_merge())
      
      tryCatch({
        plot_obj <- exclusive_plot_merge()
        
        n_genes <- max(1, min(length(parse_gene_tokens(input$exclusive_genes_to_plot_merge)), 10))
        # Adjust dimensions based on plot type and split
        # Larger canvas when the plot is split by a metadata column
        split_col <- input$exclusive_split_by_column_merge
        multiplier <- if (!is.null(split_col) && nzchar(split_col)) 1.5 else 1
        
        if (input$exclusive_plot_type_merge == "dotplot") {
          width <- max(8, n_genes * 1.2) * multiplier
          height <- 6 * multiplier
        } else if (input$exclusive_plot_type_merge == "violin") {
          ncol <- min(3, n_genes)
          nrow <- ceiling(n_genes / ncol)
          width <- ncol * 4 * multiplier
          height <- nrow * 4 * multiplier
        } else {  # feature
          ncol <- min(3, n_genes)
          nrow <- ceiling(n_genes / ncol)
          width <- ncol * 5 * multiplier
          height <- nrow * 4 * multiplier
        }
        
        # Save based on format
        if (input$exclusive_plot_format_merge == "pdf") {
          pdf(file, width = width, height = height)
          print(plot_obj)
          dev.off()
        } else if (input$exclusive_plot_format_merge == "svg") {
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
            dpi = input$exclusive_plot_dpi_merge,
            device = input$exclusive_plot_format_merge
          )
        }
        
        showNotification("Plot downloaded successfully!", type = "message")
        
      }, error = function(e) {
        showNotification(paste("Error downloading plot:", e$message), type = "error")
        print(paste("Download error:", e$message))
      })
    }
  )

  ############################# Venn Diagram Comparison #############################

  # ---- Reactive state ----

  # Initialize reactive storage for Venn diagrams
  venn_plot_rendered <- reactiveVal(NULL)

  # ---- Observers and computations ----

  # Update select inputs with available gene tables
  observe(run_safely("Interface update", persistent = FALSE, {
    updateVennSelectInputs(session, gene_table_storage(), 
                           c("venn_table_1", "venn_table_2", "venn_table_3"))
  }))

  # Enable/disable generate button based on available tables
  observe(run_safely("Interface update", persistent = FALSE, {
    tables <- gene_table_storage()
    if (length(tables) < 2) {
      shinyjs::disable("generate_venn_btn")
    } else {
      shinyjs::enable("generate_venn_btn")
    }
  }))

  # Generate Venn diagram using modular functions
  observeEvent(input$generate_venn_btn, run_safely("generate venn btn", {
    showModal(modalDialog(title = "Generating Venn Diagram", "Processing...", 
                          easyClose = FALSE, footer = NULL))
    
    # Use modular function for complete Venn generation
    result <- processVennGeneration(
      table_storage = gene_table_storage(),
      selected_tables = c(input$venn_table_1, input$venn_table_2, input$venn_table_3),
      filter_params = list(
        significant_only = c(input$significant_only_venn_1, input$significant_only_venn_2, 
                             input$significant_only_venn_3),
        log_fc_threshold = c(input$log_fc_threshold_venn_1, input$log_fc_threshold_venn_2, 
                             input$log_fc_threshold_venn_3),
        p_val_threshold = input$p_val_threshold_venn,
        use_adjusted_p = input$use_adjusted_p_venn,
        direction = input$venn_direction
      ),
      colors = c(input$venn_color_1, input$venn_color_2, input$venn_color_3)
    )
    
    removeModal()
    
    if (result$success) {
      # Store results using modular structure
      venn_plot_rendered(result$venn_plot)
      current_gene_lists(result$overlaps)
      updateSelectInput(session, "selected_gene_set", choices = names(result$overlaps))
      shinyjs::enable("download_venn_diagram")
    } else {
      showNotification(result$message, type = "error")
    }
  }))

  # ---- Plots ----

  # Render Venn diagram
  output$venn_plot <- renderPlot({
    req(venn_plot_rendered())
    grid.draw(venn_plot_rendered())
  })

  # ---- Tables ----

  # Display gene table for selected overlap
  output$venn_gene_table <- renderDT({
    req(current_gene_lists(), input$selected_gene_set)
    overlaps <- current_gene_lists()
    selected_genes <- overlaps[[input$selected_gene_set]]
    
    if (length(selected_genes) == 0) {
      return(data.frame(Gene = character(0)))
    }
    
    gene_df <- data.frame(Gene = selected_genes)
    datatable(gene_df, 
              options = list(pageLength = 15, scrollX = TRUE, dom = 'Bfrtip', 
                             buttons = c('copy', 'csv', 'excel')), 
              rownames = FALSE)
  })

  # ---- Downloads ----

  # Download handlers using modular functions
  output$download_venn_diagram <- createDownloadHandler(
    reactive_data = venn_plot_rendered,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "VennDiagram") 
    }),
    data_name = "venn_comparison",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$venn_diagram_format }),
      width = 8,
      height = 6,
      dpi = reactive({ input$venn_diagram_dpi })
    )
  )

  output$download_venn_gene_lists <- createDownloadHandler(
    reactive_data = reactive({
      gene_lists <- current_gene_lists()
      if (is.null(gene_lists)) return(NULL)
      
      # Convert gene lists to data frame
      all_genes <- data.frame()
      for (list_name in names(gene_lists)) {
        genes <- gene_lists[[list_name]]
        if (length(genes) > 0) {
          df <- data.frame(Gene = genes, Gene_Set = list_name, stringsAsFactors = FALSE)
          all_genes <- rbind(all_genes, df)
        } else {
          df <- data.frame(Gene = "No genes", Gene_Set = list_name, stringsAsFactors = FALSE)
          all_genes <- rbind(all_genes, df)
        }
      }
      return(all_genes)
    }),
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "VennDiagram") }),
    data_name = "gene_lists",
    download_type = "csv"
  )

  ############################# Helper Functions for DE Table Storage #############################

  # Function to store DE tables using modular approach (storeDETable puts the gene names in the row names and notifies)
  store_de_table <- function(table_data, table_name, description, type, parameters) {
    result <- storeDETable(gene_table_storage(), table_data, table_name, description, type, parameters)
    if (isTRUE(result$success)) gene_table_storage(result$storage)
    invisible(result)
  }

  ############################# Pie Chart for Cluster Composition #############################

  # ---- Reactive state ----
  cluster_composition_data_merge <- reactiveVal(NULL)
  cluster_pie_plot_merge <- reactiveVal(NULL)

  # ---- Observers and computations ----

  # Update cluster choices
  observe(run_safely("Interface update", persistent = FALSE, {
    req(multiple_datasets_object())
    clusters <- levels(Idents(multiple_datasets_object()))
    updateSelectInput(session, "cluster_pie_select_merge", choices = clusters)
  }))

  # Generate pie chart data when selections change
  observeEvent(c(multiple_datasets_object(), input$cluster_pie_select_merge, 
                 input$metadata_column_pie, input$normalize_by_dataset_size_merge), {
                   req(multiple_datasets_object(), input$cluster_pie_select_merge, input$metadata_column_pie)
                   
                   tryCatch({
                     obj <- multiple_datasets_object()
                     selected_cluster <- input$cluster_pie_select_merge
                     metadata_col <- input$metadata_column_pie
                     
                     # Validate metadata column exists
                     if (!metadata_col %in% colnames(obj@meta.data)) {
                       showNotification(paste("Metadata column", metadata_col, "not found!"), type = "warning")
                       cluster_composition_data_merge(NULL)
                       return()
                     }
                     
                     # Get cells in selected cluster
                     cluster_cells <- names(Idents(obj))[Idents(obj) == selected_cluster]
                     
                     if(length(cluster_cells) == 0) {
                       cluster_composition_data_merge(NULL)
                       return()
                     }
                     
                     # Count cells by metadata group
                     cluster_counts <- table(obj@meta.data[cluster_cells, metadata_col])
                     total_counts <- table(obj@meta.data[[metadata_col]])
                     
                     if(input$normalize_by_dataset_size_merge) {
                       # Calculate enrichment scores (normalized by group size)
                       enrichment_scores <- as.numeric(cluster_counts) / as.numeric(total_counts[names(cluster_counts)])
                       
                       # Renormalize to sum to 100%
                       enrichment_pcts <- enrichment_scores / sum(enrichment_scores) * 100
                       
                       pie_data <- data.frame(
                         group = names(cluster_counts),
                         count = as.numeric(cluster_counts),
                         total_in_group = as.numeric(total_counts[names(cluster_counts)]),
                         enrichment_score = round(enrichment_scores, 3),
                         percentage = round(enrichment_pcts, 1),
                         raw_pct = round(as.numeric(cluster_counts) / sum(cluster_counts) * 100, 1),
                         stringsAsFactors = FALSE
                       )
                       
                     } else {
                       # Simple raw distribution
                       pie_data <- data.frame(
                         group = names(cluster_counts),
                         count = as.numeric(cluster_counts),
                         percentage = round(as.numeric(cluster_counts) / sum(cluster_counts) * 100, 1),
                         stringsAsFactors = FALSE
                       )
                     }
                     
                     cluster_composition_data_merge(pie_data)
                     
                   }, error = function(e) {
                     showNotification(paste("Error generating pie chart data:", e$message), type = "error")
                     cluster_composition_data_merge(NULL)
                   })
                 })

  # ---- Plots ----

  # Render pie chart
  output$cluster_composition_pie_merge <- renderPlot({
    req(cluster_composition_data_merge(), input$metadata_column_pie)
    
    pie_data <- cluster_composition_data_merge()
    is_normalized <- input$normalize_by_dataset_size_merge
    metadata_col <- input$metadata_column_pie
    
    # Create base pie chart using 'group' column
    pie_plot <- ggplot(pie_data, aes(x = "", y = percentage, fill = group)) +
      geom_col(width = 1) +
      coord_polar("y", start = 0) +
      theme_void() +
      labs(
        title = input$cluster_pie_select_merge,
        subtitle = paste0("Grouped by: ", metadata_col, " | Total cells: ", sum(pie_data$count)),
        fill = metadata_col
      ) +
      geom_text(
        aes(label = paste0(percentage, "%\n(", count, " cells)")), 
        position = position_stack(vjust = 0.5), 
        color = "white", 
        size = if(is_normalized) 3.5 else 4, 
        fontface = "bold"
      )
    
    # Apply base theme
    pie_plot <- pie_plot +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.position = "bottom",
        legend.title = element_text(face = "bold")
      )
    
    # Remove legend if requested
    if(input$remove_legend_pie_merge) {
      pie_plot <- pie_plot + theme(legend.position = "none")
    }
    
    # Dark mode
    if(isTRUE(input$dark_mode_pie_chart_merge)) {
      pie_plot <- pie_plot +
        theme(
          panel.background = element_rect(fill = "black", color = "black"),
          plot.background = element_rect(fill = "black", color = "black"),
          text = element_text(color = "white"),
          plot.title = element_text(color = "white"),
          plot.subtitle = element_text(color = "white"),
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"),
          legend.background = element_rect(fill = "black", color = "black"),
          legend.key = element_rect(fill = "black", color = "black")
        )
    }
    
    cluster_pie_plot_merge(pie_plot)
    pie_plot
  })

  # ---- Downloads ----

  # Download handler for pie chart
  output$download_pie_chart_merge <- createDownloadHandler(
    reactive_data = cluster_pie_plot_merge,
    object_name_reactive = reactive({
      cluster_name <- input$cluster_pie_select_merge %||% "cluster"
      metadata_col <- input$metadata_column_pie %||% "metadata"
      normalize_text <- ifelse(input$normalize_by_dataset_size_merge, "_normalized", "_raw")
      paste0(cluster_name, "_by_", metadata_col, normalize_text)
    }),
    data_name = "cluster_composition_pie",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$pie_format_merge }),
      width = 10,
      height = 8,
      dpi = 300
    )
  )

  ########################## Volcano Plot - Multiple Datasets ##########################

  # ---- Reactive state ----

  # Reactive value for storing volcano plot
  volcano_plot_merge <- reactiveVal(NULL)

  # ---- Observers and computations ----

  # Update volcano plot source selector when DEG results are available (enrichment results left out)
  observe(run_safely("Interface update", persistent = FALSE, {
    tables <- de_tables_only(gene_table_storage())
    if (length(tables) > 0) {
      table_choices <- c("None" = "", names(tables))
      updateSelectInput(session, "volcano_deg_source_merge", choices = table_choices)
    } else {
      updateSelectInput(session, "volcano_deg_source_merge", choices = c("None" = ""))
    }
  }))

  # Generate volcano plot
  observeEvent(input$generate_volcano_btn_merge, {
    req(input$volcano_deg_source_merge, input$volcano_deg_source_merge != "")
    
    tryCatch({
      showModal(modalDialog(
        title = "Generating Volcano Plot", 
        "Processing...", 
        easyClose = FALSE, 
        footer = NULL
      ))
      
      # Get selected DEG results from storage
      tables <- gene_table_storage()
      table_entry <- tables[[input$volcano_deg_source_merge]]
      
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
        log2fc_threshold = input$volcano_log2fc_threshold_merge,
        pval_threshold = input$volcano_pval_threshold_merge,
        color_up = input$volcano_color_up_merge,
        color_down = input$volcano_color_down_merge,
        color_ns = input$volcano_color_ns_merge,
        label_top_genes = input$volcano_label_genes_merge,
        point_size = input$volcano_point_size_merge,
        point_alpha = input$volcano_point_alpha_merge
      )
      
      volcano_plot_merge(volcano)
      removeModal()
      showNotification("Volcano plot generated successfully!", type = "message")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error generating volcano plot:", e$message), type = "error")
      message("DEBUG volcano error: ", e$message)
    })
  })

  # ---- Plots ----

  # Render volcano plot
  output$volcano_plot_merge <- renderPlot({
    req(volcano_plot_merge())
    volcano_plot_merge()
  })

  # ---- Downloads ----

  # Download volcano plot
  output$download_volcano_plot_merge <- createDownloadHandler(
    reactive_data = volcano_plot_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "VolcanoPlot_Merge") 
    }),
    data_name = "volcano_plot",
    download_type = "plot",
    # reactive(): a plain input value is read once at the first download and then frozen
    plot_params = list(
      file_type = reactive({ input$volcano_plot_format_merge }),
      width = 10,
      height = 8,
      dpi = reactive({ input$volcano_plot_dpi_merge })
    )
  )

  ############################## Subseting seurat object ##############################

  # ---- Reactive state ----

  # Variable réactive pour cet onglet
  subset_seurat_merge <- reactiveVal(NULL)

  # ---- Observers and computations ----
  shinyjs::disable("download_subset_merge")

  # Mise à jour des choix d'identité de cellules pour le sous-ensemble
  observe(run_safely("Interface update", persistent = FALSE, {
    req(multiple_datasets_object())
    seurat_object <- multiple_datasets_object()
    
    tryCatch({
      # Convertir les identifiants en caractères
      cluster_idents <- Idents(seurat_object)
      cluster_choices <- as.character(unique(cluster_idents))
      cluster_choices <- cluster_choices[!is.na(cluster_choices)]
      
      if (length(cluster_choices) > 0) {
        updateSelectInput(session, "select_ident_subset_merge", 
                          choices = cluster_choices,
                          selected = cluster_choices[1])
      }
    }, error = function(e) {
      message(paste("Error updating subset choices:", e$message))
    })
  }))

  # Réinitialisation de l'objet subset_seurat_merge lors du changement de multiple_datasets_object
  observe(run_safely("Interface update", persistent = FALSE, {
    if (!is.null(multiple_datasets_object())) {
      subset_seurat_merge(multiple_datasets_object())
    }
  }))

  # Mise à jour de la sélection de sous-ensemble
  observeEvent(input$apply_subset_merge, {
    tryCatch({
      req(multiple_datasets_object())
      
      # Si aucun identifiant n'est sélectionné, retourne l'objet original
      if (length(input$select_ident_subset_merge) == 0) {
        subset_seurat_merge(multiple_datasets_object())
        shinyjs::enable("download_subset_merge")
        return()
      }
      
      subsetted_seurat <- subset(multiple_datasets_object(), idents = input$select_ident_subset_merge)
      
      # Vérifiez s'il y a des cellules dans le sous-ensemble
      if (nrow(subsetted_seurat@meta.data) == 0) {
        showNotification("No cells found with the selected identities.", type = "error")
        return()
      }
      
      subset_seurat_merge(subsetted_seurat)
      shinyjs::enable("download_subset_merge")
    }, error = function(e) {
      showNotification(paste0("Error updating subset selection: ", e$message), type = "error")
    })
  })

  observeEvent(input$apply_gene_subset_merge, {
    tryCatch({
      req(multiple_datasets_object())
      
      # Every gene must exist: the "number of genes to express" rule depends on the full list
      check <- validate_gene_input(input$gene_list_merge, multiple_datasets_object(), DefaultAssay(multiple_datasets_object()))
      if (!check$ok) { showNotification(check$message, type = "error", duration = 8); return() }
      if (length(check$missing) > 0) {
        showNotification(sprintf("Genes not found in assay '%s': %s", check$assay, paste(check$missing, collapse = ", ")), type = "error", duration = 8)
        return()
      }
      if (!is.null(check$message)) showNotification(check$message, type = "message", duration = 6)
      gene_list <- check$genes
      # Binary matrix: TRUE if cell expresses gene above threshold
      expression_matrix <- sapply(gene_list, function(gene) {
        FetchData(multiple_datasets_object(), vars = gene) >= input$expression_threshold_merge
      })
      
      n_genes_per_cell <- rowSums(expression_matrix)
      
      if (isTRUE(input$negative_gene_subset_merge)) {
        # Negative selection: keep cells NOT meeting the expression criteria
        cells_to_keep <- which(n_genes_per_cell < input$num_genes_to_express_merge)
        message(sprintf("Negative selection: keeping %d cells", length(cells_to_keep)))
      } else {
        # Positive selection: keep cells meeting the expression criteria
        cells_to_keep <- which(n_genes_per_cell >= input$num_genes_to_express_merge)
        message(sprintf("Positive selection: keeping %d cells", length(cells_to_keep)))
      }
      
      if (length(cells_to_keep) == 0) {
        showNotification("No cells found with the specified gene expression criteria.", type = "error")
        return()
      }
      
      subsetted_seurat <- subset(multiple_datasets_object(), cells = cells_to_keep)
      subset_seurat_merge(subsetted_seurat)
      shinyjs::enable("download_subset_merge")
      
    }, error = function(e) {
      showNotification(paste0("Error during gene-based subset: ", e$message), type = "error")
    })
  })

  # Mise à jour des choix de colonnes de métadonnées avec exclusions
  observe(run_safely("Interface update", persistent = FALSE, {
    req(multiple_datasets_object())
    
    # Liste complète des colonnes de métadonnées
    all_metadata_fields <- colnames(multiple_datasets_object()@meta.data)
    
    # Champs à exclure
    exclude_fields <- c( "nCount_RNA", "nCount_ATAC", "nFeature_RNA",
                         "nFeature_ATAC", "percent.mt")
    
    # Exclusion des champs qui contiennent certains motifs (optionnel)
    exclude_pattern <- "^snn_res|^pANN|^PC_|^RNA_snn|^ATAC_snn|^integrated_snn"
    exclude_fields <- c(exclude_fields, all_metadata_fields[grepl(exclude_pattern, all_metadata_fields)])
    
    # Filtrer les colonnes à afficher
    available_columns <- setdiff(all_metadata_fields, exclude_fields)
    
    # Mettre à jour le sélecteur
    updateSelectInput(session, "metadata_column_subset", choices = available_columns)
  }))

  # Gestion du bouton de subset par métadonnées
  observeEvent(input$apply_metadata_subset, {
    tryCatch({
      req(multiple_datasets_object(), input$metadata_column_subset)
      
      meta_col <- input$metadata_column_subset
      col_values <- multiple_datasets_object()@meta.data[[meta_col]]
      
      # Logique de filtrage différente selon le type de données
      if (is.numeric(col_values)) {
        # Pour les données numériques, utiliser la plage
        req(input$metadata_num_range)
        min_val <- input$metadata_num_range[1]
        max_val <- input$metadata_num_range[2]
        
        if (input$invert_metadata_selection) {
          cells_to_keep <- rownames(multiple_datasets_object()@meta.data)[col_values < min_val | col_values > max_val]
        } else {
          cells_to_keep <- rownames(multiple_datasets_object()@meta.data)[col_values >= min_val & col_values <= max_val]
        }
      } else {
        # Pour les données catégorielles, utiliser les valeurs sélectionnées
        req(input$metadata_cat_values)
        selected_values <- input$metadata_cat_values
        
        cells_to_keep <- rownames(multiple_datasets_object()@meta.data)[col_values %in% selected_values]
      }
      
      # Vérifier si des cellules correspondent au critère
      if (length(cells_to_keep) == 0) {
        showNotification("No cells found matching the metadata criteria.", type = "error")
        return()
      }
      
      # Créer le sous-ensemble
      subsetted_seurat <- subset(multiple_datasets_object(), cells = cells_to_keep)
      subset_seurat_merge(subsetted_seurat)
      
      # Activer le bouton de téléchargement
      shinyjs::enable("download_subset_merge")
      
      showNotification(paste0("Subset created with ", length(cells_to_keep), " cells."), type = "message")
    }, error = function(e) {
      showNotification(paste0("Error during metadata-based subset: ", e$message), type = "error")
    })
  })

  # ---- UI and text outputs ----

  # Interface dynamique pour les valeurs de métadonnées
  output$metadata_values_ui <- renderUI({
    req(multiple_datasets_object(), input$metadata_column_subset)
    
    # Obtenir les valeurs uniques pour la colonne sélectionnée
    col_values <- unique(multiple_datasets_object()@meta.data[[input$metadata_column_subset]])
    
    # Détecter le type de données
    if (is.numeric(col_values)) {
      # Pour les données numériques, offrir une plage
      min_val <- min(col_values, na.rm = TRUE)
      max_val <- max(col_values, na.rm = TRUE)
      
      tagList(
        sliderInput("metadata_num_range", "Value range:",
                    min = min_val, max = max_val,
                    value = c(min_val, max_val)),
        checkboxInput("invert_metadata_selection", "Invert selection", value = FALSE)
      )
    } else {
      # Pour les données catégorielles, offrir une liste de sélection
      selectInput("metadata_cat_values", "Select values:",
                  choices = sort(as.character(col_values)),
                  multiple = TRUE,
                  selected = sort(as.character(col_values))[1])
    }
  })

  # ---- Plots ----

  # Affichage de l'UMAP global
  output$global_umap_merge <- renderPlot({
    req(multiple_datasets_object())
    
    # Ne pas continuer si pas d'UMAP
    if(!"umap" %in% names(multiple_datasets_object()@reductions)) {
      return(NULL)
    }
    
    DimPlot(multiple_datasets_object(), group.by = "ident", label = TRUE, label.size = 4) +
      NoAxes() +
      NoLegend() +
      ggtitle(NULL)
  })

  # Affichage de l'UMAP du sous-ensemble
  output$subset_umap_merge <- renderPlot({
    req(subset_seurat_merge())
    
    # Ne pas continuer si pas d'UMAP
    if(!"umap" %in% names(subset_seurat_merge()@reductions)) {
      return(NULL)
    }
    
    DimPlot(subset_seurat_merge(), group.by = "ident", label = TRUE, label.size = 4) +
      NoAxes() +
      NoLegend() +
      ggtitle(NULL)
  })

  # ---- Downloads ----

  # Téléchargement de l'objet Seurat sous-ensemble
  output$download_subset_merge <- createDownloadHandler(
    reactive_data = subset_seurat_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "SubsetSeurat") 
    }),
    data_name = "subset",
    download_type = "seurat",
    show_modal = TRUE
  )

  # Download analysis parameters as CSV
  output$download_params_multi <- downloadHandler(
    filename = function() {
      paste0("multi_datasets_analysis_parameters_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      # Merge integration method into stats at export time
      cs <- qc_stats_multi()
      if (!is.null(cs)) {
        cs$integration_method <- integration_method_used()
        cs$harmony_vars       <- harmony_vars_used()
      }
      # Session parameters + one row per stored analysis (DE tables, condition comparisons, enrichment runs)
      params_df <- rbind(
        collect_analysis_params(
          module     = "multiple",
          qc_stats   = cs,
          seurat_obj = multiple_datasets_object()
        ),
        collect_stored_analyses(gene_table_storage())
      )
      write.csv(params_df, file, row.names = FALSE)
    }
  )

  output$download_report_multi <- downloadHandler(
    filename = function() {
      paste0("multi_dataset_report_", format(Sys.Date(), "%Y%m%d"), ".pdf")
    },
    content = function(file) {
      # Merge integration method into stats at export time
      cs <- qc_stats_multi()
      if (!is.null(cs)) {
        cs$integration_method <- integration_method_used()
        cs$harmony_vars       <- harmony_vars_used()
      }
      # Session parameters + one row per stored analysis (DE tables, condition comparisons, enrichment runs)
      params_df <- rbind(
        collect_analysis_params(
          module     = "multiple",
          qc_stats   = cs,
          seurat_obj = multiple_datasets_object()
        ),
        collect_stored_analyses(gene_table_storage())
      )
      result <- generate_pdf_report(params_df, file, module = "multiple")
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
  pending_table_deletion_merge <- reactiveVal(character(0))

  # ---- Observers and computations ----

  # Seurat object + session state, written by the "Save analysis" button of the Save & Session tab
  session_object_merge <- reactive({
    req(multiple_datasets_object())
    pack_session_state(multiple_datasets_object(), gene_tables = gene_table_storage(),
                       qc_stats = qc_stats_multi(), module = "multiple")
  })

  # Stored result tables, listed by their human-readable description
  session_tables_view_merge <- reactive({
    summarize_gene_tables(gene_table_storage())
  })

  observeEvent(input$delete_session_tables_merge, run_safely("delete session tables merge", {
    selected <- input$session_tables_merge_rows_selected
    if (length(selected) == 0) {
      showNotification("Select at least one table to delete.", type = "warning")
      return()
    }
    view <- session_tables_view_merge()
    pending_table_deletion_merge(view$ID[selected])
    showModal(modalDialog(
      title = "Delete stored tables",
      tags$p(sprintf("Delete %d table(s)? Venn diagrams and volcano plots will no longer be able to use them.", length(selected))),
      tags$ul(lapply(view$Table[selected], tags$li)),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_delete_session_tables_merge", "Delete", class = "btn-danger")
      )
    ))
  }))

  observeEvent(input$confirm_delete_session_tables_merge, run_safely("confirm delete session tables merge", {
    ids <- pending_table_deletion_merge()
    tables <- gene_table_storage()
    tables[intersect(ids, names(tables))] <- NULL
    gene_table_storage(tables)
    pending_table_deletion_merge(character(0))
    removeModal()
    showNotification(sprintf("%d table(s) deleted", length(ids)), type = "message")
  }))  

  # ---- UI and text outputs ----

  output$session_size_info_merge <- renderText({
    tables <- gene_table_storage()
    sprintf("%d stored table(s), %s in total. Plots are never stored.",
            length(tables), format(object.size(tables), units = "auto"))
  })

  # ---- Tables ----

  # What the object contains
  output$session_object_summary_merge <- renderTable({
    req(multiple_datasets_object())
    summarize_seurat_object(multiple_datasets_object())
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$session_tables_merge <- DT::renderDT({
    DT::datatable(
      session_tables_view_merge(),
      rownames  = FALSE,
      selection = "multiple",
      options   = list(dom = "tip", pageLength = 10,
                       # ID column is hidden: used only to delete the right entries
                       columnDefs = list(list(visible = FALSE, targets = 0)))
    )
  })

  # ---- Downloads ----

  output$download_session_merge <- createDownloadHandler(
    reactive_data = session_object_merge,
    object_name_reactive = reactive({
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets")
    }),
    data_name = "CellHub_session",
    download_type = "seurat",
    show_modal = TRUE
  )
}
