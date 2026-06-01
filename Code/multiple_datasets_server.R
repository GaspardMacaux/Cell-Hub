
############################## Script Multiple Datasets Analysis ##############################

# multiple_datasets_server.R

multiple_datasets_server <- function(input, output, session) {  
  
  # Dans multiple_datasets_server_newarch.R - Section Loading Data
  
  ############################## Loading Data ##############################
  
  # Variables for the merge dataset part
  multiple_datasets_object <- reactiveVal()
  seurat_objects <- reactiveValues()
  data_loaded <- reactiveValues()
  merged_gene_tables <- reactiveValues()
  rv_metadata <- reactiveValues(num_fields = 1)
  registered_observers <- new.env(parent = emptyenv())  # tracks which slot observers exist
  
  
  # Analysis parameter tracking for CSV export
  qc_stats_multi          <- reactiveVal(NULL)
  integration_method_used <- reactiveVal("not_run")
  harmony_vars_used       <- reactiveVal(NULL)
  
  shinyjs::disable("add_field")
  shinyjs::disable("add_metadata")
  
  output$datasets_loaded <- reactive({
    !is.null(multiple_datasets_object())
  })
  
  # Initialize metadata split column with default
  observe({
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
  })
  
  # og system pour le panneau permanent
  log_system <- createLogSystem()
  
  output$merge_loading_logs <- renderPrint({
    cat(log_system$log_value())
  })
  
  output$merge_preintegrated_logs <- renderPrint({
    cat(log_system$log_value())
  })
  
  # log system pour le modal d'intégration
  modal_log_system <- createLogSystem()
  
  output$integration_modal_logs <- renderPrint({
    cat(modal_log_system$log_value())
  })
  # Navigate to load datasets tab from metadata management
  observeEvent(input$go_to_load_merge, {
    updateTabItems(session, "tabs", "load_datasets_merge")
  })
  
  # Get cluster column - prioritize ClusterIdents over seurat_clusters over Idents
  get_cluster_column <- function(seurat_obj) {
    if ("ClusterIdents" %in% colnames(seurat_obj@meta.data)) {
      return("ClusterIdents")
    } else if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      return("seurat_clusters")
    } else {
      return("ident")  # Use active Idents as fallback
    }
  }
  
  # Get cluster values from appropriate column
  get_cluster_values <- function(seurat_obj, cluster_col = NULL) {
    if (is.null(cluster_col)) {
      cluster_col <- get_cluster_column(seurat_obj)
    }
    
    if (cluster_col == "ident") {
      return(Idents(seurat_obj))
    } else {
      return(seurat_obj@meta.data[[cluster_col]])
    }
  }
  
  
  
  ############################## Force Memory Cleanup ##############################
  
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
  
  
  
  
  outputOptions(output, 'datasets_loaded', suspendWhenHidden = FALSE)
  
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
      # Memory optimization if requested
      if (isTRUE(input$optimize_memory_preintegrated == "slim")) {
        log_system$add_log("Applying memory optimization...", "info")
        loaded_seurat <- slimSeuratObject(loaded_seurat, log_function = log_system$add_log)
      }
      logMemoryUsage(log_system$add_log, loaded_seurat, "Pre-integrated object")
      multiple_datasets_object(loaded_seurat)
      updateUIElements()
      shinyjs::enable("add_field")
      shinyjs::enable("add_metadata")
      showNotification("Pre-integrated object loaded successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading object:", conditionMessage(e)[1]), type = "error")
    })
  })
  
  # Processing individual datasets using modular functions
  observe_file_input <- function(index) {
    key <- paste0("obs_", index)
    if (isTRUE(registered_observers[[key]])) return()
    registered_observers[[key]] <- TRUE
    observeEvent(input[[paste0("merge", index)]], {
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
          # Validate file type: ZIP or H5AD for raw types, RDS or H5AD for processed
          if (dataset_type %in% c("snRNA_merge", "multiome_merge")) {
            if (!file_extension %in% c("zip", "h5ad")) {
              log_system$add_log("Error: Expected ZIP or H5AD file for raw data", "error")
              showNotification(paste("Dataset", index, ": Please upload a ZIP or H5AD file"), type = "error")
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
          if (isTRUE(input$optimize_memory_merge == "slim")) {
            log_system$add_log(paste("Optimizing memory for dataset", index, "..."), "info")
            seurat_object <- slimSeuratObject(seurat_object, log_function = log_system$add_log)
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
    }, ignoreInit = TRUE)
  }
  
  # Apply dataset names from UI text inputs to already-processed Seurat objects.
  # Called at integration time, not at upload time, so name changes by the user are captured.
  applyDatasetNamesFromInputs <- function(seurat_objects_rv, num_datasets, input) {
    for (i in seq_len(num_datasets)) {
      obj_key <- paste0("seurat_object", i)
      if (!is.null(seurat_objects_rv[[obj_key]])) {
        dataset_name <- input[[paste0("dataset_name", i)]]
        if (is.null(dataset_name) || trimws(dataset_name) == "") {
          dataset_name <- paste0("Dataset_", i)
        }
        # Overwrite integration_group with the name from the UI
        seurat_objects_rv[[obj_key]]@meta.data$integration_group <- dataset_name
        message(paste("Dataset", i, "-> integration_group:", dataset_name))
      }
    }
    return(seurat_objects_rv)
  }
  
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
      integration_method_used("fastmnn")
      logMemoryUsage(log_system$add_log, integrated_object, "FastMNN integrated object")
      
      shinyjs::enable("add_field")
      shinyjs::enable("add_metadata")
      
      removeModal()
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
  
  
  
  # Modal for dataset loading
  observeEvent(input$open_file_input_modal, {
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
                                  choices = c("LogNormalize" = "LogNormalize",
                                              "SCTransform"  = "SCTransform",
                                              "CLR"          = "CLR",
                                              "RC"           = "RC"),
                                  selected = "LogNormalize"))
          ),
          fluidRow(
            column(12,
                   tags$div(
                     style = "margin-top: 10px; background-color: white; padding: 10px; border-radius: 6px; border-left: 4px solid #11998e;",
                     radioButtons("optimize_memory_merge", NULL,
                                  choices = list(
                                    "Optimize memory (recommended)" = "slim",
                                    "Keep full object" = "full"
                                  ),
                                  selected = "slim", inline = TRUE),
                     tags$p(icon("info-circle"),
                            " Drops unused assays and scale.data to reduce RAM usage.",
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
            style = "display: flex; gap: 10px; justify-content: space-between;",
            
            # Standard Integration
            div(
              style = "flex: 1; background-color: #e3f2fd; padding: 12px; border-radius: 6px; 
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
                             style = "width: 100%; font-size: 11px; font-weight: bold;")
              )
            ),
            
            # Harmony Integration
            div(
              style = "flex: 1; background-color: #e8f5e9; padding: 12px; border-radius: 6px; 
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
                actionButton('integrate_harmony', 'Run Harmony Integration',
                             class = 'btn-success', icon = icon("play"),
                             disabled = TRUE,
                             style = "width: 100%; font-size: 11px; font-weight: bold;")
              )
            ),
            
            # fastMNN Integration
            div(
              style = "flex: 1; background-color: #ede7f6; padding: 12px; border-radius: 6px; 
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
                  style = "width: 100%; font-size: 11px; font-weight: bold;"
                )
              )
            ),
            # Simple Merge
            div(
              style = "flex: 1; background-color: #fff3e0; padding: 12px; border-radius: 6px; 
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
                             style = "width: 100%; font-size: 11px; font-weight: bold;")
              )
            )
          ),
          
          tags$div(
            style = "margin-top: 10px;",
            checkboxInput("preserve_annotations",
                          "Preserve original cluster names from individual datasets",
                          value = TRUE)
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
    observeEvent(input$num_datasets, {
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
                                accept = c('.rds', '.zip', '.h5ad', 'application/x-hdf5')))
          )
        })
      })
    })
  })
  
  # Setup file input observers
  observeEvent(input$num_datasets, {
    lapply(1:input$num_datasets, function(i) {
      observe_file_input(i)
    })
  })
  
  # Enable/disable all three buttons based on file status
  observe({
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
    
    # Enable all three buttons
    shinyjs::toggleState("integrate", condition = should_enable)  
    shinyjs::toggleState("integrate_harmony", condition = should_enable)  
    shinyjs::toggleState("simple_merge", condition = should_enable)  
    shinyjs::toggleState("integrate_fastmnn", condition = should_enable)
  })
  
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
          tags$p("Correcting for: dataset variable"),
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
      logMemoryUsage(log_system$add_log, integrated_object, "Harmony integrated object")
      
      # Enable UI controls
      shinyjs::enable("add_field")
      shinyjs::enable("add_metadata")
      
      removeModal()
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
  
  observe({
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
  })
  
  
  # Metadata management using modular functions
  observeEvent(input$add_field, {
    rv_metadata$num_fields <- rv_metadata$num_fields + 1
  })
  
  
  
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
  
  # Download handler for Seurat object with metadata
  output$download_seurat_with_metadata <- createDownloadHandler(
    reactive_data = multiple_datasets_object,
    object_name_reactive = reactive({
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets")
    }),
    data_name = "with_metadata",
    download_type = "seurat",
    show_modal = TRUE
  )
  
  
  # Function to update UI elements after data loading
  updateUIElements <- function() {
    req(multiple_datasets_object())
    tryCatch({
      # Update group_by choices
      updateSelectInput(session, "group_by_select", choices = unique(multiple_datasets_object()@meta.data$dataset))
      
      # Update gene choices using the reactive (more robust)
      gene_list <- reactive_gene_list_merge()
      updatePickerInput(session, "geneInput_merge", choices = gene_list)
      
      message(paste("UI elements updated successfully - loaded", length(gene_list) - 1, "genes"))
    }, error = function(e) {
      message(paste("Error updating UI elements:", e$message))
    })
  }
  
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
      unique_genes <- rownames(LayerData(obj, assay = assay_to_use, layer = 'counts'))
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
  # Function to update UI elements after data loading
  updateUIElements <- function() {
    updateSelectInput(session, "group_by_select", choices = unique(multiple_datasets_object()@meta.data$dataset))
    updatePickerInput(session, "geneInput_merge", choices = reactive_gene_list_merge())
  }
  
  
  ############################## Scaling and PCA reduction ##############################
  
  
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
        seurat_object_temp <- RunPCA(
          seurat_object_temp, 
          features = VariableFeatures(seurat_object_temp),
          npcs = 50,
          verbose = FALSE
        )
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
        all_genes <- rownames(seurat_object_temp)
        seurat_object_temp <- FindVariableFeatures(
          seurat_object_temp, 
          selection.method = "vst", 
          nfeatures = 3000
        )
        seurat_object_temp <- ScaleData(
          seurat_object_temp, 
          features = all_genes
        )
        seurat_object_temp <- RunPCA(
          seurat_object_temp, 
          npcs = 50
        )
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
  
  
  # Update Harmony variables choices
  observe({
    req(multiple_datasets_object())
    meta_cols <- colnames(multiple_datasets_object()@meta.data)
    updateSelectInput(session, "harmony_vars",
                      choices = meta_cols,
                      selected = "dataset")
  })
  
  
  
  # Observer forfind Clusters and display UMAP
  clustering_plot_merge <- reactiveVal()
  
  
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
      
      clustering_plot <- DimPlot(
        multiple_datasets_object(), 
        group.by = if("integration_group" %in% colnames(multiple_datasets_object()@meta.data)) {
          "integration_group"
        } else {
          "dataset"
        }
      ) + ggtitle(NULL)
      clustering_plot_merge(clustering_plot)
      
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
      
      graph_name <- if ("integrated.cca_snn" %in% available_graphs) {
        "integrated.cca_snn"
      } else if ("integrated.rpca_snn" %in% available_graphs) {
        "integrated.rpca_snn"
      } else if ("harmony_snn" %in% available_graphs) {
        "harmony_snn"
      } else {
        available_graphs[1]
      }
      
      message(paste("Using graph:", graph_name))
      metadata_backup <- seurat_integrated_temp@meta.data
      
      seurat_integrated_temp <- findClusters_reproducible(
        object     = seurat_integrated_temp,
        resolution = input$resolution_step2,
        algorithm  = as.integer(input$algorithm_select),
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
      
      plot <- DimPlot(
        multiple_datasets_object(),
        reduction = "umap",
        group.by  = "seurat_clusters",
        label     = TRUE,
        repel     = TRUE
      ) + ggtitle(NULL)
      
      if (input$remove_axes_umap_merge)   { plot <- plot + NoAxes() }
      if (input$remove_legend_umap_merge) { plot <- plot + NoLegend() }
      
      clustering_plot_merge(plot)
      
      showNotification(paste0("Clustering completed! Found ", n_clusters, " clusters."), type = "message")
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Error: ", e$message), type = "error")
      message(paste("Full error:", e))
    })
  })
  
  # Render plot with dark mode option
  output$UMAPPlot_cluster_merge <- renderPlot({
    req(clustering_plot_merge())
    plot_obj <- clustering_plot_merge()
    
    if(input$dark_mode_umap_cluster_merge) {
      plot_obj <- plot_obj + 
        theme(
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          text = element_text(color = "white"),
          panel.grid = element_blank()
        )
      
      # Make labels white if present
      if(length(plot_obj$layers) > 0) {
        for(i in seq_along(plot_obj$layers)) {
          if(inherits(plot_obj$layers[[i]]$geom, c("GeomText", "GeomTextRepel"))) {
            plot_obj$layers[[i]]$aes_params$colour <- "white"
          }
        }
      }
    }
    
    plot_obj
  })
  
  # Download handler - FIXED (format + dark mode)
  output$downloadUMAP_merge <- createDownloadHandler(
    reactive_data = reactive({
      # Get base plot
      plot <- clustering_plot_merge()
      req(plot)
      
      # Apply dark mode if enabled
      if(input$dark_mode_umap_cluster_merge) {
        plot <- plot + 
          theme(
            panel.background = element_rect(fill = "black"),
            plot.background = element_rect(fill = "black"),
            text = element_text(color = "white"),
            panel.grid = element_blank()
          )
        
        # Make labels white if present
        if(length(plot$layers) > 0) {
          for(i in seq_along(plot$layers)) {
            if(inherits(plot$layers[[i]]$geom, c("GeomText", "GeomTextRepel"))) {
              plot$layers[[i]]$aes_params$colour <- "white"
            }
          }
        }
      }
      
      return(plot)
    }),
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
  
  
  # Observer: When genes are selected in pickerInput, update all text inputs
  observeEvent(input$geneInput_merge, {
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
  })
  # OBSERVER 1: Update group_by_select choices (main axis of plots - X axis)
  observeEvent(multiple_datasets_object(), {
    seurat_object <- multiple_datasets_object()
    req(seurat_object)
    
    # Exclude technical columns
    excluded_patterns <- c(
      "percent.mt", "nCount_ATAC", "nFeature_ATAC", "nFeature_RNA", "nCount_RNA", 
      "nFeature_Spatial", "nCount_Spatial", "^RNA_snn_", "^RNA_nn_", 
      "^Spatial_snn_", "^Spatial_nn_", "original_clusters", "^pANN", "^DF", 
      "dataset_origin", "original_seurat_clusters", "original_ClusterIdents",
      "cluster_name_only", "source_format", "^integrated", "^integrated_snn",
      "^merge_method$", "^integration_method$", "^original_clusterIdents$", "^cluster_names$"
    )
    pattern <- paste(excluded_patterns, collapse = "|")
    
    metadata_fields <- colnames(seurat_object@meta.data)
    display_fields <- metadata_fields[!grepl(pattern, metadata_fields)]
    
    if (length(display_fields) == 0) {
      message("No valid metadata fields found for group_by")
      return()
    }
    
    # Determine default selection - use isolate to avoid triggering other observers
    current_selection <- isolate(input$group_by_select)
    
    if (!is.null(current_selection) && current_selection != "" && current_selection %in% display_fields) {
      default_choice <- current_selection
    } else if ("ClusterIdents" %in% display_fields) {
      default_choice <- "ClusterIdents"
    } else if ("seurat_clusters" %in% display_fields) {
      default_choice <- "seurat_clusters"
    } else {
      default_choice <- ""
    }
    
    # Add empty option at the beginning
    choices_with_empty <- c("Select grouping..." = "", display_fields)
    
    updateSelectizeInput(session, "group_by_select",
                         choices = choices_with_empty,
                         selected = default_choice,
                         server = TRUE)
    
    message(paste("Updated group_by_select. Available:", paste(head(display_fields, 5), collapse=", ")))
    
  }, priority = 10, ignoreNULL = FALSE)
  
  
  # OBSERVER 2: Update metadata_to_compare choices (for split.by visual separation)
  observeEvent(list(multiple_datasets_object(), input$group_by_select), {
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
    
  }, priority = 9, ignoreNULL = FALSE)
  
  
  # OBSERVER 3: Update cluster_order choices based on group_by_select values
  observeEvent(input$group_by_select, {
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
    
  }, priority = 8, ignoreNULL = FALSE)
  
  
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
  
  feature_plot_merge <- reactiveVal()
  # Generate feature plot
  observeEvent(input$show_feature_plot_merge, {
    tryCatch({
      req(input$gene_list_feature_merge, multiple_datasets_object())
      
      genes <- unique(trimws(unlist(strsplit(input$gene_list_feature_merge, ","))))
      seurat_object <- multiple_datasets_object()
      req(seurat_object)
      
      # Set assay
      DefaultAssay(seurat_object) <- input$viz_assay_merge
      
      # Validate genes
      present_genes <- genes[genes %in% rownames(LayerData(seurat_object, assay = "RNA", layer = "counts"))]
      missing_genes <- setdiff(genes, present_genes)
      
      if(length(missing_genes) > 0) {
        showNotification(paste("Genes not found:", paste(missing_genes, collapse = ", ")), 
                         type = "warning", duration = 5)
      }
      
      if(length(present_genes) == 0) {
        showNotification("No valid genes found", type = "error")
        return()
      }
      
      message(paste("Plotting genes:", paste(present_genes, collapse = ", ")))
      
      # Apply filters if specified
      if(!is.null(input$metadata_to_compare) && input$metadata_to_compare != "") {
        split_col <- input$metadata_to_compare
        
        if(split_col %in% colnames(seurat_object@meta.data)) {
          # Get selected values
          selected_vals <- input$split_values_filter_feature
          
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
  
  output$FeaturePlot2 <- renderPlot({
    req(feature_plot_merge())
    plot_obj <- feature_plot_merge()
    
    if(inherits(plot_obj, "gg") || inherits(plot_obj, "patchwork")) {
      if(!is.null(input$dark_mode_feature_plot_merge) && isTRUE(input$dark_mode_feature_plot_merge)) {
        plot_obj <- plot_obj & 
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
      plot_obj
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
  
  dot_plot_merge <- reactiveVal()
  
  # DotPlot observer - using modular function
  observeEvent(input$runDotPlot, {
    tryCatch({
      req(input$gene_list_dot_merge, multiple_datasets_object())
      
      # Check that group_by is selected
      if (is.null(input$group_by_select) || input$group_by_select == "") {
        showNotification("Please select a grouping column (Group By) first", type = "warning", duration = 5)
        return()
      }
      
      genes <- unique(trimws(strsplit(input$gene_list_dot_merge, ",")[[1]]))
      seurat_object <- multiple_datasets_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay_merge
      
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
      
      # STEP 3: Filter by group_by values (cluster_order) if specified
      if (!is.null(input$cluster_order_dotplot_merge) && length(input$cluster_order_dotplot_merge) > 0) {
        
        # Determine which cells to keep based on group_by column
        if (group_by == "seurat_clusters") {
          cells_to_keep <- as.character(Idents(seurat_object)) %in% input$cluster_order_dotplot_merge
        } else {
          cells_to_keep <- seurat_object@meta.data[[group_by]] %in% input$cluster_order_dotplot_merge
        }
        
        # Subset the object
        seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
        message(paste("Filtered to", ncol(seurat_object), "cells by group_by values"))
        
        # Set factor levels for ordering
        if (group_by == "seurat_clusters") {
          Idents(seurat_object) <- factor(Idents(seurat_object), levels = input$cluster_order_dotplot_merge)
        } else {
          seurat_object@meta.data[[group_by]] <- factor(seurat_object@meta.data[[group_by]], levels = input$cluster_order_dotplot_merge)
        }
      }
      
      # Check we still have cells
      if (ncol(seurat_object) == 0) {
        showNotification("No cells found with selected parameters", type = "error")
        return()
      }
      
      message(paste("Final cells:", ncol(seurat_object)))
      
      # STEP 4: Create DotPlot using modular function
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
        rotate_axis = FALSE  # No rotation for multiple datasets
      )
      
      message(paste("DotPlot created - palette:", palette_name, "| dot.scale:", dot_scale, 
                    "| split.by:", split_by_column %||% "none"))
      
      # STEP 5: Apply legend customization
      if (input$add_nolegend_dot_merge) {
        plot <- plot + guides(size = "none")
        message("Removed size legend (kept color gradient legend)")
      } else {
        message("All legends kept - showing % expressing (size) and expression gradient (color)")
      }
      
      # STEP 6: Apply axes customization
      if (input$add_noaxes_dot_merge) {
        plot <- plot + NoAxes()
        message("Removed axes")
      }
      
      # STEP 7: Invert axes if requested
      if (input$invert_axes) {
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
  
  # Render DotPlot
  output$DotPlot2 <- renderPlot({
    req(dot_plot_merge())
    dot_plot_merge()
  })
  ############################## VIOLIN PLOT ##############################
  
  vln_plot_merge <- reactiveVal()
  
  observeEvent(input$runVlnPlot, {
    tryCatch({
      req(input$gene_list_vln_merge, multiple_datasets_object())
      
      # Check that group_by is selected
      if (is.null(input$group_by_select) || input$group_by_select == "") {
        showNotification("Please select a grouping column (Group By) first", type = "warning", duration = 5)
        return()
      }
      
      genes <- unique(trimws(strsplit(input$gene_list_vln_merge, ",")[[1]]))
      seurat_object <- multiple_datasets_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay_merge
      
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
        
        # Subset the object
        seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
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
  
  # Render VlnPlot
  output$VlnPlot2 <- renderPlot({
    req(vln_plot_merge())
    vln_plot_merge()
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
  
  ridge_plot_merge <- reactiveVal()
  
  observeEvent(input$runRidgePlot, {
    tryCatch({
      req(input$gene_list_ridge_merge, multiple_datasets_object())
      
      # Check that group_by is selected
      if (is.null(input$group_by_select) || input$group_by_select == "") {
        showNotification("Please select a grouping column (Group By) first", type = "warning", duration = 5)
        return()
      }
      
      genes <- unique(trimws(strsplit(input$gene_list_ridge_merge, ",")[[1]]))
      seurat_object <- multiple_datasets_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay_merge
      
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
  
  output$Ridge_plot_merge <- renderPlot({
    req(ridge_plot_merge())
    ridge_plot_merge()
  })
  
  
  # Download handler for VlnPlot Merge - FIXED
  output$downloadVlnPlotMerge <- createDownloadHandler(
    reactive_data = vln_plot_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets") 
    }),
    data_name = reactive({ 
      paste0("VlnPlot_", paste(trimws(strsplit(input$gene_list_vln_merge, ",")[[1]]), collapse="_")) 
    }),
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$plot_format_merge }),  # ✅ Now reactive
      width = 15, 
      height = reactive({ (attr(vln_plot_merge(), "n_rows") %||% 1) * 8 }), 
      dpi = reactive({ input$dpi_input_merge })  # ✅ Now reactive
    )
  )
  
  # Download handler for FeaturePlot Merge - FIXED
  output$downloadFeaturePlotMerge <- createDownloadHandler(
    reactive_data = feature_plot_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets") 
    }),
    data_name = reactive({ 
      paste0("FeaturePlot_", paste(trimws(strsplit(input$gene_list_feature_merge, ",")[[1]]), collapse="_")) 
    }),
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$plot_format_merge }),  # ✅ Now reactive
      width = reactive({ min(12 * length(trimws(strsplit(input$gene_list_feature_merge, ",")[[1]])), 24) }), 
      height = reactive({ (attr(feature_plot_merge(), "n_rows") %||% 1) * 8 }), 
      dpi = reactive({ input$dpi_input_merge })  # ✅ Now reactive
    )
  )
  
  # Download handler for DotPlot Merge - FIXED
  output$downloadDotPlotMerge <- createDownloadHandler(
    reactive_data = dot_plot_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets") 
    }),
    data_name = reactive({ 
      paste0("DotPlot_", paste(trimws(strsplit(input$gene_list_dot_merge, ",")[[1]]), collapse="_")) 
    }),
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$plot_format_merge }),  # ✅ Now reactive
      width = reactive({ min(12 + length(trimws(strsplit(input$gene_list_dot_merge, ",")[[1]])) * 0.5, 24) }), 
      height = 8, 
      dpi = reactive({ input$dpi_input_merge })  # ✅ Now reactive
    )
  )
  
  # Download handler for RidgePlot Merge - FIXED
  output$downloadRidgePlotMerge <- createDownloadHandler(
    reactive_data = ridge_plot_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets") 
    }),
    data_name = reactive({ 
      paste0("RidgePlot_", paste(trimws(strsplit(input$gene_list_ridge_merge, ",")[[1]]), collapse="_")) 
    }),
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$plot_format_merge }),  # ✅ Now reactive
      width = reactive({ min(10 + length(trimws(strsplit(input$gene_list_ridge_merge, ",")[[1]])) * 0.5, 20) }), 
      height = 8, 
      dpi = reactive({ input$dpi_input_merge })  # ✅ Now reactive
    )
  )
  ################# Number of nuclei expressing a gene ###################
  number_of_nuclei_merge <- reactiveVal(NULL)
  
  
  # Observer pour l'analyse d'expression dans les datasets intégrés
  observeEvent(input$analyze_btn_genes_expression_merge, {
    tryCatch({
      req(input$gene_list_genes_expression_merge, multiple_datasets_object())
      
      # CRITICAL FIX: Parse the gene list string into a vector
      gene_list <- input$gene_list_genes_expression_merge
      
      # Split by comma and trim whitespace
      genes_vector <- trimws(unlist(strsplit(gene_list, ",")))
      
      # Remove empty strings
      genes_vector <- genes_vector[genes_vector != ""]
      
      if (length(genes_vector) == 0) {
        showNotification("Please enter valid gene names", type = "error")
        return()
      }
      
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
  
  # Pour le téléchargement des données
  output$download_genes_number_expression_merge <- createDownloadHandler(
    reactive_data = number_of_nuclei_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "IntegratedExpression") 
    }),
    data_name = "integrated_expression_analysis",
    download_type = "csv"
  )
  
  
  output$save_seurat_merge_2 <- createDownloadHandler(
    reactive_data = multiple_datasets_object,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "IntegratedSeurat") }),
    data_name = "integrated_object",
    download_type = "seurat",
    show_modal = TRUE
  )
  
  
  ############################## Heatmap and dual expression multi dataset ##############################
  
  heatmap_plot_multidataset <- reactiveVal()
  heatmap_n_genes_multi <- reactiveVal(0)  # Track number of genes for dynamic height
  
  # Update selectors when data changes
  observeEvent(multiple_datasets_object(), {
    req(multiple_datasets_object())
    seurat_obj <- multiple_datasets_object()
    
    # Update clusters
    clusters <- levels(Idents(seurat_obj))
    updateSelectizeInput(session, "clusters_heatmap_multi", 
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
    
    updateSelectInput(session, "split_by_heatmap_multi", 
                      choices = c("None" = "None", filtered_cols),
                      selected = "None")
  })
  
  # Update available values when split_by selection changes
  observeEvent(input$split_by_heatmap_multi, {
    req(multiple_datasets_object())
    seurat_obj <- multiple_datasets_object()
    
    split_by <- input$split_by_heatmap_multi
    
    if (is.null(split_by) || split_by == "None" || split_by == "") {
      # No split selected - clear filter choices
      updateSelectizeInput(session, "datasets_heatmap_multi",
                           choices = character(0),
                           selected = character(0),
                           server = TRUE)
    } else {
      # Update with unique values from selected split_by column
      if (split_by %in% colnames(seurat_obj@meta.data)) {
        available_values <- sort(unique(as.character(seurat_obj@meta.data[[split_by]])))
        available_values <- available_values[!is.na(available_values)]
        
        updateSelectizeInput(session, "datasets_heatmap_multi",
                             choices = available_values,
                             selected = character(0),
                             server = TRUE)
        
        message(paste("Updated filter choices with", length(available_values), 
                      "unique values from", split_by, "column"))
      }
    }
  })
  
  # Store top markers table for CSV export (NULL when using manual gene input)
  heatmap_markers_data_multi <- reactiveVal(NULL)
  
  # Generate heatmap - with split/group/merge logic
  observeEvent(input$generateHeatmapMulti, {
    showModal(modalDialog(
      title = "Generating heatmap...",
      "Please wait while generating average expression heatmap",
      easyClose = FALSE,
      footer = NULL
    ))
    
    req(multiple_datasets_object())
    seurat_object <- multiple_datasets_object()
    
    tryCatch({
      # Validate assay
      selected_assay <- input$assay_select_heatmap
      if (is.null(selected_assay)) selected_assay <- "RNA"
      if (!(selected_assay %in% names(seurat_object@assays))) selected_assay <- "RNA"
      DefaultAssay(seurat_object) <- selected_assay
      
      # Get clusters selected
      clusters_selected <- input$clusters_heatmap_multi
      if (is.null(clusters_selected) || length(clusters_selected) == 0) {
        clusters_selected <- levels(Idents(seurat_object))
      }
      
      # Get genes
      if (input$use_top10_genes_merge) {
        n_genes <- ifelse(is.null(input$n_top_genes_multi), 10, input$n_top_genes_multi)
        
        # get_top_markers now returns a list with genes + markers_table
        markers_result <- get_top_markers(
          seurat_object,
          n_genes = n_genes,
          assay = selected_assay,
          clusters_to_include = clusters_selected
        )
        genes <- markers_result$genes
        
        # Store full markers table for CSV export
        heatmap_markers_data_multi(markers_result$markers_table)
        message(paste("Using top", n_genes, "genes per cluster (for", length(clusters_selected), "selected clusters)"))
        
      } else {
        gene_text <- trimws(input$gene_select_heatmap_multi)
        if (nchar(gene_text) == 0) {
          showNotification("Please enter genes.", type = "error")
          removeModal()
          return()
        }
        gene_result <- parse_and_validate_genes(gene_text, seurat_object, selected_assay)
        genes <- gene_result$valid_genes
        
        # No per-cluster info available in manual mode
        heatmap_markers_data_multi(NULL)
        
        if (length(gene_result$missing_genes) > 0) {
          showNotification(paste(length(gene_result$missing_genes), "gene(s) not found"), 
                           type = "warning", duration = 6)
        }
        
        if (length(genes) == 0) {
          showNotification("No valid genes found in assay", type = "error")
          removeModal()
          return()
        }
      }
      
      # CRITICAL: Store number of genes for dynamic plot height
      heatmap_n_genes_multi(length(genes))
      
      # Get display options
      color_palette <- input$color_palette_heatmap_multi
      if (is.null(color_palette)) color_palette <- "viridis"
      scale_rows <- TRUE
      if (is.null(scale_rows)) scale_rows <- TRUE
      
      # Get split_by
      split_by <- input$split_by_heatmap_multi
      if (is.null(split_by) || split_by == "None" || split_by == "") {
        split_by <- NULL
      }
      
      # Get dataset filter (values to include from split_by column)
      datasets_selected <- input$datasets_heatmap_multi
      if (is.null(datasets_selected) || length(datasets_selected) == 0) {
        datasets_selected <- NULL
      }
      
      # Get grouping preference
      group_by_split <- input$group_by_dataset_heatmap_multi
      if (is.null(group_by_split)) group_by_split <- TRUE
      
      # Get merge clusters preference
      merge_clusters <- input$merge_clusters_heatmap_multi
      if (is.null(merge_clusters)) merge_clusters <- "separate"
      
      # ===== LOGIC WITH SPLIT =====
      if (!is.null(split_by)) {
        # Filter to selected split values if specified
        if (!is.null(datasets_selected) && length(datasets_selected) > 0) {
          cells_to_keep <- seurat_object@meta.data[[split_by]] %in% datasets_selected
          if (sum(cells_to_keep) == 0) {
            showNotification("No cells found for selected values", type = "error")
            removeModal()
            return()
          }
          seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
          message(paste("Filtered to", sum(cells_to_keep), "cells with values:", 
                        paste(datasets_selected, collapse = ", ")))
        }
        
        # Check merge mode
        if (merge_clusters == "merge") {
          # MERGE ALL CLUSTERS - group only by split variable
          message("Merging all clusters by split variable")
          
          # Optionally filter to selected clusters before merging
          if (!is.null(clusters_selected) && length(clusters_selected) > 0) {
            cells_to_keep <- Idents(seurat_object) %in% clusters_selected
            if (sum(cells_to_keep) > 0) {
              seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
              message(paste("Filtered to clusters:", paste(clusters_selected, collapse = ", "), "before merging"))
            }
          }
          
          avg_exp <- AverageExpression(
            seurat_object,
            features = genes,
            assays = selected_assay,
            group.by = split_by
          )
          
          subtitle <- paste("Split by:", split_by, "(all clusters merged)")
          
        } else {
          # KEEP CLUSTERS SEPARATE - create cluster___splitvalue combinations
          # Filter to selected clusters
          if (!is.null(clusters_selected) && length(clusters_selected) > 0) {
            cells_to_keep <- Idents(seurat_object) %in% clusters_selected
            if (sum(cells_to_keep) == 0) {
              showNotification("No cells found for selected clusters", type = "error")
              removeModal()
              return()
            }
            seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
            Idents(seurat_object) <- factor(Idents(seurat_object), levels = clusters_selected)
          } else {
            clusters_selected <- levels(Idents(seurat_object))
          }
          
          # Create combined categories: cluster___splitvalue
          current_idents <- as.character(Idents(seurat_object))
          split_vals <- trimws(as.character(seurat_object@meta.data[[split_by]]))
          seurat_object$temp_heatmap_group <- paste0(current_idents, "___", split_vals)
          
          # Determine order based on grouping preference
          unique_split_vals <- unique(split_vals)
          
          if (group_by_split) {
            # Group by SPLIT VALUE: all clusters within each split value
            valid_combinations <- character()
            for (split_val in sort(unique_split_vals)) {
              for (clust in clusters_selected) {
                combo <- paste0(clust, "___", split_val)
                if (combo %in% seurat_object$temp_heatmap_group) {
                  valid_combinations <- c(valid_combinations, combo)
                }
              }
            }
          } else {
            # Group by CLUSTER: all split values within each cluster
            valid_combinations <- character()
            for (clust in clusters_selected) {
              for (split_val in sort(unique_split_vals)) {
                combo <- paste0(clust, "___", split_val)
                if (combo %in% seurat_object$temp_heatmap_group) {
                  valid_combinations <- c(valid_combinations, combo)
                }
              }
            }
          }
          
          # Filter to valid combinations and set factor levels
          cells_final <- seurat_object$temp_heatmap_group %in% valid_combinations
          seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_final])
          seurat_object$temp_heatmap_group <- factor(seurat_object$temp_heatmap_group, 
                                                     levels = valid_combinations)
          
          message(paste("Final combinations:", paste(head(valid_combinations, 6), collapse = " | ")))
          
          # Compute average expression
          avg_exp <- AverageExpression(
            seurat_object,
            features = genes,
            assays = selected_assay,
            group.by = "temp_heatmap_group"
          )
          
          # Clean up temporary column
          seurat_object$temp_heatmap_group <- NULL
          
          subtitle <- paste("Split by:", split_by, 
                            if (group_by_split) "(grouped by split value)" else "(grouped by cluster)")
        }
        
      } else {
        # NO SPLIT - simple cluster grouping
        if (!is.null(clusters_selected) && length(clusters_selected) > 0) {
          seurat_subset <- subset(seurat_object, idents = clusters_selected)
          Idents(seurat_subset) <- factor(Idents(seurat_subset), levels = clusters_selected)
        } else {
          seurat_subset <- seurat_object
        }
        
        avg_exp <- AverageExpression(
          seurat_subset,
          features = genes,
          assays = selected_assay,
          group.by = "ident"
        )
        
        subtitle <- NULL
      }
      
      # Extract matrix and scale
      exp_matrix <- avg_exp[[selected_assay]]
      exp_matrix <- exp_matrix[genes[genes %in% rownames(exp_matrix)], , drop = FALSE]
      
      if (scale_rows) {
        exp_matrix <- t(scale(t(exp_matrix)))
      }
      
      # Convert to long format
      exp_df <- as.data.frame(exp_matrix)
      exp_df$Gene <- rownames(exp_df)
      exp_long <- tidyr::pivot_longer(exp_df, cols = -Gene, names_to = "Group", values_to = "Expression")
      
      # Factor genes and groups - CRITICAL: ensure unique levels
      gene_levels <- rev(unique(genes[genes %in% unique(exp_long$Gene)]))
      group_levels <- unique(colnames(exp_matrix))
      exp_long$Gene <- factor(exp_long$Gene, levels = gene_levels)
      exp_long$Group <- factor(exp_long$Group, levels = group_levels)
      
      # Create plot
      plot <- ggplot(exp_long, aes(x = Group, y = Gene, fill = Expression)) +
        geom_tile(color = "white", linewidth = 0.5) +
        labs(
          title = paste("Heatmap -", selected_assay, "assay"),
          subtitle = subtitle,
          x = if (is.null(split_by)) "Cluster" else paste("Cluster ×", split_by),
          y = "Gene",
          fill = if (scale_rows) "Scaled\nExpression" else "Average\nExpression"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
          axis.text.y = element_text(size = 10, face = "bold", hjust = 1),
          panel.grid = element_blank(),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
        )
      
      # Apply color palette
      if (color_palette == "viridis") {
        plot <- plot + scale_fill_viridis_c(option = "viridis")
      } else if (color_palette == "magma") {
        plot <- plot + scale_fill_viridis_c(option = "magma")
      } else if (color_palette == "inferno") {
        plot <- plot + scale_fill_viridis_c(option = "inferno")
      } else if (color_palette == "plasma") {
        plot <- plot + scale_fill_viridis_c(option = "plasma")
      } else if (color_palette == "RdYlBu") {
        plot <- plot + scale_fill_gradient2(
          low = "blue", mid = "yellow", high = "red",
          midpoint = if (scale_rows) 0 else median(exp_long$Expression, na.rm = TRUE)
        )
      } else if (color_palette == "RdBu") {
        plot <- plot + scale_fill_gradient2(
          low = "blue", mid = "white", high = "red",
          midpoint = if (scale_rows) 0 else median(exp_long$Expression, na.rm = TRUE)
        )
      }
      
      heatmap_plot_multidataset(plot)
      showNotification(paste0("Heatmap with ", length(genes), " genes"), type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      message(paste("Full error:", e$message))
    })
    
    removeModal()
  })
  
  # Render heatmap with dynamic height
  output$heatmap_plot_multi <- renderPlot({
    req(heatmap_plot_multidataset())
    heatmap_plot_multidataset()
  }, height = function() {
    n_genes <- heatmap_n_genes_multi()
    if (n_genes == 0) return(600)
    max(400, min(20 * n_genes, 2000))
  })
  
  # Download plot handler with dynamic height
  output$download_heatmap_multi <- createDownloadHandler(
    reactive_data = heatmap_plot_multidataset,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MultiDatasets") 
    }),
    data_name = "heatmap",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$heatmap_format }), 
      width = 10, 
      height = reactive({
        n_genes <- heatmap_n_genes_multi()
        if (n_genes == 0) return(8)
        max(6, min(0.15 * n_genes, 50))
      }), 
      dpi = reactive({ input$dpi_heatmap_multi })  
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
  
 
  # Updated selectInput to choose how to group data
  observe({
    updateSelectInput(session, "dataset_select_scatter", choices = reactive_metadata_fields())
  })
  
  
  # Mettez à jour les choix de selectInput pour les gènes/features et les clusters
  observe({
    updatePickerInput(session, "gene_select_scatter_multi1", choices = reactive_gene_list_merge())
    updatePickerInput(session, "gene_select_scatter_multi2", choices = reactive_gene_list_merge())
  })
  
  observe({
    req(multiple_datasets_object())
    clusters <- levels(Idents(multiple_datasets_object()))
    updateSelectInput(session, "cluster_selector_merge_scatter", choices = clusters, selected = clusters)
    updateTextInput(session, "text_clusters_merge_scatter", value = paste(clusters, collapse = ","))
  })
  
  
  observeEvent(input$select_all_clusters_merge_scatter, {
    req(multiple_datasets_object())
    tryCatch({
      if (input$select_all_clusters_merge_scatter) {
        shinyjs::disable("text_clusters_merge_scatter")
      } else {
        shinyjs::enable("text_clusters_merge_scatter")
      }
    }, error = function(e) {
      showNotification(paste("Erreur lors de la sélection des clusters : ", e$message), type = "error")
    })
  })
  
  
  
  ############################## Scatter Plot Multi Dataset ##############################
  
  scatter_plot_multi <- reactiveVal()
  
  # Update color_by and split_by choices when data changes
  observeEvent(multiple_datasets_object(), {
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
  })
  
  # Generate scatter plot
  observeEvent(input$generateScatter_multi, {
    req(multiple_datasets_object())
    
    showModal(modalDialog(
      title = "Processing",
      div(h4("Generating Scatter Plot...", style = "text-align: center;")),
      footer = NULL,
      easyClose = FALSE
    ))
    
    seurat_obj <- multiple_datasets_object()
    
    tryCatch({
      # Get genes
      gene1 <- trimws(input$feature1_select_multi)
      gene2 <- trimws(input$feature2_select_multi)
      
      if (gene1 == "" || gene2 == "") {
        showNotification("Please enter both gene names.", type = "error")
        removeModal()
        return()
      }
      
      # Get assay
      selected_assay <- input$assay_select_scatter_multi
      if (is.null(selected_assay)) selected_assay <- "RNA"
      
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
  })
  
  # Render scatter plot
  output$scatter_plot_multi <- renderPlot({
    req(scatter_plot_multi())
    scatter_plot_multi()
  })
  
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
  # Variable réactive pour stocker le plot actuel
  current_plot_merge <- reactiveVal()
  
  # Observer les changements dans la sélection du type de plot et mettre à jour le plot
  observe({
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
  })
  
  # Afficher le plot sélectionné
  output$selected_plot_display_merge <- renderPlot({
    req(current_plot_merge())
    current_plot_merge()
  })
  
  
  ############################## Assigning cell identity merge ##############################
  
  # Historique très simple
  cluster_history_merge <- reactiveValues(
    states = list(),
    max_history = 10
  )
  
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
  
  # Undo simple
  undo_last_action <- function() {
    if(length(cluster_history_merge$states) <= 1) return(NULL)
    cluster_history_merge$states <- cluster_history_merge$states[-length(cluster_history_merge$states)]
    if(length(cluster_history_merge$states) > 0) {
      return(cluster_history_merge$states[[length(cluster_history_merge$states)]]$clusters)
    }
    return(NULL)
  }
  
  observe({
    if (!is.null(multiple_datasets_object())) {
      updateSelectInput(session, "select_cluster_merge", choices = unique(Idents(multiple_datasets_object())))
      updateSelectInput(session, "select_color_merge", choices = levels(Idents(multiple_datasets_object())))
      get_cluster_colors_merge(multiple_datasets_object())
    }
  })
  
  # Update color picker to reflect the stored color of the selected cluster
  observeEvent(input$select_color_merge, {
    req(multiple_datasets_object(), input$select_color_merge)
    colors <- get_cluster_colors_merge(multiple_datasets_object())
    if (!is.null(colors) && input$select_color_merge %in% names(colors)) {
      updateColourInput(session, "select_cluster_merge_color", value = colors[[input$select_color_merge]])
    }
  }, ignoreInit = TRUE)
  
  
  # Ton observer existant avec juste l'ajout de la sauvegarde
  observeEvent(input$rename_single_cluster_merge_button, {
    req(input$select_cluster_merge, input$rename_single_cluster_merge, multiple_datasets_object())
    
    # SAUVEGARDER AVANT MODIFICATION
    save_cluster_state_merge(multiple_datasets_object())
    
    # Ton code existant...
    updated_seurat <- multiple_datasets_object()
    
    if (input$rename_single_cluster_merge %in% unique(Idents(updated_seurat))) {
      cells_to_merge <- which(Idents(updated_seurat) %in% c(input$select_cluster_merge, input$rename_single_cluster_merge))
      Idents(updated_seurat, cells = cells_to_merge) <- input$rename_single_cluster_merge
      showNotification(paste("Clusters merged under the name:", input$rename_single_cluster_merge), type = "message")
    } else {
      Idents(updated_seurat, cells = which(Idents(updated_seurat) == input$select_cluster_merge)) <- input$rename_single_cluster_merge
      showNotification(paste("Cluster renamed to:", input$rename_single_cluster_merge), type = "message")
    }
    
    new_ident <- setNames(levels(updated_seurat), levels(Idents(updated_seurat)))
    updated_seurat <- RenameIdents(updated_seurat, new_ident)
    
    multiple_datasets_object(updated_seurat)
    get_cluster_colors_merge(updated_seurat)
    updateSelectInput(session, "select_cluster_merge", choices = unique(Idents(updated_seurat)))
    
    updated_seurat$ClusterIdents <- Idents(updated_seurat)
    multiple_datasets_object(updated_seurat)
  })
  
  # Reactive variable for that tab
  cluster_colours_merge <- reactiveVal()
  
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
  
  # Function to apply updated colors to Seurat object and update reactive value
  apply_cluster_colours <- function(seurat_object) {
    cluster_colors <- cluster_colours_merge()
    if (!is.null(cluster_colors)) {
      seurat_object@misc$cluster_colors <- cluster_colors
    }
    return(seurat_object)
  }
  
  observeEvent(input$update_colour_merge_button, {
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
  })
  
  # Compute 3D UMAP when button is clicked
  observeEvent(input$compute_3d_umap_merge, {
    req(multiple_datasets_object())
    
    showModal(modalDialog(
      title = "Computing 3D UMAP",
      "Calculating 3D UMAP coordinates... This may take a moment.",
      easyClose = FALSE,
      footer = NULL
    ))
    
    tryCatch({
      seurat_obj <- multiple_datasets_object()
      
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
      multiple_datasets_object(seurat_obj)
      
      removeModal()
      showNotification("3D UMAP computed successfully!", type = "message")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error computing 3D UMAP:", e$message), type = "error")
      message("Error in 3D UMAP computation: ", e$message)
    })
  })
  
  
  
  # Display final UMAP with 2D/3D support
  output$umap_finale_merge <- renderPlotly({
    req(multiple_datasets_object())
    
    # Don't continue if no UMAP
    if(!"umap" %in% names(multiple_datasets_object()@reductions)) {
      return(NULL)
    }
    
    message("Generating the final UMAP for multiple datasets")
    seurat_obj <- multiple_datasets_object()
    
    # Check if 3D view is requested
    use_3d <- isTRUE(input$umap_3d_toggle_merge)
    
    # Get cluster colors
    cluster_colors <- get_cluster_colors_merge(seurat_obj)
    
    # Generate plot based on dimensionality
    if (use_3d && hasValid3DUMAP(seurat_obj)) {
      # Extract 3D UMAP coordinates
      umap_coords <- as.data.frame(Embeddings(seurat_obj, reduction = "umap3d"))
      colnames(umap_coords) <- c("UMAP_1", "UMAP_2", "UMAP_3")
      umap_coords$cluster <- Idents(seurat_obj)
      
      # Ensure all clusters have colors
      if (!all(levels(as.factor(umap_coords$cluster)) %in% names(cluster_colors))) {
        missing_clusters <- setdiff(levels(as.factor(umap_coords$cluster)), names(cluster_colors))
        default_colors <- scales::hue_pal()(length(missing_clusters))
        names(default_colors) <- missing_clusters
        cluster_colors <- c(cluster_colors, default_colors)
      }
      
      # Calculate cluster centers for labels
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
      
      # Get display options
      hide_grid <- isTRUE(input$umap_3d_hide_grid_merge)
      hide_axes <- isTRUE(input$umap_3d_hide_axes_merge)
      dark_mode <- isTRUE(input$umap_3d_dark_mode_merge)
      
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
          marker = list(size = input$pt_size_merge, opacity = 0.7),
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
          textfont = list(size = input$label_font_size_merge * 4, color = text_color),
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
            text = input$plot_title_merge, 
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
      # 2D UMAP (existing code)
      umap_data <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
      
      message("UMAP column names: ", paste(colnames(umap_data), collapse=", "))
      dim_names <- colnames(umap_data)
      if (length(dim_names) >= 2) {
        colnames(umap_data)[1:2] <- c("UMAP_1", "UMAP_2")
      } else {
        showNotification("UMAP reduction not found or has incorrect dimensions", type = "error")
        return(NULL)
      }
      
      umap_data$cluster <- Idents(seurat_obj)
      
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
        add_trace(data = umap_data, x = ~UMAP_1, y = ~UMAP_2, color = ~cluster, colors = cluster_colors,
                  type = "scatter", mode = "markers", marker = list(size = input$pt_size_merge, opacity = 0.7),
                  hoverinfo = "text", text = ~paste("Cluster:", cluster)) %>%
        layout(title = list(text = input$plot_title_merge, font = list(size = 24)),
               xaxis = list(title = "UMAP 1", zeroline = FALSE),
               yaxis = list(title = "UMAP 2", zeroline = FALSE),
               hovermode = "closest", showlegend = FALSE)
      
      for (i in 1:nrow(centers)) {
        p <- p %>% add_annotations(x = centers$UMAP_1[i], y = centers$UMAP_2[i],
                                   text = as.character(centers$cluster[i]),
                                   showarrow = FALSE,
                                   font = list(size = input$label_font_size_merge * 4, color = "black"))
      }
      
      message("2D UMAP generated for multiple datasets")
      return(p)
    }
  })
  # Observer pour Undo
  observeEvent(input$undo_cluster_merge, {
    req(multiple_datasets_object())
    
    restored_clusters <- undo_last_action()
    
    if(!is.null(restored_clusters)) {
      updated_seurat <- multiple_datasets_object()
      Idents(updated_seurat) <- restored_clusters
      updated_seurat$ClusterIdents <- Idents(updated_seurat)
      multiple_datasets_object(updated_seurat)
      get_cluster_colors_merge(updated_seurat)
      updateSelectInput(session, "select_cluster_merge", choices = unique(Idents(updated_seurat)))
      showNotification("Undid last action", type = "message")
    } else {
      showNotification("Nothing to undo", type = "warning")
    }
  })
  
  
  # Sauvegarder l'état initial
  observe({
    req(multiple_datasets_object())
    if(length(cluster_history_merge$states) == 0) {
      save_cluster_state_merge(multiple_datasets_object(), "Initial state")
    }
  })
  
  
  
  output$save_seurat_merge_3 <- createDownloadHandler(
    reactive_data = multiple_datasets_object,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "IntegratedSeurat") }),
    data_name = "integrated_object",
    download_type = "seurat",
    show_modal = TRUE
  )
  
  
  
  ############################# Differential Gene Expression Analysis  #############################
  
  # Reactive variables for this tab
  diff_genes_compare <- reactiveVal(NULL)          # Stores DE genes for a cluster vs all others
  diff_genes_compare_cluster <- reactiveVal(NULL)  # Stores DE genes for group1 vs group2
  diff_genes_compare_datasets <- reactiveVal(NULL) # Stores DE genes for dataset1 vs dataset2
  filtered_umap_plot <- reactiveVal(NULL)          # Stores the filtered UMAP plot
  gene_table_storage <- reactiveVal(list())        # Stores all DE gene tables for Venn diagrams
  current_gene_lists <- reactiveVal(NULL)          # Stores current lists for Venn diagrams
  
  # Disable download buttons at startup
  shinyjs::disable("download_markers_single_cluster_merge")
  shinyjs::disable("download_markers_multiple_clusters_merge")
  shinyjs::disable("download_diff_dataset_cluster")
  shinyjs::disable("download_venn_diagram")
  
  ############################# UMAP Filtering and Display #############################
  
  # Dropdown menu to filter by Dataset
  output$dataset_filter_ui <- renderUI({
    req(multiple_datasets_object())
    selectInput("dataset_filter",
                label = "Filter by Dataset",
                choices = c(unique(multiple_datasets_object()@meta.data$dataset)),
                selected = unique(multiple_datasets_object()@meta.data$dataset),
                multiple = TRUE)
  })
  
  
  # Update DE assay choices when merged object is loaded
  observe({
    req(multiple_datasets_object())
    available_assays <- names(multiple_datasets_object()@assays)
    default_assay    <- DefaultAssay(multiple_datasets_object())
    updateSelectInput(session, "assay_de_merge", choices = available_assays, selected = default_assay)
  })
  
  
  
  # Display the filtered UMAP plot
  output$filtered_umap_plot <- renderPlot({
    req(multiple_datasets_object(), input$dataset_filter)
    
    # Ne pas continuer si pas d'UMAP
    if(!"umap" %in% names(multiple_datasets_object()@reductions)) {
      return(NULL)
    }
    
    if (!"dataset" %in% colnames(multiple_datasets_object()@meta.data)) {
      stop("Dataset column not found in Seurat object metadata.")
    }
    
    show_labels <- input$show_labels_merge
    label_size <- ifelse(input$bold_labels_merge, 8, 5)
    
    plot <- if ("All Datasets" %in% input$dataset_filter || is.null(input$dataset_filter)) {
      DimPlot(multiple_datasets_object(),
              group.by = "ident",
              label = show_labels,
              label.size = label_size) +
        NoAxes() +
        NoLegend() +
        ggtitle(NULL)
    } else {
      valid_datasets <- input$dataset_filter %in% unique(multiple_datasets_object()@meta.data$dataset)
      subset_seurat <- subset(multiple_datasets_object(), subset = dataset %in% input$dataset_filter[valid_datasets])
      
      DimPlot(subset_seurat,
              group.by = "ident",
              label = show_labels,
              label.size = label_size) +
        NoAxes() +
        NoLegend() +
        ggtitle(NULL)
    }
    
    cluster_colors <- get_cluster_colors_merge(multiple_datasets_object())
    if (!is.null(cluster_colors)) {
      plot <- plot + scale_color_manual(values = cluster_colors)
    }
    
    if(input$dark_mode_filtered_umap) {
      plot <- plot + 
        theme(
          panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          text = element_text(color = "white"),
          panel.grid = element_blank()
        )
      
      if(show_labels) {
        for(i in seq_along(plot$layers)) {
          if(inherits(plot$layers[[i]]$geom, c("GeomText", "GeomTextRepel"))) {
            plot$layers[[i]]$aes_params$colour <- "white"
          }
        }
      }
    }
    
    filtered_umap_plot(plot)
    plot
  })
  
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
  output$save_seurat_merge_4 <- createDownloadHandler(
    reactive_data = multiple_datasets_object,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "IntegratedSeurat") }),
    data_name = "integrated_object",
    download_type = "seurat",
    show_modal = TRUE
  )
  
  
  ############################# Cluster vs All Comparison #############################
  
  # Reactive update of cluster dropdown
  observe({
    req(multiple_datasets_object())
    cluster_choices <- unique(Idents(multiple_datasets_object()))
    updateSelectInput(session, "selected_cluster", choices = cluster_choices, selected = cluster_choices[1])
  })
  
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
      
      markers <- FindMarkers(subset_seurat,
                             ident.1 = input$selected_cluster,
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
      parameters  <- list(min_pct = input$min_pct_merge, logfc_threshold = input$logfc_threshold_merge)
      
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
      df$p_val     <- format_pvalue_robust(df$p_val)
      df$p_val_adj <- format_pvalue_robust(df$p_val_adj)
      df <- df[, c("gene_link", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj")]
      colnames(df) <- c("Gene", "Log2FC", "Pct.1", "Pct.2", "P-value", "Adj. P-value")
      datatable(df, escape = FALSE,
                options = list(pageLength = 10, lengthMenu = c(10, 25, 50, 100, 200, 500),
                               order = list(list(1, "desc")), dom = "Blfrtip", scrollX = TRUE,
                               columnDefs = list(list(className = "dt-center", targets = 1:5),
                                                 list(width = "150px", targets = 0))),
                class = "cell-border stripe", rownames = FALSE)
    }, error = function(e) {
      showNotification(paste0("Error displaying gene table: ", e$message), type = "error")
    })
  })
  
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
  
  output$download_venn_gene_lists <- createDownloadHandler(
    reactive_data = current_gene_lists,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "VennDiagram") }),
    data_name = "gene_lists",
    download_type = "ods"
  )
  
  ############################# Group vs Group Comparison #############################
  
  # UI to select clusters
  output$cluster1_compare_ui <- renderUI({
    req(multiple_datasets_object())
    selectInput("cluster1_compare",
                label = "Select clusters for group 1",
                choices = unique(Idents(multiple_datasets_object())),
                selected = NULL,
                multiple = TRUE)
  })
  
  output$cluster2_compare_ui <- renderUI({
    req(multiple_datasets_object())
    # Exclude clusters already selected in group 1
    remaining_clusters <- setdiff(unique(Idents(multiple_datasets_object())), input$cluster1_compare)
    selectInput("cluster2_compare",
                label = "Select clusters for group 2",
                choices = remaining_clusters,
                selected = NULL,
                multiple = TRUE)
  })
  
  # Observer to update second group choices
  observeEvent(input$cluster1_compare, {
    req(multiple_datasets_object())
    remaining_clusters <- setdiff(unique(Idents(multiple_datasets_object())), input$cluster1_compare)
    updateSelectInput(session, "cluster2_compare", choices = remaining_clusters)
  })
  
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
      new_idents  <- as.character(Idents(seurat_obj))
      new_idents[new_idents %in% input$cluster1_compare] <- "group1"
      new_idents[new_idents %in% input$cluster2_compare] <- "group2"
      Idents(seurat_obj) <- new_idents
      
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
      df$p_val     <- format_pvalue_robust(df$p_val)
      df$p_val_adj <- format_pvalue_robust(df$p_val_adj)
      df <- df[, c("gene_link", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj", "comparison")]
      colnames(df) <- c("Gene", "Log2FC", "Pct.1", "Pct.2", "P-value", "Adj. P-value", "Comparison")
      datatable(df, escape = FALSE,
                options = list(pageLength = 10, lengthMenu = c(10, 25, 50, 100, 200, 500),
                               order = list(list(1, "desc")), dom = "Blfrtip", scrollX = TRUE,
                               columnDefs = list(list(className = "dt-center", targets = 1:6),
                                                 list(width = "150px", targets = 0))),
                class = "cell-border stripe", rownames = FALSE)
    }, error = function(e) {
      showNotification(paste0("Error displaying gene table: ", e$message), type = "error")
    })
  })
  
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
  
  # Observer to update group 2 choices when group 1 changes
  observeEvent(input$dataset1_compare, {
    req(multiple_datasets_object(), input$metadata_column_compare)
    
    metadata_values <- unique(multiple_datasets_object()@meta.data[[input$metadata_column_compare]])
    metadata_values <- sort(metadata_values[!is.na(metadata_values)])
    remaining_values <- setdiff(metadata_values, input$dataset1_compare)
    
    updateSelectInput(session, "dataset2_compare", choices = remaining_values)
  })
  
  # Observer for comparison button - COMPLETE VERSION WITH USER-SELECTED CLUSTER COLUMN
  observeEvent(input$compare_datasets_button, {
    tryCatch({
      showModal(modalDialog(
        title = "Please wait",
        "Comparing groups and clusters...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      req(input$metadata_column_compare,
          input$dataset1_compare,
          input$dataset2_compare,
          input$cluster_column_compare,
          multiple_datasets_object())
      
      if (any(input$dataset1_compare %in% input$dataset2_compare)) {
        showNotification("Groups must not have overlapping values!", type = "error")
        removeModal()
        return()
      }
      
      seurat_obj   <- multiple_datasets_object()
      metadata_col <- input$metadata_column_compare
      cluster_col  <- input$cluster_column_compare
      
      message(paste("Using metadata column:", metadata_col))
      message(paste("Using cluster column:", cluster_col))
      
      if (!metadata_col %in% colnames(seurat_obj@meta.data)) {
        showNotification(paste("Metadata column", metadata_col, "not found!"), type = "error")
        removeModal()
        return()
      }
      
      if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
        showNotification(paste("Cluster column", cluster_col, "not found!"), type = "error")
        removeModal()
        return()
      }
      
      cluster_values <- if (cluster_col == "seurat_clusters") {
        as.character(Idents(object = seurat_obj))
      } else {
        as.character(seurat_obj@meta.data[[cluster_col]])
      }
      
      seurat_obj$temp_grouping <- paste(seurat_obj@meta.data[[metadata_col]],
                                        cluster_values,
                                        sep = "_")
      
      clusters_to_analyze <- if (input$all_clusters) {
        unique(cluster_values[!is.na(cluster_values) & cluster_values != ""])
      } else {
        input$cluster_compare
      }
      
      if (!input$all_clusters && (is.null(clusters_to_analyze) || length(clusters_to_analyze) == 0)) {
        showNotification("Please select at least one cluster or check 'Compare all clusters'", type = "error")
        removeModal()
        return()
      }
      
      group1_combinations <- expand.grid(
        group   = input$dataset1_compare,
        cluster = clusters_to_analyze,
        stringsAsFactors = FALSE
      )
      group2_combinations <- expand.grid(
        group   = input$dataset2_compare,
        cluster = clusters_to_analyze,
        stringsAsFactors = FALSE
      )
      
      group1_ids <- apply(group1_combinations, 1, function(x) paste(x[1], x[2], sep = "_"))
      group2_ids <- apply(group2_combinations, 1, function(x) paste(x[1], x[2], sep = "_"))
      
      available_groups <- unique(seurat_obj$temp_grouping)
      missing_group1   <- setdiff(group1_ids, available_groups)
      missing_group2   <- setdiff(group2_ids, available_groups)
      
      if (length(missing_group1) > 0) {
        showNotification(paste("Group 1 combinations not found:", paste(missing_group1, collapse = ", ")),
                         type = "warning", duration = 10)
      }
      if (length(missing_group2) > 0) {
        showNotification(paste("Group 2 combinations not found:", paste(missing_group2, collapse = ", ")),
                         type = "warning", duration = 10)
      }
      
      group1_ids <- intersect(group1_ids, available_groups)
      group2_ids <- intersect(group2_ids, available_groups)
      
      if (length(group1_ids) == 0 || length(group2_ids) == 0) {
        showNotification("No valid group combinations found for comparison!", type = "error")
        removeModal()
        return()
      }
      
      message(paste("Comparing", length(group1_ids), "vs", length(group2_ids), "group-cluster combinations"))
      message(paste("Group 1 IDs:", paste(head(group1_ids, 3), collapse = ", "), "..."))
      message(paste("Group 2 IDs:", paste(head(group2_ids, 3), collapse = ", "), "..."))
      
      temp_res <- FindMarkers(
        object          = seurat_obj,
        ident.1         = group1_ids,
        ident.2         = group2_ids,
        group.by        = "temp_grouping",
        min.pct         = input$min_pct_compare_dataset_merge,
        logfc.threshold = input$logfc_threshold_datasets,
        return.thresh   = 0.05
      )
      
      seurat_obj$temp_grouping <- NULL
      
      if (nrow(temp_res) == 0) {
        showNotification("No differentially expressed genes found with these parameters!", type = "error")
        removeModal()
        return()
      }
      
      # Store raw numeric results for Venn diagrams BEFORE any formatting (rule #7)
      group1_text  <- paste(input$dataset1_compare, collapse = "_")
      group2_text  <- paste(input$dataset2_compare, collapse = "_")
      cluster_text <- if (input$all_clusters) "AllClusters" else paste(clusters_to_analyze, collapse = "_")
      table_name   <- paste0("Metadata_", group1_text, "_vs_", group2_text, "_", cluster_text, "_", format(Sys.time(), "%H%M%S"))
      description  <- paste0(metadata_col, ": [", group1_text, "] vs [", group2_text, "] | ",
                             cluster_col, ": ", cluster_text)
      parameters   <- list(
        metadata_col    = metadata_col,
        cluster_col     = cluster_col,
        group1          = input$dataset1_compare,
        group2          = input$dataset2_compare,
        clusters        = clusters_to_analyze,
        min_pct         = input$min_pct_compare_dataset_merge,
        logfc_threshold = input$logfc_threshold_datasets
      )
      
      store_de_table(
        table_data  = temp_res,
        table_name  = table_name,
        description = description,
        type        = "metadata_group",
        parameters  = parameters
      )
      
      # Add comparison label then store raw — formatting and filtering happen at render time
      temp_res$comparison <- sprintf(
        "Group1(%s: %s | %s: %s) vs Group2(%s: %s | %s: %s)",
        metadata_col, paste(input$dataset1_compare, collapse = ", "),
        cluster_col,  paste(if (input$all_clusters) "all" else clusters_to_analyze, collapse = ", "),
        metadata_col, paste(input$dataset2_compare, collapse = ", "),
        cluster_col,  paste(if (input$all_clusters) "all" else clusters_to_analyze, collapse = ", ")
      )
      
      diff_genes_compare_datasets(temp_res)
      shinyjs::enable("download_diff_dataset_cluster")
      
      showNotification(
        paste0("Found ", nrow(temp_res), " differentially expressed genes"),
        type = "message",
        duration = 5
      )
      
      removeModal()
      
    }, error = function(e) {
      showNotification(paste("Error during comparison:", e$message), type = "error", duration = 10)
      message(paste("Full error:", e$message))
      removeModal()
    })
  })
  
  output$diff_dataset_cluster <- renderDataTable({
    req(diff_genes_compare_datasets())
    df <- diff_genes_compare_datasets()
    threshold <- input$pval_adj_filter_datasets_merge
    if (!is.null(threshold) && !is.na(threshold)) {
      df <- df[df$p_val_adj <= threshold, , drop = FALSE]
    }
    cleaned <- clean_gene_names_for_html(rownames(df))
    df$gene_link <- paste0('<a href="#" class="gene-name" data-gene="', cleaned, '">', cleaned, '</a>')
    df$p_val     <- format_pvalue_robust(df$p_val)
    df$p_val_adj <- format_pvalue_robust(df$p_val_adj)
    display_cols <- c("gene_link", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj", "comparison")
    df <- df[, display_cols]
    colnames(df) <- c("Gene", "Log2FC", "Pct.1", "Pct.2", "P-value", "Adj. P-value", "Comparison")
    datatable(df, escape = FALSE,
              options = list(pageLength = 10, lengthMenu = c(10, 25, 50, 100, 200, 500),
                             order = list(list(1, "desc")), dom = "Blfrtip", scrollX = TRUE,
                             columnDefs = list(list(className = "dt-center", targets = 1:6),
                                               list(width = "150px", targets = 0))),
              class = "cell-border stripe", rownames = FALSE)
  })
  
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
  
  # Reactive value to store cluster composition data
  cluster_composition_multiple <- reactiveVal(NULL)
  
  
  ############################# Metadata Column Selectors for Composition Analysis #############################
  
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
  
  
  # Render cluster composition table for multiple datasets using existing function
  output$cluster_table_multiple <- renderDT({
    req(cluster_composition_multiple())
    render_cluster_composition_table(cluster_composition_multiple(), is_integrated = TRUE)
  })
  
  # Download handler for multiple datasets cluster composition using modular functions
  output$download_cluster_composition_multiple <- createDownloadHandler(
    reactive_data = reactive({ data <- cluster_composition_multiple(); data$Size_Bar <- NULL; return(data) }),
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets") }),
    data_name = "cluster_composition",
    download_type = "csv"
  )
  
  
  
  
  
  
  
  
  
  ###############Co-expression calculation#############
  
  
  
  # Reactive values for storing results
  gene_coexpression_data_multiple <- reactiveVal(NULL)
  coexpression_plot_multiple <- reactiveVal(NULL)
  
  # Update assay choices when object is loaded
  observe({
    req(multiple_datasets_object())
    updateSelectInput(session, "coexpr_assay_multiple", choices = names(multiple_datasets_object()@assays), selected = DefaultAssay(multiple_datasets_object()))
  })
  
  # Observer for co-expression analysis - FIXED
  observeEvent(input$analyze_coexpression_multiple, {
    req(multiple_datasets_object(), input$gene_text_coexpression_multiple)
    
    tryCatch({
      genes_input <- trimws(strsplit(input$gene_text_coexpression_multiple, ",")[[1]])
      genes_input <- genes_input[genes_input != ""]
      
      if (length(genes_input) < 2) {
        showNotification("Please enter at least 2 gene names", type = "error")
        return()
      }
      
      seurat_obj <- multiple_datasets_object()
      
      coexpr_results <- analyze_gene_coexpression(
        seurat_obj            = seurat_obj,
        genes                 = genes_input,
        assay_name            = input$coexpr_assay_multiple %||% "RNA",
        expression_thresholds = input$coexpr_threshold_multiple %||% 0,
        is_integrated         = TRUE
      )
      
      gene_coexpression_data_multiple(coexpr_results)
      
      if (!is.null(coexpr_results$genes_analyzed) && length(coexpr_results$genes_analyzed) > 0) {
        coexpression_plot_multiple(
          create_coexpression_plot(
            coexpr_results$data,
            coexpr_results$genes_analyzed,
            group_by        = input$coexpr_group_by_multiple %||% "cluster",
            secondary_split = if (isTRUE(input$coexpr_group_by_multiple == "dataset")) {
              input$coexpr_secondary_split_multiple_alt %||% "none"
            } else {
              input$coexpr_secondary_split_multiple %||% "none"
            }
          )
        )
      }
      
      showNotification("Co-expression analysis completed", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Comparison error:", conditionMessage(e)[1]), type = "error")
      message("Co-expression error details: ", conditionMessage(e)[1])
    })
  })
  # Render results using modular functions
  output$gene_coexpression_table_multiple <- renderDT({
    req(gene_coexpression_data_multiple())
    render_coexpression_table(gene_coexpression_data_multiple(), "gene_coexpression_table_multiple")
  })
  
  output$gene_coexpression_plot_multiple <- renderPlot({
    req(coexpression_plot_multiple())
    coexpression_plot_multiple()
  }, height = 700)
  
  # Render summary statistics using modular approach
  output$coexpression_summary_stats_multiple <- renderTable({
    req(gene_coexpression_data_multiple())
    create_coexpression_summary_stats(gene_coexpression_data_multiple()$data, gene_coexpression_data_multiple()$genes_analyzed)
  }, striped = TRUE, hover = TRUE, spacing = 'l')
  
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
  
  coexpr_de_results_multiple <- reactiveVal(NULL)
  
  # Populate cluster selector when object is available
  observe({
    req(multiple_datasets_object())
    tryCatch({
      clusters <- sort(as.character(unique(Idents(multiple_datasets_object()))))
      updateSelectInput(session, "coexpr_de_cluster_multiple", choices = clusters)
    }, error = function(e) NULL)
  })
  
  observeEvent(input$run_coexpr_de_multiple, {
    req(multiple_datasets_object(), input$gene_text_coexpression_multiple)
    tryCatch({
      genes_input <- trimws(strsplit(input$gene_text_coexpression_multiple, ",")[[1]])
      genes_input <- genes_input[nchar(genes_input) > 0]
      if (length(genes_input) < 2) {
        showNotification("Enter at least 2 genes (the LR pair)", type = "warning")
        return()
      }
      seurat_obj <- multiple_datasets_object()
      if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data))
        stop("seurat_clusters not found — run FindClusters first")

      showNotification("Running within-cluster DE...", type = "message", duration = 3)
      result <- run_coexpression_de(
        seurat_obj           = seurat_obj,
        genes                = genes_input,
        cluster              = input$coexpr_de_cluster_multiple,
        assay_name           = input$coexpr_assay_multiple %||% "RNA",
        expression_threshold = input$coexpr_threshold_multiple %||% 0
      )
      coexpr_de_results_multiple(result)
      showNotification(
        sprintf("DE complete — %d co-expressors vs %d others, %d genes",
                result$n_coexpressors, result$n_non_coexpressors, nrow(result$markers)),
        type = "message", duration = 5
      )
    }, error = function(e) {
      showNotification(paste("DE error:", conditionMessage(e)[1]), type = "error", duration = 8)
      message("coexpr_de error: ", conditionMessage(e)[1])
    })
  })
  
  output$coexpr_de_info_multiple <- renderUI({
    r <- coexpr_de_results_multiple()
    if (is.null(r)) return(NULL)
    
    bd <- r$breakdown
    breakdown_txt <- if (!is.null(bd$only_gene1)) {
      sprintf("(only %s: %d | only %s: %d | neither: %d)",
              r$genes_used[1], bd$only_gene1,
              r$genes_used[2], bd$only_gene2,
              bd$neither)
    } else {
      sprintf("(partial: %d | none: %d)", bd$partial, bd$none)
    }
    
    tags$div(
      tags$p(
        sprintf("Cluster: %s | threshold > %.2f | %d DE genes",
                r$cluster, r$threshold_used, nrow(r$markers)),
        style = "color: #aaa; font-size: 12px; margin-bottom: 2px;"
      ),
      tags$p(
        sprintf("✓ Co-expressors (both): %d cells  |  ✗ Non-co-expressors: %d cells %s",
                r$n_coexpressors, r$n_non_coexpressors, breakdown_txt),
        style = "color: #aaa; font-size: 12px; margin: 0;"
      )
    )
  })
  
  output$coexpr_de_table_multiple <- renderDT({
    req(coexpr_de_results_multiple())
    m <- coexpr_de_results_multiple()$markers
    m$avg_log2FC <- round(m$avg_log2FC, 3)
    m$pct.1      <- round(m$pct.1, 3)
    m$pct.2      <- round(m$pct.2, 3)
    m$p_val_adj  <- formatC(m$p_val_adj, format = "e", digits = 2)
    m$p_val      <- formatC(m$p_val,     format = "e", digits = 2)
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
      DT::formatStyle("avg_log2FC",
                      backgroundColor = DT::styleInterval(0, c("#c0392b22", "#1a572622")))
  })
  
  output$download_coexpr_de_multiple <- downloadHandler(
    filename = function() paste0("CoexprDE_cluster", input$coexpr_de_cluster_multiple,
                                 "_", format(Sys.time(), "%Y%m%d"), ".csv"),
    content  = function(file) {
      req(coexpr_de_results_multiple())
      write.csv(coexpr_de_results_multiple()$markers, file, row.names = FALSE)
    }
  )
  
  
  #####################Exclusive Biomarkers Discovery - Multiple Datasets###############################
  
  # Initialize reactive variable for exclusive markers
  exclusive_markers_merge <- reactiveVal()
  
  # Update dataset and cluster choices for exclusive biomarkers tab
  observe({
    req(multiple_datasets_object())
    
    # Update dataset choices if metadata has dataset column
    if ("dataset" %in% colnames(multiple_datasets_object()@meta.data)) {
      datasets <- unique(multiple_datasets_object()@meta.data$dataset)
      updateSelectInput(session, "exclusive_target_dataset_merge", 
                        choices = datasets)
    }
    
    # Update cluster choices
    updateSelectInput(session, "exclusive_target_cluster_merge", 
                      choices = levels(multiple_datasets_object()))
  })
  
  # Find exclusive biomarkers when button is clicked
  observeEvent(input$find_exclusive_markers_merge, {
    req(multiple_datasets_object(), input$exclusive_target_cluster_merge)
    
    tryCatch({
      showModal(modalDialog(
        title = "Finding Exclusive Biomarkers",
        "Analyzing gene expression patterns across integrated datasets...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      seurat_obj <- multiple_datasets_object()
      
      # Filter by dataset if specified
      if (!is.null(input$exclusive_target_dataset_merge) && 
          length(input$exclusive_target_dataset_merge) > 0 &&
          "dataset" %in% colnames(seurat_obj@meta.data)) {
        
        cells_to_keep <- colnames(seurat_obj)[seurat_obj@meta.data$dataset %in% input$exclusive_target_dataset_merge]
        seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
        
        message(paste0("Filtered to dataset(s): ", paste(input$exclusive_target_dataset_merge, collapse = ", ")))
      }
      
      # Call the function from cells_genes_expressions_newarch.R
      markers <- find_exclusive_biomarkers(
        seurat_obj = seurat_obj,
        target_cluster = input$exclusive_target_cluster_merge,
        min_pct_target = input$exclusive_min_pct_target_merge,
        max_pct_other = input$exclusive_max_pct_other_merge,
        min_log2fc = input$exclusive_min_log2fc_merge,
        detection_threshold = input$exclusive_detection_threshold_merge,
        min_mean_expr_target = input$exclusive_min_mean_expr_merge,
        statistical_test = input$exclusive_statistical_test_merge,
        max_pvalue = input$exclusive_max_pvalue_merge,
        pvalue_adjustment = "BH",
        assay_name = DefaultAssay(seurat_obj),
        top_n = NULL,
        verbose = TRUE
      )
      
      removeModal()
      
      if (nrow(markers) == 0) {
        showNotification("No exclusive biomarkers found with current parameters. Try relaxing the filters.", 
                         type = "error", 
                         duration = 10)
        exclusive_markers_merge(NULL)
        return()
      }
      
      # Store results for Venn diagrams BEFORE formatting
      markers_copy_venn <- markers
      markers_copy_venn$p_val <- ifelse("pvalue" %in% colnames(markers_copy_venn), 
                                        markers_copy_venn$pvalue, 
                                        NA)
      markers_copy_venn$p_val_adj <- ifelse("pvalue_adjusted" %in% colnames(markers_copy_venn), 
                                            markers_copy_venn$pvalue_adjusted, 
                                            NA)
      markers_copy_venn$avg_log2FC <- markers_copy_venn$log2_fold_change
      
      # Create descriptive table name
      dataset_text <- if (!is.null(input$exclusive_target_dataset_merge) && length(input$exclusive_target_dataset_merge) > 0) {
        paste0("Dataset", paste(input$exclusive_target_dataset_merge, collapse = "_"), "_")
      } else {
        "AllDatasets_"
      }
      
      cluster_text <- paste(input$exclusive_target_cluster_merge, collapse = "_")
      
      table_name <- paste0("MultipleDatasets_ExclusiveMarkers_", 
                           dataset_text,
                           "Cluster", cluster_text, 
                           "_", format(Sys.time(), "%H%M%S"))
      
      description <- paste0("Exclusive markers for cluster(s): ", 
                            paste(input$exclusive_target_cluster_merge, collapse = ", "))
      if (!is.null(input$exclusive_target_dataset_merge) && length(input$exclusive_target_dataset_merge) > 0) {
        description <- paste0(description, " in dataset(s): ", 
                              paste(input$exclusive_target_dataset_merge, collapse = ", "))
      }
      
      parameters <- list(
        target_datasets = input$exclusive_target_dataset_merge,
        target_clusters = input$exclusive_target_cluster_merge,
        min_pct_target = input$exclusive_min_pct_target_merge,
        max_pct_other = input$exclusive_max_pct_other_merge,
        min_log2fc = input$exclusive_min_log2fc_merge,
        detection_threshold = input$exclusive_detection_threshold_merge,
        statistical_test = input$exclusive_statistical_test_merge
      )
      
      result <- storeDETable(gene_table_storage(), 
                             markers_copy_venn, 
                             table_name, 
                             description, 
                             "exclusive_markers", 
                             parameters)
      if (result$success) {
        gene_table_storage(result$storage)  
      }
      
      # Format p-values for display if they exist
      if ("pvalue" %in% colnames(markers)) {
        markers$pvalue <- format_pvalue_robust(markers$pvalue)
      }
      if ("pvalue_adjusted" %in% colnames(markers)) {
        markers$pvalue_adjusted <- format_pvalue_robust(markers$pvalue_adjusted)
      }
      
      # Add gene links for HTML display
      cleaned_gene_names <- clean_gene_names_for_html(markers$gene)
      markers$gene_link <- paste0('<a href="#" class="gene-name" data-gene="',
                                  cleaned_gene_names, '">',
                                  cleaned_gene_names, '</a>')
      
      # Select and reorder columns for display
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
      
      # Store for display
      exclusive_markers_merge(markers_display)
      
      showNotification(paste("Found", nrow(markers), "exclusive biomarkers!"), 
                       type = "message")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error finding exclusive biomarkers:", conditionMessage(e)[1]),
                       type = "error", 
                       duration = 10)
      message("=== EXCLUSIVE BIOMARKERS ERROR ===\n", conditionMessage(e)[1])
      })
  })
  
  # Display results table
  output$table_exclusive_markers_merge <- renderDataTable({
    req(exclusive_markers_merge())
    df <- exclusive_markers_merge()
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
  
  #####################Exclusive Biomarkers Visualization - Multiple Datasets###############################
  
  # Reactive value to store the plot
  exclusive_plot_merge <- reactiveVal()
  
  
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
  
  
  # Generate plot when button is clicked
  observeEvent(input$generate_exclusive_plot_merge, {
    req(multiple_datasets_object(), input$exclusive_genes_to_plot_merge)
    
    tryCatch({
      # Parse genes from text input
      genes_input <- trimws(strsplit(input$exclusive_genes_to_plot_merge, ",")[[1]])
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
      seurat_obj <- multiple_datasets_object()
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
        paste("Creating", input$exclusive_plot_type_merge, "for", length(genes_input), "gene(s)..."),
        easyClose = FALSE,
        footer = NULL
      ))
      
      # Generate plot based on selected type
      plot_obj <- NULL
      # Generate the appropriate plot type
      if (input$exclusive_plot_type_merge == "dotplot") {
        # Check if split column is selected and valid
        split_column <- if (!is.null(input$exclusive_split_by_column_merge) && 
                            input$exclusive_split_by_column_merge != "" &&
                            input$exclusive_split_by_column_merge %in% colnames(seurat_obj@meta.data)) {
          input$exclusive_split_by_column_merge
        } else {
          NULL
        }
        
        # Filter by split column values if specified
        plot_seurat <- seurat_obj
        if (!is.null(split_column) && 
            !is.null(input$exclusive_split_values_filter) && 
            length(input$exclusive_split_values_filter) > 0) {
          
          cells_to_keep <- which(seurat_obj@meta.data[[split_column]] %in% input$exclusive_split_values_filter)
          
          if (length(cells_to_keep) > 0) {
            plot_seurat <- subset(seurat_obj, cells = cells_to_keep)
            message(paste("Filtered to", length(cells_to_keep), "cells"))
          } else {
            showNotification("No cells match the selected filter values", type = "warning")
            return()
          }
        }
        
        # Create DotPlot
        if (!is.null(split_column)) {
          plot_obj <- DotPlot(
            plot_seurat,
            features = genes_input,
            cols = c("lightgrey", "blue"),
            dot.scale = 8,
            split.by = split_column
          ) +
            RotatedAxis() +
            labs(
              title = paste("Expression of Exclusive Biomarkers by", split_column),
              x = "Genes",
              y = "Cluster"
            ) +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              legend.position = "right"
            )
        } else {
          plot_obj <- DotPlot(
            plot_seurat,
            features = genes_input,
            cols = c("lightgrey", "blue"),
            dot.scale = 8
          ) +
            RotatedAxis() +
            labs(
              title = "Expression of Exclusive Biomarkers across Clusters",
              x = "Genes",
              y = "Cluster"
            ) +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              legend.position = "right"
            )
        }
        
      } else if (input$exclusive_plot_type_merge == "violin") {
        # Check if split column is selected and valid
        split_column <- if (!is.null(input$exclusive_split_by_column_merge) && 
                            input$exclusive_split_by_column_merge != "" &&
                            input$exclusive_split_by_column_merge %in% colnames(seurat_obj@meta.data)) {
          input$exclusive_split_by_column_merge
        } else {
          NULL
        }
        
        # Filter by split column values if specified
        plot_seurat <- seurat_obj
        if (!is.null(split_column) && 
            !is.null(input$exclusive_split_values_filter) && 
            length(input$exclusive_split_values_filter) > 0) {
          
          cells_to_keep <- which(seurat_obj@meta.data[[split_column]] %in% input$exclusive_split_values_filter)
          
          if (length(cells_to_keep) > 0) {
            plot_seurat <- subset(seurat_obj, cells = cells_to_keep)
            message(paste("Filtered to", length(cells_to_keep), "cells"))
          } else {
            showNotification("No cells match the selected filter values", type = "warning")
            return()
          }
        }
        
        # Create VlnPlot
        if (!is.null(split_column)) {
          plot_obj <- VlnPlot(
            plot_seurat,
            features = genes_input,
            split.by = split_column,
            pt.size = 0
          ) +
            labs(
              title = paste("Expression of Exclusive Biomarkers by", split_column),
              x = "Cluster",
              y = "Expression Level"
            ) +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              legend.position = "right"
            )
        } else {
          plot_obj <- VlnPlot(
            plot_seurat,
            features = genes_input,
            pt.size = 0
          ) +
            labs(
              title = "Expression of Exclusive Biomarkers across Clusters",
              x = "Cluster",
              y = "Expression Level"
            ) +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
              axis.text.y = element_text(size = 10),
              legend.position = "right"
            )
        }
        
      } else if (input$exclusive_plot_type_merge == "feature") {
        # Check if split column is selected and valid
        split_column <- if (!is.null(input$exclusive_split_by_column_merge) && 
                            input$exclusive_split_by_column_merge != "" &&
                            input$exclusive_split_by_column_merge %in% colnames(seurat_obj@meta.data)) {
          input$exclusive_split_by_column_merge
        } else {
          NULL
        }
        
        # Filter by split column values if specified
        plot_seurat <- seurat_obj
        if (!is.null(split_column) && 
            !is.null(input$exclusive_split_values_filter) && 
            length(input$exclusive_split_values_filter) > 0) {
          
          cells_to_keep <- which(seurat_obj@meta.data[[split_column]] %in% input$exclusive_split_values_filter)
          
          if (length(cells_to_keep) > 0) {
            plot_seurat <- subset(seurat_obj, cells = cells_to_keep)
            message(paste("Filtered to", length(cells_to_keep), "cells"))
          } else {
            showNotification("No cells match the selected filter values", type = "warning")
            return()
          }
        }
        
        # Create FeaturePlot
        if (!is.null(split_column)) {
          plot_obj <- FeaturePlot(
            plot_seurat,
            features = genes_input,
            split.by = split_column,
            pt.size = 0.5,
            ncol = 2
          ) +
            labs(
              title = paste("Expression of Exclusive Biomarkers by", split_column)
            ) +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              legend.position = "right"
            )
        } else {
          plot_obj <- FeaturePlot(
            plot_seurat,
            features = genes_input,
            pt.size = 0.5,
            ncol = 2
          ) +
            labs(
              title = "Expression of Exclusive Biomarkers on UMAP"
            ) +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              legend.position = "right"
            )
        }
      }
      removeModal()
      
      if (!is.null(plot_obj)) {
        exclusive_plot_merge(plot_obj)
        showNotification("Plot generated successfully!", type = "message")
      }
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error generating plot:", e$message), type = "error", duration = 10)
      print(paste("Detailed error:", e$message))
    })
  })
  
  # Render the plot
  output$plot_exclusive_markers_merge <- renderPlot({
    req(exclusive_plot_merge())
    print(exclusive_plot_merge())
  }, height = 600)
  
  # Download handler for the plot
  output$download_exclusive_plot_merge <- downloadHandler(
    filename = function() {
      object_name <- getObjectNameForDownload(multiple_datasets_object(), default_name = "ExclusiveMarkers")
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      genes_text <- gsub(",", "_", input$exclusive_genes_to_plot_merge)
      genes_text <- gsub(" ", "", genes_text)
      genes_text <- substr(genes_text, 1, 50)  # Limit filename length
      
      paste0(object_name, "_ExclusiveMarkers_", genes_text, "_", 
             input$exclusive_plot_type_merge, "_", timestamp, ".", 
             input$exclusive_plot_format_merge)
    },
    content = function(file) {
      req(exclusive_plot_merge())
      
      tryCatch({
        plot_obj <- exclusive_plot_merge()
        
        # Determine plot dimensions based on number of genes
        genes_input <- trimws(strsplit(input$exclusive_genes_to_plot_merge, ",")[[1]])
        genes_input <- genes_input[genes_input != ""]
        n_genes <- length(genes_input)
        
        # Adjust dimensions based on plot type and split
        multiplier <- if (input$exclusive_split_by_dataset_merge) 1.5 else 1
        
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
  
  # Initialize reactive storage for Venn diagrams
  venn_plot_rendered <- reactiveVal(NULL)
  
  # Update select inputs with available gene tables
  observe({
    updateVennSelectInputs(session, gene_table_storage(), 
                           c("venn_table_1", "venn_table_2", "venn_table_3"))
  })
  
  # Enable/disable generate button based on available tables
  observe({
    tables <- gene_table_storage()
    if (length(tables) < 2) {
      shinyjs::disable("generate_venn_btn")
    } else {
      shinyjs::enable("generate_venn_btn")
    }
  })
  
  # Generate Venn diagram using modular functions
  observeEvent(input$generate_venn_btn, {
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
  })
  
  # Render Venn diagram
  output$venn_plot <- renderPlot({
    req(venn_plot_rendered())
    grid.draw(venn_plot_rendered())
  })
  
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
  
  # Download handlers using modular functions
  output$download_venn_diagram <- createDownloadHandler(
    reactive_data = venn_plot_rendered,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "VennDiagram") 
    }),
    data_name = "venn_comparison",
    download_type = "plot",
    plot_params = list(
      file_type = input$venn_diagram_format,
      width = 8,
      height = 6,
      dpi = input$venn_diagram_dpi
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
  
  # Function to store DE tables using modular approach
  store_de_table <- function(table_data, table_name, description, type, parameters) {
    # Ensure table is non-NULL and has rows
    if (is.null(table_data) || nrow(table_data) == 0) {
      showNotification("Cannot store an empty table.", type = "warning")
      return()
    }
    
    # Use modular storage function
    storeDETable(gene_table_storage(), table_data, table_name, description, type, parameters)
    
    # Update reactive storage
    current_storage <- gene_table_storage()
    current_storage[[table_name]] <- list(
      data = table_data,
      description = description,
      type = type,
      parameters = parameters,
      timestamp = Sys.time()
    )
    gene_table_storage(current_storage)
    
    showNotification(paste("Table", table_name, "stored successfully for Venn analysis"), 
                     type = "message")
  }
  ############################# Pie Chart for Cluster Composition #############################
  
  cluster_composition_data_merge <- reactiveVal(NULL)
  cluster_pie_plot_merge <- reactiveVal(NULL)
  
  # Update cluster choices
  observe({
    req(multiple_datasets_object())
    clusters <- levels(Idents(multiple_datasets_object()))
    updateSelectInput(session, "cluster_pie_select_merge", choices = clusters)
  })
  
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
  
  # Reactive value for storing volcano plot
  volcano_plot_merge <- reactiveVal(NULL)
  
  # Update volcano plot source selector when DEG results are available
  observe({
    tables <- gene_table_storage()
    if (length(tables) > 0) {
      table_choices <- c("None" = "", names(tables))
      updateSelectInput(session, "volcano_deg_source_merge", choices = table_choices)
    } else {
      updateSelectInput(session, "volcano_deg_source_merge", choices = c("None" = ""))
    }
  })
  
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
  
  # Render volcano plot
  output$volcano_plot_merge <- renderPlot({
    req(volcano_plot_merge())
    volcano_plot_merge()
  })
  
  # Download volcano plot
  output$download_volcano_plot_merge <- createDownloadHandler(
    reactive_data = volcano_plot_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(integrated_seurat_object(), default_name = "VolcanoPlot_Merge") 
    }),
    data_name = "volcano_plot",
    download_type = "plot",
    plot_params = list(
      file_type = input$volcano_plot_format_merge,
      width = 10,
      height = 8,
      dpi = input$volcano_plot_dpi_merge
    )
  )
  
  
  
  
  
  
  ############################## Subseting seurat object ##############################
  
  
  # Variable réactive pour cet onglet
  subset_seurat_merge <- reactiveVal(NULL)
  shinyjs::disable("download_subset_merge")
  
  # Mise à jour des choix d'identité de cellules pour le sous-ensemble
  observe({
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
  })
  
  # Réinitialisation de l'objet subset_seurat_merge lors du changement de multiple_datasets_object
  observe({
    if (!is.null(multiple_datasets_object())) {
      subset_seurat_merge(multiple_datasets_object())
    }
  })
  
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
      
      gene_list <- trimws(unlist(strsplit(input$gene_list_merge, ",")))
      gene_list <- gene_list[gene_list != ""]
      
      if (length(gene_list) == 0) {
        showNotification("Please enter valid gene names.", type = "error")
        return()
      }
      
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
  observe({
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
  })
  
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
      params_df <- collect_analysis_params(
        input      = input,
        module     = "multiple",
        qc_stats   = cs,
        seurat_obj = multiple_datasets_object()  # NULL if not yet integrated
      )
      write.csv(params_df, file, row.names = FALSE)
    }
  )
  output$download_report_multi <- downloadHandler(
    filename = function() {
      paste0("multi_dataset_report_", format(Sys.Date(), "%Y%m%d"), ".pdf")
    },
    content = function(file) {
      # Resolve script path locally — no external variable dependency
      script_candidates <- c(file.path(getwd(), "generate_report.py"), "/app/generate_report.py")
      script_path <- script_candidates[file.exists(script_candidates)][1]
      if (is.na(script_path)) stop("generate_report.py not found in project root or /app/")
      
      csv_tmp <- tempfile(fileext = ".csv")
      on.exit(unlink(csv_tmp), add = TRUE)
      cs <- qc_stats_multi()
      if (!is.null(cs)) {
        cs$integration_method <- integration_method_used()
        cs$harmony_vars       <- harmony_vars_used()
      }
      params_df <- collect_analysis_params(
        input      = input,
        module     = "multiple",
        qc_stats   = cs,
        seurat_obj = multiple_datasets_object()
      )
      write.csv(params_df, csv_tmp, row.names = FALSE)
      result <- system2(
        "python3",
        args   = c(shQuote(script_path),
                   "--csv",    shQuote(csv_tmp),
                   "--output", shQuote(file),
                   "--module", "multiple"),
        stdout = TRUE,
        stderr = TRUE
      )
      if (!file.exists(file)) {
        stop("PDF generation failed:\n", paste(result, collapse = "\n"))
      }
    }
  )
  
}
