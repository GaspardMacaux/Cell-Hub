############################## Script Cell Chat Analysis ##############################
# cellchat_server.R
cellchat_server <- function(input, output, session) {
  ############################## Loading Data ##############################
  
  # Initialize reactive values
  seurat_object_cellchat <- reactiveVal(NULL)
  db_cellchat <- reactiveVal(NULL)
  all_logs <- reactiveVal(character())  # Pour accumuler les logs
  db_loaded <- reactiveVal(FALSE)       # Pour tracker si DB déjà chargée
  
  
  
  # ADD THESE TWO LINES:
  cellchat_objects <- reactiveVal(list()) 
  analysis_log <- reactiveVal("")  
  
  # Loading modal function
  loadingModal <- function(text) {
    showModal(modalDialog(
      title = "Processing...",
      p(text),
      footer = NULL,
      easyClose = FALSE
    ))
  }
  
  # Add logs function
  add_logs <- function(new_logs) {
    current_logs <- all_logs()
    all_logs(c(current_logs, new_logs))
  }
  
 
  ############################## Load Seurat Object ##############################
  observeEvent(input$seurat_file_cellchat, {
    logs <- character()
    
    tryCatch({
      loadingModal("Loading Seurat object...")
      logs <- c(logs, "Starting Seurat object loading process...")
      
      loaded_data <- readRDS(input$seurat_file_cellchat$datapath)
      
      # Detect combined Seurat + CellChat structure saved from Cell-Hub
      if (is.list(loaded_data) && !inherits(loaded_data, "Seurat") &&
          !is.null(loaded_data$seurat) && !is.null(loaded_data$cellchat)) {
        
        message("Detected combined Seurat + CellChat structure in Seurat upload slot")
        
        if (!inherits(loaded_data$seurat, "Seurat")) {
          stop("Combined file does not contain a valid Seurat object under $seurat")
        }
        
        seurat_data <- loaded_data$seurat
        seurat_object_cellchat(seurat_data)
        
        # Also restore the CellChat objects and update selectors
        cellchat_objects(loaded_data$cellchat)
        update_cellchat_selectors(session, names(loaded_data$cellchat))
        
        n_analyzed <- sum(sapply(loaded_data$cellchat, function(obj) {
          !is.null(obj@net) && length(obj@net) > 0
        }))
        
        showNotification(
          sprintf("Already-analyzed object detected: Seurat + %d CellChat object(s) loaded (%d analyzed). You can go directly to visualization.",
                  length(loaded_data$cellchat), n_analyzed),
          type = "message",
          duration = 10
        )
        
      } else if (inherits(loaded_data, "Seurat")) {
        seurat_data <- loaded_data
        seurat_object_cellchat(seurat_data)
        
      } else {
        stop("File does not contain a valid Seurat object")
      }
      
      logs <- c(logs, sprintf("Seurat object loaded with %d cells", ncol(seurat_data)))
      
      output$seurat_info_cellchat <- renderPrint({
        cat(print_seurat_info(seurat_data))
      })
      
      add_logs(logs)
      removeModal()
      
    }, error = function(e) {
      logs <- c(logs, sprintf("ERROR loading Seurat: %s", e$message))
      add_logs(logs)
      removeModal()
      showModal(modalDialog(
        title = "Error",
        tags$div(p(sprintf("Error loading Seurat: %s", e$message))),
        easyClose = TRUE
      ))
    })
  })
  
  
  
  # Load previously analyzed CellChat objects
  observeEvent(input$load_analyzed_cellchat, {
    tryCatch({
      showModal(modalDialog("Loading analyzed objects...", footer = NULL))
      
      loaded_data <- readRDS(input$load_analyzed_cellchat$datapath)
      
      # Check if it's a combined Seurat + CellChat structure
      if (is.list(loaded_data) && !is.null(loaded_data$seurat) && !is.null(loaded_data$cellchat)) {
        # New format: combined structure
        message("Detected combined Seurat + CellChat structure")
        
        # Load Seurat object
        if (inherits(loaded_data$seurat, "Seurat")) {
          seurat_object_cellchat(loaded_data$seurat)
          message("Seurat object loaded successfully")
        }
        
        # Load CellChat objects
        objects <- loaded_data$cellchat
        
      } else if (is.list(loaded_data)) {
        # Old format: just CellChat objects
        message("Detected legacy format (CellChat objects only)")
        objects <- loaded_data
        
      } else {
        stop("File does not contain valid CellChat objects")
      }
      
      # Validate CellChat objects
      if (!is.list(objects) || length(objects) == 0) {
        stop("No CellChat objects found in file")
      }
      
      # Check if objects have analysis
      n_analyzed <- sum(sapply(objects, function(obj) {
        !is.null(obj@net) && length(obj@net) > 0
      }))
      
      cellchat_objects(objects)
      update_cellchat_selectors(session, names(objects))
      
      removeModal()
      
      # Show appropriate notification
      if (!is.null(loaded_data$seurat)) {
        showNotification(
          sprintf("Loaded Seurat + %d CellChat object(s) (%d analyzed)", 
                  length(objects), n_analyzed),
          type = "message",
          duration = 5
        )
      } else {
        showNotification(
          sprintf("Loaded %d CellChat object(s) (%d analyzed)", 
                  length(objects), n_analyzed),
          type = "message",
          duration = 5
        )
      }
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error loading:", e$message), type = "error")
      message(paste("Loading error:", e$message))
    })
  })
  ############################## Load Database ##############################
  observeEvent(input$species_cellchat, {
    db_loaded(FALSE)
  })
  
  
  observeEvent(input$load_db_cellchat, {
    logs <- character()
    tryCatch({
      loadingModal("Loading GaspouDB...")
      logs <- c(logs, paste("Loading GaspouDB for", input$species_cellchat))
      db_path <- if(input$species_cellchat == "mouse") {
        "databases/GaspouDB_mouse.rds"
      } else {
        "databases/GaspouDB_human.rds"
      }
      GaspouDB <- readRDS(db_path)
      db_cellchat(GaspouDB)
      db_loaded(TRUE)
      logs <- c(logs, "GaspouDB loaded successfully")
      logs <- c(logs, paste("Database contains", nrow(GaspouDB$interaction), "interactions"))
      output$db_info_cellchat <- renderPrint({
        print_db_info(GaspouDB)
      })
      add_logs(logs)
      removeModal()
    }, error = function(e) {
      logs <- c(logs, paste("ERROR loading database:", e$message))
      add_logs(logs)
      removeModal()
      showModal(modalDialog(
        title = "Error",
        tags$div(
          p(paste("Error loading GaspouDB:", e$message))
        ),
        easyClose = TRUE
      ))
    })
  })
  # Update the accumulated logs display
  output$loading_logs <- renderPrint({
    cat(paste(all_logs(), collapse = "\n"))
  })
  
  # Helper functions
  print_seurat_info <- function(seurat_obj) {
    # Use proper methods for Seurat object
    tryCatch({
      info <- list(
        ncells = ncol(seurat_obj),
        nfeatures = nrow(seurat_obj),
        assays = names(seurat_obj@assays)
      )
      
      return(sprintf(
        "Dataset loaded:\nNumber of cells: %d\nNumber of features: %d\nAvailable assays: %s",
        info$ncells,
        info$nfeatures,
        paste(info$assays, collapse = ", ")
      ))
    }, error = function(e) {
      return("Error reading Seurat object properties")
    })
  }
  
  print_db_info <- function(db) {
    paste(
      "Database loaded:",
      "\nNumber of interactions:", nrow(db$interaction),
      "\nSignaling pathways:", length(unique(db$interaction$pathway_name))
    )
  }
  
  
  ###################### TAB 2 ###################
  
  # Update column choices when Seurat object is loaded
  observe({
    req(seurat_object_cellchat())
    
    tryCatch({
      seurat_obj <- seurat_object_cellchat()
      
      # Safety check
      if (ncol(seurat_obj) == 0) {
        showNotification("Loaded Seurat object is empty!", type = "error")
        return()
      }
      
      # CRITICAL: Ensure ClusterIdents column exists with proper names
      if (!"ClusterIdents" %in% colnames(seurat_obj@meta.data)) {
        message("ClusterIdents column not found - creating it from current Idents")
        seurat_obj$ClusterIdents <- Idents(seurat_obj)
        seurat_object_cellchat(seurat_obj)
      }
      
      # Define patterns to exclude from column choices
      excluded_patterns <- c(
        "percent.mt", "nCount_ATAC", "nFeature_ATAC", "nFeature_RNA", "nCount_RNA",
        "^RNA_snn_res", "^RNA_snn", "^pANN", "^DF", "^integrated", "^integrated_snn"
      )
      
      pattern <- paste(excluded_patterns, collapse = "|")
      meta_cols <- colnames(seurat_obj@meta.data)
      filtered_cols <- meta_cols[!grepl(pattern, meta_cols)]
      
      # SMART DEFAULT SELECTION for primary column
      preferred_columns <- c("ClusterIdents", "cluster_names", "seurat_clusters", "ident")
      default_column <- NULL
      
      for (col in preferred_columns) {
        if (col %in% filtered_cols) {
          default_column <- col
          message(paste("Selected default grouping column:", col))
          break
        }
      }
      
      # If none of preferred columns exist, find first categorical column
      if (is.null(default_column)) {
        for (col in filtered_cols) {
          col_data <- seurat_obj@meta.data[[col]]
          n_unique <- length(unique(col_data))
          if (is.factor(col_data) || is.character(col_data) || n_unique < 50) {
            default_column <- col
            message(paste("No preferred column found, using first categorical column:", col))
            break
          }
        }
      }
      
      # Final fallback
      if (is.null(default_column) && length(filtered_cols) > 0) {
        default_column <- filtered_cols[1]
        message(paste("Using first available column as fallback:", default_column))
      }
      
      # Update all column selectors
      updateSelectInput(session, "group_by_cellchat", choices = filtered_cols, selected = default_column)
      updateSelectInput(session, "combine_with_column", choices = c("None" = "", filtered_cols), selected = "")
      updateSelectInput(session, "combine_with_column_2", choices = c("None" = "", filtered_cols), selected = "")
      updateSelectInput(session, "condition_filter_column", choices = filtered_cols)
      updateSelectInput(session, "subset_column_cellchat", choices = filtered_cols)
      
    }, error = function(e) {
      message("Error updating column choices:", e$message)
    })
  })
  
  
  # Preview of groups that will be created based on selected columns
  output$cellchat_groups_preview <- renderUI({
    req(seurat_object_cellchat(), input$group_by_cellchat)
    
    tryCatch({
      seurat_obj <- seurat_object_cellchat()
      
      # Collect selected columns
      selected_cols <- c(input$group_by_cellchat)
      
      if (!is.null(input$combine_with_column) && input$combine_with_column != "") {
        selected_cols <- c(selected_cols, input$combine_with_column)
      }
      
      if (!is.null(input$combine_with_column_2) && input$combine_with_column_2 != "") {
        selected_cols <- c(selected_cols, input$combine_with_column_2)
      }
      
      # Check for duplicates
      if (length(selected_cols) != length(unique(selected_cols))) {
        return(tags$p("⚠️ Duplicate columns selected", 
                      style = "color: #ff6b6b; font-size: 12px; margin: 0;"))
      }
      
      # Create preview of groups
      if (length(selected_cols) == 1) {
        # Simple grouping
        groups <- unique(seurat_obj@meta.data[[selected_cols[1]]])
        groups <- sort(groups[!is.na(groups)])
        n_groups <- length(groups)
        
        if (n_groups == 0) {
          return(tags$p("No groups found", 
                        style = "color: rgba(255,255,255,0.6); font-size: 11px; margin: 0;"))
        }
        
        preview_html <- tags$div(
          tags$p(
            tags$strong(paste(n_groups, "groups")), 
            " will be created from column: ",
            tags$code(selected_cols[1], style = "background-color: rgba(255,255,255,0.2); padding: 2px 5px; border-radius: 3px;"),
            style = "margin: 0 0 8px 0; color: white; font-size: 12px;"
          ),
          tags$p(
            paste(head(groups, 8), collapse = ", "),
            if(n_groups > 8) tags$span(paste0(" ... (", n_groups - 8, " more)"), style = "color: rgba(255,255,255,0.6);") else NULL,
            style = "margin: 0; color: rgba(255,255,255,0.85); font-size: 11px; font-family: monospace;"
          )
        )
        
      } else {
        # Combined grouping
        combined_values <- seurat_obj@meta.data[[selected_cols[1]]]
        for (col in selected_cols[-1]) {
          combined_values <- paste(combined_values, seurat_obj@meta.data[[col]], sep = "_")
        }
        
        groups <- unique(combined_values)
        groups <- sort(groups[!is.na(groups) & groups != "NA" & !grepl("NA_|_NA", groups)])
        n_groups <- length(groups)
        
        if (n_groups == 0) {
          return(tags$p("No valid groups found", 
                        style = "color: rgba(255,255,255,0.6); font-size: 11px; margin: 0;"))
        }
        
        preview_html <- tags$div(
          tags$p(
            tags$strong(paste(n_groups, "combined groups")), 
            " from: ",
            tags$code(paste(selected_cols, collapse = " + "), 
                      style = "background-color: rgba(255,255,255,0.2); padding: 2px 5px; border-radius: 3px;"),
            style = "margin: 0 0 8px 0; color: white; font-size: 12px;"
          ),
          tags$p(
            paste(head(groups, 6), collapse = ", "),
            if(n_groups > 6) tags$span(paste0(" ... (", n_groups - 6, " more)"), style = "color: rgba(255,255,255,0.6);") else NULL,
            style = "margin: 0; color: rgba(255,255,255,0.85); font-size: 11px; font-family: monospace;"
          )
        )
      }
      
      return(preview_html)
      
    }, error = function(e) {
      return(tags$p("Unable to generate preview", 
                    style = "color: rgba(255,255,255,0.5); font-size: 11px; margin: 0;"))
    })
  })
  # UI to select specific groups to keep
  output$cellchat_groups_filter_ui <- renderUI({
    req(seurat_object_cellchat(), input$group_by_cellchat)
    
    tryCatch({
      seurat_obj <- seurat_object_cellchat()
      
      # Collect selected columns
      selected_cols <- c(input$group_by_cellchat)
      
      if (!is.null(input$combine_with_column) && input$combine_with_column != "") {
        selected_cols <- c(selected_cols, input$combine_with_column)
      }
      
      if (!is.null(input$combine_with_column_2) && input$combine_with_column_2 != "") {
        selected_cols <- c(selected_cols, input$combine_with_column_2)
      }
      
      # Create groups list
      if (length(selected_cols) == 1) {
        groups <- unique(seurat_obj@meta.data[[selected_cols[1]]])
      } else {
        combined_values <- seurat_obj@meta.data[[selected_cols[1]]]
        for (col in selected_cols[-1]) {
          combined_values <- paste(combined_values, seurat_obj@meta.data[[col]], sep = "_")
        }
        groups <- unique(combined_values)
      }
      
      groups <- sort(groups[!is.na(groups) & groups != "NA" & !grepl("NA_|_NA", groups)])
      
      if (length(groups) == 0) {
        return(tags$p("No groups available", style = "color: rgba(255,255,255,0.6); font-size: 11px;"))
      }
      
      selectizeInput("selected_groups_filter",
                     label = NULL,
                     choices = groups,
                     selected = NULL,
                     multiple = TRUE,
                     options = list(
                       plugins = list('remove_button'),
                       placeholder = 'Select groups (empty = keep all)',
                       maxItems = 20
                     ))
      
    }, error = function(e) {
      return(NULL)
    })
  })
  
  
  
  
  # Update condition filter column choices when Seurat object is loaded
  observe({
    req(seurat_object_cellchat())
    
    tryCatch({
      seurat_obj <- seurat_object_cellchat()
      
      # Define patterns to exclude from condition column choices
      excluded_patterns <- c(
        "percent.mt", "nCount_ATAC", "nFeature_ATAC", "nFeature_RNA", "nCount_RNA",
        "^RNA_snn_res", "^RNA_snn", "^pANN", "^DF", "^integrated", "^integrated_snn",
        "cluster_names", "samples"
      )
      
      pattern <- paste(excluded_patterns, collapse = "|")
      meta_cols <- colnames(seurat_obj@meta.data)
      filtered_cols <- meta_cols[!grepl(pattern, meta_cols)]
      
      # Update condition filter column choices
      updateSelectInput(session, "condition_filter_column", 
                        choices = c("None" = "", filtered_cols))
      
      # Update combine with column choices
      updateSelectInput(session, "combine_with_column", 
                        choices = c("None" = "", filtered_cols))
      
    }, error = function(e) {
      print(paste("Error updating condition filter choices:", e$message))
    })
  })
  
  # Dynamic UI: show available clusters with checkboxes
  output$cluster_values_ui <- renderUI({
    req(seurat_object_cellchat())
    req(input$group_by_cellchat)
    
    tryCatch({
      seurat_obj <- seurat_object_cellchat()
      
      # Get unique clusters
      available_clusters <- unique(seurat_obj@meta.data[[input$group_by_cellchat]])
      available_clusters <- available_clusters[!is.na(available_clusters)]
      available_clusters <- sort(as.character(available_clusters))
      
      if (length(available_clusters) == 0) {
        return(tags$p("No clusters found", style = "color: red; font-size: 12px;"))
      }
      
      # Count cells per cluster
      cluster_counts <- table(seurat_obj@meta.data[[input$group_by_cellchat]])
      
      # Create choices: names = display labels, values = actual cluster names
      cluster_labels <- available_clusters
      names(cluster_labels) <- sapply(available_clusters, function(c) {
        sprintf("%s (%d cells)", c, cluster_counts[as.character(c)])
      })
      
      checkboxGroupInput(
        inputId = "selected_clusters_cellchat",
        label = NULL,
        choices = cluster_labels,
        selected = available_clusters  # All selected by default (using actual names)
      )
      
    }, error = function(e) {
      return(tags$p(paste("Error:", e$message), style = "color: red; font-size: 12px;"))
    })
  })
  
  output$condition_values_ui <- renderUI({
    req(seurat_object_cellchat())
    
    # If no column selected, show placeholder
    if (is.null(input$condition_filter_column) || input$condition_filter_column == "") {
      return(tags$p("Select a column above to see available values", 
                    style = "color: rgba(255,255,255,0.6); font-size: 12px; font-style: italic;"))
    }
    
    tryCatch({
      seurat_obj <- seurat_object_cellchat()
      selected_column <- input$condition_filter_column
      
      # Check if column exists
      if (!selected_column %in% colnames(seurat_obj@meta.data)) {
        return(tags$p("Column not found in metadata", 
                      style = "color: red; font-size: 12px;"))
      }
      
      # Get unique values from selected column
      available_values <- unique(seurat_obj@meta.data[[selected_column]])
      available_values <- available_values[!is.na(available_values)]
      available_values <- sort(as.character(available_values))
      
      if (length(available_values) == 0) {
        return(tags$p("No values found in this column", 
                      style = "color: red; font-size: 12px;"))
      }
      
      # Count cells per value
      value_counts <- table(seurat_obj@meta.data[[selected_column]])
      
      # CORRECT WAY: values first, then names
      condition_labels <- available_values
      names(condition_labels) <- sapply(available_values, function(v) {
        sprintf("%s (%d cells)", v, value_counts[as.character(v)])
      })
      
      checkboxGroupInput(
        inputId = "selected_conditions_cellchat",
        label = NULL,
        choices = condition_labels,
        selected = available_values
      )
      
    }, error = function(e) {
      return(tags$p(paste("Error:", e$message), style = "color: red; font-size: 12px;"))
    })
  })
  
  
  # Output for created objects info
  output$cellchat_objects_info <- renderPrint({
    cat(paste(analysis_logs(), collapse = "\n"))
  })
  
  # Create and analyze CellChat objects
  cellchat_objects <- reactiveVal(list())
  analysis_logs <- reactiveVal(character())
  
  add_analysis_log <- function(new_logs) {
    current_logs <- analysis_logs()
    analysis_logs(c(current_logs, new_logs))
  }
  
  # Function to create cluster names column
  create_cluster_names_column <- function(seurat_obj) {
    tryCatch({
      # If cluster_names already exists with real names, keep them
      if("cluster_names" %in% colnames(seurat_obj@meta.data)) {
        existing_names <- unique(seurat_obj@meta.data$cluster_names)
        
        # Check if names are meaningful (not just generic)
        if(!all(grepl("^Cluster_[0-9]+$", existing_names))) {
          message("Preserving existing cluster names")
          Idents(seurat_obj) <- "cluster_names"
          return(seurat_obj)
        }
      }
      
      # Create cluster names based on available data
      if("cell_type" %in% colnames(seurat_obj@meta.data)) {
        # Use cell_type if available
        seurat_obj$cluster_names <- as.character(seurat_obj$cell_type)
      } else if("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
        # Use seurat clusters
        seurat_obj$cluster_names <- paste0("Cluster_", seurat_obj$seurat_clusters)
      } else {
        # Use current identities
        seurat_obj$cluster_names <- as.character(Idents(seurat_obj))
      }
      
      # Set identities to cluster_names
      Idents(seurat_obj) <- "cluster_names"
      return(seurat_obj)
      
    }, error = function(e) {
      warning(paste("Could not create cluster names:", e$message))
      return(seurat_obj)
    })
  }
  
  
  # Function to clean CellChat object factors
  clean_cellchat_factors <- function(cellchat_obj) {
    tryCatch({
      message("=== Cleaning CellChat factors ===")
      
      # Clean identities
      if(is.factor(cellchat_obj@idents)) {
        cellchat_obj@idents <- droplevels(cellchat_obj@idents)
      }
      
      # Clean meta$labels
      if(!is.null(cellchat_obj@meta$labels)) {
        if(is.factor(cellchat_obj@meta$labels)) {
          cellchat_obj@meta$labels <- droplevels(cellchat_obj@meta$labels)
        }
      } else {
        # Create labels from idents if missing
        cellchat_obj@meta$labels <- cellchat_obj@idents
      }
      
      # Ensure consistency between idents and labels
      if(!identical(levels(cellchat_obj@idents), levels(cellchat_obj@meta$labels))) {
        message("Synchronizing idents and labels...")
        cellchat_obj@meta$labels <- factor(as.character(cellchat_obj@meta$labels), 
                                           levels = levels(cellchat_obj@idents))
      }
      
      # Clean all factor columns in meta
      for(col in names(cellchat_obj@meta)) {
        if(is.factor(cellchat_obj@meta[[col]])) {
          cellchat_obj@meta[[col]] <- droplevels(cellchat_obj@meta[[col]])
        }
      }
      
      message("Unique identities: ", length(levels(cellchat_obj@idents)))
      
      return(cellchat_obj)
      
    }, error = function(e) {
      stop(paste("Error cleaning CellChat factors:", e$message))
    })
  }
  
  # Dynamic UI: show available values when condition column is selected
  output$condition_values_ui <- renderUI({
    req(seurat_object_cellchat())
    req(input$condition_column_cellchat)
    
    if (is.null(input$condition_column_cellchat) || input$condition_column_cellchat == "") {
      return(NULL)
    }
    
    tryCatch({
      seurat_obj <- seurat_object_cellchat()
      
      # Get unique values from selected column
      available_values <- unique(seurat_obj@meta.data[[input$condition_column_cellchat]])
      available_values <- available_values[!is.na(available_values)]
      available_values <- sort(as.character(available_values))
      
      # Count cells per value
      value_counts <- table(seurat_obj@meta.data[[input$condition_column_cellchat]])
      value_labels <- sapply(available_values, function(v) {
        sprintf("%s (%d cells)", v, value_counts[v])
      })
      names(value_labels) <- available_values
      
      checkboxGroupInput(
        inputId = "selected_conditions_cellchat",
        label = "Select conditions to include:",
        choices = value_labels,
        selected = available_values  # All selected by default
      )
      
    }, error = function(e) {
      return(p("Error loading condition values", style = "color: red;"))
    })
  })
  
  output$cluster_values_ui <- renderUI({
    req(seurat_object_cellchat())
    req(input$group_by_cellchat)
    
    tryCatch({
      seurat_obj <- seurat_object_cellchat()
      
      # Get unique clusters
      available_clusters <- unique(seurat_obj@meta.data[[input$group_by_cellchat]])
      available_clusters <- available_clusters[!is.na(available_clusters)]
      available_clusters <- sort(as.character(available_clusters))
      
      if (length(available_clusters) == 0) {
        return(tags$p("No clusters found", style = "color: red; font-size: 12px;"))
      }
      
      # Count cells per cluster
      cluster_counts <- table(seurat_obj@meta.data[[input$group_by_cellchat]])
      
      # CORRECT WAY: values first, then names
      cluster_labels <- available_clusters
      names(cluster_labels) <- sapply(available_clusters, function(c) {
        sprintf("%s (%d cells)", c, cluster_counts[as.character(c)])
      })
      
      checkboxGroupInput(
        inputId = "selected_clusters_cellchat",
        label = NULL,
        choices = cluster_labels,
        selected = available_clusters
      )
      
    }, error = function(e) {
      return(tags$p(paste("Error:", e$message), style = "color: red; font-size: 12px;"))
    })
  })
  
  # Create and analyze CellChat object
  observeEvent(input$create_and_analyze_cellchat, {
    tryCatch({
      # Collect selected columns for grouping
      grouping_columns <- c(input$group_by_cellchat)
      
      if (!is.null(input$combine_with_column) && input$combine_with_column != "") {
        grouping_columns <- c(grouping_columns, input$combine_with_column)
      }
      
      if (!is.null(input$combine_with_column_2) && input$combine_with_column_2 != "") {
        grouping_columns <- c(grouping_columns, input$combine_with_column_2)
      }
      
      if (length(grouping_columns) != length(unique(grouping_columns))) {
        showNotification("Cannot use the same column multiple times!", type = "error")
        return()
      }
      
      obj_name <- trimws(input$cellchat_object_name)
      if (nchar(obj_name) == 0) {
        showNotification("Please provide an object name", type = "error")
        return()
      }
      
      current_objects <- cellchat_objects()
      if (obj_name %in% names(current_objects)) {
        showModal(modalDialog(
          title = "Object Already Exists",
          paste("An object named", obj_name, "already exists. Do you want to overwrite it?"),
          footer = tagList(
            modalButton("Cancel"),
            actionButton("confirm_overwrite", "Overwrite", class = "btn-danger")
          )
        ))
        return()
      }
      
      showModal(modalDialog(
        title = "Creating CellChat Object",
        tags$div(
          tags$h5("Please wait while CellChat object is being created and analyzed..."),
          tags$div(id = "cellchat_log",
                   style = "max-height: 300px; overflow-y: auto; background-color: #f5f5f5; padding: 10px; border-radius: 4px; font-family: monospace; font-size: 12px;")
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      
      add_analysis_log <- function(msg) {
        formatted_msg <- paste0(format(Sys.time(), "%H:%M:%S"), " - ", msg)
        add_logs(formatted_msg)
        current_logs <- all_logs()
        display_log <- paste(tail(current_logs, 30), collapse = "<br>")
        shinyjs::html("cellchat_log", display_log)
        shinyjs::runjs("document.getElementById('cellchat_log').scrollTop = document.getElementById('cellchat_log').scrollHeight;")
      }
      
      add_analysis_log("Starting CellChat object creation...")
      add_analysis_log(paste("Grouping by:", paste(grouping_columns, collapse = " + ")))
      
      result <- create_cellchat_object_multicol(
        seurat_obj          = seurat_object_cellchat(),
        grouping_columns    = grouping_columns,
        selected_groups     = input$selected_groups_filter,
        condition_column    = if (!is.null(input$condition_filter_column) && input$condition_filter_column != "") input$condition_filter_column else NULL,
        selected_conditions = input$condition_values,
        log_function        = add_analysis_log
      )
      
      if (!result$success) {
        removeModal()
        showNotification(paste("Failed to create CellChat object:", result$message), type = "error", duration = 10)
        return()
      }
      
      cellchat_obj <- result$object
      add_analysis_log("CellChat object created successfully!")
      
      add_analysis_log("Assigning database to CellChat object...")
      
      if (is.null(db_cellchat()) || !db_loaded()) {
        removeModal()
        showNotification("Please load the database first using 'Load Database' button", type = "error", duration = 10)
        return()
      }
      
      cellchat_obj@DB <- db_cellchat()
      add_analysis_log("GaspouDB database assigned successfully")
      
      add_analysis_log("Starting analysis pipeline...")
      
      analysis_result <- run_cellchat_pipeline(
        cellchat_obj = cellchat_obj,
        object_name  = obj_name,
        log_function = add_analysis_log,
        min_cells    = 10
      )
      
      if (!analysis_result$success) {
        removeModal()
        showNotification(paste("Analysis failed:", analysis_result$message), type = "error", duration = 10)
        return()
      }
      
      current_objects <- cellchat_objects()
      current_objects[[obj_name]] <- analysis_result$object
      cellchat_objects(current_objects)
      
      update_cellchat_selectors(session, names(current_objects))
      
      add_analysis_log(paste("Analysis complete! Found", analysis_result$n_interactions, "significant interactions"))
      
      removeModal()
      showNotification(
        paste("CellChat object", obj_name, "created and analyzed successfully!"),
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", conditionMessage(e)[1]), type = "error", duration = 10)
      message("Full error: ", conditionMessage(e))
    })
  })
  # Display analysis logs
  output$analysis_logs_cellchat <- renderPrint({
    cat(paste(analysis_logs(), collapse = "\n"))
  })
  
  
  #######################Comparison function##################################
  # Update comparison selectors
  observe({
    req(cellchat_objects())
    
    objects <- cellchat_objects()
    object_names <- names(objects)
    
    updateSelectInput(session, "compare_obj1_cellchat", choices = object_names)
    updateSelectInput(session, "compare_obj2_cellchat", choices = object_names)
  })
  
  # Update comparison cell type selectors based on BOTH objects
  observe({
    req(input$compare_obj1_cellchat, input$compare_obj2_cellchat)
    
    tryCatch({
      objects <- cellchat_objects()
      obj1 <- objects[[input$compare_obj1_cellchat]]
      obj2 <- objects[[input$compare_obj2_cellchat]]
      
      # Get cell types from both objects
      types1 <- if (!is.null(obj1@idents)) levels(obj1@idents) else character(0)
      types2 <- if (!is.null(obj2@idents)) levels(obj2@idents) else character(0)
      
      # Find common cell types (assuming same naming without condition suffix)
      # This extracts base names like "SatelliteCell" from "SatelliteCell_WT"
      common_types <- unique(c(types1, types2))
      
      updateSelectizeInput(session, "compare_sources_cellchat", choices = common_types)
      updateSelectizeInput(session, "compare_targets_cellchat", choices = common_types)
      
    }, error = function(e) {
      print(paste("Error updating comparison selectors:", e$message))
    })
  })
  
  # Reactive value for comparison plot
  comparison_plot <- reactiveVal()
  
  # Generate comparison plot
  observeEvent(input$generate_comparison_plot, {
    req(input$compare_obj1_cellchat, input$compare_obj2_cellchat)
    req(input$compare_sources_cellchat, input$compare_targets_cellchat)
    
    tryCatch({
      showNotification("Generating comparison plot...", type = "message", duration = 2)
      
      objects <- cellchat_objects()
      obj1 <- objects[[input$compare_obj1_cellchat]]
      obj2 <- objects[[input$compare_obj2_cellchat]]
      
      result <- generate_comparative_bubble_plot(
        cellchat_obj1 = obj1,
        cellchat_obj2 = obj2,
        obj1_name = input$compare_obj1_cellchat,
        obj2_name = input$compare_obj2_cellchat,
        sources = input$compare_sources_cellchat,
        targets = input$compare_targets_cellchat,
        threshold = input$compare_threshold
      )
      
      comparison_plot(result$plot)
      
      showNotification(
        sprintf("Comparison generated: %s (%d interactions) vs %s (%d interactions)",
                input$compare_obj1_cellchat, result$n_interactions_obj1,
                input$compare_obj2_cellchat, result$n_interactions_obj2),
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
    })
  })
  
  # Dynamic UI for comparison plot
  output$comparison_plot_ui <- renderUI({
    if (!is.null(comparison_plot())) {
      plotOutput("comparison_plot_display", 
                 height = paste0(input$compare_plot_height, "px"),
                 width = paste0(input$compare_plot_width, "px"))
    } else {
      div(class = "text-center text-muted", 
          style = "padding: 50px;",
          h4("No comparison generated yet"),
          p("Select two objects, sources, targets and click 'Generate Comparison'"))
    }
  })
  
  # Render comparison plot
  output$comparison_plot_display <- renderPlot({
    req(comparison_plot())
    print(comparison_plot())
  })
  
  # Download comparison plot
  output$download_comparison_plot <- downloadHandler(
    filename = function() {
      paste0("Comparison_", input$compare_obj1_cellchat, "_vs_", 
             input$compare_obj2_cellchat, "_", 
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      png(file, width = input$compare_plot_width, height = input$compare_plot_height, res = 300)
      print(comparison_plot())
      dev.off()
    }
  )
  
  #######################
  observe({
    req(cellchat_objects())
    req(length(cellchat_objects()) > 0)
    print("Starting initialization...")
    
    object_choices <- names(cellchat_objects())
    
    # Mise à jour unique de tous les sélecteurs d'objets
    updateSelectInput(session, "subset_choose_cellchat", choices = object_choices)
    updateSelectInput(session, "cellchat_obj_chord", choices = object_choices)
    updateSelectInput(session, "cellchat_obj_global", choices = object_choices)
    updateSelectInput(session, "cellchat_obj_specific", choices = object_choices)
  })
  
  # Observer pour mettre à jour les listes déroulantes
  observe({
    req(input$subset_choose_cellchat)
    print("Updating bubble plot inputs...")
    
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$subset_choose_cellchat]]
      req(!is.null(cellchat_obj@idents))
      cell_types <- levels(cellchat_obj@idents)
      
      # Préserver les sélections existantes
      selectedSources <- isolate(input$sources_use_cellchat)
      selectedTargets <- isolate(input$targets_use_cellchat)
      
      # Mise à jour des cell types
      updateSelectizeInput(session, "sources_use_cellchat",
                           choices = cell_types,
                           selected = intersect(selectedSources, cell_types),
                           options = list(plugins = list('remove_button'))
      )
      updateSelectizeInput(session, "targets_use_cellchat",
                           choices = cell_types,
                           selected = intersect(selectedTargets, cell_types),
                           options = list(plugins = list('remove_button'))
      )
      
      # Mise à jour des LR et pathways
      if(!is.null(cellchat_obj@net) && !is.null(cellchat_obj@net$prob)) {
        # Get L-R interactions
        lr_interactions <- NULL
        tryCatch({
          communications <- subsetCommunication(cellchat_obj)
          lr_interactions <- unique(communications$interaction_name)
        }, error = function(e) {
          message("Error getting L-R interactions: ", e$message)
          # Fallback method
          if (!is.null(cellchat_obj@LR$LRsig)) {
            lr_interactions <- unique(cellchat_obj@LR$LRsig$interaction_name)
          }
        })
        
        # Get pathways
        pathways <- NULL
        tryCatch({
          pathway_data <- subsetCommunication(cellchat_obj, slot.name = "netP")
          pathways <- unique(pathway_data$pathway_name)
        }, error = function(e) {
          message("Error getting pathways: ", e$message)
          # Fallback method
          if (!is.null(cellchat_obj@netP$pathways)) {
            pathways <- unique(unlist(lapply(cellchat_obj@netP$pathways, function(x) names(x))))
          }
        })
        
        # Update UI dropdowns
        if (!is.null(lr_interactions)) {
          updateSelectizeInput(session, "selected_lr_cellchat",
                               choices = lr_interactions,
                               selected = NULL
          )
        }
        
        if (!is.null(pathways)) {
          updateSelectizeInput(session, "signaling_pathways_bubble",
                               choices = pathways,
                               selected = NULL
          )
        }
      }
    }, error = function(e) {
      print(paste("Error updating bubble plot inputs:", e$message))
      showNotification(paste("Error updating inputs:", e$message), type = "error")
    })
  })
  
  
  
  output$download_analyzed_objects <- downloadHandler(
    filename = function() {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      paste0("CellHub_Seurat_CellChat_", timestamp, ".rds")
    },
    content = function(file) {
      req(cellchat_objects())
      
      showModal(modalDialog(
        title = "Saving Objects",
        "Saving Seurat + CellChat objects...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      # Save both Seurat and CellChat in a list
      combined_data <- list(
        seurat = seurat_object_cellchat(),
        cellchat = cellchat_objects()
      )
      
      saveRDS(combined_data, file)
      removeModal()
      
      showNotification("Seurat + CellChat objects saved successfully!", type = "message")
    }
  )
  
  
  
  # Reactive value for plot
  bubble_plot <- reactiveVal()
  
  
  
  
  # Function to extract types only for displayed interactions
  get_interaction_type_info <- function(cellchat_obj, communications) {
    # Check if the database contains the interaction_type field
    if (!is.null(cellchat_obj@DB$interaction) &&
        "interaction_type" %in% colnames(cellchat_obj@DB$interaction) &&
        "interaction_name" %in% colnames(cellchat_obj@DB$interaction)) {
      
      # Get only the displayed interaction names from the filtered communications
      displayed_interactions <- unique(communications$interaction_name)
      
      # Get all interaction types from the database
      all_interaction_info <- data.frame(
        interaction_name = cellchat_obj@DB$interaction$interaction_name,
        interaction_type = cellchat_obj@DB$interaction$interaction_type,
        stringsAsFactors = FALSE
      )
      
      # Filter to only the displayed interactions
      interaction_info <- all_interaction_info[all_interaction_info$interaction_name %in% displayed_interactions, ]
      
      # Handle any NA values
      interaction_info$interaction_type[is.na(interaction_info$interaction_type)] <- "Unknown"
      
      # Sort by type and name
      interaction_info <- interaction_info[order(interaction_info$interaction_type,
                                                 interaction_info$interaction_name), ]
      
      return(interaction_info)
    }
    
    return(NULL)
  }
  
  
  # Generate bubble plot 
  observeEvent(input$generate_plot_cellchat, {
    req(input$subset_choose_cellchat)
    
    showNotification("Generating bubble plot...", type = "message", duration = 2)
    
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$subset_choose_cellchat]]
      
      # Generate plot using helper
      result <- generate_bubble_plot(
        cellchat_obj = cellchat_obj,
        sources = input$sources_use_cellchat,
        targets = input$targets_use_cellchat,
        threshold = input$threshold_cellchat,
        flip_axes = input$flip_axes,
        exclude_intra = input$exclude_intra,
        signaling = NULL  # character vector of pathway names, e.g. c("MHC-II", "CXCL")
      )
      # Get interaction type info
      interaction_info <- get_interaction_type_info(cellchat_obj, result$communications)
      
      # Store plot and info
      bubble_plot(list(
        plot = result$plot,
        interaction_info = interaction_info
      ))
      
      # Notification
      if (!is.null(interaction_info) && "interaction_type" %in% colnames(interaction_info)) {
        type_counts <- table(interaction_info$interaction_type)
        type_summary <- paste(names(type_counts), type_counts, sep = ": ", collapse = ", ")
        showNotification(
          paste("Generated", result$n_interactions, "L-R pairs. By type:", type_summary),
          type = "message",
          duration = 5
        )
      } else {
        showNotification(
          paste("Generated", result$n_interactions, "ligand-receptor pairs"),
          type = "message"
        )
      }
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  
  # Get interaction type information
  output$interaction_info_table <- renderTable({
    req(bubble_plot())
    
    info_df <- bubble_plot()$interaction_info
    
    if (is.null(info_df) || nrow(info_df) == 0) {
      return(data.frame(
        Message = "No interaction type information available."
      ))
    }
    
    # Vérifier que les deux colonnes existent bien
    if (!("interaction_name" %in% colnames(info_df)) ||
        !("interaction_type" %in% colnames(info_df))) {
      return(data.frame(
        Message = "Expected columns not found in the data."
      ))
    }
    
    # Format the table for display - IMPORTANT: inclure explicitement les deux colonnes
    display_table <- data.frame(
      Interaction = info_df$interaction_name,
      Type = info_df$interaction_type
    )
    
    return(display_table)
  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = 'xs')
  
  # Render the plot without infinite refresh
  output$bubble_plot_cellchat <- renderPlot({
    req(bubble_plot())
    
    # Get the plot
    plot_obj <- bubble_plot()$plot
    
    # Return the plot
    plot_obj
  }, width = function() {
    # Use a default value if nothing is specified
    if(is.null(input$plot_width) || is.na(as.numeric(input$plot_width))) {
      return(1000)  # Default value
    }
    return(as.numeric(input$plot_width))
  }, height = function() {
    # Use a default value if nothing is specified
    if(is.null(input$plot_height) || is.na(as.numeric(input$plot_height))) {
      return(800)  # Default value
    }
    return(as.numeric(input$plot_height))
  })
  
  # Render the interaction info table
  output$interaction_info_table <- renderTable({
    req(bubble_plot())
    
    info_df <- bubble_plot()$interaction_info
    
    if (is.null(info_df) || nrow(info_df) == 0) {
      return(data.frame(
        Message = "No interaction type information available."
      ))
    }
    
    # Format the table for display
    display_table <- info_df %>%
      dplyr::rename(
        `Interaction` = interaction_name,
        `Type` = interaction_type
      )
    
    return(display_table)
  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = 'xs')
  
  # Download bubble plot - UTILISE LE HELPER EXISTANT
  output$download_bubble_plot <- createDownloadHandler(
    reactive_data = reactive({ bubble_plot()$plot }),  # Juste le plot
    object_name_reactive = reactive({ input$subset_choose_cellchat }),
    data_name = reactive({ 
      paste0("BubblePlot_", 
             paste(head(input$sources_use_cellchat, 2), collapse = "_"),
             "_to_",
             paste(head(input$targets_use_cellchat, 2), collapse = "_"))
    }),
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$bubble_plot_format }),
      width = reactive({ input$plot_width / 72 }),
      height = reactive({ input$plot_height / 72 }),
      dpi = 300
    )
  )
  
  # Download interaction info table - UTILISE LE HELPER
  output$download_info_table <- createDownloadHandler(
    reactive_data = reactive({ 
      info <- bubble_plot()$interaction_info
      if (is.null(info) || nrow(info) == 0) {
        return(data.frame(Message = "No interaction type information available"))
      }
      return(info)
    }),
    object_name_reactive = reactive({ input$subset_choose_cellchat }),
    data_name = "InteractionTypes",
    download_type = "csv"
  )
  
  #################Tab 3 : Circle plot##############
  observe({
    req(input$cellchat_obj_global)
    print("Updating circle plot inputs...")
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$cellchat_obj_global]]
      req(!is.null(cellchat_obj@idents))
      cell_types <- levels(cellchat_obj@idents)
      updateSelectizeInput(session, "cell_types_to_show", choices = cell_types)
      output$cellchat_obj_title_global <- renderText({
        paste("Selected object:", input$cellchat_obj_global)
      })
    }, error = function(e) {
      print(paste("Error updating circle plot inputs:", e$message))
    })
  })
  
  # Créer des reactiveVal pour stocker les plots
  global_circle_plot <- reactiveVal()
  cell_type_plots <- reactiveVal()
  FONT_CONFIG <- list(
    base = list(
      display = 72
    )
  )
  
  # Fonction pour calculer la taille de police
  get_font_scaling <- function(dpi, for_export = FALSE) {
    if(for_export) {
      return(FONT_CONFIG$base$display / dpi)
    }
    return(1)
  }
  
  # Fonction pour créer le plot global
  create_global_plot <- function(cellchat_obj, plot_type, font_size = 0.5, margin_size = 0.2,
                                 title = NULL, for_export = FALSE, dpi = NULL) {
    groupSize <- as.numeric(table(cellchat_obj@idents))
    font_scale <- get_font_scaling(dpi, for_export)
    adjusted_font_size <- font_size * font_scale
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mar = c(2, 2, 4, 2))
    par(cex = adjusted_font_size, cex.main = adjusted_font_size * 1.2)
    netVisual_circle(
      if(plot_type == "count") cellchat_obj@net$count else cellchat_obj@net$weight,
      vertex.weight = groupSize,
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = ifelse(!is.null(title),
                          paste(title, ifelse(plot_type == "count", "- Count", "- Weight")),
                          ifelse(plot_type == "count",
                                 "Number of interactions",
                                 "Interaction weights/strength")),
      margin = margin_size
    )
  }
  
  create_cell_type_plots <- function(cellchat_obj, cell_types, font_size = 0.5, margin_size = 0.15, title = NULL) {
    groupSize <- as.numeric(table(cellchat_obj@idents))
    mat <- cellchat_obj@net$weight
    
    n_plots <- length(cell_types)
    if(n_plots == 1) {
      par(mfrow=c(1,1), mar=c(2, 2, 4, 2), cex = font_size * 0.8, cex.main = font_size)
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[cell_types, ] <- mat[cell_types, ]
      netVisual_circle(mat2,
                       vertex.weight = groupSize,
                       weight.scale = TRUE,
                       edge.weight.max = max(mat),
                       title.name = ifelse(!is.null(title), paste(title, "-", cell_types), cell_types),
                       margin = margin_size)
    } else {
      grid_dim <- ceiling(sqrt(n_plots))
      layout(matrix(1:(grid_dim*grid_dim), grid_dim, grid_dim))
      for(cell_type in cell_types) {
        par(mar=c(1, 1, 3, 1), oma=c(2, 2, 2, 2), cex = font_size * 0.6, cex.main = font_size * 0.8)
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[cell_type, ] <- mat[cell_type, ]
        netVisual_circle(mat2,
                         vertex.weight = groupSize,
                         weight.scale = TRUE,
                         edge.weight.max = max(mat),
                         title.name = ifelse(!is.null(title), paste(title, "-", cell_type), cell_type),
                         margin = margin_size)
      }
    }
  }
  # Global circle plot with improved layout
  observeEvent(input$generate_global_plot, {
    req(input$cellchat_obj_global, input$plot_type_global)
    
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$cellchat_obj_global]]
      
      # Store plot function
      global_plot_function <- function() {
        groupSize <- as.numeric(table(cellchat_obj@idents))
        
        netVisual_circle(
          if(input$plot_type_global == "count") cellchat_obj@net$count else cellchat_obj@net$weight,
          vertex.weight = groupSize,
          weight.scale = TRUE,
          vertex.size = input$global_vertex_size,
          edge.width.max = input$global_edge_width,
          vertex.label.cex = input$global_label_size,
          label.edge = FALSE,
          title.name = paste(input$cellchat_obj_global, "-",
                             ifelse(input$plot_type_global == "count", 
                                    "Interaction Count", 
                                    "Interaction Strength")),
          margin = input$global_margin
        )
      }
      
      global_circle_plot(global_plot_function)
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Dynamic UI for global plot
  output$global_plot_ui <- renderUI({
    if(!is.null(global_circle_plot())) {
      plotOutput("global_plot_display", 
                 height = paste0(input$global_plot_size, "px"),
                 width = paste0(input$global_plot_size, "px"))
    } else {
      div(class = "text-center text-muted", 
          style = "padding: 50px;",
          h4("No plot generated yet"))
    }
  })
  
  # Render global plot
  output$global_plot_display <- renderPlot({
    req(global_circle_plot())
    par(mar = c(0, 0, 2, 0))
    global_circle_plot()()
  })
  
  # Generate cell-type specific plots - VERSION SIMPLIFIÉE
  observeEvent(input$generate_specific_plots, {
    req(input$cellchat_obj_specific, input$cell_types_to_show)
    
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$cellchat_obj_specific]]
      
      # Generate plots using helper
      result <- generate_specific_circle_plots(
        cellchat_obj = cellchat_obj,
        cell_types = input$cell_types_to_show,
        signal_direction = input$signal_direction,
        label_size = input$specific_label_size,
        margin = input$specific_margin,
        n_cols = input$n_cols_specific
      )
      
      # Store plot function
      cell_type_plots(result$plot_function)
      
      showNotification(
        paste("Generated", result$n_plots, "cell-type specific plots"),
        type = "message"
      )
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  
  
  # Dynamic UI for specific plots
  output$specific_plots_ui <- renderUI({
    if(!is.null(cell_type_plots())) {
      n_plots <- length(input$cell_types_to_show)
      n_cols <- min(input$n_cols_specific, n_plots)
      n_rows <- ceiling(n_plots / n_cols)
      
      total_height <- input$specific_plot_height * n_rows
      
      plotOutput("specific_plots_display", 
                 height = paste0(total_height, "px"))
    } else {
      div(class = "text-center text-muted", 
          style = "padding: 50px;",
          h4("No plots generated yet"))
    }
  })
  
  # Render specific plots
  output$specific_plots_display <- renderPlot({
    req(cell_type_plots())
    cell_type_plots()()
  })
  
  # Download handlers
  output$download_global_plot <- downloadHandler(
    filename = function() {
      paste0("GlobalCirclePlot_", input$cellchat_obj_global, "_", 
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      png(file, width = 2400, height = 2400, res = 300)
      par(mar = c(0, 0, 2, 0))
      global_circle_plot()()
      dev.off()
    }
  )
  
  output$download_specific_plots <- downloadHandler(
    filename = function() {
      paste0("CellTypeSpecificPlots_", input$cellchat_obj_specific, "_", 
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      n_plots <- length(input$cell_types_to_show)
      n_cols <- min(input$n_cols_specific, n_plots)
      n_rows <- ceiling(n_plots / n_cols)
      
      png(file, 
          width = 800 * n_cols, 
          height = 800 * n_rows, 
          res = 300)
      cell_type_plots()()
      dev.off()
    }
  )
  
  # Reactive value pour stocker le plot
  chord_plot <- reactiveVal()
  
  # Observer pour mettre à jour les inputs du chord plot
  observe({
    req(input$cellchat_obj_chord)
    print("Updating chord plot inputs...")
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$cellchat_obj_chord]]
      req(!is.null(cellchat_obj@idents))
      cell_types <- levels(cellchat_obj@idents)
      
      # Update selectizeInput instead of selectInput
      updateSelectizeInput(session, "sender_groups", 
                           choices = cell_types,
                           selected = NULL)
      updateSelectizeInput(session, "receiver_groups", 
                           choices = cell_types,
                           selected = NULL)
      
      if(!is.null(cellchat_obj@net) && !is.null(cellchat_obj@net$prob)) {
        all_comm <- subsetCommunication(cellchat_obj)
        if(!is.null(all_comm) && nrow(all_comm) > 0) {
          updateSelectizeInput(session, "selected_lr_chord",
                               choices = unique(all_comm$interaction_name))
        }
      }
      
      output$cellchat_obj_title_chord <- renderText({
        paste("Selected object:", input$cellchat_obj_chord)
      })
    }, error = function(e) {
      print(paste("Error updating chord plot inputs:", e$message))
    })
  })
  # Observer pour mettre à jour cell_types_to_show dans l'onglet specific
  observe({
    req(input$cellchat_obj_specific)
    print("Updating specific plot inputs...")
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$cellchat_obj_specific]]
      req(!is.null(cellchat_obj@idents))
      cell_types <- levels(cellchat_obj@idents)
      
      updateSelectizeInput(session, "cell_types_to_show", 
                           choices = cell_types,
                           selected = NULL)
      
      output$cellchat_obj_title_specific <- renderText({
        paste("Selected object:", input$cellchat_obj_specific)
      })
    }, error = function(e) {
      print(paste("Error updating specific plot inputs:", e$message))
    })
  })
  
  # Dynamic color inputs for chord plot — initialized from stored colors if available
  output$color_inputs_chord <- renderUI({
    tryCatch({
      req(input$cellchat_obj_chord)
      req(input$use_custom_colors)
      req(input$sender_groups)
      req(input$receiver_groups)
      
      all_groups   <- unique(c(input$sender_groups, input$receiver_groups))
      req(length(all_groups) > 0)
      
      cellchat_obj  <- cellchat_objects()[[input$cellchat_obj_chord]]
      stored_colors <- find_stored_colors(cellchat_obj, cluster_names = all_groups)
      
      n <- length(all_groups)
      fallback_pal <- if (n <= 8) {
        RColorBrewer::brewer.pal(max(3, n), "Set1")[seq_len(n)]
      } else {
        scales::hue_pal()(n)
      }
      names(fallback_pal) <- all_groups
      
      n_from_obj  <- if (!is.null(stored_colors)) sum(all_groups %in% names(stored_colors)) else 0
      info_banner <- if (n_from_obj == n) {
        tags$p(icon("circle-check"),
               sprintf(" Colors loaded from object (%d/%d)", n_from_obj, n),
               style = "color: #2e7d32; font-size: 12px; margin-bottom: 6px;")
      } else if (n_from_obj > 0) {
        tags$p(icon("triangle-exclamation"),
               sprintf(" Partial match: %d/%d colors from object", n_from_obj, n),
               style = "color: #e65100; font-size: 12px; margin-bottom: 6px;")
      } else {
        tags$p(icon("circle-info"),
               " No stored colors found — using default palette.",
               style = "color: #888; font-size: 12px; margin-bottom: 6px;")
      }
      
      color_inputs <- lapply(seq_along(all_groups), function(i) {
        group_name <- all_groups[i]
        init_color <- if (!is.null(stored_colors) && group_name %in% names(stored_colors)) {
          stored_colors[[group_name]]
        } else {
          fallback_pal[[group_name]]
        }
        colourInput(
          inputId = paste0("color_group_", i),
          label   = paste("Color for", group_name),
          value   = init_color
        )
      })
      
      div(class = "color-inputs", info_banner, color_inputs)
      
    }, error = function(e) {
      print(paste("Error in color inputs generation:", e$message))
      return(NULL)
    })
  })
  observeEvent(input$generate_chord_plot, {
    req(input$cellchat_obj_chord)
    req(length(input$sender_groups) > 0)
    req(length(input$receiver_groups) > 0)
    
    showNotification("Generating chord plot...", type = "message", duration = 2)
    
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$cellchat_obj_chord]]
      
      custom_colors <- NULL
      if (isTRUE(input$use_custom_colors)) {
        all_groups <- unique(c(input$sender_groups, input$receiver_groups))
        n          <- length(all_groups)
        
        # Same fallback palette as renderUI to stay consistent
        fallback_pal <- if (n <= 8) {
          RColorBrewer::brewer.pal(max(3, n), "Set1")[seq_len(n)]
        } else {
          scales::hue_pal()(n)
        }
        
        custom_colors <- sapply(seq_along(all_groups), function(i) {
          color_input <- input[[paste0("color_group_", i)]]
          if (is.null(color_input)) fallback_pal[i] else color_input
        })
        names(custom_colors) <- all_groups
      }
      
      result <- generate_chord_plot(
        cellchat_obj = cellchat_obj,
        senders      = input$sender_groups,
        receivers    = input$receiver_groups,
        threshold    = input$prob_threshold_chord,
        custom_colors = custom_colors,
        label_size   = input$chord_label_size
      )
      
      chord_plot(result$plot_function)
      showNotification(
        paste("Chord plot generated with", result$n_interactions, "interactions"),
        type = "message"
      )
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Dynamic UI for chord plot with adjustable size
  output$chord_plot_ui <- renderUI({
    if(!is.null(chord_plot())) {
      plotOutput("chord_plot_display", 
                 height = paste0(input$chord_plot_dims, "px"),
                 width = paste0(input$chord_plot_dims, "px"))
    } else {
      div(class = "text-center text-muted", 
          style = "padding: 50px;",
          h4("No plot generated yet"),
          p("Select parameters and click 'Generate Plot'"))
    }
  })
  
  # Render the chord plot
  output$chord_plot_display <- renderPlot({
    req(chord_plot())
    par(mar = c(0, 0, 0, 0))
    chord_plot()()  # Execute the stored function
  })
  
  # Improved download handler with better spacing
  output$download_chord_plot <- downloadHandler(
    filename = function() {
      obj_name <- input$cellchat_obj_chord
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      paste0("ChordPlot_", obj_name, "_", timestamp, ".", input$chord_plot_format)
    },
    content = function(file) {
      tryCatch({
        width_px <- input$chord_export_width
        height_px <- width_px  # Square plot
        
        if (input$chord_plot_format == "pdf") {
          pdf(file, width = width_px/72, height = height_px/72)
        } else if (input$chord_plot_format == "svg") {
          svg(file, width = width_px/72, height = height_px/72)
        } else {
          # PNG, JPEG, TIFF
          switch(input$chord_plot_format,
                 "png" = png(file, width = width_px, height = height_px, res = 300),
                 "jpeg" = jpeg(file, width = width_px, height = height_px, res = 300, quality = 100),
                 "tiff" = tiff(file, width = width_px, height = height_px, res = 300))
        }
        
        # Use larger margins for export to prevent label cutoff
        par(mar = c(2, 2, 2, 2))
        chord_plot()()  # Execute the stored function
        
        dev.off()
        showNotification("Plot downloaded successfully!", type = "message")
        
      }, error = function(e) {
        showNotification(paste("Download error:", e$message), type = "error")
      })
    }
  )
  
  
  # ---- Split chord plots ----
  split_chord_plot   <- reactiveVal(NULL)
  split_chord_lr_data <- reactiveVal(NULL)  
  
  observe({
    req(input$cellchat_obj_split_chord)
    tryCatch({
      obj <- cellchat_objects()[[input$cellchat_obj_split_chord]]
      ct  <- levels(obj@idents)
      updateSelectizeInput(session, "split_chord_senders",   choices = ct)
      updateSelectizeInput(session, "split_chord_receivers", choices = ct)
    }, error = function(e) NULL)
  })
  
  # Dynamic group definition UI for pairs mode — N groups, cell type choices from selected object
  output$split_chord_pairs_ui <- renderUI({
    req(input$split_chord_n_groups)
    n_grp <- max(2L, as.integer(input$split_chord_n_groups))
    ct <- tryCatch({
      obj <- cellchat_objects()[[input$cellchat_obj_split_chord]]
      if (!is.null(obj)) levels(obj@idents) else character(0)
    }, error = function(e) character(0))
    grp_colors <- tryCatch(
      RColorBrewer::brewer.pal(max(3L, n_grp), "Set2")[seq_len(n_grp)],
      error = function(e) scales::hue_pal()(n_grp)
    )
    ui_parts <- lapply(seq_len(n_grp), function(i) {
      tags$div(
        style = "margin-bottom: 4px;",
        tags$p(paste("Group", i),
               style = sprintf("color: %s; font-weight: bold; margin-bottom: 4px;", grp_colors[i])),
        textInput(paste0("split_chord_pair", i, "_label"), "Label:", value = paste("Group", i)),
        selectizeInput(paste0("split_chord_pair", i, "_senders"), "Sources:", choices = ct,
                       multiple = TRUE,
                       options  = list(placeholder = "Select sources", plugins = list("remove_button"))),
        selectizeInput(paste0("split_chord_pair", i, "_receivers"), "Targets:", choices = ct,
                       multiple = TRUE,
                       options  = list(placeholder = "Select targets", plugins = list("remove_button"))),
        if (i < n_grp) tags$hr(style = "border-color: rgba(255,255,255,0.2); margin: 8px 0;") else NULL
      )
    })
    do.call(tagList, ui_parts)
  })
  
  # Dynamic color inputs for split chord plot — initialized from stored colors if available
  output$split_chord_color_inputs <- renderUI({
    tryCatch({
      req(input$split_chord_custom_colors)
      all_cts <- if (isTRUE(input$split_chord_mode == "pairs")) {
        n_grp <- max(2L, as.integer(input$split_chord_n_groups))
        unique(unlist(lapply(seq_len(n_grp), function(i) {
          c(input[[paste0("split_chord_pair", i, "_senders")]],
            input[[paste0("split_chord_pair", i, "_receivers")]])
        })))
      } else {
        unique(c(input$split_chord_senders, input$split_chord_receivers))
      }
      req(length(all_cts) > 0)
      cellchat_obj  <- cellchat_objects()[[input$cellchat_obj_split_chord]]
      stored_colors <- find_stored_colors(cellchat_obj, cluster_names = all_cts)
      n <- length(all_cts)
      fallback_pal <- if (n <= 8) {
        RColorBrewer::brewer.pal(max(3, n), "Set1")[seq_len(n)]
      } else {
        scales::hue_pal()(n)
      }
      names(fallback_pal) <- all_cts
      n_from_obj  <- if (!is.null(stored_colors)) sum(all_cts %in% names(stored_colors)) else 0
      info_banner <- if (n_from_obj == n) {
        tags$p(icon("circle-check"),
               sprintf(" Colors loaded from object (%d/%d)", n_from_obj, n),
               style = "color: #2e7d32; font-size: 12px; margin-bottom: 6px;")
      } else if (n_from_obj > 0) {
        tags$p(icon("triangle-exclamation"),
               sprintf(" Partial match: %d/%d colors from object", n_from_obj, n),
               style = "color: #e65100; font-size: 12px; margin-bottom: 6px;")
      } else {
        tags$p(icon("circle-info"),
               " No stored colors found — using default palette.",
               style = "color: #888; font-size: 12px; margin-bottom: 6px;")
      }
      color_inputs <- lapply(seq_along(all_cts), function(i) {
        ct         <- all_cts[i]
        init_color <- if (!is.null(stored_colors) && ct %in% names(stored_colors)) {
          stored_colors[[ct]]
        } else {
          fallback_pal[[ct]]
        }
        colourpicker::colourInput(
          inputId = paste0("split_chord_color_", i),
          label   = paste("Color for", ct),
          value   = init_color
        )
      })
      div(info_banner, color_inputs)
    }, error = function(e) {
      message("split_chord_color_inputs error: ", conditionMessage(e)[1])
      NULL
    })
  })
  
  
  # Load LR pairs and categorize: common vs unique per group
  observeEvent(input$load_split_chord_lr, {
    req(input$cellchat_obj_split_chord)
    tryCatch({
      obj    <- cellchat_objects()[[input$cellchat_obj_split_chord]]
      thresh <- input$split_chord_threshold
      mode   <- input$split_chord_mode
      
      if (mode == "pairs") {
        n_grp <- max(2L, as.integer(input$split_chord_n_groups))
        groups_def <- lapply(seq_len(n_grp), function(i) {
          lbl <- trimws(input[[paste0("split_chord_pair", i, "_label")]])
          if (nchar(lbl) == 0) lbl <- paste("Group", i)
          list(
            label   = lbl,
            sources = input[[paste0("split_chord_pair", i, "_senders")]],
            targets = input[[paste0("split_chord_pair", i, "_receivers")]]
          )
        })
        missing_grps <- which(sapply(groups_def, function(g) {
          length(g$sources) == 0 || length(g$targets) == 0
        }))
        if (length(missing_grps) > 0)
          stop(sprintf("Groups %s: please select at least one source and one target",
                       paste(missing_grps, collapse = ", ")))
      } else {
        req(input$split_chord_senders, input$split_chord_receivers)
        split_by     <- input$split_chord_by
        sources      <- input$split_chord_senders
        targets      <- input$split_chord_receivers
        split_groups <- if (split_by == "source") sources else targets
        groups_def <- lapply(split_groups, function(grp) {
          list(label   = grp,
               sources = if (split_by == "source") grp else sources,
               targets = if (split_by == "target") grp else targets)
        })
      }
      
      split_groups <- sapply(groups_def, `[[`, "label")
      message("=== load_split_chord_lr ===")
      message("--- mode: ", mode, " | groups: ", paste(split_groups, collapse = ", "))
      
      per_group_lr <- setNames(lapply(groups_def, function(gd) {
        comms <- tryCatch(
          subsetCommunication(obj,
                              sources.use = gd$sources,
                              targets.use = gd$targets,
                              thresh      = thresh),
          error = function(e) NULL
        )
        if (is.null(comms) || nrow(comms) == 0) {
          message(sprintf("  [%s] no interactions found", gd$label))
          return(data.frame(lr_key = character(), prob = numeric()))
        }
        comms$lr_key <- paste(comms$ligand, comms$receptor, sep = " -> ")
        agg <- aggregate(prob ~ lr_key, data = comms, FUN = max)
        agg <- agg[order(agg$prob, decreasing = TRUE), ]
        message(sprintf("  [%s] %d LR pairs found", gd$label, nrow(agg)))
        agg
      }), split_groups)
      
      all_keys      <- unique(unlist(lapply(per_group_lr, function(d) d$lr_key)))
      key_to_groups <- setNames(lapply(all_keys, function(k) {
        sort(names(which(sapply(per_group_lr, function(d) k %in% d$lr_key))))
      }), all_keys)
      
      signatures <- unique(lapply(key_to_groups, identity))
      
      combo_lr <- list()
      for (sig in signatures) {
        label <- if (length(sig) == length(split_groups)) {
          "Common to all"
        } else if (length(sig) == 1) {
          paste0("Unique to ", sig[1])
        } else {
          paste0("Shared: ", paste(sig, collapse = " + "))
        }
        keys_in_sig <- names(which(sapply(key_to_groups, function(g) identical(g, sig))))
        probs <- sapply(keys_in_sig, function(k) {
          mean(sapply(sig, function(grp) {
            d <- per_group_lr[[grp]]
            v <- d$prob[d$lr_key == k]
            if (length(v) == 0) 0 else v[1]
          }))
        })
        combo_lr[[label]] <- keys_in_sig[order(probs, decreasing = TRUE)]
        message(sprintf("  [%s] %d LR pairs", label, length(keys_in_sig)))
      }
      
      sort_order <- sapply(names(combo_lr), function(label) {
        if (label == "Common to all")       return(0)
        if (startsWith(label, "Shared: "))  return(1)
        return(2)
      })
      combo_lr <- combo_lr[order(sort_order)]
      
      split_chord_lr_data(list(
        combo_lr     = combo_lr,
        per_group    = per_group_lr,
        split_groups = split_groups,
        groups_def   = groups_def
      ))
      
      total_lr <- length(all_keys)
      n_common <- length(combo_lr[["Common to all"]])
      showNotification(
        sprintf("Found %d LR pairs total (%d common to all, %d combinations)",
                total_lr, n_common, length(combo_lr)),
        type = "message", duration = 4
      )
      
    }, error = function(e) {
      showNotification(paste("Error:", conditionMessage(e)[1]), type = "error", duration = 5)
      message("load_split_chord_lr error: ", conditionMessage(e)[1])
    })
  })
  
  # Build dynamic LR selector UI: one selectize for common + one per group
  output$split_chord_lr_selector <- renderUI({
    lr_data <- split_chord_lr_data()
    if (is.null(lr_data)) {
      return(tags$p("Click 'Load LR pairs' first",
                    style = "color: #ccc; font-style: italic;"))
    }
    combo_lr  <- lr_data$combo_lr
    per_group <- lr_data$per_group
    
    ui_parts <- lapply(seq_along(combo_lr), function(i) {
      label    <- names(combo_lr)[i]
      keys     <- combo_lr[[i]]
      input_id <- paste0("split_chord_lr_combo_", i)
      if (length(keys) == 0) return(NULL)
      
      probs <- sapply(keys, function(k) {
        vals <- sapply(per_group, function(d) {
          v <- d$prob[d$lr_key == k]
          if (length(v) == 0) NA_real_ else v[1]
        })
        round(mean(vals, na.rm = TRUE), 4)
      })
      choices <- setNames(keys, sprintf("%s  [%.4f]", keys, probs))
      
      lbl_icon <- if (label == "Common to all")     icon("link") else
        if (startsWith(label, "Shared: ")) icon("code-branch") else
          icon("circle-dot")
      
      # "Select all" available for Common-to-all and Unique-to-X (not for Shared: which are partial)
      is_selectable <- label == "Common to all" || startsWith(label, "Unique to ")
      select_all_btn <- if (is_selectable) {
        btn_id    <- paste0("split_chord_select_all_", i)
        btn_label <- if (label == "Common to all") {
          tags$span(icon("check-double"), " Select all common to all",
                    style = "font-size: 12px; color: #ccc;")
        } else {
          tags$span(icon("check-double"), " Select all unique to this group",
                    style = "font-size: 12px; color: #ccc;")
        }
        tags$div(style = "margin-bottom: 4px;",
                 checkboxInput(btn_id, label = btn_label, value = FALSE))
      } else NULL
      
      tags$div(
        style = "margin-bottom: 8px;",
        select_all_btn,
        selectizeInput(
          input_id,
          label    = tags$span(lbl_icon, paste0(" ", label, ":")),
          choices  = choices,
          selected = character(0),
          multiple = TRUE,
          options  = list(plugins     = list("remove_button"),
                          placeholder = "None selected",
                          maxOptions  = 500)
        )
      )
    })
    
    do.call(tagList, Filter(Negate(is.null), ui_parts))
  })
  
  # Auto-select all unique pairs when checkbox is ticked
  observe({
    lr_data <- split_chord_lr_data()
    req(lr_data)
    combo_lr <- lr_data$combo_lr
    
    lapply(seq_along(combo_lr), function(i) {
      label <- names(combo_lr)[i]
      if (!startsWith(label, "Unique to ") && label != "Common to all") return(NULL)
      btn_id   <- paste0("split_chord_select_all_", i)
      input_id <- paste0("split_chord_lr_combo_", i)
      keys     <- combo_lr[[i]]
      checked  <- isTRUE(input[[btn_id]])
      if (checked) {
        updateSelectizeInput(session, input_id, selected = keys)
      } else {
        updateSelectizeInput(session, input_id, selected = character(0))
      }
    })
  })
  
  # Generate split chord plots
  observeEvent(input$generate_split_chord, {
    mode <- input$split_chord_mode
    req(input$cellchat_obj_split_chord)
    if (mode == "pairs") {
      n_grp <- max(2L, as.integer(input$split_chord_n_groups))
      for (i in seq_len(n_grp)) {
        req(input[[paste0("split_chord_pair", i, "_senders")]],
            input[[paste0("split_chord_pair", i, "_receivers")]])
      }
    } else {
      req(input$split_chord_senders, input$split_chord_receivers)
    }
    
    tryCatch({
      showNotification("Generating split chord plots...", type = "message", duration = 2)
      obj     <- cellchat_objects()[[input$cellchat_obj_split_chord]]
      lr_data <- split_chord_lr_data()
      
      if (is.null(lr_data)) {
        showNotification("Please click 'Load LR pairs' first", type = "warning")
        return()
      }
      
      n_combos    <- length(lr_data$combo_lr)
      selected_lr <- unique(unlist(lapply(seq_len(n_combos), function(i) {
        input[[paste0("split_chord_lr_combo_", i)]]
      })))
      selected_lr <- selected_lr[!is.null(selected_lr) & nchar(selected_lr) > 0]
      
      if (length(selected_lr) == 0) {
        showNotification("Please select at least one LR pair", type = "warning")
        return()
      }
      
      message(sprintf("generate_split_chord: %d LR pairs selected", length(selected_lr)))
      
      all_cts <- if (mode == "pairs") {
        n_grp <- max(2L, as.integer(input$split_chord_n_groups))
        unique(unlist(lapply(seq_len(n_grp), function(i) {
          c(input[[paste0("split_chord_pair", i, "_senders")]],
            input[[paste0("split_chord_pair", i, "_receivers")]])
        })))
      } else {
        unique(c(input$split_chord_senders, input$split_chord_receivers))
      }
      
      custom_colors <- NULL
      if (isTRUE(input$split_chord_custom_colors)) {
        n <- length(all_cts)
        fallback_pal <- if (n <= 8) {
          RColorBrewer::brewer.pal(max(3, n), "Set1")[seq_len(n)]
        } else {
          scales::hue_pal()(n)
        }
        custom_colors <- sapply(seq_along(all_cts), function(i) {
          col <- input[[paste0("split_chord_color_", i)]]
          if (is.null(col)) fallback_pal[i] else col
        })
        names(custom_colors) <- all_cts
      }
      
      pairs_arg <- if (mode == "pairs") lr_data$groups_def else NULL
      
      result <- generate_split_chord_plots(
        cellchat_obj     = obj,
        sources          = if (mode == "split") input$split_chord_senders   else NULL,
        targets          = if (mode == "split") input$split_chord_receivers else NULL,
        split_by         = if (mode == "split") input$split_chord_by        else "source",
        pairs            = pairs_arg,
        selected_lr_keys = selected_lr,
        threshold        = input$split_chord_threshold,
        label_distance   = input$split_chord_label_distance,
        custom_colors    = custom_colors,
        show_title       = isTRUE(input$split_chord_show_title),
        gap_inner        = input$split_chord_gap_inner,
        gap_outer        = input$split_chord_gap_outer,
        plot_size        = as.integer(input$split_chord_dims)
      )
      
      split_chord_plot(result)
      
    }, error = function(e) {
      showNotification(paste("Error:", conditionMessage(e)[1]), type = "error", duration = 6)
      message("generate_split_chord error: ", conditionMessage(e)[1])
    })
  })
  
  
  split_chord_h <- reactiveVal(600)
  split_chord_w <- reactiveVal(600)
  
  observeEvent(split_chord_plot(), {
    req(split_chord_plot())
    split_chord_h(split_chord_plot()$png_size)
    split_chord_w(split_chord_plot()$png_size * split_chord_plot()$n_plots)
  })
  
  output$split_chord_plot_display <- renderImage({
    req(split_chord_plot())
    result <- split_chord_plot()
    list(
      src     = result$output_path,
      contentType = "image/png",
      width   = result$png_size * result$n_plots,
      height  = result$png_size,
      alt     = "Split chord plot"
    )
  }, deleteFile = FALSE)
  
  output$download_split_chord <- downloadHandler(
    filename = function() {
      generateFileName(input$cellchat_obj_split_chord, "SplitChord", input$split_chord_format)
    },
    content = function(file) {
      req(split_chord_plot())
      result    <- split_chord_plot()
      n_plots   <- result$n_plots
      width_px  <- as.integer(input$split_chord_export_width)
      height_px <- width_px
      
      if (input$split_chord_format == "svg") {
        svg(file, width = (width_px / 72) * n_plots, height = width_px / 72)
        result$export_function(n_plots)
        dev.off()
      } else if (input$split_chord_format == "pdf") {
        pdf(file, width = (width_px / 72) * n_plots, height = width_px / 72)
        result$export_function(n_plots)
        dev.off()
      } else {
        switch(input$split_chord_format,
               "png"  = png(file,  width = width_px * n_plots, height = height_px, res = 300),
               "tiff" = tiff(file, width = width_px * n_plots, height = height_px, res = 300),
               "jpeg" = jpeg(file, width = width_px * n_plots, height = height_px,
                             res = 300, quality = 100)
        )
        result$plot_function()
        dev.off()
      }
    }
  )
  # Summary table: all LR pairs found, prob per group, category
  output$split_chord_lr_table_ui <- renderUI({
    lr_data <- split_chord_lr_data()
    if (is.null(lr_data)) {
      return(tags$p("Click 'Load LR pairs' first",
                    style = "color: #ccc; font-style: italic;"))
    }
    
    split_groups <- lr_data$split_groups
    per_group    <- lr_data$per_group
    combo_lr     <- lr_data$combo_lr
    groups_def   <- lr_data$groups_def  # NULL in split mode
    
    # Build human-readable column names: "Label (Src → Tgt)"
    group_col_names <- sapply(seq_along(split_groups), function(i) {
      grp <- split_groups[i]
      gd  <- if (!is.null(groups_def)) groups_def[[i]] else NULL
      if (!is.null(gd)) {
        src <- paste(gd$sources, collapse = "+")
        tgt <- paste(gd$targets, collapse = "+")
        paste0(grp, " (", src, " \u2192 ", tgt, ")")
      } else {
        grp
      }
    })
    
    # Build one row per LR key — split ligand/receptor
    all_keys <- unique(unlist(combo_lr))
    rows <- lapply(all_keys, function(k) {
      parts    <- strsplit(k, " -> ")[[1]]
      ligand   <- if (length(parts) >= 1) parts[1] else k
      receptor <- if (length(parts) >= 2) parts[2] else ""
      probs <- sapply(split_groups, function(grp) {
        d <- per_group[[grp]]
        v <- d$prob[d$lr_key == k]
        if (length(v) == 0) NA_real_ else round(v[1], 5)
      })
      category <- names(which(sapply(combo_lr, function(keys) k %in% keys)))[1]
      as.data.frame(
        t(c("Ligand"   = ligand,
            "Receptor" = receptor,
            setNames(as.character(probs), group_col_names),
            "Category" = category)),
        stringsAsFactors = FALSE
      )
    })
    
    df <- do.call(rbind, rows)
    
    priority <- ifelse(df$Category == "Common to all", 0,
                       ifelse(startsWith(df$Category, "Shared: "), 1, 2))
    max_prob <- apply(df[, group_col_names, drop = FALSE], 1, function(r) {
      vals <- suppressWarnings(as.numeric(r))
      if (all(is.na(vals))) -Inf else max(vals, na.rm = TRUE)
    })
    df <- df[order(priority, -max_prob), ]
    
    all_labels <- unique(df$Category)
    cat_colors <- setNames(scales::hue_pal()(length(all_labels)), all_labels)
    if ("Common to all" %in% all_labels) cat_colors["Common to all"] <- "#2e7d32"
    unique_labels <- all_labels[startsWith(all_labels, "Unique to ")]
    if (length(unique_labels) > 0)
      cat_colors[unique_labels] <- scales::brewer_pal(palette = "Oranges")(
        max(3, length(unique_labels)))[seq_along(unique_labels) + 1]
    
    prob_col_indices <- seq(3, 2 + length(group_col_names))
    
    DT::renderDataTable(
      DT::datatable(
        df,
        rownames   = FALSE,
        selection  = "none",
        extensions = "Buttons",
        options    = list(
          pageLength = nrow(df),
          dom        = "Bft",
          scrollX    = TRUE,
          scrollY    = "350px",
          buttons    = list(
            list(extend   = "csv",
                 filename = paste0("SplitChord_LR_", format(Sys.time(), "%Y%m%d")),
                 text     = "CSV"),
            list(extend   = "excel",
                 filename = paste0("SplitChord_LR_", format(Sys.time(), "%Y%m%d")),
                 text     = "Excel",
                 title    = paste("Split Chord — LR pairs |",
                                  paste(sapply(seq_along(split_groups), function(i) {
                                    gd <- if (!is.null(groups_def)) groups_def[[i]] else NULL
                                    if (!is.null(gd))
                                      paste0(split_groups[i], ": ",
                                             paste(gd$sources, collapse = "+"),
                                             " \u2192 ",
                                             paste(gd$targets, collapse = "+"))
                                    else split_groups[i]
                                  }), collapse = "  |  "))
            )
          ),
          columnDefs = list(
            list(className = "dt-center", targets = prob_col_indices - 1)
          )
        )
      ) |>
        DT::formatStyle(
          "Category",
          backgroundColor = DT::styleEqual(all_labels, unname(cat_colors[all_labels])),
          color      = "white",
          fontWeight = "bold"
        ) |>
        DT::formatStyle(
          group_col_names,
          backgroundColor = DT::styleInterval(0, c("#555555", "#1a5276")),
          color = "white"
        ) |>
        DT::formatStyle(
          group_col_names,
          background = DT::styleColorBar(c(0, 0.3), "#5dade2"),
          backgroundSize   = "98% 60%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ),
      server = FALSE
    )
  })
}