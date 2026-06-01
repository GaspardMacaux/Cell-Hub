############################## Differential expressed genes along the trajectory ##############################

# trajectory_server.R


trajectory_server <- function(input, output, session) {

  shinyjs::useShinyjs() # button deactivation

  ############################## Monocle Conversion and trajectory ##############################

  # Reactive values
  seurat_monocle <- reactiveVal()
  monocle_object <- reactiveVal()
  selected_cell_id <- reactiveVal(NULL)
  trajectory_plot <- reactiveVal()
  pseudotime_plot <- reactiveVal()
  # Storage for multiple pseudotime calculations
  pseudotime_history <- reactiveVal(list())
  selected_pseudotime_id <- reactiveVal(NULL)
  
  # Reactive values to store results
  pseudotime_diff_gene_results <- reactiveVal(NULL)  # Pour graph_test
  correlation_results <- reactiveVal(NULL)  # Pour correlation
  gene_trajectory_plot <- reactiveVal(NULL)
  pr_deg_ids <- reactiveVal(NULL)
  
  # Disable buttons initially
  shinyjs::disable("convertToMonocle")
  shinyjs::disable("constructGraph")
  shinyjs::disable("startRootSelection")

  ############################# Monocle Conversion Tab ##############################
  
  output$monocle_cells_count <- renderText({
    # Prefer seurat object for display, fallback to monocle object on direct load
    if (!is.null(seurat_monocle())) {
      formatC(ncol(seurat_monocle()), format = "d", big.mark = ",")
    } else if (!is.null(monocle_object())) {
      formatC(ncol(monocle_object()), format = "d", big.mark = ",")
    } else {
      "—"
    }
  })
  
  output$monocle_clusters_names <- renderText({
    # Prefer seurat object for display, fallback to monocle colData on direct load
    if (!is.null(seurat_monocle())) {
      if (!is.null(Idents(seurat_monocle()))) {
        clusters <- unique(Idents(seurat_monocle()))
        return(paste(as.character(clusters), collapse = ", "))
      } else if ("seurat_clusters" %in% colnames(seurat_monocle()@meta.data)) {
        clusters <- unique(seurat_monocle()$seurat_clusters)
        return(paste(as.character(clusters), collapse = ", "))
      }
    } else if (!is.null(monocle_object())) {
      cds <- monocle_object()
      cluster_col <- tryCatch(GetClusterColumn(cds), error = function(e) NULL)
      if (!is.null(cluster_col)) {
        clusters <- sort(unique(colData(cds)[[cluster_col]]))
        return(paste(as.character(clusters), collapse = ", "))
      }
    }
    return("N/A")
  })
  
  # Helper function to get cluster column name
  GetClusterColumn <- function(cds) {
    if ("ClusterIdents" %in% colnames(colData(cds))) {
      return("ClusterIdents")
    } else if ("seurat_clusters" %in% colnames(colData(cds))) {
      return("seurat_clusters")
    } else {
      stop("No cluster column found in object")
    }
  }
  
  # Helper function to validate pseudotime
  ValidatePseudotime <- function(cds) {
    # Check if pseudotime column exists
    if (!"pseudotime" %in% colnames(colData(cds))) {
      return(list(valid = FALSE, message = "Pseudotime not calculated. Please select a root cell first."))
    }
    
    # Check if values are valid
    pt_values <- pseudotime(cds)
    if (all(is.na(pt_values))) {
      return(list(valid = FALSE, message = "Pseudotime values are all NA. Please re-run root cell selection."))
    }
    
    # Check for sufficient valid values
    n_valid <- sum(!is.na(pt_values))
    if (n_valid < 10) {
      return(list(valid = FALSE, message = paste("Only", n_valid, "cells have valid pseudotime. This is insufficient for analysis.")))
    }
    
    return(list(valid = TRUE, values = pt_values))
  }
  
  

  

  
  # HELPER function to generate unique pseudotime ID
  GeneratePseudotimeID <- function(root_info) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    paste0("PT_", root_info, "_", timestamp)
  }
  
  # HELPER function to extract pseudotime info for display
  FormatPseudotimeInfo <- function(pt_entry) {
    paste0(
      pt_entry$label,
      " (", format(pt_entry$timestamp, "%Y-%m-%d %H:%M"), ")"
    )
  }
  
  # ---- Helper: reset all trajectory reactive state ----
  ResetTrajectoryState <- function() {
    monocle_object(NULL)
    seurat_monocle(NULL)
    trajectory_plot(NULL)
    pseudotime_plot(NULL)
    pseudotime_history(list())
    selected_pseudotime_id(NULL)
    selected_cell_id(NULL)
    
    # Reset DEG analysis results
    pseudotime_diff_gene_results(NULL)
    correlation_results(NULL)
    
    # Reset all plot outputs
    output$trajectoryPlot  <- renderPlot({ NULL })
    output$pseudotimePlot  <- renderPlot({ NULL })
    
    # Reset DEG table outputs
    output$graphTestTable    <- renderDT({ NULL })
    output$correlationTable  <- renderDT({ NULL })
    
    # Reset gene pickers in visualization tab
    updatePickerInput(session, "multi_gene_pseudotime_trajectory",
                      choices = character(0), selected = NULL)
    updateSelectInput(session, "single_gene_umap_trajectory",
                      choices = character(0), selected = NULL)
    updatePickerInput(session, "gene_picker_corr",
                      choices = character(0), selected = NULL)
    
    # Reset selectors
    updateSelectInput(session, "select_pseudotime",
                      choices = character(0), selected = NULL)
    updateSelectInput(session, "root_cluster_select",
                      choices = character(0))
    
    # Disable action buttons
    shinyjs::disable("convertToMonocle")
    shinyjs::disable("constructGraph")
    shinyjs::disable("set_root_cell")
  }
  
  # ---- Load Seurat object for Monocle conversion ----
  observeEvent(input$load_seurat_file_monocle, {
    message("Attempting to read Seurat file at: ", input$load_seurat_file_monocle$datapath)
    
    # Reset all previous state before loading
    ResetTrajectoryState()
    
    tryCatch({
      temp_seurat <- readRDS(input$load_seurat_file_monocle$datapath)
      
      if (!inherits(temp_seurat, "Seurat")) {
        showNotification("The uploaded file does not contain a valid Seurat object.", type = "error")
        return()
      }
      
      message("Seurat file successfully read.")
      seurat_monocle(temp_seurat)
      shinyjs::enable("convertToMonocle")
      showNotification(
        paste0("Seurat object loaded: ",
               formatC(ncol(temp_seurat), format = "d", big.mark = ","), " cells, ",
               length(unique(Idents(temp_seurat))), " clusters."),
        type = "message"
      )
    }, error = function(e) {
      message("Error loading Seurat object: ", e$message)
      showNotification(paste("Error loading Seurat object:", e$message), type = "error")
    })
  })
  
  # ---- Load previously saved Monocle RDS ----
  observeEvent(input$load_monocle_rds, {
    req(input$load_monocle_rds)
    
    # Reset all previous state before loading
    ResetTrajectoryState()
    
    showModal(modalDialog(
      title = "Loading Monocle Object",
      "Reading RDS file and restoring Monocle3 object...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      loaded_obj <- readRDS(input$load_monocle_rds$datapath)
      
      if (!inherits(loaded_obj, "cell_data_set")) {
        removeModal()
        showNotification(
          "The uploaded file does not contain a valid Monocle3 cell_data_set object.",
          type = "error"
        )
        return()
      }
      
      monocle_object(loaded_obj)
      
      # Detect cluster column
      cluster_column <- tryCatch(GetClusterColumn(loaded_obj), error = function(e) NULL)
      
      # Update cluster selector
      if (!is.null(cluster_column)) {
        clusters <- sort(unique(colData(loaded_obj)[[cluster_column]]))
        updateSelectInput(session, "root_cluster_select",
                          choices = clusters,
                          selected = clusters[1])
      }
      
      has_graph <- !is.null(loaded_obj@principal_graph) &&
        length(loaded_obj@principal_graph) > 0
      
      has_pseudotime <- "pseudotime" %in% colnames(colData(loaded_obj)) &&
        !all(is.na(pseudotime(loaded_obj)))
      
      # Rebuild trajectory plot if graph exists
      if (has_graph && !is.null(cluster_column)) {
        traj_plot <- plot_cells(
          loaded_obj,
          color_cells_by = cluster_column,
          cell_size = 0.5,
          label_groups_by_cluster = TRUE,
          label_leaves = FALSE,
          label_branch_points = TRUE,
          graph_label_size = 1.5
        ) +
          ggtitle("Trajectory (restored from saved object)") +
          theme_minimal()
        
        trajectory_plot(traj_plot)
        output$trajectoryPlot <- renderPlot({ trajectory_plot() })
        shinyjs::enable("set_root_cell")
      } else {
        shinyjs::enable("constructGraph")
      }
      
      # Rebuild pseudotime plot and history if pseudotime exists
      if (has_pseudotime) {
        pseudo_plot <- plot_cells(
          loaded_obj,
          color_cells_by = "pseudotime",
          cell_size = 0.5,
          label_cell_groups = FALSE,
          label_leaves = FALSE,
          label_branch_points = TRUE,
          graph_label_size = 1.5
        ) +
          ggtitle("Pseudotime (restored from saved object)") +
          theme_minimal() +
          viridis::scale_color_viridis(name = "Pseudotime", option = "plasma")
        
        pseudotime_plot(pseudo_plot)
        output$pseudotimePlot <- renderPlot({ pseudotime_plot() })
        
        # Create a minimal history entry so the picker is populated
        pt_id <- paste0("pt_restored_", format(Sys.time(), "%H%M%S"))
        entry <- list(
          label          = "Restored pseudotime",
          root_cluster   = "unknown (restored)",
          timestamp      = Sys.time(),
          monocle_object = loaded_obj
        )
        new_history <- list()
        new_history[[pt_id]] <- entry
        pseudotime_history(new_history)
        selected_pseudotime_id(pt_id)
        
        pt_choices <- sapply(new_history, function(x) x$label)
        names(pt_choices) <- names(new_history)
        updateSelectInput(session, "select_pseudotime",
                          choices = pt_choices,
                          selected = pt_id)
      }
      
      removeModal()
      showNotification(
        paste0(
          "Monocle object loaded: ",
          formatC(ncol(loaded_obj), format = "d", big.mark = ","), " cells",
          if (has_graph) " | trajectory graph present" else " | no graph",
          if (has_pseudotime) " | pseudotime present" else " | no pseudotime"
        ),
        type = "message",
        duration = 6
      )
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error loading Monocle object:", e$message), type = "error")
      message("Monocle load error: ", e$message)
    })
  })

  # Initialize gene loadings
  gene_loadings <- NULL


  initialize_gene_loadings <- function(monocle_obj, gene_loadings_matrix) {
    gene_names <- rowData(monocle_obj)$gene_short_name
    if (!is.null(gene_names)) {
      common_genes <- intersect(gene_names, rownames(gene_loadings_matrix))
      if (length(common_genes) == nrow(gene_loadings_matrix)) {
        rownames(gene_loadings_matrix) <- common_genes
        gene_loadings <<- gene_loadings_matrix
      } else {
        showNotification("Mismatch between gene names and gene loadings matrix.", type = "error")
      }
    } else {
      showNotification("Gene names not found in Monocle object.", type = "error")
    }
  }
  
  
  
  # ---- Save Monocle object as RDS ----
  output$download_monocle_object <- downloadHandler(
    filename = function() {
      paste0("monocle_object_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(monocle_object())
      
      showModal(modalDialog(
        title = "Saving Monocle Object",
        "Saving your Monocle3 object with trajectory and pseudotime...",
        footer = NULL,
        easyClose = FALSE
      ))
      
      tryCatch({
        monocle_obj <- monocle_object()
        
        # Attach export metadata in misc for traceability
        monocle_obj@metadata$cell_hub_export <- list(
          export_date    = Sys.time(),
          export_module  = "Trajectory - Monocle3",
          n_cells        = ncol(monocle_obj),
          n_genes        = nrow(monocle_obj),
          has_trajectory = !is.null(monocle_obj@principal_graph),
          pseudotime_ids = names(pseudotime_history())
        )
        
        saveRDS(monocle_obj, file = file)
        removeModal()
        showNotification("Monocle object saved successfully!", type = "message")
        
      }, error = function(e) {
        removeModal()
        showNotification(paste("Error saving object:", e$message), type = "error")
        message("Monocle save error:", e$message)
      })
    }
  )
  
 

  # FIXED: Monocle conversion
  observeEvent(input$convertToMonocle, {
    req(seurat_monocle())
    
    showModal(modalDialog(
      title = "Converting to Monocle3",
      "Converting Seurat object to Monocle format...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      seurat_obj <- seurat_monocle()
      message("Starting Seurat to Monocle conversion...")
      
      # Get expression data
      expression_matrix <- GetAssayData(seurat_obj, layer = "counts", assay = "RNA")
      message(paste("Expression matrix:", nrow(expression_matrix), "genes x", ncol(expression_matrix), "cells"))
      
      # Prepare cell metadata - ensure ClusterIdents is present
      cell_metadata <- seurat_obj@meta.data
      if (!"ClusterIdents" %in% colnames(cell_metadata)) {
        if (!is.null(Idents(seurat_obj))) {
          cell_metadata$ClusterIdents <- Idents(seurat_obj)
        } else if ("seurat_clusters" %in% colnames(cell_metadata)) {
          cell_metadata$ClusterIdents <- cell_metadata$seurat_clusters
        } else {
          stop("No cluster information found in Seurat object")
        }
      }
      
      # Prepare gene metadata
      gene_metadata <- data.frame(
        gene_short_name = rownames(expression_matrix),
        row.names = rownames(expression_matrix)
      )
      
      # Create Monocle object
      message("Creating cell_data_set object...")
      monocle_temp <- new_cell_data_set(
        expression_data = expression_matrix,
        cell_metadata = cell_metadata,
        gene_metadata = gene_metadata
      )
      
      # Transfer reductions - with dimension validation
      message("Transferring reductions...")
      
      # PCA
      if ("pca" %in% names(seurat_obj@reductions)) {
        pca_embeddings <- Embeddings(seurat_obj, "pca")
        if (nrow(pca_embeddings) == ncol(monocle_temp)) {
          reducedDims(monocle_temp)[["PCA"]] <- pca_embeddings
          message(paste("PCA transferred:", ncol(pca_embeddings), "components"))
        } else {
          warning("PCA dimensions don't match. Skipping PCA transfer.")
        }
      }
      
      # UMAP - CRITICAL for trajectory
      if ("umap" %in% names(seurat_obj@reductions)) {
        umap_embeddings <- Embeddings(seurat_obj, "umap")
        if (nrow(umap_embeddings) == ncol(monocle_temp)) {
          reducedDims(monocle_temp)[["UMAP"]] <- umap_embeddings
          message(paste("UMAP transferred:", nrow(umap_embeddings), "cells"))
        } else {
          stop("UMAP dimensions don't match cell count. Cannot proceed with trajectory analysis.")
        }
      } else {
        stop("No UMAP found in Seurat object. Run UMAP in Seurat first.")
      }
      
      # Estimate size factors
      message("Calculating size factors...")
      monocle_temp <- estimate_size_factors(monocle_temp)
      
      # Set up cluster structure WITHOUT re-clustering
      message("Setting up cluster partitions...")
      
      # Create cluster structure directly from Seurat clusters
      cell_names <- colnames(monocle_temp)
      cluster_factor <- factor(cell_metadata$ClusterIdents)
      names(cluster_factor) <- cell_names
      
      # Create single partition (all cells in one partition for trajectory)
      single_partition <- factor(rep("1", ncol(monocle_temp)))
      names(single_partition) <- cell_names
      
      # Initialize clusters slot properly
      monocle_temp@clusters <- SimpleList()
      monocle_temp@clusters[["UMAP"]] <- SimpleList(
        partitions = single_partition,
        clusters = cluster_factor
      )
      
      # Validation
      message("Validating converted object...")
      message(paste("Monocle cells:", ncol(monocle_temp)))
      message(paste("UMAP dimensions:", nrow(reducedDims(monocle_temp)[["UMAP"]]), "x", 
                    ncol(reducedDims(monocle_temp)[["UMAP"]])))
      message(paste("Clusters:", paste(unique(colData(monocle_temp)$ClusterIdents), collapse = ", ")))
      
      # Store the object
      monocle_object(monocle_temp)
      
      # Enable next button
      shinyjs::enable("constructGraph")
      removeModal()
      
      showNotification("Conversion successful! Ready to construct trajectory.", type = "message")
      
    }, error = function(e) {
      removeModal()
      message("Error during conversion:", e$message)
      showNotification(paste("Error:", conditionMessage(e)[1]), type = "error")
      })
  })
  

  # Construct trajectory graph
  observeEvent(input$constructGraph, {
    req(monocle_object())
    
    showModal(modalDialog(
      title = "Constructing Graph",
      "Creating trajectory using existing Seurat clusters...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      # Get the monocle object
      monocle_temp <- monocle_object()
      
      message("Starting graph construction...")
      message(paste("Cells:", ncol(monocle_temp), "Genes:", nrow(monocle_temp)))
      
      # Check which cluster column to use
      cluster_column <- NULL
      if ("ClusterIdents" %in% colnames(colData(monocle_temp))) {
        cluster_column <- "ClusterIdents"
      } else if ("seurat_clusters" %in% colnames(colData(monocle_temp))) {
        cluster_column <- "seurat_clusters"
      } else {
        showNotification("No cluster information found", type = "error")
        removeModal()
        return()
      }
      
      message(paste("Using cluster column:", cluster_column))
      cluster_data <- colData(monocle_temp)[[cluster_column]]
      message(paste("Unique clusters:", paste(unique(cluster_data), collapse = ", ")))
      
      # IMPORTANT: Set up clusters structure properly before graph construction
      cell_names <- colnames(monocle_temp)
      
      # Create single partition
      single_partition <- factor(rep("1", ncol(monocle_temp)))
      names(single_partition) <- cell_names
      
      # Set up clusters with names
      cluster_factor <- factor(cluster_data)
      names(cluster_factor) <- cell_names
      
      # Initialize clusters slot
      monocle_temp@clusters <- SimpleList()
      monocle_temp@clusters[["UMAP"]] <- SimpleList(
        partitions = single_partition,
        clusters = cluster_factor
      )
      
      message("Cluster structure initialized")
      
      # Construct the trajectory graph
      monocle_temp <- learn_graph(
        monocle_temp,
        use_partition = FALSE,
        close_loop = FALSE
      )
      
      message("Graph construction completed")
      
      # Update the monocle object
      monocle_object(monocle_temp)
      
      # Create trajectory plot
      base_plot <- plot_cells(
        monocle_temp,
        color_cells_by = cluster_column,
        label_groups_by_cluster = TRUE,
        label_leaves = FALSE,
        label_branch_points = TRUE,
        cell_size = 0.5
      ) + 
        ggtitle("Trajectory on Clusters") +
        theme_minimal()
      
      # Store and display plot
      trajectory_plot(base_plot)
      output$trajectoryPlot <- renderPlot({
        base_plot
      })
      
      # Update cluster choices for root selection
      clusters <- sort(unique(colData(monocle_temp)[[cluster_column]]))
      updateSelectInput(session, "root_cluster_select", 
                        choices = clusters,
                        selected = clusters[1])
      
      # Enable root selection button
      shinyjs::enable("set_root_cell")
      
      removeModal()
      showNotification("Graph construction completed! Now select a root cluster.", type = "message")
      
    }, error = function(e) {
      removeModal()
      message("Error during graph construction:", e$message)
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # MODIFIER ton observeEvent set_root_cell - À LA FIN, ajoute juste:
  observeEvent(input$set_root_cell, {
    req(monocle_object(), input$root_cluster_select)
    
    monocle_data <- monocle_object()
    
    # Determine cluster column
    cluster_column <- if("ClusterIdents" %in% colnames(colData(monocle_data))) {
      "ClusterIdents"
    } else {
      "seurat_clusters"
    }
    
    # Get cells in selected cluster
    cells_in_cluster <- which(colData(monocle_data)[[cluster_column]] == input$root_cluster_select)
    cell_barcodes <- colnames(monocle_data)[cells_in_cluster]
    
    if (length(cells_in_cluster) == 0) {
      showNotification("No cells found in selected cluster", type = "error")
      return()
    }
    
    message(paste("Selected cluster:", input$root_cluster_select, "with", length(cells_in_cluster), "cells"))
    
    # Find center cell of the cluster
    umap_coords <- reducedDims(monocle_data)[["UMAP"]][cells_in_cluster, , drop = FALSE]
    center <- colMeans(umap_coords)
    distances <- sqrt(rowSums(sweep(umap_coords, 2, center, "-")^2))
    center_idx <- which.min(distances)
    root_cell <- cell_barcodes[center_idx]
    
    selected_cell_id(root_cell)
    
    # Order cells with pseudotime
    showModal(modalDialog(
      title = "Calculating Pseudotime",
      paste("Computing pseudotime from cluster", input$root_cluster_select, "..."),
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      message(paste("Ordering cells with root:", root_cell))
      
      monocle_temp <- order_cells(
        monocle_data,
        reduction_method = "UMAP",
        root_cells = c(root_cell)
      )
      
      message("Cell ordering completed")
      
      # CRITICAL: Extract pseudotime et l'ajouter à colData
      pt_values <- pseudotime(monocle_temp)
      
      if (is.null(pt_values) || all(is.na(pt_values))) {
        removeModal()
        showNotification("Error: Could not calculate pseudotime", type = "error")
        return()
      }
      
      # AJOUT CRITIQUE: Ajouter pseudotime dans colData
      colData(monocle_temp)$pseudotime <- pt_values
      
      # Vérifier que c'est bien ajouté
      if ("pseudotime" %in% colnames(colData(monocle_temp))) {
        message("SUCCESS: Pseudotime added to colData")
      } else {
        message("ERROR: Failed to add pseudotime to colData")
      }
      
      message(paste("Pseudotime range:", round(min(pt_values, na.rm = TRUE), 2), 
                    "to", round(max(pt_values, na.rm = TRUE), 2)))
      message(paste("Cells with valid pseudotime:", sum(!is.na(pt_values))))
      
      # Stocker le pseudotime dans l'historique
      pt_id <- paste0("Root_Cluster_", input$root_cluster_select, "_", format(Sys.time(), "%H%M%S"))
      
      pt_entry <- list(
        id = pt_id,
        label = paste("Root: Cluster", input$root_cluster_select),
        root_cluster = input$root_cluster_select,
        root_cell = root_cell,
        monocle_object = monocle_temp,
        timestamp = Sys.time()
      )
      
      # Ajouter à l'historique
      current_history <- pseudotime_history()
      current_history[[pt_id]] <- pt_entry
      pseudotime_history(current_history)
      selected_pseudotime_id(pt_id)
      
      # Mettre à jour le selectInput
      pt_choices <- sapply(current_history, function(x) x$label)
      names(pt_choices) <- names(current_history)
      
      updateSelectInput(
        session,
        "select_pseudotime",
        choices = pt_choices,
        selected = pt_id
      )
      
      # Update the main monocle object
      monocle_object(monocle_temp)
      
      # Create plots
      traj_plot <- plot_cells(
        monocle_temp,
        color_cells_by = cluster_column,
        cell_size = 0.5,
        label_groups_by_cluster = TRUE,
        label_leaves = FALSE,
        label_branch_points = TRUE,
        graph_label_size = 1.5
      ) +
        ggtitle(paste("Trajectory - Root: Cluster", input$root_cluster_select)) +
        theme_minimal()
      
      pseudo_plot <- plot_cells(
        monocle_temp,
        color_cells_by = "pseudotime",
        cell_size = 0.5,
        label_cell_groups = FALSE,
        label_leaves = FALSE,
        label_branch_points = TRUE,
        graph_label_size = 1.5
      ) +
        ggtitle("Pseudotime Trajectory") +
        theme_minimal() +
        viridis::scale_color_viridis(
          name = "Pseudotime",
          option = "plasma"
        )
      
      trajectory_plot(traj_plot)
      pseudotime_plot(pseudo_plot)
      
      output$trajectoryPlot <- renderPlot({ trajectory_plot() })
      output$pseudotimePlot <- renderPlot({ pseudotime_plot() })
      
      output$root_info <- renderText({
        paste("Root cluster:", input$root_cluster_select,
              "\nNumber of cells:", length(cells_in_cluster),
              "\nRoot cell:", root_cell)
      })
      
      removeModal()
      showNotification("Pseudotime calculated successfully!", type = "message")
      
    }, error = function(e) {
      removeModal()
      message("Error in cell ordering:", e$message)
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  # Download handlers
  output$download_trajectory_umap <- downloadHandler(
    filename = function() {
      paste("trajectory_plot_", Sys.Date(), ".", input$trajectory_download_format, sep = "")
    },
    content = function(file) {
      tryCatch({
        req(trajectory_plot())
        ggsave(
          file,
          plot = trajectory_plot(),
          device = input$trajectory_download_format,
          dpi = input$trajectory_download_dpi,
          width = 10,
          height = 8
        )
      }, error = function(e) {
        showNotification(paste("Error saving plot:", e$message), type = "error")
      })
    }
  )

  output$download_pseudotime_umap <- downloadHandler(
    filename = function() {
      paste("pseudotime_plot_", Sys.Date(), ".", input$trajectory_download_format, sep = "")
    },
    content = function(file) {
      tryCatch({
        req(pseudotime_plot())
        ggsave(
          file,
          plot = pseudotime_plot(),
          device = input$trajectory_download_format,
          dpi = input$trajectory_download_dpi,
          width = 10,
          height = 8
        )
      }, error = function(e) {
        showNotification(paste("Error saving plot:", e$message), type = "error")
      })
    }
  )

  # NOUVEAU: Observer pour changer de pseudotime quand on sélectionne dans le picker
  observeEvent(input$select_pseudotime, {
    req(input$select_pseudotime, pseudotime_history())
    
    pt_id <- input$select_pseudotime
    history <- pseudotime_history()
    
    if (!pt_id %in% names(history)) return()
    
    pt_entry <- history[[pt_id]]
    
    # Mettre à jour le monocle object
    monocle_object(pt_entry$monocle_object)
    selected_pseudotime_id(pt_id)
    
    # Mettre à jour les plots
    cluster_column <- if("ClusterIdents" %in% colnames(colData(pt_entry$monocle_object))) {
      "ClusterIdents"
    } else {
      "seurat_clusters"
    }
    
    traj_plot <- plot_cells(
      pt_entry$monocle_object,
      color_cells_by = cluster_column,
      cell_size = 0.5,
      label_groups_by_cluster = TRUE,
      label_leaves = FALSE,
      label_branch_points = TRUE,
      graph_label_size = 1.5
    ) +
      ggtitle(paste("Trajectory -", pt_entry$label)) +
      theme_minimal()
    
    pseudo_plot <- plot_cells(
      pt_entry$monocle_object,
      color_cells_by = "pseudotime",
      cell_size = 0.5,
      label_cell_groups = FALSE,
      label_leaves = FALSE,
      label_branch_points = TRUE,
      graph_label_size = 1.5
    ) +
      ggtitle(paste("Pseudotime -", pt_entry$label)) +
      theme_minimal() +
      viridis::scale_color_viridis(
        name = "Pseudotime",
        option = "plasma"
      )
    
    trajectory_plot(traj_plot)
    pseudotime_plot(pseudo_plot)
    
    output$trajectoryPlot <- renderPlot({ trajectory_plot() })
    output$pseudotimePlot <- renderPlot({ pseudotime_plot() })
    
    output$root_info <- renderText({
      paste("Selected:", pt_entry$label,
            "\nRoot cell:", pt_entry$root_cell,
            "\nCalculated:", format(pt_entry$timestamp, "%Y-%m-%d %H:%M"))
    })
    
    showNotification(paste("Switched to:", pt_entry$label), type = "message", duration = 3)
  })
  
  ################## Calculate differential expression/TAB 2################
  

  
  
  
  
  CalculateDiffGenesPseudotime <- function() {
    req(monocle_object())
    
    monocle_obj <- monocle_object()
    
    # Check if trajectory is calculated
    if (is.null(monocle_obj@principal_graph) || length(monocle_obj@principal_graph) == 0) {
      showNotification("Trajectory not calculated. Complete trajectory analysis first.", type = "error")
      return(NULL)
    }
    
    # Check pseudotime exists (même si on ne l'utilise pas directement)
    if (!"pseudotime" %in% colnames(colData(monocle_obj))) {
      showNotification("Pseudotime not calculated. Please select a root cell first.", type = "error")
      return(NULL)
    }
    
    pt_values <- pseudotime(monocle_obj)
    if (all(is.na(pt_values))) {
      showNotification("Pseudotime values are all NA.", type = "error")
      return(NULL)
    }
    
    message("=== DEBUG PSEUDOTIME ===")
    message(paste("Valid values:", sum(!is.na(pt_values))))
    message(paste("Range:", round(min(pt_values, na.rm = TRUE), 2), "-", round(max(pt_values, na.rm = TRUE), 2)))
    message("========================")
    
    tryCatch({
      q_value_threshold <- input$sig_genes_cutoff
      
      message("Starting graph_test analysis (trajectory structure)...")
      message(paste("Q-value threshold:", q_value_threshold))
      
      # Run differential test
      n_cores <- min(4, parallel::detectCores() - 1)
      
      markers <- graph_test(
        monocle_obj, 
        neighbor_graph = "principal_graph", 
        cores = n_cores
      )
      
      if (is.null(markers) || nrow(markers) == 0) {
        showNotification("No results from graph_test.", type = "error")
        return(NULL)
      }
      
      # Log initial results
      message(paste("Total genes tested:", nrow(markers)))
      message(paste("Genes with q < 0.05:", sum(markers$q_value < 0.05, na.rm = TRUE)))
      message(paste("Genes with q < 0.01:", sum(markers$q_value < 0.01, na.rm = TRUE)))
      
      # Filter by q-value
      sig_genes <- markers[markers$q_value < q_value_threshold, ]
      
      if (nrow(sig_genes) == 0) {
        min_q <- min(markers$q_value, na.rm = TRUE)
        showNotification(
          paste0("No genes with q-value < ", q_value_threshold,
                 ". Minimum q-value found: ", format(min_q, scientific = TRUE)),
          type = "warning",
          duration = 10
        )
        return(NULL)
      }
      
      # Sort by morans_I
      sig_genes <- sig_genes[order(-sig_genes$morans_I), ]
      
      message(paste("Found", nrow(sig_genes), "genes with trajectory structure"))
      message(paste("Top 5 genes:", paste(head(rownames(sig_genes), 5), collapse = ", ")))
      
      # Store results
      pseudotime_diff_gene_results(sig_genes)
      
      return(sig_genes)
      
    }, error = function(e) {
      message("Error in graph_test:", e$message)
      showNotification(paste("Error:", e$message), type = "error")
      return(NULL)
    })
  }
  
  # Function 2: Correlation with pseudotime (dépend du pseudotime)
  CalculateCorrelationWithPseudotime <- function() {
    req(monocle_object())
    
    monocle_obj <- monocle_object()
    
    # Check pseudotime
    if (!"pseudotime" %in% colnames(colData(monocle_obj))) {
      showNotification("Pseudotime not calculated.", type = "error")
      return(NULL)
    }
    
    pt_values <- pseudotime(monocle_obj)
    if (all(is.na(pt_values))) {
      showNotification("Pseudotime values are all NA.", type = "error")
      return(NULL)
    }
    
    tryCatch({
      message("=== CORRELATION ANALYSIS ===")
      message("Testing correlation with pseudotime for each gene...")
      
      # Get expression matrix
      expr_matrix <- exprs(monocle_obj)
      
      # Get valid cells (with pseudotime)
      valid_cells <- !is.na(pt_values)
      pt_valid <- pt_values[valid_cells]
      expr_valid <- expr_matrix[, valid_cells]
      
      message(paste("Testing", nrow(expr_matrix), "genes across", sum(valid_cells), "cells"))
      
      # Test correlation for each gene
      results_list <- lapply(1:nrow(expr_valid), function(i) {
        gene_expr <- expr_valid[i, ]
        
        # Skip lowly expressed genes
        if (sum(gene_expr > 0) < 10) return(NULL)
        
        # Spearman correlation with exact=FALSE to avoid ties warning
        test <- cor.test(gene_expr, pt_valid, method = "spearman", exact = FALSE)
        
        data.frame(
          gene = rownames(expr_valid)[i],
          correlation = as.numeric(test$estimate),
          p_value = test$p.value,
          stringsAsFactors = FALSE
        )
      })
      
      # Combine results (remove NULL)
      results <- do.call(rbind, results_list[!sapply(results_list, is.null)])
      
      # Remove genes with NA
      results <- results[!is.na(results$correlation), ]
      
      # Adjust p-values
      results$q_value <- p.adjust(results$p_value, method = "BH")
      
      # Filter by significance
      q_threshold <- input$sig_genes_cutoff_corr
      sig_results <- results[results$q_value < q_threshold, ]
      
      message(paste("Total genes tested:", nrow(results)))
      message(paste("Significant genes (q < ", q_threshold, "):", nrow(sig_results)))
      
      if (nrow(sig_results) == 0) {
        showNotification("No genes significantly correlated with pseudotime.", type = "warning")
        return(NULL)
      }
      
      # Separate up and down
      up_genes <- sig_results[sig_results$correlation > 0, ]
      down_genes <- sig_results[sig_results$correlation < 0, ]
      
      message(paste("Up-regulated genes:", nrow(up_genes)))
      message(paste("Down-regulated genes:", nrow(down_genes)))
      
      # Sort by absolute correlation
      sig_results <- sig_results[order(-abs(sig_results$correlation)), ]
      rownames(sig_results) <- sig_results$gene
      
      message(paste("Top 5 correlated genes:", paste(head(sig_results$gene, 5), collapse = ", ")))
      message("============================")
      
      # Store results
      correlation_results(sig_results)
      
      # Update gene picker for correlation genes
      updatePickerInput(
        session, 
        "gene_picker_corr",
        choices = sig_results$gene,
        selected = head(sig_results$gene, 10)
      )
      
      return(sig_results)
      
    }, error = function(e) {
      message("Error in correlation analysis:", e$message)
      showNotification(paste("Error:", e$message), type = "error")
      return(NULL)
    })
  }
  
  # ObserveEvent 1: graph_test (lancé une seule fois, indépendant du pseudotime)
  observeEvent(input$run_graph_test, {
    req(monocle_object())
    
    showModal(modalDialog(
      title = "Analyzing Trajectory Structure",
      "Identifying genes with spatial patterns in the trajectory...",
      "This is independent of root cell selection.",
      footer = NULL,
      easyClose = FALSE
    ))
    
    results <- CalculateDiffGenesPseudotime()
    
    removeModal()
    
    if (!is.null(results) && nrow(results) > 0) {
      
      # Afficher le tableau
      output$graphTestTable <- renderDT({
        
        display_data <- as.data.frame(results)
        
        # Ajouter gene name
        if (!("gene_short_name" %in% colnames(display_data))) {
          display_data$gene_short_name <- rownames(display_data)
        }
        
        # Colonnes à afficher
        cols <- c("gene_short_name", "morans_I", "q_value", "p_value", "status")
        cols <- cols[cols %in% colnames(display_data)]
        display_data <- display_data[, cols]
        
        colnames(display_data) <- c("Gene", "Moran's I", "Q-value", "P-value", "Status")[1:ncol(display_data)]
        
        datatable(
          display_data,
          options = list(
            pageLength = 25,
            scrollX = TRUE,
            order = list(list(1, 'desc')),  # Sort by Moran's I
            dom = 'Bfrtip',
            buttons = c('copy', 'csv')
          ),
          rownames = FALSE,
          filter = 'top',
          caption = "Genes with spatial structure in trajectory (independent of root)"
        ) %>%
          formatRound(columns = "Moran's I", digits = 4) %>%
          formatSignif(columns = c("Q-value", "P-value"), digits = 3) %>%
          formatStyle(
            "Moran's I",
            color = styleInterval(c(0.3), c('black', 'darkred')),
            fontWeight = styleInterval(c(0.3), c('normal', 'bold'))
          ) %>%
          formatStyle(
            "Q-value",
            backgroundColor = styleInterval(c(0.01, 0.05), c('lightgreen', 'lightyellow', 'white'))
          )
      })
      
      showNotification(
        paste("Found", nrow(results), "genes with trajectory structure!"), 
        type = "message",
        duration = 5
      )
    }
  })
 
  
  
  
  # ObserveEvent 2: Correlation (dépend du pseudotime sélectionné)
  observeEvent(input$run_correlation_test, {
    req(monocle_object(), selected_pseudotime_id())
    
    # Get pseudotime info
    pt_id <- selected_pseudotime_id()
    history <- pseudotime_history()
    pt_entry <- history[[pt_id]]
    
    showModal(modalDialog(
      title = "Correlation Analysis",
      paste("Testing genes correlated with:", pt_entry$label),
      "Finding up/down-regulated genes along this specific pseudotime...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    results <- CalculateCorrelationWithPseudotime()
    
    removeModal()
    
    if (!is.null(results) && nrow(results) > 0) {
      
      # Save in history
      if (pt_id %in% names(history)) {
        history[[pt_id]]$correlation_results <- results
        history[[pt_id]]$correlation_timestamp <- Sys.time()
        pseudotime_history(history)
      }
      
      # Afficher le tableau
      output$correlationTable <- renderDT({
        
        display_data <- results[, c("gene", "correlation", "p_value", "q_value")]
        colnames(display_data) <- c("Gene", "Correlation", "P-value", "Q-value")
        
        # Add direction column
        display_data$Direction <- ifelse(display_data$Correlation > 0, "Up", "Down")
        display_data <- display_data[, c("Gene", "Correlation", "Direction", "P-value", "Q-value")]
        
        datatable(
          display_data,
          options = list(
            pageLength = 25,
            scrollX = TRUE,
            order = list(list(1, 'desc')),  # Sort by correlation
            dom = 'Bfrtip',
            buttons = c('copy', 'csv')
          ),
          rownames = FALSE,
          filter = 'top',
          caption = paste("Genes correlated with", pt_entry$label)
        ) %>%
          formatRound(columns = "Correlation", digits = 3) %>%
          formatSignif(columns = c("P-value", "Q-value"), digits = 3) %>%
          formatStyle(
            "Correlation",
            color = styleInterval(c(-0.3, 0.3), c('blue', 'gray', 'red')),
            fontWeight = styleInterval(c(-0.5, 0.5), c('bold', 'normal', 'bold'))
          ) %>%
          formatStyle(
            "Direction",
            backgroundColor = styleEqual(c("Up", "Down"), c('#ffcccc', '#ccccff'))
          ) %>%
          formatStyle(
            "Q-value",
            backgroundColor = styleInterval(c(0.01, 0.05), c('lightgreen', 'lightyellow', 'white'))
          )
      })
      
      showNotification(
        paste("Found", nrow(results), "genes correlated with pseudotime!"), 
        type = "message",
        duration = 5
      )
    }
  })
  
  
  # Download handler pour graph_test results
  output$download_graph_test <- createDownloadHandler(
    reactive_data = pseudotime_diff_gene_results,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(monocle_object(), default_name = "Trajectory") 
    }),
    data_name = "graph_test_trajectory_structure",
    download_type = "csv"
  )
  
  # Download handler pour correlation results
  output$download_correlation <- createDownloadHandler(
    reactive_data = correlation_results,
    object_name_reactive = reactive({
      # Ajouter le nom du root dans le nom de fichier
      pt_id <- selected_pseudotime_id()
      if (!is.null(pt_id)) {
        history <- pseudotime_history()
        pt_entry <- history[[pt_id]]
        root_name <- gsub(" ", "_", pt_entry$root_cluster)
        paste0("Trajectory_Root_", root_name)
      } else {
        "Trajectory"
      }
    }),
    data_name = "correlation_pseudotime",
    download_type = "csv"
  )

  
  
  ######################## TAB 3 / Genes visualization ###########################
  
  

  # Reactive for gene trajectory plots
  gene_trajectory_plot <- reactiveVal(NULL)
  single_gene_trajectory_plot <- reactiveVal(NULL)
  
  # Update gene pickers when results are available
  # Update gene pickers when results are available
  observe({
    # Prioriser les résultats de graph_test
    if (!is.null(pseudotime_diff_gene_results())) {
      genes <- rownames(pseudotime_diff_gene_results())
      
      # Update multi-gene picker - NO PRE-SELECTION
      updatePickerInput(
        session,
        "multi_gene_pseudotime_trajectory",
        choices = genes,
        selected = NULL  # CHANGÉ: aucune présélection
      )
      
      # Update single gene selector - NO PRE-SELECTION
      updateSelectInput(
        session,
        "single_gene_umap_trajectory",
        choices = genes,
        selected = NULL  # CHANGÉ: aucune présélection
      )
    } 
    # Sinon utiliser correlation si disponible
    else if (!is.null(correlation_results())) {
      genes <- correlation_results()$gene
      
      updatePickerInput(
        session,
        "multi_gene_pseudotime_trajectory",
        choices = genes,
        selected = NULL  # CHANGÉ: aucune présélection
      )
      
      updateSelectInput(
        session,
        "single_gene_umap_trajectory",
        choices = genes,
        selected = NULL  # CHANGÉ: aucune présélection
      )
    }
  })
  
  # Plot single gene on UMAP trajectory
  observeEvent(input$plot_single_gene_umap_trajectory, {
    req(monocle_object(), input$single_gene_umap_trajectory)
    
    showModal(modalDialog(
      title = "Plotting Gene on Trajectory",
      "Creating visualization...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      cds <- monocle_object()
      gene <- input$single_gene_umap_trajectory
      
      # Verify gene exists
      if (!gene %in% rownames(cds)) {
        removeModal()
        showNotification("Selected gene not found in dataset", type = "error")
        return()
      }
      
      # Create plot with gene expression on trajectory
      plot <- plot_cells(
        cds,
        genes = gene,
        show_trajectory_graph = input$show_trajectory_graph_single,
        label_cell_groups = FALSE,
        label_leaves = FALSE,
        label_branch_points = input$show_trajectory_graph_single,
        cell_size = input$single_gene_cell_size_trajectory,
        cell_stroke = 0.1,
        trajectory_graph_color = "black",
        trajectory_graph_segment_size = 0.75
      ) +
        ggtitle(paste("Expression of", gene, "along trajectory")) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          legend.position = "right"
        ) +
        viridis::scale_color_viridis(
          name = "Expression",
          option = "magma"
        )
      
      single_gene_trajectory_plot(plot)
      
      output$singleGeneUmapTrajectoryPlot <- renderPlot({
        single_gene_trajectory_plot()
      })
      
      removeModal()
      
      showNotification(
        paste("Successfully plotted", gene, "on trajectory"),
        type = "message",
        duration = 3
      )
      
    }, error = function(e) {
      removeModal()
      message("Error plotting single gene:", e$message)
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Update cluster/dataset filter choices
  # Update cluster filter choices (SANS datasets)
  observe({
    req(monocle_object())
    cds <- monocle_object()
    
    cluster_choices <- NULL
    
    if ("ClusterIdents" %in% colnames(colData(cds))) {
      cluster_choices <- sort(unique(as.character(colData(cds)$ClusterIdents)))
    } else if ("seurat_clusters" %in% colnames(colData(cds))) {
      cluster_choices <- sort(unique(as.character(colData(cds)$seurat_clusters)))
    }
    
    if (!is.null(cluster_choices)) {
      updatePickerInput(
        session,
        "clusters_filter_trajectory",
        choices = cluster_choices,
        selected = cluster_choices  # All selected by default
      )
    }
  })
  
  # Plot multiple genes along pseudotime
  # Plot multiple genes along pseudotime - VERSION AVEC COULEURS CLUSTERS
  observeEvent(input$plot_multi_gene_pseudotime_trajectory, {
    req(monocle_object(), input$multi_gene_pseudotime_trajectory)
    
    if (length(input$multi_gene_pseudotime_trajectory) == 0) {
      showNotification("Please select at least one gene", type = "warning")
      return()
    }
    
    if (length(input$multi_gene_pseudotime_trajectory) > 20) {
      showNotification("Maximum 20 genes allowed. Please reduce selection.", 
                       type = "warning", duration = 5)
      return()
    }
    
    showModal(modalDialog(
      title = "Generating Gene Expression Plot",
      "Creating visualization of gene expression along pseudotime...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      cds <- monocle_object()
      
      # Find cluster column name
      cluster_col <- NULL
      if ("ClusterIdents" %in% colnames(colData(cds))) {
        cluster_col <- "ClusterIdents"
      } else if ("seurat_clusters" %in% colnames(colData(cds))) {
        cluster_col <- "seurat_clusters"
      }
      
      # Apply cluster filter if selected
      if (!is.null(input$clusters_filter_trajectory) && length(input$clusters_filter_trajectory) > 0) {
        message(paste("Filtering cells by:", paste(input$clusters_filter_trajectory, collapse = ", ")))
        
        if (!is.null(cluster_col)) {
          cells_to_keep <- colData(cds)[[cluster_col]] %in% input$clusters_filter_trajectory
          
          # Get pseudotime BEFORE subsetting
          pt_values <- pseudotime(cds)
          
          # Subset the CDS
          cds_subset <- cds[, cells_to_keep]
          
          # Add pseudotime to the subset
          if (!is.null(pt_values)) {
            colData(cds_subset)$Pseudotime <- pt_values[cells_to_keep]
            
            if (!is.null(cds@principal_graph_aux) && 
                "UMAP" %in% names(cds@principal_graph_aux) &&
                "pseudotime" %in% names(cds@principal_graph_aux[["UMAP"]])) {
              cds_subset@principal_graph_aux[["UMAP"]]$pseudotime <- pt_values[cells_to_keep]
            }
            
            message(paste("Pseudotime successfully transferred to filtered cells"))
          } else {
            removeModal()
            showNotification("Pseudotime not calculated. Please calculate pseudotime first.", type = "error")
            return()
          }
          
          cds <- cds_subset
          message(paste("Filtered to", sum(cells_to_keep), "cells"))
        }
      }
      
      # Verify pseudotime exists
      has_pseudotime <- FALSE
      
      if ("Pseudotime" %in% colnames(colData(cds))) {
        has_pseudotime <- TRUE
        message("Pseudotime found in colData")
      }
      
      pt_check <- pseudotime(cds)
      if (!is.null(pt_check) && !all(is.na(pt_check))) {
        has_pseudotime <- TRUE
        message("Pseudotime found via pseudotime() function")
        if (!"Pseudotime" %in% colnames(colData(cds))) {
          colData(cds)$Pseudotime <- pt_check
          message("Added pseudotime to colData")
        }
      }
      
      if (!has_pseudotime) {
        removeModal()
        showNotification("Pseudotime not found. Please calculate pseudotime first.", type = "error")
        return()
      }
      
      # Verify genes exist
      available_genes <- rownames(cds)
      valid_genes <- input$multi_gene_pseudotime_trajectory[input$multi_gene_pseudotime_trajectory %in% available_genes]
      
      if (length(valid_genes) == 0) {
        removeModal()
        showNotification("Selected genes not found in dataset", type = "error")
        return()
      }
      
      if (length(valid_genes) < length(input$multi_gene_pseudotime_trajectory)) {
        missing_genes <- setdiff(input$multi_gene_pseudotime_trajectory, valid_genes)
        showNotification(
          paste("Warning: Some genes not found:", paste(missing_genes, collapse = ", ")), 
          type = "warning",
          duration = 5
        )
      }
      
      message(paste("Plotting", length(valid_genes), "genes along pseudotime"))
      
      # Create the plot - COLORER PAR CLUSTER
      # Create the plot - AVEC TOUTES LES AMÉLIORATIONS
      plot <- plot_genes_in_pseudotime(
        cds[valid_genes, ],
        color_cells_by = cluster_col,
        min_expr = input$multi_gene_min_expr_trajectory,
        cell_size = input$multi_gene_cell_size_trajectory,
        nrow = NULL,
        ncol = min(3, length(valid_genes)),
        panel_order = valid_genes,
        trend_formula = "~ splines::ns(Pseudotime, df=3)",
        label_by_short_name = TRUE
      )
      
      # Enhance plot - VERSION AMÉLIORÉE
      plot <- plot +
        theme_minimal(base_size = 14) +
        theme(
          # Supprimer le titre principal
          plot.title = element_blank(),
          
          # Légende à droite
          legend.position = "right",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          
          # Titres de facettes sans encadré - JUSTE LE NOM
          strip.text = element_text(size = 13, face = "bold"),
          strip.background = element_blank(),  # Pas d'encadré
          
          # Spacing
          panel.spacing = unit(1, "lines"),
          
          # Grille
          panel.grid.major = element_line(color = "grey90", size = 0.3),
          panel.grid.minor = element_blank(),
          
          # Bordures des panels
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          
          # Axes
          axis.text = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 12, face = "bold"),
          axis.line = element_line(color = "black", size = 0.5),
          
          # Fond blanc
          panel.background = element_rect(fill = "white", color = NA)
        ) +
        labs(
          x = "Pseudotime",
          y = "Expression"
        ) +
        # Format l'axe Y pour enlever les .0
        scale_y_continuous(
          labels = function(x) {
            ifelse(x == floor(x), as.character(as.integer(x)), as.character(x))
          }
        ) 

      if (length(valid_genes) == 1) {
        plot <- plot + theme(aspect.ratio = 0.6)
      }
      
      gene_trajectory_plot(plot)
      
      # Dynamic height - AJUSTÉ POUR TÉLÉCHARGEMENT
      plot_height <- 300 + (ceiling(length(valid_genes) / 3) * 200)
      
      output$multiGenePseudotimeTrajectoryPlot <- renderPlot({
        gene_trajectory_plot()
      }, height = plot_height)
      
      removeModal()
      
      showNotification(
        paste("Successfully plotted", length(valid_genes), "genes along pseudotime"), 
        type = "message",
        duration = 3
      )
      
    }, error = function(e) {
      removeModal()
      message("Error in gene trajectory plot:", e$message)
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
 
  # Download handlers using createDownloadHandler
  
  # Single gene on trajectory
  output$download_single_gene_umap_trajectory <- createDownloadHandler(
    reactive_data = single_gene_trajectory_plot,
    object_name_reactive = reactive({
      pt_id <- selected_pseudotime_id()
      gene_name <- gsub("[^A-Za-z0-9_]", "_", input$single_gene_umap_trajectory)
      
      if (!is.null(pt_id)) {
        history <- pseudotime_history()
        pt_entry <- history[[pt_id]]
        root_name <- gsub(" ", "_", pt_entry$root_cluster)
        paste0("SingleGene_", gene_name, "_root_", root_name)
      } else {
        paste0("SingleGene_", gene_name)
      }
    }),
    data_name = "trajectory_umap",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$download_format_trajectory }),
      width = 10,
      height = 8,
      dpi = reactive({ input$dpi_trajectory })
    )
  )
  
  # Multi-gene along pseudotime
  output$download_multi_gene_pseudotime_trajectory <- createDownloadHandler(
    reactive_data = gene_trajectory_plot,
    object_name_reactive = reactive({
      pt_id <- selected_pseudotime_id()
      
      n_genes <- length(input$multi_gene_pseudotime_trajectory)
      gene_part <- if(n_genes <= 3) {
        paste(input$multi_gene_pseudotime_trajectory, collapse = "_")
      } else {
        paste0(n_genes, "_genes")
      }
      
      if (!is.null(pt_id)) {
        history <- pseudotime_history()
        pt_entry <- history[[pt_id]]
        root_name <- gsub(" ", "_", pt_entry$root_cluster)
        paste0("MultiGene_", gene_part, "_root_", root_name)
      } else {
        paste0("MultiGene_", gene_part)
      }
    }),
    data_name = "pseudotime_expression",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$download_format_trajectory }),
      width = reactive({
        n_genes <- length(input$multi_gene_pseudotime_trajectory)
        # Largeur augmentée pour éviter compression
        min(15, 5 * min(3, n_genes))
      }),
      height = reactive({
        n_genes <- length(input$multi_gene_pseudotime_trajectory)
        # Hauteur augmentée
        8 + (ceiling(n_genes / 3) * 3)
      }),
      dpi = reactive({ input$dpi_trajectory })
    )
  )
  
  # SERVER CODE - 3D Trajectory Navigator with Slider
  # Add to trajectory_server_newarch.R
  
  ############################## 3D Trajectory Navigator ##############################
  
  # Reactive values for 3D trajectory
  trajectory_3d_plot <- reactiveVal()
  trajectory_top_genes <- reactiveVal()
  
  # Observer for slider changes - updates plot
  observeEvent(c(input$pseudotime_slider_3d, 
                 input$trajectory_window_size,
                 input$show_trajectory_line_3d,
                 input$dark_mode_trajectory_3d), {
                   
                   req(monocle_object())
                   
                   # Check if UMAP 3D exists
                   if(!"umap3d" %in% names(seurat_monocle()@reductions)) {
                     showNotification(
                       "3D UMAP not found. Please run 'Run UMAP 3D' first in Single Dataset module.",
                       type = "warning",
                       duration = 5
                     )
                     return()
                   }
                   
                   tryCatch({
                     # Get slider position
                     slider_pos <- input$pseudotime_slider_3d %||% 0.5
                     window_size <- input$trajectory_window_size %||% 0.05
                     show_line <- !is.null(input$show_trajectory_line_3d) && input$show_trajectory_line_3d
                     dark_mode <- !is.null(input$dark_mode_trajectory_3d) && input$dark_mode_trajectory_3d
                     
                     # Create 3D plot
                     result <- create3DTrajectoryNavigator(
                       monocle_object = monocle_object(),
                       seurat_object = seurat_monocle(),
                       reduction = "umap3d",
                       pseudotime_position = slider_pos,
                       window_size = window_size,
                       show_trajectory_line = show_line,
                       dark_mode = dark_mode
                     )
                     
                     # Store results
                     trajectory_3d_plot(result$plot)
                     trajectory_top_genes(result$top_genes)
                     
                     # Update info text
                     output$trajectory_3d_info <- renderText({
                       paste0(
                         "Cells highlighted: ", result$n_cells_highlighted, "\n",
                         "Pseudotime range: ", round(result$pseudotime_range[1], 3), 
                         " - ", round(result$pseudotime_range[2], 3)
                       )
                     })
                     
                   }, error = function(e) {
                     showNotification(paste("Error creating 3D trajectory:", e$message), type = "error")
                     print(e)
                   })
                 }, ignoreNULL = FALSE, ignoreInit = FALSE)
  
  # Render 3D plot
  output$trajectory_3d_plot <- renderPlotly({
    req(trajectory_3d_plot())
    trajectory_3d_plot()
  })
  
  # Render top genes
  output$trajectory_top_genes <- renderText({
    req(trajectory_top_genes())
    genes <- trajectory_top_genes()
    
    if(is.null(genes) || length(genes) == 0) {
      return("Not enough cells in window to calculate top genes (need at least 10)")
    }
    
    paste0("Top 15 genes at this position:\n\n", paste(genes, collapse = ", "))
  })
  
  # Download handler for 3D trajectory
  output$download_trajectory_3d <- downloadHandler(
    filename = function() {
      paste0("Trajectory3D_", Sys.Date(), ".html")
    },
    content = function(file) {
      req(trajectory_3d_plot())
      
      # Save as interactive HTML
      htmlwidgets::saveWidget(
        trajectory_3d_plot(),
        file = file,
        selfcontained = TRUE
      )
      
      showNotification("3D trajectory saved as interactive HTML", type = "message")
    }
  )
  
  # Animation control - auto-advance slider
  trajectory_animation_running <- reactiveVal(FALSE)
  
  # Play/Pause animation control
  observeEvent(input$play_trajectory_animation, {
    if(!trajectory_animation_running()) {
      # Start animation
      trajectory_animation_running(TRUE)
      updateActionButton(session, "play_trajectory_animation", label = tagList(icon("pause"), " Pause"))
    } else {
      # Stop animation
      trajectory_animation_running(FALSE)
      updateActionButton(session, "play_trajectory_animation", label = tagList(icon("play"), " Play"))
    }
  })
  
  
  # Animation loop (separate observer)
  observe({
    # Only run if animation is active
    if(trajectory_animation_running()) {
      
      isolate({
        current_pos <- input$pseudotime_slider_3d %||% 0
        
        # Increment position
        new_pos <- current_pos + 0.01
        
        # Loop back to start if at end
        if(new_pos > 1) {
          new_pos <- 0
        }
        
        # Update slider
        updateSliderInput(session, "pseudotime_slider_3d", value = new_pos)
      })
      
      # Wait before next frame (100ms = 10 fps)
      invalidateLater(100, session)
    }
  })
  
  # Reset slider to beginning
  observeEvent(input$reset_trajectory_slider, {
    updateSliderInput(session, "pseudotime_slider_3d", value = 0)
    trajectory_animation_running(FALSE)
    updateActionButton(session, "play_trajectory_animation", label = tagList(icon("play"), " Play"))
  })
  # Status for UMAP 3D
  output$umap_3d_status_trajectory <- renderText({
    req(seurat_monocle())
    
    if("umap3d" %in% names(seurat_monocle()@reductions)) {
      "✓ 3D UMAP ready"
    } else {
      "⚠ 3D UMAP not generated"
    }
  })
  
  # Generate UMAP 3D for trajectory
  observeEvent(input$run_umap_3d_trajectory, {
    req(seurat_monocle())
    
    showModal(modalDialog(
      title = "Generating 3D UMAP",
      "Please wait, this may take a few minutes...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      seurat_obj <- seurat_monocle()
      
      # Check if PCA exists
      if(!"pca" %in% names(seurat_obj@reductions)) {
        removeModal()
        showNotification("PCA not found. Please run PCA first.", type = "error")
        return()
      }
      
      # Run UMAP 3D using your existing function
      seurat_obj <- runUMAP3D_reproducible(
        object = seurat_obj,
        dims = 1:30,
        seed = 42,
        n.neighbors = 30L,
        min.dist = 0.3,
        metric = "cosine",
        verbose = TRUE
      )
      
      # Update object
      seurat_monocle(seurat_obj)
      
      removeModal()
      showNotification("✓ 3D UMAP generated successfully!", type = "message", duration = 5)
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error generating 3D UMAP:", e$message), type = "error")
      print(traceback())
    })
  })

}
