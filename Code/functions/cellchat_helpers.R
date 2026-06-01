############################## CellChat Helper Functions ##############################

# Update all CellChat object selectors in one call
update_cellchat_selectors <- function(session, object_names) {
  selector_ids <- c(
    "subset_choose_cellchat",
    "cellchat_obj_chord",
    "cellchat_obj_global",
    "cellchat_obj_specific",
    "select_objects_cellchat",
    "comp_obj1",
    "comp_obj2",
    "cellchat_obj_split_chord"  # added
  )
  for (id in selector_ids) {
    updateSelectInput(session, id, choices = object_names)
  }
}
# Update source and target cell type selectors
update_celltype_selectors <- function(session, cellchat_obj, 
                                      sources_id = "sources_use_cellchat",
                                      targets_id = "targets_use_cellchat") {
  if (is.null(cellchat_obj) || is.null(cellchat_obj@idents)) {
    return(NULL)
  }
  
  cell_types <- levels(cellchat_obj@idents)
  
  # Preserve existing selections if they're still valid
  current_sources <- session$input[[sources_id]]
  current_targets <- session$input[[targets_id]]
  
  updateSelectizeInput(
    session, sources_id,
    choices = cell_types,
    selected = intersect(current_sources, cell_types),
    options = list(plugins = list('remove_button'))
  )
  
  updateSelectizeInput(
    session, targets_id,
    choices = cell_types,
    selected = intersect(current_targets, cell_types),
    options = list(plugins = list('remove_button'))
  )
}

# Validate CellChat object before analysis or plotting
validate_cellchat_object <- function(cellchat_obj, min_groups = 2, check_db = TRUE) {
  errors <- c()
  n_groups <- 0L
  if (is.null(cellchat_obj)) {
    return(list(valid = FALSE, errors = "CellChat object is NULL", n_groups = n_groups))
  }
  if (is.null(cellchat_obj@idents)) {
    errors <- c(errors, "Object has no cell identities (@idents is NULL)")
  } else {
    n_groups <- length(levels(cellchat_obj@idents))
    if (n_groups < min_groups) {
      errors <- c(errors, sprintf("Need at least %d groups, found %d", min_groups, n_groups))
    }
  }
  if (check_db && is.null(cellchat_obj@DB)) {
    errors <- c(errors, "Database not set. Please load ligand-receptor database first")
  }
  if (is.null(cellchat_obj@data) && is.null(cellchat_obj@data.signaling)) {
    errors <- c(errors, "No expression data found in object")
  }
  return(list(valid = length(errors) == 0, errors = errors, n_groups = n_groups))

}

# Clean CellChat factors to avoid issues with unused levels
clean_cellchat_factors <- function(cellchat_obj) {
  # Clean identities
  if (!is.null(cellchat_obj@idents)) {
    cellchat_obj@idents <- droplevels(cellchat_obj@idents)
  }
  
  # Clean metadata factors
  if (!is.null(cellchat_obj@meta)) {
    for (col in names(cellchat_obj@meta)) {
      if (is.factor(cellchat_obj@meta[[col]])) {
        cellchat_obj@meta[[col]] <- droplevels(cellchat_obj@meta[[col]])
      }
    }
  }
  
  return(cellchat_obj)
}

# Prepare Seurat object for CellChat (remove NAs, clean factors)
prepare_seurat_for_cellchat <- function(seurat_obj, group_by_column) {
  tryCatch({
    message("=== Preparing Seurat object for CellChat ===")
    message("Original cells: ", ncol(seurat_obj))
    message("Grouping by: ", group_by_column)
    if (!group_by_column %in% colnames(seurat_obj@meta.data)) {
      stop(paste("Column", group_by_column, "not found in metadata"))
    }
    # Ensure RNA assay is active and layers are joined before any CellChat operation
    # (rules #1 and #2: JoinLayers mandatory before CellChat; LayerData for multi-layer objects)
    DefaultAssay(seurat_obj) <- "RNA"
    seurat_obj <- JoinLayers(seurat_obj)
    group_values <- seurat_obj@meta.data[[group_by_column]]
    group_values <- as.character(group_values)
    group_values[is.na(group_values)] <- "Unknown"
    group_values[group_values == ""]  <- "Unknown"
    group_values[group_values == "NA"] <- "Unknown"
    unique_values <- unique(group_values)
    seurat_obj@meta.data[[group_by_column]] <- factor(group_values, levels = unique_values)
    if (!"samples" %in% colnames(seurat_obj@meta.data)) {
      seurat_obj@meta.data$samples <- "sample1"
    }
    for (col in colnames(seurat_obj@meta.data)) {
      if (is.factor(seurat_obj@meta.data[[col]])) {
        seurat_obj@meta.data[[col]] <- droplevels(seurat_obj@meta.data[[col]])
      }
    }
    message("Final cells: ", ncol(seurat_obj), " (NO CELLS REMOVED)")
    message("Groups found: ", paste(unique_values, collapse = ", "))
    return(seurat_obj)
  }, error = function(e) {
    stop(paste("Error preparing Seurat object:", e$message))
  })
}

# Create single CellChat object with optional filtering and combined grouping
create_single_cellchat_object <- function(seurat_obj, 
                                          group_by,
                                          selected_clusters = NULL,
                                          condition_column = NULL,
                                          selected_conditions = NULL,
                                          combine_with_column = NULL,
                                          combined_column_name = "combined_group",
                                          log_function = message) {
  
  # Apply filters first if specified
  if (!is.null(selected_clusters) || !is.null(selected_conditions)) {
    log_function("Applying filters to Seurat object...")
    
    seurat_obj <- filter_seurat_for_cellchat(
      seurat_obj = seurat_obj,
      cluster_column = group_by,
      selected_clusters = selected_clusters,
      condition_column = condition_column,
      selected_conditions = selected_conditions,
      log_function = log_function
    )
  }
  
  # If combination is requested, create combined column
  if (!is.null(combine_with_column)) {
    log_function(sprintf("Creating combined grouping: %s + %s", 
                         group_by, combine_with_column))
    
    seurat_obj <- create_combined_grouping_column(
      seurat_obj = seurat_obj,
      primary_column = group_by,
      secondary_column = combine_with_column,
      combined_column_name = combined_column_name,
      separator = "_"
    )
    
    # Update group_by to use the combined column
    group_by <- combined_column_name
  }
  
  # Prepare object (remove NAs)
  seurat_obj <- prepare_seurat_for_cellchat(seurat_obj, group_by)
  
  # Check groups
  n_groups <- length(unique(seurat_obj@meta.data[[group_by]]))
  if (n_groups < 2) {
    return(list(success = FALSE, message = paste("Only", n_groups, "group(s)")))
  }
  
  log_function(sprintf("Creating CellChat object: %d cells, %d groups", 
                       ncol(seurat_obj), n_groups))
  
  # Extract data
  data_input <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  meta_data <- seurat_obj@meta.data
  
  # Create object
  cellchat_obj <- clean_cellchat_factors(cellchat_obj)
  
  # Propagate Seurat cluster colors into CellChat @options for downstream plots
  cellchat_obj <- store_colors_in_cellchat(
    cellchat_obj, seurat_obj,
    group_by    = group_by,
    primary_col = if (!is.null(combine_with_column)) group_by else NULL
  )
  
  return(list(success = TRUE, object = cellchat_obj))
}

# Run complete CellChat analysis pipeline for one object
run_cellchat_pipeline <- function(cellchat_obj, object_name = "object", 
                                  log_function = message, 
                                  min_cells = 10) {
  result <- list(
    success = FALSE,
    object = cellchat_obj,
    message = "",
    n_interactions = 0
  )
  
  tryCatch({
    # Validate before starting
    validation <- validate_cellchat_object(cellchat_obj, min_groups = 2, check_db = TRUE)
    if (!validation$valid) {
      result$message <- paste("Validation failed:", paste(validation$errors, collapse = "; "))
      return(result)
    }
    
    log_function(sprintf("[%s] Starting analysis pipeline...", object_name))
    
    # Step 1: Subset data
    log_function(sprintf("[%s] Subsetting data...", object_name))
    cellchat_obj <- subsetData(cellchat_obj)
    
    # Step 2: Identify over-expressed genes
    log_function(sprintf("[%s] Identifying over-expressed genes...", object_name))
    cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
    
    # Step 3: Identify over-expressed interactions
    log_function(sprintf("[%s] Identifying over-expressed interactions...", object_name))
    cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
    
    # Step 4: Compute communication probability
    log_function(sprintf("[%s] Computing communication probability...", object_name))
    cellchat_obj <- computeCommunProb(cellchat_obj)
    
    # Step 5: Filter communications
    log_function(sprintf("[%s] Filtering communications (min.cells = %d)...", object_name, min_cells))
    cellchat_obj <- filterCommunication(cellchat_obj, min.cells = min_cells)
    
    # Step 6: Compute pathway communication
    log_function(sprintf("[%s] Computing pathway communication...", object_name))
    cellchat_obj <- computeCommunProbPathway(cellchat_obj)
    
    # Step 7: Aggregate network
    log_function(sprintf("[%s] Aggregating network...", object_name))
    cellchat_obj <- aggregateNet(cellchat_obj)
    
    # Count interactions
    n_interactions <- 0
    if (!is.null(cellchat_obj@net) && !is.null(cellchat_obj@net$count)) {
      n_interactions <- sum(cellchat_obj@net$count > 0)
    }
    
    log_function(sprintf("[%s] ✓ Analysis completed! Found %d interactions", 
                         object_name, n_interactions))
    
    result$success <- TRUE
    result$object <- cellchat_obj
    result$message <- "Analysis completed successfully"
    result$n_interactions <- n_interactions
    
  }, error = function(e) {
    error_msg <- sprintf("[%s] ERROR: %s", object_name, e$message)
    log_function(error_msg)
    result$message <- error_msg
  })
  
  return(result)
}

# Get summary information about a CellChat object
get_cellchat_summary <- function(cellchat_obj) {
  if (is.null(cellchat_obj)) {
    return("No object loaded")
  }
  
  # Basic info
  n_cells <- if (!is.null(cellchat_obj@data.signaling)) ncol(cellchat_obj@data.signaling) else 0
  n_genes <- if (!is.null(cellchat_obj@data.signaling)) nrow(cellchat_obj@data.signaling) else 0
  n_groups <- if (!is.null(cellchat_obj@idents)) length(levels(cellchat_obj@idents)) else 0
  group_names <- if (!is.null(cellchat_obj@idents)) paste(levels(cellchat_obj@idents), collapse = ", ") else "None"
  
  summary_text <- sprintf(
    "Cells: %d\nGenes: %d\nGroups: %d\nGroup names: %s",
    n_cells, n_genes, n_groups, group_names
  )
  
  # Analysis status
  if (!is.null(cellchat_obj@net) && length(cellchat_obj@net) > 0) {
    n_interactions <- if (!is.null(cellchat_obj@net$count)) sum(cellchat_obj@net$count > 0) else 0
    summary_text <- paste0(
      summary_text,
      sprintf("\n\n✓ Analysis completed\nInteractions found: %d", n_interactions)
    )
  } else {
    summary_text <- paste0(summary_text, "\n\n⚠ Analysis not performed yet")
  }
  
  return(summary_text)
}

# Check if CellChat object has been analyzed
is_cellchat_analyzed <- function(cellchat_obj) {
  if (is.null(cellchat_obj)) return(FALSE)
  return(!is.null(cellchat_obj@net) && length(cellchat_obj@net) > 0)
}

############################## Plot Generation Helpers ##############################

# Validate CellChat object before plotting
validate_for_plotting <- function(cellchat_obj, check_analysis = TRUE) {
  # Basic validation
  validation <- validate_cellchat_object(cellchat_obj, min_groups = 2, check_db = FALSE)
  if (!validation$valid) {
    return(list(valid = FALSE, message = paste(validation$errors, collapse = "; ")))
  }
  
  # Check if analysis was performed
  if (check_analysis) {
    if (is.null(cellchat_obj@net) || length(cellchat_obj@net) == 0) {
      return(list(valid = FALSE, message = "Please run CellChat analysis first"))
    }
  }
  
  return(list(valid = TRUE, message = "OK"))
}

# Generate bubble plot
# Generate bubble plot
generate_bubble_plot <- function(cellchat_obj, sources, targets, threshold = 0.05,
                                 flip_axes = FALSE, exclude_intra = FALSE,
                                 signaling = NULL) {
  validation <- validate_for_plotting(cellchat_obj, check_analysis = TRUE)
  if (!validation$valid) stop(validation$message)
  if (length(sources) == 0 || length(targets) == 0) {
    stop("Please select at least one source and one target cell type")
  }
  communications <- subsetCommunication(
    cellchat_obj,
    sources.use = sources,
    targets.use = targets,
    thresh = threshold
  )
  if (is.null(communications) || nrow(communications) == 0) {
    stop("No interactions found with current parameters")
  }
  # signaling= is the correct parameter in CellChat 1.6.1; pairLR.use is broken
  plot <- netVisual_bubble(
    object = cellchat_obj,
    sources.use = sources,
    targets.use = targets,
    signaling = signaling,
    thresh = threshold,
    remove.isolate = exclude_intra,
    angle.x = 45,
    font.size = 10
  )
  if (flip_axes) plot <- plot + coord_flip()
  plot <- plot + theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10)
  )
  return(list(
    plot = plot,
    communications = communications,
    n_interactions = nrow(communications)
  ))
}

# Generate chord plot
generate_chord_plot <- function(cellchat_obj, senders, receivers, threshold = 0.05,
                                custom_colors = NULL, label_size = 0.8) {
  # Validate
  validation <- validate_for_plotting(cellchat_obj, check_analysis = TRUE)
  if (!validation$valid) {
    stop(validation$message)
  }
  
  if (length(senders) == 0 || length(receivers) == 0) {
    stop("Please select sender and receiver groups")
  }
  
  # Create plot function
  plot_function <- function() {
    netVisual_chord_gene(
      cellchat_obj,
      sources.use = senders,
      targets.use = receivers,
      thresh = threshold,
      color.use = custom_colors,
      show.legend = TRUE,
      small.gap = 3,
      big.gap = 10,
      lab.cex = label_size,
      directional = 1,
      link.target.prop = FALSE
    )
  }
  
  # Count interactions
  communications <- tryCatch({
    subsetCommunication(cellchat_obj,
                        sources.use = senders,
                        targets.use = receivers,
                        thresh = threshold)
  }, error = function(e) NULL)
  
  n_interactions <- if (!is.null(communications)) nrow(communications) else 0
  
  return(list(
    plot_function = plot_function,
    n_interactions = n_interactions
  ))
}

# Generate global circle plot
generate_global_circle_plot <- function(cellchat_obj, plot_type = "weight",
                                        vertex_size = 8, edge_width = 8,
                                        label_size = 1, margin = 0.1,
                                        title = NULL) {
  # Validate
  validation <- validate_for_plotting(cellchat_obj, check_analysis = TRUE)
  if (!validation$valid) {
    stop(validation$message)
  }
  
  # Create plot function
  plot_function <- function() {
    groupSize <- as.numeric(table(cellchat_obj@idents))
    
    netVisual_circle(
      if (plot_type == "count") cellchat_obj@net$count else cellchat_obj@net$weight,
      vertex.weight = groupSize,
      weight.scale = TRUE,
      vertex.size = vertex_size,
      edge.width.max = edge_width,
      vertex.label.cex = label_size,
      label.edge = FALSE,
      title.name = if (!is.null(title)) title else paste(plot_type, "network"),
      margin = margin
    )
  }
  
  return(list(plot_function = plot_function))
}

# Generate cell-type specific plots
generate_specific_circle_plots <- function(cellchat_obj, cell_types, 
                                           signal_direction = "outgoing",
                                           label_size = 1, margin = 0.1,
                                           n_cols = 2) {
  # Validate
  validation <- validate_for_plotting(cellchat_obj, check_analysis = TRUE)
  if (!validation$valid) {
    stop(validation$message)
  }
  
  if (length(cell_types) == 0) {
    stop("Please select at least one cell type")
  }
  
  # Create plot function
  plot_function <- function() {
    groupSize <- as.numeric(table(cellchat_obj@idents))
    n_plots <- length(cell_types)
    n_cols_actual <- min(n_cols, n_plots)
    n_rows <- ceiling(n_plots / n_cols_actual)
    
    par(mfrow = c(n_rows, n_cols_actual),
        mar = c(1, 1, 3, 1),
        oma = c(0, 0, 2, 0))
    
    # Get matrix
    mat <- if (signal_direction == "outgoing") {
      cellchat_obj@net$weight
    } else {
      t(cellchat_obj@net$weight)
    }
    
    # Plot each cell type
    for (cell_type in cell_types) {
      mat_subset <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), 
                           dimnames = dimnames(mat))
      
      if (cell_type %in% rownames(mat)) {
        mat_subset[cell_type, ] <- mat[cell_type, ]
      }
      
      netVisual_circle(
        mat_subset,
        vertex.weight = groupSize,
        weight.scale = TRUE,
        edge.width.max = max(mat),
        vertex.label.cex = label_size,
        title.name = paste(cell_type, 
                           ifelse(signal_direction == "outgoing",
                                  "(Outgoing)", "(Incoming)")),
        margin = margin
      )
    }
  }
  
  return(list(
    plot_function = plot_function,
    n_plots = length(cell_types)
  ))
}

# Get interaction type information
get_interaction_type_info <- function(cellchat_obj, communications) {
  if (is.null(cellchat_obj@DB$interaction) ||
      !("interaction_type" %in% colnames(cellchat_obj@DB$interaction))) {
    return(NULL)
  }
  
  displayed_interactions <- unique(communications$interaction_name)
  
  interaction_info <- data.frame(
    interaction_name = cellchat_obj@DB$interaction$interaction_name,
    interaction_type = cellchat_obj@DB$interaction$interaction_type,
    stringsAsFactors = FALSE
  )
  
  interaction_info <- interaction_info[
    interaction_info$interaction_name %in% displayed_interactions, 
  ]
  
  interaction_info$interaction_type[is.na(interaction_info$interaction_type)] <- "Unknown"
  interaction_info <- interaction_info[order(interaction_info$interaction_type,
                                             interaction_info$interaction_name), ]
  
  return(interaction_info)
}




#  Filter Seurat object by selected metadata values and clusters
filter_seurat_for_cellchat <- function(seurat_obj, 
                                       cluster_column,
                                       selected_clusters = NULL,
                                       condition_column = NULL,
                                       selected_conditions = NULL,
                                       log_function = message) {
  tryCatch({
    original_ncells <- ncol(seurat_obj)
    log_function(sprintf("Starting with %d cells", original_ncells))
    
    # Filter by clusters if specified
    if (!is.null(selected_clusters) && length(selected_clusters) > 0) {
      if (!cluster_column %in% colnames(seurat_obj@meta.data)) {
        stop(paste("Cluster column", cluster_column, "not found"))
      }
      
      # Convert both to character to ensure matching
      cluster_values <- as.character(seurat_obj@meta.data[[cluster_column]])
      selected_clusters_char <- as.character(selected_clusters)
      
      log_function(sprintf("Available clusters in data: %s", 
                           paste(unique(cluster_values), collapse = ", ")))
      log_function(sprintf("Selected clusters: %s", 
                           paste(selected_clusters_char, collapse = ", ")))
      
      cells_to_keep <- cluster_values %in% selected_clusters_char
      n_kept <- sum(cells_to_keep)
      
      if (n_kept == 0) {
        stop(sprintf("No cells match the selected clusters. Check cluster names."))
      }
      
      seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
      
      log_function(sprintf("After cluster filter: %d cells (%d clusters selected)", 
                           ncol(seurat_obj), length(selected_clusters_char)))
    }
    
    # Filter by conditions if specified
    if (!is.null(condition_column) && !is.null(selected_conditions) && 
        length(selected_conditions) > 0) {
      if (!condition_column %in% colnames(seurat_obj@meta.data)) {
        stop(paste("Condition column", condition_column, "not found"))
      }
      
      # Convert both to character to ensure matching
      condition_values <- as.character(seurat_obj@meta.data[[condition_column]])
      selected_conditions_char <- as.character(selected_conditions)
      
      log_function(sprintf("Available conditions in data: %s", 
                           paste(unique(condition_values), collapse = ", ")))
      log_function(sprintf("Selected conditions: %s", 
                           paste(selected_conditions_char, collapse = ", ")))
      
      cells_to_keep <- condition_values %in% selected_conditions_char
      n_kept <- sum(cells_to_keep)
      
      if (n_kept == 0) {
        stop(sprintf("No cells match the selected conditions. Check condition names."))
      }
      
      seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
      
      log_function(sprintf("After condition filter: %d cells (%d conditions selected)", 
                           ncol(seurat_obj), length(selected_conditions_char)))
    }
    
    # Check if we still have enough cells
    if (ncol(seurat_obj) < 10) {
      stop(sprintf("Too few cells remaining after filtering (%d cells)", ncol(seurat_obj)))
    }
    
    log_function(sprintf("Final: %d cells (%.1f%% of original)", 
                         ncol(seurat_obj), 
                         100 * ncol(seurat_obj) / original_ncells))
    
    return(seurat_obj)
    
  }, error = function(e) {
    stop(paste("Error filtering Seurat object:", e$message))
  })
}

# Create combined grouping column (e.g., cluster + condition)
create_combined_grouping_column <- function(seurat_obj, 
                                            primary_column, 
                                            secondary_column,
                                            combined_column_name = "combined_group",
                                            separator = "_") {
  tryCatch({
    # Check if columns exist
    if (!primary_column %in% colnames(seurat_obj@meta.data)) {
      stop(paste("Primary column", primary_column, "not found in metadata"))
    }
    if (!secondary_column %in% colnames(seurat_obj@meta.data)) {
      stop(paste("Secondary column", secondary_column, "not found in metadata"))
    }
    
    # Get values from both columns
    primary_values <- as.character(seurat_obj@meta.data[[primary_column]])
    secondary_values <- as.character(seurat_obj@meta.data[[secondary_column]])
    
    # Handle NAs
    primary_values[is.na(primary_values)] <- "Unknown"
    secondary_values[is.na(secondary_values)] <- "Unknown"
    
    # Combine with separator
    combined_values <- paste(primary_values, secondary_values, sep = separator)
    
    # Add to metadata as factor
    seurat_obj@meta.data[[combined_column_name]] <- factor(combined_values)
    
    message(sprintf("Created combined column '%s' with %d unique groups", 
                    combined_column_name, 
                    length(unique(combined_values))))
    message(sprintf("Groups: %s", 
                    paste(head(unique(combined_values), 10), collapse = ", ")))
    
    return(seurat_obj)
    
  }, error = function(e) {
    stop(paste("Error creating combined column:", e$message))
  })
}

# Create CellChat object with support for multiple grouping columns
create_cellchat_object_multicol <- function(seurat_obj, 
                                            grouping_columns,
                                            selected_groups = NULL,
                                            condition_column = NULL,
                                            selected_conditions = NULL,
                                            log_function = message) {
  
  # Validate grouping columns
  if (length(grouping_columns) == 0) {
    return(list(success = FALSE, message = "No grouping columns specified"))
  }
  
  for (col in grouping_columns) {
    if (!col %in% colnames(seurat_obj@meta.data)) {
      return(list(success = FALSE, message = paste("Column not found:", col)))
    }
  }
  
  # Filter by additional condition if specified
  if (!is.null(condition_column) && !is.null(selected_conditions) && length(selected_conditions) > 0) {
    log_function(sprintf("Filtering by %s: %s", condition_column, paste(selected_conditions, collapse = ", ")))
    cells_to_keep <- which(seurat_obj@meta.data[[condition_column]] %in% selected_conditions)
    seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
    log_function(sprintf("After filtering: %d cells remaining", ncol(seurat_obj)))
  }
  
  # Create combined grouping column if multiple columns specified
  if (length(grouping_columns) > 1) {
    log_function(sprintf("Creating combined grouping from %d columns", length(grouping_columns)))
    
    combined_column_name <- "CellChat_Combined_Groups"
    seurat_obj@meta.data[[combined_column_name]] <- seurat_obj@meta.data[[grouping_columns[1]]]
    
    for (col in grouping_columns[-1]) {
      seurat_obj@meta.data[[combined_column_name]] <- paste(
        seurat_obj@meta.data[[combined_column_name]],
        seurat_obj@meta.data[[col]],
        sep = "_"
      )
    }
    
    group_by <- combined_column_name
    
  } else {
    group_by <- grouping_columns[1]
  }
  
  # Filter specific groups if specified
  if (!is.null(selected_groups) && length(selected_groups) > 0) {
    log_function(sprintf("Keeping only %d selected groups", length(selected_groups)))
    cells_to_keep <- which(seurat_obj@meta.data[[group_by]] %in% selected_groups)
    
    if (length(cells_to_keep) == 0) {
      return(list(success = FALSE, message = "No cells match the selected groups"))
    }
    
    seurat_obj <- subset(seurat_obj, cells = cells_to_keep)
    log_function(sprintf("After group filtering: %d cells remaining", ncol(seurat_obj)))
  }
  
  # Prepare object (remove NAs from grouping column)
  seurat_obj <- prepare_seurat_for_cellchat(seurat_obj, group_by)
  
  # Check minimum requirements
  if (ncol(seurat_obj) < 50) {
    return(list(success = FALSE, message = paste("Too few cells:", ncol(seurat_obj), "cells remaining")))
  }
  
  n_groups <- length(unique(seurat_obj@meta.data[[group_by]]))
  if (n_groups < 2) {
    return(list(success = FALSE, message = paste("Need at least 2 groups, found:", n_groups)))
  }
  
  log_function(sprintf("Creating CellChat object: %d cells, %d groups", ncol(seurat_obj), n_groups))
  
  # Extract data for CellChat
  data_input <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  meta_data <- seurat_obj@meta.data
  
  # Create CellChat object
  cellchat_obj <- createCellChat(
    object = data_input,
    meta = meta_data,
    group.by = group_by
  )
  
  # Clean factors
  cellchat_obj <- clean_cellchat_factors(cellchat_obj)
  
  log_function("CellChat object created successfully")
  
  return(list(success = TRUE, object = cellchat_obj))
}

# Get available LR pairs for selected sources/targets
# Returns a named vector: names = "LIGAND_RECEPTOR (PATHWAY)", values = pathway_name
get_split_chord_lr_choices <- function(cellchat_obj, sources, targets, threshold = 0.05) {
  comms <- tryCatch(
    subsetCommunication(cellchat_obj,
                        sources.use = sources,
                        targets.use = targets,
                        thresh      = threshold),
    error = function(e) NULL
  )
  
  if (is.null(comms) || nrow(comms) == 0) {
    stop("No interactions found for selected sources/targets")
  }
  
  # Rank by mean probability
  comms <- comms[order(comms$prob, decreasing = TRUE), ]
  
  # Build named vector: label shown to user → pathway used internally
  labels   <- paste0(comms$interaction_name, " (", comms$pathway_name, ")")
  pathways <- comms$pathway_name
  
  # Deduplicate keeping highest prob per label
  seen   <- !duplicated(labels)
  labels <- labels[seen]
  pathways <- pathways[seen]
  
  result <- setNames(pathways, labels)
  return(result)
}

# Generate split chord plots: one plot per group in split_by,
# using selected LR pairs. Sectors are shared across all panels so absent
# interactions appear as gray ghosts, enabling direct visual comparison.
# Colors are resolved from: (1) custom_colors arg, (2) @options$cluster_colors
# stored at CellChat creation time, (3) deterministic brewer/Polychrome fallback.
generate_split_chord_plots <- function(
    cellchat_obj,
    sources          = NULL,
    targets          = NULL,
    split_by         = "source",
    pairs            = NULL,        # list of list(label, sources, targets) — overrides split_by logic
    selected_lr_keys = NULL,
    threshold        = 0.05,
    label_size       = 0.6,
    label_distance   = 0.3,
    custom_colors    = NULL,
    show_title       = FALSE,
    gap_inner        = 1,
    gap_outer        = 8,
    plot_size        = 1200
) {
  if (is.null(selected_lr_keys) || length(selected_lr_keys) == 0) {
    stop("Please select at least one LR pair")
  }
  
  # Determine split_groups and how to fetch comms per group
  use_pairs <- !is.null(pairs) && length(pairs) >= 2
  if (use_pairs) {
    split_groups <- sapply(pairs, `[[`, "label")
  } else {
    split_groups <- if (split_by == "source") sources else targets
  }
  
  message("=== generate_split_chord_plots ===")
  message("--- mode          : ", if (use_pairs) "pairs" else paste0("split_by=", split_by))
  message("--- split_groups  : ", paste(split_groups, collapse = ", "))
  message("--- selected LR   : ", length(selected_lr_keys), " pairs")
  message("--- threshold     : ", threshold)
  
  # --- Build per-group communication tables ---
  comms_list <- setNames(lapply(seq_along(split_groups), function(i) {
    grp   <- split_groups[i]
    s_use <- if (use_pairs) pairs[[i]]$sources else if (split_by == "source") grp else sources
    t_use <- if (use_pairs) pairs[[i]]$targets else if (split_by == "target") grp else targets
    
    comms <- tryCatch(
      subsetCommunication(cellchat_obj,
                          sources.use = s_use,
                          targets.use = t_use,
                          thresh      = threshold),
      error = function(e) {
        message(sprintf("[%s] subsetCommunication error: %s", grp, conditionMessage(e)[1]))
        NULL
      }
    )
    
    n_raw <- if (!is.null(comms)) nrow(comms) else 0
    message(sprintf("[%s] raw interactions: %d", grp, n_raw))
    
    if (!is.null(comms) && nrow(comms) > 0) {
      comms$ligand   <- as.character(comms$ligand)
      comms$receptor <- as.character(comms$receptor)
      comms$source   <- as.character(comms$source)
      comms$target   <- as.character(comms$target)
      comms$lr_pair  <- paste(comms$ligand, comms$receptor, sep = " -> ")
      comms <- comms[comms$lr_pair %in% selected_lr_keys, , drop = FALSE]
      message(sprintf("[%s] after LR filter: %d", grp, nrow(comms)))
      if (nrow(comms) > 0) {
        top_show <- head(comms[order(comms$prob, decreasing = TRUE), ], min(5, nrow(comms)))
        message(sprintf("[%s] top interactions:\n%s", grp,
                        paste(sprintf("  %s:%s -> %s:%s (%.5f)",
                                      top_show$source, top_show$ligand,
                                      top_show$receptor, top_show$target,
                                      top_show$prob), collapse = "\n")))
      }
      if (nrow(comms) == 0) return(NULL)
    }
    comms
  }), split_groups)
  
  # --- Build union of all interactions across groups ---
  union_comms <- do.call(rbind, Filter(Negate(is.null), comms_list))
  if (is.null(union_comms) || nrow(union_comms) == 0) stop("No interactions found")
  union_comms$ligand   <- as.character(union_comms$ligand)
  union_comms$receptor <- as.character(union_comms$receptor)
  union_comms$source   <- as.character(union_comms$source)
  union_comms$target   <- as.character(union_comms$target)
  union_comms$lr_key   <- paste0(union_comms$source, "|",
                                 union_comms$ligand, "->",
                                 union_comms$receptor, "|",
                                 union_comms$target)
  union_comms <- union_comms[order(union_comms$prob, decreasing = TRUE), ]
  union_comms <- union_comms[!duplicated(union_comms$lr_key), ]
  message("--- unique interactions in union: ", nrow(union_comms))
  
  # --- Resolve cell-type colors ---
  all_cts <- unique(c(union_comms$source, union_comms$target))
  n_cts   <- length(all_cts)
  
  resolve_colors <- function(color_source, cts) {
    out     <- setNames(rep("gray50", length(cts)), cts)
    matched <- cts[cts %in% names(color_source)]
    out[matched] <- color_source[matched]
    out
  }
  
  if (!is.null(custom_colors) && !is.null(names(custom_colors)) &&
      sum(all_cts %in% names(custom_colors)) > 0) {
    ct_colors <- resolve_colors(custom_colors, all_cts)
    message(sprintf("--- colors: custom_colors (%d/%d matched)", sum(all_cts %in% names(custom_colors)), n_cts))
  } else {
    stored <- tryCatch(cellchat_obj@options$cluster_colors, error = function(e) NULL)
    if (!is.null(stored) && !is.null(names(stored)) && sum(all_cts %in% names(stored)) > 0) {
      ct_colors <- resolve_colors(stored, all_cts)
      message(sprintf("--- colors: @options$cluster_colors (%d/%d matched)", sum(all_cts %in% names(stored)), n_cts))
    } else {
      pal <- if (n_cts <= 8) {
        RColorBrewer::brewer.pal(max(3, n_cts), "Set1")[seq_len(n_cts)]
      } else {
        set.seed(42)
        as.character(Polychrome::createPalette(n_cts, c("#FF0000", "#00FF00", "#0000FF")))
      }
      ct_colors <- setNames(pal, all_cts)
      message(sprintf("--- colors: fallback palette (%d colors)", n_cts))
    }
  }
  
  # --- Build sectors ---
  lig_by_src  <- split(union_comms$ligand,   union_comms$source)
  rec_by_tgt  <- split(union_comms$receptor, union_comms$target)
  lig_sectors <- unlist(lapply(names(lig_by_src), function(ct) {
    paste0(ct, "_L_", unique(lig_by_src[[ct]]))
  }))
  rec_sectors <- unlist(lapply(names(rec_by_tgt), function(ct) {
    paste0(ct, "_R_", unique(rec_by_tgt[[ct]]))
  }))
  all_sectors <- c(lig_sectors, rec_sectors)
  
  sector_gene <- sub("^.+_[LR]_", "", all_sectors)
  sector_ct   <- sub("_[LR]_[^_]+$", "", all_sectors)
  names(sector_gene) <- all_sectors
  names(sector_ct)   <- all_sectors
  
  sector_col <- sapply(all_sectors, function(s) {
    ct <- sector_ct[s]
    if (!is.na(ct) && ct %in% names(ct_colors)) ct_colors[[ct]] else "gray80"
  })
  
  sector_sizes <- setNames(rep(0, length(all_sectors)), all_sectors)
  for (j in seq_len(nrow(union_comms))) {
    ls <- paste0(union_comms$source[j], "_L_", union_comms$ligand[j])
    rs <- paste0(union_comms$target[j], "_R_", union_comms$receptor[j])
    if (ls %in% all_sectors) sector_sizes[ls] <- sector_sizes[ls] + union_comms$prob[j]
    if (rs %in% all_sectors) sector_sizes[rs] <- sector_sizes[rs] + union_comms$prob[j]
  }
  sector_sizes <- pmax(sector_sizes, 1e-6)
  
  src_cts  <- unique(sub("_L_.*$", "", lig_sectors))
  tgt_cts  <- unique(sub("_R_.*$", "", rec_sectors))
  src_lens <- sapply(src_cts, function(ct) sum(startsWith(lig_sectors, paste0(ct, "_L_"))))
  tgt_lens <- sapply(tgt_cts, function(ct) sum(startsWith(rec_sectors, paste0(ct, "_R_"))))
  gap_lig  <- as.numeric(unlist(lapply(src_lens, function(l) c(rep(gap_inner, max(l-1, 0)), gap_outer))))
  gap_rec  <- as.numeric(unlist(lapply(tgt_lens, function(l) c(rep(gap_inner, max(l-1, 0)), gap_outer))))
  gaps     <- as.numeric(c(gap_lig, gap_rec))
  
  n  <- length(split_groups)
  sg <- split_groups
  ld <- label_distance
  st <- show_title
  
  render_one_plot <- function(grp, grp_comms, grp_keys) {
    canvas_margin <- ld + 0.35
    canvas_lim    <- c(-1 - canvas_margin, 1 + canvas_margin)
    
    circlize::circos.clear()
    circlize::circos.par(
      gap.degree              = gaps,
      cell.padding            = c(0, 0, 0, 0),
      start.degree            = 90,
      clock.wise              = TRUE,
      track.margin            = c(0.005, 0.005),
      points.overflow.warning = FALSE,
      canvas.xlim             = canvas_lim,
      canvas.ylim             = canvas_lim
    )
    circlize::circos.initialize(
      factors = factor(all_sectors, levels = all_sectors),
      xlim    = cbind(rep(0, length(all_sectors)), sector_sizes[all_sectors])
    )
    circlize::circos.trackPlotRegion(
      ylim = c(0, 1), bg.col = sector_col, bg.border = "white",
      track.height = 0.06, panel.fun = function(x, y) {}
    )
    for (s in all_sectors) {
      xlim <- circlize::get.cell.meta.data("xlim", sector.index = s, track.index = 1)
      circlize::circos.text(
        x = mean(xlim), y = 1 + ld, labels = sector_gene[s],
        sector.index = s, track.index = 1,
        facing = "clockwise", niceFacing = TRUE,
        adj = c(0, 0.5), cex = label_size, col = "black", font = 3
      )
    }
    
    pos <- setNames(rep(0, length(all_sectors)), all_sectors)
    
    absent_comms <- union_comms[!union_comms$lr_key %in% grp_keys, ]
    for (j in seq_len(nrow(absent_comms))) {
      ls <- paste0(absent_comms$source[j], "_L_", absent_comms$ligand[j])
      rs <- paste0(absent_comms$target[j], "_R_", absent_comms$receptor[j])
      if (!ls %in% all_sectors || !rs %in% all_sectors) next
      prob   <- absent_comms$prob[j]
      lstart <- pos[ls]; lend <- min(lstart + prob, sector_sizes[ls])
      rstart <- pos[rs]; rend <- min(rstart + prob, sector_sizes[rs])
      if ((lend - lstart) < 1e-10 || (rend - rstart) < 1e-10) next
      circlize::circos.link(
        ls, c(lstart, lend), rs, c(rstart, rend),
        col    = adjustcolor("white", alpha.f = 0),
        border = adjustcolor("gray75", alpha.f = 0.45), lwd = 0.4
      )
      pos[ls] <- lend; pos[rs] <- rend
    }
    
    if (!is.null(grp_comms) && nrow(grp_comms) > 0) {
      grp_sorted <- grp_comms[order(grp_comms$prob, decreasing = FALSE), ]
      for (j in seq_len(nrow(grp_sorted))) {
        ls  <- paste0(grp_sorted$source[j], "_L_", grp_sorted$ligand[j])
        rs  <- paste0(grp_sorted$target[j], "_R_", grp_sorted$receptor[j])
        src <- grp_sorted$source[j]
        if (!ls %in% all_sectors || !rs %in% all_sectors) next
        arc_col    <- adjustcolor(
          if (src %in% names(ct_colors)) ct_colors[[src]] else "gray40",
          alpha.f = 0.75)
        border_col <- if (src %in% names(ct_colors)) ct_colors[[src]] else "gray40"
        lstart <- pos[ls]; lend <- min(lstart + grp_sorted$prob[j], sector_sizes[ls])
        rstart <- pos[rs]; rend <- min(rstart + grp_sorted$prob[j], sector_sizes[rs])
        if ((lend - lstart) < 1e-10 || (rend - rstart) < 1e-10) next
        circlize::circos.link(
          ls, c(lstart, lend), rs, c(rstart, rend),
          col         = arc_col, border = border_col,
          directional = 1, arr.type = "big.arrow", arr.length = 0.04, lwd = 0.5
        )
        pos[ls] <- lend; pos[rs] <- rend
      }
    }
    
    if (st) title(main = grp, cex.main = 1.8, font.main = 2, line = -1)
    
    legend(x = canvas_lim[2] - 0.05, y = canvas_lim[1] + 0.05,
           xjust = 1, yjust = 0,
           legend = c(names(ct_colors), "Absent"),
           fill   = c(unname(ct_colors), adjustcolor("gray75", alpha.f = 0.5)),
           title  = "Cell Type", cex = 0.75, bty = "n", border = NA)
    
    circlize::circos.clear()
  }
  
  # --- Render one PNG per group then stitch side-by-side ---
  tmp_files <- vapply(seq_len(n), function(i) {
    grp       <- sg[i]
    grp_comms <- comms_list[[grp]]
    grp_keys  <- if (!is.null(grp_comms) && nrow(grp_comms) > 0) {
      paste0(grp_comms$source, "|", grp_comms$ligand, "->",
             grp_comms$receptor, "|", grp_comms$target)
    } else character(0)
    
    n_present <- length(grp_keys)
    n_absent  <- sum(!union_comms$lr_key %in% grp_keys)
    message(sprintf("--- [%d/%d] %s | present: %d | absent (ghost): %d",
                    i, n, grp, n_present, n_absent))
    if (n_present > 0)
      message(sprintf("  present:\n%s", paste(" ", grp_keys, collapse = "\n")))
    if (n_absent > 0)
      message(sprintf("  ghost (first 10):\n%s",
                      paste(" ", head(union_comms$lr_key[!union_comms$lr_key %in% grp_keys], 10),
                            collapse = "\n")))
    
    tmp <- tempfile(fileext = ".png")
    png(tmp, width = plot_size, height = plot_size, res = 150)
    tryCatch(
      render_one_plot(grp, grp_comms, grp_keys),
      error = function(e) {
        message("  ERROR rendering ", grp, ": ", conditionMessage(e)[1])
        circlize::circos.clear()
        par(mar = c(2, 2, 4, 2)); plot.new()
        title(main = grp, cex.main = 1.4, font.main = 2)
        text(0.5, 0.5, paste0("Error:\n", conditionMessage(e)[1]), col = "red", cex = 0.9)
      }
    )
    dev.off()
    message(sprintf("  PNG: %d bytes", file.size(tmp)))
    tmp
  }, character(1))
  
  imgs        <- lapply(tmp_files, magick::image_read)
  combined    <- magick::image_append(do.call(c, imgs), stack = FALSE)
  output_path <- tempfile(fileext = ".png")
  magick::image_write(combined, path = output_path, format = "png")
  file.remove(tmp_files)
  message("=== combined PNG: ", file.size(output_path), " bytes")
  
  plot_function <- function() {
    img <- png::readPNG(output_path)
    par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
    plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
    rasterImage(img, 0, 0, 1, 1)
  }
  
  export_function <- function(n_plots) {
    par(mfrow = c(1, n_plots), mar = c(0, 0, 2, 0))
    for (i in seq_len(n)) {
      grp       <- sg[i]
      grp_comms <- comms_list[[grp]]
      grp_keys  <- if (!is.null(grp_comms) && nrow(grp_comms) > 0) {
        paste0(grp_comms$source, "|", grp_comms$ligand, "->",
               grp_comms$receptor, "|", grp_comms$target)
      } else character(0)
      tryCatch(
        render_one_plot(grp, grp_comms, grp_keys),
        error = function(e) {
          plot.new()
          text(0.5, 0.5, paste0(grp, "\nError: ", e$message), col = "red")
        }
      )
    }
  }
  
  return(list(
    plot_function    = plot_function,
    export_function  = export_function,
    output_path      = output_path,
    n_plots          = n,
    selected_lr_keys = selected_lr_keys,
    png_size         = plot_size
  ))
}

# Store cluster colors from a Seurat object into a CellChat object's @options slot.
# Handles combined group names (e.g., "B_cells_WT") by falling back to the primary
# group color when an exact match is not found.
# Args:
#   cellchat_obj : CellChat object to annotate
#   seurat_obj   : source Seurat object containing @misc$cluster_colors
#   group_by     : column name used as grouping (may be a combined column)
#   primary_col  : original grouping column before combination (for fallback matching)
# Returns: CellChat object with colors stored in @options$cluster_colors
store_colors_in_cellchat <- function(cellchat_obj, seurat_obj,
                                     group_by, primary_col = NULL) {
  seurat_colors <- tryCatch(seurat_obj@misc$cluster_colors, error = function(e) NULL)
  if (is.null(seurat_colors) || length(seurat_colors) == 0) return(cellchat_obj)
  
  group_levels <- levels(cellchat_obj@idents)
  if (is.null(group_levels) || length(group_levels) == 0) return(cellchat_obj)
  
  resolved <- setNames(rep(NA_character_, length(group_levels)), group_levels)
  
  # Direct match first
  direct_match <- group_levels %in% names(seurat_colors)
  resolved[direct_match] <- seurat_colors[group_levels[direct_match]]
  
  # Fallback for combined groups: match on the primary column prefix
  if (!is.null(primary_col) && any(is.na(resolved))) {
    unmatched <- group_levels[is.na(resolved)]
    for (lvl in unmatched) {
      # Combined groups are "primary_secondary" — try matching on primary prefix
      for (ct in names(seurat_colors)) {
        if (startsWith(lvl, ct)) {
          resolved[[lvl]] <- seurat_colors[[ct]]
          break
        }
      }
    }
  }
  
  # Fill remaining unmatched with a deterministic fallback palette
  still_na <- is.na(resolved)
  if (any(still_na)) {
    n_missing <- sum(still_na)
    fallback_pal <- if (n_missing <= 8) {
      RColorBrewer::brewer.pal(max(3, n_missing), "Set2")[seq_len(n_missing)]
    } else {
      set.seed(42)
      as.character(Polychrome::createPalette(n_missing, c("#FF0000", "#00FF00", "#0000FF")))
    }
    resolved[still_na] <- fallback_pal
  }
  
  cellchat_obj@options$cluster_colors <- resolved
  message(sprintf("[store_colors_in_cellchat] Stored %d colors (%d direct match, %d fallback)",
                  length(resolved), sum(direct_match), sum(still_na)))
  return(cellchat_obj)
}
