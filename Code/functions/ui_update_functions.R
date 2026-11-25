# ui_update_functions.R

updateAssayChoices <- function(session, seurat_object, assay_input_id = "viz_assay") {
  # Update assay choices based on available assays in Seurat object
  # Args:
  #   session: Shiny session object
  #   seurat_object: Seurat object
  #   assay_input_id: ID of the selectInput to update
  
  if (is.null(seurat_object)) {
    return(NULL)
  }
  
  tryCatch({
    # Get available assays
    available_assays <- names(seurat_object@assays)
    
    # Debug message
    message(paste("Available assays:", paste(available_assays, collapse = ", ")))
    
    # Update selectInput
    updateSelectInput(session, assay_input_id,
                      choices = available_assays,
                      selected = if("RNA" %in% available_assays) "RNA" else available_assays[1])
    
  }, error = function(e) {
    message(paste("Error updating assay choices:", e$message))
  })
}

updateGeneChoices <- function(session, seurat_object, selected_assay, gene_input_ids = NULL) {
  # Update gene picker choices with the features of the selected assay
  # Args:
  #   session: Shiny session object
  #   seurat_object: Seurat object
  #   selected_assay: Currently selected assay
  #   gene_input_ids: Vector of input IDs to update (if NULL, updates default ones)
  if (is.null(seurat_object) || is.null(selected_assay)) return(NULL)
  if (is.null(gene_input_ids)) gene_input_ids <- c("gene_select", "gene_select_heatmap", "gene_select_genes_analysis")
  tryCatch({
    # Features of the assay itself: LayerData(layer = "counts") reads only the first layer of a merged Seurat v5 object
    # (counts.1) and fails on assays without counts (integrated)
    gene_list <- sort(rownames(seurat_object[[selected_assay]]))
    for (input_id in gene_input_ids) updatePickerInput(session, input_id, choices = gene_list)
    message(paste("Updated", length(gene_input_ids), "gene input(s) with", length(gene_list), "genes"))
  }, error = function(e) {
    message(paste("Error updating gene choices:", e$message))
  })
}

updateClusterChoices <- function(session, seurat_object, cluster_input_ids = NULL,
                                 cluster_col = NULL) {
  # Update cluster selector choices from a metadata column (never from Idents directly).
  # Args:
  #   session          : Shiny session
  #   seurat_object    : Seurat object
  #   cluster_input_ids: Named list with keys "select" and/or "checkbox" → vectors of input IDs
  #   cluster_col      : Metadata column to use; if NULL, auto-detected via get_cluster_column()
  if (is.null(seurat_object)) return(NULL)
  if (is.null(cluster_input_ids)) {
    cluster_input_ids <- list(
      select   = c("select_cluster", "ident_1"),
      checkbox = c("ident_2")
    )
  }
  tryCatch({
    col <- if (!is.null(cluster_col) && cluster_col %in% colnames(seurat_object@meta.data)) {
      cluster_col
    } else if ("ClusterIdents" %in% colnames(seurat_object@meta.data)) {
      "ClusterIdents"
    } else if ("seurat_clusters" %in% colnames(seurat_object@meta.data)) {
      "seurat_clusters"
    } else {
      NULL
    }
    cluster_choices <- if (!is.null(col)) {
      unique(as.character(seurat_object@meta.data[[col]]))
    } else {
      as.character(unique(Idents(seurat_object)))
    }
    message(paste("updateClusterChoices — column:", ifelse(is.null(col), "Idents()", col),
                  "— clusters:", paste(cluster_choices, collapse = ", ")))
    if ("select" %in% names(cluster_input_ids)) {
      for (input_id in cluster_input_ids$select) {
        updateSelectInput(session, input_id, choices = cluster_choices)
      }
    }
    if ("checkbox" %in% names(cluster_input_ids)) {
      for (input_id in cluster_input_ids$checkbox) {
        updateCheckboxGroupInput(session, input_id, choices = cluster_choices, selected = "")
      }
    }
  }, error = function(e) {
    message(paste("Error in updateClusterChoices:", conditionMessage(e)[1]))
  })
}

updateClusterSelectors <- function(session, seurat_object, cluster_col = NULL,
                                   rename_id   = "select_cluster",
                                   color_id    = "select_color",
                                   order_ids   = NULL) {
  # Update rename + color + order selectors from the active cluster column.
  # Shared between single_dataset_server and multiple_datasets_server.
  # Args:
  #   session       : Shiny session
  #   seurat_object : Seurat object
  #   cluster_col   : Metadata column to read clusters from; auto-detected if NULL
  #   rename_id     : selectInput ID for cluster renaming
  #   color_id      : selectInput ID for color assignment
  #   order_ids     : optional character vector of additional selectInput IDs to update
  if (is.null(seurat_object)) return(NULL)
  tryCatch({
    col <- if (!is.null(cluster_col) && cluster_col %in% colnames(seurat_object@meta.data)) {
      cluster_col
    } else if ("ClusterIdents" %in% colnames(seurat_object@meta.data)) {
      "ClusterIdents"
    } else if ("seurat_clusters" %in% colnames(seurat_object@meta.data)) {
      "seurat_clusters"
    } else {
      NULL
    }
    clusters <- if (!is.null(col)) {
      unique(as.character(seurat_object@meta.data[[col]]))
    } else {
      as.character(unique(Idents(seurat_object)))
    }
    message(paste("updateClusterSelectors — column:", ifelse(is.null(col), "Idents()", col),
                  "—", length(clusters), "clusters"))
    updateSelectInput(session, rename_id, choices = clusters)
    updateSelectInput(session, color_id,  choices = clusters)
    if (!is.null(order_ids)) {
      for (id in order_ids) {
        updateSelectInput(session, id, choices = clusters)
      }
    }
  }, error = function(e) {
    message(paste("Error in updateClusterSelectors:", conditionMessage(e)[1]))
  })
}

createGeneUpdateFunction <- function(input) {
  # Returns a function that can access current input values
  function(session, new_gene, gene_text_input_ids = NULL) {
    if (is.null(new_gene) || new_gene == "") {
      return(NULL)
    }
    
    # Default gene text input IDs if not specified
    if (is.null(gene_text_input_ids)) {
      gene_text_input_ids <- list(
        gene_list_vln = input$gene_list_vln,
        gene_list_feature = input$gene_list_feature,
        gene_list_dotplot = input$gene_list_dotplot,
        gene_list_ridge_plot = input$gene_list_ridge_plot
      )
    }
    
    # Helper function to update genes
    update_genes <- function(existing_genes, new_gene) {
      if (is.null(existing_genes) || existing_genes == "") {
        return(new_gene)
      } else {
        genes <- strsplit(existing_genes, ",")[[1]]
        genes <- unique(c(trimws(genes), new_gene))
        return(paste(genes, collapse = ", "))
      }
    }
    
    # Update all specified text inputs
    for (input_id in names(gene_text_input_ids)) {
      current_value <- gene_text_input_ids[[input_id]]
      updateTextInput(session, input_id, value = update_genes(current_value, new_gene))
    }
  }
}


# Dans ui_update_functions_newarch.R, ajouter :
updateClusterTextInputs <- function(session, seurat_object, text_input_ids) {
  if (is.null(seurat_object)) return(NULL)
  tryCatch({
    cluster_list <- getClusters(seurat_object)
    cluster_text <- paste(cluster_list, collapse = ",")
    for (input_id in text_input_ids) {
      updateTextInput(session, input_id, value = cluster_text)
    }
    message(paste("Updated", length(text_input_ids), "cluster text input(s)"))
  }, error = function(e) {
    message(paste("Error updating cluster text inputs:", e$message))
  })
}

# Populate all dynamic UMAP/cluster selectors for a given module.
# Call from an observe() after object loads or changes.
# Args:
#   session    : Shiny session
#   seurat_obj : loaded Seurat object
#   module     : "single" or "merge"
update_umap_selectors <- function(session, seurat_obj, module = c("single", "merge"),
                                  selected_col = NULL) {
  module <- match.arg(module)
  if (is.null(seurat_obj)) return(invisible(NULL))
  meta_cols     <- colnames(seurat_obj@meta.data)
  priority_cols <- c("ClusterIdents", "seurat_clusters", "cell_type", "celltype",
                     "leiden", "louvain", "annotated_clusters", "clusters")
  default_col   <- if (!is.null(selected_col) && selected_col %in% meta_cols) {
    selected_col
  } else {
    col <- intersect(priority_cols, meta_cols)
    if (length(col) > 0) col[1] else meta_cols[1]
  }
  reductions    <- names(seurat_obj@reductions)
  umap_choices  <- reductions[grepl("umap", reductions, ignore.case = TRUE)]
  default_reduc <- if ("umap" %in% reductions) "umap"
  else if ("wnn.umap" %in% reductions) "wnn.umap"
  else if (length(umap_choices) > 0) umap_choices[1]
  else character(0)
  if (module == "single") {
    updateSelectInput(session, "umap_col_single",          choices = meta_cols,    selected = default_col)
    updateSelectInput(session, "cluster_col_single",       choices = meta_cols,    selected = default_col)
    updateSelectInput(session, "umap_reduction_single",    choices = umap_choices, selected = default_reduc)
    updateSelectInput(session, "cluster_reduction_single", choices = umap_choices, selected = default_reduc)
  } else {
    updateSelectInput(session, "umap_finale_col_merge",         choices = meta_cols,    selected = default_col)
    updateSelectInput(session, "filtered_umap_col_merge",       choices = meta_cols,    selected = default_col)
    updateSelectInput(session, "umap_reduction_merge",          choices = umap_choices, selected = default_reduc)
    updateSelectInput(session, "filtered_umap_reduction_merge", choices = umap_choices, selected = default_reduc)
    # Reduction selector of the Cluster Cells plot (Clustering tab), same choices as the other UMAP selectors
    updateSelectInput(session, "cluster_cells_reduction_merge", choices = umap_choices, selected = default_reduc)
  }
  message(sprintf("[update_umap_selectors] module=%s col=%s reduc=%s", module, default_col,
                  if (length(default_reduc) > 0) default_reduc else "none"))
  return(invisible(NULL))
}