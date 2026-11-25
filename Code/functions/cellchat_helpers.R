############################## CellChat Helper Functions ##############################

# Update every CellChat object selector present in the UI in one call
update_cellchat_selectors <- function(session, object_names) {
  selector_ids <- c("subset_choose_cellchat", "cellchat_obj_chord", "cellchat_obj_split_chord")
  for (id in selector_ids) updateSelectInput(session, id, choices = object_names)
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

# Prepare Seurat object for CellChat: RNA assay with joined layers, factors cleaned, and cells without a usable label
# removed. Cells with no label (NA, empty) or a placeholder label (Unknown, Unassigned, Debris, alone or as one part of a
# combined "cluster_condition" label) would otherwise form a group that sends and receives signals and enters the
# permutations. The number of removed cells is reported in the log and stored in @misc$cellchat_removed_cells.
prepare_seurat_for_cellchat <- function(seurat_obj, group_by_column, drop_labels = c("Unknown", "Unassigned", "Debris")) {
  tryCatch({
    message("=== Preparing Seurat object for CellChat ===")
    message("Original cells: ", ncol(seurat_obj))
    message("Grouping by: ", group_by_column)
    if (!group_by_column %in% colnames(seurat_obj@meta.data)) {
      stop(paste("Column", group_by_column, "not found in metadata"))
    }
    DefaultAssay(seurat_obj) <- "RNA"
    # JoinLayers only applies to Seurat v5 Assay5 objects.
    # Legacy Assay (v4) objects have no layers slot and no JoinLayers method — skip silently.
    rna_assay <- seurat_obj[["RNA"]]
    if (inherits(rna_assay, "Assay5")) {
      n_layers <- length(names(rna_assay@layers))
      if (n_layers > 1) {
        message("Joining ", n_layers, " RNA layers (Seurat v5 object)...")
        seurat_obj <- JoinLayers(seurat_obj)
      } else {
        message("RNA assay has ", n_layers, " layer(s) — no join needed.")
      }
    } else {
      message("RNA assay is Seurat v4 class (Assay) — skipping JoinLayers.")
    }
    group_values <- as.character(seurat_obj@meta.data[[group_by_column]])
    placeholder <- paste0("^(", paste(drop_labels, collapse = "|"), ")(_|$)|_(", paste(drop_labels, collapse = "|"), ")$")
    unlabeled <- is.na(group_values) | group_values %in% c("", "NA") | grepl(placeholder, group_values)
    if (all(unlabeled)) stop("No cell has a usable label in the grouping column")
    if (any(unlabeled)) {
      message(sprintf("Removing %d cells without a usable label (%s)", sum(unlabeled),
                      paste(utils::head(unique(group_values[unlabeled]), 10), collapse = ", ")))
      seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[!unlabeled])
      group_values <- group_values[!unlabeled]
    }
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
    seurat_obj@misc$cellchat_removed_cells <- sum(unlabeled)
    message("Final cells: ", ncol(seurat_obj), " (", sum(unlabeled), " cells without a usable label removed)")
    message("Groups found: ", paste(unique_values, collapse = ", "))
    return(seurat_obj)
  }, error = function(e) {
    stop(paste("Error preparing Seurat object:", conditionMessage(e)[1]))
  })
}



## Group: Run complete CellChat analysis pipeline for one object (permutations with seed 42)
run_cellchat_pipeline <- function(cellchat_obj, object_name = "object",
                                  log_function = message,
                                  min_cells = 10) {
  result <- list(success = FALSE, object = cellchat_obj, message = "", n_interactions = 0)
  validation <- validate_cellchat_object(cellchat_obj, min_groups = 2, check_db = TRUE)
  if (!validation$valid) {
    result$message <- paste("Validation failed:", paste(validation$errors, collapse = "; "))
    return(result)
  }
  # Name of the running step, reported in the error message when CellChat stops
  current_step <- "start"
  # The value of tryCatch is the result: an error handler cannot modify `result` in the calling function
  result <- tryCatch({
    log_function(sprintf("[%s] Starting analysis pipeline...", object_name))
    current_step <- "subsetData"
    log_function(sprintf("[%s] Subsetting data...", object_name))
    cellchat_obj <- subsetData(cellchat_obj)
    current_step <- "identifyOverExpressedGenes"
    log_function(sprintf("[%s] Identifying over-expressed genes...", object_name))
    cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
    current_step <- "identifyOverExpressedInteractions"
    log_function(sprintf("[%s] Identifying over-expressed interactions...", object_name))
    cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
    current_step <- "computeCommunProb"
    log_function(sprintf("[%s] Computing communication probability...", object_name))
    # Permutation test seeded with 42 (CellChat default seed.use is 1)
    cellchat_obj <- computeCommunProb(cellchat_obj, seed.use = 42)
    current_step <- "filterCommunication"
    log_function(sprintf("[%s] Filtering communications (min.cells = %d)...", object_name, min_cells))
    cellchat_obj <- filterCommunication(cellchat_obj, min.cells = min_cells)
    current_step <- "computeCommunProbPathway"
    log_function(sprintf("[%s] Computing pathway communication...", object_name))
    cellchat_obj <- computeCommunProbPathway(cellchat_obj)
    current_step <- "aggregateNet"
    log_function(sprintf("[%s] Aggregating network...", object_name))
    cellchat_obj <- aggregateNet(cellchat_obj)
    # net$count[i, j] = number of significant LR pairs from group i to group j: the total is the number of interactions
    n_interactions <- if (!is.null(cellchat_obj@net$count)) as.integer(sum(cellchat_obj@net$count)) else 0L
    n_group_pairs <- if (!is.null(cellchat_obj@net$count)) sum(cellchat_obj@net$count > 0) else 0L
    log_function(sprintf("[%s] ✓ Analysis completed! %d significant interactions between %d sender-receiver pairs",
                         object_name, n_interactions, n_group_pairs))
    list(success = TRUE, object = cellchat_obj, message = "Analysis completed successfully", n_interactions = n_interactions)
  }, error = function(e) {
    error_msg <- sprintf("[%s] ERROR at %s: %s", object_name, current_step, conditionMessage(e)[1])
    log_function(error_msg)
    list(success = FALSE, object = cellchat_obj, message = error_msg, n_interactions = 0)
  })
  return(result)
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

# Keep the size (p-value) and colour legends of a netVisual_bubble plot whole after rows of plot$data are removed:
# their breaks and labels are fixed vectors computed on all interactions, and ggplot drops breaks outside the remaining
# data range, which leaves more labels than breaks ("breaks and labels have different lengths"). The limits are set to the
# range of the original breaks, so every legend class stays shown, as in the unfiltered plot.
keep_bubble_legend_scales <- function(plot) {
  for (aesthetic in c("size", "colour")) {
    sc <- plot$scales$get_scales(aesthetic)
    if (!is.null(sc) && is.null(sc$limits) && is.numeric(sc$breaks) && length(sc$breaks) > 0) {
      sc <- sc$clone()
      sc$limits <- range(sc$breaks, na.rm = TRUE)
      plot <- suppressMessages(plot + sc)
    }
  }
  plot
}

# Grid lines of a netVisual_bubble plot (one between each LR pair and each source -> target column) are computed by
# CellChat on all interactions; after rows of plot$data are removed they would leave empty rows and columns.
# They are removed and redrawn on the remaining pairs and columns, after the points as in netVisual_bubble.
refresh_bubble_grid <- function(plot) {
  is_grid <- vapply(plot$layers, function(l) inherits(l$geom, "GeomHline") || inherits(l$geom, "GeomVline"), logical(1))
  if (!any(is_grid)) return(plot)
  first <- plot$layers[[which(is_grid)[1]]]$aes_params
  grid_col <- if (is.null(first$colour)) "grey90" else first$colour
  grid_lwd <- if (!is.null(first$linewidth)) first$linewidth else if (!is.null(first$size)) first$size else 0.1
  plot$layers <- plot$layers[!is_grid]
  n_x <- length(unique(plot$data$source.target)); n_y <- length(unique(plot$data$interaction_name_2))
  if (n_x > 1) plot <- plot + geom_vline(xintercept = seq(1.5, n_x - 0.5, 1), linewidth = grid_lwd, colour = grid_col)
  if (n_y > 1) plot <- plot + geom_hline(yintercept = seq(1.5, n_y - 0.5, 1), linewidth = grid_lwd, colour = grid_col)
  plot
}

## Group: Generate bubble plot; exclude_intra removes autocrine pairs (source == target) from the plot and the returned table
generate_bubble_plot <- function(cellchat_obj, sources, targets, threshold = 0.05,
                                 flip_axes = FALSE, exclude_intra = FALSE,
                                 signaling = NULL, genes = NULL) {
  validation <- validate_for_plotting(cellchat_obj, check_analysis = TRUE)
  if (!validation$valid) stop(validation$message)
  if (length(sources) == 0 || length(targets) == 0) stop("Please select at least one source and one target cell type")
  communications <- subsetCommunication(cellchat_obj, sources.use = sources, targets.use = targets, thresh = threshold)
  if (!is.null(communications) && isTRUE(exclude_intra)) {
    communications <- communications[as.character(communications$source) != as.character(communications$target), , drop = FALSE]
  }
  if (is.null(communications) || nrow(communications) == 0) {
    stop(if (isTRUE(exclude_intra)) "No interaction between different cell groups with current parameters" else "No interactions found with current parameters")
  }
  # Optional gene list: keeps the significant interactions whose ligand or receptor (or a complex subunit) is in the list.
  # Applied after inference, so probabilities and p-values are those computed on the whole database
  gene_match <- NULL
  if (length(genes) > 0) {
    gene_match <- match_interactions_to_genes(communications, cellchat_obj@DB$complex, genes)
    if (!any(gene_match$keep)) stop(sprintf("None of the %d genes of the list is in a significant interaction with current parameters", length(genes)))
    communications <- communications[gene_match$keep, , drop = FALSE]
  }
  # signaling= is used instead of pairLR.use (broken in CellChat 1.6.1, not re-tested on 2.1.2)
  plot <- netVisual_bubble(object = cellchat_obj, sources.use = sources, targets.use = targets, signaling = signaling,
                           thresh = threshold, angle.x = 45, font.size = 10)
  if (isTRUE(exclude_intra)) {
    # x axis of netVisual_bubble = "source -> target"; autocrine columns are dropped from the plot data
    if ("source.target" %in% colnames(plot$data)) {
      st <- as.character(plot$data$source.target)
      plot$data <- plot$data[sub(" -> .*$", "", st) != sub("^.* -> ", "", st), , drop = FALSE]
      if (is.factor(plot$data$source.target)) plot$data$source.target <- droplevels(plot$data$source.target)
    } else {
      warning("netVisual_bubble data has no 'source.target' column: intra-group interactions could not be removed from the plot")
    }
  }
  if (!is.null(gene_match)) {
    # Plot rows restricted to the interactions kept by the gene list (same approach as exclude_intra)
    if ("interaction_name" %in% colnames(plot$data)) {
      plot$data <- plot$data[as.character(plot$data$interaction_name) %in% as.character(communications$interaction_name), , drop = FALSE]
    } else if ("interaction_name_2" %in% colnames(plot$data) && "interaction_name_2" %in% colnames(communications)) {
      plot$data <- plot$data[trimws(as.character(plot$data$interaction_name_2)) %in% trimws(as.character(communications$interaction_name_2)), , drop = FALSE]
    } else {
      stop("netVisual_bubble data has no interaction name column: the gene list could not be applied to the plot")
    }
    if (nrow(plot$data) == 0) stop("No plotted interaction left after applying the gene list")
    if (is.factor(plot$data$interaction_name_2)) plot$data$interaction_name_2 <- droplevels(plot$data$interaction_name_2)
  }
  # Rows removed from plot$data: legends rebuilt on the original classes so breaks and labels stay paired
  if (isTRUE(exclude_intra) || !is.null(gene_match)) plot <- refresh_bubble_grid(keep_bubble_legend_scales(plot))
  if (flip_axes) plot <- plot + coord_flip()
  plot <- plot + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), axis.text.y = element_text(size = 10))
  return(list(plot = plot, communications = communications, n_interactions = nrow(communications),
              genes_absent = if (is.null(gene_match)) character(0) else gene_match$genes_absent))
}

# Generate gene-level chord plot (circlize-based).
# Drawn into a temporary PNG outside any Shiny render context: circlize's internal device
# handling corrupts renderPlot() (rule 8), so the server displays the file with renderImage().
# plot_function is also returned so download handlers can redraw on their own device.
# threshold is a p-value cutoff; top_n keeps the strongest interactions by communication probability (0 = all).
# genes (optional): keeps the significant interactions whose ligand or receptor (or a complex subunit) is in the list,
# before top_n, so the top interactions are chosen among those of the list; applied after inference.
generate_chord_plot <- function(cellchat_obj, senders, receivers, threshold = 0.05, top_n = 0,
                                custom_colors = NULL, label_size = 0.8, plot_size = 1200, genes = NULL) {
  validation <- validate_for_plotting(cellchat_obj, check_analysis = TRUE)
  if (!validation$valid) stop(validation$message)
  if (length(senders) == 0 || length(receivers) == 0) stop("Please select sender and receiver groups")
  communications <- tryCatch(
    subsetCommunication(cellchat_obj, sources.use = senders, targets.use = receivers, thresh = threshold),
    error = function(e) NULL
  )
  n_significant <- if (!is.null(communications)) nrow(communications) else 0
  if (n_significant == 0) stop("No interactions found with current parameters")
  chord_net <- as.data.frame(communications)
  # Optional gene list, applied to the significant interactions before the top-n selection
  gene_match <- NULL
  if (length(genes) > 0) {
    gene_match <- match_interactions_to_genes(chord_net, cellchat_obj@DB$complex, genes)
    if (!any(gene_match$keep)) stop(sprintf("None of the %d genes of the list is in a significant interaction with current parameters", length(genes)))
    chord_net <- chord_net[gene_match$keep, , drop = FALSE]
  }
  n_matching <- nrow(chord_net)
  # Permutation p-values are often 0 for every significant pair, so the strength (prob) is what reduces the plot
  if (is.numeric(top_n) && !is.na(top_n) && top_n > 0 && top_n < nrow(chord_net)) {
    chord_net <- chord_net[order(chord_net$prob, decreasing = TRUE), , drop = FALSE][seq_len(top_n), , drop = FALSE]
  }
  n_interactions <- nrow(chord_net)
  # Sectors are named by gene only: a gene used both as ligand (sender side) and as receptor (receiver side),
  # frequent with the bidirectional pairs of the merged database, would share one sector coloured by the sender,
  # so sender-to-receiver links would appear to loop back onto the sender. Such receptor sectors get a " (R)" suffix.
  shared_genes <- intersect(as.character(chord_net$ligand), as.character(chord_net$receptor))
  chord_net$receptor <- ifelse(chord_net$receptor %in% shared_genes, paste0(chord_net$receptor, " (R)"), as.character(chord_net$receptor))
  if (length(shared_genes) > 0) message("Chord plot: receptor sectors renamed to avoid ligand/receptor merging: ", paste(shared_genes, collapse = ", "))
  # Gaps scaled to the number of sectors (one per ligand x sender and per receptor x receiver) so that all gaps
  # take at most 35% of the circle; circlize fails when gaps exceed 360 degrees
  big_gap <- 10
  n_sectors <- nrow(unique(chord_net[, c("ligand", "source")])) + nrow(unique(chord_net[, c("receptor", "target")]))
  n_groups <- length(unique(c(as.character(chord_net$source), as.character(chord_net$target))))
  small_gap <- max(0.05, min(3, (360 * 0.35 - n_groups * big_gap) / n_sectors))
  message(sprintf("Chord plot: %d of %d significant interactions (%d matching the gene list), %d sectors, small gap %.2f degrees",
                  n_interactions, n_significant, n_matching, n_sectors, small_gap))
  # Draws on whatever device is currently open; circos state is reset before and after
  plot_function <- function() {
    circlize::circos.clear()
    on.exit(circlize::circos.clear(), add = TRUE)
    netVisual_chord_gene(
      cellchat_obj,
      net = chord_net,
      sources.use = senders,
      targets.use = receivers,
      thresh = threshold,
      color.use = custom_colors,
      show.legend = TRUE,
      small.gap = small_gap,
      big.gap = big_gap,
      lab.cex = label_size,
      directional = 1,
      link.target.prop = FALSE
    )
    invisible(NULL)
  }
  # Render once into a dedicated PNG device and always close that exact device, even on error
  output_path <- tempfile(pattern = "chord_", fileext = ".png")
  png(output_path, width = plot_size, height = plot_size, res = 150)
  png_dev <- dev.cur()
  draw_error <- tryCatch({
    par(mar = c(0, 0, 0, 0))
    plot_function()
    NULL
  }, error = function(e) conditionMessage(e)[1])
  if (png_dev %in% dev.list()) dev.off(png_dev)
  if (!is.null(draw_error)) {
    unlink(output_path)
    stop(paste("Chord plot rendering failed:", draw_error))
  }
  return(list(
    plot_function  = plot_function,
    output_path    = output_path,
    png_size       = plot_size,
    n_interactions = n_interactions,
    n_significant  = n_significant,
    n_matching     = n_matching,
    genes_absent   = if (is.null(gene_match)) character(0) else gene_match$genes_absent
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

## Group: Group label of each cell from one or several metadata columns joined with "_"; a cell with a missing value
## Group: (NA or empty) in any column gets NA, so it is removed instead of forming a "<cluster>_NA" group.
## Group: CellChat refuses the label "0" (createCellChat), so purely numeric labels (cluster numbers) get a "C" prefix;
## Group: preview, group filter and object creation all call this helper, so the labels stay identical everywhere.
build_cellchat_group_labels <- function(meta, columns) {
  missing_cols <- setdiff(columns, colnames(meta))
  if (length(missing_cols) > 0) stop("Column(s) not found in metadata: ", paste(missing_cols, collapse = ", "))
  parts <- lapply(columns, function(col) {
    v <- as.character(meta[[col]])
    v[!is.na(v) & !nzchar(trimws(v))] <- NA_character_
    v
  })
  labels <- do.call(paste, c(parts, sep = "_"))
  labels[Reduce(`|`, lapply(parts, is.na))] <- NA_character_
  numeric_label <- !is.na(labels) & grepl("^[0-9]+$", labels)
  labels[numeric_label] <- paste0("C", labels[numeric_label])
  # A prefixed label must not merge with an existing label of the same text (e.g. "C1" next to "1")
  clash <- intersect(unique(labels[numeric_label]), unique(labels[!numeric_label & !is.na(labels)]))
  if (length(clash) > 0) stop("Group labels ", paste(clash, collapse = ", "), " exist both as numbers and with a 'C' prefix: rename them before running CellChat")
  labels
}

# Create CellChat object with support for multiple grouping columns
## Group: Create CellChat object with support for multiple grouping columns
create_cellchat_object_multicol <- function(seurat_obj,
                                            grouping_columns,
                                            selected_groups = NULL,
                                            condition_column = NULL,
                                            selected_conditions = NULL,
                                            log_function = message) {
  if (length(grouping_columns) == 0) return(list(success = FALSE, message = "No grouping columns specified"))
  for (col in grouping_columns) {
    if (!col %in% colnames(seurat_obj@meta.data)) return(list(success = FALSE, message = paste("Column not found:", col)))
  }
  # Filter by additional condition if specified
  if (!is.null(condition_column) && !is.null(selected_conditions) && length(selected_conditions) > 0) {
    log_function(sprintf("Filtering by %s: %s", condition_column, paste(selected_conditions, collapse = ", ")))
    cells_to_keep <- which(as.character(seurat_obj@meta.data[[condition_column]]) %in% as.character(selected_conditions))
    if (length(cells_to_keep) == 0) return(list(success = FALSE, message = "No cells match the selected conditions"))
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
    log_function(sprintf("After filtering: %d cells remaining", ncol(seurat_obj)))
  }
  # Group labels: a cell missing a value in any grouping column gets NA and is removed by prepare_seurat_for_cellchat
  if (length(grouping_columns) > 1) {
    log_function(sprintf("Creating combined grouping from %d columns", length(grouping_columns)))
    group_by <- "CellChat_Combined_Groups"
  } else {
    group_by <- grouping_columns[1]
  }
  seurat_obj@meta.data[[group_by]] <- build_cellchat_group_labels(seurat_obj@meta.data, grouping_columns)
  # Filter specific groups if specified
  if (!is.null(selected_groups) && length(selected_groups) > 0) {
    log_function(sprintf("Keeping only %d selected groups", length(selected_groups)))
    cells_to_keep <- which(seurat_obj@meta.data[[group_by]] %in% selected_groups)
    if (length(cells_to_keep) == 0) return(list(success = FALSE, message = "No cells match the selected groups"))
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[cells_to_keep])
    log_function(sprintf("After group filtering: %d cells remaining", ncol(seurat_obj)))
  }
  # JoinLayers, removal of unlabeled cells (NA included) and droplevels
  seurat_obj <- prepare_seurat_for_cellchat(seurat_obj, group_by)
  if (isTRUE(seurat_obj@misc$cellchat_removed_cells > 0)) {
    log_function(sprintf("%d cells without a complete group label removed", seurat_obj@misc$cellchat_removed_cells))
  }
  if (ncol(seurat_obj) < 50) return(list(success = FALSE, message = paste("Too few cells:", ncol(seurat_obj), "cells remaining")))
  n_groups <- length(unique(seurat_obj@meta.data[[group_by]]))
  if (n_groups < 2) return(list(success = FALSE, message = paste("Need at least 2 groups, found:", n_groups)))
  log_function(sprintf("Creating CellChat object: %d cells, %d groups", ncol(seurat_obj), n_groups))
  data_input <- LayerData(seurat_obj, assay = "RNA", layer = "data")
  meta_data <- seurat_obj@meta.data
  cellchat_obj <- createCellChat(object = data_input, meta = meta_data, group.by = group_by)
  cellchat_obj <- clean_cellchat_factors(cellchat_obj)
  # Propagate Seurat cluster colors into CellChat @options for downstream plots
  cellchat_obj <- store_colors_in_cellchat(cellchat_obj, seurat_obj, group_by = group_by,
                                           primary_col = if (length(grouping_columns) > 1) grouping_columns[1] else NULL)
  log_function("CellChat object created successfully")
  return(list(success = TRUE, object = cellchat_obj))
}


## Group: Generate split chord plots: one plot per group in split_by,
## Group: using selected LR pairs. Sectors are shared across all panels so absent
## Group: interactions appear as gray ghosts, enabling direct visual comparison.
## Group: Colors are resolved from: (1) custom_colors arg, (2) @options$cluster_colors
## Group: stored at CellChat creation time, (3) deterministic brewer/Polychrome fallback.
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
  if (is.null(selected_lr_keys) || length(selected_lr_keys) == 0) stop("Please select at least one LR pair")
  # Determine split_groups and how to fetch comms per group
  use_pairs <- !is.null(pairs) && length(pairs) >= 2
  split_groups <- if (use_pairs) sapply(pairs, `[[`, "label") else if (split_by == "source") sources else targets
  if (anyDuplicated(split_groups)) stop("Group labels must be unique: ", paste(unique(split_groups[duplicated(split_groups)]), collapse = ", "))
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
      subsetCommunication(cellchat_obj, sources.use = s_use, targets.use = t_use, thresh = threshold),
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
                        paste(sprintf("  %s:%s -> %s:%s (%.5f)", top_show$source, top_show$ligand,
                                      top_show$receptor, top_show$target, top_show$prob), collapse = "\n")))
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
  union_comms$lr_key   <- paste0(union_comms$source, "|", union_comms$ligand, "->", union_comms$receptor, "|", union_comms$target)
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
  if (!is.null(custom_colors) && !is.null(names(custom_colors)) && sum(all_cts %in% names(custom_colors)) > 0) {
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
  # --- Build sectors: cell type and gene of each sector are stored, never parsed back from the sector name
  # (complex receptors such as Itga1_Itgb1 and cell types with underscores would break a regex) ---
  sector_table <- function(cts, genes) {
    genes_by_ct <- lapply(split(genes, cts), unique)
    data.frame(ct = rep(names(genes_by_ct), lengths(genes_by_ct)), gene = unlist(genes_by_ct, use.names = FALSE), stringsAsFactors = FALSE)
  }
  lig_tab <- sector_table(union_comms$source, union_comms$ligand)
  rec_tab <- sector_table(union_comms$target, union_comms$receptor)
  lig_sectors <- paste0(lig_tab$ct, "_L_", lig_tab$gene)
  rec_sectors <- paste0(rec_tab$ct, "_R_", rec_tab$gene)
  all_sectors <- c(lig_sectors, rec_sectors)
  sector_gene <- setNames(c(lig_tab$gene, rec_tab$gene), all_sectors)
  sector_ct   <- setNames(c(lig_tab$ct, rec_tab$ct), all_sectors)
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
  # One small gap between the genes of a cell type and a large gap after each cell type, in sector order
  src_cts  <- unique(lig_tab$ct)
  tgt_cts  <- unique(rec_tab$ct)
  src_lens <- as.integer(table(factor(lig_tab$ct, levels = src_cts)))
  tgt_lens <- as.integer(table(factor(rec_tab$ct, levels = tgt_cts)))
  gap_lig  <- as.numeric(unlist(lapply(src_lens, function(l) c(rep(gap_inner, max(l - 1, 0)), gap_outer))))
  gap_rec  <- as.numeric(unlist(lapply(tgt_lens, function(l) c(rep(gap_inner, max(l - 1, 0)), gap_outer))))
  gaps     <- as.numeric(c(gap_lig, gap_rec))
  n  <- length(split_groups)
  sg <- split_groups
  ld <- label_distance
  st <- show_title
  render_one_plot <- function(grp, grp_comms, grp_keys) {
    canvas_margin <- ld + 0.35
    canvas_lim    <- c(-1 - canvas_margin, 1 + canvas_margin)
    circlize::circos.clear()
    circlize::circos.par(gap.degree = gaps, cell.padding = c(0, 0, 0, 0), start.degree = 90, clock.wise = TRUE,
                         track.margin = c(0.005, 0.005), points.overflow.warning = FALSE,
                         canvas.xlim = canvas_lim, canvas.ylim = canvas_lim)
    circlize::circos.initialize(factors = factor(all_sectors, levels = all_sectors),
                                xlim = cbind(rep(0, length(all_sectors)), sector_sizes[all_sectors]))
    circlize::circos.trackPlotRegion(ylim = c(0, 1), bg.col = sector_col, bg.border = "white",
                                     track.height = 0.06, panel.fun = function(x, y) {})
    for (s in all_sectors) {
      xlim <- circlize::get.cell.meta.data("xlim", sector.index = s, track.index = 1)
      circlize::circos.text(x = mean(xlim), y = 1 + ld, labels = sector_gene[s], sector.index = s, track.index = 1,
                            facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = label_size, col = "black", font = 3)
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
      circlize::circos.link(ls, c(lstart, lend), rs, c(rstart, rend), col = adjustcolor("white", alpha.f = 0),
                            border = adjustcolor("gray75", alpha.f = 0.45), lwd = 0.4)
      pos[ls] <- lend; pos[rs] <- rend
    }
    if (!is.null(grp_comms) && nrow(grp_comms) > 0) {
      grp_sorted <- grp_comms[order(grp_comms$prob, decreasing = FALSE), ]
      for (j in seq_len(nrow(grp_sorted))) {
        ls  <- paste0(grp_sorted$source[j], "_L_", grp_sorted$ligand[j])
        rs  <- paste0(grp_sorted$target[j], "_R_", grp_sorted$receptor[j])
        src <- grp_sorted$source[j]
        if (!ls %in% all_sectors || !rs %in% all_sectors) next
        arc_col    <- adjustcolor(if (src %in% names(ct_colors)) ct_colors[[src]] else "gray40", alpha.f = 0.75)
        border_col <- if (src %in% names(ct_colors)) ct_colors[[src]] else "gray40"
        lstart <- pos[ls]; lend <- min(lstart + grp_sorted$prob[j], sector_sizes[ls])
        rstart <- pos[rs]; rend <- min(rstart + grp_sorted$prob[j], sector_sizes[rs])
        if ((lend - lstart) < 1e-10 || (rend - rstart) < 1e-10) next
        circlize::circos.link(ls, c(lstart, lend), rs, c(rstart, rend), col = arc_col, border = border_col,
                              directional = 1, arr.type = "big.arrow", arr.length = 0.04, lwd = 0.5)
        pos[ls] <- lend; pos[rs] <- rend
      }
    }
    if (st) title(main = grp, cex.main = 1.8, font.main = 2, line = -1)
    legend(x = canvas_lim[2] - 0.05, y = canvas_lim[1] + 0.05, xjust = 1, yjust = 0,
           legend = c(names(ct_colors), "Absent"), fill = c(unname(ct_colors), adjustcolor("gray75", alpha.f = 0.5)),
           title = "Cell Type", cex = 0.75, bty = "n", border = NA)
    circlize::circos.clear()
  }
  # --- Render one PNG per group then stitch side-by-side ---
  tmp_files <- vapply(seq_len(n), function(i) {
    grp       <- sg[i]
    grp_comms <- comms_list[[grp]]
    grp_keys  <- if (!is.null(grp_comms) && nrow(grp_comms) > 0) {
      paste0(grp_comms$source, "|", grp_comms$ligand, "->", grp_comms$receptor, "|", grp_comms$target)
    } else character(0)
    n_present <- length(grp_keys)
    n_absent  <- sum(!union_comms$lr_key %in% grp_keys)
    message(sprintf("--- [%d/%d] %s | present: %d | absent (ghost): %d", i, n, grp, n_present, n_absent))
    if (n_present > 0) message(sprintf("  present:\n%s", paste(" ", grp_keys, collapse = "\n")))
    if (n_absent > 0) message(sprintf("  ghost (first 10):\n%s", paste(" ", head(union_comms$lr_key[!union_comms$lr_key %in% grp_keys], 10), collapse = "\n")))
    tmp <- tempfile(fileext = ".png")
    png(tmp, width = plot_size, height = plot_size, res = 150)
    png_dev <- dev.cur()
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
    # Close this exact PNG device (never another open device)
    if (png_dev %in% dev.list()) dev.off(png_dev)
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
        paste0(grp_comms$source, "|", grp_comms$ligand, "->", grp_comms$receptor, "|", grp_comms$target)
      } else character(0)
      tryCatch(render_one_plot(grp, grp_comms, grp_keys),
               error = function(e) { plot.new(); text(0.5, 0.5, paste0(grp, "\nError: ", e$message), col = "red") })
    }
  }
  return(list(plot_function = plot_function, export_function = export_function, output_path = output_path,
              n_plots = n, selected_lr_keys = selected_lr_keys, png_size = plot_size))
}


## Group: Store cluster colors from a Seurat object into a CellChat object's @options slot.
## Group: Handles combined group names (e.g., "B_cells_WT") by falling back to the primary
## Group: group color when an exact match is not found.
## Group: Args:
## Group:   cellchat_obj : CellChat object to annotate
## Group:   seurat_obj   : source Seurat object containing @misc$cluster_colors
## Group:   group_by     : column name used as grouping (may be a combined column)
## Group:   primary_col  : original grouping column before combination (for fallback matching)
## Group: Returns: CellChat object with colors stored in @options$cluster_colors
store_colors_in_cellchat <- function(cellchat_obj, seurat_obj,
                                     group_by, primary_col = NULL) {
  group_levels <- levels(cellchat_obj@idents)
  if (is.null(group_levels) || length(group_levels) == 0) return(cellchat_obj)
  n <- length(group_levels)
  resolved <- setNames(rep(NA_character_, n), group_levels)
  # Try to get Seurat cluster colors from @misc$cluster_colors
  seurat_colors <- tryCatch(seurat_obj@misc$cluster_colors, error = function(e) NULL)
  if (!is.null(seurat_colors) && length(seurat_colors) > 0) {
    # Direct name match
    direct_match <- group_levels %in% names(seurat_colors)
    resolved[direct_match] <- seurat_colors[group_levels[direct_match]]
    # Fallback for combined groups: the longest cell type name followed by "_" ("Myonuclei_WT" -> "Myonuclei", never "Myo")
    if (!is.null(primary_col) && any(is.na(resolved))) {
      ct_names <- names(seurat_colors)[order(-nchar(names(seurat_colors)))]
      for (lvl in group_levels[is.na(resolved)]) {
        hit <- ct_names[startsWith(lvl, paste0(ct_names, "_"))]
        if (length(hit) > 0) resolved[[lvl]] <- seurat_colors[[hit[1]]]
      }
    }
    message(sprintf("[store_colors_in_cellchat] %d/%d matched from @misc$cluster_colors", sum(!is.na(resolved)), n))
  } else {
    # No Cell-Hub colors stored — reconstruct Seurat's default hue_pal palette,
    # which matches what Seurat renders on the UMAP when colors were never customized.
    message("[store_colors_in_cellchat] No @misc$cluster_colors — reconstructing Seurat default palette")
    seurat_default <- setNames(scales::hue_pal()(n), group_levels)
    default_match  <- group_levels %in% names(seurat_default)
    resolved[default_match] <- seurat_default[group_levels[default_match]]
  }
  # Fill any remaining NAs (partial miss or combined-group mismatch)
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
    message(sprintf("[store_colors_in_cellchat] %d color(s) filled with fallback palette", n_missing))
  }
  cellchat_obj@options$cluster_colors <- resolved
  message(sprintf("[store_colors_in_cellchat] Stored %d colors total", n))
  return(cellchat_obj)
}

# Source priority used by the CellHubDB curation pipeline: the highest value wins when a pair is found in several sources
CELLHUBDB_SOURCE_PRIORITY <- c(MultiNicheNet = 1, CellTalkDB = 2, cell2cell_Armingol = 3, PlantPhoneDB = 3, FlyPhoneDB2 = 3,
                               PlantCellChatDB = 4, NeuronChatDB = 6, CellPhoneDB = 7, connectomeDB2025 = 8, CellChatDB_v2 = 10)
# Sources whose pairs do not encode the ligand -> receptor direction
CELLHUBDB_NON_ORIENTED <- c("CellPhoneDB", "PlantCellChatDB")

# Row-bind tables keyed by row names (complex, cofactor), filling missing columns with "" and keeping the first definition of each row
bind_named_tables <- function(tables) {
  tables <- Filter(function(t) !is.null(t) && nrow(t) > 0, tables)
  if (length(tables) == 0) return(NULL)
  cols <- unique(unlist(lapply(tables, colnames)))
  row_ids <- unlist(lapply(tables, rownames), use.names = FALSE)
  tables <- lapply(tables, function(t) {
    t[] <- lapply(t, as.character)  # factors from CellChatDB would break rbind with "" fillers
    for (cl in setdiff(cols, colnames(t))) t[[cl]] <- rep("", nrow(t))
    t[, cols, drop = FALSE]
  })
  out <- do.call(rbind, unname(tables))
  keep <- !duplicated(row_ids)
  out <- out[keep, , drop = FALSE]
  rownames(out) <- row_ids[keep]
  out
}

# Keep the complexes used by the interaction partners and define undefined heteromers from their underscore-separated subunits
# (CellChat silently ignores a multi-subunit partner without a complex definition)
complete_complex_table <- function(complex_df, partner_names) {
  partner_names <- unique(partner_names[!is.na(partner_names)])
  if (is.null(complex_df) || nrow(complex_df) == 0) {
    complex_df <- data.frame(subunit_1 = character(0), subunit_2 = character(0), stringsAsFactors = FALSE)
  }
  complex_df <- complex_df[rownames(complex_df) %in% partner_names, , drop = FALSE]
  # Metabolite composites ("-") and named complexes with spaces need a curated definition and are never generated
  undefined <- partner_names[grepl("_", partner_names) & !grepl("[- ]", partner_names) & !(partner_names %in% rownames(complex_df))]
  if (length(undefined) == 0) return(complex_df)
  parts <- strsplit(undefined, "_")
  n_cols <- max(ncol(complex_df), max(lengths(parts)))
  for (k in seq_len(n_cols)) {
    cl <- paste0("subunit_", k)
    if (!cl %in% colnames(complex_df)) complex_df[[cl]] <- rep("", nrow(complex_df))
  }
  generated <- as.data.frame(do.call(rbind, lapply(parts, function(p) c(p, rep("", n_cols - length(p))))), stringsAsFactors = FALSE)
  colnames(generated) <- paste0("subunit_", seq_len(n_cols))
  rownames(generated) <- undefined
  for (cl in setdiff(colnames(complex_df), colnames(generated))) generated[[cl]] <- ""
  rbind(complex_df, generated[, colnames(complex_df), drop = FALSE])
}

# Merge several CellHubDB source databases of one species with the rules of curate_CellHub_DB.R:
# orientation-aware case-insensitive key, non-oriented sources aligned on oriented ones (or collapsed when reciprocal),
# highest-priority record kept, provenance/evidence aggregated, NA metadata filled from lower-priority sources.
# db_list: named list (names = source names) of CellChat-format lists (interaction, complex, cofactor, geneInfo)
merge_cellhubdb_sources <- function(db_list, priority = CELLHUBDB_SOURCE_PRIORITY, non_oriented = CELLHUBDB_NON_ORIENTED) {
  if (length(db_list) == 0) stop("No source database provided")
  if (is.null(names(db_list)) || any(names(db_list) == "")) stop("db_list must be named by source")
  inter_list <- lapply(names(db_list), function(src) {
    df <- db_list[[src]]$interaction
    if (is.null(df$source_db)) df$source_db <- src
    if (is.null(df$sources)) df$sources <- df$source_db
    df
  })
  cols <- unique(unlist(lapply(inter_list, colnames)))
  inter_list <- lapply(inter_list, function(df) {
    for (cl in setdiff(cols, colnames(df))) df[[cl]] <- NA
    df[, cols, drop = FALSE]
  })
  df <- do.call(rbind, unname(inter_list))
  rownames(df) <- NULL
  n_input <- nrow(df)
  # Orientation-aware key; non-oriented pairs are flipped onto an oriented counterpart when one exists
  lig <- tolower(df$ligand); rec <- tolower(df$receptor)
  key_fwd <- paste(lig, rec, sep = "|||"); key_rev <- paste(rec, lig, sep = "|||")
  unoriented <- df$source_db %in% non_oriented
  oriented_keys <- key_fwd[!unoriented]
  to_flip <- unoriented & !(key_fwd %in% oriented_keys) & (key_rev %in% oriented_keys)
  if (any(to_flip)) {
    old_lig <- df$ligand[to_flip]
    df$ligand[to_flip] <- df$receptor[to_flip]
    df$receptor[to_flip] <- old_lig
    key_fwd[to_flip] <- key_rev[to_flip]
  }
  # Remaining non-oriented pairs use an order-free key so that A-B and B-A collapse into one entry
  free <- unoriented & !(key_fwd %in% oriented_keys)
  key_fwd[free] <- ifelse(lig[free] < rec[free], key_fwd[free], key_rev[free])
  df$merge_key <- key_fwd
  df$priority <- ifelse(df$source_db %in% names(priority), priority[df$source_db], 0)
  df <- df[order(-df$priority, seq_len(nrow(df))), ]
  # Provenance and evidence aggregated over all entries sharing a key, in priority order
  src_by_key <- tapply(df$sources, df$merge_key, function(s) {
    paste(unique(trimws(unlist(strsplit(s[!is.na(s)], ",")))), collapse = ", ")
  })
  merged <- df[!duplicated(df$merge_key), ]
  merged$sources <- unname(src_by_key[merged$merge_key])
  if ("evidence" %in% cols) {
    ev_by_key <- tapply(df$evidence, df$merge_key, function(e) {
      e <- unique(e[!is.na(e)])
      if (length(e) == 0) NA_character_ else paste(e, collapse = "; ")
    })
    merged$evidence <- unname(ev_by_key[merged$merge_key])
  }
  # Fill missing metadata of the retained record from the highest-priority source that has it
  for (cl in intersect(c("pathway_name", "annotation", "interaction_type"), cols)) {
    has <- df[!is.na(df[[cl]]), c("merge_key", cl)]
    has <- has[!duplicated(has$merge_key), ]
    fill <- is.na(merged[[cl]])
    merged[[cl]][fill] <- has[[cl]][match(merged$merge_key[fill], has$merge_key)]
  }
  merged$interaction_name <- paste(merged$ligand, merged$receptor, sep = "_")
  merged$interaction_name_2 <- ifelse(grepl("_", merged$receptor),
                                      paste0(merged$ligand, " - (", gsub("_", "+", merged$receptor), ")"),
                                      paste(merged$ligand, "-", merged$receptor))
  merged <- merged[!duplicated(merged$interaction_name), ]
  merged$merge_key <- NULL
  merged$priority <- NULL
  rownames(merged) <- merged$interaction_name
  # Complex, cofactor and geneInfo tables pooled from the selected sources
  complex_df <- complete_complex_table(bind_named_tables(lapply(db_list, `[[`, "complex")), c(merged$ligand, merged$receptor))
  cofactor_df <- bind_named_tables(lapply(db_list, `[[`, "cofactor"))
  if (is.null(cofactor_df)) cofactor_df <- data.frame(cofactor1 = character(0), stringsAsFactors = FALSE)
  gene_info <- Filter(function(g) !is.null(g) && nrow(g) > 0, lapply(db_list, `[[`, "geneInfo"))
  gene_info <- if (length(gene_info) > 0) gene_info[[1]] else data.frame(Symbol = character(0), stringsAsFactors = FALSE)
  list(db = list(interaction = merged, complex = complex_df, cofactor = cofactor_df, geneInfo = gene_info),
       stats = list(n_input = n_input, n_merged = nrow(merged), n_multi = sum(grepl(",", merged$sources)), n_flipped = sum(to_flip)))
}

# Make a CellHubDB species database CellChat-ready whatever release it comes from: full CellChat column schema for the complex
# and cofactor tables (a one-column cofactor table collapses to a vector in subsetData), missing cofactor references set to "",
# missing pathway labels set to a placeholder. Stops when geneInfo is empty (CellChat keeps only geneInfo genes, so every
# interaction would be dropped) or when fewer than two pathways remain (computeCommunProbPathway needs at least two).
prepare_cellchat_db <- function(db, placeholder = "Unassigned") {
  ref <- CellChat::CellChatDB.mouse
  msgs <- character()
  fix_table <- function(tab, cols) {
    if (is.null(tab) || ncol(tab) == 0) tab <- data.frame(row.names = character(0))
    tab[] <- lapply(tab, as.character)
    for (cl in setdiff(cols, colnames(tab))) tab[[cl]] <- rep("", nrow(tab))
    tab <- tab[, c(cols, setdiff(colnames(tab), cols)), drop = FALSE]
    tab[is.na(tab)] <- ""
    return(tab)
  }
  db$complex <- fix_table(db$complex, union(colnames(ref$complex), grep("^subunit_", colnames(db$complex), value = TRUE)))
  db$cofactor <- fix_table(db$cofactor, union(colnames(ref$cofactor), grep("^cofactor", colnames(db$cofactor), value = TRUE)))
  inter <- db$interaction
  inter[] <- lapply(inter, function(x) if (is.factor(x)) as.character(x) else x)
  for (cl in c("agonist", "antagonist", "co_A_receptor", "co_I_receptor")) {
    if (!cl %in% colnames(inter)) inter[[cl]] <- rep("", nrow(inter))
    inter[[cl]][is.na(inter[[cl]])] <- ""
  }
  if (is.null(db$geneInfo) || !"Symbol" %in% colnames(db$geneInfo) || nrow(db$geneInfo) == 0) {
    stop("This database has no gene reference table (geneInfo), so CellChat cannot use it. Please update the CellHub Database files (v2.0.1 or later).")
  }
  missing_pw <- is.na(inter$pathway_name) | inter$pathway_name == ""
  if (any(missing_pw)) {
    inter$pathway_name[missing_pw] <- placeholder
    msgs <- c(msgs, sprintf("%d interactions without pathway labelled '%s'", sum(missing_pw), placeholder))
  }
  n_pathways <- length(unique(inter$pathway_name))
  if (n_pathways < 2) stop(sprintf("CellChat needs at least two pathways; this database has %d. Please select another database or source combination.", n_pathways))
  db$interaction <- inter
  return(list(db = db, messages = msgs))
}