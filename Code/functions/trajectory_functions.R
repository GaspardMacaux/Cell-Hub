# 3D Trajectory Navigator with Pseudotime Slider
# Add to visualization_workaround_functions.R or create new file

############################## 3D Trajectory Slider Function ##############################

create3DTrajectoryNavigator <- function(monocle_object,
                                        seurat_object = NULL,
                                        reduction = "umap3d",
                                        pseudotime_position = 0.5,
                                        window_size = 0.05,
                                        show_trajectory_line = TRUE,
                                        dark_mode = FALSE) {
  
  pt_values <- pseudotime(monocle_object)
  
  if(all(is.na(pt_values))) {
    stop("No valid pseudotime values found.")
  }
  
  pt_min <- min(pt_values, na.rm = TRUE)
  pt_max <- max(pt_values, na.rm = TRUE)
  pt_norm <- (pt_values - pt_min) / (pt_max - pt_min)
  
  if(reduction %in% names(monocle_object@principal_graph_aux)) {
    coords <- reducedDims(monocle_object)[[reduction]]
  } else if(!is.null(seurat_object) && reduction %in% names(seurat_object@reductions)) {
    coords <- Embeddings(seurat_object, reduction = reduction)
  } else {
    stop(paste("3D reduction", reduction, "not found"))
  }
  
  if(ncol(coords) < 3) {
    stop("Reduction must have 3 dimensions")
  }
  
  coords_3d <- coords[, 1:3]
  
  in_window <- abs(pt_norm - pseudotime_position) < window_size
  in_window[is.na(in_window)] <- FALSE
  
  colors <- ifelse(in_window, "#FF0000", "#CCCCCC")
  sizes <- ifelse(in_window, 5, 2)
  opacities <- ifelse(in_window, 1, 0.3)
  
  bg_color <- if(dark_mode) "#1a1a1a" else "white"
  grid_color <- if(dark_mode) "#404040" else "#e0e0e0"
  axis_color <- if(dark_mode) "#ffffff" else "#000000"
  text_color <- if(dark_mode) "#ffffff" else "#000000"
  
  hover_text <- paste0(
    "Pseudotime: ", round(pt_values, 3)
  )
  
  # Create plot - CELLS ONLY, no lines
  fig <- plot_ly() %>%
    add_trace(
      x = coords_3d[, 1],
      y = coords_3d[, 2],
      z = coords_3d[, 3],
      type = "scatter3d",
      mode = "markers",  # MARKERS ONLY - no lines!
      marker = list(
        color = colors,
        size = sizes,
        opacity = opacities,
        line = list(width = 0)
      ),
      text = hover_text,
      hoverinfo = "text",
      name = "Cells",
      showlegend = FALSE
    )
  
  # Add trajectory LINE if requested
  if(show_trajectory_line) {
    valid_pt <- !is.na(pt_values)
    order_idx <- order(pt_values[valid_pt])
    trajectory_coords <- coords_3d[valid_pt, ][order_idx, ]
    
    # Sample for smooth line
    n_points <- min(200, nrow(trajectory_coords))
    sample_idx <- seq(1, nrow(trajectory_coords), length.out = n_points)
    sample_idx <- round(sample_idx)
    trajectory_smooth <- trajectory_coords[sample_idx, ]
    
    # Add as SEPARATE trace with mode="lines" ONLY
    fig <- fig %>%
      add_trace(
        x = trajectory_smooth[, 1],
        y = trajectory_smooth[, 2],
        z = trajectory_smooth[, 3],
        type = "scatter3d",
        mode = "lines",  # LINES ONLY - no markers!
        line = list(
          color = "#FF6600",
          width = 6
        ),
        name = "Trajectory",
        showlegend = FALSE,
        hoverinfo = "skip"
      )
  }
  
  # Axes config
  axis_config <- list(
    showgrid = TRUE,
    showline = TRUE,
    zeroline = FALSE,
    showticklabels = TRUE,
    gridcolor = grid_color,
    linecolor = axis_color,
    tickcolor = axis_color,
    tickfont = list(color = text_color)
  )
  
  xaxis_config <- axis_config
  yaxis_config <- axis_config
  zaxis_config <- axis_config
  xaxis_config$title <- list(text = paste0(toupper(reduction), "_1"), font = list(color = text_color))
  yaxis_config$title <- list(text = paste0(toupper(reduction), "_2"), font = list(color = text_color))
  zaxis_config$title <- list(text = paste0(toupper(reduction), "_3"), font = list(color = text_color))
  
  fig <- fig %>%
    layout(
      scene = list(
        xaxis = xaxis_config,
        yaxis = yaxis_config,
        zaxis = zaxis_config,
        camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5)),
        bgcolor = bg_color
      ),
      paper_bgcolor = bg_color,
      plot_bgcolor = bg_color,
      title = list(
        text = paste0("Trajectory Navigator - Pseudotime: ", round(pseudotime_position, 3)),
        x = 0.5,
        font = list(color = text_color, size = 16)
      ),
      font = list(color = text_color),
      showlegend = FALSE
    )
  
  fig <- toWebGL(fig)
  
  # Calculate top genes
  cells_in_window <- which(in_window)
  top_genes <- NULL
  
  if(length(cells_in_window) >= 10) {
    expr_matrix <- exprs(monocle_object)
    window_expr <- rowMeans(as.matrix(expr_matrix[, cells_in_window]))
    top_genes <- names(sort(window_expr, decreasing = TRUE)[1:15])
  }
  
  return(list(
    plot = fig,
    top_genes = top_genes,
    n_cells_highlighted = sum(in_window),
    pseudotime_range = c(pt_min, pt_max)
  ))
}
############################## Helper: Get trajectory stats ##############################

getTrajectoryStats <- function(monocle_object) {
  # Get statistics about trajectory
  # Returns: List with pseudotime range, n_cells, etc.
  
  pt_values <- pseudotime(monocle_object)
  valid_pt <- !is.na(pt_values)
  
  return(list(
    n_cells_total = length(pt_values),
    n_cells_with_pseudotime = sum(valid_pt),
    pseudotime_min = min(pt_values, na.rm = TRUE),
    pseudotime_max = max(pt_values, na.rm = TRUE),
    pseudotime_mean = mean(pt_values, na.rm = TRUE),
    pseudotime_median = median(pt_values, na.rm = TRUE)
  ))
}