# visualization_workaround_functions.R
# Manual implementations of VlnPlot and RidgePlot to work around Seurat v5 + Docker issues

############################## Manual VlnPlot Function ##############################

createManualVlnPlot <- function(seurat_object,
                                genes,
                                group_by = "ident",
                                pt_size = 1,
                                axis_text_size = 12,
                                title_text_size = 14,
                                axis_line_width = 1,
                                add_noaxes = FALSE,
                                add_nolegend = FALSE,
                                plot_title = NULL) {
  # Create manual violin plot using ggplot2
  # Args:
  #   seurat_object: Seurat object
  #   genes: Character vector of gene names
  #   group_by: Column name for grouping (default: "ident")
  #   pt_size: Point size (0 to hide points)
  #   axis_text_size: Size of axis text
  #   title_text_size: Size of title text
  #   axis_line_width: Width of axis lines
  #   add_noaxes: Remove axes
  #   add_nolegend: Remove legend
  #   plot_title: Custom plot title
  # Returns:
  #   ggplot object
  
  # Fetch data
  if(group_by == "ident") {
    plot_data <- FetchData(seurat_object, vars = c(genes, "ident"))
  } else {
    plot_data <- FetchData(seurat_object, vars = c(genes, group_by))
    colnames(plot_data)[colnames(plot_data) == group_by] <- "ident"
  }
  
  # Handle multiple genes
  if(length(genes) > 1) {
    plot_data_long <- tidyr::pivot_longer(
      plot_data, 
      cols = all_of(genes),
      names_to = "gene",
      values_to = "expression"
    )
    
    plot <- ggplot(plot_data_long, aes(x = ident, y = expression, fill = ident)) +
      geom_violin(scale = "width", trim = TRUE) +
      facet_wrap(~gene, scales = "free_y", ncol = 2) +
      labs(x = NULL, y = "Expression Level") +
      theme_classic() +
      theme(
        axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = title_text_size),
        plot.title = element_text(size = title_text_size),
        legend.text = element_text(size = axis_text_size),
        axis.line = element_line(linewidth = axis_line_width),
        axis.ticks = element_line(linewidth = axis_line_width),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    # Add points if requested
    if(pt_size > 0) {
      plot <- plot + geom_jitter(height = 0, width = 0.2, alpha = 0.5, size = pt_size)
    }
    
  } else {
    # Single gene
    plot <- ggplot(plot_data, aes(x = ident, y = .data[[genes[1]]], fill = ident)) +
      geom_violin(scale = "width", trim = TRUE) +
      labs(x = NULL, y = "Expression Level", title = ifelse(is.null(plot_title), genes[1], plot_title)) +
      theme_classic() +
      theme(
        axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = title_text_size),
        plot.title = element_text(size = title_text_size, hjust = 0.5),
        legend.text = element_text(size = axis_text_size),
        axis.line = element_line(linewidth = axis_line_width),
        axis.ticks = element_line(linewidth = axis_line_width),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    # Add points if requested
    if(pt_size > 0) {
      plot <- plot + geom_jitter(height = 0, width = 0.2, alpha = 0.5, size = pt_size)
    }
  }
  
  # Apply NoAxes if requested
  if(add_noaxes) {
    plot <- plot + theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  }
  
  # Apply NoLegend if requested
  if(add_nolegend) {
    plot <- plot + theme(legend.position = "none")
  }
  
  return(plot)
}


############################## Manual VlnPlot for Multiple Datasets ##############################

createManualVlnPlotMultiple <- function(seurat_object,
                                        genes,
                                        group_by = "dataset",
                                        pt_size = 0.1,
                                        add_noaxes = FALSE,
                                        add_nolegend = FALSE,
                                        plot_title = NULL,
                                        ncol = 2) {
  # Create manual violin plot for multiple datasets module
  # Args:
  #   seurat_object: Seurat object
  #   genes: Character vector of gene names
  #   group_by: Column name for grouping
  #   pt_size: Point size (0 to hide points)
  #   add_noaxes: Remove axes
  #   add_nolegend: Remove legend
  #   plot_title: Custom plot title
  #   ncol: Number of columns for facets
  # Returns:
  #   ggplot object or patchwork object (for multiple genes)
  
  # Fetch data
  plot_data <- FetchData(seurat_object, vars = c(genes, group_by))
  
  # Handle multiple genes
  if(length(genes) > 1) {
    # Create individual plots for each gene
    plot_list <- lapply(genes, function(gene) {
      p <- ggplot(plot_data, aes(x = .data[[group_by]], y = .data[[gene]], fill = .data[[group_by]])) +
        geom_violin(scale = "width", trim = TRUE) +
        labs(x = NULL, y = "Expression Level", title = gene) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.margin = margin(t = 20, r = 20, b = 60, l = 20)
        )
      
      # Add points if requested
      if(pt_size > 0) {
        p <- p + geom_jitter(height = 0, width = 0.2, alpha = 0.5, size = pt_size)
      }
      
      # Apply options - ONLY if explicitly requested
      if(add_noaxes) {
        p <- p + theme(
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()
        )
      }
      if(add_nolegend) {  # Only remove if checkbox is TRUE
        p <- p + theme(legend.position = "none")
      }
      
      return(p)
    })
    
    # Combine plots
    plot <- wrap_plots(plot_list, ncol = ncol)
    attr(plot, "n_rows") <- ceiling(length(genes) / ncol)
    
  } else {
    # Single gene
    plot <- ggplot(plot_data, aes(x = .data[[group_by]], y = .data[[genes[1]]], fill = .data[[group_by]])) +
      geom_violin(scale = "width", trim = TRUE) +
      labs(x = NULL, y = "Expression Level", title = ifelse(is.null(plot_title), paste("Expression of", genes[1]), plot_title)) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(t = 20, r = 20, b = 60, l = 20)
      )
    
    # Add points if requested
    if(pt_size > 0) {
      plot <- plot + geom_jitter(height = 0, width = 0.2, alpha = 0.5, size = pt_size)
    }
    
    # Apply options - ONLY if explicitly requested
    if(add_noaxes) {
      plot <- plot + theme(
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()
      )
    }
    if(add_nolegend) {  # Only remove if checkbox is TRUE
      plot <- plot + theme(legend.position = "none")
    }
    
    attr(plot, "n_rows") <- 1
  }
  
  return(plot)
}

############################## Manual RidgePlot Function ##############################

createManualRidgePlot <- function(seurat_object,
                                  genes,
                                  group_by = "ident",
                                  axis_text_size = 12,
                                  title_text_size = 14,
                                  axis_line_width = 1,
                                  add_noaxes = FALSE,
                                  add_nolegend = FALSE) {
  # Create manual ridge plot using ggplot2 + ggridges
  # Args:
  #   seurat_object: Seurat object
  #   genes: Character vector of gene names
  #   group_by: Column name for grouping (default: "ident")
  #   axis_text_size: Size of axis text
  #   title_text_size: Size of title text
  #   axis_line_width: Width of axis lines
  #   add_noaxes: Remove axes
  #   add_nolegend: Remove legend
  # Returns:
  #   ggplot object
  
  # Load ggridges if not loaded
  if(!requireNamespace("ggridges", quietly = TRUE)) {
    stop("ggridges package is required for RidgePlot. Please install it.")
  }
  
  # Fetch data
  if(group_by == "ident") {
    plot_data <- FetchData(seurat_object, vars = c(genes, "ident"))
  } else {
    plot_data <- FetchData(seurat_object, vars = c(genes, group_by))
    colnames(plot_data)[colnames(plot_data) == group_by] <- "ident"
  }
  
  # Handle multiple genes
  if(length(genes) > 1) {
    plot_data_long <- tidyr::pivot_longer(
      plot_data, 
      cols = all_of(genes),
      names_to = "gene",
      values_to = "expression"
    )
    
    plot <- ggplot(plot_data_long, aes(x = expression, y = ident, fill = ident)) +
      ggridges::geom_density_ridges(scale = 1) +
      facet_wrap(~gene, scales = "free_x", ncol = 1) +
      labs(x = "Expression Level", y = NULL) +
      theme_classic() +
      theme(
        axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = title_text_size),
        plot.title = element_text(size = title_text_size),
        legend.text = element_text(size = axis_text_size),
        axis.line = element_line(linewidth = axis_line_width),
        axis.ticks = element_line(linewidth = axis_line_width)
      )
    
  } else {
    # Single gene
    plot <- ggplot(plot_data, aes(x = .data[[genes[1]]], y = ident, fill = ident)) +
      ggridges::geom_density_ridges(scale = 1) +
      labs(x = "Expression Level", y = NULL, title = genes[1]) +
      theme_classic() +
      theme(
        axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = title_text_size),
        plot.title = element_text(size = title_text_size, hjust = 0.5),
        legend.text = element_text(size = axis_text_size),
        axis.line = element_line(linewidth = axis_line_width),
        axis.ticks = element_line(linewidth = axis_line_width)
      )
  }
  
  # Apply NoAxes if requested
  if(add_noaxes) {
    plot <- plot + theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  }
  
  # Apply NoLegend if requested
  if(add_nolegend) {
    plot <- plot + theme(legend.position = "none")
  }
  
  return(plot)
}


############################## Manual RidgePlot for Multiple Datasets ##############################

createManualRidgePlotMultiple <- function(seurat_object,
                                          genes,
                                          group_by = "dataset",
                                          add_noaxes = FALSE,
                                          add_nolegend = FALSE) {
  # Create manual ridge plot for multiple datasets module
  # Args:
  #   seurat_object: Seurat object
  #   genes: Character vector of gene names
  #   group_by: Column name for grouping
  #   add_noaxes: Remove axes
  #   add_nolegend: Remove legend
  # Returns:
  #   ggplot object
  
  # Load ggridges if not loaded
  if(!requireNamespace("ggridges", quietly = TRUE)) {
    stop("ggridges package is required for RidgePlot. Please install it.")
  }
  
  # Fetch data
  plot_data <- FetchData(seurat_object, vars = c(genes, group_by))
  
  # Handle multiple genes
  if(length(genes) > 1) {
    plot_data_long <- tidyr::pivot_longer(
      plot_data, 
      cols = all_of(genes),
      names_to = "gene",
      values_to = "expression"
    )
    
    plot <- ggplot(plot_data_long, aes(x = expression, y = .data[[group_by]], fill = .data[[group_by]])) +
      ggridges::geom_density_ridges(scale = 1) +
      facet_wrap(~gene, scales = "free_x", ncol = 1) +
      labs(x = "Expression Level", y = NULL) +
      theme_classic()
    
  } else {
    # Single gene
    plot <- ggplot(plot_data, aes(x = .data[[genes[1]]], y = .data[[group_by]], fill = .data[[group_by]])) +
      ggridges::geom_density_ridges(scale = 1) +
      labs(x = "Expression Level", y = NULL, title = genes[1]) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  # Apply options
  if(add_noaxes) {
    plot <- plot + theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  }
  if(add_nolegend) {
    plot <- plot + theme(legend.position = "none")
  }
  
  return(plot)
}


############################## 3D FeaturePlot Function ##############################


create3DFeaturePlot <- function(seurat_object, 
                                genes, 
                                reduction = "umap3d",  # Changed from umap_3d to umap3d
                                pt_size = 3,
                                min_cutoff = NA,
                                max_cutoff = NA,
                                order = TRUE,
                                hide_grid = FALSE,
                                hide_axes = FALSE,
                                dark_mode = FALSE) {
  # Create 3D FeaturePlot using plotly with display options
  # Args:
  #   seurat_object: Seurat object with 3D reduction
  #   genes: Character vector of gene names (1-3 genes supported)
  #   reduction: Name of dimensional reduction (must have 3 dimensions)
  #   pt_size: Point size
  #   min_cutoff: Minimum expression cutoff
  #   max_cutoff: Maximum expression cutoff
  #   order: Order points by expression (highest on top)
  #   hide_grid: Hide grid lines
  #   hide_axes: Hide axes
  #   dark_mode: Enable dark mode
  # Returns:
  #   plotly object
  
  # Validate reduction has 3 dimensions
  if(!reduction %in% names(seurat_object@reductions)) {
    stop(paste("Reduction", reduction, "not found in Seurat object"))
  }
  
  coords <- Embeddings(seurat_object, reduction = reduction)
  if(ncol(coords) < 3) {
    stop(paste("Reduction", reduction, "must have at least 3 dimensions for 3D plotting"))
  }
  
  # Use first 3 dimensions
  coords_3d <- coords[, 1:3]
  colnames(coords_3d) <- c("Dim1", "Dim2", "Dim3")
  
  # Set colors based on mode
  bg_color <- if(dark_mode) "#1a1a1a" else "white"
  grid_color <- if(dark_mode) "#404040" else "#e0e0e0"
  axis_color <- if(dark_mode) "#ffffff" else "#000000"
  text_color <- if(dark_mode) "#ffffff" else "#000000"
  
  # Configure axis settings
  axis_config <- list(
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
  
  # Handle multiple genes
  if(length(genes) == 1) {
    # Single gene - color by expression
    gene_expr <- FetchData(seurat_object, vars = genes[1])[[1]]
    
    # Apply cutoffs
    if(!is.na(min_cutoff)) {
      gene_expr[gene_expr < min_cutoff] <- min_cutoff
    }
    if(!is.na(max_cutoff)) {
      gene_expr[gene_expr > max_cutoff] <- max_cutoff
    }
    
    # Order cells by expression if requested
    if(order) {
      cell_order <- order(gene_expr)
      coords_3d <- coords_3d[cell_order, ]
      gene_expr <- gene_expr[cell_order]
    }
    
    # Create axis configs with titles
    xaxis_config <- axis_config
    yaxis_config <- axis_config
    zaxis_config <- axis_config
    
    if (!hide_axes) {
      xaxis_config$title <- list(text = paste0(toupper(reduction), "_1"), font = list(color = text_color))
      yaxis_config$title <- list(text = paste0(toupper(reduction), "_2"), font = list(color = text_color))
      zaxis_config$title <- list(text = paste0(toupper(reduction), "_3"), font = list(color = text_color))
    } else {
      xaxis_config$title <- ""
      yaxis_config$title <- ""
      zaxis_config$title <- ""
    }
    
    # Create plot
    fig <- plot_ly(
      x = coords_3d[, 1],
      y = coords_3d[, 2],
      z = coords_3d[, 3],
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = pt_size,
        color = gene_expr,
        colorscale = "Viridis",
        colorbar = list(
          title = list(text = genes[1], font = list(color = text_color)),
          tickfont = list(color = text_color)
        ),
        showscale = TRUE,
        line = list(width = 0)
      ),
      text = ~paste0(
        "Cell: ", rownames(coords_3d), "<br>",
        genes[1], ": ", round(gene_expr, 2)
      ),
      hoverinfo = "text"
    ) %>%
      layout(
        scene = list(
          xaxis = xaxis_config,
          yaxis = yaxis_config,
          zaxis = zaxis_config,
          camera = list(
            eye = list(x = 1.5, y = 1.5, z = 1.5)
          ),
          bgcolor = bg_color
        ),
        paper_bgcolor = bg_color,
        plot_bgcolor = bg_color,
        title = list(
          text = paste("Expression of", genes[1]), 
          x = 0.5,
          font = list(color = text_color)
        ),
        font = list(color = text_color)
      )
    
  } else if(length(genes) == 2) {
    # Two genes - blend colors (RGB)
    gene1_expr <- FetchData(seurat_object, vars = genes[1])[[1]]
    gene2_expr <- FetchData(seurat_object, vars = genes[2])[[1]]
    
    # Normalize to 0-1 for color blending
    gene1_norm <- (gene1_expr - min(gene1_expr)) / (max(gene1_expr) - min(gene1_expr))
    gene2_norm <- (gene2_expr - min(gene2_expr)) / (max(gene2_expr) - min(gene2_expr))
    
    # Create RGB colors (gene1=red, gene2=green)
    colors <- rgb(gene1_norm, gene2_norm, 0)
    
    # Order by combined expression
    if(order) {
      combined_expr <- gene1_expr + gene2_expr
      cell_order <- order(combined_expr)
      coords_3d <- coords_3d[cell_order, ]
      colors <- colors[cell_order]
      gene1_expr <- gene1_expr[cell_order]
      gene2_expr <- gene2_expr[cell_order]
    }
    
    # Create axis configs
    xaxis_config <- axis_config
    yaxis_config <- axis_config
    zaxis_config <- axis_config
    
    if (!hide_axes) {
      xaxis_config$title <- list(text = paste0(toupper(reduction), "_1"), font = list(color = text_color))
      yaxis_config$title <- list(text = paste0(toupper(reduction), "_2"), font = list(color = text_color))
      zaxis_config$title <- list(text = paste0(toupper(reduction), "_3"), font = list(color = text_color))
    } else {
      xaxis_config$title <- ""
      yaxis_config$title <- ""
      zaxis_config$title <- ""
    }
    
    # Create plot
    fig <- plot_ly(
      x = coords_3d[, 1],
      y = coords_3d[, 2],
      z = coords_3d[, 3],
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = pt_size,
        color = colors,
        line = list(width = 0)
      ),
      text = ~paste0(
        "Cell: ", rownames(coords_3d), "<br>",
        genes[1], ": ", round(gene1_expr, 2), "<br>",
        genes[2], ": ", round(gene2_expr, 2)
      ),
      hoverinfo = "text"
    ) %>%
      layout(
        scene = list(
          xaxis = xaxis_config,
          yaxis = yaxis_config,
          zaxis = zaxis_config,
          camera = list(
            eye = list(x = 1.5, y = 1.5, z = 1.5)
          ),
          bgcolor = bg_color
        ),
        paper_bgcolor = bg_color,
        plot_bgcolor = bg_color,
        title = list(
          text = paste("Co-expression:", genes[1], "(red) +", genes[2], "(green)"),
          x = 0.5,
          font = list(color = text_color)
        ),
        font = list(color = text_color)
      )
    
  } else if(length(genes) == 3) {
    # Three genes - RGB blend
    gene1_expr <- FetchData(seurat_object, vars = genes[1])[[1]]
    gene2_expr <- FetchData(seurat_object, vars = genes[2])[[1]]
    gene3_expr <- FetchData(seurat_object, vars = genes[3])[[1]]
    
    # Normalize to 0-1
    gene1_norm <- (gene1_expr - min(gene1_expr)) / (max(gene1_expr) - min(gene1_expr))
    gene2_norm <- (gene2_expr - min(gene2_expr)) / (max(gene2_expr) - min(gene2_expr))
    gene3_norm <- (gene3_expr - min(gene3_expr)) / (max(gene3_expr) - min(gene3_expr))
    
    # Create RGB colors
    colors <- rgb(gene1_norm, gene2_norm, gene3_norm)
    
    # Order by combined expression
    if(order) {
      combined_expr <- gene1_expr + gene2_expr + gene3_expr
      cell_order <- order(combined_expr)
      coords_3d <- coords_3d[cell_order, ]
      colors <- colors[cell_order]
      gene1_expr <- gene1_expr[cell_order]
      gene2_expr <- gene2_expr[cell_order]
      gene3_expr <- gene3_expr[cell_order]
    }
    
    # Create axis configs
    xaxis_config <- axis_config
    yaxis_config <- axis_config
    zaxis_config <- axis_config
    
    if (!hide_axes) {
      xaxis_config$title <- list(text = paste0(toupper(reduction), "_1"), font = list(color = text_color))
      yaxis_config$title <- list(text = paste0(toupper(reduction), "_2"), font = list(color = text_color))
      zaxis_config$title <- list(text = paste0(toupper(reduction), "_3"), font = list(color = text_color))
    } else {
      xaxis_config$title <- ""
      yaxis_config$title <- ""
      zaxis_config$title <- ""
    }
    
    # Create plot
    fig <- plot_ly(
      x = coords_3d[, 1],
      y = coords_3d[, 2],
      z = coords_3d[, 3],
      type = "scatter3d",
      mode = "markers",
      marker = list(
        size = pt_size,
        color = colors,
        line = list(width = 0)
      ),
      text = ~paste0(
        "Cell: ", rownames(coords_3d), "<br>",
        genes[1], ": ", round(gene1_expr, 2), "<br>",
        genes[2], ": ", round(gene2_expr, 2), "<br>",
        genes[3], ": ", round(gene3_expr, 2)
      ),
      hoverinfo = "text"
    ) %>%
      layout(
        scene = list(
          xaxis = xaxis_config,
          yaxis = yaxis_config,
          zaxis = zaxis_config,
          camera = list(
            eye = list(x = 1.5, y = 1.5, z = 1.5)
          ),
          bgcolor = bg_color
        ),
        paper_bgcolor = bg_color,
        plot_bgcolor = bg_color,
        title = list(
          text = paste("Co-expression:", genes[1], "(R) +", genes[2], "(G) +", genes[3], "(B)"),
          x = 0.5,
          font = list(color = text_color)
        ),
        font = list(color = text_color)
      )
    
  } else {
    stop("Maximum 3 genes supported for 3D FeaturePlot")
  }
  
  # Enable WebGL for better performance
  fig <- toWebGL(fig)
  
  return(fig)
}

############################## DotPlot Function ##############################

create_dotplot <- function(seurat_object,
                           genes,
                           group_by = NULL,
                           split_by = NULL,
                           color_palette = "default",
                           dot_scale = 1,
                           assay = "RNA",
                           rotate_axis = TRUE,
                           axis_text_size = 12,
                           title_text_size = 14,
                           axis_line_width = 1) {
  # Create DotPlot with customizable parameters
  # Args:
  #   seurat_object: Seurat object
  #   genes: Character vector of gene names
  #   group_by: Column name for grouping (NULL uses Idents)
  #   split_by: Column name for splitting (optional)
  #   color_palette: Color palette name ("default" uses Seurat colors)
  #   dot_scale: Dot size scaling factor (default: 1)
  #   assay: Assay to use (default: "RNA")
  #   rotate_axis: Rotate x-axis labels (default: TRUE)
  #   axis_text_size: Size of axis text
  #   title_text_size: Size of title text
  #   axis_line_width: Width of axis lines
  # Returns:
  #   ggplot object (DotPlot)
  
  # Set assay
  DefaultAssay(seurat_object) <- assay
  
  # Determine which DotPlot variant to use
  use_default_colors <- (color_palette == "default")
  has_split <- !is.null(split_by) && split_by != ""
  use_seurat_clusters <- is.null(group_by) || group_by == "seurat_clusters" || group_by == ""
  
  # Resolve cols argument:
  # When split.by is active, Seurat requires an explicit named color vector
  # (one color per unique split value). Without it, Seurat tries to call
  # brewer.pal internally which fails with < 3 values, or with palette name strings.
  # When split.by is NULL, cols can be a palette name string for gradient coloring.
  if (has_split) {
    split_values <- unique(as.character(seurat_object@meta.data[[split_by]]))
    n <- length(split_values)
    
    if (use_default_colors) {
      # Seurat cannot auto-generate colors with split.by active.
      # Use Set2 as the default qualitative palette.
      # brewer.pal requires n >= 3, so clamp and slice to actual n.
      if (n <= 8) {
        colors <- RColorBrewer::brewer.pal(max(3, n), "Set2")[1:n]
      } else {
        colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n)
      }
    } else {
      colors <- switch(color_palette,
                       "RdYlBu"   = grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(n),
                       "Blues"    = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,  "Blues"))(n),
                       "Reds"     = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,  "Reds"))(n),
                       "Greens"   = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,  "Greens"))(n),
                       "Spectral" = grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n),
                       "PuOr"     = grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "PuOr"))(n),
                       "BrBG"     = grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "BrBG"))(n),
                       # Fallback: viridis
                       viridis::viridis(n)
      )
    }
    
    # Named vector required by Seurat DotPlot when split.by is active
    cols_arg <- stats::setNames(colors, split_values)
    
  } else if (!use_default_colors) {
    # No split, custom palette: pass palette name directly as gradient string
    # Seurat accepts RColorBrewer palette names in non-split mode
    cols_arg <- color_palette
  } else {
    # No split, default colors: let Seurat handle it entirely
    cols_arg <- NULL
  }
  
  # Build DotPlot arguments dynamically to avoid duplicated call blocks
  dotplot_args <- list(
    object    = seurat_object,
    features  = genes,
    dot.scale = dot_scale
  )
  
  if (!use_seurat_clusters) {
    dotplot_args$group.by <- group_by
  }
  
  if (has_split) {
    dotplot_args$split.by <- split_by
  }
  
  if (!is.null(cols_arg)) {
    dotplot_args$cols <- cols_arg
  }
  
  plot <- do.call(DotPlot, dotplot_args)
  
  # Apply theme customizations
  if (rotate_axis) {
    plot <- plot + RotatedAxis()
  }
  
  plot <- plot + theme(
    axis.text      = element_text(size = axis_text_size),
    axis.text.x    = element_text(angle = 45, hjust = 1),
    axis.title     = element_text(size = title_text_size),
    plot.title     = element_text(size = title_text_size),
    legend.text    = element_text(size = axis_text_size),
    axis.line      = element_line(linewidth = axis_line_width),
    axis.ticks     = element_line(linewidth = axis_line_width)
  )
  
  return(plot)
}