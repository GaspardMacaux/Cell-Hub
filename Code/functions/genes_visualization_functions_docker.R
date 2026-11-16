# visualization_workaround_functions.R


############################## Gene input validation ##############################
# Split a gene text field into unique tokens (pure function, no Shiny code).
# Separators: comma, semicolon, space, tab or new line (lists pasted from Excel work as is).
parse_gene_tokens <- function(gene_text) {
  gene_text <- gene_text[!is.na(gene_text)]
  if (length(gene_text) == 0) return(character(0))
  tokens <- unlist(strsplit(paste(gene_text, collapse = " "), "[,;[:space:]]+"))
  unique(tokens[nzchar(tokens)])
}

# Match the tokens of a gene field against an assay (pure function, no Shiny code).
# Exact match first, then case-insensitive match: "gapdh", "GAPDH" and "Gapdh" all return the name stored in the object.
# Returns genes (unique matched names), missing tokens, tokens / matched (same order and length, NA = not found),
# ok (at least one gene found), message (NULL when there is nothing to report) and the assay actually used.
validate_gene_input <- function(gene_text, seurat_object, assay = DefaultAssay(seurat_object)) {
  tokens <- parse_gene_tokens(gene_text)
  if (is.null(assay) || length(assay) != 1 || !assay %in% names(seurat_object@assays)) assay <- DefaultAssay(seurat_object)
  available <- rownames(seurat_object[[assay]])
  idx <- match(tokens, available)
  ci_idx <- match(tolower(tokens), tolower(available))
  idx[is.na(idx)] <- ci_idx[is.na(idx)]
  matched <- available[idx]
  present <- unique(matched[!is.na(matched)])
  missing <- tokens[is.na(matched)]
  corrected <- !is.na(matched) & matched != tokens
  messages <- character(0)
  if (length(tokens) == 0) {
    messages <- "Enter at least one gene name (separated by commas, semicolons or spaces)."
  } else if (length(present) == 0) {
    messages <- sprintf("None of the genes were found in assay '%s': %s", assay, paste(missing, collapse = ", "))
  } else {
    if (length(missing) > 0) messages <- c(messages, sprintf("Genes not found in assay '%s' (skipped): %s", assay, paste(missing, collapse = ", ")))
    if (any(corrected)) messages <- c(messages, sprintf("Gene names matched ignoring case: %s", paste(paste(tokens[corrected], "→", matched[corrected]), collapse = ", ")))
  }
  list(genes = present, missing = missing, tokens = tokens, matched = matched, ok = length(present) > 0,
       message = if (length(messages) > 0) paste(messages, collapse = " | ") else NULL, assay = assay)
}

# Keep one value per validated gene (first occurrence) from values given in token order; names are the object's gene names.
# Used by the co-expression thresholds, whose input fields are indexed by token position.
align_gene_values <- function(check, values) {
  keep <- !is.na(check$matched) & !duplicated(check$matched)
  stats::setNames(values[keep], check$matched[keep])
}
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


create3DFeaturePlot <- function(seurat_object, genes, reduction = "umap3d", pt_size = 3, min_cutoff = NA, max_cutoff = NA, order = TRUE,
                                hide_grid = FALSE, hide_axes = FALSE, dark_mode = FALSE) {
  # 3D FeaturePlot with plotly: one gene coloured by expression (Viridis), two or three genes blended as RGB channels
  # Args:
  #   seurat_object: Seurat object with a reduction of at least 3 dimensions
  #   genes: 1 to 3 gene names present in the default assay
  #   min_cutoff, max_cutoff: as in Seurat FeaturePlot, a number or "qN" (N-th percentile of each gene); NA, "NA" or "" = none
  #   order: draw cells with the highest expression last (on top)
  # Returns:
  #   plotly object (WebGL)
  if (!reduction %in% names(seurat_object@reductions)) stop(paste("Reduction", reduction, "not found in Seurat object"))
  if (length(genes) < 1 || length(genes) > 3) stop("3D FeaturePlot supports 1 to 3 genes")
  coords <- Embeddings(seurat_object, reduction = reduction)
  if (ncol(coords) < 3) stop(paste("Reduction", reduction, "must have at least 3 dimensions for 3D plotting"))
  coords <- coords[, 1:3, drop = FALSE]
  expr <- FetchData(seurat_object, vars = genes)[rownames(coords), genes, drop = FALSE]
  # Quantile cutoffs ("q10") are resolved per gene; comparing numbers with the string "q10" would turn the values into text
  resolve_cutoff <- function(cutoff, values) {
    if (is.null(cutoff) || length(cutoff) == 0 || is.na(cutoff[1]) || as.character(cutoff[1]) %in% c("", "NA")) return(NA_real_)
    if (is.character(cutoff) && grepl("^q[0-9.]+$", cutoff[1])) return(unname(stats::quantile(values, as.numeric(sub("^q", "", cutoff[1])) / 100, na.rm = TRUE)))
    return(suppressWarnings(as.numeric(cutoff[1])))
  }
  for (g in genes) {
    lo <- resolve_cutoff(min_cutoff, expr[[g]])
    hi <- resolve_cutoff(max_cutoff, expr[[g]])
    if (!is.na(lo)) expr[[g]] <- pmax(expr[[g]], lo)
    if (!is.na(hi)) expr[[g]] <- pmin(expr[[g]], hi)
  }
  if (isTRUE(order)) {
    cell_order <- order(rowSums(expr))
    coords <- coords[cell_order, , drop = FALSE]
    expr <- expr[cell_order, , drop = FALSE]
  }
  bg_color <- if (dark_mode) "#1a1a1a" else "white"
  grid_color <- if (dark_mode) "#404040" else "#e0e0e0"
  text_color <- if (dark_mode) "#ffffff" else "#000000"
  axis_config <- function(i) list(title = if (hide_axes) "" else list(text = paste0(toupper(reduction), "_", i), font = list(color = text_color)),
                                  showgrid = !hide_grid, showline = !hide_axes, zeroline = !hide_axes, showticklabels = !hide_axes,
                                  gridcolor = grid_color, gridwidth = 1, linecolor = text_color, tickcolor = text_color, tickfont = list(color = text_color))
  hover <- do.call(paste0, c(list("Cell: ", rownames(coords)), lapply(genes, function(g) paste0("<br>", g, ": ", round(expr[[g]], 2)))))
  if (length(genes) == 1) {
    marker <- list(size = pt_size, color = expr[[1]], colorscale = "Viridis", showscale = TRUE, line = list(width = 0),
                   colorbar = list(title = list(text = genes[1], font = list(color = text_color)), tickfont = list(color = text_color)))
    title <- paste("Expression of", genes[1])
  } else {
    # Each gene scaled to 0-1 on its own range; a constant gene (e.g. not detected) gives 0 instead of NaN, which rgb() rejects
    scale01 <- function(v) {
      r <- range(v, na.rm = TRUE)
      if (!all(is.finite(r)) || diff(r) == 0) return(rep(0, length(v)))
      return((v - r[1]) / diff(r))
    }
    channels <- lapply(expr, scale01)
    marker <- list(size = pt_size, color = grDevices::rgb(channels[[1]], channels[[2]], if (length(genes) == 3) channels[[3]] else 0), line = list(width = 0))
    title <- if (length(genes) == 2) paste("Co-expression:", genes[1], "(red) +", genes[2], "(green)") else paste("Co-expression:", genes[1], "(R) +", genes[2], "(G) +", genes[3], "(B)")
  }
  fig <- plotly::plot_ly(x = coords[, 1], y = coords[, 2], z = coords[, 3], type = "scatter3d", mode = "markers", marker = marker, text = hover, hoverinfo = "text")
  fig <- plotly::layout(fig, scene = list(xaxis = axis_config(1), yaxis = axis_config(2), zaxis = axis_config(3), camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5)), bgcolor = bg_color),
                        paper_bgcolor = bg_color, plot_bgcolor = bg_color, title = list(text = title, x = 0.5, font = list(color = text_color)), font = list(color = text_color))
  # WebGL rendering for large numbers of cells
  return(plotly::toWebGL(fig))
}

############################## DotPlot Function ##############################

create_dotplot <- function(seurat_object,
                           genes,
                           group_by = NULL,
                           split_by = NULL,
                           color_palette = "default",
                           dot_scale = 7,
                           assay = "RNA",
                           scale_expression = TRUE,
                           rotate_axis = TRUE,
                           axis_text_size = 12,
                           title_text_size = 14,
                           axis_line_width = 1) {
  DefaultAssay(seurat_object) <- assay
  # Rule #1: merged Seurat v5 objects must have joined layers before DotPlot
  if (inherits(seurat_object[[assay]], "Assay5") && length(Layers(seurat_object[[assay]], search = "data")) > 1) {
    message("[create_dotplot] Multiple data layers detected — calling JoinLayers()")
    seurat_object[[assay]] <- JoinLayers(seurat_object[[assay]])
  }
  # Keep only genes present in the assay and report the missing ones
  genes <- unique(genes)
  missing_genes <- setdiff(genes, rownames(seurat_object[[assay]]))
  if (length(missing_genes) > 0) message("[create_dotplot] Genes not found in assay ", assay, ": ", paste(missing_genes, collapse = ", "))
  genes <- setdiff(genes, missing_genes)
  if (length(genes) == 0) stop("None of the requested genes were found in assay '", assay, "'")
  # Resolve the grouping vector (current Idents or a metadata column), preserving factor level order
  use_seurat_idents <- is.null(group_by) || !nzchar(group_by) || group_by == "seurat_clusters"
  if (!use_seurat_idents && !group_by %in% colnames(seurat_object@meta.data)) stop("Grouping column not found in metadata: ", group_by)
  group_vec <- droplevels(as.factor(if (use_seurat_idents) Idents(seurat_object) else seurat_object@meta.data[[group_by]]))
  has_split <- !is.null(split_by) && nzchar(split_by) && split_by != "None"
  if (has_split && !split_by %in% colnames(seurat_object@meta.data)) stop("Split column not found in metadata: ", split_by)
  if (has_split) {
    # Build one identity per group x split combination so each dot is computed on the right cells
    split_vec <- droplevels(as.factor(seurat_object@meta.data[[split_by]]))
    combo_map <- data.frame(group = rep(levels(group_vec), each = nlevels(split_vec)),
                            split = rep(levels(split_vec), times = nlevels(group_vec)), stringsAsFactors = FALSE)
    combo_map$id <- paste(combo_map$group, combo_map$split, sep = " | ")
    cell_combo <- paste(group_vec, split_vec, sep = " | ")
    combo_map <- combo_map[combo_map$id %in% cell_combo, ]
    Idents(seurat_object) <- factor(cell_combo, levels = combo_map$id)
  } else {
    Idents(seurat_object) <- group_vec
  }
  # Extract averaged data from Seurat (z-score per gene across identities when scale = TRUE)
  dot_data <- DotPlot(seurat_object, features = genes, assay = assay, scale = scale_expression)$data
  n_ids <- nlevels(droplevels(dot_data$id))
  # Seurat silently disables scaling when only one identity is present
  effective_scale <- scale_expression && n_ids > 1
  if (effective_scale && n_ids < 5) {
    message(sprintf("[create_dotplot] Only %d groups — z-scores bounded to +/-%.2f, small differences are visually exaggerated",
                    n_ids, (n_ids - 1) / sqrt(n_ids)))
  }
  if (has_split) {
    dot_data$y_axis <- factor(combo_map$group[match(as.character(dot_data$id), combo_map$id)], levels = levels(group_vec))
    dot_data$split_group <- factor(combo_map$split[match(as.character(dot_data$id), combo_map$id)], levels = levels(split_vec))
  } else {
    dot_data$y_axis <- dot_data$id
  }
  # Legend title reflects the scaling mode actually applied
  legend_title <- if (effective_scale) "Avg. expression\n(z-score across groups)" else "Avg. expression\n(log1p, normalized)"
  # Symmetric limits for diverging palettes on scaled data so the palette midpoint equals the gene mean
  diverging_palettes <- c("RdYlBu", "Spectral", "PuOr", "BrBG")
  max_abs <- suppressWarnings(max(abs(dot_data$avg.exp.scaled), na.rm = TRUE))
  color_limits <- if (effective_scale && color_palette %in% diverging_palettes && is.finite(max_abs) && max_abs > 0) c(-max_abs, max_abs) else NULL
  color_scale <- switch(color_palette,
                        "RdYlBu"   = scale_color_distiller(palette = "RdYlBu",   direction = -1, limits = color_limits, name = legend_title),
                        "Spectral" = scale_color_distiller(palette = "Spectral", direction = -1, limits = color_limits, name = legend_title),
                        "PuOr"     = scale_color_distiller(palette = "PuOr",     direction = -1, limits = color_limits, name = legend_title),
                        "BrBG"     = scale_color_distiller(palette = "BrBG",     direction = -1, limits = color_limits, name = legend_title),
                        "Blues"    = scale_color_distiller(palette = "Blues",    direction = 1,  name = legend_title),
                        "Reds"     = scale_color_distiller(palette = "Reds",     direction = 1,  name = legend_title),
                        "Greens"   = scale_color_distiller(palette = "Greens",   direction = 1,  name = legend_title),
                        # default: Seurat's lightgrey-to-blue gradient
                        scale_color_gradient(low = "lightgrey", high = "blue", name = legend_title))
  # Dot radius mapped on a fixed 0-100% range so dot sizes are comparable across figures
  plot <- ggplot(dot_data, aes(x = features.plot, y = y_axis, color = avg.exp.scaled, size = pct.exp)) +
    geom_point() +
    color_scale +
    scale_radius(range = c(0, dot_scale), limits = c(0, 100), name = "% expressing") +
    labs(x = "Features", y = if (has_split) "Group" else "Identity") +
    cowplot::theme_cowplot() +
    theme(axis.text.x  = element_text(angle = if (rotate_axis) 45 else 0, hjust = if (rotate_axis) 1 else 0.5, size = axis_text_size),
          axis.text.y  = element_text(size = axis_text_size),
          axis.title   = element_text(size = title_text_size),
          plot.title   = element_text(size = title_text_size),
          legend.text  = element_text(size = axis_text_size),
          axis.line    = element_line(linewidth = axis_line_width),
          axis.ticks   = element_line(linewidth = axis_line_width))
  # One panel per split level; fixed scales keep coord_flip() compatible
  if (has_split) {
    plot <- plot + facet_wrap(~ split_group, nrow = 1) +
      theme(strip.text = element_text(size = axis_text_size, face = "bold"), strip.background = element_blank(),
            panel.spacing = unit(0.5, "lines"))
  }
  message(paste("DotPlot created - palette:", color_palette, "| dot.scale:", dot_scale,
                "| scaled:", effective_scale, "| split.by:", if (has_split) split_by else "none"))
  return(plot)
}