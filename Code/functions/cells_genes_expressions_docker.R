# cells_genes_expressions.R
# Modular functions for gene expression, co-expression, and cluster composition analysis
# Works with both single and integrated datasets


# At the top of your file to avoid namespace conflicts
filter <- dplyr::filter
select <- dplyr::select
mutate <- dplyr::mutate
case_when <- dplyr::case_when

# Clean gene names for HTML display
clean_gene_names_for_html <- function(genes) {
  # Clean gene names for HTML display
  cleaned <- gsub("[^a-zA-Z0-9_.-]", "_", genes)
  return(cleaned)
}

# Function to create visual percentage bar for expression
create_expression_percentage_bar <- function(percentage) {
  # Color gradient based on expression percentage
  if (percentage > 75) {
    color <- "#28a745"  # Green for high expression
  } else if (percentage > 50) {
    color <- "#17a2b8"  # Cyan for moderate-high
  } else if (percentage > 25) {
    color <- "#ffc107"  # Yellow for moderate
  } else if (percentage > 10) {
    color <- "#fd7e14"  # Orange for low
  } else {
    color <- "#dc3545"  # Red for very low
  }
  
  bar_html <- paste0(
    '<div style="width: 100%; height: 20px; border: 1px solid #ddd; ',
    'border-radius: 4px; overflow: hidden; background-color: #f8f9fa;">',
    '<div style="width: ', percentage, '%; height: 100%; ',
    'background: linear-gradient(90deg, ', color, ' 0%, ', 
    adjustcolor(color, alpha.f = 0.7), ' 100%); ',
    'display: flex; align-items: center; justify-content: center; ',
    'transition: width 0.3s ease;" ',
    'title="', percentage, '% of cells express this gene">',
    '<span style="font-size: 11px; color: white; font-weight: bold; ',
    'text-shadow: 1px 1px 2px rgba(0,0,0,0.5);">',
    if (percentage > 15) paste0(percentage, '%') else '',
    '</span>',
    '</div></div>'
  )
  
  return(bar_html)
}

# Function to create size visualization bar for cluster composition
create_size_bar <- function(percentage) {
  # Color gradient based on size
  if (percentage > 20) {
    color <- "#2E86AB"  # Dark blue for large clusters
  } else if (percentage > 10) {
    color <- "#A23B72"  # Purple for medium clusters
  } else if (percentage > 5) {
    color <- "#F18F01"  # Orange for small clusters
  } else {
    color <- "#C73E1D"  # Red for very small clusters
  }
  
  bar_html <- paste0(
    '<div style="width: 100%; height: 20px; border: 1px solid #ccc; ',
    'display: flex; align-items: center; background-color: #f5f5f5;">',
    '<div style="width: ', min(percentage * 4, 100), '%; height: 18px; ',
    'background-color: ', color, '; margin: 1px; ',
    'display: flex; align-items: center; justify-content: center;" ',
    'title="', percentage, '% of total cells">',
    if (percentage > 5) paste0('<span style="font-size: 10px; color: white; font-weight: bold;">', 
                               percentage, '%</span>') else '',
    '</div></div>'
  )
  
  return(bar_html)
}

# Helper function to create co-expression visualization bar
create_coexpression_bar <- function(pct_both, pct_only1, pct_only2, pct_neither) {
  bar_html <- paste0(
    '<div style="width: 150px; height: 20px; display: flex; border: 1px solid #ddd; border-radius: 3px; overflow: hidden;">',
    '<div style="width: ', pct_neither, '%; background: #e0e0e0;" title="Neither: ', pct_neither, '%"></div>',
    '<div style="width: ', pct_only2, '%; background: #87CEEB;" title="Only Gene 2: ', pct_only2, '%"></div>',
    '<div style="width: ', pct_only1, '%; background: #FFA07A;" title="Only Gene 1: ', pct_only1, '%"></div>',
    '<div style="width: ', pct_both, '%; background: #32CD32;" title="Both: ', pct_both, '%"></div>',
    '</div>'
  )
  return(bar_html)
}

# Main function to analyze gene expression with enhanced metrics
analyze_gene_expression <- function(seurat_obj,
                                    selected_genes,
                                    assay_name = "RNA",
                                    expression_threshold = 0.1,
                                    is_integrated = FALSE) {
  tryCatch({
    DefaultAssay(seurat_obj) <- assay_name
    # Split data layers (merged Seurat v5 objects) are joined so LayerData returns every cell; v3 assays (integrated)
    # and SCT assays have a single data matrix and no JoinLayers method
    if (inherits(seurat_obj[[assay_name]], "Assay5") && length(SeuratObject::Layers(seurat_obj[[assay_name]], search = "data")) > 1) {
      seurat_obj <- JoinLayers(seurat_obj, assay = assay_name)
    }
    gene_data <- LayerData(seurat_obj, assay = assay_name, layer = "data")
    available_genes <- rownames(gene_data)
    selected_genes <- intersect(selected_genes, available_genes)
    if (length(selected_genes) == 0) {
      stop("No selected genes are present in the dataset")
    }
    cluster_info <- data.frame(
      cell = colnames(seurat_obj),
      cluster = as.character(Idents(seurat_obj)),
      stringsAsFactors = FALSE
    )
    # Per-dataset breakdown needs a 'dataset' column; without it the per-cluster view is used
    if (is_integrated && !"dataset" %in% colnames(seurat_obj@meta.data)) {
      warning("No 'dataset' column in metadata: expression summarized per cluster only")
      is_integrated <- FALSE
    }
    if (is_integrated) {
      cluster_info$dataset <- as.character(seurat_obj@meta.data$dataset)
      cluster_counts_total <- cluster_info %>%
        group_by(cluster, dataset) %>%
        summarise(n = n(), .groups = 'drop')
    } else {
      cluster_counts_total <- table(cluster_info$cluster)
    }
    expression_summary_list <- lapply(selected_genes, function(gene) {
      gene_expression <- as.numeric(gene_data[gene, ])
      expressed_indices <- gene_expression > expression_threshold
      if (sum(expressed_indices) == 0) {
        return(NULL)
      }
      expressing_cells <- cluster_info[expressed_indices, ]
      if (is_integrated) {
        expression_summary <- expressing_cells %>%
          group_by(cluster, dataset) %>%
          summarise(Cells_Expressed = n(), .groups = 'drop')
        mean_expr_by_group <- sapply(1:nrow(expression_summary), function(i) {
          all_cells_in_group <- which(
            cluster_info$cluster == expression_summary$cluster[i] &
              cluster_info$dataset == expression_summary$dataset[i]
          )
          mean(gene_expression[all_cells_in_group])
        })
        expression_summary$Mean_Expression <- round(mean_expr_by_group, 3)
        expression_summary <- merge(
          expression_summary,
          cluster_counts_total,
          by = c("cluster", "dataset"),
          all.x = TRUE
        )
        names(expression_summary)[which(names(expression_summary) == "n")] <- "Total_Cells_in_Cluster"
        dataset_totals <- cluster_info %>%
          group_by(dataset) %>%
          summarise(Total_Cells_Dataset = n(), .groups = 'drop')
        expression_summary <- merge(expression_summary, dataset_totals, by = "dataset")
        expression_summary$Percentage_of_Cluster <- round(
          (expression_summary$Cells_Expressed / expression_summary$Total_Cells_in_Cluster) * 100, 2
        )
        expression_summary$Percentage_of_Dataset <- round(
          (expression_summary$Cells_Expressed / expression_summary$Total_Cells_Dataset) * 100, 2
        )
      } else {
        cluster_counts_expressed <- table(expressing_cells$cluster)
        unique_clusters <- names(cluster_counts_expressed)
        mean_expression_by_cluster <- sapply(unique_clusters, function(cluster_name) {
          all_cells_in_cluster <- which(cluster_info$cluster == cluster_name)
          mean(gene_expression[all_cells_in_cluster])
        })
        expression_summary <- data.frame(
          cluster = unique_clusters,
          Cells_Expressed = as.numeric(cluster_counts_expressed),
          Total_Cells_in_Cluster = as.numeric(cluster_counts_total[unique_clusters]),
          Mean_Expression = round(mean_expression_by_cluster, 3),
          stringsAsFactors = FALSE
        )
        expression_summary$Percentage_of_Cluster <- round(
          (expression_summary$Cells_Expressed / expression_summary$Total_Cells_in_Cluster) * 100, 2
        )
      }
      expression_summary$Gene <- gene
      expression_summary$Expression_Bar <- sapply(
        expression_summary$Percentage_of_Cluster,
        create_expression_percentage_bar
      )
      expression_summary$Cell_Summary <- paste0(
        expression_summary$Cells_Expressed, "/",
        expression_summary$Total_Cells_in_Cluster, " cells"
      )
      return(expression_summary)
    })
    valid_summaries <- Filter(Negate(is.null), expression_summary_list)
    if (length(valid_summaries) == 0) {
      stop("No cells express the selected genes above the threshold")
    }
    expression_df <- do.call(rbind, valid_summaries)
    if (is_integrated) {
      col_order <- c("Gene", "cluster", "dataset", "Cell_Summary",
                     "Cells_Expressed", "Total_Cells_in_Cluster",
                     "Percentage_of_Cluster", "Percentage_of_Dataset",
                     "Expression_Bar", "Mean_Expression")
    } else {
      col_order <- c("Gene", "cluster", "Cell_Summary",
                     "Cells_Expressed", "Total_Cells_in_Cluster",
                     "Percentage_of_Cluster", "Expression_Bar",
                     "Mean_Expression")
    }
    expression_df <- expression_df[, col_order]
    expression_df <- expression_df[order(expression_df$Gene, -expression_df$Percentage_of_Cluster), ]
    return(list(
      data = expression_df,
      threshold = expression_threshold,
      is_integrated = is_integrated
    ))
  }, error = function(e) {
    stop(paste("Error in analyze_gene_expression:", e$message))
  })
}

# Build the DT table shared by differential-expression result views.
# filter = "top" renders a range slider on numeric columns, a dropdown on factors and a free-text box on
# characters, so numeric columns (p-values included) must reach this function unformatted;
# p-values are formatted client-side with formatSignif, which keeps numeric filtering and sorting intact.
build_de_datatable <- function(df, pvalue_cols = c("P-value", "Adj. P-value"), order_col = NULL, order_dir = "desc",
                               page_length = 10, length_menu = c(10, 25, 50, 100, 200, 500), gene_col_width = "150px") {
  stopifnot(is.data.frame(df), ncol(df) >= 1)
  # Client-side significant-digit formatting is only meaningful on numeric p-value columns
  pvalue_cols <- pvalue_cols[pvalue_cols %in% colnames(df)]
  pvalue_cols <- pvalue_cols[vapply(pvalue_cols, function(col) is.numeric(df[[col]]), logical(1))]
  opts <- list(pageLength = page_length, lengthMenu = length_menu, dom = "Blfrtip", scrollX = TRUE,
               columnDefs = list(list(className = "dt-center", targets = seq_len(ncol(df) - 1)),
                                 list(width = gene_col_width, targets = 0)))
  # DataTables column indices are 0-based
  if (!is.null(order_col) && order_col %in% colnames(df)) opts$order <- list(list(match(order_col, colnames(df)) - 1, order_dir))
  dt <- DT::datatable(df, escape = FALSE, filter = "top", rownames = FALSE, class = "cell-border stripe", options = opts)
  if (length(pvalue_cols) > 0) dt <- DT::formatSignif(dt, columns = pvalue_cols, digits = 3)
  dt
}

# Enhanced function to analyze gene co-expression with better visualizations
analyze_gene_coexpression <- function(seurat_obj, genes, assay_name = "RNA",
                                      expression_thresholds = NULL, is_integrated = FALSE, threshold_layer = "data") {
  # Compute per-cell co-expression statistics for 2+ genes.
  # Every boolean comparison operates on the SAME cell barcodes — verified by
  # explicit assertions so silent index misalignment cannot happen.
  # Args:
  #   seurat_obj: Seurat object
  #   genes: Character vector or comma-separated string of gene names (>= 2)
  #   assay_name: Assay to pull expression from (default: "RNA")
  #   expression_thresholds: Single value or named vector with threshold per gene
  #   is_integrated: If TRUE, also break results down by dataset
  #   threshold_layer: "counts" (positive = at least `threshold` UMIs) or "data" (positive = log-normalized value above `threshold`)
  # Returns:
  #   List with: data (data.frame), genes_analyzed, thresholds_used, threshold_layer, is_integrated
  tryCatch({
    # Parse genes from text input
    if (is.character(genes) && length(genes) == 1) {
      genes <- unique(trimws(unlist(strsplit(genes, ","))))
    }
    genes <- genes[genes != ""]
    if (length(genes) < 2) stop("At least 2 genes are required for co-expression analysis")
    # --- Threshold setup (default: any detection, i.e. >= 1 UMI or > 0) ---
    default_threshold <- if (threshold_layer == "counts") 1 else 0
    if (is.null(expression_thresholds)) {
      expression_thresholds <- setNames(rep(default_threshold, length(genes)), genes)
    } else if (length(expression_thresholds) == 1) {
      expression_thresholds <- setNames(rep(unname(expression_thresholds), length(genes)), genes)
    } else if (is.null(names(expression_thresholds))) {
      if (length(expression_thresholds) != length(genes)) {
        stop("Number of thresholds must match number of genes or be a single value")
      }
      names(expression_thresholds) <- genes
    }
    for (gene in setdiff(genes, names(expression_thresholds))) expression_thresholds[gene] <- default_threshold
    # --- Gene availability ---
    available_genes <- rownames(seurat_obj[[assay_name]])
    missing_genes <- setdiff(genes, available_genes)
    if (length(missing_genes) > 0) {
      warning(paste("Genes not found:", paste(missing_genes, collapse = ", ")))
      genes <- intersect(genes, available_genes)
    }
    expression_thresholds <- expression_thresholds[genes]
    if (length(genes) < 2) stop("Less than 2 valid genes remaining after filtering")
    # --- Per-cell values (raw counts or normalized data), one row per cell in object order ---
    expression_data <- fetch_threshold_values(seurat_obj, genes, assay_name, threshold_layer)
    expected_cells <- colnames(seurat_obj)
    message(sprintf("[co-expression] %d cells x %d genes, thresholds on %s: %s", nrow(expression_data), length(genes), threshold_layer,
                    paste(format_threshold_label(genes, expression_thresholds, threshold_layer), collapse = ", ")))
    # --- Build per-cell binary matrix (gene x cell) ---
    gene_binary <- matrix(FALSE, nrow = length(genes), ncol = length(expected_cells), dimnames = list(genes, expected_cells))
    for (gene in genes) gene_binary[gene, ] <- call_positive(expression_data[, gene], expression_thresholds[gene], threshold_layer)
    stopifnot("[co-expression] Binary matrix colnames diverged from cell barcodes" = identical(colnames(gene_binary), expected_cells))
    # --- Metadata ---
    metadata <- data.frame(cell = expected_cells, cluster = as.character(Idents(seurat_obj)), stringsAsFactors = FALSE)
    if (is_integrated && "dataset" %in% colnames(seurat_obj@meta.data)) metadata$dataset <- seurat_obj@meta.data$dataset
    threshold_text <- format_threshold_label(genes, expression_thresholds[genes], threshold_layer)
    # --- Inner function: analyze one group of cells ---
    analyze_group <- function(group_cells, group_name, dataset_name = NULL) {
      missing_from_matrix <- setdiff(group_cells, colnames(gene_binary))
      if (length(missing_from_matrix) > 0) {
        stop(sprintf("[co-expression] %d group cells not found in binary matrix (group: %s)", length(missing_from_matrix), group_name))
      }
      group_binary_sub <- gene_binary[, group_cells, drop = FALSE]
      n_total <- length(group_cells)
      results <- data.frame(Group = group_name, N_Cells_Total = n_total, stringsAsFactors = FALSE)
      if (!is.null(dataset_name)) results$Dataset <- dataset_name
      # Per-gene positive counts
      for (gene in genes) {
        col_name <- paste0("N_", gene, "_Positive")
        results[[col_name]] <- sum(group_binary_sub[gene, ])
        pct_value <- round(100 * results[[col_name]] / n_total, 1)
        results[[paste0("Pct_", gene, "_Positive")]] <- pct_value
        results[[paste0(gene, "_Threshold")]] <- expression_thresholds[gene]
        results[[paste0(gene, "_Bar")]] <- create_expression_percentage_bar(pct_value)
      }
      # Co-expression patterns — boolean ops on the SAME cell columns
      if (length(genes) == 2) {
        vec1 <- group_binary_sub[genes[1], ]
        vec2 <- group_binary_sub[genes[2], ]
        stopifnot("[co-expression] Boolean vectors cell order mismatch" = identical(names(vec1), names(vec2)))
        results$N_Both_Positive <- sum(vec1 & vec2)
        results$N_Only_Gene1    <- sum(vec1 & !vec2)
        results$N_Only_Gene2    <- sum(!vec1 & vec2)
        results$N_Neither       <- sum(!vec1 & !vec2)
        stopifnot("[co-expression] Partition does not sum to total cell count" =
                    results$N_Both_Positive + results$N_Only_Gene1 + results$N_Only_Gene2 + results$N_Neither == n_total)
        results$Pct_Both       <- round(100 * results$N_Both_Positive / n_total, 1)
        results$Pct_Only_Gene1 <- round(100 * results$N_Only_Gene1 / n_total, 1)
        results$Pct_Only_Gene2 <- round(100 * results$N_Only_Gene2 / n_total, 1)
        results$Pct_Neither    <- round(100 * results$N_Neither / n_total, 1)
        results$Coexpression_Summary <- paste0(
          "Both (", threshold_text[1], " & ", threshold_text[2], "): ", results$N_Both_Positive, " (", results$Pct_Both, "%) | ",
          "Only ", threshold_text[1], ": ", results$N_Only_Gene1, " (", results$Pct_Only_Gene1, "%) | ",
          "Only ", threshold_text[2], ": ", results$N_Only_Gene2, " (", results$Pct_Only_Gene2, "%)")
      } else {
        # 3+ genes
        n_genes_expressed <- colSums(group_binary_sub)
        results$N_All_Positive  <- sum(n_genes_expressed == length(genes))
        results$N_Any_Positive  <- sum(n_genes_expressed > 0)
        results$N_None_Positive <- sum(n_genes_expressed == 0)
        stopifnot("[co-expression] Any + None does not sum to total" = (results$N_Any_Positive + results$N_None_Positive) == n_total)
        results$Pct_All  <- round(100 * results$N_All_Positive / n_total, 1)
        results$Pct_Any  <- round(100 * results$N_Any_Positive / n_total, 1)
        results$Pct_None <- round(100 * results$N_None_Positive / n_total, 1)
        for (i in 0:length(genes)) results[[paste0("N_Expressing_", i, "_Genes")]] <- sum(n_genes_expressed == i)
        results$Coexpression_Summary <- paste0("All genes (", paste(threshold_text, collapse = ", "), "): ", results$N_All_Positive,
                                               " (", results$Pct_All, "%) | Any gene: ", results$N_Any_Positive, " (", results$Pct_Any, "%)")
      }
      return(results)
    }
    # --- Run analysis per group ---
    if (is_integrated) {
      all_results <- metadata %>% group_by(cluster, dataset) %>% group_map(~ analyze_group(.x$cell, .y$cluster, .y$dataset)) %>% bind_rows()
      dataset_overall <- metadata %>% group_by(dataset) %>% group_map(~ analyze_group(.x$cell, paste("Overall", .y$dataset), .y$dataset)) %>% bind_rows()
      all_results <- rbind(all_results, dataset_overall)
    } else {
      cluster_results <- lapply(unique(metadata$cluster), function(cl) analyze_group(metadata$cell[metadata$cluster == cl], cl))
      all_results <- rbind(analyze_group(metadata$cell, "Overall"), do.call(rbind, cluster_results))
    }
    # Visual bars for 2-gene case
    if (length(genes) == 2) {
      all_results$Coexpression_Visual <- mapply(function(both, only1, only2, neither) create_coexpression_bar(both, only1, only2, neither),
                                                all_results$Pct_Both, all_results$Pct_Only_Gene1, all_results$Pct_Only_Gene2, all_results$Pct_Neither)
    }
    return(list(data = all_results, genes_analyzed = genes, thresholds_used = expression_thresholds,
                threshold_layer = threshold_layer, is_integrated = is_integrated))
  }, error = function(e) {
    stop(paste("Error in analyze_gene_coexpression:", e$message))
  })
}

render_coexpression_table <- function(coexpression_data, table_id) {
  df          <- coexpression_data$data
  is_integrated <- coexpression_data$is_integrated
  genes       <- coexpression_data$genes_analyzed
  
  if (length(genes) == 2) {
    # Separate Overall rows, sort clusters by Pct_Both desc, append Overall at bottom
    is_overall   <- grepl("^Overall", df$Group)
    cluster_rows <- df[!is_overall, ]
    overall_rows <- df[is_overall, ]
    cluster_rows <- cluster_rows[order(-cluster_rows$Pct_Both), ]
    display_df   <- rbind(cluster_rows, overall_rows)
    
    # Keep only the 5-6 meaningful columns
    keep_cols <- c("Group", "N_Cells_Total", "Pct_Both", "Pct_Only_Gene1", "Pct_Only_Gene2")
    if (is_integrated && "Dataset" %in% names(display_df)) {
      keep_cols <- c("Group", "Dataset", "N_Cells_Total", "Pct_Both", "Pct_Only_Gene1", "Pct_Only_Gene2")
    }
    display_df <- display_df[, intersect(keep_cols, names(display_df)), drop = FALSE]
    
    col_names <- if (is_integrated && "Dataset" %in% names(display_df)) {
      c("Cluster", "Dataset", "Total nuclei",
        "Both (%)",
        paste0("Only ", genes[1], " (%)"),
        paste0("Only ", genes[2], " (%)"))
    } else {
      c("Cluster", "Total nuclei",
        "Both (%)",
        paste0("Only ", genes[1], " (%)"),
        paste0("Only ", genes[2], " (%)"))
    }
    
    overall_group_values <- display_df$Group[grepl("^Overall", display_df$Group)]
    
    dt <- datatable(
      display_df,
      options = list(
        pageLength = 30,
        scrollX    = TRUE,
        ordering   = FALSE,   # order set in R, not client-side
        columnDefs = list(list(className = 'dt-center', targets = '_all')),
        dom        = 'Bfrtip',
        buttons    = c('copy', 'csv', 'excel')
      ),
      rownames = FALSE,
      caption  = htmltools::tags$caption(
        style = "caption-side: top; text-align: left; font-weight: bold; font-size: 14px;",
        paste0("Co-expression: ", genes[1], " & ", genes[2],
               " — % of nuclei per cluster expressing each gene")
      ),
      escape   = FALSE,
      colnames = col_names
    ) %>%
      formatStyle(
        'Group',
        target          = 'row',
        backgroundColor = styleEqual(overall_group_values, rep('#fff9c4', length(overall_group_values))),
        fontWeight      = styleEqual(overall_group_values, rep('bold',    length(overall_group_values)))
      ) %>%
      formatStyle(
        'Pct_Both',
        background         = styleColorBar(c(0, 100), '#a8d5a2'),
        backgroundSize     = '98% 70%',
        backgroundRepeat   = 'no-repeat',
        backgroundPosition = 'center',
        fontWeight         = 'bold'
      ) %>%
      formatStyle(
        'Pct_Only_Gene1',
        background         = styleColorBar(c(0, 100), '#90caf9'),
        backgroundSize     = '98% 70%',
        backgroundRepeat   = 'no-repeat',
        backgroundPosition = 'center'
      ) %>%
      formatStyle(
        'Pct_Only_Gene2',
        background         = styleColorBar(c(0, 100), '#ffcc80'),
        backgroundSize     = '98% 70%',
        backgroundRepeat   = 'no-repeat',
        backgroundPosition = 'center'
      ) %>%
      formatRound(c('Pct_Both', 'Pct_Only_Gene1', 'Pct_Only_Gene2'), digits = 1)
    
    if (is_integrated && "Dataset" %in% names(display_df)) {
      dt <- dt %>%
        formatStyle(
          'Dataset',
          backgroundColor = styleEqual(
            unique(display_df$Dataset),
            rainbow(length(unique(display_df$Dataset)), alpha = 0.3)
          )
        )
    }
    
  } else {
    # >2 genes: hide N_ columns, keep Pct_ with color bars
    hidden_cols         <- names(df)[grepl("^N_", names(df))]
    hidden_indices      <- which(names(df) %in% hidden_cols) - 1
    overall_group_values <- df$Group[grepl("^Overall", df$Group)]
    
    dt <- datatable(
      df,
      options = list(
        pageLength = 20,
        scrollX    = TRUE,
        columnDefs = list(
          list(visible = FALSE, targets = hidden_indices),
          list(className = 'dt-center', targets = '_all'),
          list(className = 'dt-left',
               targets = which(names(df) == "Coexpression_Summary") - 1)
        ),
        dom     = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      rownames = FALSE,
      caption  = paste("Co-expression Analysis -", paste(genes, collapse = ", ")),
      escape   = FALSE
    ) %>%
      formatStyle(
        'Group',
        target          = 'row',
        backgroundColor = styleEqual(overall_group_values, rep('#fff9c4', length(overall_group_values))),
        fontWeight      = styleEqual(overall_group_values, rep('bold',    length(overall_group_values)))
      )
    
    pct_cols <- names(df)[grepl("^Pct_", names(df))]
    for (col in pct_cols) {
      dt <- dt %>%
        formatStyle(
          col,
          background         = styleColorBar(c(0, 100), '#e8f4fd'),
          backgroundSize     = '98% 88%',
          backgroundRepeat   = 'no-repeat',
          backgroundPosition = 'center'
        )
    }
    
    if (is_integrated && "Dataset" %in% names(df)) {
      dt <- dt %>%
        formatStyle(
          'Dataset',
          backgroundColor = styleEqual(
            unique(df$Dataset),
            rainbow(length(unique(df$Dataset)), alpha = 0.3)
          )
        )
    }
  }
  
  return(dt)
}

## Group: Function to create cluster composition table with support for integrated datasets
create_cluster_composition_table <- function(seurat_obj, is_integrated = FALSE,
                                             metadata_column = "dataset",
                                             cluster_col = NULL) {
  # Resolve cluster column — never rely on Idents() for external objects
  col <- if (!is.null(cluster_col) && cluster_col %in% colnames(seurat_obj@meta.data)) {
    cluster_col
  } else if ("ClusterIdents" %in% colnames(seurat_obj@meta.data)) {
    "ClusterIdents"
  } else if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    "seurat_clusters"
  } else {
    NULL
  }
  cluster_values <- if (!is.null(col)) {
    as.character(seurat_obj@meta.data[[col]])
  } else {
    as.character(Idents(seurat_obj))
  }
  # QC columns absent from imported objects (e.g. no percent.mt) give NA statistics instead of an error
  meta_or_na <- function(column) if (column %in% colnames(seurat_obj@meta.data)) as.numeric(seurat_obj@meta.data[[column]]) else rep(NA_real_, ncol(seurat_obj))
  cluster_data <- data.frame(
    cluster    = cluster_values,
    nFeature   = meta_or_na("nFeature_RNA"),
    nCount     = meta_or_na("nCount_RNA"),
    percent_mt = meta_or_na("percent.mt"),
    stringsAsFactors = FALSE
  )
  if (is_integrated && metadata_column %in% colnames(seurat_obj@meta.data)) {
    cluster_data$grouping <- seurat_obj@meta.data[[metadata_column]]
  } else if (is_integrated) {
    warning(paste("Metadata column", metadata_column, "not found. Proceeding without grouping."))
    is_integrated <- FALSE
  }
  if (is_integrated) {
    cluster_stats <- cluster_data %>%
      group_by(cluster, grouping) %>%
      summarise(
        Cell_Count    = n(),
        Mean_Genes    = round(mean(nFeature, na.rm = TRUE), 1),
        Median_Genes  = round(median(nFeature, na.rm = TRUE), 1),
        Mean_UMI      = round(mean(nCount, na.rm = TRUE), 1),
        Median_UMI    = round(median(nCount, na.rm = TRUE), 1),
        Mean_MT_Percent = round(mean(percent_mt, na.rm = TRUE), 2),
        .groups = 'drop'
      )
    group_totals <- cluster_data %>%
      group_by(grouping) %>%
      summarise(Total_in_Group = n(), .groups = 'drop')
    cluster_stats <- merge(cluster_stats, group_totals, by = "grouping")
    cluster_stats$Percent_of_Group <- round(
      (cluster_stats$Cell_Count / cluster_stats$Total_in_Group) * 100, 1
    )
    colnames(cluster_stats)[colnames(cluster_stats) == "grouping"]        <- metadata_column
    colnames(cluster_stats)[colnames(cluster_stats) == "Total_in_Group"]  <- paste0("Total_in_", metadata_column)
    colnames(cluster_stats)[colnames(cluster_stats) == "Percent_of_Group"] <- paste0("Percent_of_", metadata_column)
  } else {
    cluster_stats <- cluster_data %>%
      group_by(cluster) %>%
      summarise(
        Cell_Count    = n(),
        Mean_Genes    = round(mean(nFeature, na.rm = TRUE), 1),
        Median_Genes  = round(median(nFeature, na.rm = TRUE), 1),
        Mean_UMI      = round(mean(nCount, na.rm = TRUE), 1),
        Median_UMI    = round(median(nCount, na.rm = TRUE), 1),
        Mean_MT_Percent = round(mean(percent_mt, na.rm = TRUE), 2),
        .groups = 'drop'
      )
  }
  total_cells <- nrow(cluster_data)
  cluster_stats$Percent_of_Total <- round((cluster_stats$Cell_Count / total_cells) * 100, 1)
  cluster_stats$Size_Bar <- sapply(cluster_stats$Percent_of_Total, create_size_bar)
  cluster_stats$Quality_Summary <- paste0(
    "Genes: ", cluster_stats$Mean_Genes, " | ",
    "UMI: ", cluster_stats$Mean_UMI, " | ",
    "MT: ", cluster_stats$Mean_MT_Percent, "%"
  )
  # Reorder columns based on dataset type
  if (is_integrated) {
    base_cols <- c("cluster", metadata_column, "Cell_Count", "Percent_of_Total")
    percent_col <- paste0("Percent_of_", metadata_column)
    stat_cols <- c("Size_Bar", "Mean_Genes", "Median_Genes", "Mean_UMI", "Median_UMI", "Mean_MT_Percent", "Quality_Summary")
    cluster_stats <- cluster_stats[, c(base_cols, percent_col, stat_cols)]
  } else {
    cluster_stats <- cluster_stats[, c("cluster", "Cell_Count", "Percent_of_Total", "Size_Bar",
                                       "Mean_Genes", "Median_Genes", "Mean_UMI", "Median_UMI",
                                       "Mean_MT_Percent", "Quality_Summary")]
  }
  return(cluster_stats)
}

# Enhanced render function for expression analysis table
render_expression_table <- function(expression_data, table_id) {
  
  is_integrated <- expression_data$is_integrated
  df <- expression_data$data
  
  # Define column visibility and names
  if (is_integrated) {
    hidden_cols <- c("Cells_Expressed", "Total_Cells_in_Cluster")
    col_names <- c('Gene', 'Cluster', 'Dataset', 'Cells Expressing', 
                   'Cells Expressed', 'Total Cells', '% in Cluster', 
                   '% in Dataset', 'Expression Visual', 'Mean Expression')
  } else {
    hidden_cols <- c("Cells_Expressed", "Total_Cells_in_Cluster")
    col_names <- c('Gene', 'Cluster', 'Cells Expressing',
                   'Cells Expressed', 'Total Cells', '% Expressing', 
                   'Expression Visual', 'Mean Expression')
  }
  
  # Create column indices for hiding
  all_cols <- names(df)
  hidden_indices <- which(all_cols %in% hidden_cols) - 1
  
  dt <- datatable(
    df,
    options = list(
      pageLength = 15,
      scrollX = TRUE,
      columnDefs = list(
        list(visible = FALSE, targets = hidden_indices),
        list(className = 'dt-center', targets = '_all'),
        list(width = '120px', targets = which(all_cols == "Expression_Bar") - 1),
        list(width = '100px', targets = which(all_cols == "Cell_Summary") - 1)
      ),
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel'),
      order = if (is_integrated) {
        list(list(0, 'asc'), list(6, 'desc'))
      } else {
        list(list(0, 'asc'), list(5, 'desc'))
      }
    ),
    rownames = FALSE,
    caption = paste("Gene Expression Analysis (Threshold:", expression_data$threshold, ")"),
    escape = FALSE,
    colnames = col_names
  ) %>%
    formatStyle(
      'Mean_Expression',
      backgroundColor = styleInterval(
        c(0.5, 1, 2, 5), 
        c('#f0f0f0', '#d4edda', '#fff3cd', '#f8d7da', '#d1ecf1')
      )
    )
  
  # Add styling for dataset column if integrated
  if (is_integrated && "dataset" %in% names(df)) {
    dt <- dt %>%
      formatStyle(
        'dataset',
        backgroundColor = styleEqual(
          unique(df$dataset),
          rainbow(length(unique(df$dataset)), alpha = 0.3)
        )
      )
  }
  
  return(dt)
}



# Render function for cluster composition table
render_cluster_composition_table <- function(cluster_data, is_integrated = FALSE) {
  
  # Define columns based on dataset type
  if (is_integrated) {
    col_names <- c('Cluster', 'Dataset', 'Cell Count', '% Total', '% Dataset', 'Size Visual', 
                   'Mean Genes', 'Median Genes', 'Mean UMI', 'Median UMI', 
                   'Mean MT%', 'Quality Summary')
    center_targets <- c(0, 1, 2, 3, 4, 6, 7, 8, 9, 10)
    left_targets <- c(5, 11)
  } else {
    col_names <- c('Cluster', 'Cell Count', '% Total', 'Size Visual', 
                   'Mean Genes', 'Median Genes', 'Mean UMI', 'Median UMI', 
                   'Mean MT%', 'Quality Summary')
    center_targets <- c(0, 1, 2, 4, 5, 6, 7, 8)
    left_targets <- c(3, 9)
  }
  
  dt <- datatable(
    cluster_data,
    options = list(
      pageLength = 15,
      scrollX = TRUE,
      columnDefs = list(
        list(className = 'dt-center', targets = center_targets),
        list(className = 'dt-left', targets = left_targets),
        list(width = '100px', targets = which(names(cluster_data) == "Size_Bar") - 1),
        list(width = '200px', targets = which(names(cluster_data) == "Quality_Summary") - 1)
      ),
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel')
    ),
    rownames = FALSE,
    caption = "Cluster Composition and Quality Metrics",
    escape = FALSE,
    colnames = col_names
  ) %>%
    formatStyle(
      'Cell_Count',
      background = styleColorBar(range(cluster_data$Cell_Count), '#e8f4fd'),
      backgroundSize = '80% 90%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center'
    ) %>%
    formatStyle(
      'Mean_MT_Percent',
      backgroundColor = styleInterval(c(5, 10, 20), c('#90EE90', '#FFFF99', '#FFB347', '#FF6B6B'))
    )
  
  # Add dataset coloring for integrated data
  if (is_integrated && "dataset" %in% names(cluster_data)) {
    dt <- dt %>%
      formatStyle(
        'dataset',
        backgroundColor = styleEqual(
          unique(cluster_data$dataset),
          rainbow(length(unique(cluster_data$dataset)), alpha = 0.3)
        )
      )
  }
  
  return(dt)
}

# Aggregate per-cluster-dataset rows into per-cluster totals by summing N_ columns
aggregate_coexpr_per_cluster <- function(detail_data, genes_analyzed) {
  n_cols <- grep("^N_", colnames(detail_data), value = TRUE)
  agg    <- aggregate(detail_data[, n_cols, drop = FALSE],
                      by = list(Group = detail_data$Group), FUN = sum)
  total  <- agg$N_Cells_Total
  pct    <- function(n_col) round(100 * agg[[n_col]] / total, 1)
  if (length(genes_analyzed) == 2) {
    agg$Pct_Both       <- pct("N_Both_Positive")
    agg$Pct_Only_Gene1 <- pct("N_Only_Gene1")
    agg$Pct_Only_Gene2 <- pct("N_Only_Gene2")
    agg$Pct_Neither    <- pct("N_Neither")
  } else {
    agg$Pct_All  <- pct("N_All_Positive")
    agg$Pct_Any  <- pct("N_Any_Positive")
    agg$Pct_None <- pct("N_None_Positive")
  }
  agg
}


# Co-expression stacked bar plot. Percentages are recomputed from the summed category counts (N_ columns are exact
# partitions of each cluster x dataset group), so pooling datasets or clusters never stacks two rows on the same bar.
# group_by: "cluster" or "dataset" (primary grouping, one facet per level when a breakdown is shown)
# secondary_split: "none", "dataset" (with group_by = "cluster") or "cluster" (with group_by = "dataset")
# plot_type: "bar" (stacked bars) or "heat" (tiles labelled with the percentage)
create_coexpression_plot <- function(coexpr_data, genes_analyzed, group_by = "cluster", secondary_split = "none", plot_type = "bar") {
  n_genes <- length(genes_analyzed)
  if (n_genes < 2) stop("Co-expression plot needs at least 2 genes")
  # Category count columns and labels, from the strongest co-expression pattern to none
  if (n_genes == 2) {
    count_cols <- c("N_Both_Positive", "N_Only_Gene1", "N_Only_Gene2", "N_Neither")
    category_labels <- c(paste("Both", genes_analyzed[1], "&", genes_analyzed[2]), paste("Only", genes_analyzed[1]), paste("Only", genes_analyzed[2]), "Neither")
    fill_scale <- scale_fill_manual(values = setNames(c("#2E8B57", "#4169E1", "#DC143C", "#D3D3D3"), category_labels), name = "Expression pattern")
  } else {
    counts_desc <- rev(seq(0, n_genes))
    count_cols <- paste0("N_Expressing_", counts_desc, "_Genes")
    category_labels <- paste0(counts_desc, " of ", n_genes, " genes")
    fill_scale <- scale_fill_viridis_d(name = "Genes expressed", direction = -1)
  }
  missing_cols <- setdiff(c("Group", "N_Cells_Total", count_cols), colnames(coexpr_data))
  if (length(missing_cols) > 0) stop("Co-expression table lacks columns: ", paste(missing_cols, collapse = ", "))
  # Cluster-level rows only: "Overall" (single dataset) and "Overall <dataset>" (integrated) rows are sums of them
  rows <- coexpr_data[!grepl("^Overall", coexpr_data$Group), , drop = FALSE]
  if (nrow(rows) == 0) stop("No cluster-level rows to plot")
  rows$Group <- as.character(rows$Group)
  has_dataset <- "Dataset" %in% colnames(rows) && !all(is.na(rows$Dataset))
  if (has_dataset) rows$Dataset <- as.character(rows$Dataset)
  # Without a dataset column only the per-cluster view exists
  if (!has_dataset || !group_by %in% c("cluster", "dataset")) { group_by <- "cluster"; secondary_split <- "none" }
  primary_col <- if (group_by == "dataset") "Dataset" else "Group"
  secondary_col <- if (group_by == "dataset" && identical(secondary_split, "cluster")) "Group" else if (group_by == "cluster" && identical(secondary_split, "dataset")) "Dataset" else NULL
  keys <- c(primary_col, secondary_col)
  # Sum counts over the levels not displayed (datasets pooled per cluster, or clusters pooled per dataset)
  agg <- aggregate(rows[, c("N_Cells_Total", count_cols), drop = FALSE], by = rows[, keys, drop = FALSE], FUN = sum)
  long <- do.call(rbind, lapply(seq_along(count_cols), function(i) {
    data.frame(agg[, keys, drop = FALSE], N_Total = agg$N_Cells_Total, Category = category_labels[i], N = agg[[count_cols[i]]], stringsAsFactors = FALSE)
  }))
  long$Percentage <- 100 * long$N / long$N_Total
  long$Category <- factor(long$Category, levels = category_labels)
  # Clusters ordered by their pooled share of the top category (same order as the results table), datasets alphabetically
  cluster_share <- tapply(rows[[count_cols[1]]], rows$Group, sum) / tapply(rows$N_Cells_Total, rows$Group, sum)
  if ("Group" %in% keys) long$Group <- factor(long$Group, levels = names(sort(cluster_share, decreasing = TRUE)))
  if ("Dataset" %in% keys) long$Dataset <- factor(long$Dataset, levels = sort(unique(rows$Dataset)))
  x_col <- if (is.null(secondary_col)) primary_col else secondary_col
  subtitle <- switch(paste(group_by, if (is.null(secondary_col)) "none" else "split"),
                     "cluster none" = if (has_dataset) "Per cluster, datasets pooled" else "Per cluster",
                     "cluster split" = "Per cluster, split by dataset",
                     "dataset none" = "Per dataset, clusters pooled",
                     "dataset split" = "Per dataset, split by cluster")
  plot_title <- paste("Co-expression:", paste(genes_analyzed, collapse = ", "))
  x_title <- if (x_col == "Dataset") "Dataset" else "Cluster"
  if (identical(plot_type, "heat")) {
    # One row per expression pattern (strongest at the top), one column per group, cell = % of the group's cells
    long$Category <- factor(long$Category, levels = rev(category_labels))
    p <- ggplot(long, aes(x = .data[[x_col]], y = Category, fill = Percentage)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.1f", Percentage), color = Percentage > 50), size = 3.2, show.legend = FALSE) +
      scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
      scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 100), name = "% of cells") +
      labs(title = plot_title, subtitle = subtitle, x = x_title, y = NULL) +
      theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), strip.text = element_text(face = "bold"))
  } else {
    p <- ggplot(long, aes(x = .data[[x_col]], y = Percentage, fill = Category)) +
      geom_col(position = "stack", width = 0.8) + fill_scale +
      scale_y_continuous(limits = c(0, 100.5), expand = expansion(mult = c(0, 0.02))) +
      labs(title = plot_title, subtitle = subtitle, x = x_title, y = "Percentage of cells (%)") +
      theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"))
  }
  # Breakdown: one panel per primary level; datasets stacked vertically so long cluster axes stay readable
  if (!is.null(secondary_col)) p <- p + facet_wrap(stats::as.formula(paste("~", primary_col)), scales = "free_x", ncol = if (primary_col == "Dataset") 1 else NULL)
  return(p)
}

# Within-cluster differential expression: expressors vs non-expressors.
# For 1 gene: cells expressing the gene (at or above threshold) vs cells that do not.
# For 2+ genes: cells co-expressing ALL genes vs the rest.
# Runs FindMarkers between the two groups within the selected cluster(s); several clusters are pooled into one group.
# Pure function — no Shiny reactives.
run_coexpression_de <- function(seurat_obj, genes, clusters, assay_name = "RNA",
                                expression_thresholds = NULL, threshold_layer = "data") {
  # --- Parse and validate genes ---
  if (is.character(genes) && length(genes) == 1) genes <- unique(trimws(unlist(strsplit(genes, ","))))
  genes <- genes[nzchar(genes)]
  if (length(genes) < 1) stop("At least 1 gene is required")
  # --- Threshold setup (default: any detection, i.e. >= 1 UMI or > 0) ---
  default_threshold <- if (threshold_layer == "counts") 1 else 0
  if (is.null(expression_thresholds)) {
    expression_thresholds <- setNames(rep(default_threshold, length(genes)), genes)
  } else if (length(expression_thresholds) == 1) {
    expression_thresholds <- setNames(rep(unname(expression_thresholds), length(genes)), genes)
  } else if (is.null(names(expression_thresholds))) {
    if (length(expression_thresholds) != length(genes)) stop("Number of thresholds must match number of genes or be a single value")
    names(expression_thresholds) <- genes
  }
  for (g in setdiff(genes, names(expression_thresholds))) expression_thresholds[g] <- default_threshold
  # --- Subset to selected cluster(s) ---
  DefaultAssay(seurat_obj) <- assay_name
  clusters <- as.character(clusters)
  cluster_cells <- colnames(seurat_obj)[as.character(Idents(seurat_obj)) %in% clusters]
  if (length(cluster_cells) == 0) stop(sprintf("No cells found in cluster(s): %s", paste(clusters, collapse = ", ")))
  obj_sub <- subset(seurat_obj, cells = cluster_cells)
  cluster_label <- paste(clusters, collapse = " + ")
  message(sprintf("[expr_de] Selected cluster(s): %s — %d cells", cluster_label, ncol(obj_sub)))
  # --- JoinLayers only for split Seurat v5 layers (v3 and SCT assays have one matrix per layer and no JoinLayers method) ---
  if (inherits(obj_sub[[assay_name]], "Assay5") && length(SeuratObject::Layers(obj_sub[[assay_name]], search = "data")) > 1) {
    message("[expr_de] Split data layers detected — calling JoinLayers()")
    obj_sub <- JoinLayers(obj_sub, assay = assay_name)
  }
  # --- Gene availability ---
  missing_genes <- setdiff(genes, rownames(obj_sub[[assay_name]]))
  if (length(missing_genes) > 0) stop(paste("Genes not found:", paste(missing_genes, collapse = ", ")))
  # --- Per-cell values (raw counts or normalized data) and positive calls ---
  expr_data <- fetch_threshold_values(obj_sub, genes, assay_name, threshold_layer)
  gene_above <- vapply(genes, function(g) call_positive(expr_data[, g], expression_thresholds[g], threshold_layer), logical(nrow(expr_data)))
  if (!is.matrix(gene_above)) gene_above <- matrix(gene_above, ncol = length(genes), dimnames = list(rownames(expr_data), genes))
  is_positive <- rowSums(gene_above) == length(genes)
  n_pos <- sum(is_positive)
  n_neg <- sum(!is_positive)
  threshold_text <- paste(format_threshold_label(genes, expression_thresholds[genes], threshold_layer), collapse = " & ")
  if (n_pos < 3) stop(sprintf("Only %d cell(s) with %s — FindMarkers requires at least 3 per group", n_pos, threshold_text))
  if (n_neg < 3) stop(sprintf("Only %d cell(s) without %s — FindMarkers requires at least 3 per group", n_neg, threshold_text))
  message(sprintf("[expr_de] %s: %d cells with %s vs %d others", cluster_label, n_pos, threshold_text, n_neg))
  if (length(genes) == 1) {
    breakdown <- list(positive = n_pos, negative = n_neg)
  } else if (length(genes) == 2) {
    breakdown <- list(only_gene1 = sum(gene_above[, genes[1]] & !gene_above[, genes[2]]),
                      only_gene2 = sum(!gene_above[, genes[1]] & gene_above[, genes[2]]),
                      neither    = sum(!gene_above[, genes[1]] & !gene_above[, genes[2]]))
  } else {
    n_expressed <- rowSums(gene_above)
    breakdown <- list(partial = sum(n_expressed > 0 & n_expressed < length(genes)), none = sum(n_expressed == 0))
  }
  # --- Set idents and run FindMarkers ---
  obj_sub$coexpr_group <- factor(ifelse(is_positive, "coexpressor", "non_coexpressor"), levels = c("coexpressor", "non_coexpressor"))
  Idents(obj_sub) <- "coexpr_group"
  markers <- FindMarkers(obj_sub, ident.1 = "coexpressor", ident.2 = "non_coexpressor", assay = assay_name, min.pct = 0.1, logfc.threshold = 0.25)
  markers$gene <- rownames(markers)
  rownames(markers) <- NULL
  markers <- markers[order(markers$p_val_adj, -abs(markers$avg_log2FC)), ]
  return(list(markers = markers, n_coexpressors = n_pos, n_non_coexpressors = n_neg, breakdown = breakdown, clusters = clusters,
              cluster_label = cluster_label, thresholds_used = expression_thresholds, threshold_layer = threshold_layer, genes_used = genes))
}

# Function to create co-expression summary statistics
create_coexpression_summary_stats <- function(data, genes) {
  if (length(genes) == 2) {
    # For 2 genes
    overall_data <- data[grepl("Overall", data$Group), ]
    summary_stats <- data.frame(Dataset = gsub("Overall ", "", overall_data$Group), Total_Cells = overall_data$N_Cells_Total, Both_Genes = paste0(overall_data$N_Both_Positive, " (", overall_data$Pct_Both, "%)"), Only_Gene1 = paste0(overall_data$N_Only_Gene1, " (", overall_data$Pct_Only_Gene1, "%)"), Only_Gene2 = paste0(overall_data$N_Only_Gene2, " (", overall_data$Pct_Only_Gene2, "%)"), Neither = paste0(overall_data$N_Neither, " (", overall_data$Pct_Neither, "%)"))
    names(summary_stats)[3:6] <- c("Both Genes", paste("Only", genes[1]), paste("Only", genes[2]), "Neither")
  } else {
    # For multiple genes
    overall_data <- data[grepl("Overall", data$Group), ]
    summary_stats <- data.frame(Dataset = gsub("Overall ", "", overall_data$Group), Total_Cells = overall_data$N_Cells_Total, All_Genes = paste0(overall_data$N_All_Positive, " (", overall_data$Pct_All, "%)"), Any_Gene = paste0(overall_data$N_Any_Positive, " (", overall_data$Pct_Any, "%)"), No_Genes = paste0(overall_data$N_None_Positive, " (", overall_data$Pct_None, "%)"))
  }
  return(summary_stats)
}

# Helper function for formatted percentage display
format_percentage_compact <- function(n, total, color = "blue") {
  pct <- round(100 * n / total, 1)
  sprintf(
    '<div style="display: inline-block; position: relative; width: 60px;">
    <div style="position: absolute; width: %.1f%%; height: 20px; background-color: %s; opacity: 0.3;"></div>
    <div style="position: relative; text-align: center; line-height: 20px; font-weight: bold;">%.1f%%</div>
  </div>',
    pct, color, pct
  )
}

# Format p-values for display. A p-value of exactly 0 is a floating-point underflow: the true value is below
# the smallest positive normalized double (.Machine$double.xmin, ~2.2e-308), which is not the machine precision (~2.2e-16)
format_pvalue_robust <- function(p_values) {
  vapply(p_values, function(p) {
    if (is.na(p)) return("NA")
    if (p == 0) return(sprintf("< %.2e", .Machine$double.xmin))
    if (p < 0.01) return(sprintf("%.2e", p))
    sprintf("%.4f", p)
  }, character(1))
}

#' Find exclusive biomarkers for a cluster or group of clusters
#'
#' Decision-support ranking: one-vs-rest Wilcoxon rank-sum test (Seurat FindMarkers, presto backend when installed)
#' on all genes passing a liberal detection prefilter, multiple-testing correction over that whole set, then
#' post-hoc detection-rate / fold-change filters and ranking by S = ((pct_target - pct_other) / 100) * log2FC.
#' log2FC is identical to Seurat v5 avg_log2FC on log-normalized data: log2((sum(expm1(x)) + 1) / n_x) - log2((sum(expm1(y)) + 1) / n_y).
#'
#' @param seurat_obj Seurat object with active identities set to the clusters
#' @param target_cluster Character vector of cluster identifier(s)
#' @param min_pct_target Minimum % of target cells with expression > detection_threshold (0-100, default 50)
#' @param max_pct_other Maximum % of other cells with expression > detection_threshold (0-100, default 25)
#' @param min_log2fc Minimum log2 fold change (Seurat convention, default 1.5)
#' @param detection_threshold Expression value above which a gene counts as detected (default 0)
#' @param min_mean_expr_target Minimum mean log-normalized expression in the target (default 0.5)
#' @param statistical_test "wilcox" or "none" (default "wilcox")
#' @param max_pvalue Maximum adjusted p-value (default 0.05, ignored when statistical_test = "none")
#' @param pvalue_adjustment "BH", "fdr", "bonferroni" or "none" (default "BH")
#' @param assay_name Assay to use (default "RNA")
#' @param top_n Number of top markers to return (NULL = all)
#' @param prefilter_min_pct Liberal prefilter (% detected in target) defining the tested gene set; it does not
#'   depend on the strict filters, so the multiple-testing correction is not restricted to post-selected genes (default 10)
#' @param verbose Print progress messages
#' @return Data frame of exclusive biomarker candidates sorted by specificity score
find_exclusive_biomarkers <- function(seurat_obj, target_cluster, min_pct_target = 50, max_pct_other = 25, min_log2fc = 1.5,
                                      detection_threshold = 0, min_mean_expr_target = 0.5, statistical_test = "wilcox",
                                      max_pvalue = 0.05, pvalue_adjustment = "BH", assay_name = "RNA", top_n = NULL,
                                      prefilter_min_pct = 10, verbose = TRUE) {
  current_idents <- as.character(Idents(seurat_obj))
  target_cluster <- as.character(target_cluster)
  missing_clusters <- setdiff(target_cluster, unique(current_idents))
  if (length(missing_clusters) > 0) {
    stop(paste0("Target cluster(s) not found in active identities: ", paste(missing_clusters, collapse = ", "),
                "\nAvailable clusters: ", paste(unique(current_idents), collapse = ", ")))
  }
  if (min_pct_target < 0 || min_pct_target > 100) stop("min_pct_target must be between 0 and 100")
  if (max_pct_other < 0 || max_pct_other > 100) stop("max_pct_other must be between 0 and 100")
  if (!statistical_test %in% c("wilcox", "none")) stop("statistical_test must be 'wilcox' or 'none'")
  if (!pvalue_adjustment %in% c("bonferroni", "BH", "fdr", "none")) stop("pvalue_adjustment must be 'bonferroni', 'BH', 'fdr', or 'none'")
  DefaultAssay(seurat_obj) <- assay_name
  # Split Seurat v5 layers must be joined before reading the data layer (LayerData fails on multiple layers)
  if (inherits(seurat_obj[[assay_name]], "Assay5")) seurat_obj <- JoinLayers(seurat_obj, assay = assay_name)
  gene_data <- LayerData(seurat_obj, assay = assay_name, layer = "data")
  gene_data <- gene_data[!is.na(rownames(gene_data)), , drop = FALSE]
  is_target <- current_idents %in% target_cluster
  n_target <- sum(is_target)
  n_other <- sum(!is_target)
  if (n_target < 3 || n_other < 3) stop("Target and other groups need at least 3 cells each")
  if (verbose) message(sprintf("=== Exclusive biomarkers: %d target vs %d other cells, %d genes ===", n_target, n_other, nrow(gene_data)))
  # Liberal prefilter on target detection only: defines the tested gene set independently of the strict filters
  pct_target_all <- 100 * Matrix::rowSums(gene_data[, is_target, drop = FALSE] > detection_threshold) / n_target
  genes_to_test <- rownames(gene_data)[pct_target_all >= min(prefilter_min_pct, min_pct_target)]
  if (length(genes_to_test) == 0) {
    warning("No gene detected in enough target cells. Try lowering min_pct_target.")
    return(data.frame())
  }
  if (verbose) message(sprintf("Genes passing the %s%% detection prefilter: %d", min(prefilter_min_pct, min_pct_target), length(genes_to_test)))
  x_target <- gene_data[genes_to_test, is_target, drop = FALSE]
  x_other <- gene_data[genes_to_test, !is_target, drop = FALSE]
  pct_detected_target <- 100 * Matrix::rowSums(x_target > detection_threshold) / n_target
  pct_detected_other <- 100 * Matrix::rowSums(x_other > detection_threshold) / n_other
  mean_expr_target <- Matrix::rowMeans(x_target)
  mean_expr_other <- Matrix::rowMeans(x_other)
  # Seurat v5 fold change (FoldChange.Assay, slot "data"): group means in linear (expm1) space with a pseudocount of 1 added
  # to the summed expression; replaces the former ratio of log-scale means with a 1e-10 pseudocount, which exploded when
  # expression elsewhere was ~0
  log_linear_mean <- function(m) {
    if (inherits(m, "dgCMatrix")) m@x <- expm1(m@x) else m <- expm1(as.matrix(m))
    log2((Matrix::rowSums(m) + 1) / ncol(m))
  }
  log2fc <- log_linear_mean(x_target) - log_linear_mean(x_other)
  detection_difference <- pct_detected_target - pct_detected_other
  specificity_score <- (detection_difference / 100) * log2fc
  # Wilcoxon on every prefiltered gene, then correction over the whole tested set (test before filtering)
  pvalues <- rep(NA_real_, length(genes_to_test))
  pvalues_adjusted <- rep(NA_real_, length(genes_to_test))
  if (statistical_test == "wilcox") {
    if (verbose) message("Running Wilcoxon rank-sum tests (FindMarkers) on all prefiltered genes...")
    de <- FindMarkers(seurat_obj, ident.1 = target_cluster, ident.2 = setdiff(unique(current_idents), target_cluster),
                      assay = assay_name, features = genes_to_test, test.use = "wilcox", logfc.threshold = 0, min.pct = 0,
                      only.pos = FALSE, verbose = FALSE)
    pvalues <- de$p_val[match(genes_to_test, rownames(de))]
    pvalues_adjusted <- if (pvalue_adjustment == "none") pvalues else stats::p.adjust(pvalues, method = pvalue_adjustment)
  }
  results <- data.frame(gene = genes_to_test,
                        n_cells_detected_target = round(pct_detected_target * n_target / 100),
                        pct_detected_target = round(pct_detected_target, 2),
                        n_cells_detected_other = round(pct_detected_other * n_other / 100),
                        pct_detected_other = round(pct_detected_other, 2),
                        detection_difference = round(detection_difference, 2),
                        mean_expr_target = round(mean_expr_target, 4), mean_expr_other = round(mean_expr_other, 4),
                        log2_fold_change = round(log2fc, 3), fold_change = round(2^log2fc, 2),
                        specificity_score = round(specificity_score, 4), stringsAsFactors = FALSE)
  if (statistical_test != "none") {
    results$pvalue <- pvalues
    results$pvalue_adjusted <- pvalues_adjusted
  }
  keep <- results$pct_detected_target >= min_pct_target & results$pct_detected_other <= max_pct_other &
    results$log2_fold_change >= min_log2fc & results$mean_expr_target >= min_mean_expr_target
  # A gene without a valid test result is never reported as significant
  if (statistical_test != "none") keep <- keep & !is.na(results$pvalue_adjusted) & results$pvalue_adjusted <= max_pvalue
  out <- results[keep, , drop = FALSE]
  out <- out[order(-out$specificity_score, -out$log2_fold_change), , drop = FALSE]
  if (!is.null(top_n) && nrow(out) > top_n) out <- out[seq_len(top_n), , drop = FALSE]
  # Medians only for the reported candidates (dense conversion of the full matrix is avoided)
  if (nrow(out) > 0) {
    out$median_expr_target <- round(apply(as.matrix(gene_data[out$gene, is_target, drop = FALSE]), 1, stats::median), 4)
    out$median_expr_other <- round(apply(as.matrix(gene_data[out$gene, !is_target, drop = FALSE]), 1, stats::median), 4)
  }
  rownames(out) <- NULL
  if (verbose) {
    message(sprintf("=== %d exclusive biomarker candidates ===", nrow(out)))
    if (nrow(out) == 0) message("No gene passed the filters: lower min_pct_target or min_log2fc, or raise max_pct_other")
  }
  out
}

############################## Color Palette Helper ##############################

# Generate color palettes for split.by visualization
generateSplitColorPalette <- function(n_colors, palette_name = "default") {
  # Generate color palette for DotPlot/VlnPlot split visualization
  # Args:
  #   n_colors: Number of colors needed
  #   palette_name: Name of palette to use
  # Returns:
  #   Character vector of hex colors
  
  if (n_colors <= 1) {
    return(c("lightgrey", "blue"))
  }
  
  # Add 1 for the base color (usually lightgrey or low expression)
  n_total <- n_colors + 1
  
  palette <- switch(palette_name,
                    "viridis" = {
                      viridis::viridis(n_total, option = "D")
                    },
                    "plasma" = {
                      viridis::viridis(n_total, option = "C")
                    },
                    "inferno" = {
                      viridis::viridis(n_total, option = "B")
                    },
                    "magma" = {
                      viridis::viridis(n_total, option = "A")
                    },
                    "rocket" = {
                      viridis::viridis(n_total, option = "H")
                    },
                    "turbo" = {
                      viridis::turbo(n_total)
                    },
                    "blues" = {
                      if (n_total <= 9) {
                        RColorBrewer::brewer.pal(max(3, n_total), "Blues")
                      } else {
                        colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(n_total)
                      }
                    },
                    "reds" = {
                      if (n_total <= 9) {
                        RColorBrewer::brewer.pal(max(3, n_total), "Reds")
                      } else {
                        colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))(n_total)
                      }
                    },
                    "greens" = {
                      if (n_total <= 9) {
                        RColorBrewer::brewer.pal(max(3, n_total), "Greens")
                      } else {
                        colorRampPalette(RColorBrewer::brewer.pal(9, "Greens"))(n_total)
                      }
                    },
                    "rdylbu" = {
                      if (n_total <= 11) {
                        RColorBrewer::brewer.pal(max(3, n_total), "RdYlBu")
                      } else {
                        colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(n_total)
                      }
                    },
                    "spectral" = {
                      if (n_total <= 11) {
                        RColorBrewer::brewer.pal(max(3, n_total), "Spectral")
                      } else {
                        colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(n_total)
                      }
                    },
                    "set1" = {
                      if (n_total <= 9) {
                        RColorBrewer::brewer.pal(max(3, n_total), "Set1")
                      } else {
                        colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n_total)
                      }
                    },
                    "set2" = {
                      if (n_total <= 8) {
                        RColorBrewer::brewer.pal(max(3, n_total), "Set2")
                      } else {
                        colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_total)
                      }
                    },
                    "paired" = {
                      if (n_total <= 12) {
                        RColorBrewer::brewer.pal(max(3, n_total), "Paired")
                      } else {
                        colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_total)
                      }
                    },
                    "default" = {
                      # Professional default palette
                      if (n_total <= 2) {
                        c("lightgrey", "#3182bd")
                      } else if (n_total <= 9) {
                        c("lightgrey", RColorBrewer::brewer.pal(max(3, n_total - 1), "Set2"))
                      } else {
                        c("lightgrey", colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_total - 1))
                      }
                    }
  )
  
  return(palette)
}


# ==============================================================================
# Pseudobulk differential expression (sample-level DESeq2)
# Cells are not independent replicates: cells from the same animal share its biology, so a
# cell-level test (Wilcoxon) treats thousands of cells as thousands of replicates and gives
# p-values that are far too small (pseudoreplication). Pseudobulk sums the raw counts of each
# sample (animal) and tests the samples, so the between-animal variability enters the test.
# ==============================================================================

# Number of cells per sample and condition; `condition` is "Group1", "Group2" or NA for each cell
summarize_pseudobulk_design <- function(samples, condition, min_cells = 10) {
  keep <- !is.na(condition) & !is.na(samples) & samples != ""
  if (!any(keep)) return(data.frame(sample = character(0), condition = character(0), n_cells = integer(0), used = logical(0)))
  tab <- as.data.frame(table(sample = as.character(samples[keep]), condition = condition[keep]), stringsAsFactors = FALSE)
  tab <- tab[tab$Freq > 0, , drop = FALSE]
  design <- data.frame(sample = tab$sample, condition = tab$condition, n_cells = as.integer(tab$Freq), stringsAsFactors = FALSE)
  design$used <- design$n_cells >= min_cells
  design[order(design$condition, design$sample), , drop = FALSE]
}

# NULL when the design can be tested, otherwise the message explaining why it cannot
check_pseudobulk_design <- function(design, min_samples = 2) {
  if (nrow(design) == 0) return("No cells match the selected groups and clusters.")
  shared <- unique(design$sample[duplicated(design$sample)])
  if (length(shared) > 0) return(paste0("Sample(s) present in both groups: ", paste(shared, collapse = ", "), ". The sample column must identify biological replicates (e.g. one animal) that each belong to a single group."))
  n_used <- vapply(c("Group1", "Group2"), function(g) sum(design$used & design$condition == g), integer(1))
  if (any(n_used < min_samples)) return(sprintf("Pseudobulk needs at least %d samples per group with enough cells (Group 1: %d, Group 2: %d). With one sample per group the between-sample variance cannot be estimated: use the Wilcoxon test and report it as descriptive.", min_samples, n_used[1], n_used[2]))
  NULL
}

# Sample-level DESeq2 test on a genes x cells count matrix; returns FindMarkers-compatible columns
pseudobulk_deseq2 <- function(counts, samples, condition, min_cells = 10, min_samples = 2, min_count = 10) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Package DESeq2 is not installed.")
  design <- summarize_pseudobulk_design(samples, condition, min_cells)
  problem <- check_pseudobulk_design(design, min_samples)
  if (!is.null(problem)) stop(problem)
  used <- design[design$used, , drop = FALSE]
  cells <- which(!is.na(condition) & as.character(samples) %in% used$sample)
  counts <- counts[, cells, drop = FALSE]
  sample_f <- factor(as.character(samples[cells]), levels = used$sample)
  cond <- condition[cells]
  # DESeq2 models raw integer counts: normalized or corrected values (data, SCT) are refused
  values <- if (inherits(counts, "dgCMatrix")) counts@x else as.vector(counts)
  if (any(values < 0) || any(abs(values - round(values)) > 1e-8)) stop("The counts layer does not contain raw integer counts: pseudobulk needs the raw RNA counts.")
  # Sum of the raw counts of each sample, via a sparse cells x samples indicator matrix
  indicator <- Matrix::sparseMatrix(i = seq_along(sample_f), j = as.integer(sample_f), x = 1, dims = c(length(sample_f), nlevels(sample_f)), dimnames = list(NULL, levels(sample_f)))
  pb <- as.matrix(counts %*% indicator)
  storage.mode(pb) <- "integer"
  # Genes kept when they reach min_count in at least as many samples as the smallest group
  smallest_group <- min(table(used$condition))
  pb <- pb[rowSums(pb >= min_count) >= smallest_group, , drop = FALSE]
  if (nrow(pb) < 10) stop("Too few genes pass the expression filter for a pseudobulk test.")
  col_data <- data.frame(condition = factor(used$condition, levels = c("Group2", "Group1")), row.names = used$sample)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = pb, colData = col_data, design = ~ condition)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- as.data.frame(DESeq2::results(dds, contrast = c("condition", "Group1", "Group2")))
  # Genes with NA p-value (all-zero or Cook's outliers) cannot be ranked; NA padj (independent filtering) becomes 1
  res <- res[!is.na(res$pvalue), , drop = FALSE]
  res$padj[is.na(res$padj)] <- 1
  # Cell-level detection rates, as in FindMarkers (fraction of cells with at least one count)
  genes <- rownames(res)
  g1 <- cond == "Group1"
  pct1 <- Matrix::rowMeans(counts[genes, g1, drop = FALSE] > 0)
  pct2 <- Matrix::rowMeans(counts[genes, !g1, drop = FALSE] > 0)
  out <- data.frame(avg_log2FC = res$log2FoldChange, pct.1 = round(pct1, 3), pct.2 = round(pct2, 3), p_val = res$pvalue, p_val_adj = res$padj,
                    baseMean = res$baseMean, lfcSE = res$lfcSE, row.names = genes)
  out <- out[order(out$p_val), , drop = FALSE]
  attr(out, "design") <- design
  out
}

# Pseudobulk DESeq2 on a Seurat object; `condition` gives "Group1" / "Group2" / NA for each cell (column order)
run_pseudobulk_de <- function(seurat_obj, condition, sample_col, assay = "RNA", min_cells = 10, min_samples = 2, min_count = 10) {
  if (!sample_col %in% colnames(seurat_obj@meta.data)) stop(paste("Sample column", sample_col, "not found."))
  if (!assay %in% names(seurat_obj@assays)) stop(paste("Assay", assay, "not found: pseudobulk needs the raw counts of the RNA assay."))
  cells <- colnames(seurat_obj)[!is.na(condition)]
  sub <- subset(seurat_obj, cells = cells)
  # Split counts layers (one per dataset) are joined so that all cells share one genes x cells matrix
  if (length(SeuratObject::Layers(sub[[assay]], search = "counts")) > 1) sub[[assay]] <- SeuratObject::JoinLayers(sub[[assay]])
  counts <- SeuratObject::LayerData(sub, assay = assay, layer = "counts")
  cond <- condition[match(colnames(sub), colnames(seurat_obj))]
  pseudobulk_deseq2(counts, sub@meta.data[[sample_col]], cond, min_cells = min_cells, min_samples = min_samples, min_count = min_count)
}

# Post-test display filters: BH is computed on all tested genes, then the user thresholds are applied
filter_pseudobulk_results <- function(res, logfc_threshold = 0, min_pct = 0, max_p = 0.05) {
  keep <- abs(res$avg_log2FC) >= logfc_threshold & pmax(res$pct.1, res$pct.2) >= min_pct & res$p_val < max_p
  out <- res[keep, , drop = FALSE]
  attr(out, "design") <- attr(res, "design")
  out
}

# "Group1" / "Group2" / NA label of each cell for a metadata-based comparison restricted to clusters (NULL = all clusters)
compare_cell_condition <- function(seurat_obj, metadata_col, cluster_col, group1, group2, clusters = NULL) {
  cluster_values <- if (cluster_col == "seurat_clusters") as.character(Seurat::Idents(seurat_obj)) else as.character(seurat_obj@meta.data[[cluster_col]])
  meta_values <- as.character(seurat_obj@meta.data[[metadata_col]])
  in_clusters <- if (is.null(clusters)) !is.na(cluster_values) & cluster_values != "" else cluster_values %in% clusters
  out <- rep(NA_character_, length(meta_values))
  out[in_clusters & meta_values %in% group1] <- "Group1"
  out[in_clusters & meta_values %in% group2] <- "Group2"
  out
}

# Log-normalized RNA assay for differential expression on integrated or SCTransform objects. The integrated assay holds
# batch-corrected values of the integration features only (~2,000 genes) and SCT assays with several models require
# PrepSCTFindMarkers, so DE always reads RNA. An RNA assay without a 'data' layer (SCTransform workflow) is log-normalized
# here, on this copy only; split Seurat v5 layers are joined afterwards so that FindMarkers reads one matrix.
# Returns list(obj, assay, normalized_here).
prepare_rna_assay_for_de <- function(obj, assay = "RNA") {
  if (!assay %in% names(obj@assays)) stop(sprintf("Assay '%s' not found: differential expression needs the RNA counts", assay))
  normalized_here <- length(SeuratObject::Layers(obj[[assay]], search = "data")) == 0
  if (normalized_here) obj <- Seurat::NormalizeData(obj, assay = assay, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  if (inherits(obj[[assay]], "Assay5") && length(SeuratObject::Layers(obj[[assay]], search = "data")) > 1) obj <- SeuratObject::JoinLayers(obj, assay = assay)
  SeuratObject::DefaultAssay(obj) <- assay
  return(list(obj = obj, assay = assay, normalized_here = normalized_here))
}

# Cells passing the QC filter: nFeature_RNA and percent.mt with strict inequalities (no UMI filter: nCount_RNA is
# strongly correlated with nFeature_RNA). Shared by the single- and multiple-dataset modules so that both apply the
# same rule. Returns a logical vector, one value per row of md.
qc_pass_cells <- function(md, min_features, max_features, max_mt) {
  keep <- md$nFeature_RNA > min_features & md$nFeature_RNA < max_features & md$percent.mt < max_mt
  keep[is.na(keep)] <- FALSE
  return(keep)
}
#' Average expression heatmap of clusters, optionally crossed with a second metadata axis (e.g. genotype)
#'
#' Shared by the single- and multiple-dataset modules. Genes come from a manual list, the top markers of each cluster
#' (vs all other clusters) or the top DE genes between two groups of split_by within each cluster.
#' @param cluster_values Cluster label of every cell (the column shown in the cluster selector), in cell order
#' @param clusters Clusters displayed, in display order (empty = all)
#' @param split_by Metadata column of the second axis (NULL = clusters only)
#' @param gene_source "manual", "cluster_markers" or "condition_de"
#' @param order_by "cluster" (one panel per cluster, split values inside) or "split" (one panel per split value)
#' @return List with plot, genes, markers_table (NULL in manual mode), skipped (clusters without DE result),
#'   gene_warning (message from the gene list check) and n_columns
build_grouped_expression_heatmap <- function(seurat_object, cluster_values, clusters = NULL, split_by = NULL, gene_source = "manual",
                                             gene_text = NULL, n_top = 10, group1 = NULL, group2 = NULL, direction = "both",
                                             order_by = "cluster", scale_mode = "zscore", color_palette = "viridis", assay = "RNA") {
  if (length(cluster_values) != ncol(seurat_object)) stop("cluster_values must have one value per cell")
  if (!gene_source %in% c("manual", "cluster_markers", "condition_de")) stop("Unknown gene source: ", gene_source)
  cluster_values <- as.character(cluster_values)
  all_clusters <- unique(cluster_values[!is.na(cluster_values) & cluster_values != ""])
  if (length(clusters) == 0) clusters <- all_clusters
  clusters <- clusters[clusters %in% all_clusters]
  if (length(clusters) == 0) stop("None of the selected clusters exists in the current object (renamed?): reselect them")
  if (!is.null(split_by) && !split_by %in% colnames(seurat_object@meta.data)) stop("Column ", split_by, " not found in metadata")
  if (gene_source == "condition_de") {
    if (is.null(split_by)) stop("Choose in 'Split by' the column holding the conditions (e.g. genotype)")
    if (length(group1) == 0 || length(group2) == 0) stop("Select at least one value for Group 1 and for Group 2")
    if (length(intersect(group1, group2)) > 0) stop("A value cannot be in both groups: ", paste(intersect(group1, group2), collapse = ", "))
  }
  n_top <- max(1L, as.integer(n_top))
  # Identities = the displayed cluster column, so marker and DE searches use the same clusters as the heatmap
  Idents(seurat_object) <- factor(stats::setNames(cluster_values, colnames(seurat_object)), levels = all_clusters)
  # Cells displayed: selected clusters and, with a split column, the chosen values (both groups in condition mode)
  keep <- cluster_values %in% clusters
  split_levels <- NULL
  if (!is.null(split_by)) {
    split_vals <- trimws(as.character(seurat_object@meta.data[[split_by]]))
    split_levels <- if (gene_source == "condition_de") c(group1, group2) else sort(unique(split_vals[keep & !is.na(split_vals) & split_vals != ""]))
    keep <- keep & split_vals %in% split_levels
  }
  if (!any(keep)) stop("No cells left for the selected clusters and groups")
  obj <- prepare_rna_assay_for_de(subset(seurat_object, cells = colnames(seurat_object)[keep]), assay)$obj
  # Gene selection
  markers_table <- NULL; skipped <- character(); gene_warning <- NULL; subtitle <- NULL
  if (gene_source == "cluster_markers") {
    # Markers of each cluster against all other clusters of the full object
    res <- get_top_markers(seurat_object, n_genes = n_top, assay = assay, clusters_to_include = clusters)
    genes <- res$genes; markers_table <- res$markers_table
    subtitle <- sprintf("Top %d markers per cluster (vs all other clusters)", n_top)
  } else if (gene_source == "condition_de") {
    res <- get_top_condition_de_genes(obj, split_by, group1, group2, clusters = clusters, n_genes = n_top, direction = direction, assay = assay)
    genes <- res$genes; markers_table <- res$markers_table; skipped <- res$skipped
    subtitle <- sprintf("Top %d DE genes per cluster, %s vs %s (cell-level Wilcoxon, adj. p < 0.05)", n_top, paste(group1, collapse = "+"), paste(group2, collapse = "+"))
  } else {
    check <- validate_gene_input(gene_text, obj, assay)
    if (!check$ok) stop(check$message)
    genes <- check$genes; gene_warning <- check$message
  }
  if (length(genes) == 0) stop("No gene to display")
  # Column layout: one column per cluster, or per cluster x split value with panels on the primary axis
  obj_clusters <- as.character(Idents(obj))
  if (is.null(split_by)) {
    col_id <- obj_clusters
    ids <- clusters[clusters %in% col_id]
    col_info <- data.frame(id = ids, x = ids, facet = NA_character_, stringsAsFactors = FALSE)
    x_title <- "Cluster"
  } else {
    col_id <- paste(obj_clusters, trimws(as.character(obj@meta.data[[split_by]])), sep = "___")
    # expand.grid varies its first column fastest: split values inside each cluster, or clusters inside each split value
    grid <- if (order_by == "cluster") expand.grid(split = split_levels, cluster = clusters, stringsAsFactors = FALSE) else expand.grid(cluster = clusters, split = split_levels, stringsAsFactors = FALSE)
    grid$id <- paste(grid$cluster, grid$split, sep = "___")
    grid <- grid[grid$id %in% col_id, , drop = FALSE]
    col_info <- if (order_by == "cluster") data.frame(id = grid$id, x = grid$split, facet = grid$cluster, stringsAsFactors = FALSE) else data.frame(id = grid$id, x = grid$cluster, facet = grid$split, stringsAsFactors = FALSE)
    x_title <- if (order_by == "cluster") split_by else "Cluster"
  }
  avg <- compute_group_average_expression(obj, genes, factor(col_id, levels = col_info$id), assay = assay)
  plot <- plot_grouped_expression_heatmap(avg, col_info, scale_mode = scale_mode, color_palette = color_palette,
                                          title = paste("Heatmap -", assay, "assay"), subtitle = subtitle, x_title = x_title)
  return(list(plot = plot, genes = genes, markers_table = markers_table, skipped = skipped, gene_warning = gene_warning, n_columns = nrow(col_info)))
}
# Per-cell values used to call a gene positive: raw UMI counts ("counts") or log-normalized expression ("data").
# Counts are the robust choice against ambient RNA: one ambient UMI in a small nucleus gives a high log-normalized value
# (log1p(1 / 500 * 1e4) = 3.0 for a 500-UMI nucleus), so a normalized threshold cannot separate it from true expression.
# Raw UMI counts are always read from the RNA assay: the integrated assay has no counts and SCT counts are corrected values.
fetch_threshold_values <- function(seurat_obj, genes, assay_name = "RNA", threshold_layer = "data") {
  if (!threshold_layer %in% c("counts", "data")) stop("threshold_layer must be 'counts' or 'data'")
  if (threshold_layer == "counts") {
    if ("RNA" %in% names(seurat_obj@assays)) assay_name <- "RNA"
    if (length(SeuratObject::Layers(seurat_obj[[assay_name]], search = "counts")) == 0) {
      stop(sprintf("The %s assay has no counts layer: use normalized-expression thresholds", assay_name))
    }
    missing_raw <- setdiff(genes, rownames(seurat_obj[[assay_name]]))
    if (length(missing_raw) > 0) stop(sprintf("Genes absent from the %s counts: %s", assay_name, paste(missing_raw, collapse = ", ")))
  }
  DefaultAssay(seurat_obj) <- assay_name
  # FetchData reads every split layer (merged Seurat v5 objects), one row per cell in object order
  values <- FetchData(seurat_obj, vars = genes, layer = threshold_layer)
  stopifnot("[co-expression] FetchData rows do not match the object cells" = identical(rownames(values), colnames(seurat_obj)),
            "[co-expression] FetchData is missing gene columns" = all(genes %in% colnames(values)))
  values
}

# TRUE when a cell counts as positive: at least `threshold` UMIs (counts), or expression above `threshold` (data)
call_positive <- function(values, threshold, threshold_layer = "data") {
  if (threshold_layer == "counts") values >= threshold else values > threshold
}

# Text of a threshold for tables and summaries, e.g. "Pax7 >= 2 UMI" or "Pax7 > 0.5"
format_threshold_label <- function(gene, threshold, threshold_layer = "data") {
  if (threshold_layer == "counts") sprintf("%s >= %s UMI", gene, format(threshold)) else sprintf("%s > %s", gene, format(threshold))
}

#' Background (ambient RNA) level of each gene in cells that should not express it
#'
#' Ambient RNA reaches every droplet, so its level is measured in background cells chosen by the user (a cell type known
#' not to express the genes). The threshold is set just above a high quantile of the background distribution, so that at
#' most (1 - quantile) of the background cells are still called positive.
#' @param background_cells Cell barcodes of the background population
#' @param quantile Quantile of the background distribution (default 0.95)
#' @return data.frame with gene, threshold, background_quantile_value and pct_background_positive (% of background cells
#'   still called positive at that threshold)
estimate_coexpression_thresholds <- function(seurat_obj, genes, background_cells, assay_name = "RNA", threshold_layer = "counts", quantile = 0.95) {
  background_cells <- intersect(background_cells, colnames(seurat_obj))
  if (length(background_cells) < 20) stop(sprintf("Only %d background cells: select a larger background population (>= 20 cells)", length(background_cells)))
  if (quantile <= 0 || quantile >= 1) stop("quantile must be between 0 and 1")
  values <- fetch_threshold_values(subset(seurat_obj, cells = background_cells), genes, assay_name, threshold_layer)
  out <- do.call(rbind, lapply(genes, function(g) {
    q <- as.numeric(stats::quantile(values[[g]], probs = quantile, names = FALSE, type = 7))
    # Counts: minimum number of UMIs strictly above the background quantile; data: the quantile itself (positive = above it)
    thr <- if (threshold_layer == "counts") floor(q) + 1 else round(q, 2)
    data.frame(gene = g, threshold = thr, background_quantile_value = round(q, 3),
               pct_background_positive = round(100 * mean(call_positive(values[[g]], thr, threshold_layer)), 1), stringsAsFactors = FALSE)
  }))
  out
}


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