# cells_genes_expressions_newarch.R
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
    # Join layers before extraction — merged objects have fragmented layers
    # that cause LayerData to fail on the "data" slot (rules #1 and #2)
    seurat_obj <- JoinLayers(seurat_obj)
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
    if (is_integrated && "dataset" %in% colnames(seurat_obj@meta.data)) {
      cluster_info$dataset <- seurat_obj@meta.data$dataset
    }
    if (is_integrated) {
      cluster_counts_total <- cluster_info %>%
        group_by(cluster, dataset) %>%
        summarise(n = n(), .groups = 'drop')
    } else {
      cluster_counts_total <- table(cluster_info$cluster)
    }
    expression_summary_list <- lapply(selected_genes, function(gene) {
      gene_expression <- gene_data[gene, ]
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


# Enhanced function to analyze gene co-expression with better visualizations
analyze_gene_coexpression <- function(seurat_obj, genes, assay_name = "RNA",
                                      expression_thresholds = NULL, is_integrated = FALSE) {
  tryCatch({
    if (is.character(genes) && length(genes) == 1) {
      genes <- unique(trimws(unlist(strsplit(genes, ","))))
    }
    genes <- genes[genes != ""]
    if (length(genes) < 2) stop("At least 2 genes are required for co-expression analysis")
    
    if (is.null(expression_thresholds)) {
      expression_thresholds <- rep(0, length(genes))
      names(expression_thresholds) <- genes
    } else if (length(expression_thresholds) == 1) {
      single_threshold <- expression_thresholds
      expression_thresholds <- rep(single_threshold, length(genes))
      names(expression_thresholds) <- genes
    } else if (is.null(names(expression_thresholds))) {
      if (length(expression_thresholds) != length(genes)) {
        stop("Number of thresholds must match number of genes or be a single value")
      }
      names(expression_thresholds) <- genes
    }
    missing_thresholds <- setdiff(genes, names(expression_thresholds))
    for (gene in missing_thresholds) expression_thresholds[gene] <- 0
    
    available_genes <- rownames(seurat_obj[[assay_name]])
    missing_genes <- setdiff(genes, available_genes)
    if (length(missing_genes) > 0) {
      warning(paste("Genes not found:", paste(missing_genes, collapse = ", ")))
      genes <- intersect(genes, available_genes)
      expression_thresholds <- expression_thresholds[genes]
    }
    if (length(genes) < 2) stop("Less than 2 valid genes remaining after filtering")
    
    DefaultAssay(seurat_obj) <- assay_name
    if (inherits(seurat_obj[[assay_name]], "Assay5")) {
      seurat_obj <- JoinLayers(seurat_obj, assay = assay_name)
    }
    expression_data <- FetchData(seurat_obj, vars = genes)
    
    gene_binary <- matrix(FALSE, nrow = length(genes), ncol = ncol(seurat_obj))
    rownames(gene_binary) <- genes
    colnames(gene_binary) <- colnames(seurat_obj)
    for (gene in genes) {
      gene_binary[gene, ] <- expression_data[, gene] > expression_thresholds[gene]
    }
    
    metadata <- data.frame(
      cell    = colnames(seurat_obj),
      cluster = as.character(Idents(seurat_obj)),
      stringsAsFactors = FALSE
    )
    if (is_integrated && "dataset" %in% colnames(seurat_obj@meta.data)) {
      metadata$dataset <- seurat_obj@meta.data$dataset
    }
    
    analyze_group <- function(group_cells, group_name, dataset_name = NULL) {
      group_cells <- intersect(group_cells, colnames(gene_binary))
      if (length(group_cells) == 0) return(NULL)
      group_binary <- gene_binary[, group_cells, drop = FALSE]
      n_total <- length(group_cells)
      results <- data.frame(Group = group_name, N_Cells_Total = n_total, stringsAsFactors = FALSE)
      if (!is.null(dataset_name)) results$Dataset <- dataset_name
      
      for (gene in genes) {
        n_pos <- sum(group_binary[gene, ])
        pct   <- round(100 * n_pos / n_total, 1)
        results[[paste0("N_",   gene, "_Positive")]] <- n_pos
        results[[paste0("Pct_", gene, "_Positive")]] <- pct
        results[[paste0(gene,   "_Threshold")]]       <- expression_thresholds[gene]
        results[[paste0(gene,   "_Bar")]]             <- create_expression_percentage_bar(pct)
      }
      
      if (length(genes) == 2) {
        b1 <- group_binary[genes[1], ]
        b2 <- group_binary[genes[2], ]
        n_both    <- sum( b1 &  b2)
        n_only1   <- sum( b1 & !b2)
        n_only2   <- sum(!b1 &  b2)
        n_neither <- sum(!b1 & !b2)
        results$N_Both_Positive <- n_both
        results$N_Only_Gene1    <- n_only1
        results$N_Only_Gene2    <- n_only2
        results$N_Neither       <- n_neither
        results$Pct_Both       <- round(100 * n_both    / n_total, 1)
        results$Pct_Only_Gene1 <- round(100 * n_only1   / n_total, 1)
        results$Pct_Only_Gene2 <- round(100 * n_only2   / n_total, 1)
        results$Pct_Neither    <- round(100 * n_neither / n_total, 1)
        results$Coexpression_Summary <- paste0(
          "Both (>", expression_thresholds[genes[1]], " & >", expression_thresholds[genes[2]], "): ",
          n_both,  " (", results$Pct_Both,       "%) | ",
          "Only ", genes[1], ": ", n_only1, " (", results$Pct_Only_Gene1, "%) | ",
          "Only ", genes[2], ": ", n_only2, " (", results$Pct_Only_Gene2, "%)"
        )
      } else {
        n_genes_expressed       <- colSums(group_binary)
        results$N_All_Positive  <- sum(n_genes_expressed == length(genes))
        results$N_Any_Positive  <- sum(n_genes_expressed > 0)
        results$N_None_Positive <- sum(n_genes_expressed == 0)
        results$Pct_All  <- round(100 * results$N_All_Positive  / n_total, 1)
        results$Pct_Any  <- round(100 * results$N_Any_Positive  / n_total, 1)
        results$Pct_None <- round(100 * results$N_None_Positive / n_total, 1)
        for (i in 0:length(genes)) {
          results[[paste0("N_Expressing_", i, "_Genes")]] <- sum(n_genes_expressed == i)
        }
        threshold_info <- paste(genes, ">", expression_thresholds[genes], collapse = ", ")
        results$Coexpression_Summary <- paste0(
          "All (", threshold_info, "): ", results$N_All_Positive, " (", results$Pct_All, "%) | ",
          "Any: ", results$N_Any_Positive, " (", results$Pct_Any, "%)"
        )
      }
      return(results)
    }
    
    if (is_integrated) {
      cluster_results <- metadata %>%
        group_by(cluster, dataset) %>%
        group_map(function(.x, .y) {
          analyze_group(
            group_cells  = .x$cell,
            group_name   = as.character(.y$cluster),
            dataset_name = as.character(.y$dataset)
          )
        })
      cluster_results <- Filter(Negate(is.null), cluster_results)
      
      overall_results <- metadata %>%
        group_by(dataset) %>%
        group_map(function(.x, .y) {
          analyze_group(
            group_cells  = .x$cell,
            group_name   = paste("Overall", as.character(.y$dataset)),
            dataset_name = as.character(.y$dataset)
          )
        })
      overall_results <- Filter(Negate(is.null), overall_results)
      
      all_results <- bind_rows(c(cluster_results, overall_results))
      
    } else {
      all_clusters <- unique(metadata$cluster)
      cluster_results <- lapply(all_clusters, function(cl) {
        analyze_group(
          group_cells = metadata$cell[metadata$cluster == cl],
          group_name  = cl
        )
      })
      cluster_results <- Filter(Negate(is.null), cluster_results)
      overall_result  <- analyze_group(metadata$cell, "Overall")
      all_results <- do.call(rbind, c(list(overall_result), cluster_results))
    }
    
    if (length(genes) == 2) {
      all_results$Coexpression_Visual <- mapply(
        function(both, only1, only2, neither) create_coexpression_bar(both, only1, only2, neither),
        all_results$Pct_Both,
        all_results$Pct_Only_Gene1,
        all_results$Pct_Only_Gene2,
        all_results$Pct_Neither
      )
    }
    
    return(list(
      data            = all_results,
      genes_analyzed  = genes,
      thresholds_used = expression_thresholds,
      is_integrated   = is_integrated
    ))
    
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

# Function to create cluster composition table with support for integrated datasets
create_cluster_composition_table <- function(seurat_obj, is_integrated = FALSE, metadata_column = "dataset") {
  # Get cluster information
  cluster_data <- data.frame(
    cluster = as.character(Idents(seurat_obj)),
    nFeature = seurat_obj$nFeature_RNA,
    nCount = seurat_obj$nCount_RNA,
    percent_mt = seurat_obj$percent.mt,
    stringsAsFactors = FALSE
  )
  
  # Add metadata grouping info for integrated data
  if (is_integrated && metadata_column %in% colnames(seurat_obj@meta.data)) {
    cluster_data$grouping <- seurat_obj@meta.data[[metadata_column]]
  } else if (is_integrated) {
    warning(paste("Metadata column", metadata_column, "not found. Proceeding without grouping."))
    is_integrated <- FALSE
  }
  
  # Calculate statistics
  if (is_integrated) {
    # Group by cluster and metadata column
    cluster_stats <- cluster_data %>%
      group_by(cluster, grouping) %>%
      summarise(
        Cell_Count = n(),
        Mean_Genes = round(mean(nFeature, na.rm = TRUE), 1),
        Median_Genes = round(median(nFeature, na.rm = TRUE), 1),
        Mean_UMI = round(mean(nCount, na.rm = TRUE), 1),
        Median_UMI = round(median(nCount, na.rm = TRUE), 1),
        Mean_MT_Percent = round(mean(percent_mt, na.rm = TRUE), 2),
        .groups = 'drop'
      )
    
    # Calculate percentages within each metadata group
    group_totals <- cluster_data %>%
      group_by(grouping) %>%
      summarise(Total_in_Group = n(), .groups = 'drop')
    
    cluster_stats <- merge(cluster_stats, group_totals, by = "grouping")
    cluster_stats$Percent_of_Group <- round(
      (cluster_stats$Cell_Count / cluster_stats$Total_in_Group) * 100, 1
    )
    
    # Rename grouping column to reflect actual metadata column name
    colnames(cluster_stats)[colnames(cluster_stats) == "grouping"] <- metadata_column
    colnames(cluster_stats)[colnames(cluster_stats) == "Total_in_Group"] <- paste0("Total_in_", metadata_column)
    colnames(cluster_stats)[colnames(cluster_stats) == "Percent_of_Group"] <- paste0("Percent_of_", metadata_column)
    
  } else {
    # Standard single dataset processing
    cluster_stats <- cluster_data %>%
      group_by(cluster) %>%
      summarise(
        Cell_Count = n(),
        Mean_Genes = round(mean(nFeature, na.rm = TRUE), 1),
        Median_Genes = round(median(nFeature, na.rm = TRUE), 1),
        Mean_UMI = round(mean(nCount, na.rm = TRUE), 1),
        Median_UMI = round(median(nCount, na.rm = TRUE), 1),
        Mean_MT_Percent = round(mean(percent_mt, na.rm = TRUE), 2),
        .groups = 'drop'
      )
  }
  
  # Calculate total cells for percentage
  total_cells <- nrow(cluster_data)
  cluster_stats$Percent_of_Total <- round((cluster_stats$Cell_Count / total_cells) * 100, 1)
  
  # Create visual representation of cluster sizes
  cluster_stats$Size_Bar <- sapply(cluster_stats$Percent_of_Total, create_size_bar)
  
  # Create quality summary
  cluster_stats$Quality_Summary <- paste0(
    "Genes: ", cluster_stats$Mean_Genes, " | ",
    "UMI: ", cluster_stats$Mean_UMI, " | ",
    "MT: ", cluster_stats$Mean_MT_Percent, "%"
  )
  
  # Reorder columns based on dataset type
  if (is_integrated) {
    # Dynamic column ordering based on metadata_column
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

# Enhanced render function for co-expression table
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


# Function to create co-expression plot
create_coexpression_plot <- function(coexpr_data, genes_analyzed,
                                     group_by = "cluster", secondary_split = "none") { 
  tryCatch({
    if (length(genes_analyzed) == 2) {
      plot_data <- coexpr_data %>%
        filter(Group != "Overall") %>%
        select(Group, Pct_Both, Pct_Only_Gene1, Pct_Only_Gene2, Pct_Neither) %>%
        pivot_longer(cols = starts_with("Pct_"),
                     names_to = "Category",
                     values_to = "Percentage") %>%
        mutate(
          Category = case_when(
            Category == "Pct_Both"       ~ paste("Both", genes_analyzed[1], "&", genes_analyzed[2]),
            Category == "Pct_Only_Gene1" ~ paste("Only", genes_analyzed[1]),
            Category == "Pct_Only_Gene2" ~ paste("Only", genes_analyzed[2]),
            Category == "Pct_Neither"    ~ "Neither"
          ),
          Category = factor(Category, levels = c(
            paste("Both", genes_analyzed[1], "&", genes_analyzed[2]),
            paste("Only", genes_analyzed[1]),
            paste("Only", genes_analyzed[2]),
            "Neither"
          ))
        )
      p <- ggplot(plot_data, aes(x = Group, y = Percentage, fill = Category)) +
        geom_col(position = "stack") +
        scale_fill_manual(values = c("#2E8B57", "#4169E1", "#DC143C", "#D3D3D3")) +
        labs(
          title = paste("Co-expression Analysis:", genes_analyzed[1], "vs", genes_analyzed[2]),
          x = "Group",
          y = "Percentage of Cells (%)",
          fill = "Expression Pattern"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      # For multiple genes: percentage computed per group using N_Cells_Total
      # N_Expressing_*_Genes columns count cells per bucket; divide by group total, not global sum
      plot_data <- coexpr_data %>%
        filter(Group != "Overall") %>%
        select(Group, N_Cells_Total, starts_with("N_Expressing_")) %>%
        pivot_longer(cols = starts_with("N_Expressing_"),
                     names_to = "N_Genes",
                     values_to = "N_Cells") %>%
        mutate(
          N_Genes = as.numeric(gsub("N_Expressing_(.*)_Genes", "\\1", N_Genes)),
          Percentage = N_Cells / N_Cells_Total * 100
        )
      p <- ggplot(plot_data, aes(x = Group, y = Percentage, fill = factor(N_Genes))) +
        geom_col(position = "stack") +
        scale_fill_viridis_d(name = "Number of\nGenes Expressed") +
        labs(
          title = paste("Co-expression Analysis:", paste(genes_analyzed, collapse = ", ")),
          x = "Group",
          y = "Percentage of Cells (%)"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    return(p)
  }, error = function(e) {
    stop(paste("Error creating co-expression plot:", e$message))
  })
}

# Within-cluster DE: co-expressing cells (all genes > threshold) vs non-co-expressing
run_coexpression_de <- function(seurat_obj, genes, cluster, assay_name = "RNA",
                                expression_threshold = 0, min_pct = 0.1,
                                logfc_threshold = 0.1) {
  if (inherits(seurat_obj[[assay_name]], "Assay5")) {
    tryCatch(seurat_obj <- JoinLayers(seurat_obj, assay = assay_name), error = function(e) NULL)
  }
  DefaultAssay(seurat_obj) <- assay_name
  expr <- FetchData(seurat_obj, vars = genes)
  
  cluster_cells <- if (identical(cluster, "__all__")) {
    colnames(seurat_obj)
  } else {
    colnames(seurat_obj)[as.character(Idents(seurat_obj)) == as.character(cluster)]
  }
  if (length(cluster_cells) == 0)
    stop(sprintf("No cells found in cluster '%s'", cluster))
  
  expr_cluster  <- expr[cluster_cells, , drop = FALSE]
  expressed     <- expr_cluster > expression_threshold
  n_genes       <- length(genes)
  
  coexpr_mask           <- rowSums(expressed) == n_genes
  coexpressor_cells     <- cluster_cells[coexpr_mask]
  non_coexpressor_cells <- cluster_cells[!coexpr_mask]
  
  # Detailed breakdown of non-co-expressors
  breakdown <- if (n_genes == 2) {
    only_g1  <- cluster_cells[expressed[, 1] & !expressed[, 2]]
    only_g2  <- cluster_cells[!expressed[, 1] & expressed[, 2]]
    neither  <- cluster_cells[!expressed[, 1] & !expressed[, 2]]
    list(
      only_gene1 = length(only_g1),
      only_gene2 = length(only_g2),
      neither    = length(neither)
    )
  } else {
    partial <- cluster_cells[rowSums(expressed) > 0 & !coexpr_mask]
    none    <- cluster_cells[rowSums(expressed) == 0]
    list(partial = length(partial), none = length(none))
  }
  
  if (length(coexpressor_cells) < 3)
    stop(sprintf("Only %d co-expressing cells (min 3). Lower threshold or pick another cluster.",
                 length(coexpressor_cells)))
  if (length(non_coexpressor_cells) < 3)
    stop(sprintf("Only %d non-co-expressing cells (min 3).", length(non_coexpressor_cells)))
  
  seurat_obj$coexpr_de_group <- NA_character_
  seurat_obj$coexpr_de_group[coexpressor_cells]     <- "coexpressor"
  seurat_obj$coexpr_de_group[non_coexpressor_cells] <- "non_coexpressor"
  
  markers <- FindMarkers(
    seurat_obj,
    group.by        = "coexpr_de_group",
    ident.1         = "coexpressor",
    ident.2         = "non_coexpressor",
    assay           = assay_name,
    test.use        = "wilcox",
    only.pos        = FALSE,
    min.pct         = min_pct,
    logfc.threshold = logfc_threshold
  )
  markers$gene <- rownames(markers)
  markers      <- markers[order(markers$avg_log2FC, decreasing = TRUE), ]
  
  list(
    markers            = markers,
    genes_used         = genes,
    cluster            = cluster,
    n_coexpressors     = length(coexpressor_cells),
    n_non_coexpressors = length(non_coexpressor_cells),
    breakdown          = breakdown,
    threshold_used     = expression_threshold
  )
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


# Wrapper functions for backward compatibility
create_single_cluster_composition_table <- function(seurat_obj) {
  create_cluster_composition_table(seurat_obj, is_integrated = FALSE)
}

analyze_gene_coexpression_single <- function(seurat_obj, gene_text, expression_threshold = 0) {
  analyze_gene_coexpression(
    seurat_obj            = seurat_obj,
    genes                 = gene_text,
    assay_name            = DefaultAssay(seurat_obj),
    expression_thresholds = expression_threshold,
    is_integrated         = FALSE
  )
}

create_coexpression_plot_single <- function(results_data, genes) {
  create_coexpression_plot(results_data, genes)
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

# Format p-values for display, handling extreme values properly
# Supports p-values down to machine precision limit (~10^-308)
# @param p_values: numeric vector of p-values
# @return: character vector with properly formatted p-values
format_pvalue_robust <- function(p_values) {
  sapply(p_values, function(p) {
    # Handle NA values
    if (is.na(p)) {
      return("NA")
    }
    
    # Handle exact zero - this means underflow below machine precision
    # For double precision, this is below ~2.225e-308
    if (p == 0) {
      return("< 2.23e-308")
    }
    
    # For p-values < 0.01, always use scientific notation
    if (p < 0.01) {
      return(sprintf("%.2e", p))
    }
    
    # For larger p-values (>= 0.01), use fixed decimal notation
    return(sprintf("%.4f", p))
  })
}

#' Find exclusive biomarkers for a specific cluster or group of clusters using robust statistical approach
#'
#' Identifies genes that are predominantly expressed in target cluster(s)
#' using a combination of detection rate, expression magnitude, and statistical testing
#'
#' @param seurat_obj Seurat object containing clustered data
#' @param target_cluster Character vector of cluster identifier(s) to find biomarkers for
#' @param min_pct_target Minimum percentage of cells with detectable expression in target cluster (0-100, default: 50)
#' @param max_pct_other Maximum percentage of cells with detectable expression in other clusters (0-100, default: 25)
#' @param min_log2fc Minimum log2 fold change between target and other clusters (default: 1.5)
#' @param detection_threshold Minimum expression value to consider a gene as detected (default: 0, any expression > 0)
#' @param min_mean_expr_target Minimum mean expression level in target cluster (default: 0.5)
#' @param statistical_test Type of test to use: "wilcox" (Wilcoxon rank-sum) or "none" (default: "wilcox")
#' @param max_pvalue Maximum adjusted p-value to consider significant (default: 0.05, only used if statistical_test != "none")
#' @param pvalue_adjustment Method for p-value adjustment: "bonferroni", "BH" (Benjamini-Hochberg), "fdr", "none" (default: "BH")
#' @param assay_name Name of assay to use (default: "RNA")
#' @param top_n Number of top markers to return (default: NULL for all passing filters)
#' @param verbose Print detailed progress messages (default: TRUE)
#'
#' @return Data frame with exclusive biomarker candidates sorted by specificity score
#' @export
find_exclusive_biomarkers <- function(seurat_obj,
                                      target_cluster,
                                      min_pct_target = 50,
                                      max_pct_other = 25,
                                      min_log2fc = 1.5,
                                      detection_threshold = 0,
                                      min_mean_expr_target = 0.5,
                                      statistical_test = "wilcox",
                                      max_pvalue = 0.05,
                                      pvalue_adjustment = "BH",
                                      assay_name = "RNA",
                                      top_n = NULL,
                                      verbose = TRUE) {
  
  current_idents <- as.character(Idents(seurat_obj))
  target_cluster <- as.character(target_cluster)
  
  missing_clusters <- setdiff(target_cluster, unique(current_idents))
  if (length(missing_clusters) > 0) {
    stop(paste0("Target cluster(s) not found in active identities: ",
                paste(missing_clusters, collapse = ", "),
                "\nAvailable clusters: ",
                paste(unique(current_idents), collapse = ", ")))
  }
  
  if (min_pct_target < 0 || min_pct_target > 100) stop("min_pct_target must be between 0 and 100")
  if (max_pct_other < 0 || max_pct_other > 100)   stop("max_pct_other must be between 0 and 100")
  if (!statistical_test %in% c("wilcox", "none"))  stop("statistical_test must be 'wilcox' or 'none'")
  if (!pvalue_adjustment %in% c("bonferroni", "BH", "fdr", "none")) {
    stop("pvalue_adjustment must be 'bonferroni', 'BH', 'fdr', or 'none'")
  }
  
  if (verbose) {
    message("=== Finding exclusive biomarkers ===")
    if (length(target_cluster) == 1) {
      message(paste0("Target cluster: ", target_cluster))
    } else {
      message(paste0("Target cluster group: ", paste(target_cluster, collapse = ", ")))
    }
    message(paste0("Minimum % detected in target: ", min_pct_target, "%"))
    message(paste0("Maximum % detected in others: ", max_pct_other, "%"))
    message(paste0("Detection threshold: ", detection_threshold))
    message(paste0("Minimum mean expression in target: ", min_mean_expr_target))
    message(paste0("Minimum log2FC: ", min_log2fc))
    message(paste0("Statistical test: ", statistical_test))
    if (statistical_test != "none") {
      message(paste0("P-value adjustment: ", pvalue_adjustment))
      message(paste0("Maximum adjusted p-value: ", max_pvalue))
    }
  }
  
  DefaultAssay(seurat_obj) <- assay_name
  
  # JoinLayers silently ignored for legacy Assay objects (no method available).
  # GetAssayData with layer= works for both Assay and Assay5 in SeuratObject v5.
  tryCatch(
    seurat_obj <- JoinLayers(seurat_obj, assay = assay_name),
    error = function(e) NULL
  )
  gene_data <- GetAssayData(seurat_obj, assay = assay_name, layer = "data")
  
  cluster_info <- data.frame(
    cell    = colnames(seurat_obj),
    cluster = current_idents,
    stringsAsFactors = FALSE
  )
  
  target_cells <- cluster_info$cell[cluster_info$cluster %in% target_cluster]
  other_cells  <- cluster_info$cell[!cluster_info$cluster %in% target_cluster]
  n_target <- length(target_cells)
  n_other  <- length(other_cells)
  
  if (verbose) {
    message(paste0("Analyzing ", n_target, " cells in target cluster(s) vs ",
                   n_other, " cells in other clusters"))
  }
  
  all_genes     <- rownames(gene_data)
  all_genes     <- all_genes[!is.na(all_genes)]  # drop unnamed genes
  gene_data     <- gene_data[all_genes, , drop = FALSE]
  n_genes_total <- length(all_genes)
  if (verbose) message(paste0("Total genes in dataset: ", n_genes_total))
  
  if (verbose) message("Pre-filtering genes based on target cluster expression...")
  
  target_data            <- gene_data[, target_cells, drop = FALSE]
  mean_expr_in_target    <- Matrix::rowMeans(target_data)
  n_detected_in_target   <- Matrix::rowSums(target_data > detection_threshold)
  pct_detected_in_target <- (n_detected_in_target / n_target) * 100
  
  prefilter_min_pct  <- max(0, min_pct_target - 10)
  prefilter_min_expr <- max(0, min_mean_expr_target * 0.5)
  
  genes_pass_target_filter <- (
    pct_detected_in_target >= prefilter_min_pct &
      mean_expr_in_target    >= prefilter_min_expr
  )
  
  genes_to_test <- all_genes[genes_pass_target_filter]
  
  if (verbose) {
    message(paste0("After target cluster pre-filtering: ", length(genes_to_test), " genes (",
                   round(100 * length(genes_to_test) / n_genes_total, 1), "% of total)"))
    message(paste0("Excluded ", n_genes_total - length(genes_to_test),
                   " genes with low expression in target cluster(s)"))
  }
  
  if (length(genes_to_test) == 0) {
    warning("No genes passed the pre-filtering step. Try relaxing min_pct_target or min_mean_expr_target.")
    return(data.frame())
  }
  
  gene_data_filtered <- gene_data[genes_to_test, , drop = FALSE]
  pct_detected_target <- pct_detected_in_target[genes_pass_target_filter]
  mean_expr_target    <- mean_expr_in_target[genes_pass_target_filter]
  n_detected_target   <- n_detected_in_target[genes_pass_target_filter]
  
  if (verbose) message("Calculating expression statistics for other clusters...")
  other_expr_matrix  <- gene_data_filtered[, other_cells, drop = FALSE] > detection_threshold
  n_detected_other   <- Matrix::rowSums(other_expr_matrix)
  pct_detected_other <- (n_detected_other / n_other) * 100
  mean_expr_other    <- Matrix::rowMeans(gene_data_filtered[, other_cells, drop = FALSE])
  
  if (verbose) message("Calculating median expression...")
  median_expr_target   <- apply(as.matrix(gene_data_filtered[, target_cells, drop = FALSE]), 1, median)
  median_expr_other    <- apply(as.matrix(gene_data_filtered[, other_cells,  drop = FALSE]), 1, median)
  log2fc               <- log2((mean_expr_target + 1e-10) / (mean_expr_other + 1e-10))
  log2fc[is.nan(log2fc) | is.infinite(log2fc)] <- NA_real_
  detection_difference <- pct_detected_target - pct_detected_other
  specificity_score    <- (detection_difference / 100) * log2fc
  
  basic_filter <- (
    pct_detected_target >= min_pct_target &
      pct_detected_other  <= max_pct_other  &
      !is.na(log2fc) & log2fc >= min_log2fc &
      mean_expr_target    >= min_mean_expr_target
  )
  
  genes_for_testing <- genes_to_test[basic_filter]
  if (verbose) {
    message(paste0("After basic filtering: ", length(genes_for_testing), " genes pass criteria"))
    message(paste0("Will perform statistical tests on these ", length(genes_for_testing), " genes only"))
  }
  
  pvalues <- rep(NA, length(genes_to_test))
  if (statistical_test == "wilcox" && length(genes_for_testing) > 0) {
    if (verbose) message("Performing Wilcoxon rank-sum tests on filtered gene set...")
    indices_to_test <- which(basic_filter)
    pvalues_tested <- sapply(seq_along(indices_to_test), function(i) {
      if (i %% 500 == 0 && verbose) {
        message(paste0("  Tested ", i, "/", length(indices_to_test), " genes..."))
      }
      idx           <- indices_to_test[i]
      gene          <- genes_to_test[idx]
      target_values <- as.numeric(gene_data_filtered[gene, target_cells])
      other_values  <- as.numeric(gene_data_filtered[gene, other_cells])
      if (var(target_values) > 0 && var(other_values) > 0) {
        tryCatch({
          test_result <- wilcox.test(target_values, other_values, alternative = "greater")
          return(test_result$p.value)
        }, error = function(e) return(NA))
      } else {
        return(NA)
      }
    })
    pvalues[indices_to_test] <- pvalues_tested
    if (pvalue_adjustment != "none") {
      if (verbose) message(paste0("Adjusting p-values using ", pvalue_adjustment, " method..."))
      valid_pvalues <- !is.na(pvalues)
      if (sum(valid_pvalues) > 0) {
        pvalues_adjusted <- rep(NA, length(pvalues))
        pvalues_adjusted[valid_pvalues] <- p.adjust(pvalues[valid_pvalues], method = pvalue_adjustment)
      } else {
        pvalues_adjusted <- pvalues
      }
    } else {
      pvalues_adjusted <- pvalues
    }
  } else {
    pvalues_adjusted <- rep(NA, length(genes_to_test))
  }
  
  results <- data.frame(
    gene                   = genes_to_test,
    n_cells_detected_target = n_detected_target,
    pct_detected_target    = round(pct_detected_target, 2),
    n_cells_detected_other = n_detected_other,
    pct_detected_other     = round(pct_detected_other, 2),
    detection_difference   = round(detection_difference, 2),
    mean_expr_target       = round(mean_expr_target, 4),
    mean_expr_other        = round(mean_expr_other, 4),
    median_expr_target     = round(median_expr_target, 4),
    median_expr_other      = round(median_expr_other, 4),
    log2_fold_change       = round(log2fc, 3),
    fold_change            = round(2^log2fc, 2),
    specificity_score      = round(specificity_score, 4),
    stringsAsFactors = FALSE
  )
  
  if (statistical_test != "none") {
    results$pvalue          <- pvalues
    results$pvalue_adjusted <- pvalues_adjusted
  }
  
  filter_conditions <- (
    !is.na(results$log2_fold_change) &
      results$pct_detected_target >= min_pct_target &
      results$pct_detected_other  <= max_pct_other  &
      results$log2_fold_change    >= min_log2fc      &
      results$mean_expr_target    >= min_mean_expr_target
  )
  
  if (statistical_test != "none") {
    filter_conditions <- filter_conditions & (results$pvalue_adjusted <= max_pvalue | is.na(results$pvalue_adjusted))
  }
  
  exclusive_markers <- results[filter_conditions, ]
  exclusive_markers <- exclusive_markers[order(-exclusive_markers$specificity_score,
                                               -exclusive_markers$log2_fold_change), ]
  if (!is.null(top_n) && nrow(exclusive_markers) > top_n) {
    exclusive_markers <- exclusive_markers[1:top_n, ]
  }
  rownames(exclusive_markers) <- NULL
  
  if (verbose) {
    message(paste0("=== Found ", nrow(exclusive_markers), " exclusive biomarkers ==="))
    if (nrow(exclusive_markers) > 0) {
      message("\nTop 5 markers:")
      top5 <- head(exclusive_markers, 5)
      for (i in 1:nrow(top5)) {
        msg <- sprintf("  %d. %s - %.1f%% target vs %.1f%% others (log2FC=%.2f, FC=%.1fx)",
                       i, top5$gene[i],
                       top5$pct_detected_target[i],
                       top5$pct_detected_other[i],
                       top5$log2_fold_change[i],
                       top5$fold_change[i])
        if (statistical_test != "none" && !is.na(top5$pvalue_adjusted[i])) {
          msg <- paste0(msg, sprintf(", adj.p=%.2e", top5$pvalue_adjusted[i]))
        }
        message(msg)
      }
    } else {
      message("\nNo genes passed the filtering criteria.")
      message("Consider relaxing the parameters:")
      message("  - Decrease min_pct_target")
      message("  - Increase max_pct_other")
      message("  - Decrease min_log2fc")
      message("  - Decrease min_mean_expr_target")
      if (statistical_test != "none") message("  - Increase max_pvalue")
    }
  }
  
  return(exclusive_markers)
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


