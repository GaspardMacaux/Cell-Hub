# volcano_plot_functions.R

generateVolcanoPlot <- function(deg_results, log2fc_threshold = 0.5, pval_threshold = 0.05,
                                color_up = "#00ba38", color_down = "#f8766d", color_ns = "#619cff",
                                label_top_genes = 10, point_size = 2, point_alpha = 0.6) {
  if (is.null(deg_results) || nrow(deg_results) == 0) stop("DEG results are empty or NULL")
  required_cols <- c("avg_log2FC", "p_val_adj")
  if (!all(required_cols %in% colnames(deg_results))) stop(paste("DEG results must contain columns:", paste(required_cols, collapse = ", ")))
  plot_data <- deg_results
  plot_data$gene <- rownames(plot_data)
  plot_data$p_val_adj <- suppressWarnings(as.numeric(plot_data$p_val_adj))
  plot_data$avg_log2FC <- suppressWarnings(as.numeric(plot_data$avg_log2FC))
  plot_data <- plot_data[!is.na(plot_data$p_val_adj) & !is.na(plot_data$avg_log2FC), ]
  if (nrow(plot_data) == 0) stop("No valid data after removing NA values (tables without adjusted p-values cannot be drawn)")
  # p-values of 0 give Inf: drawn just above the largest finite value
  plot_data$neg_log10_pval <- -log10(plot_data$p_val_adj)
  finite_vals <- plot_data$neg_log10_pval[is.finite(plot_data$neg_log10_pval)]
  plot_data$neg_log10_pval[is.infinite(plot_data$neg_log10_pval)] <- if (length(finite_vals) > 0) max(finite_vals) + 50 else 350
  plot_data$regulation <- "NS"
  plot_data$regulation[plot_data$avg_log2FC > log2fc_threshold & plot_data$p_val_adj < pval_threshold] <- "Up"
  plot_data$regulation[plot_data$avg_log2FC < -log2fc_threshold & plot_data$p_val_adj < pval_threshold] <- "Down"
  plot_data$regulation <- factor(plot_data$regulation, levels = c("Up", "Down", "NS"))
  plot_data$label <- ""
  if (label_top_genes > 0) {
    top_up <- plot_data[plot_data$regulation == "Up", ]
    if (nrow(top_up) > 0) {
      top_up_genes <- head(top_up$gene[order(-top_up$avg_log2FC)], label_top_genes)
      plot_data$label[plot_data$gene %in% top_up_genes] <- plot_data$gene[plot_data$gene %in% top_up_genes]
    }
    top_down <- plot_data[plot_data$regulation == "Down", ]
    if (nrow(top_down) > 0) {
      top_down_genes <- head(top_down$gene[order(top_down$avg_log2FC)], label_top_genes)
      plot_data$label[plot_data$gene %in% top_down_genes] <- plot_data$gene[plot_data$gene %in% top_down_genes]
    }
  }
  counts <- table(plot_data$regulation)
  # Named labels and drop = FALSE: positional labels were shifted when a category was empty (NS points shown as "Down (0)")
  legend_labels <- c(Up = paste0("Up (", counts[["Up"]], ")"), Down = paste0("Down (", counts[["Down"]], ")"), NS = paste0("NS (", counts[["NS"]], ")"))
  p <- ggplot(plot_data, aes(x = avg_log2FC, y = neg_log10_pval, color = regulation)) +
    geom_point(size = point_size, alpha = point_alpha) +
    scale_color_manual(values = c(Up = color_up, Down = color_down, NS = color_ns), labels = legend_labels, breaks = c("Up", "Down", "NS"), drop = FALSE) +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "gray40", linewidth = 0.5) +
    labs(x = "log2 Fold Change", y = "-log10 (Adjusted p-value)", color = "Regulation", title = "Volcano Plot - Differential Gene Expression") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), legend.position = "right", panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  if (any(plot_data$label != "")) {
    p <- p + ggrepel::geom_text_repel(data = plot_data[plot_data$label != "", ], aes(label = label), size = 3, box.padding = 0.5,
                                      point.padding = 0.3, max.overlaps = 20, color = "black")
  }
  return(p)
}

