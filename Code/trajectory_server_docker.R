############################## Trajectory analysis (Monocle 3) ##############################

# trajectory_server.R
# Reactive logic of the trajectory module; computations and plots are in functions/trajectory_functions.R

trajectory_server <- function(input, output, session) {
  
  shinyjs::useShinyjs() # button deactivation
  
  ############################## Reactive state ##############################
  
  seurat_monocle <- reactiveVal(NULL)          # Seurat object to convert
  monocle_object <- reactiveVal(NULL)          # current cell_data_set (graph and selected pseudotime)
  pseudotime_history <- reactiveVal(list())    # one entry per root
  selected_pseudotime_id <- reactiveVal(NULL)
  trajectory_plot <- reactiveVal(NULL)
  pseudotime_plot <- reactiveVal(NULL)
  root_info_text <- reactiveVal(NULL)
  pseudotime_diff_gene_results <- reactiveVal(NULL)  # graph_test
  correlation_results <- reactiveVal(NULL)           # correlation with pseudotime
  correlation_caption <- reactiveVal(NULL)
  single_gene_trajectory_plot <- reactiveVal(NULL)
  gene_trajectory_plot <- reactiveVal(NULL)
  gene_trajectory_n <- reactiveVal(0)
  
  ############################## Download buttons state ##############################
  # Each download button stays disabled until the result it exports exists
  bind_download_state("download_trajectory_umap", function() req(trajectory_plot()))
  bind_download_state("download_pseudotime_umap", function() req(pseudotime_plot()))
  bind_download_state("download_graph_test", function() req(pseudotime_diff_gene_results()))
  bind_download_state("download_correlation", function() req(correlation_results()))
  bind_download_state("download_single_gene_umap_trajectory", function() req(single_gene_trajectory_plot()))
  bind_download_state("download_multi_gene_pseudotime_trajectory", function() req(gene_trajectory_plot()))
  
  # Buttons are enabled step by step
  shinyjs::disable("constructGraph")
  shinyjs::disable("set_root_cell")
  
  # Selected pseudotime entry (NULL when none)
  selected_pseudotime_entry <- function() {
    id <- selected_pseudotime_id()
    if (is.null(id)) return(NULL)
    pseudotime_history()[[id]]
  }
  
  # Clear every result of a previous object before loading a new one
  reset_trajectory_state <- function() {
    monocle_object(NULL); seurat_monocle(NULL)
    trajectory_plot(NULL); pseudotime_plot(NULL); root_info_text(NULL)
    pseudotime_history(list()); selected_pseudotime_id(NULL)
    pseudotime_diff_gene_results(NULL); correlation_results(NULL); correlation_caption(NULL)
    single_gene_trajectory_plot(NULL); gene_trajectory_plot(NULL); gene_trajectory_n(0)
    updatePickerInput(session, "multi_gene_pseudotime_trajectory", choices = character(0), selected = NULL)
    updateSelectInput(session, "single_gene_umap_trajectory", choices = character(0), selected = NULL)
    updatePickerInput(session, "gene_picker_corr", choices = character(0), selected = NULL)
    updateSelectInput(session, "select_pseudotime", choices = character(0), selected = NULL)
    updateSelectInput(session, "root_cluster_select", choices = character(0))
    shinyjs::disable("constructGraph"); shinyjs::disable("set_root_cell")
    return(invisible(NULL))
  }
  
  # Store a pseudotime calculation, select it and refresh the plots
  add_pseudotime_entry <- function(entry) {
    history <- pseudotime_history()
    history[[entry$id]] <- entry
    pseudotime_history(history)
    selected_pseudotime_id(entry$id)
    # A new root has no correlation result yet: the table of the previous root is cleared
    correlation_results(NULL); correlation_caption(NULL)
    updateSelectInput(session, "select_pseudotime", choices = pseudotime_choices(history), selected = entry$id)
  }
  
  ############################## Tab 1: load and convert ##############################
  
  output$monocle_cells_count <- renderText({
    obj <- seurat_monocle() %||% monocle_object()
    if (is.null(obj)) "—" else formatC(ncol(obj), format = "d", big.mark = ",")
  })
  
  output$monocle_clusters_names <- renderText({
    if (!is.null(seurat_monocle())) return(paste(as.character(unique(Idents(seurat_monocle()))), collapse = ", "))
    clusters <- if (!is.null(monocle_object())) get_cds_clusters(monocle_object()) else NULL
    if (is.null(clusters)) "N/A" else paste(clusters, collapse = ", ")
  })
  
  # Seurat object checked and converted to Monocle 3 (clusters kept, no re-clustering);
  # the previous trajectory is cleared only after a successful conversion
  load_seurat_for_trajectory <- function(seurat_obj, source_label) {
    temp_seurat <- prepare_seurat_for_monocle(seurat_obj)
    showModal(modalDialog(title = "Converting to Monocle3", "Converting Seurat object to Monocle format...", footer = NULL, easyClose = FALSE))
    cds <- tryCatch(seurat_to_monocle(temp_seurat), finally = removeModal())
    reset_trajectory_state()
    seurat_monocle(temp_seurat)
    monocle_object(cds)
    shinyjs::enable("constructGraph")
    showNotification(sprintf("Seurat object (%s) converted to Monocle3: %s cells, %d clusters. Ready to construct the trajectory.", source_label,
                             formatC(ncol(cds), format = "d", big.mark = ","), length(get_cds_clusters(cds))), type = "message", duration = 6)
    return(invisible(cds))
  }
  
  # Monocle object saved earlier (graph and pseudotime restored when present)
  restore_monocle_object <- function(loaded) {
    reset_trajectory_state()
    monocle_object(loaded)
    clusters <- get_cds_clusters(loaded)
    if (!is.null(clusters)) updateSelectInput(session, "root_cluster_select", choices = clusters, selected = clusters[1])
    has_graph <- has_trajectory_graph(loaded)
    has_pt <- has_valid_pseudotime(loaded)
    if (has_graph && !is.null(clusters)) {
      trajectory_plot(plot_trajectory_clusters(loaded, "Trajectory (restored from saved object)"))
      shinyjs::enable("set_root_cell")
    } else {
      shinyjs::enable("constructGraph")
    }
    if (has_pt) {
      pseudotime_plot(plot_trajectory_pseudotime(loaded, "Pseudotime (restored from saved object)"))
      entry <- new_pseudotime_entry(loaded, "unknown (restored)", NA_character_, label = "Restored pseudotime")
      add_pseudotime_entry(entry)
    }
    showNotification(paste0("Monocle object loaded: ", formatC(ncol(loaded), format = "d", big.mark = ","), " cells",
                            if (has_graph) " | trajectory graph present" else " | no graph",
                            if (has_pt) " | pseudotime present" else " | no pseudotime"), type = "message", duration = 6)
    return(invisible(loaded))
  }
  
  # One upload for every RDS: a saved Monocle object is restored, a Seurat object is converted automatically
  observeEvent(input$load_trajectory_rds, run_safely("Load RDS (trajectory)", {
    req(input$load_trajectory_rds)
    showModal(modalDialog(title = "Loading RDS file", "Reading the RDS file...", footer = NULL, easyClose = FALSE))
    loaded <- tryCatch(readRDS(input$load_trajectory_rds$datapath), finally = removeModal())
    if (inherits(loaded, "cell_data_set")) {
      restore_monocle_object(loaded)
    } else if (inherits(loaded, "Seurat")) {
      load_seurat_for_trajectory(loaded, "file")
    } else {
      stop("The file contains neither a Seurat object nor a Monocle3 cell_data_set")
    }
  }))
  
  # Seurat object already loaded in the Single or Multiple module (snapshot, no live link), converted automatically
  observeEvent(input$import_session_object_trajectory, run_safely("Import session object (trajectory)", {
    source_key <- input$trajectory_import_source
    source_label <- if (identical(source_key, "single")) "Single dataset" else "Multiple datasets"
    shared_rv <- switch(source_key, single = session$userData$single_dataset_object, multiple = session$userData$multiple_datasets_object, NULL)
    temp_seurat <- if (is.function(shared_rv)) shared_rv() else NULL
    # Checked before converting, so a failed import does not wipe an existing trajectory
    if (is.null(temp_seurat)) stop(sprintf("No Seurat object loaded in the %s module", source_label))
    load_seurat_for_trajectory(temp_seurat, paste(source_label, "module"))
  }))
  
  # Save the Monocle object (graph, pseudotime and export information)
  output$download_monocle_object <- downloadHandler(
    filename = function() paste0("monocle_object_", Sys.Date(), ".rds"),
    content = function(file) {
      req(monocle_object())
      showModal(modalDialog(title = "Saving Monocle Object", "Saving your Monocle3 object with trajectory and pseudotime...", footer = NULL, easyClose = FALSE))
      on.exit(removeModal(), add = TRUE)
      obj <- monocle_object()
      obj@metadata$cell_hub_export <- list(export_date = Sys.time(), export_module = "Trajectory - Monocle3", n_cells = ncol(obj), n_genes = nrow(obj),
                                           has_trajectory = has_trajectory_graph(obj), pseudotime_ids = names(pseudotime_history()))
      saveRDS(obj, file = file)
    }
  )
  
  
  ############################## Tab 1: graph and pseudotime ##############################
  
  # Principal graph (seed 42)
  observeEvent(input$constructGraph, run_safely("Trajectory graph", {
    req(monocle_object())
    showModal(modalDialog(title = "Constructing Graph", "Creating trajectory using existing Seurat clusters...", footer = NULL, easyClose = FALSE))
    cds <- tryCatch(learn_trajectory_graph(monocle_object(), seed = 42), finally = removeModal())
    monocle_object(cds)
    trajectory_plot(plot_trajectory_clusters(cds))
    clusters <- get_cds_clusters(cds)
    updateSelectInput(session, "root_cluster_select", choices = clusters, selected = clusters[1])
    shinyjs::enable("set_root_cell")
    showNotification("Graph construction completed! Now select a root cluster.", type = "message")
  }))
  
  # Pseudotime from the centre cell of the chosen root cluster
  observeEvent(input$set_root_cell, run_safely("Pseudotime", {
    req(monocle_object(), input$root_cluster_select)
    root <- find_root_cell(monocle_object(), input$root_cluster_select)
    showModal(modalDialog(title = "Calculating Pseudotime", paste("Computing pseudotime from cluster", input$root_cluster_select, "..."), footer = NULL, easyClose = FALSE))
    cds <- tryCatch(order_cells_from_root(monocle_object(), root$cell), finally = removeModal())
    entry <- new_pseudotime_entry(cds, input$root_cluster_select, root$cell)
    add_pseudotime_entry(entry)
    monocle_object(cds)
    trajectory_plot(plot_trajectory_clusters(cds, paste("Trajectory - Root: Cluster", input$root_cluster_select)))
    pseudotime_plot(plot_trajectory_pseudotime(cds))
    root_info_text(paste("Root cluster:", input$root_cluster_select, "\nNumber of cells:", root$n_cells, "\nRoot cell:", root$cell))
    showNotification("Pseudotime calculated successfully!", type = "message")
  }))
  
  # Switch between the pseudotime calculations of the history
  observeEvent(input$select_pseudotime, run_safely("Select pseudotime", persistent = FALSE, {
    req(input$select_pseudotime)
    entry <- pseudotime_history()[[input$select_pseudotime]]
    req(entry)
    # The picker is also updated when a new pseudotime is added (already selected and plotted): only a user switch is redrawn
    if (!identical(selected_pseudotime_id(), entry$id)) {
      monocle_object(entry$monocle_object)
      selected_pseudotime_id(entry$id)
      trajectory_plot(plot_trajectory_clusters(entry$monocle_object, paste("Trajectory -", entry$label)))
      pseudotime_plot(plot_trajectory_pseudotime(entry$monocle_object, paste("Pseudotime -", entry$label)))
      root_info_text(paste("Selected:", entry$label, "\nRoot cell:", entry$root_cell, "\nCalculated:", format(entry$timestamp, "%Y-%m-%d %H:%M")))
      # Correlation table and its download follow the selected root (NULL when not computed for this root)
      correlation_results(entry$correlation_results)
      correlation_caption(if (!is.null(entry$correlation_results)) paste("Genes correlated with", entry$label) else NULL)
      showNotification(paste("Switched to:", entry$label), type = "message", duration = 3)
    }
  }))
  
  output$trajectoryPlot <- renderPlot({ req(trajectory_plot()); trajectory_plot() })
  # Pseudotime plot with its display options, shared by the screen and the download
  pseudotime_plot_styled <- reactive({
    req(pseudotime_plot())
    style_trajectory_plot(pseudotime_plot(), hide_grid = isTRUE(input$pseudotime_hide_grid),
                          hide_axes = isTRUE(input$pseudotime_hide_axes), hide_title = isTRUE(input$pseudotime_hide_title))
  })
  output$pseudotimePlot <- renderPlot({ pseudotime_plot_styled() })  
  output$root_info <- renderText({ root_info_text() %||% "" })
  
  # Trajectory and pseudotime plots saved at their on-screen size (renderPlot draws at 72 px per inch),
  # so point, node and label sizes match the display; dpi only changes the sharpness
  save_trajectory_plot <- function(plot, file, output_id, format = "png", dpi = 300) {
    width_px <- session$clientData[[paste0("output_", output_id, "_width")]]
    height_px <- session$clientData[[paste0("output_", output_id, "_height")]]
    width_in <- if (is.numeric(width_px) && width_px > 0) width_px / 72 else 10
    height_in <- if (is.numeric(height_px) && height_px > 0) height_px / 72 else 8
    ggplot2::ggsave(file, plot = plot, device = format, dpi = dpi, width = width_in, height = height_in, units = "in", bg = "white")
    return(invisible(file))
  }
  output$download_trajectory_umap <- downloadHandler(
    filename = function() paste0("trajectory_plot_", Sys.Date(), ".", input$trajectory_download_format_main %||% "png"),
    content = function(file) {
      req(trajectory_plot())
      save_trajectory_plot(trajectory_plot(), file, "trajectoryPlot", input$trajectory_download_format_main %||% "png", input$trajectory_download_dpi %||% 300)
    }
  )
  output$download_pseudotime_umap <- downloadHandler(
    filename = function() paste0("pseudotime_plot_", Sys.Date(), ".", input$pseudotime_download_format %||% "png"),
    content = function(file) {
      save_trajectory_plot(pseudotime_plot_styled(), file, "pseudotimePlot", input$pseudotime_download_format %||% "png", input$pseudotime_download_dpi %||% 300)
    }
  )
  ############################## Tab 2: genes along the trajectory ##############################
  
  # Genes with structure along the principal graph (Moran's I), independent of the root
  observeEvent(input$run_graph_test, run_safely("Trajectory graph test", {
    req(monocle_object())
    showModal(modalDialog(title = "Analyzing Trajectory Structure", "Identifying genes with spatial patterns in the trajectory...",
                          "This is independent of root cell selection.", footer = NULL, easyClose = FALSE))
    res <- tryCatch(run_trajectory_graph_test(monocle_object(), input$sig_genes_cutoff %||% 0.05), finally = removeModal())
    if (nrow(res$significant) == 0) {
      pseudotime_diff_gene_results(NULL)
      showNotification(sprintf("No gene with q-value < %g among %d tested (minimum q-value %s).", input$sig_genes_cutoff,
                               res$n_tested, format(res$min_q, scientific = TRUE, digits = 3)), type = "warning", duration = 10)
    } else {
      pseudotime_diff_gene_results(res$significant)
      showNotification(sprintf("Found %d genes with trajectory structure (%d tested).", nrow(res$significant), res$n_tested), type = "message", duration = 5)
    }
  }))
  
  # Genes correlated with the selected pseudotime
  observeEvent(input$run_correlation_test, run_safely("Pseudotime correlation", {
    req(monocle_object())
    entry <- selected_pseudotime_entry()
    if (is.null(entry)) stop("Calculate a pseudotime first (select a root cluster)")
    showModal(modalDialog(title = "Correlation Analysis", paste("Testing genes correlated with:", entry$label),
                          "Finding up/down-regulated genes along this specific pseudotime...", footer = NULL, easyClose = FALSE))
    res <- tryCatch(run_pseudotime_correlation(monocle_object(), input$sig_genes_cutoff_corr %||% 0.05), finally = removeModal())
    if (nrow(res$significant) == 0) {
      correlation_results(NULL)
      showNotification(sprintf("No gene significantly correlated with pseudotime (%d tested).", res$n_tested), type = "warning")
    } else {
      correlation_results(res$significant)
      correlation_caption(paste("Genes correlated with", entry$label))
      # Kept with the pseudotime entry, so switching roots keeps each result
      history <- pseudotime_history()
      history[[entry$id]]$correlation_results <- res$significant
      history[[entry$id]]$correlation_timestamp <- Sys.time()
      pseudotime_history(history)
      showNotification(sprintf("Found %d genes correlated with pseudotime (%d tested).", nrow(res$significant), res$n_tested), type = "message", duration = 5)
    }
  }))
  
  output$graphTestTable <- DT::renderDT({ req(pseudotime_diff_gene_results()); graph_test_datatable(pseudotime_diff_gene_results()) })
  output$correlationTable <- DT::renderDT({ req(correlation_results()); correlation_datatable(correlation_results(), correlation_caption()) })
  
  output$download_graph_test <- createDownloadHandler(
    reactive_data = pseudotime_diff_gene_results,
    object_name_reactive = reactive({ getObjectNameForDownload(monocle_object(), default_name = "Trajectory") }),
    data_name = "graph_test_trajectory_structure",
    download_type = "csv"
  )
  output$download_correlation <- createDownloadHandler(
    reactive_data = correlation_results,
    object_name_reactive = reactive({
      tag <- pseudotime_root_tag(selected_pseudotime_entry())
      if (is.null(tag)) "Trajectory" else paste0("Trajectory_Root_", tag)
    }),
    data_name = "correlation_pseudotime",
    download_type = "csv"
  )
  
  ############################## Tab 3: gene visualization ##############################
  
  # Gene choices: graph_test genes, otherwise correlated genes (no pre-selection)
  observe(run_safely("Interface update", persistent = FALSE, {
    genes <- if (!is.null(pseudotime_diff_gene_results())) rownames(pseudotime_diff_gene_results()) else if (!is.null(correlation_results())) correlation_results()$gene else NULL
    req(genes)
    updatePickerInput(session, "multi_gene_pseudotime_trajectory", choices = genes, selected = NULL)
    updateSelectInput(session, "single_gene_umap_trajectory", choices = genes, selected = NULL)
  }))
  
  # Cluster filter of the multi-gene plot (all selected by default; a valid user selection is kept across root changes)
  observe(run_safely("Interface update", persistent = FALSE, {
    req(monocle_object())
    clusters <- get_cds_clusters(monocle_object())
    req(clusters)
    current <- isolate(input$clusters_filter_trajectory)
    selected <- if (length(current) > 0 && all(current %in% clusters)) current else clusters
    updatePickerInput(session, "clusters_filter_trajectory", choices = clusters, selected = selected)
  }))
  
  # One gene on the trajectory UMAP
  observeEvent(input$plot_single_gene_umap_trajectory, run_safely("Gene on trajectory", {
    req(monocle_object(), input$single_gene_umap_trajectory)
    single_gene_trajectory_plot(plot_gene_on_trajectory(monocle_object(), input$single_gene_umap_trajectory,
                                                        show_graph = isTRUE(input$show_trajectory_graph_single),
                                                        cell_size = input$single_gene_cell_size_trajectory %||% 0.5))
  }))
  
  # Several genes along pseudotime, restricted to the selected clusters
  observeEvent(input$plot_multi_gene_pseudotime_trajectory, run_safely("Genes along pseudotime", {
    req(monocle_object())
    genes <- input$multi_gene_pseudotime_trajectory
    if (length(genes) == 0) stop("Select at least one gene")
    if (length(genes) > 20) stop("Maximum 20 genes: reduce the selection")
    showModal(modalDialog(title = "Generating Gene Expression Plot", "Creating visualization of gene expression along pseudotime...", footer = NULL, easyClose = FALSE))
    res <- tryCatch({
      cds <- subset_cds_by_clusters(monocle_object(), input$clusters_filter_trajectory)
      plot_genes_along_pseudotime(cds, genes, min_expr = input$multi_gene_min_expr_trajectory %||% 0.1,
                                  cell_size = input$multi_gene_cell_size_trajectory %||% 0.75, trend_size = input$trend_line_size_trajectory %||% 1)
    }, finally = removeModal())
    missing <- setdiff(genes, res$genes)
    if (length(missing) > 0) showNotification(paste("Genes not found:", paste(missing, collapse = ", ")), type = "warning", duration = 5)
    gene_trajectory_plot(res$plot)
    gene_trajectory_n(length(res$genes))
  }))
  
  output$singleGeneUmapTrajectoryPlot <- renderPlot({ req(single_gene_trajectory_plot()); single_gene_trajectory_plot() })
  # Height grows with the number of gene panels (3 per row)
  output$multiGenePseudotimeTrajectoryPlot <- renderPlot({ req(gene_trajectory_plot()); gene_trajectory_plot() },
                                                         height = function() 300 + ceiling(max(1, gene_trajectory_n()) / 3) * 200)
  
  output$download_single_gene_umap_trajectory <- createDownloadHandler(
    reactive_data = single_gene_trajectory_plot,
    object_name_reactive = reactive({
      gene <- gsub("[^A-Za-z0-9_]", "_", input$single_gene_umap_trajectory %||% "gene")
      tag <- pseudotime_root_tag(selected_pseudotime_entry())
      if (is.null(tag)) paste0("SingleGene_", gene) else paste0("SingleGene_", gene, "_root_", tag)
    }),
    data_name = "trajectory_umap",
    download_type = "plot",
    plot_params = list(file_type = reactive({ input$download_format_trajectory }), width = 10, height = 8, dpi = reactive({ input$dpi_trajectory }))
  )
  output$download_multi_gene_pseudotime_trajectory <- createDownloadHandler(
    reactive_data = gene_trajectory_plot,
    object_name_reactive = reactive({
      genes <- input$multi_gene_pseudotime_trajectory
      gene_part <- if (length(genes) <= 3) paste(genes, collapse = "_") else paste0(length(genes), "_genes")
      tag <- pseudotime_root_tag(selected_pseudotime_entry())
      if (is.null(tag)) paste0("MultiGene_", gene_part) else paste0("MultiGene_", gene_part, "_root_", tag)
    }),
    data_name = "pseudotime_expression",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$download_format_trajectory }),
      width = reactive({ min(15, 5 * min(3, max(1, gene_trajectory_n()))) }),
      height = reactive({ 8 + ceiling(max(1, gene_trajectory_n()) / 3) * 3 }),
      dpi = reactive({ input$dpi_trajectory })
    )
  )
}