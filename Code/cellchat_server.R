############################## Script Cell Chat Analysis ##############################
# cellchat_server.R
cellchat_server <- function(input, output, session) {
  ############################## Loading Data ##############################
  
  # Initialize reactive values
  seurat_object_cellchat <- reactiveVal(NULL)
  db_cellchat <- reactiveVal(NULL)
  all_logs <- reactiveVal(character())  # Pour accumuler les logs
  db_loaded <- reactiveVal(FALSE)       # Pour tracker si DB déjà chargée
  overwrite_confirmed <- reactiveVal(FALSE)
  cellchat_objects <- reactiveVal(list())
  
  ############################## Download buttons state ##############################
  # Each download button stays disabled until the result it exports exists
  bind_download_state("download_analyzed_objects", function() cellchat_objects())
  bind_download_state("download_bubble_plot", function() req(bubble_plot()$plot))
  bind_download_state("download_info_table", function() req(bubble_plot()))
  bind_download_state("download_chord_plot", function() req(chord_plot()))
  bind_download_state("download_split_chord", function() req(split_chord_plot()))
  
  # Loading modal function
  loadingModal <- function(text) {
    showModal(modalDialog(
      title = "Processing...",
      p(text),
      footer = NULL,
      easyClose = FALSE
    ))
  }
  
  # Add logs function
  add_logs <- function(new_logs) {
    current_logs <- all_logs()
    all_logs(c(current_logs, new_logs))
  }
  
 
  ############################## Load Seurat Object ##############################
  # One upload for every Cell-Hub RDS: a Seurat object, a saved analysis (Seurat + CellChat objects),
  # a legacy list of CellChat objects or a single CellChat object
  observeEvent(input$seurat_file_cellchat, run_safely("seurat file cellchat", {
    logs <- character()
    
    tryCatch({
      loadingModal("Loading Seurat object...")
      logs <- c(logs, "Starting Seurat object loading process...")
      
      loaded_data <- readRDS(input$seurat_file_cellchat$datapath)
      seurat_data <- NULL
      objects <- NULL
      
      if (inherits(loaded_data, "Seurat")) {
        seurat_data <- loaded_data
      } else if (inherits(loaded_data, "CellChat")) {
        message("Detected a single CellChat object")
        objects <- stats::setNames(list(loaded_data), tools::file_path_sans_ext(input$seurat_file_cellchat$name))
      } else if (is.list(loaded_data) && !is.null(loaded_data$seurat) && !is.null(loaded_data$cellchat)) {
        message("Detected combined Seurat + CellChat structure in Seurat upload slot")
        if (!inherits(loaded_data$seurat, "Seurat")) {
          stop("Combined file does not contain a valid Seurat object under $seurat")
        }
        seurat_data <- loaded_data$seurat
        objects <- loaded_data$cellchat
      } else if (is.list(loaded_data) && length(loaded_data) > 0 &&
                 all(vapply(loaded_data, inherits, logical(1), what = "CellChat"))) {
        message("Detected legacy format (CellChat objects only)")
        objects <- loaded_data
      } else {
        stop("File contains neither a Seurat object nor CellChat objects saved by Cell-Hub")
      }
      
      if (!is.null(seurat_data)) {
        seurat_object_cellchat(seurat_data)
        logs <- c(logs, sprintf("Seurat object loaded with %d cells", ncol(seurat_data)))
        output$seurat_info_cellchat <- renderPrint({
          cat(print_seurat_info(seurat_data))
        })
      }
      
      # Saved CellChat objects are restored and the visualization selectors updated
      if (!is.null(objects)) {
        if (!is.list(objects) || length(objects) == 0) stop("No CellChat objects found in file")
        if (is.null(names(objects))) names(objects) <- paste0("CellChat_", seq_along(objects))
        n_analyzed <- sum(sapply(objects, function(obj) {
          !is.null(obj@net) && length(obj@net) > 0
        }))
        cellchat_objects(objects)
        update_cellchat_selectors(session, names(objects))
        logs <- c(logs, sprintf("%d CellChat object(s) restored (%d analyzed)", length(objects), n_analyzed))
        showNotification(
          if (!is.null(seurat_data)) {
            sprintf("Already-analyzed object detected: Seurat + %d CellChat object(s) loaded (%d analyzed). You can go directly to visualization.",
                    length(objects), n_analyzed)
          } else {
            sprintf("Loaded %d CellChat object(s) (%d analyzed)", length(objects), n_analyzed)
          },
          type = "message",
          duration = 10
        )
      }
      
      add_logs(logs)
      removeModal()
      
    }, error = function(e) {
      logs <- c(logs, sprintf("ERROR loading Seurat: %s", e$message))
      add_logs(logs)
      removeModal()
      showModal(modalDialog(
        title = "Error",
        tags$div(p(sprintf("Error loading Seurat: %s", e$message))),
        easyClose = TRUE
      ))
    })
  }))
  
  ############################## Import Seurat Object from current session ##############################
  # Snapshot of the object currently loaded in the Single or Multiple module (no live link:
  # upstream edits after import are not propagated, so existing CellChat objects stay consistent).
  # JoinLayers/droplevels are handled downstream by prepare_seurat_for_cellchat()
  observeEvent(input$import_session_object_cellchat, run_safely("import session object cellchat", {
    logs <- character()
    tryCatch({
      source_key   <- input$cellchat_import_source
      source_label <- if (identical(source_key, "single")) "Single dataset" else "Multiple datasets"
      shared_rv <- switch(source_key,
                          single   = session$userData$single_dataset_object,
                          multiple = session$userData$multiple_datasets_object,
                          NULL)
      seurat_data <- if (is.function(shared_rv)) shared_rv() else NULL
      if (is.null(seurat_data) || !inherits(seurat_data, "Seurat")) stop(sprintf("No Seurat object loaded in the %s module.", source_label))
      if (ncol(seurat_data) == 0) stop(sprintf("The Seurat object in the %s module is empty.", source_label))
      if (!"RNA" %in% names(seurat_data@assays)) stop("No RNA assay found: CellChat requires log-normalized RNA data.")
      # CellChat reads the RNA 'data' layer: a Seurat v5 object normalized only with SCTransform has none
      if (inherits(seurat_data[["RNA"]], "Assay5") && length(Layers(seurat_data[["RNA"]], search = "data")) == 0) {
        stop("RNA assay is not log-normalized (SCTransform-only objects have no RNA 'data' layer). Run LogNormalize on the RNA assay before CellChat.")
      }
      loadingModal(sprintf("Importing Seurat object from the %s module...", source_label))
      logs <- c(logs, sprintf("Importing Seurat object from %s module...", source_label))
      seurat_object_cellchat(seurat_data)
      logs <- c(logs, sprintf("Seurat object imported with %d cells", ncol(seurat_data)))
      output$seurat_info_cellchat <- renderPrint({
        cat(print_seurat_info(seurat_data))
      })
      add_logs(logs)
      removeModal()
      showNotification(sprintf("Seurat object imported from the %s module (%d cells)", source_label, ncol(seurat_data)),
                       type = "message", duration = 5)
    }, error = function(e) {
      logs <- c(logs, sprintf("ERROR importing Seurat: %s", conditionMessage(e)[1]))
      add_logs(logs)
      removeModal()
      showModal(modalDialog(
        title = "Error",
        tags$div(p(sprintf("Error importing Seurat object: %s", conditionMessage(e)[1]))),
        easyClose = TRUE
      ))
    })
  }))
  
 
  ############################## Load Database ##############################
  # Update available database checkboxes and the recommended confidence tiers when species changes
  observeEvent(input$species_cellchat, run_safely("species cellchat", {
    sp <- input$species_cellchat
    # Recommended tiers; the confirmed-roles restriction exists only for species whose database carries role_corroborated
    updateCheckboxGroupInput(session, "db_filter_tiers_cellchat", selected = cellhubdb_recommended_tiers(sp))
    updateCheckboxInput(session, "db_medium_roles_cellchat", value = cellhubdb_recommended_roles_only(sp))
    shinyjs::toggleState("db_medium_roles_cellchat", condition = sp %in% cellhubdb_species_with_roles())
    sp_dir <- file.path("databases", sp)
    if (!dir.exists(sp_dir)) {
      updateCheckboxGroupInput(session, "db_sources_cellchat", choices = NULL, selected = NULL)
      return()
    }
    # Individual source files (merged file excluded): CellHubDB_mouse_CellChatDB_v2.rds -> CellChatDB_v2
    source_files <- list.files(sp_dir, pattern = "CellHubDB_.*\\.rds$")
    source_files <- source_files[!grepl("_merged\\.rds$", source_files)]
    source_names <- sub("\\.rds$", "", sub(paste0("CellHubDB_", sp, "_"), "", source_files))
    # Checkbox labels: display name of the source and number of interactions; values stay the file labels
    choices <- setNames(source_names, sapply(source_names, function(src) {
      db <- readRDS(file.path(sp_dir, paste0("CellHubDB_", sp, "_", src, ".rds")))
      paste0(cellhubdb_source_label(src), " (", nrow(db$interaction), ")")
    }))
    updateCheckboxGroupInput(session, "db_sources_cellchat", choices = choices, selected = choices)
  }), ignoreInit = FALSE)
  
  # A loaded database no longer matches the settings once they change: it must be loaded again
  observeEvent(list(input$db_filter_tiers_cellchat, input$db_medium_roles_cellchat, input$db_mode_cellchat, input$db_sources_cellchat), run_safely("database settings cellchat", {
    if (isTRUE(db_loaded())) {
      db_loaded(FALSE)
      db_cellchat(NULL)
      showNotification("Database settings changed: click 'Load Database' again to apply them.", type = "warning", duration = 6)
    }
  }), ignoreInit = TRUE, ignoreNULL = FALSE)
  
  # Load database based on user selection, then keep the selected confidence tiers
  observeEvent(input$load_db_cellchat, run_safely("load db cellchat", {
    logs <- character()
    tryCatch({
      sp <- input$species_cellchat
      mode <- input$db_mode_cellchat
      sp_dir <- file.path("databases", sp)
      if (!dir.exists(sp_dir)) {
        stop(paste("No databases found for species:", sp))
      }
      tiers <- input$db_filter_tiers_cellchat
      if (length(tiers) == 0) stop("Select at least one confidence tier.")
      # Confirmed-roles restriction only where the database carries the flag (human, mouse)
      roles_only <- isTRUE(input$db_medium_roles_cellchat) && sp %in% cellhubdb_species_with_roles()
      loadingModal("Loading CellHub Database...")
      merged_path <- file.path(sp_dir, paste0("CellHubDB_", sp, "_merged.rds"))
      if (mode == "merged") {
        # Load pre-built merged database
        if (!file.exists(merged_path)) {
          stop(paste("Merged database not found:", merged_path))
        }
        CellHubDB <- readRDS(merged_path)
        logs <- c(logs, paste("Loaded merged CellHub Database for", sp))
        all_sources <- unique(unlist(strsplit(CellHubDB$interaction$sources, ", ")))
        logs <- c(logs, paste("Sources:", paste(cellhubdb_source_label(sort(all_sources)), collapse = ", ")))
        CellHubDB$interaction <- annotate_cellhubdb_tier(CellHubDB$interaction)
      } else {
        # Load and merge selected individual sources
        selected <- input$db_sources_cellchat
        if (is.null(selected) || length(selected) == 0) {
          stop("No databases selected. Please select at least one database.")
        }
        db_list <- setNames(lapply(selected, function(src) {
          readRDS(file.path(sp_dir, paste0("CellHubDB_", sp, "_", src, ".rds")))
        }), selected)
        for (src in selected) logs <- c(logs, paste("  +", src, ":", nrow(db_list[[src]]$interaction), "interactions"))
        if (length(selected) == 1) {
          # Single source: the per-source file is already a complete CellChat database
          CellHubDB <- db_list[[1]]
          logs <- c(logs, paste("Loaded single database:", selected))
        } else {
          # Multiple sources: same merge rules as the pre-built merged file (curate_CellHub_DB.R)
          logs <- c(logs, paste("Merging", length(selected), "databases..."))
          merge_result <- merge_cellhubdb_sources(db_list)
          CellHubDB <- merge_result$db
          logs <- c(logs, sprintf("Merged: %d -> %d interactions (%d supported by 2+ sources, %d non-oriented pairs aligned)",
                                  merge_result$stats$n_input, merge_result$stats$n_merged,
                                  merge_result$stats$n_multi, merge_result$stats$n_flipped))
        }
        # Per-source files carry no tier: curated sources and confirmed roles joined from the merged release file,
        # tier recomputed on the selected sources
        merged_interaction <- if (file.exists(merged_path)) readRDS(merged_path)$interaction else NULL
        if (!is.null(merged_interaction)) {
          n_unmatched <- sum(is.na(match_interaction_rows(CellHubDB$interaction, merged_interaction)))
          if (n_unmatched > 0) logs <- c(logs, sprintf("%d interactions not found in the merged release file: no confidence tier", n_unmatched))
        }
        CellHubDB$interaction <- annotate_cellhubdb_tier(CellHubDB$interaction, merged_interaction, selected_sources = selected)
        logs <- c(logs, "Confidence tiers recomputed on the selected databases")
      }
      # Confidence tiers kept
      tier_counts <- cellhubdb_tier_counts(CellHubDB$interaction)
      if (!is.null(tier_counts)) logs <- c(logs, paste("Tiers:", paste(names(tier_counts), tier_counts, collapse = " | ")))
      filtered <- filter_cellhubdb_tiers(CellHubDB, tiers, medium_roles_only = roles_only)
      CellHubDB <- filtered$db
      tier_label <- paste0(paste(filtered$tiers, collapse = " + "), if (filtered$medium_roles_only) ", medium with confirmed roles only" else "")
      logs <- c(logs, sprintf("Tiers kept (%s): %d -> %d interactions", tier_label, filtered$n_before, filtered$n_after))
      # Schema, gene reference and pathway labels checked for every release and source combination before CellChat uses them
      prepared <- prepare_cellchat_db(CellHubDB)
      CellHubDB <- prepared$db
      logs <- c(logs, prepared$messages)
      # Tier selection travels with the database and the CellChat objects built on it
      attr(CellHubDB, "cellhub_tiers") <- tier_label
      db_cellchat(CellHubDB)
      db_loaded(TRUE)
      logs <- c(logs, paste("Database contains", nrow(CellHubDB$interaction), "interactions and", nrow(CellHubDB$complex), "complexes"))
      output$db_info_cellchat <- renderPrint({
        cat(print_db_info(CellHubDB))
      })
      add_logs(logs)
      removeModal()
    }, error = function(e) {
      logs <- c(logs, paste("ERROR loading database:", conditionMessage(e)[1]))
      add_logs(logs)
      removeModal()
      showModal(modalDialog(
        title = "Error",
        tags$div(
          p(paste("Error loading CellHub Database:", conditionMessage(e)[1]))
        ),
        easyClose = TRUE
      ))
    })
  }))
  
  # Update the accumulated logs display
  output$loading_logs <- renderPrint({
    cat(paste(all_logs(), collapse = "\n"))
  })
  
  ############################## Genes of interest (named gene lists) ##############################
  gene_lists_cellchat <- reactiveVal(list())
  gene_query_result <- reactiveVal(NULL)
  # A newly loaded (or unloaded) database invalidates the previous query
  observeEvent(db_cellchat(), gene_query_result(NULL), ignoreNULL = FALSE)
  
  # Genes of the list chosen in a plot selector ("" = no list)
  selected_gene_list <- function(name) {
    if (is.null(name) || !nzchar(name)) return(character(0))
    genes <- gene_lists_cellchat()[[name]]
    if (is.null(genes)) stop("Gene list not found: ", name)
    genes
  }
  
  # Plot selectors follow the saved lists and keep their choice while it still exists
  observe(run_safely("Gene list selectors", persistent = FALSE, {
    choices <- gene_list_choices(gene_lists_cellchat())
    for (id in c("bubble_gene_list_cellchat", "chord_gene_list_cellchat")) {
      current <- isolate(input[[id]])
      updateSelectInput(session, id, choices = choices, selected = if (!is.null(current) && current %in% choices) current else "")
    }
  }))
  
  # Add or replace a named gene list; the fields are kept so that the list can be edited and added again under the same name
  observeEvent(input$add_gene_list_cellchat, run_safely("add gene list", {
    if (is.null(input$gene_list_name_cellchat) || !nzchar(trimws(input$gene_list_name_cellchat))) {
      stop("Type a name in 'List name' before clicking 'Add list' (an existing name replaces that list)")
    }
    res <- add_gene_list(gene_lists_cellchat(), input$gene_list_name_cellchat, input$gene_list_text_cellchat)
    gene_lists_cellchat(res$lists)
    updateSelectInput(session, "gene_list_selected_cellchat", choices = gene_list_choices(res$lists, none_label = NULL), selected = res$name)
    showNotification(sprintf("%s gene list '%s': %d genes", if (res$replaced) "Replaced" else "Added", res$name, length(res$genes)), type = "message")
  }))
  
  observeEvent(input$delete_gene_list_cellchat, run_safely("delete gene list", {
    name <- input$gene_list_selected_cellchat
    if (is.null(name) || !nzchar(name)) stop("Select a list to delete")
    lists <- gene_lists_cellchat()
    lists[[name]] <- NULL
    gene_lists_cellchat(lists)
    updateSelectInput(session, "gene_list_selected_cellchat", choices = gene_list_choices(lists, none_label = NULL))
    showNotification(sprintf("Gene list '%s' deleted", name), type = "message")
  }))
  
  output$gene_lists_table_cellchat <- renderTable({
    lists <- gene_lists_cellchat()
    req(length(lists) > 0)
    gene_lists_table(lists)
  }, striped = TRUE, spacing = "xs", width = "100%")
  
  observeEvent(input$run_gene_query_cellchat, run_safely("gene query cellchat", {
    db <- db_cellchat()
    if (is.null(db) || !isTRUE(db_loaded())) stop("Load a database first: the query searches the loaded database, with its confidence tiers")
    name <- input$gene_list_selected_cellchat
    if (is.null(name) || !nzchar(name)) stop("Add a gene list first, then select it")
    res <- query_cellhubdb_genes(db, selected_gene_list(name))
    res$list_name <- name
    gene_query_result(res)
    msg <- sprintf("%s: %d of %d genes found in %d of %d interactions", name, length(res$genes_found), res$n_genes, nrow(res$table), res$n_db)
    if (length(res$genes_absent) > 0) msg <- paste0(msg, ". Not in the loaded database: ", paste(head(res$genes_absent, 20), collapse = ", "), if (length(res$genes_absent) > 20) ", ..." else "")
    showNotification(msg, type = if (nrow(res$table) == 0) "warning" else "message", duration = 8)
  }))
  
  output$gene_query_summary_cellchat <- renderText({
    res <- gene_query_result(); req(res)
    sprintf("%s: %d of %d genes found | %d interactions of the loaded database (%d in total)%s", res$list_name, length(res$genes_found), res$n_genes, nrow(res$table), res$n_db,
            if (length(res$genes_absent) > 0) paste0(" | not found: ", paste(head(res$genes_absent, 20), collapse = ", "), if (length(res$genes_absent) > 20) ", ..." else "") else "")
  })
  
  # Client-side table so that the CSV / Excel buttons export every row, not only the visible page
  output$gene_query_table_cellchat <- DT::renderDataTable({
    res <- gene_query_result(); req(res, nrow(res$table) > 0)
    file_stub <- paste0("CellHubDB_gene_query_", gsub("[^A-Za-z0-9]+", "_", res$list_name), "_", input$species_cellchat, "_", format(Sys.Date(), "%Y%m%d"))
    DT::datatable(res$table, rownames = FALSE, selection = "none", extensions = "Buttons",
                  options = list(pageLength = 25, dom = "Bfrtip", scrollX = TRUE,
                                 buttons = list(list(extend = "csv", filename = file_stub, text = "CSV"),
                                                list(extend = "excel", filename = file_stub, text = "Excel"))))
  }, server = FALSE)
  
  # Helper functions
  print_seurat_info <- function(seurat_obj) {
    # Use proper methods for Seurat object
    tryCatch({
      info <- list(
        ncells = ncol(seurat_obj),
        nfeatures = nrow(seurat_obj),
        assays = names(seurat_obj@assays)
      )
      
      return(sprintf(
        "Dataset loaded:\nNumber of cells: %d\nNumber of features: %d\nAvailable assays: %s",
        info$ncells,
        info$nfeatures,
        paste(info$assays, collapse = ", ")
      ))
    }, error = function(e) {
      return("Error reading Seurat object properties")
    })
  }
  
  print_db_info <- function(db) {
    inter <- db$interaction
    n_interactions <- nrow(inter)
    n_pathways <- length(unique(inter$pathway_name[!is.na(inter$pathway_name)]))
    # Individual source names (the sources column lists every source of a pair, comma-separated), shown with their display names
    source_names <- if ("sources" %in% colnames(inter)) sort(unique(unlist(strsplit(inter$sources[!is.na(inter$sources)], "\\s*,\\s*")))) else character(0)
    source_labels <- cellhubdb_source_label(source_names)
    # Count by interaction type
    type_counts <- table(inter$interaction_type, useNA = "no")
    type_str <- paste(names(type_counts), type_counts, sep = ": ", collapse = " | ")
    # Count multi-source entries
    n_multi <- if ("sources" %in% colnames(inter)) sum(grepl(",", inter$sources), na.rm = TRUE) else 0
    tier_counts <- cellhubdb_tier_counts(inter)
    paste0(
      "CellHub Database loaded:",
      sprintf("\n  Interactions: %d", n_interactions),
      sprintf("\n  Sources: %d (%s)", length(source_labels), paste(source_labels, collapse = ", ")),
      sprintf("\n  Pathways: %d", n_pathways),
      if (nchar(type_str) > 0) sprintf("\n  Types: %s", type_str) else "",
      if (n_multi > 0) sprintf("\n  Cross-validated (multi-source): %d (%.0f%%)", n_multi, 100 * n_multi / n_interactions) else "",
      if (!is.null(tier_counts)) sprintf("\n  Confidence tiers: %s", paste(names(tier_counts), tier_counts, collapse = " | ")) else ""
    )
  }
  
  ###################### TAB 2 ###################
  
  # Update column choices when Seurat object is loaded
  observe({
    req(seurat_object_cellchat())
    
    tryCatch({
      seurat_obj <- seurat_object_cellchat()
      
      # Safety check
      if (ncol(seurat_obj) == 0) {
        showNotification("Loaded Seurat object is empty!", type = "error")
        return()
      }
      
      # CRITICAL: Ensure ClusterIdents column exists with proper names
      if (!"ClusterIdents" %in% colnames(seurat_obj@meta.data)) {
        message("ClusterIdents column not found - creating it from current Idents")
        seurat_obj$ClusterIdents <- Idents(seurat_obj)
        seurat_object_cellchat(seurat_obj)
      }
      
      # Define patterns to exclude from column choices
      excluded_patterns <- c(
        "percent.mt", "nCount_ATAC", "nFeature_ATAC", "nFeature_RNA", "nCount_RNA",
        "^RNA_snn_res", "^RNA_snn", "^pANN", "^DF", "^integrated", "^integrated_snn"
      )
      
      pattern <- paste(excluded_patterns, collapse = "|")
      meta_cols <- colnames(seurat_obj@meta.data)
      filtered_cols <- meta_cols[!grepl(pattern, meta_cols)]
      
      # SMART DEFAULT SELECTION for primary column
      preferred_columns <- c("ClusterIdents", "cluster_names", "seurat_clusters", "ident")
      default_column <- NULL
      
      for (col in preferred_columns) {
        if (col %in% filtered_cols) {
          default_column <- col
          message(paste("Selected default grouping column:", col))
          break
        }
      }
      
      # If none of preferred columns exist, find first categorical column
      if (is.null(default_column)) {
        for (col in filtered_cols) {
          col_data <- seurat_obj@meta.data[[col]]
          n_unique <- length(unique(col_data))
          if (is.factor(col_data) || is.character(col_data) || n_unique < 50) {
            default_column <- col
            message(paste("No preferred column found, using first categorical column:", col))
            break
          }
        }
      }
      
      # Final fallback
      if (is.null(default_column) && length(filtered_cols) > 0) {
        default_column <- filtered_cols[1]
        message(paste("Using first available column as fallback:", default_column))
      }
      
      # Update all column selectors
      updateSelectInput(session, "group_by_cellchat", choices = filtered_cols, selected = default_column)
      updateSelectInput(session, "combine_with_column", choices = c("None" = "", filtered_cols), selected = "")
      updateSelectInput(session, "combine_with_column_2", choices = c("None" = "", filtered_cols), selected = "")
      updateSelectInput(session, "condition_filter_column", choices = c("None" = "", filtered_cols), selected = "")

    }, error = function(e) {
      showNotification(paste("Error updating column choices:", conditionMessage(e)[1]), type = "warning", duration = 10)
    })
  })
  
  
  # Preview of groups that will be created based on selected columns (same labels as create_cellchat_object_multicol)
  output$cellchat_groups_preview <- renderUI({
    req(seurat_object_cellchat(), input$group_by_cellchat)
    tryCatch({
      seurat_obj <- seurat_object_cellchat()
      selected_cols <- c(input$group_by_cellchat)
      if (!is.null(input$combine_with_column) && input$combine_with_column != "") selected_cols <- c(selected_cols, input$combine_with_column)
      if (!is.null(input$combine_with_column_2) && input$combine_with_column_2 != "") selected_cols <- c(selected_cols, input$combine_with_column_2)
      if (length(selected_cols) != length(unique(selected_cols))) {
        return(tags$p("⚠️ Duplicate columns selected", style = "color: #ff6b6b; font-size: 12px; margin: 0;"))
      }
      labels <- build_cellchat_group_labels(seurat_obj@meta.data, selected_cols)
      groups <- sort(unique(labels[!is.na(labels)]))
      n_groups <- length(groups)
      n_missing <- sum(is.na(labels))
      if (n_groups == 0) return(tags$p("No valid groups found", style = "color: rgba(255,255,255,0.6); font-size: 11px; margin: 0;"))
      n_show <- if (length(selected_cols) == 1) 8 else 6
      tags$div(
        tags$p(tags$strong(paste(n_groups, if (length(selected_cols) == 1) "groups" else "combined groups")),
               if (length(selected_cols) == 1) " will be created from column: " else " from: ",
               tags$code(paste(selected_cols, collapse = " + "), style = "background-color: rgba(255,255,255,0.2); padding: 2px 5px; border-radius: 3px;"),
               style = "margin: 0 0 8px 0; color: white; font-size: 12px;"),
        tags$p(paste(head(groups, n_show), collapse = ", "),
               if (n_groups > n_show) tags$span(paste0(" ... (", n_groups - n_show, " more)"), style = "color: rgba(255,255,255,0.6);") else NULL,
               style = "margin: 0; color: rgba(255,255,255,0.85); font-size: 11px; font-family: monospace;"),
        # Cells with a missing value in any selected column are excluded from the analysis
        if (n_missing > 0) tags$p(sprintf("%d cells with a missing value will be excluded", n_missing),
                                  style = "margin: 6px 0 0 0; color: #ffcc80; font-size: 11px;") else NULL
      )
    }, error = function(e) {
      tags$p("Unable to generate preview", style = "color: rgba(255,255,255,0.5); font-size: 11px; margin: 0;")
    })
  })
  
  # UI to select specific groups to keep (same labels as create_cellchat_object_multicol)
  output$cellchat_groups_filter_ui <- renderUI({
    req(seurat_object_cellchat(), input$group_by_cellchat)
    tryCatch({
      seurat_obj <- seurat_object_cellchat()
      selected_cols <- c(input$group_by_cellchat)
      if (!is.null(input$combine_with_column) && input$combine_with_column != "") selected_cols <- c(selected_cols, input$combine_with_column)
      if (!is.null(input$combine_with_column_2) && input$combine_with_column_2 != "") selected_cols <- c(selected_cols, input$combine_with_column_2)
      selected_cols <- unique(selected_cols)
      labels <- build_cellchat_group_labels(seurat_obj@meta.data, selected_cols)
      groups <- sort(unique(labels[!is.na(labels)]))
      if (length(groups) == 0) return(tags$p("No groups available", style = "color: rgba(255,255,255,0.6); font-size: 11px;"))
      selectizeInput("selected_groups_filter", label = NULL, choices = groups, selected = NULL, multiple = TRUE,
                     options = list(plugins = list('remove_button'), placeholder = 'Select groups (empty = keep all)', maxItems = 20))
    }, error = function(e) NULL)
  })
  
  
  
  
 
  output$condition_values_ui <- renderUI({
    req(seurat_object_cellchat())
    
    # If no column selected, show placeholder
    if (is.null(input$condition_filter_column) || input$condition_filter_column == "") {
      return(tags$p("Select a column above to see available values", 
                    style = "color: rgba(255,255,255,0.6); font-size: 12px; font-style: italic;"))
    }
    
    tryCatch({
      seurat_obj <- seurat_object_cellchat()
      selected_column <- input$condition_filter_column
      
      # Check if column exists
      if (!selected_column %in% colnames(seurat_obj@meta.data)) {
        return(tags$p("Column not found in metadata", 
                      style = "color: red; font-size: 12px;"))
      }
      
      # Get unique values from selected column
      available_values <- unique(seurat_obj@meta.data[[selected_column]])
      available_values <- available_values[!is.na(available_values)]
      available_values <- sort(as.character(available_values))
      
      if (length(available_values) == 0) {
        return(tags$p("No values found in this column", 
                      style = "color: red; font-size: 12px;"))
      }
      
      # Count cells per value
      value_counts <- table(seurat_obj@meta.data[[selected_column]])
      
      # CORRECT WAY: values first, then names
      condition_labels <- available_values
      names(condition_labels) <- sapply(available_values, function(v) {
        sprintf("%s (%d cells)", v, value_counts[as.character(v)])
      })
      
      checkboxGroupInput(
        inputId = "selected_conditions_cellchat",
        label = NULL,
        choices = condition_labels,
        selected = available_values
      )
      
    }, error = function(e) {
      return(tags$p(paste("Error:", e$message), style = "color: red; font-size: 12px;"))
    })
  })

 
 

  # Create and analyze CellChat object
  observeEvent(input$create_and_analyze_cellchat, {
    tryCatch({
      # No dataset loaded (or the import failed): stop before any other check, with a message the user can act on
      if (is.null(seurat_object_cellchat())) {
        showNotification("No dataset loaded: import a Seurat object in the 'Load data' tab first.", type = "error", duration = 10)
        return()
      }
      
      # Collect selected columns for grouping
      grouping_columns <- c(input$group_by_cellchat)
      
      if (!is.null(input$combine_with_column) && input$combine_with_column != "") {
        grouping_columns <- c(grouping_columns, input$combine_with_column)
      }
      
      if (!is.null(input$combine_with_column_2) && input$combine_with_column_2 != "") {
        grouping_columns <- c(grouping_columns, input$combine_with_column_2)
      }
      
      if (length(grouping_columns) != length(unique(grouping_columns))) {
        showNotification("Cannot use the same column multiple times!", type = "error")
        return()
      }
      
      # Condition filter: the checkbox group rendered by condition_values_ui is selected_conditions_cellchat
      condition_column <- if (!is.null(input$condition_filter_column) && input$condition_filter_column != "") input$condition_filter_column else NULL
      selected_conditions <- input$selected_conditions_cellchat
      # An active filter column with every value unchecked would silently analyze all cells — block it instead
      if (!is.null(condition_column) && length(selected_conditions) == 0) {
        showNotification(paste0("No value selected for filter column '", condition_column, "'. Select at least one value or set the filter column to None."),
                         type = "error", duration = 10)
        return()
      }
      
      obj_name <- trimws(input$cellchat_object_name)
      if (nchar(obj_name) == 0) {
        showNotification("Please provide an object name", type = "error")
        return()
      }
      
      current_objects <- cellchat_objects()
      if (obj_name %in% names(current_objects)) {
        if (isTRUE(overwrite_confirmed())) {
          overwrite_confirmed(FALSE)
        } else {
          showModal(modalDialog(
            title = "Object Already Exists",
            paste("An object named", obj_name, "already exists. Do you want to overwrite it?"),
            footer = tagList(
              modalButton("Cancel"),
              actionButton("confirm_overwrite", "Overwrite", class = "btn-danger")
            )
          ))
          return()
        }
      }
      
      # The database is checked before the (long) object creation, not after it
      if (is.null(db_cellchat()) || !db_loaded()) {
        showNotification("Please load the database first using 'Load Database' button", type = "error", duration = 10)
        return()
      }
      showModal(modalDialog(
        title = "Creating CellChat Object",
        tags$div(
          tags$h5("Please wait while CellChat object is being created and analyzed..."),
          tags$div(id = "cellchat_log",
                   style = "max-height: 300px; overflow-y: auto; background-color: #f5f5f5; padding: 10px; border-radius: 4px; font-family: monospace; font-size: 12px;")
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      
      add_analysis_log <- function(msg) {
        formatted_msg <- paste0(format(Sys.time(), "%H:%M:%S"), " - ", msg)
        add_logs(formatted_msg)
        current_logs <- all_logs()
        display_log <- paste(tail(current_logs, 30), collapse = "<br>")
        shinyjs::html("cellchat_log", display_log)
        shinyjs::runjs("document.getElementById('cellchat_log').scrollTop = document.getElementById('cellchat_log').scrollHeight;")
      }
      
      add_analysis_log("Starting CellChat object creation...")
      add_analysis_log(paste("Grouping by:", paste(grouping_columns, collapse = " + ")))
      
      result <- create_cellchat_object_multicol(
        seurat_obj          = seurat_object_cellchat(),
        grouping_columns    = grouping_columns,
        selected_groups     = input$selected_groups_filter,
        condition_column    = condition_column,
        selected_conditions = selected_conditions,
        log_function        = add_analysis_log
      )
      
      if (!result$success) {
        removeModal()
        showNotification(paste("Failed to create CellChat object:", result$message), type = "error", duration = 10)
        return()
      }
      
      cellchat_obj <- result$object
      add_analysis_log("CellChat object created successfully!")
      
      add_analysis_log("Assigning database to CellChat object...")
      
      cellchat_obj@DB <- db_cellchat()
      add_analysis_log("CellHub Database assigned successfully")
      
      add_analysis_log("Starting analysis pipeline...")
      
      analysis_result <- run_cellchat_pipeline(
        cellchat_obj = cellchat_obj,
        object_name  = obj_name,
        log_function = add_analysis_log,
        min_cells    = 10
      )
      
      if (!analysis_result$success) {
        removeModal()
        showNotification(paste("Analysis failed:", analysis_result$message), type = "error", duration = 10)
        return()
      }
      
      current_objects <- cellchat_objects()
      current_objects[[obj_name]] <- analysis_result$object
      cellchat_objects(current_objects)
      
      update_cellchat_selectors(session, names(current_objects))
      
      add_analysis_log(paste("Analysis complete! Found", analysis_result$n_interactions, "significant interactions"))
      
      removeModal()
      showNotification(
        paste("CellChat object", obj_name, "created and analyzed successfully!"),
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error:", conditionMessage(e)[1]), type = "error", duration = 10)
      message("Full error: ", conditionMessage(e))
    })
  })
  
  # Display analysis logs (accumulated in all_logs by the creation/analysis observer)
  output$analysis_logs_cellchat <- renderPrint({
    cat(paste(all_logs(), collapse = "\n"))
  })
  
  observeEvent(input$confirm_overwrite, run_safely("confirm overwrite", {
    removeModal()
    overwrite_confirmed(TRUE)
    shinyjs::click("create_and_analyze_cellchat")
  }))
  
  
 
  # Fill the object selectors of the visualization tabs
  observe(run_safely("Interface update", persistent = FALSE, {
    req(cellchat_objects())
    req(length(cellchat_objects()) > 0)
    object_choices <- names(cellchat_objects())
    updateSelectInput(session, "subset_choose_cellchat", choices = object_choices)
    updateSelectInput(session, "cellchat_obj_chord", choices = object_choices)
  }))
  
  # Source and target selectors of the bubble plot, existing selections kept when still valid
  observe(run_safely("Interface update", persistent = FALSE, {
    req(input$subset_choose_cellchat)
    cellchat_obj <- cellchat_objects()[[input$subset_choose_cellchat]]
    req(cellchat_obj, !is.null(cellchat_obj@idents))
    cell_types <- levels(cellchat_obj@idents)
    selected_sources <- isolate(input$sources_use_cellchat)
    selected_targets <- isolate(input$targets_use_cellchat)
    updateSelectizeInput(session, "sources_use_cellchat", choices = cell_types, selected = intersect(selected_sources, cell_types),
                         options = list(plugins = list('remove_button')))
    updateSelectizeInput(session, "targets_use_cellchat", choices = cell_types, selected = intersect(selected_targets, cell_types),
                         options = list(plugins = list('remove_button')))
  }))
  
  
  output$download_analyzed_objects <- downloadHandler(
    filename = function() paste0("CellHub_Seurat_CellChat_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"),
    content = function(file) {
      req(cellchat_objects())
      showModal(modalDialog(title = "Saving Objects", "Saving Seurat + CellChat objects...", easyClose = FALSE, footer = NULL))
      # The modal is closed whatever happens during saving
      on.exit(removeModal(), add = TRUE)
      combined_data <- list(seurat = seurat_object_cellchat(), cellchat = cellchat_objects())
      saveRDS(combined_data, file)
      showNotification("Seurat + CellChat objects saved successfully!", type = "message")
    }
  )
  
  # Reactive value for plot
  bubble_plot <- reactiveVal()
  
  # Generate bubble plot (optionally restricted to a gene list, applied after inference)
  observeEvent(input$generate_plot_cellchat, run_safely("generate plot cellchat", {
    req(input$subset_choose_cellchat)
    
    showNotification("Generating bubble plot...", type = "message", duration = 2)
    
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$subset_choose_cellchat]]
      genes <- selected_gene_list(input$bubble_gene_list_cellchat)
      
      # Generate plot using helper
      result <- generate_bubble_plot(
        cellchat_obj = cellchat_obj,
        sources = input$sources_use_cellchat,
        targets = input$targets_use_cellchat,
        threshold = input$threshold_cellchat,
        flip_axes = input$flip_axes,
        exclude_intra = input$exclude_intra,
        signaling = NULL,  # character vector of pathway names, e.g. c("MHC-II", "CXCL")
        genes = genes
      )
      # Get interaction type info
      interaction_info <- get_interaction_type_info(cellchat_obj, result$communications)
      
      # Store plot and info
      bubble_plot(list(
        plot = result$plot,
        interaction_info = interaction_info
      ))
      
      # Genes of the list absent from the significant interactions
      if (length(result$genes_absent) > 0) {
        showNotification(paste0(length(result$genes_absent), " of ", length(genes), " genes in no significant interaction: ",
                                paste(head(result$genes_absent, 20), collapse = ", "), if (length(result$genes_absent) > 20) ", ..." else ""),
                         type = "warning", duration = 8)
      }
      
      # Notification
      if (!is.null(interaction_info) && "interaction_type" %in% colnames(interaction_info)) {
        type_counts <- table(interaction_info$interaction_type)
        type_summary <- paste(names(type_counts), type_counts, sep = ": ", collapse = ", ")
        showNotification(
          paste("Generated", result$n_interactions, "L-R pairs. By type:", type_summary),
          type = "message",
          duration = 5
        )
      } else {
        showNotification(
          paste("Generated", result$n_interactions, "ligand-receptor pairs"),
          type = "message"
        )
      }
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  }))
  
  # Render the plot without infinite refresh
  output$bubble_plot_cellchat <- renderPlot({
    req(bubble_plot())
    
    # Get the plot
    plot_obj <- bubble_plot()$plot
    
    # Return the plot
    plot_obj
  }, width = function() {
    # Use a default value if nothing is specified
    if(is.null(input$plot_width) || is.na(as.numeric(input$plot_width))) {
      return(1000)  # Default value
    }
    return(as.numeric(input$plot_width))
  }, height = function() {
    # Use a default value if nothing is specified
    if(is.null(input$plot_height) || is.na(as.numeric(input$plot_height))) {
      return(800)  # Default value
    }
    return(as.numeric(input$plot_height))
  })
  
  # Download bubble plot - UTILISE LE HELPER EXISTANT
  output$download_bubble_plot <- createDownloadHandler(
    reactive_data = reactive({ bubble_plot()$plot }),  # Juste le plot
    object_name_reactive = reactive({ input$subset_choose_cellchat }),
    data_name = reactive({ 
      paste0("BubblePlot_", 
             paste(head(input$sources_use_cellchat, 2), collapse = "_"),
             "_to_",
             paste(head(input$targets_use_cellchat, 2), collapse = "_"))
    }),
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$bubble_plot_format }),
      width = reactive({ input$plot_width / 72 }),
      height = reactive({ input$plot_height / 72 }),
      dpi = 300
    )
  )
  
  # Download interaction info table - UTILISE LE HELPER
  output$download_info_table <- createDownloadHandler(
    reactive_data = reactive({ 
      info <- bubble_plot()$interaction_info
      if (is.null(info) || nrow(info) == 0) {
        return(data.frame(Message = "No interaction type information available"))
      }
      return(info)
    }),
    object_name_reactive = reactive({ input$subset_choose_cellchat }),
    data_name = "InteractionTypes",
    download_type = "csv"
  )
  
 
  
  # Reactive value pour stocker le plot
  chord_plot <- reactiveVal()
  
  # Update chord plot cell-type selectors from the selected object
  observe({
    req(input$cellchat_obj_chord)
    cellchat_obj <- cellchat_objects()[[input$cellchat_obj_chord]]
    req(cellchat_obj, !is.null(cellchat_obj@idents))
    cell_types <- levels(cellchat_obj@idents)
    updateSelectizeInput(session, "sender_groups", choices = cell_types, selected = NULL)
    updateSelectizeInput(session, "receiver_groups", choices = cell_types, selected = NULL)
  })
  
  # Dynamic color inputs for chord plot — initialized from stored colors if available
  output$color_inputs_chord <- renderUI({
    tryCatch({
      req(input$cellchat_obj_chord)
      req(input$use_custom_colors)
      req(input$sender_groups)
      req(input$receiver_groups)
      
      all_groups   <- unique(c(input$sender_groups, input$receiver_groups))
      req(length(all_groups) > 0)
      
      cellchat_obj  <- cellchat_objects()[[input$cellchat_obj_chord]]
      stored_colors <- find_stored_colors(cellchat_obj, cluster_names = all_groups)
      
      n <- length(all_groups)
      fallback_pal <- if (n <= 8) {
        RColorBrewer::brewer.pal(max(3, n), "Set1")[seq_len(n)]
      } else {
        scales::hue_pal()(n)
      }
      names(fallback_pal) <- all_groups
      
      n_from_obj  <- if (!is.null(stored_colors)) sum(all_groups %in% names(stored_colors)) else 0
      info_banner <- if (n_from_obj == n) {
        tags$p(icon("circle-check"),
               sprintf(" Colors loaded from object (%d/%d)", n_from_obj, n),
               style = "color: #2e7d32; font-size: 12px; margin-bottom: 6px;")
      } else if (n_from_obj > 0) {
        tags$p(icon("triangle-exclamation"),
               sprintf(" Partial match: %d/%d colors from object", n_from_obj, n),
               style = "color: #e65100; font-size: 12px; margin-bottom: 6px;")
      } else {
        tags$p(icon("circle-info"),
               " No stored colors found — using default palette.",
               style = "color: #888; font-size: 12px; margin-bottom: 6px;")
      }
      
      color_inputs <- lapply(seq_along(all_groups), function(i) {
        group_name <- all_groups[i]
        init_color <- if (!is.null(stored_colors) && group_name %in% names(stored_colors)) {
          stored_colors[[group_name]]
        } else {
          fallback_pal[[group_name]]
        }
        colourInput(
          inputId = paste0("color_group_", i),
          label   = paste("Color for", group_name),
          value   = init_color
        )
      })
      
      div(class = "color-inputs", info_banner, color_inputs)
      
    }, error = function(e) {
      print(paste("Error in color inputs generation:", e$message))
      return(NULL)
    })
  })
  
  observeEvent(input$generate_chord_plot, run_safely("generate chord plot", {
    req(input$cellchat_obj_chord)
    req(length(input$sender_groups) > 0)
    req(length(input$receiver_groups) > 0)
    
    showNotification("Generating chord plot...", type = "message", duration = 2)
    
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$cellchat_obj_chord]]
      genes <- selected_gene_list(input$chord_gene_list_cellchat)
      
      custom_colors <- NULL
      if (isTRUE(input$use_custom_colors)) {
        all_groups <- unique(c(input$sender_groups, input$receiver_groups))
        n          <- length(all_groups)
        
        # Same fallback palette as renderUI to stay consistent
        fallback_pal <- if (n <= 8) {
          RColorBrewer::brewer.pal(max(3, n), "Set1")[seq_len(n)]
        } else {
          scales::hue_pal()(n)
        }
        
        custom_colors <- sapply(seq_along(all_groups), function(i) {
          color_input <- input[[paste0("color_group_", i)]]
          if (is.null(color_input)) fallback_pal[i] else color_input
        })
        names(custom_colors) <- all_groups
      }
      
      # Render at >= 1200 px so the image stays sharp when the display size is increased
      png_size <- max(1200L, as.integer(input$chord_plot_dims %||% 800))
      
      result <- generate_chord_plot(
        cellchat_obj  = cellchat_obj,
        senders       = input$sender_groups,
        receivers     = input$receiver_groups,
        threshold     = input$prob_threshold_chord,
        top_n         = input$chord_top_n %||% 0,
        custom_colors = custom_colors,
        label_size    = input$chord_label_size,
        plot_size     = png_size,
        genes         = genes
      )
      
      # Delete the previous temporary PNG before replacing it
      previous <- chord_plot()
      if (!is.null(previous) && !is.null(previous$output_path) && file.exists(previous$output_path)) unlink(previous$output_path)
      chord_plot(result)
      # Genes of the list absent from the significant interactions
      if (length(result$genes_absent) > 0) {
        showNotification(paste0(length(result$genes_absent), " of ", length(genes), " genes in no significant interaction: ",
                                paste(head(result$genes_absent, 20), collapse = ", "), if (length(result$genes_absent) > 20) ", ..." else ""),
                         type = "warning", duration = 8)
      }
      showNotification(
        paste("Chord plot generated with", result$n_interactions, "interactions",
              if (length(genes) > 0) paste0("(", result$n_matching, " of ", result$n_significant, " significant match the gene list)") else "",
              if (result$n_interactions < result$n_matching) paste0("(top ", result$n_interactions, " of ", result$n_matching, ")") else ""),
        type = "message"
      )
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  }))
  
  # Dynamic UI for chord plot with adjustable size
  output$chord_plot_ui <- renderUI({
    if(!is.null(chord_plot())) {
      imageOutput("chord_plot_display", 
                  height = paste0(input$chord_plot_dims, "px"),
                  width = paste0(input$chord_plot_dims, "px"))
    } else {
      div(class = "text-center text-muted", 
          style = "padding: 50px;",
          h4("No plot generated yet"),
          p("Select parameters and click 'Generate Plot'"))
    }
  })
  
  # Display the pre-rendered PNG: circlize must never run inside renderPlot (rule 8)
  output$chord_plot_display <- renderImage({
    req(chord_plot())
    result <- chord_plot()
    dims <- as.integer(input$chord_plot_dims %||% result$png_size)
    list(
      src         = result$output_path,
      contentType = "image/png",
      width       = dims,
      height      = dims,
      alt         = "Chord plot"
    )
  }, deleteFile = FALSE)
  
  # Download handler: redraws on its own device (outside renderPlot) and always closes it
  # Download handler: redraws on its own device (outside renderPlot) and always closes it.
  # Every format uses the physical size of the displayed PNG (png_size px at 150 dpi) and the same margins,
  # so text and sectors keep the same proportions; the export width only sets the pixel count of raster files
  output$download_chord_plot <- downloadHandler(
    filename = function() {
      obj_name <- input$cellchat_obj_chord
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      paste0("ChordPlot_", obj_name, "_", timestamp, ".", input$chord_plot_format)
    },
    content = function(file) {
      req(chord_plot())
      tryCatch({
        size_in <- (chord_plot()$png_size %||% 1200) / 150
        raster_res <- max(72, as.numeric(input$chord_export_dpi %||% 300))
        width_px <- round(size_in * raster_res)
        if (input$chord_plot_format == "pdf") {
          pdf(file, width = size_in, height = size_in)
        } else if (input$chord_plot_format == "svg") {
          svg(file, width = size_in, height = size_in)
        } else {
          # PNG, JPEG, TIFF: same physical size, resolution derived from the requested width
          switch(input$chord_plot_format,
                 "png" = png(file, width = width_px, height = width_px, res = raster_res),
                 "jpeg" = jpeg(file, width = width_px, height = width_px, res = raster_res, quality = 100),
                 "tiff" = tiff(file, width = width_px, height = width_px, res = raster_res, compression = "lzw"))
        }
        export_dev <- dev.cur()
        # Close the export device even if drawing fails, so later plots are not drawn into it
        on.exit(if (export_dev %in% dev.list()) dev.off(export_dev), add = TRUE)
        # Same margins as the displayed PNG
        par(mar = c(0, 0, 0, 0))
        chord_plot()$plot_function()
        
        showNotification("Plot downloaded successfully!", type = "message")
        
      }, error = function(e) {
        showNotification(paste("Download error:", e$message), type = "error")
      })
    }
  )
  # ---- Split chord plots ----
  split_chord_plot   <- reactiveVal(NULL)
  split_chord_lr_data <- reactiveVal(NULL)  
  
  observe({
    req(input$cellchat_obj_split_chord)
    tryCatch({
      obj <- cellchat_objects()[[input$cellchat_obj_split_chord]]
      ct  <- levels(obj@idents)
      updateSelectizeInput(session, "split_chord_senders",   choices = ct)
      updateSelectizeInput(session, "split_chord_receivers", choices = ct)
    }, error = function(e) NULL)
  })
  
  # Dynamic group definition UI for pairs mode — N groups, cell type choices from selected object
  output$split_chord_pairs_ui <- renderUI({
    req(input$split_chord_n_groups)
    n_grp <- max(2L, as.integer(input$split_chord_n_groups))
    ct <- tryCatch({
      obj <- cellchat_objects()[[input$cellchat_obj_split_chord]]
      if (!is.null(obj)) levels(obj@idents) else character(0)
    }, error = function(e) character(0))
    grp_colors <- tryCatch(
      RColorBrewer::brewer.pal(max(3L, n_grp), "Set2")[seq_len(n_grp)],
      error = function(e) scales::hue_pal()(n_grp)
    )
    ui_parts <- lapply(seq_len(n_grp), function(i) {
      tags$div(
        style = "margin-bottom: 4px;",
        tags$p(paste("Group", i),
               style = sprintf("color: %s; font-weight: bold; margin-bottom: 4px;", grp_colors[i])),
        textInput(paste0("split_chord_pair", i, "_label"), "Label:", value = paste("Group", i)),
        selectizeInput(paste0("split_chord_pair", i, "_senders"), "Sources:", choices = ct,
                       multiple = TRUE,
                       options  = list(placeholder = "Select sources", plugins = list("remove_button"))),
        selectizeInput(paste0("split_chord_pair", i, "_receivers"), "Targets:", choices = ct,
                       multiple = TRUE,
                       options  = list(placeholder = "Select targets", plugins = list("remove_button"))),
        if (i < n_grp) tags$hr(style = "border-color: rgba(255,255,255,0.2); margin: 8px 0;") else NULL
      )
    })
    do.call(tagList, ui_parts)
  })
  
  # Dynamic color inputs for split chord plot — initialized from stored colors if available
  output$split_chord_color_inputs <- renderUI({
    tryCatch({
      req(input$split_chord_custom_colors)
      all_cts <- if (isTRUE(input$split_chord_mode == "pairs")) {
        n_grp <- max(2L, as.integer(input$split_chord_n_groups))
        unique(unlist(lapply(seq_len(n_grp), function(i) {
          c(input[[paste0("split_chord_pair", i, "_senders")]],
            input[[paste0("split_chord_pair", i, "_receivers")]])
        })))
      } else {
        unique(c(input$split_chord_senders, input$split_chord_receivers))
      }
      req(length(all_cts) > 0)
      cellchat_obj  <- cellchat_objects()[[input$cellchat_obj_split_chord]]
      stored_colors <- find_stored_colors(cellchat_obj, cluster_names = all_cts)
      n <- length(all_cts)
      fallback_pal <- if (n <= 8) {
        RColorBrewer::brewer.pal(max(3, n), "Set1")[seq_len(n)]
      } else {
        scales::hue_pal()(n)
      }
      names(fallback_pal) <- all_cts
      n_from_obj  <- if (!is.null(stored_colors)) sum(all_cts %in% names(stored_colors)) else 0
      info_banner <- if (n_from_obj == n) {
        tags$p(icon("circle-check"),
               sprintf(" Colors loaded from object (%d/%d)", n_from_obj, n),
               style = "color: #2e7d32; font-size: 12px; margin-bottom: 6px;")
      } else if (n_from_obj > 0) {
        tags$p(icon("triangle-exclamation"),
               sprintf(" Partial match: %d/%d colors from object", n_from_obj, n),
               style = "color: #e65100; font-size: 12px; margin-bottom: 6px;")
      } else {
        tags$p(icon("circle-info"),
               " No stored colors found — using default palette.",
               style = "color: #888; font-size: 12px; margin-bottom: 6px;")
      }
      color_inputs <- lapply(seq_along(all_cts), function(i) {
        ct         <- all_cts[i]
        init_color <- if (!is.null(stored_colors) && ct %in% names(stored_colors)) {
          stored_colors[[ct]]
        } else {
          fallback_pal[[ct]]
        }
        colourpicker::colourInput(
          inputId = paste0("split_chord_color_", i),
          label   = paste("Color for", ct),
          value   = init_color
        )
      })
      div(info_banner, color_inputs)
    }, error = function(e) {
      message("split_chord_color_inputs error: ", conditionMessage(e)[1])
      NULL
    })
  })
  
  
  # Load LR pairs and categorize: common vs unique per group
  observeEvent(input$load_split_chord_lr, {
    req(input$cellchat_obj_split_chord)
    tryCatch({
      obj    <- cellchat_objects()[[input$cellchat_obj_split_chord]]
      thresh <- input$split_chord_threshold
      mode   <- input$split_chord_mode
      
      if (mode == "pairs") {
        n_grp <- max(2L, as.integer(input$split_chord_n_groups))
        groups_def <- lapply(seq_len(n_grp), function(i) {
          lbl <- trimws(input[[paste0("split_chord_pair", i, "_label")]])
          if (nchar(lbl) == 0) lbl <- paste("Group", i)
          list(
            label   = lbl,
            sources = input[[paste0("split_chord_pair", i, "_senders")]],
            targets = input[[paste0("split_chord_pair", i, "_receivers")]]
          )
        })
        missing_grps <- which(sapply(groups_def, function(g) {
          length(g$sources) == 0 || length(g$targets) == 0
        }))
        if (length(missing_grps) > 0)
          stop(sprintf("Groups %s: please select at least one source and one target",
                       paste(missing_grps, collapse = ", ")))
      } else {
        req(input$split_chord_senders, input$split_chord_receivers)
        split_by     <- input$split_chord_by
        sources      <- input$split_chord_senders
        targets      <- input$split_chord_receivers
        split_groups <- if (split_by == "source") sources else targets
        groups_def <- lapply(split_groups, function(grp) {
          list(label   = grp,
               sources = if (split_by == "source") grp else sources,
               targets = if (split_by == "target") grp else targets)
        })
      }
      
      split_groups <- sapply(groups_def, `[[`, "label")
      message("=== load_split_chord_lr ===")
      message("--- mode: ", mode, " | groups: ", paste(split_groups, collapse = ", "))
      
      per_group_lr <- setNames(lapply(groups_def, function(gd) {
        comms <- tryCatch(
          subsetCommunication(obj,
                              sources.use = gd$sources,
                              targets.use = gd$targets,
                              thresh      = thresh),
          error = function(e) NULL
        )
        if (is.null(comms) || nrow(comms) == 0) {
          message(sprintf("  [%s] no interactions found", gd$label))
          return(data.frame(lr_key = character(), prob = numeric()))
        }
        comms$lr_key <- paste(comms$ligand, comms$receptor, sep = " -> ")
        agg <- aggregate(prob ~ lr_key, data = comms, FUN = max)
        agg <- agg[order(agg$prob, decreasing = TRUE), ]
        message(sprintf("  [%s] %d LR pairs found", gd$label, nrow(agg)))
        agg
      }), split_groups)
      
      all_keys      <- unique(unlist(lapply(per_group_lr, function(d) d$lr_key)))
      key_to_groups <- setNames(lapply(all_keys, function(k) {
        sort(names(which(sapply(per_group_lr, function(d) k %in% d$lr_key))))
      }), all_keys)
      
      signatures <- unique(lapply(key_to_groups, identity))
      
      combo_lr <- list()
      for (sig in signatures) {
        label <- if (length(sig) == length(split_groups)) {
          "Common to all"
        } else if (length(sig) == 1) {
          paste0("Unique to ", sig[1])
        } else {
          paste0("Shared: ", paste(sig, collapse = " + "))
        }
        keys_in_sig <- names(which(sapply(key_to_groups, function(g) identical(g, sig))))
        probs <- sapply(keys_in_sig, function(k) {
          mean(sapply(sig, function(grp) {
            d <- per_group_lr[[grp]]
            v <- d$prob[d$lr_key == k]
            if (length(v) == 0) 0 else v[1]
          }))
        })
        combo_lr[[label]] <- keys_in_sig[order(probs, decreasing = TRUE)]
        message(sprintf("  [%s] %d LR pairs", label, length(keys_in_sig)))
      }
      
      sort_order <- sapply(names(combo_lr), function(label) {
        if (label == "Common to all")       return(0)
        if (startsWith(label, "Shared: "))  return(1)
        return(2)
      })
      combo_lr <- combo_lr[order(sort_order)]
      
      split_chord_lr_data(list(
        combo_lr     = combo_lr,
        per_group    = per_group_lr,
        split_groups = split_groups,
        groups_def   = groups_def
      ))
      
      total_lr <- length(all_keys)
      n_common <- length(combo_lr[["Common to all"]])
      showNotification(
        sprintf("Found %d LR pairs total (%d common to all, %d combinations)",
                total_lr, n_common, length(combo_lr)),
        type = "message", duration = 4
      )
      
    }, error = function(e) {
      showNotification(paste("Error:", conditionMessage(e)[1]), type = "error", duration = 5)
      message("load_split_chord_lr error: ", conditionMessage(e)[1])
    })
  })
  
  # Build dynamic LR selector UI: one selectize for common + one per group
  output$split_chord_lr_selector <- renderUI({
    lr_data <- split_chord_lr_data()
    if (is.null(lr_data)) {
      return(tags$p("Click 'Load LR pairs' first",
                    style = "color: #ccc; font-style: italic;"))
    }
    combo_lr  <- lr_data$combo_lr
    per_group <- lr_data$per_group
    
    ui_parts <- lapply(seq_along(combo_lr), function(i) {
      label    <- names(combo_lr)[i]
      keys     <- combo_lr[[i]]
      input_id <- paste0("split_chord_lr_combo_", i)
      if (length(keys) == 0) return(NULL)
      
      probs <- sapply(keys, function(k) {
        vals <- sapply(per_group, function(d) {
          v <- d$prob[d$lr_key == k]
          if (length(v) == 0) NA_real_ else v[1]
        })
        round(mean(vals, na.rm = TRUE), 4)
      })
      choices <- setNames(keys, sprintf("%s  [%.4f]", keys, probs))
      
      lbl_icon <- if (label == "Common to all")     icon("link") else
        if (startsWith(label, "Shared: ")) icon("code-branch") else
          icon("circle-dot")
      
      btn_id    <- paste0("split_chord_select_all_", i)
      btn_label <- if (label == "Common to all") {
        tags$span(icon("check-double"), " Include all common",
                  style = "font-size: 12px; color: #ccc;")
      } else if (startsWith(label, "Unique to ")) {
        tags$span(icon("check-double"), " Include all unique",
                  style = "font-size: 12px; color: #ccc;")
      } else {
        tags$span(icon("check-double"), " Include shared",
                  style = "font-size: 12px; color: #ccc;")
      }
      
      tags$div(
        style = "margin-bottom: 8px;",
        tags$div(style = "margin-bottom: 4px;",
                 checkboxInput(btn_id, label = btn_label, value = TRUE)),
        selectizeInput(
          input_id,
          label    = tags$span(lbl_icon, paste0(" ", label, ":")),
          choices  = choices,
          selected = keys,
          multiple = TRUE,
          options  = list(plugins     = list("remove_button"),
                          placeholder = "None selected",
                          maxOptions  = 500)
        )
      )
    })
    
    do.call(tagList, Filter(Negate(is.null), ui_parts))
  })
    
  # Last known state of each "Include all" checkbox, so only the box that changed updates its own selector
  split_chord_checkbox_state <- new.env(parent = emptyenv())
  observe(run_safely("Interface update", persistent = FALSE, {
    lr_data <- split_chord_lr_data()
    req(lr_data)
    combo_lr <- lr_data$combo_lr
    # New LR pairs: selectors were just rendered with every pair selected and every checkbox ticked
    if (!identical(split_chord_checkbox_state$combo_lr, combo_lr)) {
      split_chord_checkbox_state$combo_lr <- combo_lr
      split_chord_checkbox_state$checked <- rep(TRUE, length(combo_lr))
    }
    for (i in seq_along(combo_lr)) {
      checked <- input[[paste0("split_chord_select_all_", i)]]
      if (is.null(checked) || identical(isTRUE(checked), split_chord_checkbox_state$checked[i])) next
      split_chord_checkbox_state$checked[i] <- isTRUE(checked)
      updateSelectizeInput(session, paste0("split_chord_lr_combo_", i), selected = if (isTRUE(checked)) combo_lr[[i]] else character(0))
    }
  }))
  
  # Generate split chord plots
  observeEvent(input$generate_split_chord, run_safely("generate split chord", {
    mode <- input$split_chord_mode
    req(input$cellchat_obj_split_chord)
    if (mode == "pairs") {
      n_grp <- max(2L, as.integer(input$split_chord_n_groups))
      for (i in seq_len(n_grp)) {
        req(input[[paste0("split_chord_pair", i, "_senders")]],
            input[[paste0("split_chord_pair", i, "_receivers")]])
      }
    } else {
      req(input$split_chord_senders, input$split_chord_receivers)
    }
    
    tryCatch({
      showNotification("Generating split chord plots...", type = "message", duration = 2)
      obj     <- cellchat_objects()[[input$cellchat_obj_split_chord]]
      lr_data <- split_chord_lr_data()
      
      if (is.null(lr_data)) {
        showNotification("Please click 'Load LR pairs' first", type = "warning")
        return()
      }
      
      n_combos    <- length(lr_data$combo_lr)
      selected_lr <- unique(unlist(lapply(seq_len(n_combos), function(i) {
        input[[paste0("split_chord_lr_combo_", i)]]
      })))
      selected_lr <- selected_lr[!is.null(selected_lr) & nchar(selected_lr) > 0]
      
      if (length(selected_lr) == 0) {
        showNotification("Please select at least one LR pair", type = "warning")
        return()
      }
      
      message(sprintf("generate_split_chord: %d LR pairs selected", length(selected_lr)))
      
      all_cts <- if (mode == "pairs") {
        n_grp <- max(2L, as.integer(input$split_chord_n_groups))
        unique(unlist(lapply(seq_len(n_grp), function(i) {
          c(input[[paste0("split_chord_pair", i, "_senders")]],
            input[[paste0("split_chord_pair", i, "_receivers")]])
        })))
      } else {
        unique(c(input$split_chord_senders, input$split_chord_receivers))
      }
      
      custom_colors <- NULL
      if (isTRUE(input$split_chord_custom_colors)) {
        n <- length(all_cts)
        fallback_pal <- if (n <= 8) {
          RColorBrewer::brewer.pal(max(3, n), "Set1")[seq_len(n)]
        } else {
          scales::hue_pal()(n)
        }
        custom_colors <- sapply(seq_along(all_cts), function(i) {
          col <- input[[paste0("split_chord_color_", i)]]
          if (is.null(col)) fallback_pal[i] else col
        })
        names(custom_colors) <- all_cts
      }
      
      pairs_arg <- if (mode == "pairs") lr_data$groups_def else NULL
      
      result <- generate_split_chord_plots(
        cellchat_obj     = obj,
        sources          = if (mode == "split") input$split_chord_senders   else NULL,
        targets          = if (mode == "split") input$split_chord_receivers else NULL,
        split_by         = if (mode == "split") input$split_chord_by        else "source",
        pairs            = pairs_arg,
        selected_lr_keys = selected_lr,
        threshold        = input$split_chord_threshold,
        label_distance   = input$split_chord_label_distance,
        custom_colors    = custom_colors,
        show_title       = isTRUE(input$split_chord_show_title),
        gap_inner        = input$split_chord_gap_inner,
        gap_outer        = input$split_chord_gap_outer,
        plot_size        = as.integer(input$split_chord_dims)
      )
      
      split_chord_plot(result)
      
    }, error = function(e) {
      showNotification(paste("Error:", conditionMessage(e)[1]), type = "error", duration = 6)
      message("generate_split_chord error: ", conditionMessage(e)[1])
    })
  }))
  
  
  split_chord_h <- reactiveVal(600)
  split_chord_w <- reactiveVal(600)
  
  observeEvent(split_chord_plot(), run_safely("Interface update", persistent = FALSE, {
    req(split_chord_plot())
    split_chord_h(split_chord_plot()$png_size)
    split_chord_w(split_chord_plot()$png_size * split_chord_plot()$n_plots)
  }))
  
  output$split_chord_plot_display <- renderImage({
    req(split_chord_plot())
    result <- split_chord_plot()
    list(
      src     = result$output_path,
      contentType = "image/png",
      width   = result$png_size * result$n_plots,
      height  = result$png_size,
      alt     = "Split chord plot"
    )
  }, deleteFile = FALSE)
  
  output$download_split_chord <- downloadHandler(
    filename = function() generateFileName(input$cellchat_obj_split_chord, "SplitChord", input$split_chord_format),
    content = function(file) {
      req(split_chord_plot())
      result    <- split_chord_plot()
      n_plots   <- result$n_plots
      width_px  <- as.integer(input$split_chord_export_width %||% 1200)
      if (is.na(width_px) || width_px < 100) width_px <- 1200L
      dpi       <- as.integer(input$split_chord_export_dpi %||% 300)
      width_in  <- 8 * n_plots
      height_in <- 8
      switch(input$split_chord_format,
             "svg"  = svg(file, width = width_in, height = height_in),
             "pdf"  = pdf(file, width = width_in, height = height_in, pointsize = 12),
             "tiff" = tiff(file, width = width_px * n_plots, height = width_px, res = dpi, units = "px", compression = "lzw"),
             "jpeg" = jpeg(file, width = width_px * n_plots, height = width_px, res = dpi, units = "px", quality = 100),
             png(file, width = width_px * n_plots, height = width_px, res = dpi, units = "px"))
      export_dev <- dev.cur()
      # The export device is closed even when drawing fails
      on.exit(if (export_dev %in% dev.list()) dev.off(export_dev), add = TRUE)
      result$export_function(n_plots)
    }
  )
  
  # Summary table: all LR pairs found, prob per group, category
  output$split_chord_lr_table_ui <- renderUI({
    lr_data <- split_chord_lr_data()
    if (is.null(lr_data)) {
      return(tags$p("Click 'Load LR pairs' first",
                    style = "color: #ccc; font-style: italic;"))
    }
    
    split_groups <- lr_data$split_groups
    per_group    <- lr_data$per_group
    combo_lr     <- lr_data$combo_lr
    groups_def   <- lr_data$groups_def  # NULL in split mode
    
    # Build human-readable column names: "Label (Src → Tgt)"
    group_col_names <- sapply(seq_along(split_groups), function(i) {
      grp <- split_groups[i]
      gd  <- if (!is.null(groups_def)) groups_def[[i]] else NULL
      if (!is.null(gd)) {
        src <- paste(gd$sources, collapse = "+")
        tgt <- paste(gd$targets, collapse = "+")
        paste0(grp, " (", src, " \u2192 ", tgt, ")")
      } else {
        grp
      }
    })
    
    # Build one row per LR key — split ligand/receptor
    all_keys <- unique(unlist(combo_lr))
    rows <- lapply(all_keys, function(k) {
      parts    <- strsplit(k, " -> ")[[1]]
      ligand   <- if (length(parts) >= 1) parts[1] else k
      receptor <- if (length(parts) >= 2) parts[2] else ""
      probs <- sapply(split_groups, function(grp) {
        d <- per_group[[grp]]
        v <- d$prob[d$lr_key == k]
        if (length(v) == 0) NA_real_ else round(v[1], 5)
      })
      category <- names(which(sapply(combo_lr, function(keys) k %in% keys)))[1]
      as.data.frame(
        t(c("Ligand"   = ligand,
            "Receptor" = receptor,
            setNames(as.character(probs), group_col_names),
            "Category" = category)),
        stringsAsFactors = FALSE
      )
    })
    
    df <- do.call(rbind, rows)
    
    priority <- ifelse(df$Category == "Common to all", 0,
                       ifelse(startsWith(df$Category, "Shared: "), 1, 2))
    max_prob <- apply(df[, group_col_names, drop = FALSE], 1, function(r) {
      vals <- suppressWarnings(as.numeric(r))
      if (all(is.na(vals))) -Inf else max(vals, na.rm = TRUE)
    })
    df <- df[order(priority, -max_prob), ]
    
    all_labels <- unique(df$Category)
    cat_colors <- setNames(scales::hue_pal()(length(all_labels)), all_labels)
    if ("Common to all" %in% all_labels) cat_colors["Common to all"] <- "#2e7d32"
    unique_labels <- all_labels[startsWith(all_labels, "Unique to ")]
    if (length(unique_labels) > 0)
      cat_colors[unique_labels] <- scales::brewer_pal(palette = "Oranges")(
        max(3, length(unique_labels)))[seq_along(unique_labels) + 1]
    
    prob_col_indices <- seq(3, 2 + length(group_col_names))
    
    DT::renderDataTable(
      DT::datatable(
        df,
        rownames   = FALSE,
        selection  = "none",
        extensions = "Buttons",
        options    = list(
          pageLength = nrow(df),
          dom        = "Bft",
          scrollX    = TRUE,
          scrollY    = "350px",
          buttons    = list(
            list(extend   = "csv",
                 filename = paste0("SplitChord_LR_", format(Sys.time(), "%Y%m%d")),
                 text     = "CSV"),
            list(extend   = "excel",
                 filename = paste0("SplitChord_LR_", format(Sys.time(), "%Y%m%d")),
                 text     = "Excel",
                 title    = paste("Split Chord — LR pairs |",
                                  paste(sapply(seq_along(split_groups), function(i) {
                                    gd <- if (!is.null(groups_def)) groups_def[[i]] else NULL
                                    if (!is.null(gd))
                                      paste0(split_groups[i], ": ",
                                             paste(gd$sources, collapse = "+"),
                                             " \u2192 ",
                                             paste(gd$targets, collapse = "+"))
                                    else split_groups[i]
                                  }), collapse = "  |  "))
            )
          ),
          columnDefs = list(
            list(className = "dt-center", targets = prob_col_indices - 1)
          )
        )
      ) |>
        DT::formatStyle(
          "Category",
          backgroundColor = DT::styleEqual(all_labels, unname(cat_colors[all_labels])),
          color      = "white",
          fontWeight = "bold"
        ) |>
        DT::formatStyle(
          group_col_names,
          backgroundColor = DT::styleInterval(0, c("#555555", "#1a5276")),
          color = "white"
        ) |>
        DT::formatStyle(
          group_col_names,
          background = DT::styleColorBar(c(0, 0.3), "#5dade2"),
          backgroundSize   = "98% 60%",
          backgroundRepeat = "no-repeat",
          backgroundPosition = "center"
        ),
      server = FALSE
    )
  })
}
