# ==============================================================================
# enrichment_functions.R
# Pure helper functions for pathway enrichment (no Shiny reactive code).
# Gene Ontology for mouse, human and rat through their OrgDb packages (offline, gene symbols), or user gene sets
# from a GMT file for any species. Over-representation analysis (ORA, one-sided hypergeometric test) of a gene list against
# a universe, and gene set enrichment analysis (GSEA, fgsea multilevel) of a ranked gene list, seed 42.
# clusterProfiler is not used: its Bioconductor 3.20 dependency chain (enrichplot -> ggtree) does not build with ggplot2 4.x.
# The ORA statistics follow clusterProfiler::enricher (same universe, set-size filter, GeneRatio / BgRatio, BH correction).
# Reactive logic lives in single_dataset_server.R (IDs *_single) and multiple_datasets_server.R (IDs *_merge).
# ==============================================================================

############################## Gene sets ##############################

# OrgDb object of a species with built-in Gene Ontology annotation (mouse, human, rat)
enrichment_orgdb <- function(species) {
  pkg <- switch(species, mouse = "org.Mm.eg.db", human = "org.Hs.eg.db", rat = "org.Rn.eg.db",
                stop("No built-in Gene Ontology annotation for '", species, "': load a GMT file"))
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package ", pkg, " is not installed in this image")
  return(getExportedValue(pkg, pkg))
}

# Text describing the built-in annotation of a species for the parameters table: OrgDb package version and GO release date.
# The GO date is the first non-empty GO*SOURCEDATE field (GOSOURCEDATE preferred) of the OrgDb metadata, then of GO.db
# (in Bioconductor 3.20 GOSOURCEDATE is empty and the date is in GOEGSOURCEDATE).
orgdb_description <- function(species) {
  orgdb <- enrichment_orgdb(species)
  pkg <- switch(species, mouse = "org.Mm.eg.db", human = "org.Hs.eg.db", rat = "org.Rn.eg.db")
  go_date_from <- function(db) {
    md <- tryCatch(AnnotationDbi::metadata(db), error = function(e) NULL)
    if (is.null(md)) return(NA_character_)
    md <- md[grepl("^GO.*SOURCEDATE$", md$name) & !is.na(md$value) & nzchar(trimws(md$value)), , drop = FALSE]
    if (nrow(md) == 0) return(NA_character_)
    md <- md[order(md$name != "GOSOURCEDATE"), , drop = FALSE]
    return(md$value[1])
  }
  go_date <- go_date_from(orgdb)
  if (is.na(go_date) && requireNamespace("GO.db", quietly = TRUE)) go_date <- go_date_from(GO.db::GO.db)
  return(sprintf("Gene Ontology, %s (%s %s, GO source date %s)", species, pkg, as.character(utils::packageVersion(pkg)), if (is.na(go_date)) "unknown" else go_date))
}
# Gene sets of a GMT file (one set per line: name <tab> description <tab> gene <tab> gene ...) as a two-column
# data.frame (term, gene). Duplicated pairs and empty genes are removed.
read_gmt_gene_sets <- function(path) {
  lines <- readLines(path, warn = FALSE)
  parts <- strsplit(lines[nzchar(trimws(lines))], "\t", fixed = TRUE)
  parts <- parts[lengths(parts) >= 3]
  if (length(parts) == 0) stop("No gene set found: a GMT file has one set per line, name <tab> description <tab> genes separated by tabs")
  term2gene <- data.frame(term = rep(trimws(vapply(parts, `[`, character(1), 1)), lengths(parts) - 2),
                          gene = trimws(unlist(lapply(parts, function(p) p[-(1:2)]), use.names = FALSE)), stringsAsFactors = FALSE)
  term2gene <- unique(term2gene[nzchar(term2gene$gene), , drop = FALSE])
  return(term2gene)
}

# Gene sets restricted to the dataset genes. Genes are matched exactly, then case-insensitively (e.g. "MYH7" -> "Myh7").
# Case-insensitive matching is not an orthology mapping: gene sets of another species should be converted beforehand.
# Returns the restricted sets and the fraction of set genes found in the dataset.
match_gene_sets_to_dataset <- function(term2gene, dataset_genes) {
  exact <- term2gene$gene %in% dataset_genes
  lookup <- stats::setNames(dataset_genes, toupper(dataset_genes))
  lookup <- lookup[!duplicated(names(lookup))]
  mapped <- ifelse(exact, term2gene$gene, unname(lookup[toupper(term2gene$gene)]))
  out <- unique(data.frame(term = term2gene$term, gene = mapped, stringsAsFactors = FALSE)[!is.na(mapped), , drop = FALSE])
  return(list(term2gene = out, fraction_matched = mean(!is.na(mapped)), n_sets = length(unique(out$term))))
}

# Gene Ontology gene sets of the given genes (gene symbols) from the species OrgDb, as a data.frame (term, gene).
# GOALL includes the annotations inherited from child terms (true-path rule), as clusterProfiler::enrichGO does.
go_gene_sets <- function(genes, species, ontology = "BP") {
  ann <- suppressMessages(AnnotationDbi::select(enrichment_orgdb(species), keys = unique(genes), keytype = "SYMBOL", columns = c("GOALL", "ONTOLOGYALL")))
  ann <- ann[!is.na(ann$GOALL) & ann$ONTOLOGYALL == ontology, , drop = FALSE]
  return(unique(data.frame(term = ann$GOALL, gene = ann$SYMBOL, stringsAsFactors = FALSE)))
}

# Readable names of gene set IDs: GO term names from GO.db for GO IDs, the ID itself otherwise (GMT set names)
gene_set_names <- function(ids) {
  out <- stats::setNames(ids, ids)
  go <- ids[grepl("^GO:[0-9]{7}$", ids)]
  if (length(go) > 0 && requireNamespace("GO.db", quietly = TRUE)) {
    # A failed GO.db lookup leaves the GO IDs as names instead of stopping the analysis
    terms <- tryCatch(suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = go, keytype = "GOID", columns = "TERM")), error = function(e) NULL)
    if (!is.null(terms)) {
      terms <- terms[!is.na(terms$TERM) & !duplicated(terms$GOID), , drop = FALSE]
      out[terms$GOID] <- terms$TERM
    }
  }
  return(unname(out))
}

# Gene sets for a test: GO of the species when term2gene is NULL, otherwise the user sets; restricted to the universe
# and to sets of min_size-max_size universe genes. Returns a named list of gene vectors (fgsea pathway format).
prepare_gene_sets <- function(universe, species, ontology, term2gene, min_size, max_size) {
  t2g <- if (is.null(term2gene)) go_gene_sets(universe, species, ontology) else term2gene[term2gene$gene %in% universe, c("term", "gene"), drop = FALSE]
  if (nrow(t2g) == 0) stop("No gene of the universe is annotated in the selected gene sets")
  sets <- split(t2g$gene, t2g$term)
  sets <- lapply(sets, unique)
  sizes <- lengths(sets)
  return(list(sets = sets[sizes >= min_size & sizes <= max_size], annotated = unique(t2g$gene)))
}

############################## Inputs ##############################

# Genes detected (count > 0) in at least min_pct of the cells of an assay: ORA universe when the tested genes are unknown
detected_genes <- function(obj, assay = NULL, min_pct = 0.01) {
  assay <- assay %||% SeuratObject::DefaultAssay(obj)
  counts <- tryCatch(SeuratObject::LayerData(obj, assay = assay, layer = "counts"), error = function(e) NULL)
  if (is.null(counts)) counts <- SeuratObject::LayerData(obj, assay = assay, layer = "data")
  return(rownames(counts)[Matrix::rowMeans(counts > 0) >= min_pct])
}

# Differential expression table (FindMarkers or pseudobulk columns, genes in a 'gene' column or as row names) to a clean data.frame.
# rank_stat (GSEA ranking) is kept when present (DESeq2 Wald statistic), otherwise it is the log2 fold change.
standardize_de_table <- function(de) {
  de <- as.data.frame(de)
  if (!"gene" %in% colnames(de)) de$gene <- rownames(de)
  if (!all(c("avg_log2FC", "p_val_adj") %in% colnames(de))) stop("The table needs avg_log2FC and p_val_adj columns")
  de$avg_log2FC <- as.numeric(de$avg_log2FC); de$p_val_adj <- as.numeric(de$p_val_adj)
  de$rank_stat <- if ("rank_stat" %in% colnames(de)) as.numeric(de$rank_stat) else de$avg_log2FC
  de <- de[!is.na(de$gene) & !duplicated(de$gene) & is.finite(de$avg_log2FC), c("gene", "avg_log2FC", "p_val_adj", "rank_stat"), drop = FALSE]
  rownames(de) <- NULL
  return(de)
}

# Significant genes of a DE table for ORA: adjusted p-value below the cutoff and fold change beyond the threshold
# in the chosen direction ("up" or "down")
select_ora_genes <- function(de, direction = "up", max_p_adj = 0.05, min_abs_log2fc = 0.25) {
  keep <- !is.na(de$p_val_adj) & de$p_val_adj <= max_p_adj
  keep <- keep & if (identical(direction, "down")) de$avg_log2FC <= -min_abs_log2fc else de$avg_log2FC >= min_abs_log2fc
  return(de$gene[keep])
}

# Genes ranked for GSEA by rank_stat (decreasing). Ties are broken at random with a fixed seed rather than by gene name,
# which would place gene families with similar names (Rpl, Rps, Nduf...) next to each other within tied values.
# The fraction of genes sharing their value with another gene is returned as an attribute for the parameters table.
rank_genes_for_gsea <- function(de, seed = 42) {
  de <- de[is.finite(de$rank_stat), , drop = FALSE]
  if (nrow(de) == 0) stop("No gene with a finite ranking statistic")
  set.seed(seed)
  tiebreak <- sample.int(nrow(de))
  de <- de[order(-de$rank_stat, tiebreak), , drop = FALSE]
  ranked <- stats::setNames(de$rank_stat, de$gene)
  attr(ranked, "fraction_tied") <- mean(duplicated(de$rank_stat) | duplicated(de$rank_stat, fromLast = TRUE))
  return(ranked)
}

# Metadata columns usable to group cells (factor or character, 2 to 200 values)
enrichment_grouping_columns <- function(md) {
  return(colnames(md)[vapply(md, function(v) (is.factor(v) || is.character(v)) && length(unique(stats::na.omit(v))) %in% 2:200, logical(1))])
}

# Values of a grouping column (factor levels in their order, otherwise sorted)
enrichment_group_values <- function(v) {
  return(if (is.factor(v)) levels(droplevels(v)) else sort(unique(as.character(stats::na.omit(v)))))
}

# Choices of the stored-table selector: differential expression tables labelled by their description
# (stored enrichment results are not gene tables and are left out)
stored_table_choices <- function(tables) {
  tables <- tables[vapply(tables, function(e) !identical(e$type, "enrichment"), logical(1))]
  if (length(tables) == 0) return(character(0))
  return(stats::setNames(names(tables), vapply(names(tables), function(n) tables[[n]]$description %||% n, character(1))))
}

# Cells of the two groups of a comparison. group2 empty = all other cells with a value; restrict_values limits both groups
# to cells of the active identities (e.g. one cluster).
select_comparison_cells <- function(obj, group_col, group1, group2 = character(0), restrict_values = character(0)) {
  if (is.null(group_col) || !group_col %in% colnames(obj@meta.data)) stop("Select a grouping column")
  if (length(group1) == 0) stop("Select at least one value for group 1")
  if (length(intersect(group1, group2)) > 0) stop("A value cannot be in both groups")
  in_scope <- if (length(restrict_values) > 0) as.character(Seurat::Idents(obj)) %in% restrict_values else rep(TRUE, ncol(obj))
  v <- as.character(obj@meta.data[[group_col]])
  cells1 <- colnames(obj)[in_scope & !is.na(v) & v %in% group1]
  cells2 <- colnames(obj)[in_scope & !is.na(v) & if (length(group2) > 0) v %in% group2 else !v %in% group1]
  if (length(cells1) < 3 || length(cells2) < 3) stop(sprintf("Too few cells: %d in group 1, %d in group 2", length(cells1), length(cells2)))
  return(list(cells1 = cells1, cells2 = cells2))
}

############################## Tests ##############################

# Over-representation analysis: one-sided hypergeometric test per gene set, Benjamini-Hochberg over the tested sets.
# Same definitions as clusterProfiler::enricher / enrichGO: N = universe genes annotated in the gene sets, n = query genes
# among them, M = set size within the universe (sets of min_size-max_size genes), k = query genes in the set; only sets with
# k >= 1 are tested. Gene Ontology of the species OrgDb when term2gene is NULL, otherwise the user gene sets.
# Every tested set is returned; filtering is done at display time.
run_ora <- function(genes, universe, species = "mouse", ontology = "BP", term2gene = NULL, min_size = 10, max_size = 500) {
  universe <- unique(universe)
  gs <- prepare_gene_sets(universe, species, ontology, term2gene, min_size, max_size)
  genes <- intersect(unique(genes), gs$annotated)
  if (length(genes) < 3) stop(sprintf("Only %d significant gene(s) annotated in the universe: relax the thresholds", length(genes)))
  N <- length(gs$annotated); n <- length(genes)
  overlap <- lapply(gs$sets, intersect, genes)
  k <- lengths(overlap); keep <- k > 0
  if (!any(keep)) return(list(table = data.frame(), n_genes = n, n_universe = N))
  overlap <- overlap[keep]; k <- k[keep]; M <- lengths(gs$sets[keep])
  pvalue <- stats::phyper(k - 1, M, N - M, n, lower.tail = FALSE)
  table <- data.frame(ID = names(overlap), Description = gene_set_names(names(overlap)), GeneRatio = paste0(k, "/", n), BgRatio = paste0(M, "/", N),
                      pvalue = pvalue, p.adjust = stats::p.adjust(pvalue, method = "BH"),
                      geneID = vapply(overlap, function(g) paste(sort(g), collapse = "/"), character(1)), Count = unname(k),
                      gene_ratio = unname(k / n), fold_enrichment = unname((k / n) / (M / N)), stringsAsFactors = FALSE)
  table <- table[order(table$p.adjust, -table$gene_ratio, table$ID), , drop = FALSE]
  rownames(table) <- NULL
  return(list(table = table, n_genes = n, n_universe = N))
}

# ORA run separately on the up- and down-regulated genes ("both"), or on one direction. Each direction is corrected on its own
# (Benjamini-Hochberg over the sets tested for that list). A direction with too few significant genes is skipped with a note.
# Returns list(table (with a 'direction' column), n_genes (named by direction), n_universe, notes).
run_ora_directions <- function(de, direction, universe, max_p_adj, min_log2fc, species, ontology, term2gene, min_size, max_size) {
  directions <- if (identical(direction, "both")) c("up", "down") else direction
  labels <- c(up = "Up in group 1", down = "Down in group 1")
  tables <- list(); n_genes <- integer(0); notes <- character(0)
  for (d in directions) {
    res <- tryCatch(run_ora(select_ora_genes(de, d, max_p_adj, min_log2fc), universe, species, ontology, term2gene, min_size, max_size),
                    error = function(e) { notes <<- c(notes, sprintf("%s: %s", labels[[d]], conditionMessage(e))); NULL })
    if (is.null(res)) next
    n_genes[[d]] <- res$n_genes
    if (nrow(res$table) > 0) tables[[d]] <- cbind(direction = labels[[d]], res$table, stringsAsFactors = FALSE)
  }
  if (length(n_genes) == 0) stop(paste(notes, collapse = "; "))
  table <- if (length(tables) == 0) data.frame() else do.call(rbind, unname(tables))
  rownames(table) <- NULL
  return(list(table = table, n_genes = n_genes, n_universe = length(universe), notes = notes))
}

# Gene set enrichment analysis with fgsea multilevel (eps = 0 for exact small p-values, 10,000 simple permutations,
# single process so the seed fully determines the result), Benjamini-Hochberg over the tested sets.
# Gene sets restricted to the ranked genes and to min_size-max_size of them. Every tested set is returned.
run_gsea <- function(ranked, species = "mouse", ontology = "BP", term2gene = NULL, min_size = 10, max_size = 500, seed = 42) {
  ranked <- ranked[is.finite(ranked) & !is.na(names(ranked)) & !duplicated(names(ranked))]
  if (length(ranked) < 50) stop(sprintf("Only %d ranked genes: GSEA needs the full list of tested genes", length(ranked)))
  ranked <- sort(ranked, decreasing = TRUE)
  gs <- prepare_gene_sets(names(ranked), species, ontology, term2gene, min_size, max_size)
  if (length(gs$sets) == 0) stop(sprintf("No gene set of %d-%d ranked genes", min_size, max_size))
  set.seed(seed)
  # Ties in log2 fold changes are expected (many genes near 0) and only affect the order within tied genes
  res <- withCallingHandlers(fgsea::fgsea(pathways = gs$sets, stats = ranked, minSize = min_size, maxSize = max_size, eps = 0, nPermSimple = 10000, nproc = 1),
                             warning = function(w) if (grepl("ties in the preranked stats", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning"))
  res <- as.data.frame(res[!is.na(res$pval), ])
  if (nrow(res) == 0) return(list(table = data.frame(), n_ranked = length(ranked)))
  table <- data.frame(ID = res$pathway, Description = gene_set_names(res$pathway), setSize = res$size, enrichmentScore = res$ES, NES = res$NES,
                      pvalue = res$pval, p.adjust = stats::p.adjust(res$pval, method = "BH"),
                      core_enrichment = vapply(res$leadingEdge, paste, character(1), collapse = "/"), stringsAsFactors = FALSE)
  table <- table[order(table$p.adjust, -abs(table$NES), table$ID), , drop = FALSE]
  rownames(table) <- NULL
  return(list(table = table, n_ranked = length(ranked)))
}

# Full enrichment run from the tab settings (no Shiny code: inputs are passed as a list).
# settings: source ("compare"/"stored"), group_col, group1, group2, restrict_values, min_pct, stored_table, gene_sets
# ("mouse"/"human"/"rat"/"gmt"), ontology, min_size, max_size, methods, direction ("up"/"down"/"both"), ora_max_p_adj, ora_min_log2fc,
# de_method ("wilcox", default, or "pseudobulk"), sample_col and pb_min_cells (pseudobulk only).
# Differential expression always uses the log-normalized RNA assay (prepare_rna_assay_for_de, cells_genes_expressions.R).
# Wilcoxon: cell-level FindMarkers, ORA genes by Bonferroni-adjusted p, GSEA ranked by log2 fold change (descriptive for
# condition comparisons). Pseudobulk: DESeq2 on raw counts summed per sample (run_pseudobulk_de, cells_genes_expressions.R),
# ORA genes by Benjamini-Hochberg adjusted p, GSEA ranked by the Wald statistic.
# stored_tables: the DE table storage (list of entries with $data and $description). gmt: list(term2gene, file, fraction_matched, md5) or NULL.
# Returns list(ora, gsea, label, params, gsea_skipped, gene_set_label); gsea_skipped is TRUE when GSEA was asked on a stored table.
run_enrichment_analysis <- function(obj, settings, stored_tables = NULL, gmt = NULL) {
  methods <- settings$methods
  if (length(methods) == 0) stop("Select ORA, GSEA or both")
  use_gmt <- identical(settings$gene_sets, "gmt")
  if (use_gmt && is.null(gmt)) stop("Load a GMT file first")
  term2gene <- if (use_gmt) gmt$term2gene else NULL
  min_pct <- settings$min_pct %||% 0.05; gsea_skipped <- FALSE
  de_method <- if (identical(settings$de_method, "pseudobulk")) "pseudobulk" else "wilcox"
  prep <- prepare_rna_assay_for_de(obj, "RNA")
  obj <- prep$obj; assay <- prep$assay
  assay_desc <- paste0(assay, " assay, log-normalized data", if (prep$normalized_here) " (LogNormalize applied for this analysis: the RNA assay had no data layer)" else "")
  ranking_desc <- "log2 fold change"
  if (identical(settings$source, "stored")) {
    entry <- stored_tables[[settings$stored_table %||% ""]]
    if (is.null(entry)) stop("Select a stored table")
    de <- standardize_de_table(entry$data)
    universe <- union(detected_genes(obj, assay = assay, min_pct = 0.01), de$gene)
    label <- entry$description %||% settings$stored_table
    source_desc <- "Stored differential expression table; universe = genes detected in >= 1% of all cells plus the genes of the table"
    # A stored table holds only the genes that passed its thresholds: no full ranking for GSEA
    if ("GSEA" %in% methods) { gsea_skipped <- TRUE; methods <- setdiff(methods, "GSEA") }
    if (length(methods) == 0) stop("GSEA needs all tested genes ranked: it runs only on a new comparison")
  } else {
    cells <- select_comparison_cells(obj, settings$group_col, settings$group1, settings$group2, settings$restrict_values)
    label <- sprintf("%s: %s vs %s%s", settings$group_col, paste(settings$group1, collapse = "+"),
                     if (length(settings$group2) > 0) paste(settings$group2, collapse = "+") else "all other cells",
                     if (length(settings$restrict_values) > 0) paste0(" (clusters ", paste(settings$restrict_values, collapse = ", "), ")") else "")
    if (identical(de_method, "pseudobulk")) {
      sample_col <- settings$sample_col
      if (is.null(sample_col) || !sample_col %in% colnames(obj@meta.data)) stop("Select the sample column (one value = one animal / biological replicate)")
      condition <- rep(NA_character_, ncol(obj))
      condition[colnames(obj) %in% cells$cells1] <- "Group1"
      condition[colnames(obj) %in% cells$cells2] <- "Group2"
      pb <- run_pseudobulk_de(obj, condition, sample_col, assay = assay, min_cells = settings$pb_min_cells %||% 10, min_samples = 2)
      design <- attr(pb, "design")
      n_samples <- c(sum(design$used & design$condition == "Group1"), sum(design$used & design$condition == "Group2"))
      # Wald statistic of DESeq2 (log2 fold change / standard error) ranks the genes for GSEA
      pb$rank_stat <- pb$avg_log2FC / pb$lfcSE
      de <- standardize_de_table(pb)
      ranking_desc <- "DESeq2 Wald statistic"
      source_desc <- sprintf("New comparison, pseudobulk DESeq2 on raw RNA counts summed per sample (column %s), %d vs %d samples (%d vs %d cells, samples with >= %d cells); ORA genes by Benjamini-Hochberg adjusted p",
                             sample_col, n_samples[1], n_samples[2], length(cells$cells1), length(cells$cells2), settings$pb_min_cells %||% 10)
    } else {
      de <- standardize_de_table(Seurat::FindMarkers(obj, ident.1 = cells$cells1, ident.2 = cells$cells2, assay = assay, logfc.threshold = 0, min.pct = min_pct, verbose = FALSE))
      source_desc <- sprintf("New comparison, Wilcoxon (FindMarkers, logfc.threshold = 0, min.pct = %g), %d vs %d cells; ORA genes by Bonferroni adjusted p (cell-level test: descriptive for condition comparisons)",
                             min_pct, length(cells$cells1), length(cells$cells2))
    }
    universe <- de$gene
  }
  species <- if (use_gmt) NA_character_ else settings$gene_sets
  ontology <- if (use_gmt) NA_character_ else settings$ontology %||% "BP"
  min_size <- settings$min_size %||% 10; max_size <- settings$max_size %||% 500
  max_p_adj <- settings$ora_max_p_adj %||% 0.05; min_log2fc <- settings$ora_min_log2fc %||% 0.25; direction <- settings$direction %||% "up"
  ora <- if ("ORA" %in% methods) run_ora_directions(de, direction, universe, max_p_adj, min_log2fc, species, ontology, term2gene, min_size, max_size) else NULL
  ranked <- if ("GSEA" %in% methods) rank_genes_for_gsea(de, seed = 42) else NULL
  gsea <- if (!is.null(ranked)) run_gsea(ranked, species, ontology, term2gene, min_size, max_size, seed = 42) else NULL
  gene_sets_desc <- if (use_gmt) sprintf("GMT file: %s%s", gmt$file, if (!is.null(gmt$md5)) paste0(" (md5 ", gmt$md5, ")") else "") else orgdb_description(species)
  params <- data.frame(parameter = c("Comparison", "Genes", "Assay", "Gene sets", "Ontology", "Set size", "ORA genes", "ORA universe", "GSEA ranking", "Multiple testing", "Seed", "Software"),
                       value = c(label, source_desc, assay_desc, gene_sets_desc,
                                 if (use_gmt) "-" else ontology, sprintf("%d-%d genes", min_size, max_size),
                                 if (is.null(ora)) "-" else sprintf("%s (adjusted p <= %g, |log2FC| >= %g; each direction tested and corrected separately)", paste(sprintf("%s: %d genes", names(ora$n_genes), ora$n_genes), collapse = ", "), max_p_adj, min_log2fc),
                                 if (is.null(ora)) "-" else sprintf("%d genes", ora$n_universe),
                                 if (is.null(gsea)) "-" else sprintf("%d genes ranked by %s (preranked GSEA, null from random gene sets of the same size; %.1f%% tied values, ties broken at random with seed 42)", gsea$n_ranked, ranking_desc, 100 * attr(ranked, "fraction_tied")),
                                 "Benjamini-Hochberg: ORA over the sets containing at least one query gene (per direction), GSEA over all tested sets",
                                 "42 (GSEA ranking ties and permutations)",
                                 sprintf("ORA: hypergeometric test (stats::phyper); GSEA: fgsea %s (multilevel, eps = 0, nPermSimple = 10000); GO.db %s%s",
                                         as.character(utils::packageVersion("fgsea")), if (requireNamespace("GO.db", quietly = TRUE)) as.character(utils::packageVersion("GO.db")) else "not installed",
                                         if (identical(de_method, "pseudobulk") && !identical(settings$source, "stored")) paste0("; DESeq2 ", as.character(utils::packageVersion("DESeq2"))) else "")),
                       stringsAsFactors = FALSE)
  if (use_gmt && !is.null(gmt$fraction_matched)) {
    params <- rbind(params, data.frame(parameter = "Gene sets matched", value = sprintf("%.1f%% of set genes found in the dataset (case-insensitive symbol matching, not an orthology mapping)", 100 * gmt$fraction_matched), stringsAsFactors = FALSE))
  }
  gene_set_label <- if (use_gmt) paste("Gene sets:", gmt$file) else sprintf("Gene Ontology %s, %s", c(BP = "biological process", MF = "molecular function", CC = "cellular component")[[ontology]], species)
  if (!is.null(ora) && length(ora$notes) > 0) params <- rbind(params, data.frame(parameter = "ORA skipped", value = paste(ora$notes, collapse = "; "), stringsAsFactors = FALSE))
  return(list(ora = ora, gsea = gsea, label = label, params = params, gsea_skipped = gsea_skipped, gene_set_label = gene_set_label))
}

# Entry of the stored-table list (same structure as storeDETable) holding an enrichment run: ORA and GSEA rows in one
# table (column 'method'), the parameters table as a named list, and the number of significant sets per method.
# Stored with type "enrichment" so it is saved with the session and summarized in the parameter export and PDF report.
enrichment_storage_entry <- function(res, max_p_adj = 0.05) {
  cols <- c("method", "direction", "ID", "Description", "GeneRatio", "BgRatio", "fold_enrichment", "Count", "setSize",
            "enrichmentScore", "NES", "pvalue", "p.adjust", "geneID", "core_enrichment")
  pieces <- list()
  for (m in c("ORA", "GSEA")) {
    tab <- res[[tolower(m)]]$table
    if (is.null(tab) || nrow(tab) == 0) next
    tab$method <- m
    # Same columns for both methods so the rows can be bound (NA where a column does not apply)
    for (k in setdiff(cols, colnames(tab))) tab[[k]] <- NA
    pieces[[m]] <- tab[, cols, drop = FALSE]
  }
  if (length(pieces) == 0) stop("No tested term to store")
  data <- do.call(rbind, unname(pieces))
  rownames(data) <- make.unique(paste(data$method, data$direction %||% "", data$ID, sep = "_"))
  n_sig <- vapply(c(ORA = "ORA", GSEA = "GSEA"), function(m) sum(data$method == m & !is.na(data$p.adjust) & data$p.adjust <= max_p_adj), numeric(1))
  return(list(name = NULL, description = paste0("Enrichment: ", res$label, " (", res$gene_set_label, ")"), data = data, type = "enrichment",
              gene_list = character(0), timestamp = Sys.time(), parameters = as.list(stats::setNames(res$params$value, res$params$parameter)),
              n_significant = n_sig, max_p_adj = max_p_adj))
}

# Number of sets below the adjusted p-value cutoff in each result, for the completion message
enrichment_done_message <- function(res, max_p_adj = 0.05) {
  n_sig <- function(r) if (is.null(r) || nrow(r$table) == 0) 0 else sum(r$table$p.adjust <= max_p_adj, na.rm = TRUE)
  return(paste0("Enrichment done: ", if (is.null(res$ora)) "" else sprintf("ORA %d significant sets ", n_sig(res$ora)),
                if (is.null(res$gsea)) "" else sprintf("GSEA %d significant sets", n_sig(res$gsea))))
}

# Text summary of a run (label and parameters)
enrichment_summary_text <- function(res) {
  return(paste0(res$label, "\n", paste(res$params$parameter, res$params$value, sep = ": ", collapse = "\n")))
}

# File name stem of the downloads
enrichment_file_stub <- function(label) {
  return(paste0("enrichment_", gsub("[^A-Za-z0-9]+", "_", label %||% "results"), "_", Sys.Date()))
}

############################## Tables ##############################

# Columns shown and exported for each method
enrichment_table_columns <- function(method = "ORA") {
  if (identical(method, "GSEA")) return(c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "core_enrichment"))
  return(c("direction", "ID", "Description", "GeneRatio", "BgRatio", "fold_enrichment", "Count", "pvalue", "p.adjust", "geneID"))
}

# Signed significance of each term: -log10(adjusted p) with a + sign when the set goes up in group 1 and a - sign when it
# goes down (ORA: gene list of the row; GSEA: sign of the NES). Sorting on it puts the most significant up terms first and
# the most significant down terms last.
enrichment_signed_score <- function(tab, method = "ORA") {
  sgn <- if (identical(method, "GSEA")) ifelse(!is.na(tab$NES) & tab$NES < 0, -1, 1) else ifelse(tab$direction %in% "Down in group 1", -1, 1)
  return(round(sgn * -log10(pmax(tab$p.adjust, .Machine$double.xmin)), 2))
}

# Colour scale of the signed score: white for terms above the adjusted p-value cutoff, then four shades of red (up) or blue
# (down) at cutoff, +2, +5 and +10 units of -log10(adjusted p). Returns the 8 breaks and the 9 background / 9 text colours.
enrichment_score_scale <- function(max_p_adj = 0.05) {
  cut <- -log10(max_p_adj) - 1e-9
  steps <- c(cut, cut + 2, cut + 5, cut + 10)
  return(list(breaks = c(-rev(steps), steps),
              fill = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
              text = c("#FFFFFF", "#FFFFFF", "#000000", "#000000", "#000000", "#000000", "#000000", "#FFFFFF", "#FFFFFF")))
}

# Table of one method with the signed score after the term name, sorted by it (most up first, most down last);
# only terms below the adjusted p-value cutoff unless show_all
enrichment_result_table <- function(tab, method = "ORA", max_p_adj = 0.05, show_all = FALSE) {
  tab <- tab[, intersect(enrichment_table_columns(method), colnames(tab)), drop = FALSE]
  if (!show_all) tab <- tab[!is.na(tab$p.adjust) & tab$p.adjust <= max_p_adj, , drop = FALSE]
  score <- enrichment_signed_score(tab, method)
  pos <- match("Description", colnames(tab))
  tab <- cbind(tab[, seq_len(pos), drop = FALSE], signed_log10_padj = score, tab[, -seq_len(pos), drop = FALSE])
  tab <- tab[order(-tab$signed_log10_padj), , drop = FALSE]
  rownames(tab) <- NULL
  return(tab)
}

# Result table for display: rows shaded by the signed score (darker red = more significant up, darker blue = more
# significant down in group 1, white = not significant); p-values kept numeric and formatted by DT only
enrichment_datatable <- function(tab, method = "ORA", max_p_adj = 0.05, show_all = FALSE) {
  tab <- enrichment_result_table(tab, method, max_p_adj, show_all)
  sc <- enrichment_score_scale(max_p_adj)
  dt <- DT::datatable(tab, rownames = FALSE, filter = "top",
                      colnames = ifelse(colnames(tab) == "signed_log10_padj", "Signed -log10(adj. p)", colnames(tab)),
                      options = list(pageLength = 15, scrollX = TRUE))
  dt <- DT::formatStyle(dt, "signed_log10_padj", target = "row", backgroundColor = DT::styleInterval(sc$breaks, sc$fill), color = DT::styleInterval(sc$breaks, sc$text))
  return(DT::formatSignif(dt, intersect(c("pvalue", "p.adjust", "NES", "enrichmentScore", "fold_enrichment"), colnames(tab)), digits = 3))
}

# Excel export of one method: every tested term sorted by the signed score, rows shaded with the same scale as the app,
# numeric p-values in scientific format, header frozen and filterable, and the parameters of the run in a second sheet
write_enrichment_xlsx <- function(res, method = "ORA", file, max_p_adj = 0.05) {
  tab <- res[[tolower(method)]]$table
  if (is.null(tab) || nrow(tab) == 0) stop("No result to export")
  tab <- enrichment_result_table(tab, method, max_p_adj, show_all = TRUE)
  sc <- enrichment_score_scale(max_p_adj)
  bin <- findInterval(tab$signed_log10_padj, sc$breaks, left.open = TRUE) + 1
  wb <- openxlsx::createWorkbook()
  sheet <- paste(method, "results")
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, tab, headerStyle = openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9D9D9", border = "bottom"), withFilter = TRUE)
  cols <- seq_len(ncol(tab))
  for (b in setdiff(unique(bin), 5)) {
    openxlsx::addStyle(wb, sheet, openxlsx::createStyle(fgFill = sc$fill[b], fontColour = sc$text[b]), rows = which(bin == b) + 1, cols = cols, gridExpand = TRUE, stack = TRUE)
  }
  num_cols <- which(colnames(tab) %in% c("pvalue", "p.adjust"))
  if (length(num_cols) > 0) openxlsx::addStyle(wb, sheet, openxlsx::createStyle(numFmt = "0.00E+00"), rows = seq_len(nrow(tab)) + 1, cols = num_cols, gridExpand = TRUE, stack = TRUE)
  openxlsx::freezePane(wb, sheet, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet, cols = cols, widths = ifelse(colnames(tab) %in% c("Description", "geneID", "core_enrichment"), 45, 14))
  openxlsx::addWorksheet(wb, "Parameters")
  openxlsx::writeData(wb, "Parameters", res$params, headerStyle = openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9D9D9"))
  legend <- data.frame(parameter = c("signed_log10_padj", "Row colours"),
                       value = c("-log10(adjusted p), positive when the set goes up in group 1, negative when it goes down",
                                 sprintf("White = adjusted p > %g; red = up, blue = down in group 1, darker = more significant", max_p_adj)))
  openxlsx::writeData(wb, "Parameters", legend, startRow = nrow(res$params) + 3, colNames = FALSE)
  openxlsx::setColWidths(wb, "Parameters", cols = 1:2, widths = c(25, 100))
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  return(invisible(file))
}

# Plot of a result with errors (e.g. no significant term) turned into a message shown in the plot area
enrichment_plot_or_message <- function(table, method, max_p_adj, n_terms, title, subtitle = NULL) {
  return(tryCatch(plot_enrichment_terms(table, method, max_p_adj, n_terms, title, subtitle),
                  error = function(e) ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = conditionMessage(e), size = 5, colour = "grey40") + ggplot2::theme_void()))
}

############################## Plots ##############################

# Dot plot of the most significant terms (adjusted p-value below the cutoff), most significant at the top.
# Term labels carry the set identifier (e.g. GO:0007517). Colour = direction (red = up, blue = down in group 1),
# shade = -log10 adjusted p-value (darker = more significant).
# ORA: x = fold enrichment (fraction of the significant genes in the term / fraction of all tested genes in the term),
#      size = significant genes in the term; one panel per direction (up / down in group 1), n_terms per panel.
# GSEA: x = normalized enrichment score (NES > 0: genes of the term at the top of the ranking, i.e. up in group 1;
#       NES < 0: down in group 1), size = term size.
plot_enrichment_terms <- function(table, method = "ORA", max_p_adj = 0.05, n_terms = 15, title = NULL, subtitle = NULL) {
  if (is.null(table) || nrow(table) == 0) stop("No term was tested")
  # Sets whose p-value fgsea could not compute (NA) are left out
  top <- table[!is.na(table$p.adjust) & table$p.adjust <= max_p_adj, , drop = FALSE]
  if (nrow(top) == 0) stop(sprintf("No term with adjusted p-value <= %g", max_p_adj))
  if (identical(method, "GSEA")) top$direction <- "All tested genes"
  if (!"direction" %in% colnames(top)) top$direction <- "Significant genes"
  top$direction <- factor(top$direction, levels = intersect(c("Up in group 1", "Down in group 1", "Significant genes", "All tested genes"), unique(top$direction)))
  # Top terms of each panel
  top <- do.call(rbind, lapply(split(top, top$direction, drop = TRUE), function(d) d[order(d$p.adjust, d$pvalue %||% d$p.adjust), , drop = FALSE][seq_len(min(n_terms, nrow(d))), , drop = FALSE]))
  if (identical(method, "GSEA")) {
    top$x <- top$NES; top$size <- top$setSize; x_label <- "Normalized enrichment score (NES)"; size_label <- "Genes in term"
  } else {
    top$x <- top$fold_enrichment %||% top$gene_ratio; top$size <- top$Count; x_label <- "Fold enrichment"; size_label <- "Significant genes"
  }
  top$score <- -log10(pmax(top$p.adjust, .Machine$double.xmin))
  # Colour = direction (red up / blue down in group 1; GSEA: sign of the NES), shade = significance
  top$tone <- if (identical(method, "GSEA")) ifelse(top$NES >= 0, "Up in group 1", "Down in group 1") else ifelse(top$direction == "Down in group 1", "Down in group 1", "Up in group 1")
  # Identifier added to the term name when it differs from it (GO terms; GMT sets have the name as identifier)
  label <- stringr::str_wrap(ifelse(top$ID == top$Description, top$Description, paste0(top$Description, " (", top$ID, ")")), 55)
  # Panel-specific y key so a term found in both directions keeps its own position in each panel
  top$key <- factor(paste0(as.integer(top$direction), "___", label), levels = rev(paste0(as.integer(top$direction), "___", label)))
  p <- ggplot2::ggplot(top, ggplot2::aes(x = x, y = key, size = size, colour = tone, alpha = score)) + ggplot2::geom_point() +
    ggplot2::scale_y_discrete(labels = function(k) sub("^[0-9]+___", "", k)) +
    ggplot2::scale_colour_manual(values = c(`Up in group 1` = "#B2182B", `Down in group 1` = "#2166AC"), name = "Direction",
                                 guide = if (identical(method, "GSEA")) ggplot2::guide_legend(override.aes = list(size = 4), order = 1) else "none") +
    ggplot2::scale_alpha_continuous(range = c(0.3, 1), name = "-log10\nadjusted p", guide = ggplot2::guide_legend(override.aes = list(colour = "grey20", size = 4), order = 2)) +
    ggplot2::scale_size_continuous(name = size_label, range = c(3, 8)) +
    ggplot2::labs(x = x_label, y = NULL, title = title, subtitle = subtitle) + ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey30"),
                   strip.text = ggplot2::element_text(face = "bold", size = 12), plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  if (nlevels(droplevels(top$direction)) > 1) p <- p + ggplot2::facet_wrap(~direction, ncol = 2, scales = "free_y")
  p <- if (identical(method, "GSEA")) p + ggplot2::geom_vline(xintercept = 0, colour = "grey60", linetype = "dashed") else p + ggplot2::geom_vline(xintercept = 1, colour = "grey60", linetype = "dashed") + ggplot2::expand_limits(x = 1)
  return(p)
}