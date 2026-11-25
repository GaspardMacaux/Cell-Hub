# cellhubdb_filter_functions.R
# Confidence-tier filter applied to a CellHub Database (CellChat-native list: interaction, complex, cofactor, geneInfo)
# before CellChat runs. Pure functions (no reactive code).

# Species whose databases carry the role_corroborated flag (gene-level: the ligand gene is used as a ligand and the
# receptor gene as a receptor by 2 or more curated sources). Only human and mouse have enough curated sources for it.
cellhubdb_species_with_roles <- function() c("human", "mouse")

# Recommended tiers per species. Human and mouse: high + medium, medium restricted to confirmed ligand / receptor roles
# (the recommended set of the CellHubDB paper). Other species rely on few sources, mostly orthology-inferred or predicted:
# every tier is kept and results are read as hypotheses.
cellhubdb_recommended_tiers <- function(species) {
  if (species %in% cellhubdb_species_with_roles()) c("high", "medium") else c("high", "medium", "low")
}

# Whether medium interactions are restricted to confirmed roles by default for this species
cellhubdb_recommended_roles_only <- function(species) {
  species %in% cellhubdb_species_with_roles()
}

# Row of `reference` matching each row of `interaction`: case-insensitive ligand / receptor pair, then the reversed pair
# (merge_cellhubdb_sources flips non-oriented pairs and rebuilds interaction_name, so names alone do not match), then interaction_name
match_interaction_rows <- function(interaction, reference) {
  key <- function(a, b) paste(tolower(a), tolower(b), sep = "|||")
  ref_key <- key(reference$ligand, reference$receptor)
  idx <- match(key(interaction$ligand, interaction$receptor), ref_key)
  miss <- is.na(idx)
  idx[miss] <- match(key(interaction$receptor[miss], interaction$ligand[miss]), ref_key)
  miss <- is.na(idx)
  if (any(miss) && "interaction_name" %in% colnames(reference)) idx[miss] <- match(interaction$interaction_name[miss], reference$interaction_name)
  idx
}

#' Confidence tier of every interaction
#'
#' The merged release file carries confidence_tier, curated_sources and role_corroborated. Per-source files and databases
#' merged in the app do not: both annotation columns are joined from the merged file of the species, and the tier is
#' recomputed on the selected sources only (high = curated support in >= 2 selected sources, medium = 1, low = 0), so a
#' pair is never "high" because of a source the user left out. role_corroborated keeps its release value (gene-level flag).
#' @param merged_interaction Interaction table of the merged release file (NULL when the table already has the columns)
#' @param selected_sources Selected source names; NULL keeps the release tier
#' @return The interaction table with a lower-case confidence_tier column (NA when unknown)
annotate_cellhubdb_tier <- function(interaction, merged_interaction = NULL, selected_sources = NULL) {
  if (!is.null(merged_interaction)) {
    missing_cols <- setdiff(intersect(c("curated_sources", "role_corroborated"), colnames(merged_interaction)), colnames(interaction))
    if (length(missing_cols) > 0) {
      idx <- match_interaction_rows(interaction, merged_interaction)
      for (col in missing_cols) interaction[[col]] <- merged_interaction[[col]][idx]
    }
  }
  if (!is.null(selected_sources) && "curated_sources" %in% colnames(interaction)) {
    curated <- strsplit(ifelse(is.na(interaction$curated_sources), "", interaction$curated_sources), "\\s*[,;|]\\s*")
    n_sel <- vapply(curated, function(v) sum(v[nzchar(v)] %in% selected_sources), integer(1))
    interaction$confidence_tier <- ifelse(is.na(interaction$curated_sources), NA_character_, ifelse(n_sel >= 2, "high", ifelse(n_sel == 1, "medium", "low")))
  }
  if ("confidence_tier" %in% colnames(interaction)) interaction$confidence_tier <- tolower(as.character(interaction$confidence_tier))
  interaction
}

#' Keep the interactions of the selected confidence tiers
#' @param tiers Tiers kept among "high", "medium", "low"; all three keep every interaction (including unknown tiers)
#' @param medium_roles_only Keep only the medium interactions whose ligand / receptor roles are confirmed (role_corroborated);
#'   ignored when the database has no role_corroborated column
#' @return list(db, n_before, n_after, tiers, medium_roles_only)
filter_cellhubdb_tiers <- function(db, tiers, medium_roles_only = FALSE) {
  tiers <- tolower(tiers)
  if (length(tiers) == 0) stop("Select at least one confidence tier.")
  inter <- db$interaction
  n_before <- nrow(inter)
  has_roles <- "role_corroborated" %in% colnames(inter)
  roles_applied <- isTRUE(medium_roles_only) && has_roles && "medium" %in% tiers
  if (!all(c("high", "medium", "low") %in% tiers) || roles_applied) {
    if (!"confidence_tier" %in% colnames(inter)) stop("This database has no confidence tiers (CellHub Database 2.1 or later needed): keep all three tiers.")
    keep <- !is.na(inter$confidence_tier) & inter$confidence_tier %in% tiers
    if (all(c("high", "medium", "low") %in% tiers)) keep <- keep | is.na(inter$confidence_tier)
    if (roles_applied) keep <- keep & (is.na(inter$confidence_tier) | inter$confidence_tier != "medium" | as.logical(inter$role_corroborated) %in% TRUE)
    keep[is.na(keep)] <- FALSE
    if (!any(keep)) stop(sprintf("No interaction of tier %s in this database: select more tiers.", paste(tiers, collapse = " + ")))
    db$interaction <- inter[keep, , drop = FALSE]
  }
  list(db = db, n_before = n_before, n_after = nrow(db$interaction), tiers = tiers, medium_roles_only = roles_applied)
}

# Number of interactions per confidence tier, medium split by confirmed roles when the database has the flag
cellhubdb_tier_counts <- function(interaction) {
  if (!"confidence_tier" %in% colnames(interaction)) return(NULL)
  tier <- tolower(as.character(interaction$confidence_tier))
  counts <- c(high = sum(tier == "high", na.rm = TRUE), medium = sum(tier == "medium", na.rm = TRUE))
  if ("role_corroborated" %in% colnames(interaction)) {
    counts <- c(counts, "medium with confirmed roles" = sum(tier == "medium" & as.logical(interaction$role_corroborated) %in% TRUE, na.rm = TRUE))
  }
  c(counts, low = sum(tier == "low", na.rm = TRUE))
}

# Tier selection shown under the database selection of the CellChat module
# (input ids: db_filter_tiers_cellchat, db_medium_roles_cellchat)
cellhubdb_filter_panel <- function() {
  tags$div(
    style = "margin-top: 10px;",
    checkboxGroupInput("db_filter_tiers_cellchat", "Confidence tiers to keep:",
                       choices = c("High" = "high", "Medium" = "medium", "Low" = "low"), selected = c("high", "medium"), inline = TRUE),
    checkboxInput("db_medium_roles_cellchat", "Medium: confirmed ligand / receptor roles only", value = TRUE),
    tags$p(tags$b("High:"), " supported by 2 or more curated sources. ", tags$b("Medium:"), " 1 curated source; ",
           tags$i("confirmed roles"), " = the ligand gene is used as a ligand and the receptor gene as a receptor by 2 or more curated sources (human and mouse). ",
           tags$b("Low:"), " orthology-inferred or predicted only. ",
           "Default: high + medium with confirmed roles for human and mouse, all tiers for other species.",
           style = "color: rgba(255,255,255,0.75); font-size: 11px; margin: -4px 0 6px 0; line-height: 1.35;")
  )
}

# Genes typed or pasted by the user: separated by spaces, new lines, commas, semicolons or tabs; duplicates removed (case-insensitive)
parse_gene_list <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) return(character(0))
  genes <- trimws(unlist(strsplit(text, "[\\s,;]+", perl = TRUE)))
  genes <- genes[nzchar(genes)]
  genes[!duplicated(toupper(genes))]
}

# Genes of each partner (ligand or receptor): the name itself plus its complex subunits (complex table of the database),
# or its "_"-separated parts when the complex is not in the table
partner_gene_sets <- function(partners, complex = NULL) {
  has_cpx <- !is.null(complex) && nrow(complex) > 0
  lapply(as.character(partners), function(p) {
    if (has_cpx && p %in% rownames(complex)) {
      g <- unlist(complex[p, ], use.names = FALSE)
      return(unique(c(p, g[!is.na(g) & nzchar(g)])))
    }
    unique(c(p, strsplit(p, "_", fixed = TRUE)[[1]]))
  })
}

#' Rows whose ligand or receptor (or one of their complex subunits) is in a gene list
#' @param interaction Database interaction table or CellChat communications (columns ligand and receptor)
#' @param complex Complex table of the database (rownames = complex names), or NULL
#' @param genes Gene symbols, matched case-insensitively
#' @return list(keep = logical per row, matched = list genes found in each row (";"-separated), genes_found, genes_absent)
match_interactions_to_genes <- function(interaction, complex, genes) {
  genes_up <- toupper(genes)
  if (is.null(interaction) || nrow(interaction) == 0) return(list(keep = logical(0), matched = character(0), genes_found = character(0), genes_absent = genes))
  lig <- partner_gene_sets(interaction$ligand, complex)
  rec <- partner_gene_sets(interaction$receptor, complex)
  matched <- vapply(seq_along(lig), function(i) {
    hit <- genes_up[genes_up %in% toupper(c(lig[[i]], rec[[i]]))]
    paste(genes[match(hit, genes_up)], collapse = ";")
  }, character(1))
  keep <- nzchar(matched)
  found_up <- unique(toupper(unlist(strsplit(matched[keep], ";", fixed = TRUE))))
  list(keep = keep, matched = matched, genes_found = genes[genes_up %in% found_up], genes_absent = genes[!genes_up %in% found_up])
}

#' Interactions of a loaded CellHub Database that involve a gene list (query of the database, no CellChat run)
#' @param db CellChat-native database (interaction, complex, cofactor, geneInfo), after the tier filter
#' @param genes Gene symbols (parse_gene_list output)
#' @return list(table, genes_found, genes_absent, n_genes, n_db)
query_cellhubdb_genes <- function(db, genes) {
  if (length(genes) == 0) stop("Enter at least one gene.")
  inter <- db$interaction
  if (is.null(inter) || nrow(inter) == 0) stop("The loaded database has no interaction.")
  m <- match_interactions_to_genes(inter, db$complex, genes)
  cols <- intersect(c("interaction_name", "ligand", "receptor", "pathway_name", "annotation", "confidence_tier", "role_corroborated", "sources", "evidence"), colnames(inter))
  table <- data.frame(query_genes = m$matched[m$keep], inter[m$keep, cols, drop = FALSE], stringsAsFactors = FALSE, check.names = FALSE)
  rownames(table) <- NULL
  list(table = table, genes_found = m$genes_found, genes_absent = m$genes_absent, n_genes = length(genes), n_db = nrow(inter))
}

# Choices of a gene list selector: "Name (n genes)" labels, values = list names; none_label adds an empty first choice
gene_list_choices <- function(gene_lists, none_label = "None (all interactions)") {
  choices <- if (length(gene_lists) > 0) setNames(names(gene_lists), sprintf("%s (%d genes)", names(gene_lists), lengths(gene_lists))) else character(0)
  if (!is.null(none_label)) choices <- c(setNames("", none_label), choices)
  choices
}

#' Add (or replace) a named gene list
#' @param gene_lists Named list of character vectors (current lists)
#' @param name List name typed by the user; an existing name is replaced
#' @param text Genes pasted by the user (parse_gene_list rules)
#' @return list(lists, name, genes, replaced)
add_gene_list <- function(gene_lists, name, text) {
  name <- if (is.null(name)) "" else trimws(name)
  if (!nzchar(name)) stop("Give the list a name")
  genes <- parse_gene_list(text)
  if (length(genes) == 0) stop("Paste at least one gene symbol")
  replaced <- name %in% names(gene_lists)
  gene_lists[[name]] <- genes
  list(lists = gene_lists, name = name, genes = genes, replaced = replaced)
}

# One row per gene list for the summary table: name, number of genes, first genes
gene_lists_table <- function(gene_lists) {
  if (length(gene_lists) == 0) return(data.frame(List = character(0), Genes = integer(0), First_genes = character(0)))
  data.frame(List = names(gene_lists), Genes = lengths(gene_lists),
             First_genes = vapply(gene_lists, function(g) paste0(paste(head(g, 8), collapse = ", "), if (length(g) > 8) ", ..." else ""), character(1)),
             row.names = NULL, stringsAsFactors = FALSE, check.names = FALSE)
}

# Display name of a CellHub Database source in the interface. The release files keep their labels (sources column,
# per-source file names); the interface adds the name used by the source authors where the two differ
cellhubdb_source_label <- function(src) {
  labels <- c(MultiNicheNet = "NicheNet v2 (MultiNicheNet)")
  out <- unname(labels[src])
  ifelse(is.na(out), src, out)
}