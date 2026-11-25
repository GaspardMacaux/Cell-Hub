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
