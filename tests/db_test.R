# db_test.R — CellHub Database checks inside the Cell-Hub image, every species and every source file, with the app's own functions.
# Usage (inside the image): Rscript /tests/db_test.R /srv/shiny-server/mon_app /out
args <- commandArgs(trailingOnly = TRUE)
app_dir <- if (length(args) >= 1) args[1] else "/srv/shiny-server/mon_app"
out_dir <- if (length(args) >= 2) args[2] else "/out"
suppressMessages({ library(shiny); library(CellChat) })
# App helpers: dev names or _docker names, whichever the image ships
src <- function(stem) {
  f <- c(file.path(app_dir, "functions", paste0(stem, "_docker.R")), file.path(app_dir, "functions", paste0(stem, ".R")))
  f <- f[file.exists(f)][1]
  if (is.na(f)) stop("helper file not found: ", stem)
  source(f, local = globalenv())
}
`%||%` <- function(a, b) if (is.null(a)) b else a
src("cellchat_helpers"); src("cellhubdb_filter_functions")
db_root <- file.path(app_dir, "databases")
species <- sort(list.dirs(db_root, full.names = FALSE, recursive = FALSE))
rows <- list(); missing_genes <- list()
add <- function(sp, check, status, detail = "", n = NA) rows[[length(rows) + 1]] <<- data.frame(species = sp, check = check, status = status, n = n, detail = substr(detail, 1, 600), stringsAsFactors = FALSE)
safe <- function(sp, check, expr) tryCatch(expr, error = function(e) { add(sp, check, "KO", conditionMessage(e)); NULL })

# Content checks of one CellChat-format database; returns nothing, records rows
check_db <- function(db, sp, label) {
  it <- db$interaction
  add(sp, paste(label, "interactions"), if (nrow(it) > 0) "OK" else "KO", n = nrow(it))
  empty <- sum(is.na(it$ligand) | it$ligand == "" | is.na(it$receptor) | it$receptor == "")
  add(sp, paste(label, "empty ligand / receptor"), if (empty == 0) "OK" else "KO", n = empty)
  dup <- sum(duplicated(it$interaction_name))
  add(sp, paste(label, "duplicated interaction_name"), if (dup == 0) "OK" else "KO", n = dup)
  no_pw <- sum(is.na(it$pathway_name) | it$pathway_name == "")
  add(sp, paste(label, "interactions without pathway"), if (no_pw == 0) "OK" else "WARN", n = no_pw)
  # Multi-subunit partners must be defined in the complex table (CellChat silently drops undefined ones)
  partners <- unique(c(it$ligand, it$receptor))
  multi <- partners[grepl("_", partners)]
  undefined <- setdiff(multi, rownames(db$complex))
  add(sp, paste(label, "multi-subunit partners without complex definition"), if (length(undefined) == 0) "OK" else "KO", n = length(undefined), paste(head(undefined, 15), collapse = ", "))
  # Every partner gene (complex subunits expanded) must be in geneInfo, or CellChat drops the interaction
  sub_cols <- grep("^subunit", colnames(db$complex), value = TRUE)
  # Non-protein ligands (metabolites, small molecules) are not genes: they are expected to be absent from geneInfo
  non_protein <- if ("interaction_type" %in% colnames(it)) unique(it$ligand[it$interaction_type %in% "Non-protein Signaling"]) else character(0)
  partners <- setdiff(partners, non_protein)
  genes <- unique(c(partners[!partners %in% rownames(db$complex)], unlist(db$complex[intersect(partners, rownames(db$complex)), sub_cols])))
  genes <- genes[!is.na(genes) & genes != ""]
  miss <- setdiff(genes, db$geneInfo$Symbol)
  pct <- round(100 * length(miss) / max(1, length(genes)), 2)
  add(sp, paste(label, "partner genes absent from geneInfo"), if (length(miss) == 0) "OK" else if (pct < 5) "WARN" else "KO", n = length(miss), sprintf("%.2f%% of %d genes: %s", pct, length(genes), paste(head(miss, 15), collapse = ", ")))
  if (length(miss) > 0) missing_genes[[length(missing_genes) + 1]] <<- data.frame(species = sp, database = label, gene = miss, stringsAsFactors = FALSE)
}

for (sp in species) {
  sp_dir <- file.path(db_root, sp)
  merged_path <- file.path(sp_dir, paste0("CellHubDB_", sp, "_merged.rds"))
  t0 <- Sys.time()
  merged <- safe(sp, "merged file loads", readRDS(merged_path))
  if (is.null(merged)) next
  add(sp, "merged file loads", "OK", sprintf("%.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  check_db(merged, sp, "merged")
  # Confidence tiers as the app reads them, and every tier combination offered in the interface
  merged$interaction <- annotate_cellhubdb_tier(merged$interaction)
  tc <- cellhubdb_tier_counts(merged$interaction)
  add(sp, "confidence tiers present", if (is.null(tc)) "KO" else "OK", detail = if (is.null(tc)) "no confidence_tier column" else paste(names(tc), tc, collapse = " | "))
  na_tier <- if ("confidence_tier" %in% colnames(merged$interaction)) sum(is.na(merged$interaction$confidence_tier)) else NA
  add(sp, "interactions without tier", if (isTRUE(na_tier == 0)) "OK" else "KO", n = na_tier)
  combos <- list(list(t = c("high", "medium"), r = TRUE), list(t = c("high", "medium"), r = FALSE), list(t = c("high", "medium", "low"), r = FALSE),
                 list(t = "high", r = FALSE), list(t = "medium", r = TRUE), list(t = "medium", r = FALSE), list(t = "low", r = FALSE))
  for (cb in combos) {
    tiers <- cb$t
    label <- paste("tiers", paste(tiers, collapse = "+"), if (cb$r) "(medium roles confirmed)" else "")
    res <- tryCatch(filter_cellhubdb_tiers(merged, tiers, medium_roles_only = cb$r), error = function(e) e)
    if (inherits(res, "error")) { add(sp, label, "INFO", conditionMessage(res)); next }
    prep <- tryCatch(prepare_cellchat_db(res$db), error = function(e) e)
    # A failure is a defect only for the default combination of the species (the one the interface pre-selects)
    is_default <- setequal(tiers, cellhubdb_recommended_tiers(sp)) && identical(cb$r, cellhubdb_recommended_roles_only(sp))
    if (inherits(prep, "error")) add(sp, paste(label, "-> CellChat-ready"), if (is_default) "KO" else "WARN", conditionMessage(prep), n = res$n_after)
    else add(sp, paste(label, "-> CellChat-ready"), "OK", sprintf("%d interactions, %d pathways", nrow(prep$db$interaction), length(unique(prep$db$interaction$pathway_name))), n = res$n_after)
  }
  # Per-source files: each alone, then all merged with the app's merge, tiers recomputed on the selected sources
  src_files <- list.files(sp_dir, pattern = paste0("^CellHubDB_", sp, "_.*\\.rds$"))
  src_files <- src_files[!grepl("_merged\\.rds$", src_files)]
  src_names <- sub("\\.rds$", "", sub(paste0("^CellHubDB_", sp, "_"), "", src_files))
  db_list <- list()
  for (s in src_names) {
    d <- safe(sp, paste("source", s, "loads"), readRDS(file.path(sp_dir, paste0("CellHubDB_", sp, "_", s, ".rds"))))
    if (is.null(d)) next
    db_list[[s]] <- d
    check_db(d, sp, paste("source", s))
    prep <- tryCatch(prepare_cellchat_db(d), error = function(e) e)
    add(sp, paste("source", s, "-> CellChat-ready"), if (inherits(prep, "error")) "WARN" else "OK", if (inherits(prep, "error")) conditionMessage(prep) else "", n = nrow(d$interaction))
  }
  if (length(db_list) >= 2) {
    mr <- safe(sp, "app merge of all sources", merge_cellhubdb_sources(db_list))
    if (!is.null(mr)) {
      add(sp, "app merge of all sources", if (abs(mr$stats$n_merged - nrow(merged$interaction)) <= 0.01 * nrow(merged$interaction)) "OK" else "WARN",
          sprintf("app merge %d vs release merged %d", mr$stats$n_merged, nrow(merged$interaction)), n = mr$stats$n_merged)
      unmatched <- sum(is.na(match_interaction_rows(mr$db$interaction, merged$interaction)))
      add(sp, "app-merged pairs not found in the release merged file", if (unmatched == 0) "OK" else "WARN", n = unmatched)
      ann <- annotate_cellhubdb_tier(mr$db$interaction, merged$interaction, selected_sources = names(db_list))
      tc2 <- cellhubdb_tier_counts(ann)
      add(sp, "tiers recomputed on all sources vs release tiers", if (!is.null(tc2) && all(abs(tc2 - tc) <= 0.02 * sum(tc))) "OK" else "WARN",
          paste("recomputed:", paste(names(tc2), tc2, collapse = " | "), "/ release:", paste(names(tc), tc, collapse = " | ")))
      cur <- unique(trimws(unlist(strsplit(stats::na.omit(merged$interaction$curated_sources), "\\s*[,;|]\\s*"))))
      cur <- cur[nzchar(cur)]
      add(sp, "curated_sources names match the per-source files", if (all(cur %in% names(db_list))) "OK" else "KO",
          paste("not a file name:", paste(setdiff(cur, names(db_list)), collapse = ", ")))
    }
  }
}

# Reviewer example: EGFR partners in human, with the tiers and corroborated roles
hum <- tryCatch(readRDS(file.path(db_root, "human", "CellHubDB_human_merged.rds")), error = function(e) NULL)
if (!is.null(hum)) {
  it <- annotate_cellhubdb_tier(hum$interaction)
  e <- it[it$ligand == "EGFR" | it$receptor == "EGFR" | grepl("(^|_)EGFR(_|$)", it$ligand) | grepl("(^|_)EGFR(_|$)", it$receptor), ]
  add("human", "EGFR partners (tier high/medium/low)", "INFO", paste(names(table(e$confidence_tier)), table(e$confidence_tier), collapse = " | "), n = nrow(e))
  for (g in c("CALM2", "ARF4", "GSTP1", "HSP90AA1")) {
    r <- e[e$ligand == g | e$receptor == g, ]
    add("human", paste("EGFR -", g), "INFO", if (nrow(r) == 0) "absent" else paste0(r$confidence_tier, ", role_corroborated=", r$role_corroborated %||% NA, ", sources=", r$sources, collapse = "; "))
  }
}
res <- do.call(rbind, rows)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(res, file.path(out_dir, "db_results.csv"), row.names = FALSE)
if (length(missing_genes) > 0) write.csv(do.call(rbind, missing_genes), file.path(out_dir, "db_partner_genes_missing.csv"), row.names = FALSE)
cat(sprintf("DB checks: %d rows, %d KO, %d WARN\n", nrow(res), sum(res$status == "KO"), sum(res$status == "WARN")))
