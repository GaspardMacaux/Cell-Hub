# db_deep_check.R — CellHub Database 2.1: every number the manuscript can cite, recomputed from the release files.
# Usage: Rscript db_deep_check.R <release_dir> <out_dir>
#   release_dir: folder holding CellHubDB.rds and app_databases/<species>/CellHubDB_<species>_{merged|<source>}.rds
# Output (out_dir): deep_summary.md (per-species tables to compare with the manuscript), deep_metrics.csv (long format),
#   and detail tables: missing_symbols.csv, complex_issues.csv, symbol_convention_issues.csv, reciprocal_pairs.csv,
#   duplicates.csv, egfr_human.csv
args <- commandArgs(trailingOnly = TRUE)
release_dir <- if (length(args) >= 1) args[1] else "/release"
out_dir <- if (length(args) >= 2) args[2] else "/out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a
app_root <- file.path(release_dir, "app_databases")
release <- readRDS(file.path(release_dir, "CellHubDB.rds"))
species <- sort(list.dirs(app_root, full.names = FALSE, recursive = FALSE))
metrics <- list(); details <- list()
put <- function(sp, metric, value, scope = "merged") metrics[[length(metrics) + 1]] <<- data.frame(species = sp, scope = scope, metric = metric, value = as.character(value), stringsAsFactors = FALSE)
keep <- function(name, df) if (!is.null(df) && nrow(df) > 0) details[[name]] <<- rbind(details[[name]], df)
col_or <- function(df, col, default = NA) if (col %in% colnames(df)) df[[col]] else rep(default, nrow(df))
# Symbol nomenclature per species, as described in the manuscript (HGNC upper case, MGI capitalized, ZFIN lower case)
upper_species <- c("human", "chimp", "macaque", "marmoset", "cow", "dog", "horse", "pig", "sheep", "chicken")
rodent_species <- c("mouse", "rat")
lower_species <- c("zebrafish")

# Partner genes of an interaction table: complex names expanded to their subunits
partner_genes <- function(partners, complex) {
  sub_cols <- grep("^subunit", colnames(complex), value = TRUE)
  in_cx <- partners %in% rownames(complex)
  genes <- c(partners[!in_cx], unlist(complex[partners[in_cx], sub_cols, drop = FALSE]))
  unique(genes[!is.na(genes) & genes != ""])
}

# Interactions whose every partner gene is in geneInfo (CellChat can score them)
scorable <- function(it, complex, symbols) {
  sub_cols <- grep("^subunit", colnames(complex), value = TRUE)
  ok_partner <- function(p) {
    g <- if (p %in% rownames(complex)) unlist(complex[p, sub_cols]) else p
    g <- g[!is.na(g) & g != ""]
    length(g) > 0 && all(g %in% symbols)
  }
  partners <- unique(c(it$ligand, it$receptor))
  ok <- vapply(partners, ok_partner, logical(1))
  ok[it$ligand] & ok[it$receptor]
}

analyze_db <- function(db, sp, scope) {
  it <- db$interaction; cx <- db$complex; symbols <- db$geneInfo$Symbol
  n <- nrow(it)
  put(sp, "interactions", n, scope)
  put(sp, "pathways (distinct, non-empty)", length(unique(na.omit(it$pathway_name[it$pathway_name != ""]))), scope)
  put(sp, "% without pathway", round(100 * mean(is.na(it$pathway_name) | it$pathway_name == ""), 1), scope)
  # Interaction type: the curated interaction_type column (four canonical categories); the CellChat annotation column is
  # reported separately because standalone CellChat users filter on it (subsetDB(search = "Secreted Signaling"))
  ann <- if ("interaction_type" %in% colnames(it)) it$interaction_type else col_or(it, "annotation")
  canonical <- c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor", "Non-protein Signaling")
  put(sp, "% CellChat annotation column with a canonical category", round(100 * mean(col_or(it, "annotation") %in% canonical), 1), scope)
  for (a in sort(unique(na.omit(ann)))) put(sp, paste0("type: ", a), sum(ann == a, na.rm = TRUE), scope)
  put(sp, "type: missing", sum(is.na(ann) | ann == ""), scope)
  # Non-protein ligands (hormones, ions, metabolites): not genes, expected outside geneInfo
  non_protein <- grepl("non-protein", ann, ignore.case = TRUE)
  put(sp, "non-protein signaling interactions", sum(non_protein), scope)
  # Scorable interactions, protein interactions only
  sc <- scorable(it, cx, symbols)
  put(sp, "% scorable (all partner genes in geneInfo), protein interactions", round(100 * mean(sc[!non_protein]), 1), scope)
  genes_protein <- partner_genes(unique(c(it$ligand[!non_protein], it$receptor[!non_protein])), cx)
  miss <- setdiff(genes_protein, symbols)
  put(sp, "partner genes (protein interactions)", length(genes_protein), scope)
  put(sp, "partner genes absent from geneInfo", length(miss), scope)
  if (length(miss) > 0) {
    kind <- ifelse(grepl("^LOC[0-9]+$", miss), "NCBI LOC id", ifelse(grepl("^ENS[A-Z]*G[0-9]+$", miss), "Ensembl id",
                   ifelse(grepl("/", miss), "two ids joined by '/'", "symbol")))
    keep("missing_symbols", data.frame(species = sp, scope = scope, gene = miss, kind = kind, stringsAsFactors = FALSE))
    for (k in unique(kind)) put(sp, paste0("absent from geneInfo: ", k), sum(kind == k), scope)
  }
  # Duplicates: identifiers, and the same oriented pair listed twice
  put(sp, "duplicated interaction_name", sum(duplicated(it$interaction_name)), scope)
  key <- paste(toupper(it$ligand), toupper(it$receptor), sep = "->")
  dup_pairs <- key[duplicated(key)]
  put(sp, "duplicated oriented ligand->receptor pairs", length(dup_pairs), scope)
  if (length(dup_pairs) > 0) keep("duplicates", data.frame(species = sp, scope = scope, pair = unique(dup_pairs), stringsAsFactors = FALSE))
  # Reciprocal pairs (A->B and B->A both present): legitimate for directed signaling, reported for review
  rev_key <- paste(toupper(it$receptor), toupper(it$ligand), sep = "->")
  recip <- key %in% rev_key & toupper(it$ligand) != toupper(it$receptor)
  put(sp, "interactions in reciprocal pairs", sum(recip), scope)
  if (any(recip)) keep("reciprocal_pairs", data.frame(species = sp, scope = scope, ligand = it$ligand[recip], receptor = it$receptor[recip],
                                                      sources = col_or(it, "sources")[recip], stringsAsFactors = FALSE))
  # Complexes: undefined multi-subunit partners, definitions with subunits outside geneInfo, empty definitions
  partners <- unique(c(it$ligand, it$receptor))
  used_cx <- intersect(partners, rownames(cx))
  undefined <- setdiff(partners[grepl("_", partners)], rownames(cx))
  sub_cols <- grep("^subunit", colnames(cx), value = TRUE)
  cx_used <- cx[used_cx, sub_cols, drop = FALSE]
  n_sub <- apply(cx_used, 1, function(r) sum(!is.na(r) & r != ""))
  bad_sub <- apply(cx_used, 1, function(r) { r <- r[!is.na(r) & r != ""]; paste(setdiff(r, symbols), collapse = ", ") })
  put(sp, "complexes used", length(used_cx), scope)
  put(sp, "multi-subunit partners without definition", length(undefined), scope)
  put(sp, "complexes with < 2 subunits", sum(n_sub < 2), scope)
  put(sp, "complexes with subunits absent from geneInfo", sum(bad_sub != ""), scope)
  cx_issues <- rbind(if (length(undefined)) data.frame(species = sp, scope = scope, complex = undefined, issue = "no definition", detail = "", stringsAsFactors = FALSE),
                     if (any(n_sub < 2)) data.frame(species = sp, scope = scope, complex = names(n_sub)[n_sub < 2], issue = "fewer than 2 subunits", detail = "", stringsAsFactors = FALSE),
                     if (any(bad_sub != "")) data.frame(species = sp, scope = scope, complex = names(bad_sub)[bad_sub != ""], issue = "subunit absent from geneInfo", detail = bad_sub[bad_sub != ""], stringsAsFactors = FALSE))
  keep("complex_issues", cx_issues)
  # Cofactors referenced by interactions must be defined
  cof_cols <- intersect(c("agonist", "antagonist", "co_A_receptor", "co_I_receptor"), colnames(it))
  cof_used <- unique(unlist(it[, cof_cols]))
  cof_used <- cof_used[!is.na(cof_used) & cof_used != ""]
  put(sp, "cofactor names used but undefined", length(setdiff(cof_used, rownames(db$cofactor %||% data.frame()))), scope)
  # Nomenclature: case convention and symbols typical of another species
  genes_all <- partner_genes(partners, cx)
  bad_case <- if (sp %in% upper_species) genes_all[grepl("[a-z]", genes_all) & !grepl("^(LOC|ENS|orf)", genes_all) & !grepl("orf", genes_all)]
    else if (sp %in% rodent_species) genes_all[grepl("^[A-Z0-9-]{3,}$", genes_all) & grepl("[A-Z]{2,}", genes_all)]
    else if (sp %in% lower_species) genes_all[grepl("[A-Z]", genes_all) & !grepl("^(LOC|ENS|si:|zgc:|CABZ|BX|CR|CU|FP|AL)", genes_all)]
    else character(0)
  foreign <- c(if (!sp %in% c("mouse", "rat")) genes_all[grepl("^H2-", genes_all)],
               if (!sp %in% c("human", "chimp", "macaque", "marmoset")) genes_all[grepl("^HLA-", genes_all)],
               genes_all[grepl("^Gm[0-9]+$", genes_all)])
  put(sp, "symbols breaking the species case convention", length(bad_case), scope)
  put(sp, "symbols typical of another species (H2-, HLA-, Gm)", length(foreign), scope)
  keep("symbol_convention_issues", rbind(if (length(bad_case)) data.frame(species = sp, scope = scope, symbol = bad_case, issue = "case convention", stringsAsFactors = FALSE),
                                         if (length(foreign)) data.frame(species = sp, scope = scope, symbol = foreign, issue = "other species", stringsAsFactors = FALSE)))
  invisible(it)
}

for (sp in species) {
  sp_dir <- file.path(app_root, sp)
  merged <- readRDS(file.path(sp_dir, paste0("CellHubDB_", sp, "_merged.rds")))
  it <- analyze_db(merged, sp, "merged")
  # Release object (CellHubDB.rds) identical to the app's merged file
  rel <- release[[sp]]
  put(sp, "release CellHubDB.rds identical to app merged file",
      if (is.null(rel)) "species missing from CellHubDB.rds" else isTRUE(all.equal(rel$interaction, merged$interaction, check.attributes = FALSE)))
  # Tiers, curated sources and role corroboration
  tier <- col_or(it, "confidence_tier")
  for (t in c("high", "medium", "low")) put(sp, paste0("tier ", t), sum(tier == t, na.rm = TRUE))
  put(sp, "tier missing", sum(is.na(tier)))
  rc <- col_or(it, "role_corroborated")
  if (!all(is.na(rc))) {
    put(sp, "medium with roles corroborated", sum(tier == "medium" & rc %in% TRUE, na.rm = TRUE))
    put(sp, "medium without corroborated roles", sum(tier == "medium" & !(rc %in% TRUE), na.rm = TRUE))
  }
  # Sources: entries per source file, and merged interactions crediting each source
  src_files <- list.files(sp_dir, pattern = paste0("^CellHubDB_", sp, "_.*\\.rds$"))
  src_files <- src_files[!grepl("_merged\\.rds$", src_files)]
  sources_col <- strsplit(as.character(col_or(it, "sources", "")), ",\\s*")
  put(sp, "number of sources", length(src_files))
  put(sp, "interactions from a single source", sum(lengths(sources_col) == 1))
  for (f in src_files) {
    s <- sub("\\.rds$", "", sub(paste0("^CellHubDB_", sp, "_"), "", f))
    d <- readRDS(file.path(sp_dir, f))
    analyze_db(d, sp, paste("source", s))
    put(sp, paste0("merged interactions crediting ", s), sum(vapply(sources_col, function(x) s %in% x, logical(1))))
  }
  if (sp == "human") {
    eg <- it[grepl("(^|_)EGFR(_|$)", it$ligand) | grepl("(^|_)EGFR(_|$)", it$receptor), , drop = FALSE]
    keep("egfr_human", data.frame(ligand = eg$ligand, receptor = eg$receptor, tier = col_or(eg, "confidence_tier"), role_corroborated = col_or(eg, "role_corroborated"),
                                  curated_sources = col_or(eg, "curated_sources"), sources = col_or(eg, "sources"), stringsAsFactors = FALSE))
  }
  message("done: ", sp)
}

m <- do.call(rbind, metrics)
write.csv(m, file.path(out_dir, "deep_metrics.csv"), row.names = FALSE)
for (nm in names(details)) write.csv(details[[nm]], file.path(out_dir, paste0(nm, ".csv")), row.names = FALSE)

# Markdown summary: one table of key numbers per species (merged), then the per-source entries
mm <- m[m$scope == "merged", ]
wide_metric <- function(metric) vapply(species, function(sp) { v <- mm$value[mm$species == sp & mm$metric == metric]; if (length(v)) v[1] else "" }, character(1))
key_metrics <- c("interactions", "number of sources", "tier high", "tier medium", "tier low", "medium with roles corroborated",
                 "% scorable (all partner genes in geneInfo), protein interactions", "non-protein signaling interactions", "pathways (distinct, non-empty)",
                 "% without pathway", "type: ECM-Receptor", "type: missing", "% CellChat annotation column with a canonical category", "complexes used", "multi-subunit partners without definition", "complexes with subunits absent from geneInfo",
                 "duplicated oriented ligand->receptor pairs", "interactions in reciprocal pairs", "symbols breaking the species case convention",
                 "symbols typical of another species (H2-, HLA-, Gm)", "release CellHubDB.rds identical to app merged file")
tab <- data.frame(metric = key_metrics, do.call(cbind, lapply(species, function(sp) vapply(key_metrics, function(k) { v <- mm$value[mm$species == sp & mm$metric == k]; if (length(v)) v[1] else "" }, character(1)))), check.names = FALSE)
colnames(tab) <- c("metric", species)
md <- c("# CellHubDB 2.1 — numbers recomputed from the release", "", paste("Release folder:", release_dir), paste("Run:", format(Sys.time())), "",
        sprintf("Total merged interactions, all species: **%s**", format(sum(as.numeric(wide_metric("interactions"))), big.mark = ",")), "",
        "## Merged database per species", "", paste("|", paste(colnames(tab), collapse = " | "), "|"), paste("|", paste(rep("---", ncol(tab)), collapse = " | "), "|"),
        apply(tab, 1, function(r) paste("|", paste(r, collapse = " | "), "|")), "",
        "## Entries per source file (source-level, before cross-source deduplication)", "")
src_rows <- m[grepl("^source ", m$scope) & m$metric == "interactions", ]
md <- c(md, "| species | source | entries | merged interactions crediting it |", "| --- | --- | --- | --- |",
        apply(src_rows, 1, function(r) {
          s <- sub("^source ", "", r[["scope"]])
          cr <- mm$value[mm$species == r[["species"]] & mm$metric == paste0("merged interactions crediting ", s)]
          sprintf("| %s | %s | %s | %s |", r[["species"]], s, r[["value"]], if (length(cr)) cr else "")
        }),
        "", sprintf("Total source-level entries: **%s**", format(sum(as.numeric(src_rows$value)), big.mark = ",")), "",
        "Detail tables: missing_symbols.csv, complex_issues.csv, symbol_convention_issues.csv, reciprocal_pairs.csv, duplicates.csv, egfr_human.csv")
writeLines(md, file.path(out_dir, "deep_summary.md"))
message("Summary: ", file.path(out_dir, "deep_summary.md"))
