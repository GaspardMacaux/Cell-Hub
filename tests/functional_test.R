# functional_test.R — Cell-Hub functional test: the app's own functions are run on the demo data and compared with an
# independent reference computed here with plain Seurat / CellChat / harmony calls, plus biological sanity checks.
# Runs inside the Cell-Hub image (same package versions as the app), without the interface.
# Usage: Rscript functional_test.R <app_dir> <demo_dir> <out_dir>
# Output: <out_dir>/functional_results.csv (section, check, status, value, expected, detail) and functional_test.log
# Status: OK = identical to the reference or biologically expected; WARN = plausible but outside the expected range;
#         KO = the app differs from the reference or fails; SKIP = package or file not available.
args <- commandArgs(trailingOnly = TRUE)
app_dir <- if (length(args) >= 1) args[1] else "/srv/shiny-server/mon_app"
demo_dir <- if (length(args) >= 2) args[2] else "/demo"
out_dir <- if (length(args) >= 3) args[3] else "/out"
suppressPackageStartupMessages({ library(Seurat); library(shiny); library(Matrix); library(dplyr); library(ggplot2) })
options(warn = 1, future.globals.maxSize = 8 * 1024^3)
`%||%` <- function(a, b) if (is.null(a)) b else a
quiet <- function(...) invisible(NULL)
t_start <- Sys.time()

############################## Result recording ##############################

results <- list()
rec <- function(section, check, status, value = "", expected = "", detail = "") {
  results[[length(results) + 1]] <<- data.frame(section = section, check = check, status = status, value = paste(value, collapse = " "),
                                                expected = paste(expected, collapse = " "), detail = substr(paste(detail, collapse = " "), 1, 800),
                                                stringsAsFactors = FALSE)
  message(sprintf("[%s] %s / %s: %s %s", status, section, check, paste(value, collapse = " "), if (nzchar(paste(detail, collapse = ""))) paste0("(", substr(paste(detail, collapse = " "), 1, 200), ")") else ""))
}
save_results <- function() {
  if (length(results) > 0) utils::write.csv(do.call(rbind, results), file.path(out_dir, "functional_results.csv"), row.names = FALSE)
}
# Runs one section; an error stops that section only and is recorded as KO
run_section <- function(name, fn) {
  t0 <- Sys.time()
  message("\n==================== ", name, " ====================")
  tryCatch(fn(), error = function(e) rec(name, "section aborted", "KO", detail = conditionMessage(e)))
  rec(name, "duration (min)", "INFO", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  save_results()
}
# OK when the condition holds, otherwise the given failure status
check <- function(section, name, ok, value = "", expected = "", detail = "", fail = "KO") rec(section, name, if (isTRUE(ok)) "OK" else fail, value, expected, detail)

############################## Small statistics ##############################

# Adjusted Rand index of two partitions (1 = identical partitions up to label names)
ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); s2 <- function(x) sum(choose(x, 2))
  sa <- s2(rowSums(tab)); sb <- s2(colSums(tab)); expected <- sa * sb / choose(n, 2)
  (s2(tab) - expected) / ((sa + sb) / 2 - expected)
}
# Batch mixing: mean fraction of the k nearest neighbours of each cell that come from another batch, and the value expected
# under perfect mixing (batch proportions). Uncorrected data with a batch effect sit well below the expected value.
batch_mixing <- function(emb, batch, k = 20) {
  batch <- as.character(batch)
  nn <- RANN::nn2(emb, k = k + 1)$nn.idx[, -1, drop = FALSE]
  observed <- mean(rowMeans(matrix(batch[nn], nrow = nrow(nn)) != batch))
  sizes <- table(batch)
  expected <- mean(1 - (as.numeric(sizes[batch]) - 1) / (length(batch) - 1))
  c(observed = observed, expected = expected)
}
max_abs_diff <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)), na.rm = TRUE)
# Cluster with the highest mean log-normalized expression of a gene
top_cluster <- function(obj, gene, clusters) {
  if (!gene %in% rownames(obj)) return(NA_character_)
  x <- LayerData(obj, assay = "RNA", layer = "data")[gene, ]
  m <- tapply(as.numeric(x), as.character(clusters), mean)
  names(m)[which.max(m)]
}
read_h5_counts <- function(path) {
  x <- Read10X_h5(path)
  if (is.list(x)) x <- x[["Gene Expression"]]
  x
}

############################## App helpers ##############################

# Helper files of the app: _docker names (image) or dev names, anywhere under the app folder (server and UI files excluded)
source_helper <- function(stem) {
  files <- list.files(app_dir, pattern = paste0("^", stem, "(_docker)?\\.R$"), recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("server|_ui", basename(files))]
  if (length(files) == 0) { rec("setup", paste("source", stem), "KO", detail = "file not found"); return(invisible(FALSE)) }
  f <- files[order(!grepl("_docker", files))][1]
  ok <- tryCatch({ source(f, local = globalenv()); TRUE }, error = function(e) { rec("setup", paste("source", stem), "KO", detail = conditionMessage(e)); FALSE })
  if (ok) rec("setup", paste("source", stem), "OK", basename(f))
  invisible(ok)
}
for (stem in c("dimensional_reduction_functions", "cells_genes_expressions", "integration_functions", "data_loading_functions",
               "cellchat_helpers", "cellhubdb_filter_functions", "spatial_functions", "workspace_management")) source_helper(stem)
# The app attaches these packages at start-up (its functions call RunHarmony etc. without a namespace prefix)
for (p in c("harmony", "batchelor")) if (requireNamespace(p, quietly = TRUE)) suppressPackageStartupMessages(library(p, character.only = TRUE))
for (p in c("CellChat", "harmony", "batchelor", "presto", "RANN")) rec("setup", paste("package", p), if (requireNamespace(p, quietly = TRUE)) "INFO" else "SKIP",
                                                                  if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "not installed")
v3_file <- file.path(demo_dir, "pbmc_1k_v3_filtered_feature_bc_matrix.h5")
v2_file <- file.path(demo_dir, "pbmc_1k_v2_filtered_feature_bc_matrix.h5")
visium_zip <- file.path(demo_dir, "visium_mouse_brain.zip")
qc <- list(min_features_merge = 200, max_features_merge = 6000, max_mt_percent_merge = 20)
pbmc_genes <- c(T = "CD3E", B = "MS4A1", Mono = "CD14", NK = "NKG7")

############################## 1. Loading and QC ##############################

ref_raw <- NULL; app_v3 <- NULL; app_v2 <- NULL
run_section("1 loading + QC", function() {
  s <- "1 loading + QC"
  counts <- read_h5_counts(v3_file)
  raw <- CreateSeuratObject(counts, project = "PBMC_v3", min.cells = 3, min.features = qc$min_features_merge)
  raw[["percent.mt"]] <- PercentageFeatureSet(raw, pattern = "^MT-")
  ref_raw <<- raw
  md <- raw@meta.data
  # Reference rule written here: strict inequalities on genes and mitochondrial %, as documented in the parameter export
  keep_ref <- md$nFeature_RNA > qc$min_features_merge & md$nFeature_RNA < qc$max_features_merge & md$percent.mt < qc$max_mt_percent_merge
  keep_app <- qc_pass_cells(md, qc$min_features_merge, qc$max_features_merge, qc$max_mt_percent_merge)
  check(s, "qc_pass_cells = reference rule", identical(unname(keep_ref), unname(keep_app)), sum(keep_app), sum(keep_ref),
        sprintf("%d cells differ", sum(keep_ref != keep_app)))
  # Multiple-dataset loader (Cell Ranger .h5): cells kept and log-normalized values compared with the reference
  app_v3 <<- preprocessRawDataset(v3_file, "snRNA_merge", "human", "PBMC_v3", qc, "LogNormalize", log_function = quiet)
  app_v2 <<- preprocessRawDataset(v2_file, "snRNA_merge", "human", "PBMC_v2", qc, "LogNormalize", log_function = quiet)
  ref_cells <- colnames(ref_raw)[keep_ref]
  check(s, "loader: same cells as the reference", setequal(colnames(app_v3), ref_cells), ncol(app_v3), length(ref_cells),
        sprintf("only app %d, only reference %d", length(setdiff(colnames(app_v3), ref_cells)), length(setdiff(ref_cells, colnames(app_v3)))))
  ref_norm <- NormalizeData(subset(ref_raw, cells = ref_cells), verbose = FALSE)
  genes <- intersect(rownames(app_v3), rownames(ref_norm))
  common <- intersect(colnames(app_v3), ref_cells)
  d <- max_abs_diff(LayerData(app_v3, "RNA", layer = "data")[genes, common], LayerData(ref_norm, "RNA", layer = "data")[genes, common])
  check(s, "loader: log-normalized values = NormalizeData", d < 1e-8, signif(d, 3), "< 1e-8")
  check(s, "loader: counts unchanged", max_abs_diff(LayerData(app_v3, "RNA", layer = "counts")[genes, common], counts[genes, common]) == 0, "", "0 difference")
  check(s, "PBMC v3 cells after QC in plausible range", ncol(app_v3) >= 900 && ncol(app_v3) <= 1250, ncol(app_v3), "900-1250", fail = "WARN")
  check(s, "PBMC v2 cells after QC in plausible range", ncol(app_v2) >= 800 && ncol(app_v2) <= 1250, ncol(app_v2), "800-1250", fail = "WARN")
})

############################## 2. Single-dataset pipeline vs plain Seurat ##############################

app_obj <- NULL
run_section("2 single pipeline", function() {
  s <- "2 single pipeline"
  base <- subset(ref_raw, cells = colnames(ref_raw)[qc_pass_cells(ref_raw@meta.data, 200, 6000, 20)])
  base <- NormalizeData(base, verbose = FALSE)
  base <- FindVariableFeatures(base, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  base <- ScaleData(base, features = VariableFeatures(base), verbose = FALSE)
  # App wrappers (seed 42)
  app_path <- function(o) {
    o <- runPCA_reproducible(o, npcs = 50, seed = 42)
    o <- findNeighbors_reproducible(o, dims = 1:30, seed = 42)
    o <- findClusters_reproducible(o, resolution = 0.5, algorithm = 1, seed = 42)
    runUMAP_reproducible(o, dims = 1:30, seed = 42)
  }
  app_obj <<- app_path(base)
  # Plain Seurat with the same parameters
  ref <- RunPCA(base, npcs = 50, seed.use = 42, verbose = FALSE)
  ref <- FindNeighbors(ref, dims = 1:30, k.param = 20, verbose = FALSE)
  ref <- FindClusters(ref, resolution = 0.5, algorithm = 1, random.seed = 42, verbose = FALSE)
  ref <- RunUMAP(ref, dims = 1:30, seed.use = 42, n.neighbors = 30L, min.dist = 0.3, metric = "cosine", verbose = FALSE)
  d_pca <- max_abs_diff(Embeddings(app_obj, "pca"), Embeddings(ref, "pca"))
  check(s, "PCA = Seurat RunPCA", d_pca < 1e-8, signif(d_pca, 3), "< 1e-8")
  a <- ari(app_obj$seurat_clusters, ref$seurat_clusters)
  check(s, "clusters = Seurat FindClusters (ARI)", a == 1, round(a, 4), "1")
  d_umap <- max_abs_diff(Embeddings(app_obj, "umap"), Embeddings(ref, "umap"))
  check(s, "UMAP = Seurat RunUMAP", d_umap < 1e-6, signif(d_umap, 3), "< 1e-6")
  # Reproducibility: the app path run a second time gives identical results
  again <- app_path(base)
  check(s, "reproducible: clusters identical on re-run", identical(as.character(again$seurat_clusters), as.character(app_obj$seurat_clusters)))
  check(s, "reproducible: UMAP identical on re-run", identical(Embeddings(again, "umap"), Embeddings(app_obj, "umap")))
  n_cl <- length(unique(app_obj$seurat_clusters))
  check(s, "number of clusters plausible for PBMC (res 0.5)", n_cl >= 6 && n_cl <= 16, n_cl, "6-16", fail = "WARN")
  # Biology: T, B, monocyte and NK markers peak in different clusters
  tops <- vapply(pbmc_genes, function(g) top_cluster(app_obj, g, app_obj$seurat_clusters), character(1))
  check(s, "markers peak in distinct clusters (T/B/Mono/NK)", length(unique(tops)) == length(tops), paste(names(tops), tops, sep = ":", collapse = " "),
        "4 distinct clusters", fail = "WARN")
})

############################## 3. Differential expression vs Seurat ##############################

run_section("3 exclusive markers", function() {
  s <- "3 exclusive markers"
  obj <- app_obj
  Idents(obj) <- "seurat_clusters"
  target <- top_cluster(obj, "CD14", obj$seurat_clusters)
  ex <- find_exclusive_biomarkers(obj, target_cluster = target, verbose = FALSE)
  check(s, "markers found for the CD14 cluster", nrow(ex) > 0, nrow(ex), "> 0")
  if (nrow(ex) == 0) return(invisible())
  fm <- FindMarkers(obj, ident.1 = target, features = ex$gene, logfc.threshold = 0, min.pct = 0, verbose = FALSE)
  fm <- fm[ex$gene, ]
  d_all <- abs(ex$log2_fold_change - fm$avg_log2FC); d_fc <- max(d_all, na.rm = TRUE); w <- which.max(d_all)
  # Exact agreement expected; a difference below 0.01 log2 unit is reported as WARN (numerical, no effect on any filter)
  # find_exclusive_biomarkers rounds log2FC to 3 decimals: differences up to 0.0005 are that rounding
  rec(s, "log2FC = Seurat avg_log2FC (3-decimal rounding)", if (d_fc <= 5e-4 + 1e-9) "OK" else if (d_fc < 0.01) "WARN" else "KO", signif(d_fc, 3), "<= 0.0005",
      sprintf("largest difference: %s app %.6f vs Seurat %.6f | median difference %.2g | genes differing > 1e-6: %d of %d", ex$gene[w], ex$log2_fold_change[w],
              fm$avg_log2FC[w], stats::median(d_all, na.rm = TRUE), sum(d_all > 1e-6, na.rm = TRUE), length(d_all)))
  # Seurat rounds pct to 3 decimals
  d_pct <- max_abs_diff(ex$pct_detected_target, 100 * fm$pct.1)
  check(s, "% detected target = Seurat pct.1", d_pct <= 0.1, signif(d_pct, 3), "<= 0.1")
  if ("pvalue" %in% colnames(ex)) {
    ok <- ex$pvalue > 0 & fm$p_val > 0
    d_p <- if (any(ok)) max_abs_diff(log10(ex$pvalue[ok]), log10(fm$p_val[ok])) else 0
    check(s, "p-values = Seurat Wilcoxon", d_p < 1e-6, signif(d_p, 3), "< 1e-6 (log10)")
    check(s, "p-values stored as numbers (Venn rule)", is.numeric(ex$pvalue) && is.numeric(ex$pvalue_adjusted))
  }
  # Every returned gene satisfies the default filters
  ok_filters <- all(ex$pct_detected_target >= 50 & ex$pct_detected_other <= 25 & ex$log2_fold_change >= 1.5)
  check(s, "returned genes satisfy the filters (50% / 25% / log2FC 1.5)", ok_filters)
  check(s, "monocyte genes among the markers", any(c("CD14", "LYZ", "S100A8", "S100A9", "FCN1") %in% ex$gene),
        paste(head(ex$gene, 10), collapse = ", "), "CD14 / LYZ / S100A8 / S100A9 / FCN1", fail = "WARN")
})

############################## 4. CellChat with CellHubDB vs plain CellChat ##############################

run_section("4 CellChat", function() {
  s <- "4 CellChat"
  if (!requireNamespace("CellChat", quietly = TRUE)) { rec(s, "CellChat", "SKIP", detail = "not installed"); return(invisible()) }
  suppressPackageStartupMessages(library(CellChat))
  merged <- readRDS(file.path(app_dir, "databases", "human", "CellHubDB_human_merged.rds"))
  db_version <- attr(merged, "version") %||% "0"
  check(s, "database version", isTRUE(utils::compareVersion(db_version, "2.1.2") >= 0), db_version, ">= 2.1.2 (canonical annotation + complete geneInfo)")
  # CellChat drops single-gene partners absent from geneInfo (extractGeneSubset): every partner gene must be listed
  it0 <- merged$interaction; partners <- unique(c(it0$ligand, it0$receptor)); single <- setdiff(partners, rownames(merged$complex))
  check(s, "every single-gene partner in geneInfo", all(single %in% merged$geneInfo$Symbol), sum(!single %in% merged$geneInfo$Symbol), "0 missing")
  merged$interaction <- annotate_cellhubdb_tier(merged$interaction)
  f <- filter_cellhubdb_tiers(merged, c("high", "medium"), medium_roles_only = TRUE)
  it <- merged$interaction
  ref_keep <- it$interaction_name[(it$confidence_tier %in% "high") | (it$confidence_tier %in% "medium" & as.logical(it$role_corroborated) %in% TRUE)]
  check(s, "tier filter = reference subset (high + medium confirmed roles)", setequal(f$db$interaction$interaction_name, ref_keep), f$n_after, length(ref_keep))
  db <- prepare_cellchat_db(f$db)$db
  check(s, "annotation column canonical (subsetDB-compatible)", mean(db$interaction$annotation %in% c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor", "Non-protein Signaling")) > 0.5,
        round(100 * mean(db$interaction$annotation %in% c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor", "Non-protein Signaling")), 1), "> 50 %")
  # App path: numeric seurat_clusters as group labels (the case that failed with "Cell labels cannot contain 0")
  cc <- create_cellchat_object_multicol(app_obj, "seurat_clusters", log_function = quiet)
  check(s, "CellChat object from numeric clusters", isTRUE(cc$success), "", "created", cc$message %||% "")
  if (!isTRUE(cc$success)) return(invisible())
  cc$object@DB <- db
  r <- run_cellchat_pipeline(cc$object, "app", log_function = quiet)
  check(s, "app pipeline completes", isTRUE(r$success), r$n_interactions, "> 0 interactions", r$message)
  if (!isTRUE(r$success)) return(invisible())
  # Plain CellChat on the same matrix, labels (same level order) and database
  run_plain <- function(database) {
    lab <- paste0("C", as.character(app_obj$seurat_clusters))
    meta <- data.frame(group = factor(lab, levels = levels(r$object@idents)), row.names = colnames(app_obj))
    m <- createCellChat(object = LayerData(app_obj, assay = "RNA", layer = "data"), meta = meta, group.by = "group")
    m@DB <- database
    m <- subsetData(m); m <- identifyOverExpressedGenes(m); m <- identifyOverExpressedInteractions(m)
    m <- computeCommunProb(m, seed.use = 42); m <- filterCommunication(m, min.cells = 10)
    m <- computeCommunProbPathway(m); aggregateNet(m)
  }
  ref <- run_plain(db)
  p_app <- r$object@net$prob; p_ref <- ref@net$prob
  same_dims <- identical(dimnames(p_app), dimnames(p_ref))
  check(s, "same groups and LR pairs as plain CellChat", same_dims, paste(dim(p_app), collapse = "x"), paste(dim(p_ref), collapse = "x"))
  if (same_dims) {
    check(s, "communication probabilities = plain CellChat", max_abs_diff(p_app, p_ref) < 1e-12, signif(max_abs_diff(p_app, p_ref), 3), "< 1e-12")
    check(s, "significant interactions = plain CellChat", identical(unname(as.matrix(r$object@net$count)), unname(as.matrix(ref@net$count))), sum(r$object@net$count), sum(ref@net$count))
  }
  pw <- r$object@netP$pathways
  expected_pw <- c("MIF", "MHC-I", "MHC-II", "GALECTIN", "CD99", "APP", "ANNEXIN", "CLEC", "CD45", "ITGB2", "LCK", "BAFF", "CD22", "CCL")
  check(s, "PBMC pathways found", sum(expected_pw %in% pw) >= 3, paste(intersect(expected_pw, pw), collapse = ", "), ">= 3 of the usual PBMC pathways",
        paste("top:", paste(head(pw, 15), collapse = ", ")), fail = "WARN")
  # Agreement with the official CellChatDB.human: share of its significant ligand-receptor pairs also found with CellHubDB
  off <- run_plain(CellChat::CellChatDB.human)
  key <- function(o) { d <- subsetCommunication(o); unique(toupper(paste(d$ligand, d$receptor))) }
  k_off <- key(off); k_hub <- key(r$object)
  rec_frac <- if (length(k_off) > 0) mean(k_off %in% k_hub) else NA
  check(s, "official CellChatDB pairs recovered by CellHubDB", isTRUE(rec_frac >= 0.6), sprintf("%.0f %% of %d", 100 * rec_frac, length(k_off)), ">= 60 %",
        paste("missing (first 15):", paste(head(setdiff(k_off, k_hub), 15), collapse = ", ")), fail = "WARN")
})

############################## 5. Integration vs plain harmony, batch mixing ##############################

run_section("5 integration", function() {
  s <- "5 integration"
  lst <- list(app_v3, app_v2)
  m0 <- performSimpleMerge(lst)
  check(s, "simple merge keeps every cell", ncol(m0) == ncol(app_v3) + ncol(app_v2), ncol(m0), ncol(app_v3) + ncol(app_v2))
  check(s, "simple merge keeps both datasets", setequal(unique(m0$dataset), c("PBMC_v3", "PBMC_v2")), paste(unique(m0$dataset), collapse = ", "), "PBMC_v3, PBMC_v2")
  pca_mix <- batch_mixing(Embeddings(m0, "pca")[, 1:30], m0$dataset)
  rec(s, "batch mixing without correction (PCA)", "INFO", round(pca_mix["observed"], 3), paste("perfect mixing =", round(pca_mix["expected"], 3)))
  if (requireNamespace("harmony", quietly = TRUE)) {
    h <- performHarmonyIntegration(lst)
    # Plain harmony on the same merge, same steps and seed
    ref <- performSimpleMerge(lst)
    ref <- NormalizeData(ref, verbose = FALSE)
    ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
    ref <- ScaleData(ref, features = VariableFeatures(ref), verbose = FALSE)
    ref <- RunPCA(ref, npcs = 50, seed.use = 42, features = VariableFeatures(ref), verbose = FALSE)
    set.seed(42)
    ref <- if (utils::packageVersion("harmony") >= "1.0.0") {
      harmony::RunHarmony(ref, group.by.vars = "dataset", dims.use = 1:30, max_iter = 20, verbose = FALSE, plot_convergence = FALSE)
    } else {
      harmony::RunHarmony(ref, group.by.vars = "dataset", dims.use = 1:30, max.iter.harmony = 20, verbose = FALSE, plot_convergence = FALSE)
    }
    d_h <- max_abs_diff(Embeddings(h, "harmony"), Embeddings(ref, "harmony")[colnames(h), ])
    check(s, "Harmony embedding = plain RunHarmony", d_h < 1e-6, signif(d_h, 3), "< 1e-6")
    h_mix <- batch_mixing(Embeddings(h, "harmony")[, 1:min(30, ncol(Embeddings(h, "harmony")))], h$dataset)
    check(s, "Harmony improves batch mixing", h_mix["observed"] > pca_mix["observed"] + 0.05, round(h_mix["observed"], 3),
          sprintf("> %.3f (uncorrected + 0.05)", pca_mix["observed"]), fail = "WARN")
    # Biology kept after correction: T, B and monocyte markers still peak in distinct clusters
    hc <- findNeighbors_reproducible(h, dims = 1:30, reduction = "harmony", seed = 42)
    hc <- findClusters_reproducible(hc, resolution = 0.5, seed = 42)
    tops <- vapply(pbmc_genes[1:3], function(g) top_cluster(hc, g, hc$seurat_clusters), character(1))
    check(s, "Harmony keeps T / B / Mono separated", length(unique(tops)) == 3, paste(names(tops), tops, sep = ":", collapse = " "), "3 distinct clusters", fail = "WARN")
  } else rec(s, "Harmony", "SKIP", detail = "harmony not installed")
  if (requireNamespace("batchelor", quietly = TRUE)) {
    fm <- performFastMNNIntegration(lst)
    mnn_mix <- batch_mixing(Embeddings(fm, "integrated.mnn")[, 1:30], fm$dataset)
    check(s, "fastMNN improves batch mixing", mnn_mix["observed"] > pca_mix["observed"] + 0.05, round(mnn_mix["observed"], 3),
          sprintf("> %.3f (uncorrected + 0.05)", pca_mix["observed"]), fail = "WARN")
  } else rec(s, "fastMNN", "SKIP", detail = "batchelor not installed")
  # RPCA: same embedding as a plain Seurat v5 IntegrateLayers run on the same split object, and better batch mixing
  rp <- performRPCAIntegration(lst)
  ref_rp <- performSimpleMerge(lst)
  ref_rp[["RNA"]] <- split(ref_rp[["RNA"]], f = as.character(ref_rp$dataset))
  ref_rp <- NormalizeData(ref_rp, verbose = FALSE)
  ref_rp <- FindVariableFeatures(ref_rp, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  ref_rp <- ScaleData(ref_rp, features = VariableFeatures(ref_rp), verbose = FALSE)
  ref_rp <- RunPCA(ref_rp, npcs = 50, seed.use = 42, features = VariableFeatures(ref_rp), verbose = FALSE)
  set.seed(42)
  ref_rp <- IntegrateLayers(ref_rp, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", dims = 1:30, k.anchor = 5,
                            k.weight = rp@misc$integration_params$k_weight, verbose = FALSE)
  d_rp <- max_abs_diff(Embeddings(rp, "integrated.rpca"), Embeddings(ref_rp, "integrated.rpca")[colnames(rp), ])
  check(s, "RPCA embedding = plain IntegrateLayers(RPCAIntegration)", d_rp < 1e-6, signif(d_rp, 3), "< 1e-6")
  check(s, "RPCA layers joined back", length(Layers(rp[["RNA"]], search = "counts")) == 1, paste(Layers(rp[["RNA"]]), collapse = ", "), "one counts layer")
  rp_mix <- batch_mixing(Embeddings(rp, "integrated.rpca")[, 1:30], rp$dataset)
  check(s, "RPCA improves batch mixing", rp_mix["observed"] > pca_mix["observed"] + 0.05, round(rp_mix["observed"], 3),
        sprintf("> %.3f (uncorrected + 0.05)", pca_mix["observed"]), fail = "WARN")
  cca <- performStandardIntegration(lst)
  DefaultAssay(cca) <- "integrated"
  cca <- ScaleData(cca, verbose = FALSE)
  cca <- runPCA_reproducible(cca, npcs = 30, seed = 42)
  cca_mix <- batch_mixing(Embeddings(cca, "pca")[, 1:30], cca$dataset)
  check(s, "Seurat CCA improves batch mixing", cca_mix["observed"] > pca_mix["observed"] + 0.05, round(cca_mix["observed"], 3),
        sprintf("> %.3f (uncorrected + 0.05)", pca_mix["observed"]), fail = "WARN")
})

############################## 6. Spatial: loading and co-localization vs direct computation ##############################

run_section("6 spatial", function() {
  s <- "6 spatial"
  tmp <- tempfile(pattern = "visium_")
  staged <- prepare_visium_directory(visium_zip, tmp, log_function = quiet)
  sp <- load_visium_object(staged, species = "mouse")
  sp@misc$spot_size_um <- staged$spot_size_um
  check(s, "spot size of standard Visium", isTRUE(staged$spot_size_um == 55), staged$spot_size_um, "55")
  h5 <- list.files(tmp, pattern = "filtered_feature_bc_matrix\\.h5$", recursive = TRUE, full.names = TRUE)[1]
  ref_counts <- read_h5_counts(h5)
  common <- intersect(colnames(sp), colnames(ref_counts))
  check(s, "every loaded spot is in the Space Ranger matrix", length(common) == ncol(sp), ncol(sp), paste(ncol(ref_counts), "in the matrix"))
  genes <- intersect(rownames(sp), rownames(ref_counts))
  check(s, "loaded counts = Space Ranger matrix", max_abs_diff(LayerData(sp, "Spatial", layer = "counts")[genes, common], ref_counts[genes, common]) == 0, "", "0 difference")
  pos <- get_spot_positions_um(sp, 55)
  # Visium spots are printed on a 100 um centre-to-centre grid: the micrometre conversion must give that pitch
  check(s, "grid pitch ~100 um (calibration)", pos$pitch_um >= 95 && pos$pitch_um <= 105, round(pos$pitch_um, 1), "95-105",
        sprintf("um per pixel %.4f, spot_diameter_fullres %.2f px", pos$um_per_px, sp@images[[1]]@scale.factors$spot))
  sp <- NormalizeData(sp, assay = "Spatial", verbose = FALSE)
  g <- c("Mbp", "Plp1")
  if (!all(g %in% rownames(sp))) { rec(s, "co-localization", "SKIP", detail = "Mbp / Plp1 absent"); return(invisible()) }
  res <- analyze_gene_colocalization(sp, g, max_distance_um = 100, spot_size_um = 55, assay = "Spatial", bidirectional = TRUE)
  # Direct computation: positive spots, then every pair distance (no neighbour-search library)
  xy <- pos$xy; rownames(xy) <- pos$barcodes
  cts <- LayerData(sp, "Spatial", layer = "counts")[g, pos$barcodes]
  p1 <- which(cts[1, ] > 0); p2 <- which(cts[2, ] > 0)
  check(s, "positive spots counted like the reference (counts > 0)", res$gene1_total_spots == length(p1) && res$gene2_total_spots == length(p2),
        paste(res$gene1_total_spots, res$gene2_total_spots), paste(length(p1), length(p2)))
  near <- function(a, b) {
    d2 <- outer(xy[a, 1], xy[b, 1], "-")^2 + outer(xy[a, 2], xy[b, 2], "-")^2
    sum(rowSums(d2 <= 100^2 * (1 + 1e-9)) > 0)
  }
  n12 <- near(p1, p2); n21 <- near(p2, p1)
  check(s, "Mbp spots near Plp1 = direct computation", res$gene1_near_gene2 == n12, res$gene1_near_gene2, n12)
  check(s, "Plp1 spots near Mbp = direct computation", isTRUE(res$gene2_near_gene1 == n21), res$gene2_near_gene1, n21)
  check(s, "Mbp / Plp1 strongly co-localized (white matter)", isTRUE(res$gene1_coloc_rate >= 80), res$gene1_coloc_rate, ">= 80 %", fail = "WARN")
  unlink(tmp, recursive = TRUE)
})

rec("total", "duration (min)", "INFO", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1))
save_results()
res_df <- do.call(rbind, results)
message(sprintf("\nFunctional test: %d checks | OK %d | WARN %d | KO %d | SKIP %d", sum(res_df$status != "INFO"), sum(res_df$status == "OK"),
                sum(res_df$status == "WARN"), sum(res_df$status == "KO"), sum(res_df$status == "SKIP")))
