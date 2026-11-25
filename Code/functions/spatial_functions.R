# ==============================================================================
# spatial_functions.R
# Pure helper functions for the spatial transcriptomics module (10x Visium / Visium HD).
# NO Shiny reactive code here: no reactiveVal(), observe(), output$, input$, session$.
# All reactive logic lives in spatial_transcriptomic_server.R.
# No stochastic operation is performed here (seed-42 wrappers live in dimensional_reduction_functions.R).
# Optional dependencies checked at run time: nanoparquet (Visium HD parquet positions), RANN (co-localization),
# vips / vipsheader command-line tools (fast downsampling of uploaded custom images, magick fallback otherwise).
# ==============================================================================

############################## Data loading (Space Ranger ZIP) ##############################

# Mitochondrial gene regex for a species (mouse "mt-", human "MT-", rat "Mt-"; unknown species match all three)
get_mt_pattern_spatial <- function(species = "mouse") {
  pattern <- switch(species %||% "mouse", "mouse" = "^mt-|^Mt-", "human" = "^MT-", "rat" = "^Mt-", "^mt-|^Mt-|^MT-")
  return(pattern)
}

# Unzip a Space Ranger output and stage it into a clean folder with the layout Seurat expects:
#   <out_dir>/staged/filtered_feature_bc_matrix.h5
#   <out_dir>/staged/spatial/{scalefactors_json.json, tissue_positions*.csv, tissue_*_image.png}
# Extraction (<out_dir>/raw) and staging are separate, so Load10X_Spatial never sees a stray binned_outputs/ tree.
# Visium HD: every file comes from the square_XXXum folder matching bin_size_um, so matrix, positions and scale factors agree.
# Parquet positions are converted to a headered CSV with nanoparquet, so loading never depends on arrow.
# Returns a list describing the staged data, including spot_size_um (55 um Visium spot, bin size for Visium HD).
prepare_visium_directory <- function(zip_path, out_dir, bin_size_um = 8, log_function = message) {
  raw_dir <- file.path(out_dir, "raw"); staged_dir <- file.path(out_dir, "staged"); spatial_dir <- file.path(staged_dir, "spatial")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE); dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)
  # The system unzip handles zip64 members >= 4 GB that R's internal unzip can truncate (large Visium HD matrices)
  unzip_bin <- Sys.which("unzip")
  if (nzchar(unzip_bin)) utils::unzip(zip_path, exdir = raw_dir, unzip = unzip_bin) else utils::unzip(zip_path, exdir = raw_dir)
  all_files <- list.files(raw_dir, full.names = TRUE, recursive = TRUE)
  # Ignore macOS resource-fork junk and molecule-level files that would otherwise match the patterns
  all_files <- all_files[!grepl("__MACOSX|/\\._|molecule_info\\.h5$", all_files)]
  if (length(all_files) == 0) stop("ZIP archive is empty")
  log_function(paste("Extracted", length(all_files), "files"))
  # Visium HD binned outputs: restrict the search to a single bin folder
  hd_dirs <- unique(sub("^(.*/square_[0-9]+um)/.*$", "\\1", grep("/square_[0-9]+um/", all_files, value = TRUE)))
  if (length(hd_dirs) > 0) {
    bin_label <- sprintf("square_%03dum", as.integer(bin_size_um))
    chosen <- hd_dirs[basename(hd_dirs) == bin_label]
    if (length(chosen) == 0) stop(sprintf("Visium HD bin %s not found in the ZIP. Available: %s", bin_label, paste(basename(hd_dirs), collapse = ", ")))
    log_function(paste("Visium HD binned output detected - using", bin_label))
    all_files <- all_files[startsWith(all_files, paste0(chosen[1], "/"))]
  }
  # First file matching the patterns, tried in order of preference
  pick <- function(patterns) {
    for (p in patterns) {
      hit <- grep(p, all_files, value = TRUE, ignore.case = TRUE)
      if (length(hit)) return(hit[1])
    }
    return(NULL)
  }
  # Prefer the filtered matrix explicitly so a raw_feature_bc_matrix.h5 in the same ZIP is never picked
  h5_src <- pick(c("filtered_feature_bc_matrix\\.h5$", "(matrix|expression).*\\.h5$", "\\.h5$"))
  json_src <- pick(c("scalefactors_json\\.json$", "scalefactors.*\\.json$"))
  parquet_src <- pick(c("tissue_positions.*\\.parquet$", "positions.*\\.parquet$"))
  csv_src <- pick(c("tissue_positions.*\\.csv$", "positions.*\\.csv$"))
  hires_src <- pick("tissue_hires_image.*\\.png$")
  lowres_src <- pick("tissue_lowres_image.*\\.png$")
  log_function(sprintf("Detected - H5: %s | JSON: %s | Parquet: %s | CSV: %s | hires: %s | lowres: %s",
                       !is.null(h5_src), !is.null(json_src), !is.null(parquet_src), !is.null(csv_src), !is.null(hires_src), !is.null(lowres_src)))
  if (is.null(h5_src)) stop("No .h5 expression matrix found in the ZIP")
  if (is.null(json_src)) stop("scalefactors_json.json is missing - spots cannot be aligned to the image without it")
  if (is.null(parquet_src) && is.null(csv_src)) stop("No tissue_positions file (.parquet or .csv) found in the ZIP")
  if (is.null(hires_src) && is.null(lowres_src)) stop("No tissue image found (tissue_lowres_image.png or tissue_hires_image.png)")
  # Read the scale factors before staging: Visium HD files carry bin_size_um, standard Visium spots are 55 um in diameter
  scalefactors <- tryCatch(jsonlite::fromJSON(json_src), error = function(e) list())
  is_hd <- !is.null(scalefactors$bin_size_um)
  spot_size_um <- if (is_hd) as.numeric(scalefactors$bin_size_um) else 55
  # A ZIP without square_XXXum folders holds a single Visium HD bin: refuse a different selection before the long
  # parquet conversion and matrix loading, instead of silently loading another resolution
  if (is_hd && abs(spot_size_um - as.numeric(bin_size_um)) > 1e-6) {
    stop(sprintf("This ZIP contains Visium HD %g um bins, but %g um was selected. Select %g um, or zip the binned_outputs/square_%03dum folder (or the whole binned_outputs folder).",
                 spot_size_um, as.numeric(bin_size_um), spot_size_um, as.integer(bin_size_um)))
  }
  # Move rather than copy: extraction and staging share the temp filesystem, so no large file is duplicated on disk
  stage <- function(src, dest) {
    if (is.null(src)) return(invisible(NULL))
    if (!(suppressWarnings(file.rename(src, dest)) || file.copy(src, dest, overwrite = TRUE))) stop("Could not stage ", basename(src))
    return(dest)
  }
  h5_path <- stage(h5_src, file.path(staged_dir, "filtered_feature_bc_matrix.h5"))
  stage(json_src, file.path(spatial_dir, "scalefactors_json.json"))
  if (!is.null(parquet_src)) {
    # Convert parquet to the Space Ranger >= 2.0 headered CSV layout (barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres)
    if (!requireNamespace("nanoparquet", quietly = TRUE)) stop("Package 'nanoparquet' is required to read tissue_positions.parquet - run install.packages('nanoparquet')")
    log_function("Converting tissue_positions.parquet to CSV with nanoparquet")
    positions <- as.data.frame(nanoparquet::read_parquet(parquet_src))
    data.table::fwrite(positions, file.path(spatial_dir, "tissue_positions.csv"))
    log_function(paste("Positions converted:", nrow(positions), "rows"))
  } else {
    # Space Ranger >= 2.0 writes tissue_positions.csv WITH a header; < 2.0 writes tissue_positions_list.csv WITHOUT one.
    # Seurat picks the parser from the file name, so the name must match the actual header state.
    first_line <- readLines(csv_src, n = 1, warn = FALSE)
    csv_name <- if (grepl("barcode", first_line, ignore.case = TRUE)) "tissue_positions.csv" else "tissue_positions_list.csv"
    stage(csv_src, file.path(spatial_dir, csv_name))
  }
  stage(hires_src, file.path(spatial_dir, "tissue_hires_image.png"))
  stage(lowres_src, file.path(spatial_dir, "tissue_lowres_image.png"))
  log_function(paste("Physical spot / bin size:", spot_size_um, "um"))
  return(list(data_dir = staged_dir, spatial_dir = spatial_dir, h5_path = h5_path, positions_format = if (!is.null(parquet_src)) "parquet->csv" else "csv",
              has_hires = !is.null(hires_src), has_lowres = !is.null(lowres_src), is_hd = is_hd, spot_size_um = spot_size_um))
}

# Read a Visium / Visium HD tissue image, preferring lowres and falling back to hires.
# When only the hires image exists, the lowres scale factor is redirected to the hires one,
# because Seurat's spatial plots map fullres coordinates using scale.factors$lowres by default.
read_visium_image <- function(spatial_dir, assay = "Spatial", slice = "slice1", filter.matrix = TRUE, log_function = message) {
  if (file.exists(file.path(spatial_dir, "tissue_lowres_image.png"))) {
    log_function("Reading lowres tissue image")
    img <- Seurat::Read10X_Image(image.dir = spatial_dir, image.name = "tissue_lowres_image.png", assay = assay, slice = slice, filter.matrix = filter.matrix)
    return(img)
  }
  if (!file.exists(file.path(spatial_dir, "tissue_hires_image.png"))) stop("No tissue image found in ", spatial_dir)
  log_function("Lowres image missing - reading hires tissue image instead")
  img <- Seurat::Read10X_Image(image.dir = spatial_dir, image.name = "tissue_hires_image.png", assay = assay, slice = slice, filter.matrix = filter.matrix)
  if (is.null(img@scale.factors$hires)) stop("scalefactors_json.json lacks tissue_hires_scalef - cannot align hires image")
  # Redirect the default lowres scale factor to the hires one so spots align with the loaded image
  img@scale.factors$lowres <- img@scale.factors$hires
  return(img)
}

# Keep every Space Ranger tissue image of a staged spatial folder as compressed PNG bytes with its scale factor,
# so the displayed image can be switched later (set_spatial_image_resolution) without keeping two decoded arrays in memory.
# Returns a list named by resolution ("lowres", "hires"), each element list(png = raw, scalef = numeric).
collect_tissue_images <- function(spatial_dir) {
  scalefactors <- tryCatch(jsonlite::fromJSON(file.path(spatial_dir, "scalefactors_json.json")), error = function(e) list())
  image_files <- c(lowres = "tissue_lowres_image.png", hires = "tissue_hires_image.png")
  image_scalefs <- c(lowres = scalefactors$tissue_lowres_scalef %||% NA_real_, hires = scalefactors$tissue_hires_scalef %||% NA_real_)
  images <- list()
  for (res in names(image_files)) {
    path <- file.path(spatial_dir, image_files[[res]])
    if (file.exists(path) && is.finite(image_scalefs[[res]])) images[[res]] <- list(png = readBin(path, what = "raw", n = file.info(path)$size), scalef = as.numeric(image_scalefs[[res]]))
  }
  return(images)
}

# Build a Seurat spatial object from a staged Visium directory (see prepare_visium_directory).
# Approach 1: Load10X_Spatial with a pre-built image. Approach 2: manual H5 + image attach (barcode suffix mismatches).
# Adds in_tissue (all retained spots are in tissue) and percent.mt metadata.
load_visium_object <- function(staged, species = "mouse", assay = "Spatial", slice = "slice1", log_function = message) {
  load_errors <- character(0)
  sp_img <- read_visium_image(staged$spatial_dir, assay = assay, slice = slice, log_function = log_function)
  # Approach 1: standard Seurat loader with the image read above
  spatial_data <- tryCatch({
    log_function("Attempting Load10X_Spatial with pre-built image...")
    Load10X_Spatial(data.dir = staged$data_dir, filename = basename(staged$h5_path), assay = assay, slice = slice, filter.matrix = TRUE, to.upper = FALSE, image = sp_img)
  }, error = function(e) {
    load_errors <<- c(load_errors, paste("Load10X_Spatial:", conditionMessage(e)))
    log_function(paste("Standard loading failed:", conditionMessage(e)))
    NULL
  })
  # Approach 2: manual H5 loading, then attach the image on the shared barcodes
  if (is.null(spatial_data)) {
    spatial_data <- tryCatch({
      log_function("Attempting manual H5 loading...")
      h5_data <- Read10X_h5(staged$h5_path)
      # Multi-modal H5 files return a list - keep the gene expression matrix
      if (is.list(h5_data)) h5_data <- h5_data[["Gene Expression"]] %||% h5_data[[1]]
      obj <- CreateSeuratObject(counts = h5_data, assay = assay, project = "SpatialData")
      img_bc <- Cells(sp_img)
      common <- intersect(img_bc, colnames(obj))
      if (length(common) == 0) {
        # Align barcodes that differ only by the "-1" suffix
        obj_stripped <- sub("-1$", "", colnames(obj))
        idx <- match(obj_stripped, sub("-1$", "", img_bc))
        if (all(is.na(idx))) stop("No shared barcodes between tissue positions and expression matrix")
        new_names <- colnames(obj)
        new_names[!is.na(idx)] <- img_bc[idx[!is.na(idx)]]
        obj <- RenameCells(obj, new.names = new_names)
        common <- intersect(img_bc, colnames(obj))
      }
      log_function(paste("Shared barcodes:", length(common)))
      obj <- subset(obj, cells = common)
      img_sub <- sp_img[common]
      DefaultAssay(img_sub) <- assay
      obj[[slice]] <- img_sub
      DefaultAssay(obj) <- assay
      obj
    }, error = function(e) {
      load_errors <<- c(load_errors, paste("Manual:", conditionMessage(e)))
      log_function(paste("Manual loading failed:", conditionMessage(e)))
      NULL
    })
  }
  if (is.null(spatial_data)) stop("All loading methods failed: ", paste(load_errors, collapse = " | "))
  if (ncol(spatial_data) == 0) stop("No spots found in the loaded data")
  # Image was read with filter.matrix = TRUE, so every retained spot is in tissue (required by the QC tissue filter)
  spatial_data$in_tissue <- 1L
  spatial_data[["percent.mt"]] <- PercentageFeatureSet(spatial_data, pattern = get_mt_pattern_spatial(species))
  # Store both Space Ranger images (compressed) for switching; the decoded one is the image read above (lowres when present)
  spatial_data@misc$tissue_images <- collect_tissue_images(staged$spatial_dir)
  spatial_data@misc$active_tissue_image <- if (file.exists(file.path(staged$spatial_dir, "tissue_lowres_image.png"))) "lowres" else "hires"
  log_function(paste("Stored tissue images:", paste(names(spatial_data@misc$tissue_images), collapse = ", "), "| displayed:", spatial_data@misc$active_tissue_image))
  log_function(paste("Spatial data loaded:", ncol(spatial_data), "spots,", nrow(spatial_data), "genes"))
  return(spatial_data)
}

############################## Data loading (saved RDS object) ##############################

# Spatial-specific Seurat loader.
# ALWAYS calls UpdateSeuratObject() (spatial images need it even with Assay5),
# converts v3/v4 Assay -> Assay5 (never SCTAssay) and joins genuinely split layers.
# Returns: Seurat object ready for spatial analysis.
loadSpatialSeuratObject <- function(rds_path, species = "mouse", log_function = message) {
  log_function("Loading spatial Seurat object from RDS...")
  loaded_data <- readRDS(rds_path)
  # Handle combined Seurat + CellChat list files
  if (is.list(loaded_data) && !is.null(loaded_data$seurat)) {
    log_function("Detected combined Seurat + CellChat file, extracting Seurat object")
    loaded_obj <- loaded_data$seurat
  } else {
    loaded_obj <- loaded_data
  }
  rm(loaded_data); gc(verbose = FALSE)
  if (!inherits(loaded_obj, "Seurat")) stop("The loaded file is not a valid Seurat object")
  # VisiumV1/V2 image structures need UpdateSeuratObject even when assays are already Assay5
  # (loadSeuratObject() skips it for Assay5 objects, which breaks SpatialDimPlot: "data aligns to image slice1")
  tryCatch({
    loaded_obj <- UpdateSeuratObject(loaded_obj)
    log_function("UpdateSeuratObject completed")
  }, error = function(e) log_function(paste("UpdateSeuratObject warning:", conditionMessage(e)[1])))
  for (assay_name in names(loaded_obj@assays)) {
    assay_obj <- loaded_obj[[assay_name]]
    # SCTAssay extends Assay: converting it to Assay5 would drop the SCT models, so it is left untouched
    if (inherits(assay_obj, "Assay") && !inherits(assay_obj, "SCTAssay") && !inherits(assay_obj, "Assay5")) {
      tryCatch({
        log_function(paste("Converting", assay_name, "from Assay to Assay5"))
        loaded_obj[[assay_name]] <- as(assay_obj, "Assay5")
      }, error = function(e) log_function(paste("Failed to convert", assay_name, ":", conditionMessage(e)[1])))
    }
    assay_obj <- loaded_obj[[assay_name]]
    # Join only genuinely split layers (counts.1, data.sampleA...), not the standard counts/data/scale.data trio
    if (inherits(assay_obj, "Assay5") && any(grepl("^(counts|data)\\.", names(assay_obj@layers)))) {
      tryCatch({
        log_function(paste("Joining split layers in", assay_name))
        loaded_obj[[assay_name]] <- JoinLayers(assay_obj)
      }, error = function(e) log_function(paste("JoinLayers failed for", assay_name, ":", conditionMessage(e)[1])))
    }
  }
  DefaultAssay(loaded_obj) <- if ("Spatial" %in% names(loaded_obj@assays)) "Spatial" else names(loaded_obj@assays)[1]
  if (!"percent.mt" %in% colnames(loaded_obj@meta.data)) {
    log_function("Adding mitochondrial percentage...")
    loaded_obj[["percent.mt"]] <- PercentageFeatureSet(loaded_obj, pattern = get_mt_pattern_spatial(species))
  }
  log_function(paste("Loaded:", ncol(loaded_obj), "spots,", nrow(loaded_obj), "genes"))
  log_function(paste("Assays:", paste(names(loaded_obj@assays), collapse = ", ")))
  log_function(paste("Reductions:", paste(names(loaded_obj@reductions), collapse = ", ")))
  return(loaded_obj)
}

############################## Coordinates and physical calibration ##############################

# Spot centroids in full-resolution pixels with unambiguous column names (barcode, imagerow, imagecol).
# VisiumV1 exposes imagerow / imagecol directly. VisiumV2 exposes centroids as x / y (barcodes in a "cell" column) and the
# meaning of x depends on the Seurat version: up to Seurat 5.3, x = imagerow (vertical axis); newer versions store x = imagecol
# (horizontal axis) and flag it with coords_x_orientation = "horizontal" (UpdateSeuratObject converts older objects and sets the flag).
# scale = NULL requests unscaled (full-resolution) coordinates for both classes.
# @misc$spot_coords_transposed = TRUE swaps the two axes: fallback used at load time only when the stored coordinates cannot be
# corrected by transpose_spot_coordinates().
get_spot_coordinates_spatial <- function(obj, image = NULL) {
  if (length(obj@images) == 0) stop("No spatial image in object")
  image <- image %||% names(obj@images)[1]
  img <- obj@images[[image]]
  coords <- as.data.frame(Seurat::GetTissueCoordinates(img, scale = NULL))
  if (all(c("imagerow", "imagecol") %in% colnames(coords))) {
    out <- data.frame(barcode = rownames(coords), imagerow = coords$imagerow, imagecol = coords$imagecol, stringsAsFactors = FALSE)
  } else if (all(c("x", "y") %in% colnames(coords))) {
    barcodes <- if ("cell" %in% colnames(coords)) coords$cell else rownames(coords)
    # Horizontal x axis (current Seurat convention) only when the flag says so; objects without the flag use the former convention
    x_horizontal <- methods::.hasSlot(img, "coords_x_orientation") && identical(as.character(methods::slot(img, "coords_x_orientation")), "horizontal")
    out <- if (x_horizontal) {
      data.frame(barcode = barcodes, imagerow = coords$y, imagecol = coords$x, stringsAsFactors = FALSE)
    } else {
      data.frame(barcode = barcodes, imagerow = coords$x, imagecol = coords$y, stringsAsFactors = FALSE)
    }
  } else {
    stop("Unrecognized tissue coordinate columns: ", paste(colnames(coords), collapse = ", "))
  }
  # Correction set at load time when the coordinates were found transposed relative to the tissue image
  if (isTRUE(obj@misc$spot_coords_transposed)) out[, c("imagerow", "imagecol")] <- out[, c("imagecol", "imagerow")]
  # Keep only spots present in the object (guards against image centroids out of sync with the cells)
  out <- out[out$barcode %in% colnames(obj), , drop = FALSE]
  rownames(out) <- out$barcode
  return(out)
}

# Availability check used before spatial cluster plots; relies on get_spot_coordinates_spatial for barcode handling
check_spatial_coordinates <- function(obj) {
  if (length(obj@images) == 0) return(list(available = FALSE, message = "No spatial images found"))
  coords <- tryCatch(get_spot_coordinates_spatial(obj), error = function(e) e)
  if (inherits(coords, "error")) return(list(available = FALSE, message = paste("Error accessing coordinates:", conditionMessage(coords))))
  if (nrow(coords) == 0) return(list(available = FALSE, message = "No matching coordinates for cells"))
  return(list(available = TRUE, message = paste("Coordinates available for", nrow(coords), "spots"), coords = coords))
}

# Micrometres per full-resolution pixel. Space Ranger's spot_diameter_fullres spans a theoretical 65 um spot for standard
# Visium (10x converts pixels with 65 / spot_diameter_fullres; the printed spot itself is 55 um), and one bin for Visium HD.
# With 55 / spot_diameter_fullres every distance was 15 % too short (grid pitch measured at 84 um instead of 100 um).
get_um_per_pixel_spatial <- function(obj, spot_size_um, image = NULL) {
  image <- image %||% names(obj@images)[1]
  spot_px <- obj@images[[image]]@scale.factors$spot
  if (is.null(spot_px) || !is.finite(spot_px) || spot_px <= 0) stop("spot_diameter_fullres is missing from the image scale factors")
  reference_um <- if (isTRUE(spot_size_um >= 50)) 65 else spot_size_um
  return(reference_um / spot_px)
}

# Spot centroids in micrometres (x = image column, y = image row) plus calibration diagnostics.
# pitch_um is the median nearest-neighbour distance: ~100 um for Visium, ~bin size for Visium HD.
get_spot_positions_um <- function(obj, spot_size_um, image = NULL) {
  coords <- get_spot_coordinates_spatial(obj, image)
  um_per_px <- get_um_per_pixel_spatial(obj, spot_size_um, image)
  xy <- cbind(x = coords$imagecol, y = coords$imagerow) * um_per_px
  rownames(xy) <- coords$barcode
  pitch_um <- if (nrow(xy) > 1) stats::median(RANN::nn2(xy, k = 2)$nn.dists[, 2]) else NA_real_
  return(list(barcodes = coords$barcode, xy = xy, um_per_px = um_per_px, pitch_um = pitch_um))
}

# Physical spot / bin size (um) of a spatial object, for micrometre conversions.
# Standard Visium spots (55 um) sit on a 100 um grid, so the grid pitch is ~1.8 spot diameters; Visium HD bins touch
# (pitch = 1 bin diameter). For Visium HD the bin size comes from a valid stored value (2, 8 or 16 um), else from the
# names Space Ranger / Seurat give to barcodes, images or assays (s_008um_..., slice1.008um, Spatial.008um).
# Returns list(size_um, technology, source); size_um is NA for a Visium HD object whose bin size cannot be determined.
infer_spot_size_um <- function(obj) {
  stored <- suppressWarnings(as.numeric(obj@misc$spot_size_um))[1]
  coords <- get_spot_coordinates_spatial(obj)
  spot_px <- obj@images[[1]]@scale.factors$spot
  if (nrow(coords) < 2 || is.null(spot_px) || !is.finite(spot_px) || spot_px <= 0) {
    return(list(size_um = if (isTRUE(stored > 0)) stored else NA_real_, technology = "unknown", source = "no grid information"))
  }
  # Pitch measured on the first 5,000 spots (deterministic, fast on Visium HD)
  xy <- cbind(coords$imagecol, coords$imagerow)
  query <- xy[seq_len(min(nrow(xy), 5000)), , drop = FALSE]
  pitch_px <- stats::median(RANN::nn2(xy, query, k = 2)$nn.dists[, 2])
  if (pitch_px / spot_px > 1.4) return(list(size_um = 55, technology = "Visium", source = "spot grid"))
  if (isTRUE(stored %in% c(2, 8, 16))) return(list(size_um = stored, technology = "Visium HD", source = "stored value"))
  tags <- c(utils::head(colnames(obj), 20), names(obj@images), names(obj@assays))
  hit <- regmatches(tags, regexpr("0(02|08|16)um", tags))
  if (length(hit) > 0) return(list(size_um = as.numeric(sub("um$", "", hit[1])), technology = "Visium HD", source = paste("name", hit[1])))
  return(list(size_um = NA_real_, technology = "Visium HD", source = "unknown"))
}

# Fraction of spots lying on tissue in the tissue image, for spot centroids given in full-resolution pixels.
# tissue: result of get_tissue_image_spatial(). Tissue = pixels clearly darker than the background (median grey level of the
# image minus 15, the image being mostly background); spots outside the image count as off tissue.
spots_on_tissue_fraction <- function(tissue, col, row) {
  gray <- magick::image_data(magick::image_convert(tissue$img, colorspace = "gray"), channels = "gray")
  g <- matrix(as.integer(gray[1, , ]), nrow = dim(gray)[2])
  threshold <- stats::median(g) - 15
  px <- floor(col * tissue$scalef) + 1; py <- floor(row * tissue$scalef) + 1
  inside <- px >= 1 & px <= nrow(g) & py >= 1 & py <= ncol(g)
  on <- rep(FALSE, length(px))
  on[inside] <- g[cbind(px[inside], py[inside])] < threshold
  return(mean(on))
}

# Check that the spot coordinates are not transposed relative to the tissue image (x / y convention lost between Seurat
# versions or when an object is saved and reloaded): compares the fraction of spots on tissue with the coordinates as read
# and with x and y swapped. transposed = TRUE when the swapped version fits clearly better (+15 points of spots on tissue).
# The check ignores any previous correction stored in @misc$spot_coords_transposed.
detect_spot_transposition <- function(obj) {
  obj@misc$spot_coords_transposed <- NULL
  coords <- get_spot_coordinates_spatial(obj)
  tissue <- get_tissue_image_spatial(obj)
  as_read <- spots_on_tissue_fraction(tissue, coords$imagecol, coords$imagerow)
  swapped <- spots_on_tissue_fraction(tissue, coords$imagerow, coords$imagecol)
  return(list(transposed = swapped > as_read + 0.15, on_tissue_as_read = as_read, on_tissue_swapped = swapped))
}

# Swap the two axes of the stored spot centroids of every tissue image (VisiumV1 imagerow / imagecol, VisiumV2 centroid
# x / y and bounding box), for objects whose coordinates are transposed relative to the tissue image. Correcting the stored
# coordinates keeps Seurat plots (SpatialFeaturePlot, SpatialDimPlot) and the Cell-Hub views aligned on the same image.
transpose_spot_coordinates <- function(obj) {
  for (image in names(obj@images)) {
    img <- obj@images[[image]]
    if (methods::.hasSlot(img, "coordinates") && all(c("imagerow", "imagecol") %in% colnames(img@coordinates))) {
      co <- img@coordinates
      co[, c("imagerow", "imagecol")] <- co[, c("imagecol", "imagerow")]
      img@coordinates <- co
    } else if (methods::.hasSlot(img, "boundaries") && !is.null(img@boundaries$centroids)) {
      cents <- img@boundaries$centroids
      xy <- methods::slot(cents, "coords")
      xy[, c(1, 2)] <- xy[, c(2, 1)]
      methods::slot(cents, "coords", check = FALSE) <- xy
      if (methods::.hasSlot(cents, "bbox")) {
        bb <- methods::slot(cents, "bbox")
        bb[c(1, 2), ] <- bb[c(2, 1), ]
        methods::slot(cents, "bbox", check = FALSE) <- bb
      }
      img@boundaries$centroids <- cents
    } else {
      stop("Unsupported tissue image class: ", class(img)[1])
    }
    obj@images[[image]] <- img
  }
  return(obj)
}

############################## Tissue images: stored images and custom upload (magick / VIPS) ##############################

# TRUE when the vips command-line tool is installed and runs.
# Sys.which is checked first: system2() captures the shell's "vips: not found" message, which would otherwise look like valid output.
ensure_vips <- function() {
  if (!nzchar(Sys.which("vips"))) return(FALSE)
  out <- suppressWarnings(tryCatch(system2("vips", "--version", stdout = TRUE, stderr = TRUE), error = function(e) character()))
  return(is.null(attr(out, "status")) && length(out) > 0)
}

# Image width and height in pixels, read from the header only (vipsheader), with a magick fallback that loads the image.
# Returns c(NA, NA) when no method can read the file.
vips_tiff_size <- function(path) {
  # 1) vipsheader -a lists every header field, including "width: N" and "height: N"
  out <- suppressWarnings(tryCatch(system2("vipsheader", c("-a", shQuote(path)), stdout = TRUE), error = function(e) character()))
  if (length(out)) {
    w <- suppressWarnings(as.numeric(sub(".*width:\\s*", "", out[grep("width:", out, ignore.case = TRUE)[1]])))
    h <- suppressWarnings(as.numeric(sub(".*height:\\s*", "", out[grep("height:", out, ignore.case = TRUE)[1]])))
    if (is.finite(w) && is.finite(h) && w > 0 && h > 0) return(c(w, h))
  }
  # 2) Plain vipsheader, typical line: "in.tif: 100000x80000 uchar, 3 bands, srgb, tiffload"
  out2 <- suppressWarnings(tryCatch(system2("vipsheader", shQuote(path), stdout = TRUE), error = function(e) character()))
  if (length(out2)) {
    m <- regmatches(out2, regexpr("\\b(\\d+)x(\\d+)\\b", out2))
    if (length(m) && nzchar(m[1])) {
      nums <- as.numeric(strsplit(m[1], "x")[[1]])
      if (length(nums) == 2 && all(is.finite(nums)) && all(nums > 0)) return(nums)
    }
  }
  # 3) Fallback (slow on huge files): read the image with magick
  img <- tryCatch(magick::image_read(path), error = function(e) NULL)
  if (!is.null(img)) {
    info <- magick::image_info(img)
    w <- as.numeric(info$width[1]); h <- as.numeric(info$height[1])
    if (is.finite(w) && is.finite(h) && w > 0 && h > 0) return(c(w, h))
  }
  return(c(NA_real_, NA_real_))
}

# Display one of the Space Ranger images stored at load time (see collect_tissue_images) in the first slice.
# Only one image is decoded at a time; the other stays compressed in @misc$tissue_images.
# scale.factors$lowres is set to the displayed image's factor because Seurat plots map full-resolution
# coordinates with it by default (image.scale = "lowres"), so spots stay aligned with either image.
set_spatial_image_resolution <- function(obj, resolution = "lowres", log_function = message) {
  if (length(obj@images) == 0) stop("No spatial image slot in this object")
  stored <- obj@misc$tissue_images
  if (is.null(stored[[resolution]])) {
    stop(sprintf("No %s image stored in this object (stored: %s). Include tissue_%s_image.png in the ZIP and reload it.",
                 resolution, if (length(stored)) paste(names(stored), collapse = ", ") else "none", resolution))
  }
  slice <- names(obj@images)[1]
  img <- png::readPNG(stored[[resolution]]$png)
  # Seurat expects a height x width x 3 array: drop the alpha channel if the PNG has one
  if (length(dim(img)) == 3 && dim(img)[3] > 3) img <- img[, , 1:3]
  img_obj <- obj@images[[slice]]
  img_obj@image <- img
  img_obj@scale.factors$lowres <- stored[[resolution]]$scalef
  # Restore the genuine hires factor (replace_spatial_image overwrites it with a custom ratio)
  if (!is.null(stored$hires)) img_obj@scale.factors$hires <- stored$hires$scalef
  obj@images[[slice]] <- img_obj
  obj@misc$active_tissue_image <- resolution
  log_function(sprintf("Tissue image set to %s: %d x %d px (scale factor %.4f)", resolution, dim(img)[2], dim(img)[1], stored[[resolution]]$scalef))
  return(obj)
}

# Replace the tissue image of the first slice with the FULL-RESOLUTION image Space Ranger aligned the spots to.
# The image is downsampled to max_side pixels and lowres/hires scale factors are set to the downsampling ratio,
# so spot coordinates (stored in full-resolution pixels) stay aligned. Any other image frame cannot be aligned this way.
replace_spatial_image <- function(obj, img_path, max_side = 2000, log_function = message) {
  if (length(obj@images) == 0) stop("No spatial image slot to replace")
  slice <- names(obj@images)[1]
  tmp_png <- tempfile(fileext = ".png"); on.exit(unlink(tmp_png), add = TRUE)
  if (ensure_vips()) {
    # VIPS streams the file, so whole-slide TIFFs are downsampled without being loaded in memory
    full_dims <- vips_tiff_size(img_path)
    out <- system2("vips", c("thumbnail", shQuote(img_path), shQuote(tmp_png), max_side, "--height", max_side, "--size", "down"), stdout = TRUE, stderr = TRUE)
    if (!file.exists(tmp_png)) stop("VIPS thumbnail failed: ", paste(out, collapse = " "))
    img <- magick::image_read(tmp_png)
  } else {
    img <- magick::image_read(img_path)
    info0 <- magick::image_info(img); full_dims <- c(as.numeric(info0$width[1]), as.numeric(info0$height[1]))
    if (max(full_dims) > max_side) img <- magick::image_resize(img, sprintf("%dx%d", max_side, max_side))
  }
  if (any(!is.finite(full_dims))) stop("Could not read image dimensions")
  coords <- get_spot_coordinates_spatial(obj, slice)
  # Spots are stored in full-resolution pixels: they must fall inside the uploaded frame
  if (max(coords$imagecol) > full_dims[1] || max(coords$imagerow) > full_dims[2]) {
    stop(sprintf("Spots extend to %.0f x %.0f px but the image is %.0f x %.0f px - upload the full-resolution image used by Space Ranger",
                 max(coords$imagecol), max(coords$imagerow), full_dims[1], full_dims[2]))
  }
  info <- magick::image_info(img)
  ratio <- as.numeric(info$width[1]) / full_dims[1]
  img_obj <- obj@images[[slice]]
  # magick bitmap -> height x width x 3 array scaled to [0, 1], the layout Seurat stores in @image
  img_obj@image <- as.integer(magick::image_data(img, channels = "rgb")) / 255
  # Seurat maps full-resolution coordinates onto the image with these factors (image.scale = "lowres" by default)
  img_obj@scale.factors$lowres <- ratio
  img_obj@scale.factors$hires <- ratio
  obj@images[[slice]] <- img_obj
  # Flag the custom image so the lowres / hires selector shows no stored image as active
  obj@misc$active_tissue_image <- "custom"
  log_function(sprintf("Tissue image replaced: %d x %d px (scale factor %.4f)", info$width[1], info$height[1], ratio))
  return(obj)
}

############################## Tissue sections on one capture area ##############################

# Separate tissue sections lying on the same capture area: spots closer than max_gap_pitch times the grid pitch
# (median nearest-neighbour distance) are connected, and each connected group of spots is one section.
# Groups smaller than min_spots are labelled "Debris" (tissue fragments, isolated spots). Sections are numbered
# S1, S2, ... from left to right (then top to bottom) by the position of their centre on the image. Deterministic.
# Returns the section of every spot (named by barcode), a per-section summary and the grid pitch in pixels.
detect_tissue_sections <- function(sp_obj, max_gap_pitch = 2, min_spots = 50) {
  if (!requireNamespace("RANN", quietly = TRUE) || !requireNamespace("igraph", quietly = TRUE)) stop("Packages 'RANN' and 'igraph' are required")
  max_gap_pitch <- as.numeric(max_gap_pitch); min_spots <- as.integer(min_spots)
  if (!is.finite(max_gap_pitch) || max_gap_pitch < 1) stop("The maximum gap must be at least 1 grid pitch")
  coords <- get_spot_coordinates_spatial(sp_obj)
  n <- nrow(coords)
  if (n < 3) stop("Too few spots to detect sections")
  xy <- cbind(coords$imagecol, coords$imagerow)
  pitch <- stats::median(RANN::nn2(xy, xy, k = 2)$nn.dists[, 2])
  # Upper bound on grid points inside the search disc, so radius searches are not truncated
  k <- min(n, ceiling(4 * (max_gap_pitch + 1)^2) + 1)
  nn <- RANN::nn2(xy, xy, k = k, searchtype = "radius", radius = max_gap_pitch * pitch * (1 + 1e-6))
  keep <- as.vector(nn$nn.idx) > 0
  adjacency <- Matrix::sparseMatrix(i = rep(seq_len(n), k)[keep], j = as.vector(nn$nn.idx)[keep], x = 1, dims = c(n, n))
  graph <- igraph::graph_from_adjacency_matrix(adjacency, mode = "max", diag = FALSE)
  component <- igraph::components(graph)$membership
  sizes <- tabulate(component)
  big <- which(sizes >= min_spots)
  if (length(big) == 0) stop(sprintf("No group of at least %d contiguous spots: lower the minimum section size", min_spots))
  # Sections ordered left to right, then top to bottom, by their centre
  centre_x <- tapply(xy[, 1], component, mean)[as.character(big)]
  centre_y <- tapply(xy[, 2], component, mean)[as.character(big)]
  big <- big[order(round(centre_x / (10 * pitch)), centre_y)]
  labels <- stats::setNames(paste0("S", seq_along(big)), big)
  section <- ifelse(component %in% big, labels[as.character(component)], "Debris")
  section <- factor(section, levels = c(paste0("S", seq_along(big)), if (any(!component %in% big)) "Debris"))
  names(section) <- coords$barcode
  summary <- data.frame(section = levels(section), n_spots = as.integer(table(section)),
                        centre_x = round(as.numeric(tapply(xy[, 1], section, mean))), centre_y = round(as.numeric(tapply(xy[, 2], section, mean))),
                        stringsAsFactors = FALSE)
  return(list(section = section, summary = summary, pitch_px = pitch, params = list(max_gap_pitch = max_gap_pitch, min_spots = min_spots)))
}

# Store the detected sections in the spot metadata: section_detected keeps the detection labels (S1, S2, ..., Debris),
# section holds the names given by the user (initially the same); sample is reset to the section name and condition to NA
# until the user assigns them. Spots absent from the coordinates get NA.
add_tissue_sections_to_object <- function(sp_obj, sections) {
  sp_obj$section_detected <- factor(as.character(sections$section)[match(colnames(sp_obj), names(sections$section))], levels = levels(sections$section))
  sp_obj$section <- sp_obj$section_detected
  sp_obj$sample <- as.character(sp_obj$section)
  sp_obj$condition <- NA_character_
  sp_obj@misc$tissue_sections <- c(sections$params, list(summary = sections$summary, pitch_px = sections$pitch_px))
  return(sp_obj)
}

# Apply the user's assignment: assignment is a data.frame with one row per detected section (S1, S2, ..., Debris) and columns
# section, name, sample, condition. Giving two sections the same name merges them (e.g. a section split by a tear).
# Sections named "Debris" (or left empty) get the Debris label. The assignment can be changed and applied again.
assign_tissue_sections <- function(sp_obj, assignment) {
  if (!"section_detected" %in% colnames(sp_obj@meta.data)) stop("Detect the tissue sections first")
  assignment$name <- trimws(as.character(assignment$name))
  assignment$name[is.na(assignment$name) | assignment$name == ""] <- "Debris"
  current <- as.character(sp_obj$section_detected)
  idx <- match(current, assignment$section)
  if (anyNA(idx[!is.na(current)])) stop("Some sections of the object are missing from the assignment table")
  clean <- function(v) { v <- trimws(as.character(v)); v[is.na(v) | v == ""] <- NA_character_; v }
  new_names <- assignment$name[idx]
  sp_obj$section <- factor(new_names, levels = unique(c(setdiff(assignment$name, "Debris"), if (any(new_names == "Debris", na.rm = TRUE)) "Debris")))
  sample_values <- clean(assignment$sample)[idx]
  sp_obj$sample <- ifelse(is.na(sample_values), new_names, sample_values)
  sp_obj$condition <- clean(assignment$condition)[idx]
  sp_obj@misc$tissue_sections$assignment <- assignment
  return(sp_obj)
}

# Spots, samples and conditions per section, for the summary table
summarize_tissue_sections <- function(sp_obj) {
  if (!"section" %in% colnames(sp_obj@meta.data)) return(NULL)
  md <- sp_obj@meta.data
  out <- do.call(rbind, lapply(levels(droplevels(md$section)), function(s) {
    rows <- md$section == s & !is.na(md$section)
    data.frame(section = s, n_spots = sum(rows), sample = paste(unique(stats::na.omit(md$sample[rows])), collapse = ", "),
               condition = paste(unique(stats::na.omit(md$condition[rows])), collapse = ", "),
               median_umi = round(stats::median(md$nCount_Spatial[rows])), stringsAsFactors = FALSE)
  }))
  return(out)
}

# Tissue section of each spot (NULL when sections have not been detected), used to keep spot graphs inside sections
spot_sections <- function(sp_obj, barcodes) {
  if (!"section" %in% colnames(sp_obj@meta.data)) return(NULL)
  return(as.character(sp_obj@meta.data[barcodes, "section"]))
}

############################## Quality control ##############################

# Data-driven QC suggestions, used only to help the user when no spot passes their thresholds.
# Lower bounds sit at half the first quartile of non-empty spots; no upper bound on genes/UMIs,
# because high-UMI spots in spatial data reflect dense tissue rather than doublets.
analyze_data_distribution <- function(obj) {
  md <- obj@meta.data
  features <- md$nFeature_Spatial; counts <- md$nCount_Spatial; mt <- md$percent.mt
  features[is.na(features)] <- 0; counts[is.na(counts)] <- 0; mt[is.na(mt)] <- 0
  suggested_filters <- list(
    min_features = max(1, round(stats::quantile(features[features > 0], 0.25, na.rm = TRUE) * 0.5)), max_features = max(features),
    min_counts = max(1, round(stats::quantile(counts[counts > 0], 0.25, na.rm = TRUE) * 0.5)), max_counts = max(counts),
    max_mt = min(round(stats::quantile(mt, 0.95, na.rm = TRUE)), 50))
  would_pass <- features >= suggested_filters$min_features & counts >= suggested_filters$min_counts & mt <= suggested_filters$max_mt
  return(list(suggested_filters = suggested_filters, spots_passing = sum(would_pass, na.rm = TRUE)))
}

############################## Clusters ##############################

# Metadata column holding spatial clusters (annotated names take precedence).
# Spatial-specific: multiple_datasets_server.R has its own get_cluster_column() with a different priority order.
get_cluster_column_spatial <- function(obj) {
  cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) "annotated_clusters" else "seurat_clusters"
  return(cluster_col)
}

# Cluster labels of the current spatial cluster column (factor levels, or sorted values for character columns)
get_cluster_names_spatial <- function(obj) {
  col <- get_cluster_column_spatial(obj)
  if (!col %in% colnames(obj@meta.data)) return(character(0))
  vals <- obj@meta.data[[col]]
  cluster_names <- if (is.factor(vals)) levels(droplevels(vals)) else sort(unique(as.character(vals)))
  return(cluster_names)
}

# Per-spot cluster labels (named by barcode), ordered cluster names and the column used.
# Objects without clusters get a single "No_clusters" group so co-localization still runs.
get_cluster_info <- function(spatial_obj) {
  cluster_col <- get_cluster_column_spatial(spatial_obj)
  if (!cluster_col %in% colnames(spatial_obj@meta.data)) {
    spots_clusters <- setNames(rep("No_clusters", ncol(spatial_obj)), colnames(spatial_obj))
    return(list(spots = spots_clusters, names = "No_clusters", column_used = NULL))
  }
  spots_clusters <- setNames(as.character(spatial_obj@meta.data[[cluster_col]]), colnames(spatial_obj))
  return(list(spots = spots_clusters, names = get_cluster_names_spatial(spatial_obj), column_used = cluster_col))
}

############################## Plot themes ##############################

# Black background theme for the annotation plots (dark mode)
.dark_theme_bg <- function(base_size = 14) {
  dark_theme <- theme_void(base_size = base_size) +
    theme(plot.background = element_rect(fill = "black", color = NA), panel.background = element_rect(fill = "black", color = NA),
          legend.background = element_rect(fill = "black", color = NA), legend.key = element_rect(fill = "black", color = NA),
          legend.text = element_text(color = "white"), legend.title = element_text(color = "white"),
          plot.title = element_text(color = "white", hjust = 0.5))
  return(dark_theme)
}

############################## Interactive tissue viewer ##############################
# Windows are c(xmin, xmax, ymin, ymax) in full-resolution image pixels, y running top-down (image rows).
# The viewer crops the tissue image stored in the object to the current window and redraws the spots in it (server-side zoom).

# Keep the viewing window inside image bounds and enforce a minimum size (fraction of the image, at least 16 px)
clamp_view <- function(win, w, h, min_frac = 0.02) {
  if (is.null(win)) return(c(0, w, 0, h))
  xmin <- max(0, min(win[1], win[2])); xmax <- min(w, max(win[1], win[2]))
  ymin <- max(0, min(win[3], win[4])); ymax <- min(h, max(win[3], win[4]))
  min_w <- max(16, w * min_frac); min_h <- max(16, h * min_frac)
  if ((xmax - xmin) < min_w) { cx <- (xmin + xmax) / 2; xmin <- cx - min_w / 2; xmax <- cx + min_w / 2 }
  if ((ymax - ymin) < min_h) { cy <- (ymin + ymax) / 2; ymin <- cy - min_h / 2; ymax <- cy + min_h / 2 }
  xmin <- max(0, xmin); ymin <- max(0, ymin); xmax <- min(w, xmax); ymax <- min(h, ymax)
  return(c(xmin, xmax, ymin, ymax))
}

# Window centred on (x, y) whose width and height are the current ones divided by factor (factor > 1 zooms in)
zoom_center <- function(win, w, h, x, y, factor) {
  if (is.null(win)) win <- c(0, w, 0, h)
  cur_w <- (win[2] - win[1]) / factor
  cur_h <- (win[4] - win[3]) / factor
  return(c(x - cur_w / 2, x + cur_w / 2, y - cur_h / 2, y + cur_h / 2))
}

# Crop a magick image to a window given in image pixels c(xmin, xmax, ymin, ymax) and downscale it so its longest side fits display_size
prepare_display_image <- function(img_full, view_window, display_size) {
  crop_geo <- sprintf("%dx%d+%d+%d", as.integer(view_window[2] - view_window[1]), as.integer(view_window[4] - view_window[3]),
                      as.integer(view_window[1]), as.integer(view_window[3]))
  img_crop <- magick::image_crop(img_full, geometry = crop_geo)
  info <- magick::image_info(img_crop)
  if (max(info$width, info$height) > display_size) {
    scale <- display_size / max(info$width, info$height)
    img_crop <- magick::image_scale(img_crop, geometry = sprintf("%dx%d", as.integer(info$width * scale), as.integer(info$height * scale)))
  }
  return(img_crop)
}

# Move a window by (dx, dy) full-resolution pixels, keeping its size and keeping it inside the image (panning, zoom out at the edges).
# A window larger than the image is reduced to the image size.
pan_view <- function(win, w, h, dx = 0, dy = 0) {
  ww <- min(win[2] - win[1], w); hh <- min(win[4] - win[3], h)
  xmin <- min(max(0, win[1] + dx), w - ww); ymin <- min(max(0, win[3] + dy), h - hh)
  return(c(xmin, xmin + ww, ymin, ymin + hh))
}

# Best tissue image available in the object, as a magick image, with the factor mapping full-resolution pixels to its pixels.
# Priority: stored hires PNG, stored lowres PNG (@misc$tissue_images, see collect_tissue_images), then the image displayed by Seurat.
# width_full / height_full give the image extent in full-resolution pixels (the coordinate system of the spots).
get_tissue_image_spatial <- function(obj) {
  if (length(obj@images) == 0) stop("No spatial image in this object")
  stored <- obj@misc$tissue_images
  res <- intersect(c("hires", "lowres"), names(stored))
  if (length(res) > 0) {
    img <- magick::image_read(stored[[res[1]]]$png); scalef <- stored[[res[1]]]$scalef; source <- res[1]
  } else {
    # Displayed image: scale.factors$lowres always describes it (see read_visium_image / set_spatial_image_resolution).
    # The height x width x 3 array is encoded once as an in-memory PNG: as.raster() on a 6000 px image builds a matrix of
    # ~36 million colour strings and takes minutes, the PNG round trip takes seconds.
    img_obj <- obj@images[[1]]
    img <- magick::image_read(png::writePNG(img_obj@image)); scalef <- img_obj@scale.factors$lowres; source <- "object"
  }
  if (is.null(scalef) || !is.finite(scalef) || scalef <= 0) stop("Missing scale factor for the tissue image")
  info <- magick::image_info(img)
  width_px <- as.numeric(info$width[1]); height_px <- as.numeric(info$height[1])
  return(list(img = img, scalef = scalef, source = source, width_px = width_px, height_px = height_px,
              width_full = width_px / scalef, height_full = height_px / scalef))
}

# Crop the tissue image to a window given in full-resolution pixels and downscale it so its longest side fits display_size.
# The crop is snapped to whole image pixels, so its extent (returned in full-resolution pixels) can slightly exceed the window.
# Returns NULL when the window lies outside the image.
crop_tissue_image <- function(tissue, view, display_size = 1500) {
  s <- tissue$scalef
  xmin <- max(0, floor(view[1] * s)); xmax <- min(tissue$width_px, ceiling(view[2] * s))
  ymin <- max(0, floor(view[3] * s)); ymax <- min(tissue$height_px, ceiling(view[4] * s))
  if (xmax - xmin < 1 || ymax - ymin < 1) return(NULL)
  img <- prepare_display_image(tissue$img, c(xmin, xmax, ymin, ymax), display_size)
  return(list(img = img, extent = c(xmin, xmax, ymin, ymax) / s))
}

# Tissue viewer plot: cropped tissue image with the spots of the window drawn at their physical size for the current zoom.
# spots: data.frame with x, y (full-resolution pixels, y = image row) and value (factor = clusters, numeric = gene or QC metric),
# or NULL to draw the tissue image alone.
# view: window c(xmin, xmax, ymin, ymax); crop: result of crop_tissue_image() or NULL.
# panel_px: plot size in points (1 pt = 1 px in Shiny at 72 dpi), used to turn the spot diameter into a ggplot point size.
# The y axis is the negated image row, so the tissue keeps its orientation without a reversed scale (brush coordinates stay simple).
build_tissue_viewer_plot <- function(spots, view, crop = NULL, value_label = "", cluster_colors = NULL, highlight = NULL,
                                     spot_diameter_px = 10, spot_shape = 16, size_factor = 1, spot_alpha = 0.8, image_alpha = 1,
                                     show_labels = FALSE, show_legend = TRUE, panel_px = c(1000, 750)) {
  draw_spots <- !is.null(spots)
  # Window corners as a blank layer: the x / y scales exist even when no spot is drawn (tissue-only view, brush zoom still works)
  p <- ggplot(data.frame(x = view[1:2], y_plot = -view[3:4]), aes(x = x, y = y_plot)) + geom_blank()
  if (!is.null(crop)) {
    img <- crop$img
    # Lower opacity = fade the tissue towards white
    if (image_alpha < 1) img <- magick::image_colorize(img, opacity = round((1 - image_alpha) * 100), color = "white")
    p <- p + annotation_raster(grDevices::as.raster(img), xmin = crop$extent[1], xmax = crop$extent[2],
                               ymin = -crop$extent[4], ymax = -crop$extent[3], interpolate = TRUE)
  }
  if (draw_spots) {
    spots$y_plot <- -spots$y
    is_discrete <- is.factor(spots$value)
    # Continuous values: highest values drawn last so they stay visible on top
    if (!is_discrete) spots <- spots[order(spots$value, na.last = FALSE), , drop = FALSE]
    # Points per full-resolution pixel: coord_fixed keeps the aspect ratio, so the tighter dimension sets the scale
    usable_w <- panel_px[1] * (if (show_legend) 0.82 else 0.95); usable_h <- panel_px[2] * 0.95
    pt_per_px <- min(usable_w / (view[2] - view[1]), usable_h / (view[4] - view[3]))
    point_size <- min(max(spot_diameter_px * pt_per_px * size_factor / ggplot2::.pt, 0.05), 25)
    if (is_discrete) {
      if (!is.null(highlight) && highlight %in% levels(spots$value)) {
        # Highlighted cluster in color, every other spot in light grey
        in_hl <- spots$value == highlight
        p <- p + geom_point(data = spots[!in_hl, , drop = FALSE], color = "#BDBDBD", alpha = spot_alpha * 0.4, size = point_size, shape = spot_shape) +
          geom_point(data = spots[in_hl, , drop = FALSE], aes(color = value), alpha = spot_alpha, size = point_size, shape = spot_shape)
      } else {
        p <- p + geom_point(data = spots, aes(color = value), alpha = spot_alpha, size = point_size, shape = spot_shape)
      }
      if (!is.null(cluster_colors)) p <- p + scale_color_manual(values = cluster_colors, name = value_label)
      p <- p + guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))
      if (show_labels && nrow(spots) > 0) {
        # Label each cluster at the median position of its spots in the window
        centers <- stats::aggregate(cbind(x, y_plot) ~ value, data = spots, FUN = stats::median)
        p <- p + geom_label(data = centers, aes(x = x, y = y_plot, label = value), inherit.aes = FALSE, size = 3.5, alpha = 0.8)
      }
    } else {
      p <- p + geom_point(data = spots, aes(color = value), alpha = spot_alpha, size = point_size, shape = spot_shape) +
        scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(11, "Spectral")), name = value_label)
    }
  }
  p <- p + coord_fixed(xlim = c(view[1], view[2]), ylim = c(-view[4], -view[3]), expand = FALSE) +
    theme_void() + theme(legend.position = if (show_legend && draw_spots) "right" else "none", plot.background = element_rect(fill = "white", colour = NA))
  return(p)
}

############################## Gene co-localization ##############################

# HTML bar showing a distance as a fraction of the maximum search distance (DT tables, escape = FALSE)
create_distance_bar <- function(distance_um, max_distance_um) {
  if (is.na(distance_um) || is.na(max_distance_um) || max_distance_um == 0) return("<span style='color: grey;'>N/A</span>")
  pct <- round(100 * distance_um / max_distance_um, 1)
  pct <- min(max(pct, 0), 100)
  # Green for short distances, orange for distances close to the threshold
  color <- if (pct < 25) "#28a745" else if (pct < 50) "#85c88a" else if (pct < 75) "#ffc107" else "#fd7e14"
  bar_html <- paste0(
    '<div style="width: 100%; height: 20px; border: 1px solid #ddd; ',
    'border-radius: 4px; overflow: hidden; background-color: #f8f9fa;">',
    '<div style="width: ', pct, '%; height: 100%; ',
    'background: linear-gradient(90deg, ', color, ' 0%, ',
    adjustcolor(color, alpha.f = 0.7), ' 100%); ',
    'display: flex; align-items: center; justify-content: center; ',
    'transition: width 0.3s ease;" ',
    'title="', distance_um, ' µm (', pct, '% of max distance)">',
    '<span style="font-size: 10px; color: white; font-weight: bold; ',
    'text-shadow: 1px 1px 2px rgba(0,0,0,0.5);">',
    distance_um, ' µm',
    '</span>',
    '</div></div>'
  )
  return(bar_html)
}

# Inter-cluster statistics for co-localized pairs (integer spot indices), with the percentage of reference spots involved.
# gene1_pct = % of reference-gene-positive spots of ref_cluster having >= 1 target neighbour in target_cluster
# (the column keeps this name in both directions because the server tables expect it).
calculate_inter_cluster_stats <- function(pairs, ref_spots, spot_clusters) {
  if (is.null(pairs) || nrow(pairs) == 0) return(data.frame())
  dt <- data.table::as.data.table(pairs)
  out <- dt[, list(n_pairs = .N, n_ref_spots = data.table::uniqueN(ref_spot), mean_distance_um = round(mean(distance_um), 1),
                   median_distance_um = round(stats::median(distance_um), 1), min_distance_um = round(min(distance_um), 1),
                   max_distance_um = round(max(distance_um), 1)), by = list(ref_cluster, target_cluster)]
  # Denominator: reference-gene-positive spots in the reference cluster
  ref_totals <- table(spot_clusters[ref_spots])
  out <- as.data.frame(out)
  out$gene1_pct <- round(100 * out$n_ref_spots / as.numeric(ref_totals[as.character(out$ref_cluster)]), 1)
  out <- out[order(-out$n_pairs), , drop = FALSE]
  return(data.frame(ref_cluster = as.character(out$ref_cluster), target_cluster = as.character(out$target_cluster), n_pairs = out$n_pairs,
                    gene1_pct = out$gene1_pct, mean_distance_um = out$mean_distance_um, median_distance_um = out$median_distance_um,
                    min_distance_um = out$min_distance_um, max_distance_um = out$max_distance_um,
                    same_cluster = as.character(out$ref_cluster) == as.character(out$target_cluster), stringsAsFactors = FALSE))
}

# Per-cluster summary of co-localized spots; each reference spot is counted once per category
# (a spot with neighbours both inside and outside its cluster counts once in the total).
# spot_clusters is a factor over all spots; pairs_12 / pairs_21 use integer spot indices (NULL when no pair).
calculate_cluster_summary <- function(g1_spots, g2_spots, pairs_12, pairs_21, spot_clusters) {
  cluster_names <- levels(spot_clusters)
  count_by_cluster <- function(spots) as.numeric(table(spot_clusters[spots]))
  # Unique reference spots per cluster among the pairs selected by keep()
  count_ref <- function(pairs, keep) {
    if (is.null(pairs) || nrow(pairs) == 0) return(rep(0, length(cluster_names)))
    return(count_by_cluster(unique(pairs$ref_spot[keep(pairs)])))
  }
  df <- data.frame(cluster = cluster_names, gene1_spots = count_by_cluster(g1_spots), gene2_spots = count_by_cluster(g2_spots), stringsAsFactors = FALSE)
  safe_pct <- function(num, den) ifelse(den > 0, round(100 * num / den, 1), 0)
  for (tag in c("gene1", "gene2")) {
    pairs <- if (tag == "gene1") pairs_12 else pairs_21
    same <- count_ref(pairs, function(p) p$same_cluster %in% TRUE)
    different <- count_ref(pairs, function(p) p$same_cluster %in% FALSE)
    total <- count_ref(pairs, function(p) rep(TRUE, nrow(p)))
    den <- df[[paste0(tag, "_spots")]]
    df[[paste0(tag, "_coloc_same_cluster")]] <- same
    df[[paste0(tag, "_coloc_different_cluster")]] <- different
    df[[paste0(tag, "_total_coloc")]] <- total
    df[[paste0(tag, "_coloc_rate_same")]] <- safe_pct(same, den)
    df[[paste0(tag, "_coloc_rate_diff")]] <- safe_pct(different, den)
    df[[paste0(tag, "_total_coloc_rate")]] <- safe_pct(total, den)
  }
  return(df)
}

# Distance-based co-localization between gene-positive spots, in micrometres.
# Full-resolution pixel coordinates are converted with spot_size_um / spot_diameter_fullres (55 um Visium, bin size Visium HD).
# 2 genes: directional pair statistics (gene1 -> gene2, optionally gene2 -> gene1).
# 3+ genes: a spot of gene g counts when every other gene has a positive spot within max_distance_um.
# include_same_spot = TRUE counts co-expression within the same spot (distance 0) as co-localization.
# Pair-level tables above max_pairs_returned rows are not returned (Excel row limit, memory on Visium HD).
# The assay must cover every spot: a sketched assay (random subset) is refused, as for all spatial views.
# n_permutations > 0 (2 genes): the positive spots of the target gene are placed at random among all spots, keeping their number,
# to give the rate expected without spatial association, the enrichment (observed / expected) and an empirical one-sided
# p-value (1 + #null >= observed) / (1 + n_permutations). The seed makes the permutations reproducible.
analyze_gene_colocalization <- function(sp_obj, gene_list, max_distance_um = 50, spot_size_um = 55, assay = "Spatial",
                                        bidirectional = TRUE, include_same_spot = TRUE, max_pairs_returned = 1e6,
                                        n_permutations = 0, seed = 42) {
  gene_list <- unique(gene_list); max_distance_um <- as.numeric(max_distance_um); spot_size_um <- as.numeric(spot_size_um)
  n_permutations <- suppressWarnings(as.integer(n_permutations))
  if (length(n_permutations) != 1 || is.na(n_permutations) || n_permutations < 0) n_permutations <- 0L
  stopifnot(length(gene_list) >= 2, is.finite(max_distance_um), max_distance_um >= 0, is.finite(spot_size_um), spot_size_um > 0)
  if (!requireNamespace("RANN", quietly = TRUE)) stop("Package 'RANN' is required for the neighbour search")
  if (length(sp_obj@images) == 0) stop("No spatial image / coordinates.")
  if (!assay %in% names(sp_obj@assays)) stop("Assay not found: ", assay)
  missing_genes <- setdiff(gene_list, rownames(sp_obj[[assay]]))
  if (length(missing_genes) > 0) stop("Genes not found: ", paste(missing_genes, collapse = ", "))
  # Spot positions in micrometres and grid pitch (calibration check)
  pos <- get_spot_positions_um(sp_obj, spot_size_um)
  barcodes <- pos$barcodes; xy <- pos$xy; pitch_um <- pos$pitch_um
  n_missing_spots <- length(setdiff(barcodes, colnames(sp_obj[[assay]])))
  if (n_missing_spots > 0) stop(sprintf("Assay '%s' covers only %d of %d spots (sketched subset?) - use a full-dataset assay such as Spatial or SCT",
                                        assay, length(barcodes) - n_missing_spots, length(barcodes)))
  expr <- tryCatch(SeuratObject::LayerData(sp_obj, assay = assay, layer = "data")[gene_list, barcodes, drop = FALSE],
                   error = function(e) stop(sprintf("Assay '%s' has no usable 'data' layer (normalize it first): %s", assay, conditionMessage(e))))
  cl_info <- get_cluster_info(sp_obj)
  spot_clusters <- factor(cl_info$spots[barcodes], levels = cl_info$names)
  # Integer indices (into barcodes) of the spots expressing each gene
  positive <- lapply(setNames(gene_list, gene_list), function(g) which(as.numeric(expr[g, ]) > 0))
  radius <- max_distance_um + 1e-6
  # Upper bound on grid points inside the search disc, so radius searches are not truncated
  k_max <- if (is.finite(pitch_um) && pitch_um > 0) ceiling(4 * (max_distance_um / pitch_um + 1)^2) + 1 else length(barcodes)
  # All (reference, target) pairs within the radius as integer spot indices
  find_pairs <- function(ref, tgt) {
    if (length(ref) == 0 || length(tgt) == 0) return(NULL)
    k <- min(length(tgt), k_max); truncated <- FALSE
    # Chunk size adapted to k so each nn2 result stays around 2e7 cells (bounded memory on Visium HD)
    chunk_size <- max(1000L, as.integer(floor(2e7 / k)))
    chunks <- split(seq_along(ref), ceiling(seq_along(ref) / chunk_size))
    res <- lapply(chunks, function(ix) {
      nn <- RANN::nn2(xy[tgt, , drop = FALSE], xy[ref[ix], , drop = FALSE], k = k, searchtype = "radius", radius = radius)
      hit <- nn$nn.idx > 0
      if (k < length(tgt) && any(rowSums(hit) == k)) truncated <<- TRUE
      return(data.frame(ref_spot = ref[ix][row(hit)[hit]], target_spot = tgt[nn$nn.idx[hit]], distance_um = nn$nn.dists[hit]))
    })
    if (truncated) warning("Neighbour search hit its k bound: pair counts may be underestimated")
    pairs <- do.call(rbind, res)
    if (!include_same_spot) pairs <- pairs[pairs$ref_spot != pairs$target_spot, , drop = FALSE]
    if (nrow(pairs) == 0) return(NULL)
    pairs$ref_cluster <- spot_clusters[pairs$ref_spot]; pairs$target_cluster <- spot_clusters[pairs$target_spot]
    pairs$same_cluster <- as.integer(pairs$ref_cluster) == as.integer(pairs$target_cluster)
    return(pairs)
  }
  # TRUE for each reference spot having at least one target spot within the radius
  has_neighbor <- function(ref, tgt) {
    if (length(ref) == 0) return(logical(0))
    if (length(tgt) == 0) return(rep(FALSE, length(ref)))
    # Without same-spot counting, the second neighbour is needed when the nearest one is the reference spot itself
    k <- min(length(tgt), if (include_same_spot) 1 else 2)
    nn <- RANN::nn2(xy[tgt, , drop = FALSE], xy[ref, , drop = FALSE], k = k, searchtype = "radius", radius = radius)
    ok <- nn$nn.idx > 0
    if (!include_same_spot) ok <- ok & matrix(tgt[pmax(nn$nn.idx, 1)], nrow = nrow(ok)) != ref
    return(rowSums(ok) > 0)
  }
  params <- list(max_distance_um = max_distance_um, spot_size_um = spot_size_um, um_per_px = pos$um_per_px,
                 pitch_um = round(pitch_um, 1), assay = assay, include_same_spot = include_same_spot, n_permutations = n_permutations)
  safe_pct <- function(num, den) if (den > 0) round(100 * num / den, 1) else NA_real_
  # Two genes: directional pair statistics
  if (length(gene_list) == 2) {
    g1 <- gene_list[1]; g2 <- gene_list[2]; p1 <- positive[[g1]]; p2 <- positive[[g2]]
    message(sprintf("=== Co-localization: %s vs %s | %g um | %.3f um/px | pitch %.1f um ===", g1, g2, max_distance_um, pos$um_per_px, pitch_um))
    run_21 <- isTRUE(bidirectional)
    pairs_12 <- find_pairs(p1, p2)
    pairs_21 <- if (run_21) find_pairs(p2, p1) else NULL
    near_12 <- if (is.null(pairs_12)) integer(0) else unique(pairs_12$ref_spot)
    near_21 <- if (is.null(pairs_21)) integer(0) else unique(pairs_21$ref_spot)
    # Permutation test: target-gene positive spots relocated at random among all spots, reference spots fixed
    permutation <- NULL
    if (n_permutations > 0 && length(p1) > 0 && length(p2) > 0) {
      set.seed(seed)
      n_spots <- length(barcodes)
      null_rates <- function(ref, n_tgt) vapply(seq_len(n_permutations), function(b) mean(has_neighbor(ref, sample.int(n_spots, n_tgt))), numeric(1))
      summarize_null <- function(observed, null) {
        expected <- mean(null)
        return(list(observed_pct = round(100 * observed, 1), expected_pct = round(100 * expected, 1),
                    enrichment = if (expected > 0) round(observed / expected, 2) else NA_real_,
                    p_value = (1 + sum(null >= observed)) / (n_permutations + 1)))
      }
      permutation <- list(n_permutations = n_permutations, seed = seed,
                          gene1_to_gene2 = summarize_null(length(near_12) / length(p1), null_rates(p1, length(p2))))
      if (run_21) permutation$gene2_to_gene1 <- summarize_null(length(near_21) / length(p2), null_rates(p2, length(p1)))
      message(sprintf("Permutations (%d): %s -> %s p = %.4g", n_permutations, g1, g2, permutation$gene1_to_gene2$p_value))
    }
    idx <- seq_along(barcodes)
    spatial_df <- data.frame(spot_id = barcodes, x = xy[, "x"], y = xy[, "y"], gene1_expr = as.numeric(expr[g1, ]), gene2_expr = as.numeric(expr[g2, ]),
                             gene1_positive = idx %in% p1, gene2_positive = idx %in% p2, gene1_has_neighbor = idx %in% near_12,
                             gene2_has_neighbor = idx %in% near_21, cluster = as.character(spot_clusters), stringsAsFactors = FALSE)
    # Pair tables exported with barcodes; NULL when too large to keep in memory or write to Excel
    export_pairs <- function(pairs) {
      if (is.null(pairs)) return(data.frame())
      if (nrow(pairs) > max_pairs_returned) return(NULL)
      return(data.frame(ref_spot = barcodes[pairs$ref_spot], ref_cluster = as.character(pairs$ref_cluster), target_spot = barcodes[pairs$target_spot],
                        target_cluster = as.character(pairs$target_cluster), same_cluster = pairs$same_cluster, distance_um = round(pairs$distance_um, 1),
                        stringsAsFactors = FALSE))
    }
    return(list(
      analysis_type = "bidirectional",
      parameters = c(list(gene1 = g1, gene2 = g2, bidirectional = run_21), params),
      gene1_total_spots = length(p1), gene2_total_spots = length(p2),
      gene1_near_gene2 = length(near_12), gene1_coloc_rate = safe_pct(length(near_12), length(p1)),
      gene2_near_gene1 = if (run_21) length(near_21) else NA, gene2_coloc_rate = if (run_21) safe_pct(length(near_21), length(p2)) else NA,
      permutation = permutation,
      cluster_pairs_gene1_to_gene2 = export_pairs(pairs_12), cluster_pairs_gene2_to_gene1 = export_pairs(pairs_21),
      cluster_summary = calculate_cluster_summary(p1, p2, pairs_12, pairs_21, spot_clusters),
      inter_cluster_gene1_to_gene2 = calculate_inter_cluster_stats(pairs_12, p1, spot_clusters),
      inter_cluster_gene2_to_gene1 = if (run_21) calculate_inter_cluster_stats(pairs_21, p2, spot_clusters) else data.frame(),
      spatial_data = spatial_df))
  }
  # Three or more genes: simultaneous co-localization (no permutation test)
  message(sprintf("=== Multi-gene co-localization: %s | %g um | pitch %.1f um ===", paste(gene_list, collapse = " + "), max_distance_um, pitch_um))
  coloc_counts <- sapply(gene_list, function(g) {
    ok <- Reduce(`&`, lapply(setdiff(gene_list, g), function(h) has_neighbor(positive[[g]], positive[[h]])), rep(TRUE, length(positive[[g]])))
    return(sum(ok))
  })
  totals <- sapply(gene_list, function(g) length(positive[[g]]))
  return(list(analysis_type = "simultaneous", parameters = c(list(gene_list = gene_list), params), total_spots_per_gene = totals,
              colocalized_spots_per_gene = coloc_counts, colocalization_rates = ifelse(totals > 0, round(100 * coloc_counts / totals, 1), 0)))
}

# Co-localization category of each spot of a two-gene analysis, for display in the tissue viewer.
# spatial_df and params come from analyze_gene_colocalization(). Levels are in drawing order (later levels on top):
# Neither, gene 1 only, gene 2 only, gene 1 next to gene 2, gene 2 next to gene 1, both genes in the same spot.
# partners_only = TRUE shows positive spots without a partner within the distance as Neither.
# Without the bidirectional analysis, gene 2 spots are never marked "next to gene 1".
colocalization_categories <- function(spatial_df, params, partners_only = FALSE) {
  g1 <- params$gene1; g2 <- params$gene2
  levels <- c("Neither", paste(g1, "only"), paste(g2, "only"), paste(g1, "next to", g2), paste(g2, "next to", g1), "Both in the same spot")
  category <- rep(levels[1], nrow(spatial_df))
  if (!partners_only) {
    category[spatial_df$gene1_positive] <- levels[2]
    category[spatial_df$gene2_positive] <- levels[3]
  }
  category[spatial_df$gene1_positive & spatial_df$gene1_has_neighbor] <- levels[4]
  category[spatial_df$gene2_positive & spatial_df$gene2_has_neighbor] <- levels[5]
  category[spatial_df$gene1_positive & spatial_df$gene2_positive] <- levels[6]
  colors <- stats::setNames(c("#D9D9D9", "#9ECAE1", "#FCBBA1", "#08519C", "#A50F15", "#6A3D9A"), levels)
  return(list(category = stats::setNames(category, spatial_df$spot_id), levels = levels, colors = colors))
}

# Static co-localization map: tissue image with the spots of a two-gene analysis coloured by category
# (see colocalization_categories). spots: data.frame with barcode, x, y in full-resolution pixels (y = image row).
# Spots expressing neither gene are not drawn, so the tissue stays visible; without a tissue image they are drawn in light
# grey to show the outline. point_size is a fixed ggplot size: on a whole-slide view Visium HD bins would be sub-pixel.
# The legend gives the number of spots of each category.
build_colocalization_map <- function(spots, categories, crop = NULL, view, point_size = 0.8, title = NULL, subtitle = NULL) {
  value <- unname(categories$category[spots$barcode])
  value[is.na(value)] <- categories$levels[1]
  spots$value <- factor(value, levels = categories$levels)
  spots$y_plot <- -spots$y
  counts <- table(spots$value)
  shown <- spots[as.integer(spots$value) > 1, , drop = FALSE]
  # Later levels drawn on top: spots with a partner stay visible over the single-gene spots
  shown <- shown[order(shown$value), , drop = FALSE]
  p <- ggplot(data.frame(x = view[1:2], y_plot = -view[3:4]), aes(x = x, y = y_plot)) + geom_blank()
  if (!is.null(crop)) {
    p <- p + annotation_raster(grDevices::as.raster(crop$img), xmin = crop$extent[1], xmax = crop$extent[2],
                               ymin = -crop$extent[4], ymax = -crop$extent[3], interpolate = TRUE)
  } else {
    p <- p + geom_point(data = spots[as.integer(spots$value) == 1, , drop = FALSE], color = "#E5E5E5", size = 0.1, shape = 16)
  }
  present <- categories$levels[-1][counts[categories$levels[-1]] > 0]
  p <- p + geom_point(data = shown, aes(color = value), size = point_size, shape = 16, alpha = 0.9) +
    scale_color_manual(values = categories$colors, breaks = present, labels = function(b) sprintf("%s (%s)", b, formatC(as.integer(counts[b]), big.mark = ",", format = "d")),
                       name = NULL, drop = TRUE) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1)) +
    coord_fixed(xlim = view[1:2], ylim = c(-view[4], -view[3]), expand = FALSE) +
    labs(title = title, subtitle = subtitle) + theme_void() +
    theme(legend.position = "right", legend.text = element_text(size = 11), plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5, colour = "grey30"), plot.background = element_rect(fill = "white", colour = NA))
  return(p)
}

# Spots with a partner of the other gene within the distance, per cluster and direction (two-gene analysis).
# clusters: cluster label of every spot, named by barcode. Returns cluster, direction, n_spots (pairs only, empty rows dropped).
colocalization_by_cluster <- function(spatial_df, params, clusters) {
  cl <- as.character(clusters[spatial_df$spot_id])
  cl[is.na(cl)] <- "No cluster"
  rows <- list(data.frame(cluster = cl[spatial_df$gene1_positive & spatial_df$gene1_has_neighbor], direction = paste(params$gene1, "next to", params$gene2)),
               data.frame(cluster = cl[spatial_df$gene2_positive & spatial_df$gene2_has_neighbor], direction = paste(params$gene2, "next to", params$gene1)))
  d <- do.call(rbind, rows)
  if (nrow(d) == 0) return(data.frame(cluster = character(0), direction = character(0), n_spots = integer(0)))
  out <- as.data.frame(table(cluster = d$cluster, direction = d$direction), stringsAsFactors = FALSE)
  names(out)[3] <- "n_spots"
  out <- out[out$n_spots > 0, , drop = FALSE]
  # Cluster order of the object when the labels are its levels
  levels_cl <- if (is.factor(clusters)) intersect(levels(clusters), out$cluster) else sort(unique(out$cluster))
  out$cluster <- factor(out$cluster, levels = unique(c(levels_cl, out$cluster)))
  return(out)
}

############################## Spatially variable genes (Moran's I) ##############################

# Symmetric binary k-nearest-neighbour graph between spots (sparse spots x spots matrix, zero diagonal).
# Positions can be in pixels or micrometres: only the neighbour order and relative distances matter.
# Neighbours farther than max_pitch times the grid pitch (median nearest-neighbour distance) are dropped, so spots of
# separate tissue pieces, or across holes and tears, are never linked; a spot may then have fewer than k neighbours.
# groups (optional, one value per spot, e.g. the tissue section): spots of different groups are never linked.
build_spatial_knn_weights <- function(xy, k = 6, max_pitch = 3, groups = NULL) {
  n <- nrow(xy); k <- as.integer(k)
  if (is.na(k) || k < 1) stop("The number of neighbours must be a positive integer")
  if (n < k + 2) stop(sprintf("At least %d spots are needed for a %d-neighbour graph", k + 2, k))
  nn <- RANN::nn2(xy, xy, k = k + 1)
  pitch <- stats::median(nn$nn.dists[, 2])
  # Column-major nn.idx: entry (r, c) belongs to spot r
  i <- rep(seq_len(n), k + 1); j <- as.vector(nn$nn.idx)
  keep <- as.vector(nn$nn.dists) <= max_pitch * pitch
  if (!is.null(groups)) {
    groups <- as.character(groups); groups[is.na(groups)] <- "<none>"
    keep <- keep & groups[i] == groups[j]
  }
  a <- Matrix::sparseMatrix(i = i[keep], j = j[keep], x = 1, dims = c(n, n))
  w <- methods::as(((a + Matrix::t(a)) > 0) * 1, "CsparseMatrix")
  Matrix::diag(w) <- 0
  return(Matrix::drop0(w))
}

# Moran's I of each gene (rows of expr, genes x spots) on the spot graph w, with its z-score under the normality
# assumption (Cliff & Ord): E[I] = -1/(n-1). One-sided p-value for positive spatial autocorrelation.
# Genes are processed in chunks so the dense block stays small on Visium HD.
compute_morans_i <- function(expr, w, chunk_size = 100) {
  n <- ncol(expr)
  s0 <- sum(w); s1 <- 0.5 * sum((w + Matrix::t(w))^2); s2 <- sum((Matrix::rowSums(w) + Matrix::colSums(w))^2)
  expected <- -1 / (n - 1)
  variance <- (n^2 * s1 - n * s2 + 3 * s0^2) / ((n^2 - 1) * s0^2) - expected^2
  chunks <- split(seq_len(nrow(expr)), ceiling(seq_len(nrow(expr)) / chunk_size))
  morans <- unlist(lapply(chunks, function(ix) {
    z <- as.matrix(Matrix::t(expr[ix, , drop = FALSE]))
    z <- sweep(z, 2, colMeans(z))
    return((n / s0) * colSums(z * as.matrix(w %*% z)) / colSums(z^2))
  }), use.names = FALSE)
  z_score <- (morans - expected) / sqrt(variance)
  return(data.frame(gene = rownames(expr), morans_I = morans, z_score = z_score, p_val = stats::pnorm(z_score, lower.tail = FALSE), stringsAsFactors = FALSE))
}

# Spatially variable genes: Moran's I on the normalized 'data' layer over a k-nearest-neighbour spot graph,
# Benjamini-Hochberg correction over the genes tested. Deterministic (no random step).
# Tested genes: the variable features of the assay ("variable") or the n_features most expressed genes ("expressed");
# genes detected in fewer than min_pct of the spots are skipped. The assay must cover every spot (no sketched subset).
# Each gene is also linked to the cluster where its mean expression is highest: top_cluster_share is that mean divided by
# the sum of the per-cluster means (independent of cluster sizes; 1 / number of clusters = evenly spread, 1 = single cluster).
# P-values are returned as raw numerics (formatting is done at display time).
find_spatially_variable_genes <- function(sp_obj, assay = "Spatial", feature_source = "variable", n_features = 2000, k = 6, min_pct = 0.01, log_function = message) {
  if (!requireNamespace("RANN", quietly = TRUE)) stop("Package 'RANN' is required for the neighbour graph")
  if (!assay %in% names(sp_obj@assays)) stop("Assay not found: ", assay)
  n_features <- as.integer(n_features)
  if (is.na(n_features) || n_features < 1) stop("The number of genes to test must be a positive integer")
  coords <- get_spot_coordinates_spatial(sp_obj)
  barcodes <- coords$barcode
  n_missing <- length(setdiff(barcodes, colnames(sp_obj[[assay]])))
  if (n_missing > 0) stop(sprintf("Assay '%s' covers only %d of %d spots (sketched subset?) - use a full-dataset assay such as Spatial or SCT", assay, length(barcodes) - n_missing, length(barcodes)))
  data_layer <- tryCatch(SeuratObject::LayerData(sp_obj, assay = assay, layer = "data"), error = function(e) NULL)
  if (is.null(data_layer)) stop(sprintf("Assay '%s' has no 'data' layer: normalize the data first", assay))
  expr_all <- data_layer[, barcodes, drop = FALSE]
  pct <- Matrix::rowMeans(expr_all > 0)
  candidates <- names(pct)[pct >= min_pct]
  if (length(candidates) == 0) stop("No gene is detected in the minimum fraction of spots")
  features <- if (identical(feature_source, "variable")) intersect(SeuratObject::VariableFeatures(sp_obj, assay = assay), candidates) else character(0)
  if (length(features) == 0) {
    if (identical(feature_source, "variable")) log_function("No variable features for this assay: the most expressed genes are tested instead")
    mean_expr <- Matrix::rowMeans(expr_all[candidates, , drop = FALSE])
    features <- names(sort(mean_expr, decreasing = TRUE))
  }
  features <- utils::head(features, n_features)
  log_function(sprintf("Moran's I: %d genes, %d spots, %d-neighbour graph", length(features), length(barcodes), as.integer(k)))
  w <- build_spatial_knn_weights(cbind(coords$imagecol, coords$imagerow), k = k, groups = spot_sections(sp_obj, coords$barcode))
  res <- compute_morans_i(expr_all[features, , drop = FALSE], w)
  # Constant genes give 0/0: they cannot be tested
  res <- res[is.finite(res$morans_I), , drop = FALSE]
  res$p_val_adj <- stats::p.adjust(res$p_val, method = "BH")
  res$pct_spots <- round(as.numeric(pct[res$gene]), 3)
  res$mean_expr <- as.numeric(Matrix::rowMeans(expr_all[res$gene, , drop = FALSE]))
  # Cluster with the highest mean expression, through a sparse spots x clusters indicator (no dense genes x spots matrix)
  cl <- get_cluster_info(sp_obj)
  spot_clusters <- factor(cl$spots[barcodes], levels = cl$names)
  if (nlevels(spot_clusters) > 1 && !anyNA(spot_clusters)) {
    indicator <- Matrix::sparseMatrix(i = seq_along(spot_clusters), j = as.integer(spot_clusters), x = 1, dims = c(length(spot_clusters), nlevels(spot_clusters)))
    cluster_means <- as.matrix(expr_all[res$gene, , drop = FALSE] %*% indicator) %*% diag(1 / as.numeric(table(spot_clusters)), nrow = nlevels(spot_clusters))
    top <- max.col(cluster_means, ties.method = "first")
    res$top_cluster <- levels(spot_clusters)[top]
    res$top_cluster_share <- round(cluster_means[cbind(seq_len(nrow(res)), top)] / rowSums(cluster_means), 3)
  } else {
    res$top_cluster <- NA_character_
    res$top_cluster_share <- NA_real_
  }
  res <- res[order(-res$morans_I), , drop = FALSE]
  rownames(res) <- NULL
  attr(res, "params") <- list(method = "Moran's I (kNN graph, normality assumption, BH)", assay = assay, feature_source = feature_source,
                              n_tested = nrow(res), k = as.integer(k), min_pct = min_pct, n_spots = length(barcodes),
                              cluster_column = cl$column_used, n_clusters = nlevels(spot_clusters))
  return(res)
}
