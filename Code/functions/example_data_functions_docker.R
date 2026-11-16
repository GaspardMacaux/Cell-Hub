# example_data_functions.R
# Demonstration data linked from the loading pages. Nothing is shipped in the image: the files are deposited on Zenodo
# (public 10x Genomics datasets, CC BY 4.0, redistributed unchanged) and downloaded by the user. UI builder only (no reactive code).

############################## Demonstration data ##############################

# DOI of the Zenodo record holding the demonstration files (reserved on Zenodo before publication)
demo_data_doi <- "10.5281/zenodo.24198709"
demo_data_url <- paste0("https://doi.org/", demo_data_doi)

# Files of the record used by each module
demo_files <- list(
  pbmc_v3 = "pbmc_1k_v3_filtered_feature_bc_matrix.h5",
  pbmc_v2 = "pbmc_1k_v2_filtered_feature_bc_matrix.h5",
  visium = "visium_mouse_brain.zip"
)

############################## Example data box ##############################

# "Try with demonstration data" box; module is one of "single", "multiple", "cellchat", "trajectory", "spatial";
# collapsed = TRUE shows only the title until the user clicks it
example_data_box <- function(module, collapsed = FALSE) {
  record <- tags$a("Cell-Hub demonstration data (Zenodo)", href = demo_data_url, target = "_blank", rel = "noopener")
  file_tag <- function(key) tags$code(demo_files[[key]])
  steps <- switch(module,
                  single = tagList(
                    tags$li("Download ", file_tag("pbmc_v3"), " from ", record, " and load it above as raw data (species: human)."),
                    tags$li("1,000 human blood cells: QC, clustering, markers and annotation (SingleR) run in a few minutes.")),
                  multiple = tagList(
                    tags$li("Download ", file_tag("pbmc_v2"), " and ", file_tag("pbmc_v3"), " from ", record, " and load them as two datasets (type scRNA-seq, species human)."),
                    tags$li("The same donor sequenced with two chemistries forms a technical batch: compare CCA, Harmony, fastMNN and simple merge.")),
                  cellchat = tagList(
                    tags$li("Process ", file_tag("pbmc_v3"), " (", record, ") in Single Dataset Analysis up to annotation, then use ", tags$em("Use current object"), " (species: human).")),
                  trajectory = tagList(
                    tags$li("Process ", file_tag("pbmc_v3"), " (", record, ") in Single Dataset Analysis up to clustering, then import the current object.")),
                  spatial = tagList(
                    tags$li("Download ", file_tag("visium"), " from ", record, " and load the ZIP as it is (species: mouse)."),
                    tags$li("Visium section of adult mouse brain: clusters, tissue viewer, spatially variable genes and co-localization.")),
                  stop("Unknown module: ", module))
  tags$div(
    style = "margin-top: 20px; background-color: #f1f8e9; padding: 15px 20px; border-radius: 8px; border-left: 4px solid #7cb342;",
    tags$div(style = "cursor: pointer;", onclick = "$(this).next().slideToggle();",
             tags$h4(icon("flask"), " No data at hand? Try the demonstration data", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"),
                     style = "margin: 0; color: #33691e; font-weight: bold;")),
    tags$div(style = paste0("margin-top: 12px;", if (collapsed) " display: none;" else ""),
             tags$ul(steps, style = "margin-bottom: 6px;"),
             tags$p("Public 10x Genomics datasets (CC BY 4.0) redistributed unchanged on Zenodo; downloaded by the user, not shipped with Cell-Hub.",
                    style = "font-size: 11px; color: #6c757d; margin: 0;"))
  )
}