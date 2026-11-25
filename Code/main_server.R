############################## Server ##############################

############################## Library ##############################

# Interface
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(shinyFiles)
library(fontawesome)
library(jsonlite)
library(tools)
library(R.utils)

# Graphics and tables
library(plotly)
library(tidyverse)   
library(patchwork)
library(Matrix)
library(ggrepel)
library(DT)
library(cowplot)
library(data.table)
library(httr)
library(rhdf5)
library(hdf5r)
library(colourpicker)
library(uwot)
library(circlize)
library(callr)
library(viridis)

# Single cell
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(monocle3)
library(SingleCellExperiment)
library(harmony)
library(clustree)
library(EnsDb.Hsapiens.v86)
library(biomaRt)
library(DoubletFinder)
library(CellChat)
library(VennDiagram)
library(grid)
library(openxlsx)
library(glmGamPoi)

# Image handling
library(tiff)
library(png)
library(jpeg)
library(readODS)
library(schard)


#New library
library(nanoparquet)
library(SingleR)
library(celldex)

############################## Source Files ##############################
# Increase memory limits for large dataset uploads
options(shiny.maxRequestSize = 100000*1024^2)  # ~100 GB max per file upload
options(future.globals.maxSize = 6000000 * 1024^2)  # ~6 TB max for future globals
Sys.setenv('R_MAX_VSIZE'=20000000000000)


# Source server files:
source("single_dataset_server.R")
source("multiple_datasets_server.R")
source("cellchat_server.R")
source("spatial_transcriptomic_server.R")
source("trajectory_server.R")
#Source function:
source("functions/data_loading_functions.R")
source("functions/workspace_management.R")
source("functions/download_functions.R")
source("functions/venn_diagram_functions.R")
source("functions/ui_update_functions.R")
source("functions/cells_genes_expressions.R")
source("functions/integration_functions.R")
source("functions/cellchat_helpers.R")
source("functions/dimensional_reduction_functions.R")
source("functions/heatmap_scatterplot_functions.R") 
source("functions/volcano_plot_functions.R")
source("functions/genes_visualization_functions.R")
source("functions/trajectory_functions.R")
source("functions/spatial_functions.R")
source("functions/annotation_functions.R")
source("functions/enrichment_functions.R")
source("functions/error_handling_functions.R")
source("functions/example_data_functions.R")
source("functions/cellhubdb_filter_functions.R")

#Source UI:
source("main_ui.R")


server <- function(input, output, session) {
single_dataset_server(input, output, session)
multiple_datasets_server(input, output, session)
trajectory_server(input, output, session)
cellchat_server(input, output, session)
spatial_transcriptomic_server(input, output, session)
}

# Launch the Shiny application
options(shiny.launch.browser = TRUE)

shinyApp(ui = ui, server = server)
