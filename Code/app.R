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


source("main_server_docker.R")



# Exécuter l'application
options(shiny.launch.browser = TRUE)
# Lancer l'application avec ui et server chargés
shinyApp(ui = ui, server = server)
