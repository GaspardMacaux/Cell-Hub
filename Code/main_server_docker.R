############################## Server ##############################


############################## Source Files ##############################
# Importer les données
options(shiny.maxRequestSize = 300000*1024^2)  # Augmente à 30GB par fichier
options(future.globals.maxSize = 6000000 * 1024^2)  # 100GB
Sys.setenv('R_MAX_VSIZE'=20000000000000)


# Source server files:
source("single_dataset_server_docker.R")
source("multiple_datasets_server_docker.R")
source("cellchat_server_docker.R")
source("spatial_transcriptomic_server_docker.R")
source("trajectory_server_docker.R")

#Source function:
source("functions/data_loading_functions_docker.R")
source("functions/workspace_management_docker.R")
source("functions/download_functions_docker.R")
source("functions/venn_diagram_functions_docker.R")
source("functions/ui_update_functions_docker.R")
source("functions/cells_genes_expressions_docker.R")
source("functions/integration_functions_docker.R")
source("functions/cellchat_helpers_docker.R")
source("functions/dimensional_reduction_functions_docker.R")
source("functions/heatmap_scatterplot_functions_docker.R") 
source("functions/volcano_plot_functions_docker.R")
source("functions/genes_visualization_functions_docker.R")
source("functions/trajectory_functions_docker.R")
source("functions/spatial_functions_docker.R")
source("functions/annotation_functions_docker.R")
source("functions/enrichment_functions_docker.R")
source("functions/error_handling_functions_docker.R")
source("functions/example_data_functions_docker.R")
source("functions/cellhubdb_filter_functions_docker.R")


#Source UI:
source("main_ui_docker.R")

server <- function(input, output, session) {
  single_dataset_server(input, output, session)
  multiple_datasets_server(input, output, session)
  trajectory_server(input, output, session)
  cellchat_server(input, output, session)
  spatial_transcriptomic_server (input, output, session)
}


# Exécuter l'application
options(shiny.launch.browser = TRUE)

shinyApp(ui = ui, server = server)
