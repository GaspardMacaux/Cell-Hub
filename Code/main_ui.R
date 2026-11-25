############################## User Interface ##############################

ui <- dashboardPage(
  
  ############################## HTML et CSS ##############################
  
  dashboardHeader(title = tags$img(src = "header.png",   height = "53px", width = "218px")),
  dashboardSidebar(
    tags$head(
      tags$style(HTML("
                       
                        /* Aligner le contenu de la sidebar avec le contenu principal */
                        .sidebar .sidebar-menu {
                          margin-top: 5px;
                        }

                        /* Augmenter la largeur de la bulle d'information */
                        .popover {
                          width: 300px;
                        }

                        /* Style pour ajouter l'image au header */
                        .main-header {
                          background-image: url('www/header.png');
                          background-repeat: no-repeat;
                          background-position: center center;
                        }

                        .info-box {
                          width: 100%;
                        }
/* Style pour le plot carré réactif */
                .reactive-square-plot {
                    width: 100%;
                    height: auto;
                    aspect-ratio: 1 / 1;
                }
                      ")),
      tags$script(src = "my_script.js")
    ),
    useShinyjs(),
    
    ############################## Sidebar Menu ##############################
    sidebarMenu(
      menuItem("Introduction", tabName = "introduction", icon = icon("mug-hot")),
      menuItem(
        "Single Dataset Analysis", icon = icon("virus-covid"),
        menuItem("Load dataset", tabName = "load_dataset", icon = icon("download")),
        menuItem("Data cleanup & Variable features", tabName = "qc", icon = icon("bug")),
        menuItem("Dimensional reduction", tabName = "dimensional_reduction", icon = icon("home")),
        menuItem("Clustering", tabName = "clustering", icon = icon("brain")),
        menuItem("Doublet Detection", tabName = "doublet_detection",icon = icon("dna")),
        menuItem("Plot Gene expressions", tabName = "plots_genes_expressions", icon = icon("flask")),
        menuItem("Assigning cell types identity", tabName = "assigning_cell_type_identity", icon = icon("id-card")),
        menuItem("Clusters comparison", tabName = "cluster_comparison", icon = icon("earth-americas")),
        menuItem("Pathway Enrichment", tabName = "enrichment_single", icon = icon("diagram-project")),
        menuItem("DEG Visualizations", tabName = "deg_visualizations", icon = icon("volcano")),
        menuItem("Exclusive Biomarkers", tabName = "exclusive_biomarkers", icon = icon("star")),  
        menuItem("Heatmaps", tabName = "heatmaps_dualexpression", icon = icon("thermometer-half")),
        menuItem("Co-expression Analysis", tabName = "coexpression_single", icon = icon("project-diagram")),
        menuItem("Subset", tabName = "subset", icon = icon("project-diagram")),
        menuItem("Save & Session", tabName = "session_single", icon = icon("floppy-disk"))
      ),
      menuItem(
        "Multiple Datasets Analysis", icon = icon("viruses"),
        menuItem("Load datasets", tabName = "load_datasets_merge", icon = icon("download")),
        menuItem("Metadata Management", tabName = "metadata_management_merge", icon = icon("tags")),  
        menuItem("Clustering", tabName = "clustering_merge", icon = icon("brain")),
        menuItem("Plot Gene expressions", tabName = "plots_genes_expressions_merge", icon = icon("flask")),
        menuItem("Assigning cell types identity", tabName = "assigning_cell_type_identity_merge", icon = icon("id-card")),
        menuItem("Clusters comparison", tabName = "cluster_comparison_merge", icon = icon("earth-americas")),
        menuItem("Pathway Enrichment", tabName = "enrichment_merge", icon = icon("diagram-project")),
        menuItem("Exclusive Biomarkers", tabName = "exclusive_biomarkers_merge", icon = icon("star")), 
        menuItem("Heatmaps", tabName = "heatmaps_dualexpression_merge", icon = icon("thermometer-half")),
        menuItem("Co-expression Analysis", tabName = "coexpression_merge", icon = icon("project-diagram")),
        menuItem("DEG Visualizations", tabName = "deg_visualizations_merge", icon = icon("volcano")),
        menuItem("Subset", tabName = "subset_merge", icon = icon("project-diagram")),
        menuItem("Save & Session", tabName = "session_merge", icon = icon("floppy-disk"))
      ),
      menuItem(
        "Cell Chat", icon = icon("envelope"),
        menuItem("Load dataset", tabName = "load_data_cellchat", icon = icon("download")),
        menuItem("Ligand-Receptor", tabName = "ligand_receptor_cellchat", icon = icon("link")),
        menuItem("Circle Plot", tabName = "circle_plot_cellchat", icon = icon("circle"))
        
      ),
      
      
      menuItem(
        "Monocle", icon = icon("magnifying-glass"),
        menuItem("Load dataset", tabName = "load_trajectory", icon = icon("download")),
        menuItem("Trajectory construction", tabName = "trajectory_construction", icon = icon("route")),
        menuItem("Differentialy expressed genes", tabName = "differentialy_expressed_genes_trajectory", icon = icon("route")),
        menuItem("Gene expression visualisation", tabName = "genes_visualization_trajectory", icon = icon("route"))
        
      ),
      menuItem(
        "Transcriptomic Spatial", tabName = "transcriptomic_spatial", icon = icon("map"),
        menuItem("Load Spatial dataset", tabName = "load_spatial_dataset", icon = icon("download")),
        menuItem("Tissue Sections", tabName = "spatial_sections", icon = icon("object-ungroup")),
        menuItem("Normalisation Step", tabName = "spatial_normalisation", icon = icon("chart-line")),
        menuItem("Clustering", tabName = "spatial_clustering", icon = icon("circle")),
        menuItem("Interactive Tissue Viewer", tabName = "spatial_interactive_visualization", icon = icon("magnifying-glass")),
        menuItem("Gene expression visualisation", tabName = "spatial_gene_expression_visualisation", icon = icon("database")),
        menuItem("Cluster Annotation", tabName = "spatial_cluster_annotation", icon = icon("tags")),
        menuItem("Co-localization", tabName = "spatial_colocalization", icon = icon("circle-dot")),
        menuItem("Spatially Variable Genes", tabName = "spatial_svg", icon = icon("braille")),
        menuItem("Marker Analysis", tabName = "spatial_marker_analysis", icon = icon("dna"))
        
      ),
      menuItem("Acknowledgement & Licence", tabName = "acknowledgement", icon = icon("heart"))
    )
  ),
  
  ############################## Single/Introduction ##############################
  dashboardBody(
    includeCSS("www/custom_styles.css"),
    
    tabItems(
      tabItem(tabName = "introduction",
              tags$div(
                class = "header-violet",
                style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 50px 30px; margin: -15px -15px 30px -15px; text-align: center; color: white;",
                tags$h1(icon("flask"), " Cell-Hub", style = "font-size: 48px; font-weight: bold; margin-bottom: 15px;"),
                tags$h3("Single-Cell & Spatial Transcriptomics Analysis Platform", style = "font-weight: 300; opacity: 0.95;")
              ),
              fluidRow(
                column(3,
                       tags$div(
                         class = "info-card",
                         style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 25px; border-radius: 10px; color: white; text-align: center; box-shadow: 0 4px 6px rgba(0,0,0,0.1);",
                         icon("database", style = "font-size: 3em; margin-bottom: 10px;"),
                         tags$h3("Single-Cell", style = "margin: 10px 0; font-size: 18px;"),
                         tags$p("Seurat-based analysis", style = "opacity: 0.9; margin: 0;")
                       )
                ),
                column(3,
                       tags$div(
                         class = "info-card",
                         style = "background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); padding: 25px; border-radius: 10px; color: white; text-align: center; box-shadow: 0 4px 6px rgba(0,0,0,0.1);",
                         icon("project-diagram", style = "font-size: 3em; margin-bottom: 10px;"),
                         tags$h3("Cell-Cell", style = "margin: 10px 0; font-size: 18px;"),
                         tags$p("Communication analysis", style = "opacity: 0.9; margin: 0;")
                       )
                ),
                column(3,
                       tags$div(
                         class = "info-card",
                         style = "background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); padding: 25px; border-radius: 10px; color: white; text-align: center; box-shadow: 0 4px 6px rgba(0,0,0,0.1);",
                         icon("route", style = "font-size: 3em; margin-bottom: 10px;"),
                         tags$h3("Trajectory", style = "margin: 10px 0; font-size: 18px;"),
                         tags$p("Monocle pseudotime", style = "opacity: 0.9; margin: 0;")
                       )
                ),
                column(3,
                       tags$div(
                         class = "info-card",
                         style = "background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); padding: 25px; border-radius: 10px; color: white; text-align: center; box-shadow: 0 4px 6px rgba(0,0,0,0.1);",
                         icon("map-marked-alt", style = "font-size: 3em; margin-bottom: 10px;"),
                         tags$h3("Spatial", style = "margin: 10px 0; font-size: 18px;"),
                         tags$p("Visium & Visium HD", style = "opacity: 0.9; margin: 0;")
                       )
                )
              ),
              tags$div(
                style = "background: white; padding: 40px; margin-top: 30px; border-radius: 10px; box-shadow: 0 2px 8px rgba(0,0,0,0.1);",
                tags$h3(icon("sitemap"), " Analysis Workflows", style = "color: #667eea; margin-bottom: 40px; text-align: center;"),
                tags$h4("Single-Cell RNA-seq Analysis", style = "text-align: center; color: #2c3e50; margin-bottom: 25px; font-weight: bold;"),
                tags$div(
                  style = "display: grid; grid-template-columns: repeat(2, 1fr); gap: 30px; max-width: 1000px; margin: 0 auto 40px auto; position: relative;",
                  tags$div(
                    style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 25px; border-radius: 15px;",
                    tags$h4("Preprocessing", style = "text-align: center; margin-bottom: 15px; font-weight: bold; color: white;"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "barcodes.tsv / features.tsv / matrix.mtx"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "If multiple datasets, perform integration steps"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Preprocessing & QC"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Normalization steps"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; text-align: center; font-size: 14px;", "Find neighbors & Determine clusters")
                  ),
                  tags$div(
                    style = "background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); padding: 25px; border-radius: 15px;",
                    tags$h4("Characterization", style = "text-align: center; margin-bottom: 15px; font-weight: bold; color: white;"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Assigning cell types identity"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Differentially expressed genes"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Visualize gene expresions"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Compare genes expresions between clusters"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; text-align: center; font-size: 14px;", "Exclusive biomarker detection")
                  ),
                  tags$div(
                    style = "background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); padding: 25px; border-radius: 15px;",
                    tags$h4("Cell-Cell Communication", style = "text-align: center; margin-bottom: 15px; font-weight: bold; color: white;"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Load Seurat object"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Load CellHub Database"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Create Cellchat object"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; text-align: center; font-size: 14px;", "Calculate DE genes between cell types"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; text-align: center; font-size: 14px;", "Visualize ligands & receptors")
                  ),
                  tags$div(
                    style ="background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); padding: 25px; border-radius: 15px;",
                    tags$h4("Trajectory", style = "text-align: center; margin-bottom: 15px; font-weight: bold; color: white;"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Convert Seurat object to Monocle"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Determine trajectory & Pseudotime"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "DE genes along Pseudotime"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; text-align: center; font-size: 14px;", "Visualize genes expression on UMAP"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; text-align: center; font-size: 14px;", "Visualize genes expression along pseudotime")
                  ),
                  tags$div(
                    style = "position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); font-size: 40px; color: #666; opacity: 0.3;",
                    "↻"
                  )
                ),
                tags$hr(style = "margin: 50px 0; border: none; border-top: 2px solid #e0e0e0;"),
                tags$h4("Spatial Transcriptomics Analysis", style = "text-align: center; color: #2c3e50; margin-bottom: 25px; font-weight: bold;"),
                tags$div(
                  style = "max-width: 600px; margin: 0 auto;",
                  tags$div(
                    style = "background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); padding: 30px; border-radius: 15px;",
                    tags$h4(icon("map-marked-alt"), " Spatial Workflow", style = "text-align: center; margin-bottom: 20px; font-weight: bold; color: white;"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Load Visium / Visium HD data"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Spatial clustering & annotation"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; margin-bottom: 10px; text-align: center; font-size: 14px;", "Spatial marker identification"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; text-align: center; font-size: 14px;", "Spatial genes identification"),
                    tags$div(style = "text-align: center; font-size: 20px; margin: 8px 0; color: white;", "↓"),
                    tags$div(style = "background: rgba(255,255,255,0.25); color: white; padding: 12px; border-radius: 8px; text-align: center; font-size: 14px;", "Gene co-localization analysis")
                  )
                )
              ),
              tags$div(
                style = "text-align: center; margin-top: 50px; padding: 20px; color: #888;",
                tags$p(icon("code"), " Built with R, R studio & Shiny,"),
                tags$p(icon("heart"), " Developed for the single-cell community")
              )
      )
      ,
      
      
      ############################################ Single Dataset Analysis ###############################################
      
      ###############################################Load Dataset#############################################
      tabItem(
        tabName = "load_dataset",
        
        tags$div(
          class = "header-violet",
          tags$h2("Single Dataset Analysis", class = "header-title"),
          tags$p("Load and analyze your single-cell RNA sequencing data", class = "header-subtitle")
        ),
        tags$div(
          style = "background-color: #e8eaf6; padding: 20px; border-radius: 8px; border-left: 4px solid #667eea; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Single Dataset Analysis", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Single-cell RNA sequencing analyzes gene expression at the individual cell level, revealing cellular heterogeneity and identifying distinct cell populations."),
            tags$h5(icon("star"), " Supported Formats:", style = "color: #667eea; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("ZIP files:"), " Containing barcodes.tsv.gz, matrix.mtx.gz, and features.tsv.gz"),
              tags$li(tags$strong("H5 files:"), " 10X HDF5 format (filtered_feature_bc_matrix.h5)")
            ),
            tags$h5(icon("flask"), " Dataset Types:", style = "color: #667eea; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("snRNA-seq:"), " Single-nucleus RNA sequencing"),
              tags$li(tags$strong("Multiome:"), " Combined gene expression and chromatin accessibility")
            ),
            tags$div(
              style = "background-color: #fff9c4; padding: 12px; border-radius: 6px; border-left: 3px solid #fbc02d; margin-top: 15px;",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " For Multiome data, RNA expression will be used for downstream analysis", style = "margin: 0;")
            )
          )
        ),
        
        tags$div(
          style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; border: 2px solid #dee2e6; margin-bottom: 20px;",
          tags$h4(icon("cog"), " Dataset Parameters", style = "color: #495057; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
          fluidRow(
            column(6,
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #667eea;",
                     tags$label("Species:", style = "font-weight: bold; color: #495057;"),
                     selectInput("species_choice", NULL, choices = c("Mouse" = "mouse", "Human" = "human", "Rat" = "rat"), selected = "mouse"),
                     tags$p(icon("info-circle"), " Select organism for mitochondrial gene detection", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                   )
            ),
            column(6,
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #764ba2;",
                     tags$label("Dataset Type:", style = "font-weight: bold; color: #495057;"),
                     radioButtons("dataset_type", NULL, choices = list("snRNA-seq" = "snRNA", "Multiome" = "multiome"), selected = "snRNA"),
                     tags$p(icon("info-circle"), " Choose technology used for sequencing", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                   )
            )
          ),
          tags$div(
            style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; border: 2px solid #dee2e6; margin-bottom: 20px;",
            tags$h4(icon("memory"), " Memory Optimization", style = "color: #495057; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
            tags$div(
              style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #11998e;",
              radioButtons("optimize_memory_single", NULL,
                           choices = list(
                             "Safe: drop RNA scale.data only (recommended)" = "slim",
                             "Aggressive: also drop SCT / integrated assays" = "aggressive",
                             "Keep full object" = "full"
                           ),
                           selected = "slim", inline = FALSE),
              tags$p(icon("info-circle"),
                     " Safe: removes RNA scale.data only (recomputed before every PCA); nothing is lost. Aggressive: also removes SCT and integrated assays; their corrected values are lost and re-running PCA or clustering then uses uncorrected RNA data. Counts, normalized data, reductions and metadata are always kept.",
                     style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
            )
          ),
        ),
        
        tags$h4("Load Your Data", style = "color: #495057; margin-bottom: 15px; font-weight: bold;"),
        
        fluidRow(
          column(6,
                 tags$div(
                   class = "card-violet",
                   style = "min-height: 400px; display: flex; flex-direction: column;",
                   tags$div(
                     style = "text-align: center; margin-bottom: 20px;",
                     tags$div(
                       style = "background-color: rgba(255,255,255,0.2); width: 80px; height: 80px; border-radius: 50%; display: inline-flex; align-items: center; justify-content: center; margin-bottom: 15px;",
                       icon("file-archive", style = "font-size: 36px;")
                     ),
                     tags$h4("Load Raw 10X Data", style = "margin: 0; font-weight: bold;")
                   ),
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.1); padding: 15px; border-radius: 6px; margin-bottom: 20px; flex-grow: 1;",
                     tags$p(tags$strong("Supported formats:"), style = "margin-bottom: 10px;"),
                     tags$ul(
                       style = "margin: 0 0 15px 0; padding-left: 20px; font-size: 14px;",
                       tags$li("ZIP files (barcodes, matrix, features)"),
                       tags$li("H5 files (10X HDF5 format)"),
                       tags$li("H5AD files (AnnData / Scanpy)")
                     ),
                     tags$p(tags$strong("Perfect for:"), style = "margin-bottom: 10px;"),
                     tags$ul(
                       style = "margin: 0; padding-left: 20px; font-size: 14px;",
                       tags$li("Fresh CellRanger output"),
                       tags$li("Complete quality control pipeline"),
                       tags$li("Full preprocessing workflow")
                     )
                   ),
                   tags$div(
                     style = "margin-top: auto;",
                     fileInput('file', 'Select file with raw data',
                               accept = c('.zip', '.h5', '.h5ad', '.hdf5', 'application/x-hdf5')),
                     tags$p(class = "text-muted",
                            style = "color: rgba(255,255,255,0.8) !important; margin-top: -10px; font-size: 12px;",
                            "Supported: 10X ZIP, 10X HDF5 (.h5), or AnnData (.h5ad)")
                   )
                 )
          ),
          column(6,
                 tags$div(
                   class = "card-green",
                   style = "min-height: 400px; display: flex; flex-direction: column;",
                   tags$div(
                     style = "text-align: center; margin-bottom: 20px;",
                     tags$div(
                       style = "background-color: rgba(255,255,255,0.2); width: 80px; height: 80px; border-radius: 50%; display: inline-flex; align-items: center; justify-content: center; margin-bottom: 15px;",
                       icon("database", style = "font-size: 36px;")
                     ),
                     tags$h4("Load Processed Data", style = "margin: 0; font-weight: bold;")
                   ),
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.1); padding: 15px; border-radius: 6px; margin-bottom: 20px; flex-grow: 1;",
                     tags$p(tags$strong("Requirements:"), style = "margin-bottom: 10px;"),
                     tags$ul(
                       style = "margin: 0 0 15px 0; padding-left: 20px; font-size: 14px;",
                       tags$li("Seurat object (.rds file)"),
                       tags$li("AnnData object (.h5ad file - via schard)"),
                       tags$li("Normalized expression data"),
                       tags$li("Optional: Clustering results")
                     ),
                     tags$p(tags$strong("Perfect for:"), style = "margin-bottom: 10px;"),
                     tags$ul(
                       style = "margin: 0; padding-left: 20px; font-size: 14px;",
                       tags$li("Previously analyzed datasets"),
                       tags$li("Skip quality control steps"),
                       tags$li("Continue from saved analysis")
                     )
                   ),
                   tags$div(
                     style = "margin-top: auto;",
                     fileInput("load_seurat_file", "Select processed object (.rds or .h5ad)",
                               accept = c(".rds", ".h5ad", ".hdf5", "application/x-hdf5")),
                     tags$p(class = "text-muted", style = "color: rgba(255,255,255,0.8) !important; margin-top: -10px; font-size: 12px;", "Upload a pre-processed Seurat object saved as RDS file")
                   )
                 )
          )
        ),
        tags$div(
          style = "margin-top: 20px; background-color: #f8f9fa; padding: 20px; 
         border-radius: 8px; border: 2px solid #dee2e6;",
          tags$h4(
            icon("terminal"), " Single Dataset Loading Logs", 
            style = "color: #495057; margin-top: 0; margin-bottom: 15px; font-weight: bold;"
          ),
          tags$pre(
            style = "background-color: #1e272e; color: #00ff00; padding: 15px; 
           border-radius: 6px; font-family: 'Courier New', monospace; 
           font-size: 11px; max-height: 400px; overflow-y: auto; 
           margin: 0; border-left: 4px solid #667eea; white-space: pre-wrap;",
            verbatimTextOutput("single_loading_logs")
          )
        ),
        example_data_box("single"),
        fluidRow(
          style = "margin-top: 30px;",
          column(3, tags$div(style = "background-color: #e8eaf6; padding: 15px; border-radius: 6px; border-left: 4px solid #667eea; height: 100%;", tags$h6(icon("check-circle"), strong(" File Structure"), style = "color: #5e35b1; margin-top: 0;"), tags$p("Ensure proper 10X Genomics folder structure", style = "font-size: 12px; margin: 0; color: #555;"))),
          column(3, tags$div(style = "background-color: #f3e5f5; padding: 15px; border-radius: 6px; border-left: 4px solid #9c27b0; height: 100%;", tags$h6(icon("memory"), strong(" File Size"), style = "color: #6a1b9a; margin-top: 0;"), tags$p("Typical datasets: 100MB-5GB compressed", style = "font-size: 12px; margin: 0; color: #555;"))),
          column(3, tags$div(style = "background-color: #fff3e0; padding: 15px; border-radius: 6px; border-left: 4px solid #ff9800; height: 100%;", tags$h6(icon("clock"), strong(" Processing Time"), style = "color: #e65100; margin-top: 0;"), tags$p("Loading takes 30s-3min depending on size", style = "font-size: 12px; margin: 0; color: #555;"))),
          column(3, tags$div(style = "background-color: #e3f2fd; padding: 15px; border-radius: 6px; border-left: 4px solid #2196f3; height: 100%;", tags$h6(icon("save"), strong(" Save Progress"), style = "color: #1976d2; margin-top: 0;"), tags$p("Export processed data as RDS for later use", style = "font-size: 12px; margin: 0; color: #555;")))
        )
      ) ,
      
      ############################## Single/QC metrics and normalization  ##############################
      tabItem(
        tabName = "qc",
        
        tags$div(
          class = "header-violet",
          tags$h2("Quality Control & Filtering", class = "header-title"),
          tags$p("Assess data quality and filter low-quality cells", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8f5e9; padding: 20px; border-radius: 8px; border-left: 4px solid #11998e; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " Quality Control Overview", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Quality control is essential to remove low-quality cells and technical artifacts before downstream analysis."),
            tags$h5(icon("star"), " Key Metrics:", style = "color: #11998e; margin-top: 15px;"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(strong("Feature Selection:"), "Set the range of unique genes detected per cell. Generally, cells with very few genes may be empty droplets or poor quality cells, while those with too many genes might be doublets."),
              tags$li(strong("Mitochondrial Filtering:"), "Set the maximum percentage of mitochondrial genes. High mitochondrial content often indicates cell stress or death."),
              tags$li(strong("Visual Inspection:"), "Use the violin plots and scatter plots to evaluate the distribution of these metrics across your cells."),
              tags$li(strong("Apply Filters:"), "Click 'Apply QC Filters' to remove cells that don't meet the criteria."),
              tags$li(strong("Normalization:"), "Finally, normalize the filtered data to account for technical variations between cells.")
            ),
            tags$div(
              style = "background-color: #fff9c4; padding: 12px; border-radius: 6px; border-left: 3px solid #fbc02d; margin-top: 15px;",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Low feature count = empty droplets, High MT% = dying cells, Very high features = doublets", style = "margin: 0;")
            ),
            tags$p(
              style = "margin-top: 15px; font-size: 12px; font-style: italic;",
              "Note: The default parameters are general guidelines. Optimal values may vary depending on your specific experimental context."
            )
          )
        ),
        
        sidebarLayout(
          sidebarPanel(
            uiOutput("stats_cards_output"),
            
            div(
              style = "display: inline-block; width: 80%;",
              sliderInput("nFeature_range", "Unique genes detected in each cell", min = 0, max = 8000, value = c(200, 3500))
            ),
            div(
              style = "display: inline-block; width: 18%;",
              actionButton("exclamation1", label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-html` = "true", `data-content` = "The number of unique genes detected in each cell.<ul><li>Low-quality cells or empty droplets will often have very few genes</li><li>Cell doublets or multiplets may exhibit an aberrantly high gene count.</li><li>Be careful, in most cases it's better to stay between 200 and 2500</li></ul>")
            ),
            div(
              style = "display: inline-block; width: 80%;",
              sliderInput("percent.mt_max", "Maximum value for mitochondrial genome:", value = 5, min = 0, max = 100)
            ),
            div(
              style = "display: inline-block; width: 18%;",
              actionButton("exclamation3", label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-content` = "It is normal to have mitochondrial DNA contamination for single cell, but be careful with this value for single nuclei")
            ),
            div(
              style = "display: inline-block; width: 80%;",
              sliderInput(inputId = "scale_factor", label = "Scale factor", value = 10000, min = 500, max = 50000, step = 500)
            ),
            div(
              style = "display: inline-block; width: 18%;",
              actionButton("exclamation4", label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-content` = "Feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result")
            ),
            actionButton("show_plots", "Feature Scatter Plots", class = "btn-white-violet btn-full-width"),
            tags$br(),
            actionButton("apply_qc", "Apply QC Filters", class = "btn-white-violet btn-full-width"),
            tags$br(),
            numericInput("num_var_features", "Number of Variable Features:", value = 2000, min = 1),
            tags$br(),
            selectInput("normalization_method_single", 
                        "Normalization Method:", 
                        choices = c("LogNormalize (standard)" = "LogNormalize", 
                                    "SCTransform (recommended)" = "SCTransform",
                                    "CLR (centered log-ratio)" = "CLR", 
                                    "RC (relative counts)" = "RC"),
                        selected = "LogNormalize"),
            tags$br(),
            actionButton("normalize_data", "Normalize Data", class = "btn-gradient-green btn-full-width")
          ),
          mainPanel(
            tags$div(
              class = "box-gradient-green",
              tags$h4(icon("chart-bar"), " QC Metrics Distribution", class = "box-title"),
              checkboxInput("qc_show_points_single", "Show cells as points", value = TRUE),
              tags$div(
                class = "plot-container",
                style = "min-height: 500px;",
                plotOutput("vlnplot", height = "500px")
              )
            ),
            fluidRow(
              column(6,
                     tags$div(
                       class = "box-gradient-green box-thin",
                       tags$h4(icon("braille"), " UMI vs Mitochondrial %", class = "box-title"),
                       tags$div(
                         class = "plot-container-thin",
                         style = "min-height: 400px;",
                         plotOutput("scatter_plot1", height = "400px")
                       )
                     )
              ),
              column(6,
                     tags$div(
                       class = "box-gradient-green box-thin",
                       tags$h4(icon("braille"), " UMI vs Feature Count", class = "box-title"),
                       tags$div(
                         class = "plot-container-thin",
                         style = "min-height: 400px;",
                         plotOutput("scatter_plot2", height = "400px")
                       )
                     )
              )
            ),
            tags$div(
              class = "box-gradient-violet",
              tags$h4(icon("dna"), " Variable Features", class = "box-title"),
              tags$div(
                class = "plot-container",
                style = "min-height: 400px;",
                plotOutput("variable_feature_plot", height = "400px")
              )
            )
          )
        )
      ),
      
      ############################## Single/Scaling, PCA and elbow plot ##############################
      tabItem(
        tabName = "dimensional_reduction",
        
        tags$div(
          class = "header-violet",
          tags$h2("Dimensional Reduction", class = "header-title"),
          tags$p("Scale data and perform PCA to reduce dimensionality", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8f5e9; padding: 20px; border-radius: 8px; border-left: 4px solid #11998e; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Dimensional Reduction", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Dimensional reduction simplifies your data while preserving important biological signals."),
            tags$h5(icon("star"), " Key Steps:", style = "color: #11998e; margin-top: 15px;"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(strong("Data Scaling:"), "Normalize the expression of each gene to give equal weight to all genes, regardless of their expression level. This prevents highly-expressed genes from dominating the analysis."),
              tags$li(strong("Principal Component Analysis (PCA):"), "PCA identifies the main sources of variation in your data, transforming thousands of gene expression measurements into a smaller set of meaningful components."),
              tags$li(strong("Elbow Plot Analysis:"), "The elbow plot helps determine how many principal components to use in downstream analysis by showing the percentage of variance explained by each component.")
            ),
            tags$div(
              style = "background-color: #fff9c4; padding: 12px; border-radius: 6px; border-left: 3px solid #fbc02d; margin-top: 15px;",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " The 'elbow' in the plot indicates the optimal number of dimensions. Typically, 15-30 components are sufficient.", style = "margin: 0;")
            ),
            tags$p(
              style = "margin-top: 15px; font-size: 12px; font-style: italic;",
              "Note: More components capture more variation but increase computational time."
            )
          )
        ),
        
        fluidRow(
          column(4,
                 tags$div(
                   style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; border: 2px solid #dee2e6; min-height: 400px; display: flex; flex-direction: column;",
                   tags$h4(icon("cogs"), " Analysis Workflow", style = "color: #495057; margin-top: 0; margin-bottom: 20px; font-weight: bold;"),
                   tags$div(
                     style = "flex-grow: 1;",
                     tags$div(
                       style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #667eea; margin-bottom: 15px;",
                       tags$h5(icon("balance-scale"), " Step 1: Scale & PCA", style = "color: #667eea; margin-top: 0;"),
                       tags$p("Normalize gene expression and identify principal components", style = "font-size: 13px; margin: 0; color: #555;")
                     ),
                     tags$div(
                       style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #11998e; margin-bottom: 15px;",
                       tags$h5(icon("chart-line"), " Step 2: Elbow Plot", style = "color: #11998e; margin-top: 0;"),
                       tags$p("Determine optimal number of dimensions for clustering", style = "font-size: 13px; margin: 0; color: #555;")
                     ),
                     tags$div(
                       style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #f093fb;",
                       tags$h5(icon("info-circle"), " Next Steps", style = "color: #f093fb; margin-top: 0;"),
                       tags$p("After PCA, proceed to clustering analysis", style = "font-size: 13px; margin: 0; color: #555;")
                     )
                   ),
                   tags$div(
                     style = "margin-top: auto;",
                     actionButton("scale_button", tagList(icon("balance-scale"), " Scale Data & Run PCA"), class = "btn-gradient-violet btn-full-width"),
                     tags$br(),
                     actionButton("run_elbow", tagList(icon("chart-line"), " Show Elbow Plot"), class = "btn-gradient-green btn-full-width")
                   )
                 )
          ),
          column(8,
                 tags$div(
                   class = "box-gradient-green",
                   tags$h4(icon("table"), " PCA Results", class = "box-title"),
                   tags$div(
                     class = "plot-container",
                     style = "min-height: 200px; max-height: 200px; overflow-y: auto;",
                     verbatimTextOutput("pca_results")
                   )
                 ),
                 tags$div(
                   class = "box-gradient-violet",
                   tags$h4(icon("chart-area"), " PCA Loading Plot", class = "box-title"),
                   tags$div(
                     class = "plot-container",
                     style = "min-height: 300px;",
                     plotOutput("loading_plot", height = "300px")
                   )
                 )
          )
        ),
        
        fluidRow(
          style = "margin-top: 20px;",
          column(6,
                 tags$div(
                   class = "box-gradient-green",
                   tags$h4(icon("project-diagram"), " PCA Dimension Plot", class = "box-title"),
                   tags$div(
                     class = "plot-container",
                     style = "min-height: 400px;",
                     plotOutput("dim_plot", height = "400px")
                   )
                 )
          ),
          column(6,
                 tags$div(
                   class = "box-gradient-violet",
                   tags$h4(icon("chart-line"), " Elbow Plot", class = "box-title"),
                   tags$div(
                     class = "plot-container",
                     style = "min-height: 400px;",
                     plotOutput("elbow_plot", height = "400px")
                   )
                 )
          )
        )
      ),
      #######################################Single Dataset Analysis/Clustering####################################
      tabItem(
        tabName = "clustering",
        
        tags$div(
          class = "header-violet",
          tags$h2("Cell Clustering", class = "header-title"),
          tags$p("Group cells based on expression similarities", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8f5e9; padding: 20px; border-radius: 8px; border-left: 4px solid #11998e; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Cell Clustering", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("This step groups cells based on their gene expression similarities to identify distinct cell populations."),
            tags$h5(icon("star"), " Key Steps:", style = "color: #11998e; margin-top: 15px;"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(strong("Neighbor Finding:"), "First, we identify cells with similar gene expression profiles. The number of dimensions used affects how these similarities are calculated."),
              tags$li(strong("Clustering Resolution:"), "Then, we group similar cells together. The resolution parameter controls how finely the cells are grouped - higher values create more clusters."),
              tags$li(strong("Algorithm Selection:"), "Different clustering algorithms are available, each with its own approach to grouping cells.")
            ),
            tags$div(
              style = "background-color: #fff9c4; padding: 12px; border-radius: 6px; border-left: 3px solid #fbc02d; margin-top: 15px;",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Finding the right balance between too many and too few clusters is key for meaningful biological interpretation.", style = "margin: 0;")
            )
          )
        ),
        
        tags$div(
          style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; border: 2px solid #dee2e6; margin-bottom: 20px;",
          tags$h4(icon("cogs"), " Clustering Parameters", style = "color: #495057; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
          fluidRow(
            column(3,
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #667eea;",
                     tags$label("Number of Dimensions:", style = "font-weight: bold; color: #495057;"),
                     numericInput("dimension_1", NULL, value = 15, min = 1),
                     actionButton("run_neighbors", tagList(icon("project-diagram"), " Find Neighbors"), class = "btn-gradient-violet btn-full-width"),
                     actionButton("compute_3d_umap", tagList(icon("cube"), " Compute 3D UMAP"), class = "btn-white-violet btn-full-width", style = "margin-top: 8px;"),
                     tags$p(icon("info-circle"), " Use PCA dimensions", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                   )
            ),
            column(3,
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #11998e;",
                     tags$label("Resolution:", style = "font-weight: bold; color: #495057;"),
                     numericInput("resolution_step1", NULL, min = 0.1, max = 2, step = 0.1, value = 0.5),
                     actionButton("run_clustering", tagList(icon("sitemap"), " Find Clusters"), class = "btn-gradient-green btn-full-width"),
                     tags$p(icon("info-circle"), " Higher = more clusters", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                   )
            ),
            column(3,
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #764ba2;",
                     tags$label("Algorithm:", style = "font-weight: bold; color: #495057;"),
                     selectInput("algorithm_select_single", NULL, choices = list("Original Louvain" = 1, "Louvain with Multilevel Refinement" = 2, "SLM Algorithm" = 3)),
                     tags$label("Export Format:", style = "font-weight: bold; color: #495057; margin-top: 10px;"),
                     selectInput("umap_export_format", NULL, choices = c("TIFF" = "tiff", "PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg", "SVG" = "svg"), selected = "tiff"),
                     tags$p(icon("info-circle"), " Clustering method", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                   )
            ),
            column(3,
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #38ef7d;",
                     tags$label("Resolution (DPI):", style = "font-weight: bold; color: #495057;"),
                     numericInput("dpi_umap", NULL, value = 300, min = 72, max = 1200, step = 72),
                     downloadButton("downloadUMAP", "Download UMAP", class = "btn-white-green", style = "width: 100%; margin-top: 44px;"),
                     tags$p(icon("info-circle"), " Image quality", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                   )
            )
          )
        ),
        tags$div(
          class = "box-gradient-green",
          tags$h4(icon("project-diagram"), " Clustering Results", class = "box-title"),
          tags$div(
            style = "background-color: white; padding: 15px; border-radius: 6px; margin-bottom: 15px; color: #495057;",
            fluidRow(
              column(4, radioButtons("cluster_label_mode_single", "Cluster names:",
                                     choices = c("On the UMAP" = "plot", "Side legend" = "legend", "Both" = "both", "None" = "none"),
                                     selected = "legend", inline = TRUE)),
              column(2, numericInput("cluster_label_size_single", "Name size:", value = 5, min = 1, max = 20, step = 0.5)),
              column(2, checkboxInput("cluster_bold_labels_single", "Bold names", value = FALSE)),
              column(2, checkboxInput("cluster_dark_mode_single", tagList(icon("moon"), " Black background"), value = FALSE)),
              column(2, checkboxInput("remove_axes", "Remove Axes", FALSE))
            ),
            tags$p(icon("info-circle"), " Display options apply instantly and to the downloaded image; no need to re-run clustering.",
                   style = "font-size: 11px; color: #6c757d; margin: 0;")
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 600px;",
            plotOutput("clustering_plot", height = "600px")
          )
        ),
        # 3D UMAP panel, shown once "Compute 3D UMAP" has been run
        conditionalPanel(
          condition = "output.has_umap3d_single == true",
          tags$div(
            style = "background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800; margin-top: 20px;",
            tags$h4(icon("cube"), " 3D UMAP", style = "color: #ff9800; margin-top: 0; font-weight: bold;"),
            fluidRow(
              column(2, checkboxInput("umap3d_dark_clust_single", tagList(icon("moon"), " Dark mode"), value = FALSE)),
              column(2, checkboxInput("umap3d_hide_grid_clust_single", tagList(icon("border-none"), " Hide grid"), value = FALSE)),
              column(2, checkboxInput("umap3d_hide_axes_clust_single", tagList(icon("arrows-alt"), " Hide axes"), value = FALSE)),
              column(3, numericInput("umap3d_pt_size_clust_single", "Point size:", value = 2, min = 0.5, max = 20, step = 0.5)),
              column(3, tags$p(icon("info-circle"), " Click and drag to rotate", style = "font-size: 12px; color: #6c757d; margin-top: 10px;"))
            ),
            plotlyOutput("clustering_umap3d_single", height = "650px")
          )
        )
      ),
      ############################## Single/DoubletFinder ##############################
      tabItem(
        tabName = "doublet_detection",
        
        tags$div(
          class = "header-violet",
          tags$h2("Doublet Detection", class = "header-title"),
          tags$p("Identify and remove technical artifacts from your data", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8f5e9; padding: 20px; border-radius: 8px; border-left: 4px solid #11998e; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Doublet Detection", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Doublets are artifacts where two cells are captured together, leading to mixed expression profiles. DoubletFinder helps identify and remove these artifacts to ensure data quality."),
            tags$h5(icon("star"), " Key Parameters:", style = "color: #11998e; margin-top: 15px;"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(strong("Expected Rate:"), "Set based on your experimental protocol, typically 5-10% for 10x Genomics data."),
              tags$li(strong("Parameter Settings:"), "pN and pK values control the algorithm's sensitivity and specificity."),
              tags$li(strong("Visualization:"), "Results are shown on a UMAP plot where doublets are highlighted.")
            ),
            tags$div(
              style = "background-color: #fff9c4; padding: 12px; border-radius: 6px; border-left: 3px solid #fbc02d; margin-top: 15px;",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), "To use Doublet Finder, you must have normalized and scaled the data.", style = "margin: 0;")
            )
          )
        ),
        
        tags$div(
          style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; border: 2px solid #dee2e6; margin-bottom: 20px;",
          tags$h4(icon("sliders-h"), " Detection Parameters", style = "color: #495057; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
          fluidRow(
            column(3,
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #667eea;",
                     tags$label("Expected Doublet Rate (%):", style = "font-weight: bold; color: #495057;"),
                     numericInput("doublet_rate", NULL, value = 7.5, min = 0, max = 100, step = 0.5),
                     tags$p(icon("info-circle"), " Typical: 5-10% for 10x", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                   )
            ),
            column(3,
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #11998e;",
                     tags$label("Number of PCs:", style = "font-weight: bold; color: #495057;"),
                     numericInput("pc_use", NULL, value = 10, min = 1, max = 50),
                     tags$p(icon("info-circle"), " Principal components to use", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                   )
            ),
            column(3,
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #764ba2;",
                     tags$label("pN value:", style = "font-weight: bold; color: #495057;"),
                     numericInput("pN_value", NULL, value = 0.25, min = 0, max = 1, step = 0.05),
                     tags$p(icon("info-circle"), " Proportion of artificial doublets", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                   )
            ),
            column(3,
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #38ef7d;",
                     tags$label("pK selection:", style = "font-weight: bold; color: #495057;"),
                     radioButtons("doublet_pK_mode", NULL,
                                  choices = c("Auto (BCmvn)" = "auto", "Manual" = "manual"),
                                  selected = "auto", inline = TRUE),
                     conditionalPanel(
                       condition = "input.doublet_pK_mode == 'manual'",
                       numericInput("pK_value", NULL, value = 0.09, min = 0.0005, max = 0.3, step = 0.005)
                     ),
                     tags$p(icon("info-circle"), " PC neighborhood size (auto = maximum BCmvn)", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                   )
            )
          ),
          fluidRow(
            style = "margin-top: 15px;",
            column(6,
                   actionButton("run_doubletfinder", tagList(icon("search"), " Detect Doublets"), class = "btn-gradient-violet btn-full-width")
            ),
            column(6,
                   actionButton("remove_doublets", tagList(icon("trash-alt"), " Remove Doublets"), class = "btn-gradient-green btn-full-width")
            )
          )
        ),
        
        fluidRow(
          column(6,
                 tags$div(
                   class = "box-gradient-violet",
                   tags$h4(icon("project-diagram"), " Doublet UMAP", class = "box-title"),
                   tags$div(
                     class = "plot-container",
                     style = "min-height: 500px;",
                     plotOutput("doublet_umap", height = "500px")
                   )
                 )
          ),
          column(6,
                 tags$div(
                   class = "box-gradient-green",
                   tags$h4(icon("table"), " Detection Statistics", class = "box-title"),
                   tags$div(
                     class = "plot-container",
                     style = "min-height: 500px; display: flex; flex-direction: column;",
                     DTOutput("doublet_stats"),
                     tags$div(
                       style = "margin-top: auto; padding-top: 15px;",
                       downloadButton("download_doublet_results", "Download Results", class = "btn-white-green", style = "width: 100%;")
                     )
                   )
                 )
          )
        )
      ),
      
      
      ############################## Single/Visualization of expressed genes ##############################
      tabItem(
        tabName = "plots_genes_expressions",
        
        tags$div(
          class = "header-violet",
          tags$h2("Gene Expression Visualization", class = "header-title"),
          tags$p("Explore gene expression patterns across your data", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8eaf6; padding: 20px; border-radius: 8px; border-left: 4px solid #667eea; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Gene Expression Visualization", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Choose from different visualization methods to understand gene expression patterns."),
            tags$h5(icon("star"), " Available Visualizations:", style = "color: #667eea; margin-top: 15px;"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(strong("Feature Plot:"), "Shows gene expression on UMAP projection. Useful for seeing spatial expression patterns."),
              tags$li(strong("Violin Plot:"), "Displays expression distribution across clusters. Good for comparing expression levels."),
              tags$li(strong("Dot Plot:"), "Shows both expression level and percentage of expressing cells. Perfect for comparing multiple genes."),
              tags$li(strong("Ridge Plot:"), "Alternative to violin plots, better for visualizing many clusters.")
            ),
            tags$div(
              style = "background-color: #fff9c4; padding: 12px; border-radius: 6px; border-left: 3px solid #fbc02d; margin-top: 15px;",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Use Feature Plots for spatial patterns, Violin/Dot plots for quantitative comparisons", style = "margin: 0;")
            )
          )
        ),
        
        tags$div(
          style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; border: 2px solid #dee2e6; margin-bottom: 20px;",
          tags$h4(icon("dna"), " Gene Selection & Global Parameters", style = "color: #495057; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
          fluidRow(
            column(4,
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #667eea;",
                     tags$label("Select Genes:", style = "font-weight: bold; color: #495057;"),
                     pickerInput("gene_select", NULL, choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                     tags$label("Title Size:", style = "font-weight: bold; color: #495057; margin-top: 10px;"),
                     numericInput("title_text_size", NULL, value = 14, min = 8, max = 32),
                     tags$p(icon("info-circle"), " Search and select genes", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                   )
            ),
            column(4,
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #11998e;",
                     tags$label("Assay:", style = "font-weight: bold; color: #495057;"),
                     selectInput("viz_assay", NULL, choices = NULL, selected = "RNA"),
                     tags$label("Axis Text Size:", style = "font-weight: bold; color: #495057; margin-top: 10px;"),
                     numericInput("axis_text_size", NULL, value = 12, min = 6, max = 30),
                     tags$label("Axis Line Width:", style = "font-weight: bold; color: #495057; margin-top: 10px;"),
                     numericInput("axis_line_width", NULL, value = 1, min = 0.5, max = 3, step = 0.1),
                     checkboxInput("bold_axis_text", "Bold text", value = FALSE),
                     tags$p(icon("info-circle"), " Assay and axis styling", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                   )
            ),
            column(4,
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #764ba2;",
                     tags$label("Export Settings:", style = "font-weight: bold; color: #495057;"),
                     selectInput("plot_format", NULL, choices = c("PNG" = "png", "JPEG" = "jpeg", "TIFF" = "tiff", "SVG" = "svg", "PDF" = "pdf"), selected = "png"),
                     tags$label("DPI:", style = "font-weight: bold; color: #495057; margin-top: 10px;"),
                     numericInput("dpi_plot", NULL, value = 300, min = 72, step = 72),
                     tags$p(icon("info-circle"), " Export settings", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                   )
            )
          )
        ),
        
        tags$div(
          class = "box-gradient-green",
          tags$h4(icon("map-marked-alt"), " Feature Plot", class = "box-title"),
          fluidRow(
            column(12,
                   textInput("gene_list_feature", "Selected genes for FeaturePlot (separated by , ; or spaces):", value = "", width = "100%")
            )
          ),
          fluidRow(
            column(3,
                   actionButton("show_feature", tagList(icon("eye"), " Display Plot"), class = "btn-white-violet btn-full-width"),
                   checkboxInput("add_noaxes_feature", "Remove Axes", FALSE),
                   checkboxInput("enable_3d_feature", "Display in 3D", FALSE)
            ),
            column(3,
                   selectInput("min_cutoff", "Minimum Cutoff:", choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30", "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80", "q90" = "q90", "q99" = "q99"), selected = NA),
                   checkboxInput("add_nolegend_feature", "Remove Legend", FALSE),
                   checkboxInput("hide_grid_feature", "Hide grid (3D)", FALSE)
            ),
            column(3,
                   selectInput("max_cutoff", "Maximum Cutoff:", choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30", "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80", "q90" = "q90", "q99" = "q99"), selected = NA),
                   checkboxInput("show_coexpression", "Show co-expression", value = FALSE),
                   checkboxInput("dark_mode_feature", "Dark mode", FALSE)
            ),
            column(3,
                   downloadButton("downloadFeaturePlot", "Download", class = "btn-white-violet", style = "width: 100%;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px; margin-top: 15px;",
            conditionalPanel(
              condition = "!input.enable_3d_feature",
              plotOutput("feature_plot", height = "600px")
            ),
            conditionalPanel(
              condition = "input.enable_3d_feature",
              plotly::plotlyOutput("feature_plot_3d", height = "600px")
            )                  )
        ),
        
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("braille"), " Dot Plot", class = "box-title"),
          fluidRow(
            column(12,
                   textInput("gene_list_dotplot", "Selected genes (separated by , ; or spaces):", value = "", width = "100%"),
                   selectizeInput("cluster_order_dot", "Order of Clusters:", 
                                  choices = NULL, multiple = TRUE, width = "100%",
                                  options = list(plugins = list('remove_button', 'drag_drop')))
            )
          ),
          fluidRow(
            column(4,
                   actionButton("show_dot", tagList(icon("eye"), " Display"), class = "btn-white-violet btn-full-width"),
                   checkboxInput("add_noaxes_dot", "Remove Axes", FALSE),
                   checkboxInput("add_nolegend_dot", "Remove Legend", FALSE)
            ),
            column(4,
                   checkboxInput("invert_axes_single", "Invert Axes", value = TRUE),
                   checkboxInput("dotplot_scale_single", "Scale expression (z-score across groups)", value = TRUE),
                   sliderInput("dot_scale_single", "Dot size scale:", 
                               min = 0.5, max = 15, value = 7, step = 0.5)
            ),
            column(4,
                   selectInput("color_palette_dotplot_single", "Color palette:",
                               choices = c("Default (Seurat)" = "default",
                                           "RdYlBu", "Blues", "Reds", "Greens", 
                                           "Spectral", "PuOr", "BrBG"),
                               selected = "default"),
                   downloadButton("downloadDotPlot", "Download", class = "btn-white-violet", style = "width: 100%;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px; margin-top: 15px;",
            plotOutput("dot_plot", height = "500px")
          )
        ),
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("chart-bar"), " Violin Plot", class = "box-title"),
          fluidRow(
            column(12,
                   textInput("gene_list_vln", "Selected genes (separated by , ; or spaces):", value = "", width = "100%")
            )
          ),
          fluidRow(
            column(12,
                   selectInput("cluster_order_vln", "Order of Clusters:", choices = NULL, multiple = TRUE, width = "100%")
            )
          ),
          fluidRow(
            column(2,
                   actionButton("show_vln", tagList(icon("eye"), " Display"), class = "btn-white-green btn-full-width")
            ),
            column(2,
                   checkboxInput("hide_vln_points", "Hide Points", FALSE)
            ),
            column(2,
                   checkboxInput("add_noaxes_vln", "Remove Axes", FALSE)
            ),
            column(2,
                   checkboxInput("add_nolegend_vln", "Remove Legend", FALSE)
            ),
            column(4,
                   downloadButton("downloadVlnPlot", "Download", class = "btn-white-green", style = "width: 100%;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px; margin-top: 15px;",
            plotOutput("vln_plot", height = "500px")
          )
        ),
        
        
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("mountain"), " Ridge Plot", class = "box-title"),
          fluidRow(
            column(12,
                   textInput("gene_list_ridge_plot", "Selected genes (separated by , ; or spaces):", value = "", width = "100%")
            )
          ),
          fluidRow(
            column(3,
                   actionButton("show_ridge", tagList(icon("eye"), " Display"), class = "btn-white-green btn-full-width")
            ),
            column(3,
                   checkboxInput("add_noaxes_ridge", "Remove Axes", FALSE)
            ),
            column(3,
                   checkboxInput("add_nolegend_ridge", "Remove Legend", FALSE)
            ),
            column(3,
                   downloadButton("download_ridge_plot", "Download", class = "btn-white-green", style = "width: 100%;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px; margin-top: 15px;",
            plotOutput("ridge_plot", height = "500px")
          )
        ),
        
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("table"), " Cell Expression Analysis", class = "box-title"),
          fluidRow(
            column(4,
                   pickerInput("gene_select_genes_analysis", "Select Genes:", choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE))
            ),
            column(4,
                   numericInput("expression_threshold_single", "Minimum Expression Threshold:", value = 1, min = 0, max = 10),
                   actionButton("analyze_btn", tagList(icon("calculator"), " Analyze Expression"), class = "btn-white-violet btn-full-width")
            ),
            column(4,
                   downloadButton("download_genes_number_expresion", "Download Results", class = "btn-white-violet", style = "width: 100%;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 300px; margin-top: 15px;",
            dataTableOutput("expression_summary")
          )
        )
      ),
      
      ############################## Single/Heatmaps ##############################
      tabItem(
        tabName = "heatmaps_dualexpression",
        
        tags$div(
          class = "header-violet",
          tags$h2("Heatmaps", class = "header-title"),
          tags$p("Advanced expression analysis with heatmaps", class = "header-subtitle")
        ),
        tags$div(
          style = "background-color: #e8f5e9; padding: 20px; border-radius: 8px; border-left: 4px solid #2c5f2d; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(
              icon("info-circle"), 
              " About Heatmap Visualization", 
              tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;")
            )
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p(
              "Visualize average gene expression across multiple genes and clusters simultaneously.",
              style = "margin-bottom: 20px; font-size: 1.05em;"
            ),
            tags$div(
              style = "padding: 15px; background-color: white; border-radius: 6px;",
              tags$h5(icon("th"), " Heatmap Analysis", style = "color: #2c5f2d; margin-top: 0; margin-bottom: 10px;"),
              tags$p(
                tags$strong("What it shows:"), " Each column represents the mean expression of all cells within a cluster (or cluster-dataset combination). Each row is one gene.",
                style = "margin-bottom: 8px;"
              ),
              tags$ul(
                style = "margin-left: 20px; margin-bottom: 0;",
                tags$li(tags$strong("Z-score scaling:"), " Shows relative expression (high/low) for each gene across clusters. Yellow = high, purple = low."),
                tags$li(tags$strong("Split by:"), " Compare the same clusters across different conditions (e.g., treated vs untreated)."),
                tags$li(tags$strong("Cluster ordering:"), " Drag and reorder clusters in the selection field to group related cell types together.")
              )
            )
          )
        ),
        tags$div(
          class = "box-gradient-green",
          tags$h4(icon("th"), " Heatmap Analysis", class = "box-title"),
          fluidRow(
            column(4,
                   tags$h5(icon("cog"), " Analysis Settings", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   selectInput("split_by_heatmap_single", "Split by (second axis, optional):", choices = NULL),
                   conditionalPanel(
                     condition = "input.split_by_heatmap_single != 'None' && input.split_by_heatmap_single != ''",
                     radioButtons("heatmap_order_single", "Panels:", choices = c("One panel per cluster" = "cluster", "One panel per split value" = "split"), selected = "cluster", inline = TRUE)
                   ),
                   tags$hr(style = "margin: 10px 0;"),
                   radioButtons("heatmap_scale_single", "Values:", choices = c("Row z-score" = "zscore", "Average expression (log)" = "none"), selected = "zscore", inline = TRUE),
                   selectInput("color_palette_heatmap_single", "Color Palette:",
                               choices = c("Viridis" = "viridis", "Magma" = "magma", "Inferno" = "inferno",
                                           "Plasma" = "plasma", "Red-Yellow-Blue" = "RdYlBu",
                                           "Red-White-Blue" = "RdBu", "Blue-White-Red" = "BlueRed",
                                           "Yellow-Red" = "YellowRed"),
                               selected = "viridis")
            ),
            column(4,
                   tags$h5(icon("layer-group"), " Selection", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   selectizeInput("clusters_heatmap_single", 
                                  "Select Clusters (in order):",
                                  choices = NULL, multiple = TRUE,
                                  options = list(plugins = list('remove_button', 'drag_drop'), persist = FALSE)),
                   tags$hr(style = "margin: 10px 0;"),
                   radioButtons("gene_source_heatmap_single", "Genes:",
                                choices = c("Manual list" = "manual", "Top markers per cluster" = "cluster_markers",
                                            "Top DE genes between two groups (per cluster)" = "condition_de"),
                                selected = "manual"),
                   conditionalPanel(
                     condition = "input.gene_source_heatmap_single == 'manual'",
                     textInput("gene_select_heatmap", "Select Genes (separated by , ; or spaces):")
                   ),
                   conditionalPanel(
                     condition = "input.gene_source_heatmap_single == 'condition_de'",
                     selectizeInput("heatmap_group1_single", "Group 1 (values of the split column):", choices = NULL, multiple = TRUE),
                     selectizeInput("heatmap_group2_single", "Group 2:", choices = NULL, multiple = TRUE),
                     selectInput("heatmap_de_direction_single", "Direction:", choices = c("Up and down" = "both", "Higher in Group 1" = "up", "Higher in Group 2" = "down"), selected = "both"),
                     tags$small("Genes ranked by a cell-level Wilcoxon test within each cluster: use it to choose genes to show, and the pseudobulk DE for statistics.", style = "color: #666;")
                   ),
                   conditionalPanel(
                     condition = "input.gene_source_heatmap_single != 'manual'",
                     numericInput("n_top_genes_single", "Number of genes per cluster:", value = 10, min = 1, max = 50, step = 1),
                     downloadButton("download_heatmap_markers_single", "Export genes CSV", class = "btn-sm btn-info")
                   ),
                   tags$hr(style = "margin: 10px 0;"),
                   actionButton("generateCustomHeatmap", tagList(icon("fire"), " Generate Heatmap"), 
                                class = "btn-white-green", style = "width: 100%;")
            ),
            column(4,
                   tags$h5(icon("download"), " Export Options", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   numericInput("dpi_heatmap_single", "Resolution (DPI):", value = 300, min = 72, step = 72),
                   selectInput("format_heatmap_single", "Format:", 
                               choices = c("PNG" = "png", "PDF" = "pdf", "TIFF" = "tiff", "JPEG" = "jpeg"), selected = "png"),
                   downloadButton("download_heatmap_single", "Download Heatmap", class = "btn-white-green", style = "width: 100%;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "margin-top: 15px; overflow-y: auto; max-height: 1000px; border: 1px solid #ddd;",
            plotOutput("heatmap_single")
          )
        )
      ),
      ######################################### Co-expression Analysis - Multiple Datasets #####################################################"
      tabItem(
        tabName = "coexpression_single",
        tags$div(
          class = "header-violet",
          tags$h2("Co-expression Analysis", class = "header-title"),
          tags$p("Gene co-expression statistics and expression-based differential analysis", class = "header-subtitle")
        ),
        # --- Scatter plot ---
        tags$div(
          class = "box-gradient-violet",
          tags$h4(icon("braille"), " Gene Co-expression Scatter", class = "box-title"),
          fluidRow(
            column(3,
                   tags$h5(icon("dna"), " Gene Selection", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   textInput("feature1_select", "First Gene:", ""),
                   textInput("feature2_select", "Second Gene:", "")
            ),
            column(3,
                   tags$h5(icon("palette"), " Display", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   selectInput("color_by_scatter_single", "Color by:", choices = c("Cluster" = "cluster"), selected = "cluster"),
                   selectInput("split_by_scatter_single", "Split by (optional):", choices = c("None" = "None"), selected = "None")
            ),
            column(3,
                   tags$h5(icon("filter"), " Filters", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   textInput("scatter_text_clusters", "Clusters (optional):", placeholder = "e.g., IIa, IIb"),
                   tags$hr(style = "margin: 10px 0;"),
                   numericInput("threshold_gene1_scatter_single", "Gene 1 Threshold:", value = 0, min = 0, step = 0.1),
                   numericInput("threshold_gene2_scatter_single", "Gene 2 Threshold:", value = 0, min = 0, step = 0.1)
            ),
            column(3,
                   tags$h5(icon("download"), " Export", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   selectInput("format_scatter_single", "Format:", choices = c("PNG" = "png", "JPEG" = "jpeg", "TIFF" = "tiff", "SVG" = "svg", "PDF" = "pdf"), selected = "png"),
                   numericInput("dpi_scatter_single", "Resolution (DPI):", value = 300, min = 72, step = 72),
                   downloadButton("download_scatter_single", "Download Plot", class = "btn-white-green", style = "width: 100%;")
            )
          ),
          fluidRow(
            column(12,
                   tags$hr(style = "margin: 15px 0;"),
                   actionButton("generateScatter_single", tagList(icon("braille"), " Generate Plot"), class = "btn-white-green", style = "width: 100%; height: 45px; font-size: 16px; font-weight: bold;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px; margin-top: 15px;",
            plotOutput("scatter_plot_single", height = "500px")
          )
        ),
        # --- Co-expression stats + DE ---
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("project-diagram"), " Gene Co-expression Analysis", class = "box-title"),
          fluidRow(
            column(width = 4,
                   textAreaInput("gene_text_coexpression_single", "Enter genes (separated by , ; or spaces):", value = "", placeholder = "e.g., Pax7, Myod1", height = "80px", width = "100%"),
                   selectInput("coexpr_plot_type", "Visualization:", choices = c("Stacked Bar" = "bar", "Heatmap" = "heat"), selected = "bar")
            ),
            column(width = 4,
                   radioButtons("coexpr_threshold_unit_single", "Thresholds on:", choices = c("Raw UMI counts (robust to ambient RNA)" = "counts", "Normalized expression (log)" = "data"), selected = "counts"),
                   uiOutput("coexpr_thresholds_ui_single")
            ),
            column(width = 4,
                   selectInput("coexpr_background_clusters_single", "Background cluster(s): cells that should not express these genes", choices = NULL, multiple = TRUE),
                   actionButton("estimate_coexpr_thresholds_single", tagList(icon("magic"), " Estimate thresholds from background"), class = "btn-white-violet btn-full-width"),
                   tags$small("Sets each threshold just above the 95th percentile of the gene in the background cells (ambient RNA level)", style = "color: #666;"),
                   uiOutput("coexpr_background_info_single")
            )
          ),
          tags$hr(),
          tags$div(
            style = "background-color: white; padding: 15px; border-radius: 8px;",
            tabsetPanel(
              id = "coexpr_main_tabs_single", type = "pills",
              tabPanel(
                tagList(icon("chart-bar"), " Co-expression Statistics"), value = "coexpr_stats_s",
                tags$div(
                  style = "background-color: #f5f5f5; padding: 12px 15px; border-radius: 6px; margin: 15px 0; border-left: 3px solid #6a4c93;",
                  tags$p(tags$b("What this does:"), " For each cluster, counts how many cells express each gene individually and how many co-express all genes simultaneously. Reports exact percentages per group.",
                         style = "margin: 0; font-size: 13px;")
                ),
                fluidRow(
                  column(width = 4,
                         selectInput("coexpr_groupby_col_single", "Group by column:", choices = NULL),
                         tags$small("Metadata column to restrict cells by", style = "color: #666;")                  ),
                  column(width = 4,
                         selectInput("coexpr_filter_groups_single", "Restrict to group(s):", choices = NULL, multiple = TRUE),
                         tags$small("Leave empty to use all cells", style = "color: #666;")
                  )
                ),
                fluidRow(
                  column(width = 4, actionButton("analyze_coexpression_single", tagList(icon("calculator"), " Analyze Co-expression"), class = "btn-white-violet btn-full-width", style = "margin-bottom: 10px;")),
                  column(width = 4, downloadButton("download_coexpression_table_single", "Download Table", style = "width: 100%;")),
                  column(width = 4, downloadButton("download_coexpression_plot_single", "Download Plot", style = "width: 100%;"))
                ),
                tags$hr(),
                tabsetPanel(
                  tabPanel("Results Table", DTOutput("gene_coexpression_table_single")),
                  tabPanel("Visualization", plotOutput("gene_coexpression_plot_single", height = "500px"))
                )
              ),
              tabPanel(
                tagList(icon("dna"), " DE: Expressors vs Rest"), value = "coexpr_de_s",
                tags$div(
                  style = "background-color: #f5f5f5; padding: 12px 15px; border-radius: 6px; margin: 15px 0; border-left: 3px solid #6a4c93;",
                  tags$p(tags$b("What this does:"), " Within the selected cluster(s), cells are split into two groups based on your gene(s) and thresholds, then differential expression (FindMarkers) identifies what distinguishes the two populations.",
                         style = "margin: 0 0 6px 0; font-size: 13px;"),
                  tags$p(tags$b("1 gene:"), " Cells expressing the gene (> threshold) vs cells that do not.",
                         style = "margin: 0 0 3px 0; font-size: 12px; color: #555;"),
                  tags$p(tags$b("2+ genes:"), " Cells co-expressing ALL genes (> thresholds) vs the rest (partial expressors + negatives).",
                         style = "margin: 0; font-size: 12px; color: #555;")
                ),
                fluidRow(
                  column(4,
                         selectInput("coexpr_de_cluster_single", "Cluster(s) to analyze:", choices = NULL, multiple = TRUE),
                         tags$p("Select one or more clusters. Cells are pooled if multiple clusters are selected.",
                                style = "color: #888; font-size: 11px; margin-top: -6px;"),
                         tags$br(),
                         actionButton("run_coexpr_de_single", tagList(icon("dna"), " Run DE"), class = "btn-white-violet btn-full-width"),
                         tags$br(), tags$br(),
                         uiOutput("coexpr_de_info_single"),
                         tags$br(),
                         downloadButton("download_coexpr_de_single", "Download DE table", style = "width: 100%;")
                  ),
                  column(8, DTOutput("coexpr_de_table_single"))
                )
              )
            ) )
        )
      ),
      ############################## Assigning cell type##############################
      tabItem(
        tabName = "assigning_cell_type_identity",
        
        tags$div(
          class = "header-violet",
          tags$h2("Cell Type Assignment", class = "header-title"),
          tags$p("Annotate clusters with biological identities", class = "header-subtitle")
        ),
        tags$div(
          style = "background-color: #e8f5e9; padding: 20px; border-radius: 8px; border-left: 4px solid #11998e; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Cell Type Assignment", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Give each cluster a biological identity. Three routes are available and can be combined: automated label transfer from your own annotated reference, automated annotation against built-in references (SingleR), and manual renaming based on marker genes. Automated labels are suggestions to be confirmed, not a final annotation."),
            tags$h5(icon("route"), " Annotation methods", style = "color: #11998e; margin-top: 15px;"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(strong("Your reference dataset (recommended):"), " upload an annotated Seurat object (.rds) and choose the column holding the cell types. Seurat anchors (FindTransferAnchors + TransferData) project the reference labels onto each cell and give a per-cell prediction score. Use a reference of the same species, tissue and technology (single-nucleus vs single-cell)."),
              tags$li(strong("Built-in references (SingleR):"), " each cell is correlated with bulk expression profiles of sorted cell types from celldex (mouse: 2 references, human: 5). Labels are broad and mostly immune: suitable for blood/immune cells or as a first pass."),
              tags$li(strong("Manual renaming:"), " in the Manual renaming tab, select a cluster and type its name. Typing an existing name merges the two clusters.")
            ),
            tags$h5(icon("list-check"), " Workflow for automated annotation", style = "color: #11998e; margin-top: 15px;"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(strong("Run:"), " each cell receives a label; each cluster (current cluster column) receives the majority label of its cells. Nothing is renamed at this step."),
              tags$li(strong("Review:"), " the cluster summary gives the agreement (% of cells carrying the majority label), the unassigned fraction (cells below the score threshold or pruned), the median score/delta and a High/Low confidence flag. Red rows must be checked with marker genes."),
              tags$li(strong("Apply:"), " 'Apply Labels to Clusters' renames the clusters (ClusterIdents). Low-confidence clusters can keep their current name, and clusters sharing a label can be numbered (Myonuclei_1, Myonuclei_2)."),
              tags$li(strong("Correct:"), " rename any cluster manually afterwards, or click 'Undo' (Manual renaming tab) to restore the previous names.")
            ),
            tags$p(icon("database"), " Per-cell predictions are stored as metadata columns (RefTransfer_* and SingleR_*: label, pruned_label, score, delta, cluster_label). They can be used as grouping variables in the other tabs and downloaded as CSV (Summary / Per cell).",
                   style = "font-size: 13px; margin-top: 10px;"),
            tags$div(
              style = "background-color: #ffebee; padding: 12px; border-radius: 6px; border-left: 3px solid #e53935; margin-top: 15px;",
              tags$p(icon("triangle-exclamation"), tags$strong(" Limitations of automated annotation"), style = "margin: 0 0 6px 0;"),
              tags$ul(
                style = "margin: 0; padding-left: 18px; line-height: 1.6;",
                tags$li("A method can only return labels present in its reference. A population absent from the reference receives the closest label, often with a high score (e.g. skeletal muscle myonuclei labelled 'Cardiomyocytes' and FAPs 'Fibroblasts' with MouseRNAseqData)."),
                tags$li("A high prediction score means that the reference neighbours agree, not that the cell type is correct."),
                tags$li("Rare, transitional or disease-specific states are usually merged into the nearest reference type.")
              )
            ),
            tags$h5(icon("palette"), " Display", style = "color: #11998e; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(strong("Plot Settings:"), " choose the colouring column, reduction, title, label size and point size of the final UMAP."),
              tags$li(strong("Cluster Colors:"), " assign a colour to each cluster; colours are stored in the object and reused across plots.")
            ),
            tags$div(
              style = "background-color: #fff9c4; padding: 12px; border-radius: 6px; border-left: 3px solid #fbc02d; margin-top: 15px;",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " check each label with known markers: display plots from the gene visualization section by selecting the plot type in the Alternative Visualizations box.", style = "margin: 0;")
            )
          )
        ),
        tags$div(
          style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; border: 2px solid #dee2e6; margin-bottom: 20px;",
          tags$h4(icon("cogs"), " Cluster Identity Assignment", style = "color: #495057; margin-top: 0; margin-bottom: 8px; font-weight: bold;"),
          tags$p("Name clusters manually or from an automated annotation. Automated predictions are first stored as new metadata columns; cluster names only change when you click 'Apply Labels'.",
                 style = "font-size: 13px; color: #6c757d; margin-bottom: 15px;"),
          tabsetPanel(
            id = "annotation_method_tabs_single", selected = "manual",
            ############ Tab 1: manual renaming ############
            tabPanel(
              tagList(icon("pen"), " Manual renaming"), value = "manual",
              tags$div(
                style = "margin-top: 15px;",
                fluidRow(
                  column(width = 6,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #667eea;",
                           tags$h5("Rename Clusters", style = "color: #667eea; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           fluidRow(
                             column(7,
                                    selectInput("select_cluster", "Select cluster:", choices = NULL),
                                    textInput("rename_single_cluster", "New name:")
                             ),
                             column(5,
                                    actionButton("rename_single_cluster_button", tagList(icon("check"), " Apply Name"), class = "btn-gradient-violet btn-full-width", style = "margin-top: 25px;"),
                                    actionButton("undo_cluster_single", tagList(icon("undo"), " Undo"), class = "btn-white-violet btn-full-width", style = "margin-top: 40px;")
                             )
                           ),
                           tags$p(icon("info-circle"), " Undo also reverts the last 'Apply Labels' of the reference or SingleR tabs.", style = "font-size: 11px; color: #6c757d; margin: 0;")
                         )
                  ),
                  column(width = 6,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #11998e; font-size: 13px; color: #555;",
                           tags$h5(icon("lightbulb"), " How it works", style = "color: #11998e; margin-top: 0; margin-bottom: 10px; font-weight: bold;"),
                           tags$ul(
                             style = "padding-left: 18px; line-height: 1.7; margin: 0;",
                             tags$li("Select a cluster, type its name, click 'Apply Name'."),
                             tags$li("Typing a name that already exists merges the two clusters."),
                             tags$li("Names are written to the ClusterIdents column and used by all other tabs.")
                           )
                         )
                  )
                )
              )
            ),
            ############ Tab 2: label transfer from the user's reference ############
            tabPanel(
              tagList(icon("database"), " Your reference dataset"), value = "transfer",
              tags$div(
                style = "margin-top: 15px;",
                tags$div(
                  style = "background-color: #fff3cd; padding: 10px 15px; border-radius: 6px; border-left: 4px solid #f0ad4e; margin-bottom: 15px; font-size: 13px;",
                  tags$div(
                    style = "cursor: pointer;", onclick = "$(this).next().slideToggle();",
                    icon("triangle-exclamation"), tags$strong(" Labels can only come from your reference: a population missing from it gets the closest reference label. "),
                    tags$span("Show limitations", style = "text-decoration: underline; color: #8a6d3b;")
                  ),
                  tags$ul(
                    style = "display: none; margin: 8px 0 0 0; padding-left: 18px; line-height: 1.6;",
                    tags$li("A missing population is often assigned with a high prediction score: the score measures agreement between neighbouring reference cells, not how well the cell matches the reference."),
                    tags$li("Use a reference of the same species, tissue and technology (single-nucleus vs single-cell data differ, e.g. myonuclei are absent from single-cell atlases)."),
                    tags$li("Check every cluster with marker genes before applying labels.")
                  )
                ),
                fluidRow(
                  column(width = 4,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #667eea; min-height: 360px;",
                           tags$h5(tags$span("1", style = "display: inline-block; width: 24px; height: 24px; border-radius: 50%; background: #667eea; color: white; text-align: center; line-height: 24px; font-size: 13px; margin-right: 6px;"),
                                   "Load the reference", style = "color: #667eea; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           fileInput("annotation_reference_file_single", "Annotated Seurat object (.rds):", accept = ".rds"),
                           selectInput("annotation_reference_label_single", "Column containing the cell types:", choices = NULL),
                           tags$p(icon("circle-question"), " The metadata column whose values (cell types) will be transferred to your cells.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           uiOutput("annotation_reference_info_single")
                         )
                  ),
                  column(width = 4,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #11998e; min-height: 360px;",
                           tags$h5(tags$span("2", style = "display: inline-block; width: 24px; height: 24px; border-radius: 50%; background: #11998e; color: white; text-align: center; line-height: 24px; font-size: 13px; margin-right: 6px;"),
                                   "Transfer the labels", style = "color: #11998e; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           numericInput("transfer_min_score_single", "Min. prediction score per cell (0-1):", value = 0.5, min = 0, max = 1, step = 0.05),
                           tags$p(icon("circle-question"), " Each cell gets the reference label with the highest score (0-1). Cells below this value are marked 'Unassigned'. Higher = stricter; lower it (e.g. 0.3) if many cells end up unassigned with a reference that has many close subtypes. Re-run after changing it.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           numericInput("transfer_max_cells_single", "Max. reference cells per label (0 = all):", value = 500, min = 0, max = 10000, step = 100),
                           tags$p(icon("circle-question"), " Each reference cell type is randomly reduced to this number of cells (seed 42): faster, and abundant types do not outweigh rare ones. 0 keeps every cell (slower, more memory).", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           numericInput("transfer_dims_single", "PCA dimensions:", value = 30, min = 5, max = 100, step = 5),
                           tags$p(icon("circle-question"), " Number of principal components used to match your cells to the reference. 30 suits most data; 40-50 for references with many fine subtypes.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           actionButton("run_transfer_single", tagList(icon("play"), " Transfer Labels"), class = "btn-gradient-violet btn-full-width"),
                           uiOutput("transfer_status_single")
                         )
                  ),
                  column(width = 4,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #764ba2; min-height: 360px;",
                           tags$h5(tags$span("3", style = "display: inline-block; width: 24px; height: 24px; border-radius: 50%; background: #764ba2; color: white; text-align: center; line-height: 24px; font-size: 13px; margin-right: 6px;"),
                                   "Review & apply", style = "color: #764ba2; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           numericInput("transfer_min_agreement_single", "Min. agreement for a confident cluster (%):", value = 50, min = 10, max = 100, step = 5),
                           tags$p(icon("circle-question"), " A cluster is confident when at least this % of its cells share the same label and at most 30% of its cells are unassigned. Changes the summary immediately, no re-run needed.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           checkboxInput("transfer_keep_low_single", "Keep current name for low-confidence clusters", value = TRUE),
                           tags$p(icon("circle-question"), " Checked: only confident (green) clusters are renamed. Unchecked: every cluster takes its majority label.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           checkboxInput("transfer_unique_single", "Number clusters sharing a label (Myonuclei_1, Myonuclei_2...)", value = TRUE),
                           tags$p(icon("circle-question"), " Unchecked: clusters with the same label are merged into one.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           actionButton("apply_transfer_labels_single", tagList(icon("tags"), " Apply Labels to Clusters"), class = "btn-gradient-violet btn-full-width"),
                           fluidRow(
                             style = "margin-top: 10px;",
                             column(6, downloadButton("download_transfer_clusters_single", "Summary", class = "btn-white-violet", style = "width: 100%;")),
                             column(6, downloadButton("download_transfer_cells_single", "Per cell", class = "btn-white-violet", style = "width: 100%;"))
                           )
                         )
                  )
                ),
                tags$div(
                  style = "background-color: white; padding: 15px; border-radius: 6px; margin-top: 15px;",
                  tags$h5(icon("table"), " Cluster summary", style = "color: #495057; margin-top: 0; font-weight: bold;"),
                  tags$p("Agreement: % of cells carrying the majority label. Unassigned: % of cells below the minimum score. Score: median prediction score. Delta: gap between the best and second-best label. Red rows need checking with marker genes.",
                         style = "font-size: 12px; color: #6c757d;"),
                  DTOutput("transfer_summary_single"),
                  fluidRow(
                    style = "margin-top: 15px;",
                    column(width = 6, plotOutput("transfer_umap_labels_single", height = "500px")),
                    column(width = 6, plotOutput("transfer_umap_score_single", height = "500px"))
                  )
                )
              )
            ),
            ############ Tab 3: SingleR with built-in celldex references ############
            tabPanel(
              tagList(icon("book"), " Built-in references (SingleR)"), value = "singler",
              tags$div(
                style = "margin-top: 15px;",
                tags$div(
                  style = "background-color: #fff3cd; padding: 10px 15px; border-radius: 6px; border-left: 4px solid #f0ad4e; margin-bottom: 15px; font-size: 13px;",
                  tags$div(
                    style = "cursor: pointer;", onclick = "$(this).next().slideToggle();",
                    icon("triangle-exclamation"), tags$strong(" Broad, mainly immune reference types: cell types absent from the reference get the closest label. "),
                    tags$span("Show limitations", style = "text-decoration: underline; color: #8a6d3b;")
                  ),
                  tags$ul(
                    style = "display: none; margin: 8px 0 0 0; padding-left: 18px; line-height: 1.6;",
                    tags$li("The celldex references are bulk RNA-seq or microarray profiles of sorted cells (mouse: 2 references; human: 5)."),
                    tags$li("A cell type absent from the reference is given the closest available label (e.g. skeletal muscle myonuclei are labelled 'Cardiomyocytes', FAPs 'Fibroblasts')."),
                    tags$li("Suitable for blood/immune cells or as a first pass. For other tissues, prefer a reference dataset of the same tissue (previous tab), and always confirm labels with marker genes.")
                  )
                ),
                fluidRow(
                  column(width = 4,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #667eea; min-height: 360px;",
                           tags$h5(tags$span("1", style = "display: inline-block; width: 24px; height: 24px; border-radius: 50%; background: #667eea; color: white; text-align: center; line-height: 24px; font-size: 13px; margin-right: 6px;"),
                                   "Choose the reference", style = "color: #667eea; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           selectInput("singler_species_single", "Species:", choices = c("Mouse" = "mouse", "Human" = "human"), selected = "mouse"),
                           selectInput("singler_reference_single", "Reference dataset:", choices = NULL),
                           tags$p(icon("circle-question"), " Bulk expression profiles of sorted cell types (celldex). Pick the one closest to your tissue; most are immune-centred.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           radioButtons("singler_label_level_single", "Label resolution:",
                                        choices = c("Broad cell types" = "label.main", "Fine subtypes" = "label.fine"), selected = "label.main"),
                           tags$p(icon("circle-question"), " Broad: main types (e.g. 'T cells', 'Fibroblasts'), more robust. Fine: subtypes (e.g. 'T cells (CD4+, naive)'), more precise but more often wrong or ambiguous.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;")
                         )
                  ),
                  column(width = 4,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #11998e; min-height: 360px;",
                           tags$h5(tags$span("2", style = "display: inline-block; width: 24px; height: 24px; border-radius: 50%; background: #11998e; color: white; text-align: center; line-height: 24px; font-size: 13px; margin-right: 6px;"),
                                   "Run SingleR", style = "color: #11998e; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           radioButtons("singler_mode_single", "Annotation mode:",
                                        choices = c("Per cell + cluster vote" = "cell", "Per cluster (faster)" = "cluster"), selected = "cell"),
                           tags$p(icon("circle-question"), " Per cell: every cell is annotated, then each cluster takes its majority label; gives agreement and per-cell confidence. Per cluster: the average profile of each cluster is annotated once; much faster, but agreement is always 100% so mixed clusters cannot be detected.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           actionButton("run_singler_single", tagList(icon("play"), " Run SingleR"), class = "btn-gradient-violet btn-full-width"),
                           uiOutput("singler_status_single")
                         )
                  ),
                  column(width = 4,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #764ba2; min-height: 360px;",
                           tags$h5(tags$span("3", style = "display: inline-block; width: 24px; height: 24px; border-radius: 50%; background: #764ba2; color: white; text-align: center; line-height: 24px; font-size: 13px; margin-right: 6px;"),
                                   "Review & apply", style = "color: #764ba2; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           numericInput("singler_min_agreement_single", "Min. agreement for a confident cluster (%):", value = 50, min = 10, max = 100, step = 5),
                           tags$p(icon("circle-question"), " A cluster is confident when at least this % of its cells share the same label and at most 30% of its cells are unassigned (pruned by SingleR). Changes the summary immediately, no re-run needed.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           checkboxInput("singler_keep_low_single", "Keep current name for low-confidence clusters", value = TRUE),
                           tags$p(icon("circle-question"), " Checked: only confident (green) clusters are renamed. Unchecked: every cluster takes its majority label.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           checkboxInput("singler_unique_single", "Number clusters sharing a label (Neurons_1, Neurons_2...)", value = TRUE),
                           tags$p(icon("circle-question"), " Unchecked: clusters with the same label are merged into one.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           actionButton("apply_singler_labels_single", tagList(icon("tags"), " Apply Labels to Clusters"), class = "btn-gradient-violet btn-full-width"),
                           fluidRow(
                             style = "margin-top: 10px;",
                             column(6, downloadButton("download_singler_clusters_single", "Summary", class = "btn-white-violet", style = "width: 100%;")),
                             column(6, downloadButton("download_singler_cells_single", "Per cell", class = "btn-white-violet", style = "width: 100%;"))
                           )
                         )
                  )
                ),
                tags$div(
                  style = "background-color: white; padding: 15px; border-radius: 6px; margin-top: 15px;",
                  tags$h5(icon("table"), " Cluster summary", style = "color: #495057; margin-top: 0; font-weight: bold;"),
                  tags$p("Agreement: % of cells carrying the majority label. Unassigned: % of cells whose label was pruned as unreliable. Delta: gap between the best and second-best label score. Red rows need checking with marker genes.",
                         style = "font-size: 12px; color: #6c757d;"),
                  DTOutput("singler_summary_single"),
                  fluidRow(
                    style = "margin-top: 15px;",
                    column(width = 6, plotOutput("singler_umap_labels_single", height = "500px")),
                    column(width = 6, plotOutput("singler_umap_delta_single", height = "500px"))
                  )
                )
              )
            )
          ),
          # Display banner for the final UMAP, visible whatever the annotation tab
          tags$div(
            style = "background-color: white; padding: 15px; border-radius: 6px; margin-top: 20px; border-top: 4px solid #11998e;",
            tags$h5(icon("sliders-h"), " Plot Settings", style = "color: #11998e; margin-top: 0; margin-bottom: 10px; font-weight: bold;"),
            fluidRow(
              column(3, selectInput("cluster_col_single", label = "Color by:", choices = NULL)),
              column(3, selectInput("cluster_reduction_single", label = "UMAP reduction:", choices = NULL)),
              column(2, textInput("plot_title", "Plot title:", value = "")),
              column(2, numericInput("label_font_size", "Label size:", value = 5, min = 1, max = 20, step = 0.5)),
              column(2, numericInput("pt_size", "Point size:", value = 1, min = 0.1, max = 3, step = 0.1))
            ),
            tags$hr(style = "margin: 5px 0 10px 0;"),
            tags$h5(icon("palette"), " Cluster Colors", style = "color: #764ba2; margin-top: 0; margin-bottom: 10px; font-weight: bold;"),
            fluidRow(
              column(4, selectInput("cluster_select", "Select cluster:", choices = NULL)),
              column(4, colourInput("cluster_colour", "Choose color:", value = "red")),
              column(4, actionButton("update_colour", tagList(icon("palette"), " Update Color"), class = "btn-gradient-violet btn-full-width", style = "margin-top: 25px;"))
            )
          )
        ),
        tags$div(
          class = "box-gradient-green",
          tags$h4(icon("project-diagram"), " Interactive UMAP", class = "box-title"),
          tags$div(
            class = "plot-container",
            style = "min-height: 600px;",
            plotlyOutput("umap_finale", height = "600px")
          )
        ),
        
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("chart-bar"), " Alternative Visualizations", class = "box-title"),
          fluidRow(
            column(width = 3,
                   selectInput("plot_type_select", "Plot type:", choices = c("FeaturePlot", "VlnPlot", "DotPlot", "RidgePlot"))
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px; margin-top: 15px;",
            plotOutput("selected_plot_display", height = "500px")
          )
        )
      ),
      
      
      ############################## Single/Find markers for a specific cluster ##############################
      tabItem(
        tabName = "cluster_comparison",
        
        tags$div(
          class = "header-violet",
          tags$h2("Cluster Comparison", class = "header-title"),
          tags$p("Identify biomarkers and compare cell populations", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8eaf6; padding: 20px; border-radius: 8px; border-left: 4px solid #667eea; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " Cluster Biomarker Analysis", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Find genes that characterize clusters. Exclusive markers (genes almost absent elsewhere) are in the Exclusive Biomarkers tab; condition comparisons (e.g. WT vs KO) are in the Multiple Datasets module."),
            tags$h5(icon("star"), " Analysis Types:", style = "color: #667eea; margin-top: 15px;"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(strong("Cluster vs all others:"), " markers of one cluster compared with all other cells."),
              tags$li(strong("Cluster vs cluster:"), " genes distinguishing two clusters or groups of clusters.")
            )
          )
        ),
        
        tags$div(
          class = "box-gradient-green",
          tags$h4(icon("project-diagram"), " Cluster Overview", class = "box-title"),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px;",
            plotOutput("umap_plot", height = "500px")
          ),
          fluidRow(
            style = "margin-top: 15px;",
            column(width = 3, selectInput("umap_col_single", label = "Color by:", choices = NULL)),
            column(width = 3, selectInput("umap_reduction_single", label = "UMAP reduction:", choices = NULL)),
            column(width = 3, numericInput("dpi_umap_cluster", "Plot resolution (DPI):", value = 300, min = 72, step = 72)),
            column(width = 3, selectInput("umap_cluster_format", "File format:", choices = c("TIFF" = "tiff", "PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg", "SVG" = "svg"), selected = "tiff"))
          ),
          fluidRow(
            column(width = 3, checkboxInput("show_labels", "Show cluster labels", value = TRUE)),
            column(width = 3, offset = 6, downloadButton("downloadUMAPCluster", "Download UMAP", class = "btn-white-violet", style = "width: 100%;"))
          )
        ),
        
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("dna"), " Differential Expression", class = "box-title"),
          fluidRow(
            column(width = 4, selectInput("assay_de_single", "Assay for DE analyses:", choices = "RNA", selected = "RNA")),
            column(width = 8, tags$p(icon("info-circle"), " log2FC > 0 means higher expression in the first group (target cluster or Group 1). Test thresholds apply when the analysis is run; the adjusted p-value filter only changes the table display.",
                                     style = "color: white; margin-top: 25px; font-size: 13px;"))
          ),
          tags$div(
            class = "de-tabs",
            tabsetPanel(
              id = "de_tabs_single",
              type = "tabs",
              tabPanel(
                title = tagList(icon("dna"), " Cluster vs all others"),
                value = "one_vs_all",
                tags$div(
                  class = "de-tab-help",
                  tags$p(icon("lightbulb"), " Markers of one cluster: genes enriched in the selected cluster compared with all other cells.", style = "margin: 0;")
                ),
                fluidRow(
                  column(width = 4,
                         tags$h5("1. Groups", class = "de-step-title"),
                         selectInput("target_cluster_global_single", "Target cluster:", choices = NULL)
                  ),
                  column(width = 4,
                         tags$h5("2. Test parameters (re-run to apply)", class = "de-step-title"),
                         numericInput("min_pct_single", "Min fraction of cells expressing (0-1):", value = 0.01, min = 0, max = 1, step = 0.01),
                         numericInput("logfc_threshold_single", "Min |log2FC|:", value = 0.1, min = 0, max = 5, step = 0.1),
                         actionButton("find_markers_global_single", tagList(icon("play"), " Run analysis"), class = "btn-white-violet btn-full-width")
                  ),
                  column(width = 4,
                         tags$h5("3. Display & export", class = "de-step-title"),
                         numericInput("pval_adj_filter_global_single", "Max adj. p-value (display):", value = 0.05, min = 0, max = 1, step = 0.01),
                         downloadButton("download_markers_global_single", "Download results", class = "btn-white-violet", style = "width: 100%;")
                  )
                ),
                tags$div(class = "plot-container", style = "margin-top: 15px;", DTOutput("table_global_single"))
              ),
              tabPanel(
                title = tagList(icon("exchange-alt"), " Cluster vs cluster"),
                value = "pairwise",
                tags$div(
                  class = "de-tab-help",
                  tags$p(icon("lightbulb"), " Genes distinguishing two clusters or groups of clusters, e.g. two closely related cell types or states.", style = "margin: 0;")
                ),
                fluidRow(
                  column(width = 4,
                         tags$h5("1. Groups", class = "de-step-title"),
                         selectInput("cluster1_pairwise_single", "Group 1 clusters:", choices = NULL, multiple = TRUE),
                         selectInput("cluster2_pairwise_single", "Group 2 clusters:", choices = NULL, multiple = TRUE)
                  ),
                  column(width = 4,
                         tags$h5("2. Test parameters (re-run to apply)", class = "de-step-title"),
                         numericInput("min_pct_pairwise_single", "Min fraction of cells expressing (0-1):", value = 0.01, min = 0, max = 1, step = 0.01),
                         numericInput("logfc_threshold_pairwise_single", "Min |log2FC|:", value = 0.1, min = 0, max = 5, step = 0.1),
                         actionButton("compare_markers_pairwise_single", tagList(icon("play"), " Run analysis"), class = "btn-white-violet btn-full-width")
                  ),
                  column(width = 4,
                         tags$h5("3. Display & export", class = "de-step-title"),
                         numericInput("pval_adj_filter_pairwise_single", "Max adj. p-value (display):", value = 0.05, min = 0, max = 1, step = 0.01),
                         downloadButton("download_markers_pairwise_single", "Download results", class = "btn-white-violet", style = "width: 100%;")
                  )
                ),
                tags$div(class = "plot-container", style = "margin-top: 15px;", DTOutput("table_pairwise_single"))
              )
            )
          )
        ),
        
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("chart-pie"), " Cluster Composition Analysis", class = "box-title"),
          fluidRow(
            column(6,
                   actionButton("generate_cluster_composition_single", tagList(icon("calculator"), " Generate Cluster Composition"), class = "btn-white-green btn-full-width")
            ),
            column(6,
                   downloadButton("download_cluster_composition_single", "Download Table", class = "btn-white-green", style = "width: 100%;")
            )
          ),
          tags$br(),
          tags$div(
            class = "plot-container",
            DTOutput("cluster_composition_single")
          )
        )
      ),
      ############################## Single/Pathway enrichment ##############################
      tabItem(
        tabName = "enrichment_single",
        
        tags$div(
          class = "header-violet",
          tags$h2("Pathway Enrichment", class = "header-title"),
          tags$p("Find the Gene Ontology terms or custom gene sets enriched between two groups of cells", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8eaf6; padding: 20px; border-radius: 8px; border-left: 4px solid #667eea; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " Enrichment Methods", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Pathway enrichment asks which biological processes or pathways are over-represented among the genes that differ between two groups of cells (for example a knockout and its control within one cluster). A gene set is a list of genes sharing a function, such as a Gene Ontology term (e.g. ", tags$i("muscle contraction, GO:0006936"), ")."),
            tags$h5(icon("list-check"), " 1. Choose the genes", style = "color: #667eea; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(strong("New comparison:"), " Cell-Hub compares group 1 with group 2 (Wilcoxon test on every gene detected in at least 5% of the cells of either group, no fold-change threshold). Both methods below can then run."),
              tags$li(strong("Stored table:"), " a differential expression table saved from another page (cluster markers, condition comparison, pseudobulk test). It contains only the genes that passed its thresholds, so only ORA can run.")
            ),
            tags$h5(icon("chart-bar"), " 2. Two complementary methods", style = "color: #667eea; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(strong("Over-representation analysis (ORA)."), " Takes the significant genes (adjusted p-value and |log2FC| thresholds) and asks, for each gene set, whether more of its genes are in this list than expected by chance among all tested genes (hypergeometric test). Up- and down-regulated genes are tested separately, so each set is reported as enriched among the genes that go up, or that go down, in group 1. ",
                      tags$br(), tags$i("Read: Fold enrichment = how many times more genes of the set than expected (1 = chance); GeneRatio = set genes found / significant genes; BgRatio = set genes / all tested genes; geneID = the significant genes of the set.")),
              tags$li(strong("Gene set enrichment analysis (GSEA)."), " Uses all tested genes ranked from most up to most down in group 1, without any threshold, and asks whether the genes of a set are concentrated at the top or the bottom of the ranking (fgsea, 42 as seed). It detects coordinated moderate changes that ORA misses. ",
                      tags$br(), tags$i("Read: NES > 0 = the set is shifted towards genes up in group 1, NES < 0 = towards genes down in group 1; core_enrichment (leading edge) = the genes that drive the result."))
            ),
            tags$p("In the plots and tables, ", tags$span("red = up in group 1", style = "color: #B2182B; font-weight: bold;"), " and ", tags$span("blue = down in group 1", style = "color: #2166AC; font-weight: bold;"), "; darker = more significant. Tables list the terms below the adjusted p-value cutoff (tick the box to see every tested term, in white when not significant) and are sorted by the column Signed -log10(adj. p) (positive = up, negative = down): the most significant up terms come first and the most significant down terms last. A set can appear both up and down in ORA when some of its genes increase and others decrease; GSEA then gives the dominant direction."),
            tags$h5(icon("database"), " 3. Gene sets", style = "color: #667eea; margin-top: 15px;"),
            tags$p("Gene Ontology (biological process, molecular function or cellular component) is built into the application for mouse, human and rat and works offline; its version is written in the Parameters tab. For another species, or for other collections (MSigDB hallmarks, KEGG, Reactome, your own lists), load a GMT file whose gene symbols match those of the dataset."),
            tags$h5(icon("exclamation-triangle"), " 4. Interpret with care", style = "color: #667eea; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8; color: #555;",
              tags$li("P-values are corrected for the number of gene sets tested (Benjamini-Hochberg); keep terms with adjusted p-value < 0.05."),
              tags$li("Gene Ontology terms overlap: several terms often reflect the same genes. Check the gene column before concluding that several processes change."),
              tags$li("A very high fold enrichment carried by 3 to 5 genes is specific but fragile."),
              tags$li("Cell-level tests treat every cell as an independent replicate, which inflates significance when conditions come from few animals. With several animals or donors per condition, run the pseudobulk test (Multiple datasets, Clusters comparison), store its table and use it here as a stored table.")
            )
          )
        ),
        
        # Settings
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("sliders-h"), " Enrichment Settings", class = "box-title"),
          fluidRow(
            column(4,
                   radioButtons("enrichment_source_single", "Genes to analyse:",
                                choices = c("New comparison between two groups (ORA and GSEA)" = "compare", "Stored differential expression table (ORA only)" = "stored"),
                                selected = "compare"),
                   conditionalPanel(
                     condition = "input.enrichment_source_single == 'compare'",
                     selectInput("enrichment_group_col_single", "Group cells by:", choices = NULL),
                     selectizeInput("enrichment_group1_single", "Group 1:", choices = NULL, multiple = TRUE),
                     selectizeInput("enrichment_group2_single", "Group 2 (empty = all other cells):", choices = NULL, multiple = TRUE),
                     selectizeInput("enrichment_restrict_single", "Restrict to cells of (optional, e.g. one cluster):", choices = NULL, multiple = TRUE),
                     numericInput("enrichment_min_pct_single", "Genes tested: detected in at least this fraction of cells of either group:", value = 0.05, min = 0, max = 1, step = 0.01)
                   ),
                   conditionalPanel(
                     condition = "input.enrichment_source_single == 'stored'",
                     selectInput("enrichment_stored_table_single", "Stored table:", choices = NULL),
                     tags$p(icon("info-circle"), " Tables are stored from the differential expression pages (same list as the Venn diagrams).", style = "font-size: 12px; color: #666;")
                   )
            ),
            column(4,
                   radioButtons("enrichment_gene_sets_single", "Gene sets:",
                                choices = c("Gene Ontology — mouse" = "mouse", "Gene Ontology — human" = "human", "Gene Ontology — rat" = "rat", "My gene sets (GMT file)" = "gmt"),
                                selected = "mouse"),
                   conditionalPanel(condition = "input.enrichment_gene_sets_single != 'gmt'",
                                    selectInput("enrichment_ontology_single", "Ontology:", choices = c("Biological process" = "BP", "Molecular function" = "MF", "Cellular component" = "CC"), selected = "BP")),
                   conditionalPanel(condition = "input.enrichment_gene_sets_single == 'gmt'",
                                    fileInput("enrichment_gmt_file_single", "GMT file:", accept = c(".gmt", ".txt")),
                                    tags$p(icon("info-circle"), " One set per line: name, description, then genes, separated by tabs (MSigDB format). Matching ignores case but does not convert between species.",
                                           style = "font-size: 12px; color: #666;"),
                                    verbatimTextOutput("enrichment_gmt_info_single")),
                   fluidRow(column(6, numericInput("enrichment_min_size_single", "Min. genes per set:", value = 10, min = 2)),
                            column(6, numericInput("enrichment_max_size_single", "Max. genes per set:", value = 500, min = 10)))
            ),
            column(4,
                   checkboxGroupInput("enrichment_methods_single", "Methods:", choices = c("Over-representation (ORA)" = "ORA", "Gene set enrichment (GSEA)" = "GSEA"), selected = c("ORA", "GSEA")),
                   selectInput("enrichment_direction_single", "ORA genes:", choices = c("Up and down (tested separately)" = "both", "Up-regulated in group 1" = "up", "Down-regulated in group 1" = "down"), selected = "both"),
                   fluidRow(column(6, numericInput("enrichment_ora_max_p_adj_single", "Adjusted p-value threshold:", value = 0.05, min = 0, max = 1, step = 0.01)),
                            column(6, numericInput("enrichment_ora_min_log2fc_single", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05))),
                   numericInput("enrichment_show_max_p_adj_single", "Show terms with adjusted p-value below:", value = 0.05, min = 0, max = 1, step = 0.01),
                   checkboxInput("enrichment_show_all_single", "Tables: also show non-significant terms", value = FALSE)
            )
          ),
          fluidRow(
            column(12,
                   actionButton("enrichment_run_single", tagList(icon("play"), " Run Enrichment"), class = "btn-white-violet btn-full-width", style = "margin-top: 10px; margin-bottom: 15px;"),
                   tags$p(icon("box-archive"), " Each run is stored with the other result tables: it is saved with the session and its parameters appear in the report. Delete unwanted runs in Save & Session.",
                          style = "font-size: 12px; color: #666; margin-top: -8px;")
            )
          ),
          tags$div(
            class = "plot-container",
            verbatimTextOutput("enrichment_summary_single")
          )
        ),
        
        # Results
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("diagram-project"), " Enrichment Results", class = "box-title"),
          tags$div(
            class = "de-tabs",
            tabsetPanel(
              tabPanel("Over-representation (ORA)",
                       tags$div(class = "plot-container", style = "min-height: 550px; margin-top: 15px;", plotOutput("enrichment_ora_plot_single", height = "550px")),
                       fluidRow(
                         column(4, numericInput("enrichment_ora_n_terms_single", "Terms shown:", value = 15, min = 1, max = 50)),
                         column(4, downloadButton("download_enrichment_ora_plot_single", "Download Plot (PNG)", class = "btn-white-green", style = "width: 100%; margin-top: 25px;")),
                         column(4, downloadButton("download_enrichment_ora_table_single", "Download Table (Excel)", class = "btn-white-green", style = "width: 100%; margin-top: 25px;"))
                       ),
                       tags$div(class = "plot-container", style = "margin-top: 15px;", DT::DTOutput("enrichment_ora_table_single"))),
              tabPanel("Gene set enrichment (GSEA)",
                       tags$div(class = "plot-container", style = "min-height: 550px; margin-top: 15px;", plotOutput("enrichment_gsea_plot_single", height = "550px")),
                       fluidRow(
                         column(4, numericInput("enrichment_gsea_n_terms_single", "Terms shown:", value = 15, min = 1, max = 50)),
                         column(4, downloadButton("download_enrichment_gsea_plot_single", "Download Plot (PNG)", class = "btn-white-green", style = "width: 100%; margin-top: 25px;")),
                         column(4, downloadButton("download_enrichment_gsea_table_single", "Download Table (Excel)", class = "btn-white-green", style = "width: 100%; margin-top: 25px;"))
                       ),
                       tags$div(class = "plot-container", style = "margin-top: 15px;", DT::DTOutput("enrichment_gsea_table_single"))),
              tabPanel("Parameters",
                       tags$div(class = "plot-container", style = "margin-top: 15px;", tableOutput("enrichment_params_table_single")),
                       downloadButton("download_enrichment_params_single", "Download Parameters (CSV)", class = "btn-white-green", style = "margin-top: 10px;"))
            )
          )
        )
      ),
      ############################## Single/Exclusive Biomarkers ##############################
      tabItem(
        tabName = "exclusive_biomarkers",
        
        tags$div(
          class = "header-violet",
          tags$h2("Exclusive Biomarkers Discovery", class = "header-title"),
          tags$p("Identify genes predominantly expressed in specific clusters", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8eaf6; padding: 20px; border-radius: 8px; border-left: 4px solid #667eea; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " What are Exclusive Biomarkers?", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Exclusive biomarkers are genes that show strong, specific expression in one cluster with minimal expression elsewhere."),
            tags$h5(icon("flask"), " Key Features:", style = "color: #d97706; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(strong("Detection Rate:"), " Percentage of cells expressing the gene in target vs other clusters"),
              tags$li(strong("Expression Magnitude:"), " Mean/median expression levels to ensure biological relevance"),
              tags$li(strong("Fold Change:"), " How much more the gene is expressed in target cluster"),
              tags$li(strong("Statistical Validation:"), " Optional Wilcoxon test for robust marker identification"),
              tags$li(strong("Specificity Score:"), " Combined metric weighing both detection difference and fold change")
            )
          )
        ),
        tags$div(
          class = "box-gradient-green",
          tags$h4(icon("sliders-h"), " Analysis Parameters", class = "box-title"),
          fluidRow(
            column(width = 3,
                   selectInput("exclusive_groupby_col_single",
                               "Group by column:",
                               choices = NULL),
                   tags$small("Metadata column to split cells by (e.g. seurat_clusters, condition, genotype)", style = "color: #666;")
            ),
            column(width = 3,
                   selectInput("exclusive_target_cluster_single",
                               "Target group(s):",
                               choices = NULL,
                               multiple = TRUE),
                   tags$small("Group(s) to find exclusive markers for", style = "color: #666;")
            ),
            column(width = 3,
                   selectInput("exclusive_compare_clusters_single",
                               "Compare against group(s):",
                               choices = NULL,
                               multiple = TRUE),
                   tags$small("Leave empty to compare against all other groups", style = "color: #666;")
            ),
            column(width = 3,
                   numericInput("exclusive_min_pct_target_single", 
                                "Min % in Target:", 
                                value = 50, min = 0, max = 100, step = 5),
                   tags$small("Minimum % of cells expressing gene in target. 50% (default) favours robust markers; use 30% for genes captured at low rates (transcription factors, single-nucleus data, Visium HD bins).", style = "color: #666;")            )
          ),
          fluidRow(
            style = "margin-top: 10px;",
            column(width = 3,
                   numericInput("exclusive_max_pct_other_single", 
                                "Max % in Others:", 
                                value = 25, min = 0, max = 100, step = 5),
                   tags$small("Maximum % of cells expressing gene in other groups", style = "color: #666;")
            ),
            column(width = 3,
                   numericInput("exclusive_min_log2fc_single", 
                                "Min Log2FC:", 
                                value = 1.5, min = 0, max = 5, step = 0.1),
                   tags$small("Minimum fold change (log2 scale)", style = "color: #666;")
            ),
            column(width = 3,
                   numericInput("exclusive_detection_threshold_single", 
                                "Detection Threshold:", 
                                value = 0, min = 0, max = 5, step = 0.1),
                   tags$small("Minimum expression to consider gene as detected", style = "color: #666;")
            ),
            column(width = 3,
                   numericInput("exclusive_min_mean_expr_single", 
                                "Min Mean Expression:", 
                                value = 0.5, min = 0, max = 10, step = 0.1),
                   tags$small("Minimum mean expression in target", style = "color: #666;")
            )
          ),
          fluidRow(
            style = "margin-top: 10px;",
            column(width = 3,
                   selectInput("exclusive_statistical_test_single", 
                               "Statistical Test:", 
                               choices = c("Wilcoxon" = "wilcox", "None" = "none"),
                               selected = "wilcox"),
                   tags$small("Use Wilcoxon test for validation", style = "color: #666;")
            ),
            column(width = 3,
                   numericInput("exclusive_max_pvalue_single", 
                                "Max P-value:", 
                                value = 0.05, min = 0.001, max = 1, step = 0.01),
                   tags$small("Maximum adjusted p-value", style = "color: #666;")
            )
          ),
          fluidRow(
            style = "margin-top: 15px;",
            column(width = 6,
                   actionButton("find_exclusive_markers_single", 
                                tagList(icon("rocket"), " Find Exclusive Biomarkers"), 
                                class = "btn-white-orange btn-full-width")
            ),
            column(width = 6,
                   downloadButton("download_exclusive_markers_single", 
                                  "Download Results", 
                                  class = "btn-white-orange", 
                                  style = "width: 100%;")
            )
          )
        ),
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("table"), " Exclusive Biomarkers Results", class = "box-title"),
          tags$div(
            class = "plot-container",
            DTOutput("table_exclusive_markers_single")
          )
        ),
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("chart-bar"), " Expression Visualization", class = "box-title"),
          tags$div(
            style = "background-color: #d1fae5; padding: 15px; border-radius: 5px; margin-bottom: 15px;",
            tags$p(
              icon("info-circle"), 
              " Enter gene names from the results table above to visualize their expression pattern across clusters.",
              style = "margin: 0; color: #065f46;"
            )
          ),
          fluidRow(
            column(width = 5,
                   textInput("exclusive_genes_to_plot_single", 
                             "Gene Names (separated by , ; or spaces):", 
                             placeholder = "e.g., Gene1, Gene2, Gene3"),
                   tags$br(),
                   selectInput("exclusive_plot_type_single", 
                               "Plot Type:", 
                               choices = c("Dot Plot" = "dotplot", 
                                           "Violin Plot" = "violin",
                                           "Feature Plot" = "feature"),
                               selected = "dotplot")
            ),
            column(width = 4,
                   tags$br(),
                   tags$br(),
                   actionButton("generate_exclusive_plot_single", 
                                tagList(icon("chart-line"), " Generate Plot"), 
                                class = "btn-white-green",
                                style = "width: 100%;")
            ),
            column(width = 3,
                   numericInput("exclusive_plot_dpi_single", "Resolution (DPI):", 
                                value = 300, min = 72, max = 600, step = 72),
                   selectInput("exclusive_plot_format_single", "Format:", 
                               choices = c("TIFF" = "tiff", "PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg", "SVG" = "svg"), 
                               selected = "tiff")
            )
          ),
          fluidRow(
            style = "margin-top: 15px;",
            column(width = 12,
                   downloadButton("download_exclusive_plot_single", 
                                  "Download Plot", 
                                  class = "btn-white-green", 
                                  style = "width: 100%;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px;",
            plotOutput("plot_exclusive_markers_single", height = "600px")
          )
        )
      ),
      
      ############################## Single/DEG Visualizations ##############################
      tabItem(
        tabName = "deg_visualizations",
        
        tags$div(
          class = "header-violet",
          tags$h2("DEG Visualizations", class = "header-title"),
          tags$p("Explore differential expression results through Venn diagrams and volcano plots", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8eaf6; padding: 20px; border-radius: 8px; border-left: 4px solid #667eea; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " DEG Visualization Tools", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Visualize and compare differential expression results from your analyses."),
            tags$h5(icon("chart-bar"), " Available Visualizations:", style = "color: #667eea; margin-top: 15px;"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(strong("Venn Diagrams:"), "Compare gene lists across different comparisons to identify overlaps and unique markers."),
              tags$li(strong("Volcano Plots:"), "Visualize the relationship between fold change and statistical significance.")
            )
          )
        ),
        
        # Venn Diagram Section
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("circle-notch"), " Venn Diagram Comparison", class = "box-title"),
          fluidRow(
            column(4,
                   selectInput("venn_table_1_single", "Select first gene list:", choices = c("None" = ""), selected = ""),
                   checkboxInput("significant_only_venn_1_single", "Include only significant genes", value = TRUE),
                   numericInput("log_fc_threshold_venn_1_single", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05)
            ),
            column(4,
                   selectInput("venn_table_2_single", "Select second gene list:", choices = c("None" = ""), selected = ""),
                   checkboxInput("significant_only_venn_2_single", "Include only significant genes", value = TRUE),
                   numericInput("log_fc_threshold_venn_2_single", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05)
            ),
            column(4,
                   selectInput("venn_table_3_single", "Select third gene list (optional):", choices = c("None" = ""), selected = ""),
                   checkboxInput("significant_only_venn_3_single", "Include only significant genes", value = TRUE),
                   numericInput("log_fc_threshold_venn_3_single", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05)
            )
          ),
          fluidRow(
            column(4,
                   selectInput("venn_direction_single", "Gene selection criteria:", 
                               choices = c("Up-regulated (avg_log2FC > threshold)" = "up", 
                                           "Down-regulated (avg_log2FC < -threshold)" = "down", 
                                           "Both directions (|avg_log2FC| > threshold)" = "both"), 
                               selected = "up")
            ),
            column(4,
                   numericInput("p_val_threshold_venn_single", "p-value threshold:", value = 0.05, min = 0, max = 1, step = 0.01)
            ),
            column(4,
                   checkboxInput("use_adjusted_p_venn_single", "Use adjusted p-values", value = TRUE),
                   actionButton("generate_venn_btn_single", tagList(icon("circle-notch"), " Generate Venn Diagram"), 
                                class = "btn-white-violet btn-full-width", style = "margin-top: 23px;")
            )
          ),
          fluidRow(
            column(4,
                   colourInput("venn_color_1_single", "Color for first set:", value = "#56B4E9")
            ),
            column(4,
                   colourInput("venn_color_2_single", "Color for second set:", value = "#E69F00")
            ),
            column(4,
                   colourInput("venn_color_3_single", "Color for third set (if used):", value = "#009E73")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 400px; margin-top: 15px;",
            plotOutput("venn_plot_single", height = "400px")
          ),
          fluidRow(
            column(6,
                   selectInput("venn_diagram_format_single", "Export format:", 
                               choices = c("PNG" = "png", "TIFF" = "tiff", "PDF" = "pdf", "JPEG" = "jpeg"), 
                               selected = "png"),
                   numericInput("venn_diagram_dpi_single", "DPI:", value = 300, min = 72, step = 72)
            ),
            column(6,
                   downloadButton("download_venn_diagram_single", "Download Venn Diagram", 
                                  class = "btn-white-violet", style = "width: 100%; margin-bottom: 10px;"),
                   downloadButton("download_venn_gene_lists_single", "Download Gene Lists", 
                                  class = "btn-white-violet", style = "width: 100%;")
            )
          ),
          tags$div(
            class = "plot-container",  
            style = "margin-top: 15px;",
            selectInput("selected_gene_set_single", "View gene list:", choices = NULL),
            DTOutput("venn_gene_table_single")
          )
        ),                tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("volcano"), " Volcano Plot", class = "box-title"),
          fluidRow(
            column(4,
                   selectInput("volcano_deg_source_single", "Select DEG results:", 
                               choices = c("None" = ""), selected = ""),
                   numericInput("volcano_log2fc_threshold_single", "Log2FC threshold:", 
                                value = 0.5, min = 0, max = 5, step = 0.1)
            ),
            column(4,
                   numericInput("volcano_pval_threshold_single", "Adjusted p-value threshold:", 
                                value = 0.05, min = 0, max = 1, step = 0.01),
                   numericInput("volcano_label_genes_single", "Number of genes to label:", 
                                value = 10, min = 0, max = 50, step = 1)
            ),
            column(4,
                   numericInput("volcano_point_size_single", "Point size:", 
                                value = 2, min = 0.5, max = 5, step = 0.5),
                   numericInput("volcano_point_alpha_single", "Point transparency:", 
                                value = 0.6, min = 0.1, max = 1, step = 0.1)
            )
          ),
          fluidRow(
            column(4,
                   colourInput("volcano_color_up_single", "Color for upregulated genes:", value = "#00ba38")
            ),
            column(4,
                   colourInput("volcano_color_down_single", "Color for downregulated genes:", value = "#f8766d")
            ),
            column(4,
                   colourInput("volcano_color_ns_single", "Color for non-significant genes:", value = "#619cff")
            )
          ),
          fluidRow(
            column(12,
                   actionButton("generate_volcano_btn_single", 
                                tagList(icon("volcano"), " Generate Volcano Plot"), 
                                class = "btn-white-green btn-full-width", 
                                style = "margin-top: 10px; margin-bottom: 15px;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px; margin-top: 15px;",
            plotOutput("volcano_plot_single", height = "500px")
          ),
          fluidRow(
            column(6,
                   selectInput("volcano_plot_format_single", "Export format:", 
                               choices = c("PNG" = "png", "TIFF" = "tiff", "PDF" = "pdf", "JPEG" = "jpeg", "SVG" = "svg"), 
                               selected = "png"),
                   numericInput("volcano_plot_dpi_single", "DPI:", value = 300, min = 72, step = 72)
            ),
            column(6,
                   downloadButton("download_volcano_plot_single", "Download Volcano Plot", 
                                  class = "btn-white-green", style = "width: 100%; margin-top: 25px;")
            )
          )
        )
      ),
      
      
      
      ############################## Single/Subsetting ##############################
      tabItem(
        tabName = "subset",
        
        tags$div(
          class = "header-violet",
          tags$h2("Data Subsetting", class = "header-title"),
          tags$p("Create focused datasets from specific cell populations", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8eaf6; padding: 20px; border-radius: 8px; border-left: 4px solid #667eea; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " Data Subsetting Tools", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Select specific cell populations either by cluster identity or gene expression patterns. All previous analyses are preserved in the subset."),
            tags$h5(icon("star"), " Subsetting Methods:", style = "color: #667eea; margin-top: 15px;"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(strong("Cluster-based:"), "Select specific clusters to retain"),
              tags$li(strong("Expression-based:"), "Filter cells based on expression of key genes"),
              tags$li(strong("Results:"), "Compare original and subset UMAPs")
            )
          )
        ),
        
        tags$div(
          class = "box-gradient-green",
          tags$h4(icon("project-diagram"), " Original Dataset", class = "box-title"),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px;",
            plotOutput("global_umap", height = "500px")
          )
        ),
        
        fluidRow(
          style = "margin-top: 20px;",
          column(width = 6,
                 tags$div(
                   class = "box-gradient-violet",
                   tags$h4(icon("filter"), " Subset by Clusters", class = "box-title"),
                   tags$div(
                     class = "plot-container",
                     selectInput("select_ident_subset", "Select clusters:", choices = NULL, multiple = TRUE),
                     actionButton("apply_subset", tagList(icon("cut"), " Create Subset"), class = "btn-white-violet btn-full-width")
                   )
                 )
          ),
          column(width = 6,
                 tags$div(
                   class = "box-gradient-green",
                   tags$h4(icon("dna"), " Subset by Expression", class = "box-title"),
                   tags$div(
                     class = "plot-container",
                     numericInput("expression_threshold", "Expression threshold:", value = 0.1),
                     textInput("gene_list_subset", "Genes (separated by , ; or spaces):", value = ""),
                     numericInput("num_genes_to_express", "Minimum expressed genes:", value = 1, min = 1),
                     checkboxInput("negative_gene_subset", "Negative selection (exclude expressing cells)", value = FALSE),
                     actionButton("apply_gene_subset", tagList(icon("cut"), " Create Subset"), class = "btn-white-green btn-full-width")
                   )
                 )
          )
        ),
        
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("eye"), " Subset Preview", class = "box-title"),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px;",
            plotOutput("subset_umap", height = "500px")
          ),
          tags$div(
            style = "margin-top: 15px; text-align: center;",
            downloadButton("download_subset_seurat", "Save Subset as RDS", class = "btn-white-violet", style = "padding: 10px 30px; font-size: 16px;")
          )
        )

      ),
      ############################## Single/Session ##############################
      tabItem(
        tabName = "session_single",
        tags$div(
          class = "header-violet",
          tags$h2("Save & Session", class = "header-title"),
          tags$p("Save the whole analysis in one file, export the analysis parameters, see what the object contains and manage stored result tables", class = "header-subtitle")
        ),
        fluidRow(
          column(5,
                 tags$div(
                   class = "box-gradient-violet",
                   tags$h4(icon("floppy-disk"), " Save analysis", class = "box-title"),
                   tags$div(
                     style = "background: white; border-radius: 6px; padding: 15px; color: #2c3e50;",
                     tags$p("One .rds file containing the Seurat object (clusters, cell type names, colors, metadata, reductions, analysis parameters) and all stored result tables. Plots are not stored: they are regenerated from the object."),
                     tags$p(icon("rotate-left"), " To resume: ", tags$b("Load dataset → processed object → select this file."), style = "font-size: 12px; color: #6c757d;"),
                     downloadButton("download_session_single", tagList(icon("floppy-disk"), " Save analysis (.rds)"),
                                    class = "btn-gradient-violet", style = "width: 100%;")
                   )
                 ),
                 tags$div(
                   class = "box-gradient-green",
                   style = "margin-top: 20px;",
                   tags$h4(icon("file-export"), " Export report", class = "box-title"),
                   tags$div(
                     style = "background: white; border-radius: 6px; padding: 15px; color: #2c3e50;",
                     tags$p("Every analysis parameter actually applied (QC thresholds, normalization, PCs, resolution, seeds, DoubletFinder, annotation), as a CSV table or a PDF report.",
                            style = "font-size: 12px; color: #6c757d;"),
                     fluidRow(
                       column(6, downloadButton("download_params_single", tagList(icon("file-csv"), " Parameters (CSV)"), class = "btn-gradient-green btn-full-width")),
                       column(6, downloadButton("download_report_single", tagList(icon("file-pdf"), " PDF report"), class = "btn-gradient-green btn-full-width"))
                     )
                   )
                 ),
                 tags$div(
                   class = "box-gradient-green",
                   style = "margin-top: 20px;",
                   tags$h4(icon("circle-info"), " Object contents", class = "box-title"),
                   tags$div(
                     style = "background: white; border-radius: 6px; padding: 10px; color: #2c3e50;",
                     tableOutput("session_object_summary_single")
                   )
                 )
          ),
          column(7,
                 tags$div(
                   class = "box-gradient-green",
                   tags$h4(icon("table"), " Stored result tables", class = "box-title"),
                   tags$div(
                     style = "background: white; border-radius: 6px; padding: 15px; color: #2c3e50;",
                     tags$p("Tables created by cluster comparisons and exclusive biomarker searches (used by Venn diagrams and volcano plots). Select rows, then delete them if they are no longer needed.",
                            style = "font-size: 12px; color: #6c757d;"),
                     DT::DTOutput("session_tables_single"),
                     tags$div(
                       style = "margin-top: 10px; display: flex; justify-content: space-between; align-items: center; gap: 10px; flex-wrap: wrap;",
                       tags$span(textOutput("session_size_info_single", inline = TRUE), style = "font-size: 12px; color: #6c757d;"),
                       actionButton("delete_session_tables_single", tagList(icon("trash-alt"), " Delete selected"),
                                    class = "btn-gradient-violet")
                     )
                   )
                 )
          )
        )
      ),
      ############################## Multiple/Loading Data ##############################
      tabItem(
        tabName = "load_datasets_merge",
        tags$div(
          style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
             padding: 30px; margin: -15px -15px 30px -15px; 
             border-radius: 0; color: white; text-align: center;",
          tags$h2("Multiple Dataset Analysis", 
                  style = "margin: 0 0 10px 0; font-weight: bold; font-size: 32px;"),
          tags$p("Integrate and analyze multiple single-cell datasets together", 
                 style = "margin: 0; font-size: 18px; opacity: 0.9;")
        ),
        tags$div(
          style = "background-color: #e8f4f8; padding: 20px; border-radius: 8px; 
             border-left: 4px solid #667eea; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(
              icon("info-circle"), " About Dataset Integration",
              tags$span(icon("chevron-down"), 
                        style = "float: right; font-size: 14px;")
            )
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Integrate multiple single-cell datasets to perform comparative analysis 
              across conditions, tissues, or time points while correcting for batch effects."),
            
            tags$h5(icon("star"), " Key Capabilities:", 
                    style = "color: #667eea; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Mix & Match:"), " Combine raw 10X data with pre-processed Seurat objects"),
              tags$li(tags$strong("Five Methods:"), " Seurat anchors (CCA), Seurat reciprocal PCA (RPCA), Harmony, fastMNN, or Simple Merge (no batch correction)"),
              tags$li(tags$strong("Batch Correction:"), " Automatically correct for technical variation"),
              tags$li(tags$strong("Quality Control:"), " Unified QC parameters across all datasets")
            )
          )
        ),
        tags$h4("Choose Your Loading Method", 
                style = "color: #495057; margin-bottom: 15px; font-weight: bold;"),
        fluidRow(
          column(6,
                 tags$div(
                   style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                 padding: 25px; border-radius: 10px; min-height: 450px;
                 box-shadow: 0 4px 6px rgba(0,0,0,0.1); color: white;
                 display: flex; flex-direction: column;",
                   tags$div(
                     style = "text-align: center; margin-bottom: 20px;",
                     tags$div(
                       style = "background-color: rgba(255,255,255,0.2); 
                     width: 80px; height: 80px; border-radius: 50%; 
                     display: inline-flex; align-items: center; 
                     justify-content: center; margin-bottom: 15px;",
                       icon("layer-group", style = "font-size: 36px;")
                     ),
                     tags$h4("Load Multiple Datasets", 
                             style = "margin: 0; font-weight: bold;")
                   ),
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.1); 
                   padding: 15px; border-radius: 6px; margin-bottom: 20px; flex-grow: 1;",
                     tags$p(tags$strong("Perfect for:"), style = "margin-bottom: 10px;"),
                     tags$ul(
                       style = "margin: 0 0 15px 0; padding-left: 20px; font-size: 14px;",
                       tags$li("Comparing multiple conditions"),
                       tags$li("Time-series experiments"),
                       tags$li("Cross-tissue comparisons"),
                       tags$li("Batch effect correction")
                     ),
                     tags$p(tags$strong("Supported formats:"), style = "margin-bottom: 10px;"),
                     tags$ul(
                       style = "margin: 0 0 15px 0; padding-left: 20px; font-size: 14px;",
                       tags$li("Raw 10X data (ZIP files)"),
                       tags$li("H5 files (10X HDF5 format)"),
                       tags$li("Pre-processed Seurat objects (.rds)"),
                       tags$li("Mix different types in same analysis!")
                     ),
                     tags$div(
                       style = "background-color: rgba(255,255,255,0.15); 
                     padding: 10px; border-radius: 4px; margin-top: 12px;",
                       icon("cogs"), 
                       tags$strong(" Choose integration method after loading"),
                       style = "font-size: 13px;"
                     )
                   ),
                   tags$div(
                     style = "margin-top: auto;",
                     actionButton("open_file_input_modal",
                                  tagList(icon("upload"), " Load Multiple Datasets"),
                                  class = "btn-lg",
                                  style = "width: 100%; background-color: white; 
                               color: #667eea; border: none; 
                               font-weight: bold; padding: 15px;
                               border-radius: 8px; font-size: 16px;
                               box-shadow: 0 2px 4px rgba(0,0,0,0.1);")
                   ),
                   tags$div(
                     style = "text-align: center; margin-top: 15px; 
                   font-size: 12px; opacity: 0.9;",
                     icon("clock"), " Processing time: 5-15 minutes depending on size"
                   )
                 )
          ),
          column(6,
                 tags$div(
                   style = "background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%); 
                 padding: 25px; border-radius: 10px; min-height: 450px;
                 box-shadow: 0 4px 6px rgba(0,0,0,0.1); color: white;
                 display: flex; flex-direction: column;",
                   tags$div(
                     style = "text-align: center; margin-bottom: 20px;",
                     tags$div(
                       style = "background-color: rgba(255,255,255,0.2); 
                     width: 80px; height: 80px; border-radius: 50%; 
                     display: inline-flex; align-items: center; 
                     justify-content: center; margin-bottom: 15px;",
                       icon("database", style = "font-size: 36px;")
                     ),
                     tags$h4("Load Pre-integrated Object", 
                             style = "margin: 0; font-weight: bold;")
                   ),
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.1); 
                   padding: 15px; border-radius: 6px; margin-bottom: 20px; flex-grow: 1;",
                     tags$p(tags$strong("Perfect for:"), style = "margin-bottom: 10px;"),
                     tags$ul(
                       style = "margin: 0 0 15px 0; padding-left: 20px; font-size: 14px;",
                       tags$li("Previously integrated datasets"),
                       tags$li("Continuing analysis from saved checkpoint"),
                       tags$li("Sharing integrated data between users"),
                       tags$li("Skip integration and start analyzing")
                     ),
                     tags$p(tags$strong("Requirements:"), style = "margin-bottom: 10px;"),
                     tags$ul(
                       style = "margin: 0 0 15px 0; padding-left: 20px; font-size: 14px;",
                       tags$li("Seurat object with multiple datasets (.rds)"),
                       tags$li("Already integrated or merged"),
                       tags$li("Contains 'dataset' metadata column"),
                       tags$li("Normalized and scaled data ready")
                     ),
                     tags$div(
                       style = "background-color: rgba(255,255,255,0.15); 
                     padding: 10px; border-radius: 4px; margin-top: 12px;",
                       icon("zap"), 
                       tags$strong(" Instant loading - no integration needed"),
                       style = "font-size: 13px;"
                     )
                   ),
                   tags$div(
                     style = "background-color: white; padding: 15px; 
                   border-radius: 8px; margin-top: auto;",
                     tags$label("Select Pre-integrated Seurat Object:", 
                                style = "font-weight: bold; color: #495057; 
                              margin-bottom: 10px; display: block;"),
                     fileInput("load_seurat_file_merge", NULL,
                               accept = c('.rds', '.h5ad', 'application/x-hdf5'),
                               buttonLabel = tagList(icon("upload"), " Browse RDS/H5AD"),
                               placeholder = "No file selected"),
                     tags$div(
                       style = "font-weight: bold; color: #495057;",
                       radioButtons("optimize_memory_preintegrated", NULL,
                                    choices = list(
                                      "Safe (recommended)" = "slim",
                                      "Aggressive" = "aggressive",
                                      "Keep full object" = "full"
                                    ),
                                    selected = "slim", inline = TRUE)
                     )
                   ),
                   tags$div(
                     style = "text-align: center; margin-top: 15px; 
                   font-size: 12px; opacity: 0.9;",
                     icon("bolt"), " Loads in seconds"
                   )
                 )
          )
        ),
        tags$div(
          style = "margin-top: 20px;",
          tags$div(
            style = "background-color: #1e272e; padding: 15px; border-radius: 8px;
                     border: 2px solid #dee2e6; max-height: 300px; overflow-y: auto;",
            tags$div(
              style = "display: flex; align-items: center; gap: 8px; margin-bottom: 10px;",
              icon("terminal", style = "color: #667eea;"),
              tags$strong("Loading Logs", style = "color: white;")
            ),
            verbatimTextOutput("merge_preintegrated_logs")
          )
        ),
        conditionalPanel(
          condition = "output.datasets_loaded",
          tags$div(
            style = "margin-top: 30px;",
            tags$h4("Integration Methods Available", 
                    style = "color: #495057; margin-bottom: 15px; font-weight: bold;"),
            # One card per method; descriptions match performDataIntegration() in integration_functions.R
            fluidRow(
              column(3,
                     tags$div(
                       style = "background-color: #e3f2fd; padding: 15px; border-radius: 8px; border-left: 4px solid #2196F3; height: 100%;",
                       tags$h6(icon("project-diagram"), strong(" Seurat anchors (CCA)"), style = "color: #1976D2; margin-top: 0;"),
                       tags$p("FindIntegrationAnchors + IntegrateData on 2,000 anchor features (dims 1:30, k parameters adapted to the smallest dataset). Creates an 'integrated' assay. Best for 2-3 datasets; memory-intensive.",
                              style = "font-size: 12px; margin: 0; color: #555;")
                     )
              ),
              column(3,
                     tags$div(
                       style = "background-color: #f3e5f5; padding: 15px; border-radius: 8px; border-left: 4px solid #9c27b0; height: 100%;",
                       tags$h6(icon("wave-square"), strong(" Harmony"), style = "color: #7b1fa2; margin-top: 0;"),
                       tags$p("PCA (50 PCs) on the merged data, then Harmony on the 'dataset' variable (dims 1:30, 20 iterations, seed 42). Corrects the embedding only, expression values are unchanged. Fast and scalable to many datasets.",
                              style = "font-size: 12px; margin: 0; color: #555;")
                     )
              ),
              column(3,
                     tags$div(
                       style = "background-color: #ede7f6; padding: 15px; border-radius: 8px; border-left: 4px solid #673AB7; height: 100%;",
                       tags$h6(icon("link"), strong(" fastMNN"), style = "color: #4527A0; margin-top: 0;"),
                       tags$p("batchelor::fastMNN on log-normalized variable genes (d = 50, k = 20, seed 42). Creates the 'integrated.mnn' reduction. Tolerates cell types present in only some datasets.",
                              style = "font-size: 12px; margin: 0; color: #555;")
                     )
              ),
              column(3,
                     tags$div(
                       style = "background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800; height: 100%;",
                       tags$h6(icon("layer-group"), strong(" Simple Merge"), style = "color: #e65100; margin-top: 0;"),
                       tags$p("Concatenation without batch correction (PCA 30 PCs, UMAP on PCs 1-20). Use for datasets processed together, or to visualize the uncorrected batch effect.",
                              style = "font-size: 12px; margin: 0; color: #555;")
                     )
              )
            )
          )
        ),
      # Demonstration data instructions (PBMC v3 + v2: a technical batch to compare the integration methods)
      example_data_box("multiple")),
    
      
      ############################## Multiple/Metadata Management ##############################
      tabItem(
        tabName = "metadata_management_merge",
        
        # Header
        tags$div(
          style = "background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); 
                       padding: 30px; margin: -15px -15px 30px -15px; 
                       border-radius: 0; color: white; text-align: center;",
          tags$h2("Metadata Management", 
                  style = "margin: 0 0 10px 0; font-weight: bold; font-size: 32px;"),
          tags$p("Add and manage custom metadata columns for your integrated datasets", 
                 style = "margin: 0; font-size: 18px; opacity: 0.9;")
        ),
        
        # Info box
        tags$div(
          style = "background-color: #e8f4f8; padding: 20px; border-radius: 8px; 
                       border-left: 4px solid #f093fb; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(
              icon("info-circle"), " About Metadata Columns",
              tags$span(icon("chevron-down"), 
                        style = "float: right; font-size: 14px;")
            )
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Metadata columns allow you to annotate your datasets with experimental information 
                          that can be used for grouping, coloring, and filtering your analyses."),
            
            tags$h5(icon("lightbulb"), " Common Use Cases:", 
                    style = "color: #f093fb; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Experimental conditions:"), " Control vs Treatment, Young vs Old"),
              tags$li(tags$strong("Time points:"), " Day0, Day3, Day7"),
              tags$li(tags$strong("Biological replicates:"), " Rep1, Rep2, Rep3"),
              tags$li(tags$strong("Technical batches:"), " Batch_A, Batch_B"),
              tags$li(tags$strong("Tissue types:"), " Muscle, Heart, Brain")
            ),
            tags$div(
              style = "background-color: #fff3cd; padding: 12px; border-radius: 6px; 
                         border-left: 3px solid #ffc107; margin-top: 15px;",
              tags$p(
                icon("exclamation-triangle"), 
                tags$strong(" Important:"), 
                " You must load and integrate datasets first before adding metadata columns.",
                style = "margin: 0;"
              )
            )
          )
        ),
        
        # Check if data is loaded
        conditionalPanel(
          condition = "!output.datasets_loaded",
          tags$div(
            style = "background-color: #f8d7da; padding: 30px; border-radius: 8px; 
                         border-left: 4px solid #dc3545; text-align: center;",
            icon("exclamation-circle", style = "font-size: 48px; color: #dc3545; margin-bottom: 15px;"),
            tags$h4("No Datasets Loaded", style = "color: #721c24; margin-bottom: 10px;"),
            tags$p("Please load and integrate your datasets first in the 'Load datasets' tab.", 
                   style = "color: #721c24; margin-bottom: 20px;"),
            
          )
        ),
        
        # Main content - only shown when data is loaded
        conditionalPanel(
          condition = "output.datasets_loaded",
          
          # Display current datasets
          tags$div(
            style = "background-color: #ffffff; padding: 20px; border-radius: 8px; 
                         border: 2px solid #dee2e6; margin-bottom: 25px;",
            tags$h4(
              icon("database"), " Loaded Datasets",
              style = "color: #495057; margin-top: 0; margin-bottom: 15px; font-weight: bold;"
            ),
            uiOutput("loaded_datasets_display_metadata"),
            tags$hr(style = "margin: 20px 0;"),
            tags$h5(
              icon("columns"), " Existing Metadata Columns",
              style = "color: #495057; margin-bottom: 10px; font-weight: bold;"
            ),
            uiOutput("existing_metadata_display")
          ),
          
          # Metadata creation interface
          tags$div(
            style = "background-color: #ffffff; padding: 25px; border-radius: 8px; 
                         border: 2px solid #dee2e6; margin-bottom: 25px;",
            tags$h4(
              icon("edit"), " Create New Metadata Columns",
              style = "color: #495057; margin-top: 0; margin-bottom: 20px; font-weight: bold;"
            ),
            
            # NEW: Split column selector
            uiOutput("metadata_split_selector"),
            
            # Dynamic metadata inputs
            uiOutput("metadata_inputs"),
            
            # Action buttons
            # Action buttons
            tags$div(
              style = "margin-top: 20px; padding-top: 20px; border-top: 2px solid #dee2e6;",
              fluidRow(
                column(6,
                       actionButton(
                         "add_field",
                         tagList(icon("plus-circle"), " Add New Metadata Field"),
                         class = "btn-info btn-lg",
                         style = "width: 100%; font-weight: bold;"
                       )
                ),
                column(6,
                       actionButton(
                         "add_metadata",
                         tagList(icon("check-circle"), " Apply Metadata to Datasets"),
                         class = "btn-success btn-lg",
                         style = "width: 100%; font-weight: bold;"
                       )
                )
              ),
              tags$div(
                style = "text-align: center; margin-top: 15px; color: #6c757d; font-size: 13px;",
                icon("info-circle"), 
                " Fill in field names and values for each group, then click 'Apply' to save changes. Save the updated object from the 'Save & Session' tab."
              )
            )
          )
        )
      ),
      
      ############################## Multiple/Scaling and PCA reduction ##############################
      tabItem(
        tabName = "clustering_merge",
        tags$div(
          class = "header-violet",
          tags$h2("Clustering - Multiple Datasets", class = "header-title"),
          tags$p("Cluster integrated datasets and visualize combined populations", class = "header-subtitle")
        ),
        
        tags$div(
          class = "info-box-turquoise",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Integrated Clustering", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("After integration, clustering identifies shared cell populations across all datasets while accounting for batch effects."),
            tags$h5(icon("star"), " Key Steps:"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Scaling & PCA:"), " Reduce dimensionality and identify main sources of variation across integrated data."),
              tags$li(tags$strong("Find Neighbors:"), " Calculate cell-cell similarities using the integrated reduction."),
              tags$li(tags$strong("Clustering:"), " Group cells into clusters representing distinct populations across all datasets.")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Use the elbow plot to determine optimal number of dimensions for clustering")
            )
          )
        ),
        tags$div(
          class = "box-gradient-green",
          tags$h4(icon("compress"), " Scaling & PCA", class = "box-title"),
          fluidRow(
            column(6,
                   actionButton("runScalePCA", tagList(icon("play"), " Run Scaling, PCA and Elbow Plot"), class = "btn-white-violet btn-full-width")
            ),
            column(6,
                   actionButton("force_cleanup_merge", tagList(icon("broom"), " Free Memory"), class = "btn-white-violet btn-full-width")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 400px; margin-top: 15px;",
            plotOutput("elbow_plot2", height = "400px")
          )
        ),
        
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("sitemap"), " Cluster Cells", class = "box-title"),
          tags$div(
            style = "background-color: rgba(255,255,255,0.9); padding: 20px; border-radius: 8px; margin-bottom: 15px;",
            fluidRow(
              column(width = 3,
                     tags$div(
                       style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #667eea;",
                       tags$label("Number of Dimensions:", style = "font-weight: bold; color: #495057;"),
                       numericInput("dimension_2", NULL, value = 15, min = 1),
                       actionButton("runFindNeighbors", tagList(icon("project-diagram"), " Find Neighbors & UMAP"), class = "btn-gradient-violet btn-full-width"),
                       actionButton("compute_3d_umap_merge", tagList(icon("cube"), " Compute 3D UMAP"), class = "btn-white-violet btn-full-width", style = "margin-top: 8px;"), 
                       tags$p(icon("info-circle"), " Use integrated dimensions", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                     )
              ),
              column(width = 3,
                     tags$div(
                       style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #11998e;",
                       tags$label("Resolution:", style = "font-weight: bold; color: #495057;"),
                       numericInput("resolution_step2", NULL, min = 0.01, step = 0.1, value = 0.5),
                       actionButton("runFindClusters", tagList(icon("sitemap"), " Find Clusters"), class = "btn-gradient-green btn-full-width"),
                       tags$p(icon("info-circle"), " Higher = more clusters", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                     )
              ),
              column(width = 3,
                     tags$div(
                       style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #764ba2;",
                       tags$label("Algorithm:", style = "font-weight: bold; color: #495057;"),
                       selectInput("algorithm_select_merge", NULL, choices = list("Original Louvain" = 1, "Louvain with Multilevel Refinement" = 2, "SLM Algorithm" = 3)),
                       tags$label("Export Format:", style = "font-weight: bold; color: #495057; margin-top: 10px;"),
                       selectInput("umap_merge_format", NULL, choices = c("TIFF" = "tiff", "PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg", "SVG" = "svg"), selected = "tiff"),
                       tags$p(icon("info-circle"), " Clustering method", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                     )
              ),
              column(width = 3,
                     tags$div(
                       style = "background-color: white; padding: 15px; border-radius: 6px; border-left: 4px solid #38ef7d;",
                       tags$label("Resolution (DPI):", style = "font-weight: bold; color: #495057;"),
                       numericInput("dpi_umap_merge", NULL, value = 300, min = 72, max = 1200, step = 72),
                       downloadButton("downloadUMAP_merge", "Download UMAP", class = "btn-white-green", style = "width: 100%; margin-top: 10px;"),
                       tags$p(icon("info-circle"), " Image quality", style = "font-size: 11px; color: #6c757d; margin: 5px 0 0 0;")
                     )
              )
            )
          ),
          tags$div(
            style = "background-color: white; padding: 15px; border-radius: 6px; margin-bottom: 15px; color: #495057;",
            fluidRow(
              column(2, selectInput("cluster_cells_reduction_merge", "Reduction:", choices = NULL)),
              column(4, radioButtons("cluster_label_mode_merge", "Cluster names:",
                                     choices = c("On the UMAP" = "plot", "Side legend" = "legend", "Both" = "both", "None" = "none"),
                                     selected = "both", inline = TRUE)),
              column(2, numericInput("cluster_label_size_merge", "Name size:", value = 5, min = 1, max = 20, step = 0.5)),
              column(2, checkboxInput("cluster_bold_labels_merge", "Bold names", value = FALSE),
                     checkboxInput("dark_mode_umap_cluster_merge", tagList(icon("moon"), " Black background"), value = FALSE)),
              column(2, checkboxInput("remove_axes_umap_merge", "Remove Axes", FALSE))
            ),
            tags$p(icon("info-circle"), " Display options apply instantly and to the downloaded image; no need to re-run clustering.",
                   style = "font-size: 11px; color: #6c757d; margin: 0;")
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 600px;",
            plotOutput("UMAPPlot_cluster_merge", height = "600px")
          )
        ),
        # 3D UMAP panel, shown once "Compute 3D UMAP" has been run
        conditionalPanel(
          condition = "output.has_umap3d_merge == true",
          tags$div(
            style = "background-color: #fff3e0; padding: 15px; border-radius: 8px; border-left: 4px solid #ff9800; margin-top: 20px;",
            tags$h4(icon("cube"), " 3D UMAP", style = "color: #ff9800; margin-top: 0; font-weight: bold;"),
            fluidRow(
              column(2, checkboxInput("umap3d_dark_clust_merge", tagList(icon("moon"), " Dark mode"), value = FALSE)),
              column(2, checkboxInput("umap3d_hide_grid_clust_merge", tagList(icon("border-none"), " Hide grid"), value = FALSE)),
              column(2, checkboxInput("umap3d_hide_axes_clust_merge", tagList(icon("arrows-alt"), " Hide axes"), value = FALSE)),
              column(3, numericInput("umap3d_pt_size_clust_merge", "Point size:", value = 2, min = 0.5, max = 20, step = 0.5)),
              column(3, tags$p(icon("info-circle"), " Click and drag to rotate", style = "font-size: 12px; color: #6c757d; margin-top: 10px;"))
            ),
            plotlyOutput("clustering_umap3d_merge", height = "650px")
          )
        )
      ),
      ############################## Multiple/Visualize genes expressions ##############################
      tabItem(
        tabName = "plots_genes_expressions_merge",
        
        tags$div(
          class = "header-violet",
          tags$h2("Gene Expression Visualization - Multiple Datasets", class = "header-title"),
          tags$p("Explore and compare gene expression across integrated datasets", class = "header-subtitle")
        ),
        fluidRow(
          column(12,
                 tags$div(
                   class = "info-box-turquoise",
                   tags$div(
                     class = "info-toggle",
                     onclick = "$(this).next('.info-content-hidden').slideToggle(); $(this).find('.chevron').toggleClass('fa-chevron-down fa-chevron-up');",
                     tags$h4(
                       icon("question-circle"), 
                       " Understanding Visualization Options",
                       tags$span(class = "chevron", icon("chevron-down"), style = "float: right; font-size: 14px;")
                     )
                   ),
                   tags$div(
                     class = "info-content-hidden",
                     style = "display: none; margin-top: 15px;",
                     
                     tags$p(
                       "This interface allows you to visualize gene expression by grouping cells based on metadata columns from your integrated Seurat object (",
                       tags$code("seurat_object@meta.data"), ").",
                       style = "margin-bottom: 20px;"
                     ),
                     tags$div(
                       style = "display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-bottom: 20px;",
                       tags$div(
                         tags$h5(icon("layer-group"), " Group By - Choose Your Primary Grouping:", style = "color: #00acc1; margin-top: 0;"),
                         tags$p("Select which metadata column defines your main cell groups:", style = "margin-bottom: 10px; font-size: 13px;"),
                         tags$ul(
                           style = "line-height: 1.6; font-size: 13px;",
                           tags$li(
                             tags$strong(tags$code("dataset"), ":"), 
                             " Groups cells by dataset origin",
                             tags$ul(
                               style = "font-size: 12px;",
                               tags$li("Compare expression between different samples"),
                               tags$li("Example: Dataset1 vs Dataset2 vs Dataset3")
                             )
                           ),
                           tags$li(
                             tags$strong(tags$code("seurat_clusters"), ":"), 
                             " Groups cells by Seurat clustering results",
                             tags$ul(
                               style = "font-size: 12px;",
                               tags$li("Cluster IDs from ", tags$code("FindClusters()"), " (0, 1, 2, 3...)"),
                               tags$li("Shows all clusters across all datasets")
                             )
                           ),
                           tags$li(
                             tags$strong(tags$code("orig.ident"), ":"), 
                             " Groups cells by original sample identity",
                             tags$ul(
                               style = "font-size: 12px;",
                               tags$li("Original sample names before integration")
                             )
                           ),
                           tags$li(
                             tags$strong(tags$code("ClusterIdents"), ":"), 
                             " Groups cells by annotated identities",
                             tags$ul(
                               style = "font-size: 12px;",
                               tags$li("Your manual cell type annotations")
                             )
                           ),
                           tags$li(
                             tags$strong("Other columns:"), 
                             " Any metadata column",
                             tags$ul(
                               style = "font-size: 12px;",
                               tags$li(tags$code("species"), ", ", tags$code("dataset_type"), ", ", tags$code("treatment"), ", etc.")
                             )
                           )
                         )
                       ),
                       tags$div(
                         tags$h5(icon("list"), " Select clusters to remove:", style = "color: #00acc1; margin-top: 0;"),
                         tags$p(
                           "This dropdown is dynamically populated based on your selections:",
                           style = "margin-bottom: 10px; font-size: 13px;"
                         ),
                         tags$ul(
                           style = "line-height: 1.6; font-size: 13px;",
                           tags$li(
                             tags$strong("If Group By = dataset:"),
                             tags$ul(
                               style = "font-size: 12px;",
                               tags$li("Without Split by: Shows datasets"),
                               tags$li("With Split by: Shows split column values")
                             )
                           ),
                           tags$li(
                             tags$strong("If Group By = seurat_clusters:"), 
                             " Shows cluster IDs to filter display"
                           ),
                           tags$li(
                             tags$strong("If Group By = ClusterIdents:"), 
                             " Shows cell type names to filter"
                           ),
                           tags$li(
                             tags$strong("If Group By = other:"), 
                             " Shows unique values from that column"
                           )
                         ),
                         
                         tags$div(
                           class = "tip-box",
                           style = "margin-top: 15px;",
                           tags$p(
                             icon("lightbulb"), 
                             tags$strong(" Quick Start:"), 
                             tags$br(),
                             "1. Group By = ", tags$code("dataset"),
                             tags$br(),
                             "2. Add Split by = ", tags$code("seurat_clusters"),
                             tags$br(),
                             "3. Use ", tags$em("Compare selected clusters"),
                             style = "margin: 0; line-height: 1.8; font-size: 12px;"
                           )
                         )
                       )
                     ),
                     tags$h5(icon("sitemap"), " Advanced Options (only when Group By = dataset):", style = "color: #00acc1; margin-top: 20px; margin-bottom: 15px;"),
                     tags$div(
                       style = "display: grid; grid-template-columns: 1fr 1fr; gap: 20px;",
                       tags$div(
                         style = "background-color: #f0f9ff; padding: 15px; border-radius: 6px; border-left: 3px solid #00acc1;",
                         tags$h5(icon("filter"), " Split by - Second Grouping Layer", style = "margin-top: 0; color: #0288d1; font-size: 14px;"),
                         tags$p(
                           "Subdivide each dataset by another metadata column.", 
                           style = "margin: 5px 0 10px 0; font-size: 12px;"
                         ),
                         
                         tags$div(
                           style = "background-color: white; padding: 8px; border-radius: 4px; margin-top: 10px;",
                           tags$strong("Example 1:", style = "font-size: 12px;"),
                           tags$ul(
                             style = "margin: 3px 0 0 0; font-size: 11px;",
                             tags$li("Group By: ", tags$code("dataset")),
                             tags$li("Split by: ", tags$code("seurat_clusters")),
                             tags$li(tags$em("→ Dataset1_0, Dataset1_1, Dataset2_0..."))
                           )
                         ),
                         tags$div(
                           style = "background-color: white; padding: 8px; border-radius: 4px; margin-top: 8px;",
                           tags$strong("Example 2:", style = "font-size: 12px;"),
                           tags$ul(
                             style = "margin: 3px 0 0 0; font-size: 11px;",
                             tags$li("Group By: ", tags$code("dataset")),
                             tags$li("Split by: ", tags$code("ClusterIdents")),
                             tags$li(tags$em("→ Dataset1_Tcells, Dataset1_Bcells..."))
                           )
                         )
                       ),
                       tags$div(
                         style = "background-color: #f0f9ff; padding: 15px; border-radius: 6px; border-left: 3px solid #00acc1;",
                         tags$h5(icon("code-branch"), " Comparison Mode", style = "margin-top: 0; color: #0288d1; font-size: 14px;"),
                         tags$p("Choose how to display the subdivided groups:", style = "margin: 5px 0 10px 0; font-size: 12px;"),
                         tags$div(
                           style = "margin-bottom: 10px;",
                           tags$strong("Split by clusters:", style = "font-size: 12px;"),
                           tags$div(
                             style = "font-size: 11px; margin-top: 3px;",
                             "Show ", tags$em("all combinations"), " of Dataset × Split"
                           ),
                           tags$div(
                             style = "background-color: white; padding: 6px; border-radius: 4px; margin-top: 5px; font-family: monospace; font-size: 10px; color: #555;",
                             "Dataset1_0 | Dataset1_1 | Dataset2_0 | Dataset2_1"
                           )
                         ),
                         tags$div(
                           tags$strong("Compare selected clusters:", style = "font-size: 12px;"),
                           tags$div(
                             style = "font-size: 11px; margin-top: 3px;",
                             "Filter to ", tags$em("specific values"), ", compare across datasets"
                           ),
                           tags$div(
                             style = "background-color: white; padding: 6px; border-radius: 4px; margin-top: 5px; font-size: 10px;",
                             tags$div("1. Select values (e.g., cluster 0, 3)"),
                             tags$div("2. Keep only matching cells", style = "margin-top: 2px;"),
                             tags$div("3. Compare between datasets:", style = "margin-top: 2px;"),
                             tags$div(
                               style = "font-family: monospace; padding-left: 10px; margin-top: 2px; color: #555;",
                               "Dataset1 | Dataset2 | Dataset3"
                             )
                           )
                         )
                       )
                     )
                   )
                 )
          )
        ),
        tags$div(
          class = "box-gradient-green",
          tags$h4(icon("dna"), " Gene Selection & Global Parameters", class = "box-title"),
          fluidRow(
            column(4,
                   selectizeInput("group_by_select",  "Group By (X-axis):",  choices = c("Select grouping..." = ""), options = list( placeholder = 'Select a column', onInitialize = I('function() { this.setValue(""); }')))
                   
            ),
            column(4,
                   selectInput("viz_assay_merge", "Select Assay:", choices = c("RNA", "integrated"), selected = "RNA")
            ),
            column(4,
                   numericInput("dpi_input_merge", "Images resolution for download:", value = 300, min = 72, step = 72)
            )
          ),
          fluidRow(
            column(4,
                   selectizeInput("metadata_to_compare", "Split By (Optional):", choices = c("None" = ""),options = list( placeholder = 'Select to split (optional)', onInitialize = I('function() { this.setValue(""); }'))),
                   uiOutput("split_values_filter_ui_gene_viz")
            ),
            column(4,
                   selectInput("plot_format_merge", "Select Image Format:", choices = c("PNG" = "png", "JPEG" = "jpeg", "TIFF" = "tiff", "SVG" = "svg", "PDF" = "pdf"), selected = "png")
            ),
            column(4,
                   numericInput("axis_text_size_merge", "Axis text size:", value = 12, min = 8, max = 24, step = 1),
                   numericInput("title_text_size_merge", "Title text size:", value = 14, min = 10, max = 28, step = 1),
                   numericInput("axis_line_width_merge", "Line width:", value = 0.5, min = 0.1, max = 2, step = 0.1),
                   checkboxInput("bold_axis_text_merge", "Bold text", value = FALSE),
            )
          ),
          fluidRow(
            column(4,
                   pickerInput("geneInput_merge", "Select Genes:", choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE))
            )
          )
        ),
        
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("map-marked-alt"), " Feature Plot", class = "box-title"),
          textInput("gene_list_feature_merge", "Selected genes for FeaturePlot (separated by , ; or spaces):", value = "", width = "100%"),
          fluidRow(
            column(3,
                   actionButton("show_feature_plot_merge", tagList(icon("eye"), " Run Feature Plot"), class = "btn-white-violet btn-full-width"),
                   checkboxInput("add_nolegend_feature_merge", "Remove Legend", value = FALSE),
                   checkboxInput("add_noaxes_feature_merge", "Remove Axes", value = FALSE)
            ),
            column(3,
                   selectInput("min_cutoff_feature_merge", "Minimum Cutoff:", 
                               choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30", "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80", "q90" = "q90", "q99" = "q99"), 
                               selected = NA),
                   checkboxInput("show_coexpression_merge", "Show Co-expression", value = FALSE),
                   checkboxInput("enable_3d_feature_merge", "Display in 3D", FALSE)
            ),
            column(3,
                   selectInput("max_cutoff_feature_merge", "Maximum Cutoff:", 
                               choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30", "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80", "q90" = "q90", "q99" = "q99"), 
                               selected = NA),
                   checkboxInput("hide_grid_feature_merge", "Hide grid (3D)", FALSE),
                   checkboxInput("dark_mode_feature_plot_merge", "Dark mode", value = FALSE)
            ),
            column(3,
                   downloadButton("downloadFeaturePlotMerge", "Download (2D only)", class = "btn-white-violet", style = "width: 100%;"),
                   helpText("3D download not supported", style = "font-size: 10px; margin-top: 5px; color: #888;")
            )
          ),
          tags$div(
            class = "plot-output-container",
            conditionalPanel(
              condition = "!input.enable_3d_feature_merge",
              plotOutput("FeaturePlot2", height = "auto")
            ),
            conditionalPanel(
              condition = "input.enable_3d_feature_merge",
              plotly::plotlyOutput("FeaturePlot2_3d", height = "600px")
            )                          )
        ),
        
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("braille"), " Dot Plot", class = "box-title"),
          textInput("gene_list_dot_merge", "Selected genes for DotPlot (separated by , ; or spaces):", value = "", width = "100%"),
          fluidRow(
            column(4,
                   actionButton("runDotPlot", tagList(icon("eye"), " Generate DotPlot"), class = "btn-white-violet btn-full-width"),
                   checkboxInput("add_nolegend_dot_merge", "Remove Legend", value = FALSE),
                   checkboxInput("add_noaxes_dot_merge", "Remove Axes", value = FALSE),
                   checkboxInput("invert_axes_merge", "Invert Axes", value = TRUE)
            ),
            column(4,
                   selectizeInput("cluster_order_dotplot_merge", "Select clusters to show:", 
                                  choices = NULL, multiple = TRUE,
                                  options = list(plugins = list('remove_button'))),
                   checkboxInput("dotplot_scale_merge", "Scale expression (z-score across groups)", value = TRUE)
            ),
            column(4,
                   selectInput("color_palette_dotplot", "Color palette:",
                               choices = c("Default (Seurat)" = "default",
                                           "RdYlBu", "Blues", "Reds", "Greens", 
                                           "Spectral", "PuOr", "BrBG"),
                               selected = "default"),
                   downloadButton("downloadDotPlotMerge", "Download", class = "btn-white-violet", style = "width: 100%;"),
                   sliderInput("dot_scale_dotplot", "Dot size scale:", 
                               min = 0.5, max = 15, value = 7, step = 0.5)
            )
          ),
          tags$div(
            class = "plot-output-container",  
            plotOutput("DotPlot2") 
          )
        ),
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("chart-bar"), " Violin Plot", class = "box-title"),
          textInput("gene_list_vln_merge", "Selected genes for VlnPlot (separated by , ; or spaces):", value = "", width = "100%"),
          fluidRow(
            column(4,
                   actionButton("runVlnPlot", tagList(icon("eye"), " Run Vln Plot"), class = "btn-white-green btn-full-width"),
                   checkboxInput("add_nolegend_vln_merge", "Remove Legend", value = FALSE),
                   checkboxInput("add_noaxes_vln_merge", "Remove Axes", value = FALSE)
            ),
            column(4,
                   selectInput("cluster_order_vln_merge", "Select clusters to show:", choices = NULL, multiple = TRUE),
                   checkboxInput("hide_vln_points_merge", "Hide points", value = FALSE)
            ),
            column(4,
                   downloadButton("downloadVlnPlotMerge", "Download", class = "btn-white-green", style = "width: 100%;")
            )
          ),
          tags$div(
            class = "plot-output-container",
            plotOutput("VlnPlot2", height = "auto")
          )
        ),
        
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("mountain"), " Ridge Plot", class = "box-title"),
          textInput("gene_list_ridge_merge", "Selected genes for Ridge Plot (separated by , ; or spaces):", value = "", width = "100%"),
          fluidRow(
            column(4,
                   actionButton("runRidgePlot", tagList(icon("eye"), " Run Ridge Plot"), class = "btn-white-green btn-full-width"),
                   checkboxInput("add_nolegend_ridge_merge", "Remove Legend", value = FALSE),
                   checkboxInput("add_noaxes_ridge_merge", "Remove Axes", value = FALSE)
            ),
            column(4,
                   downloadButton("downloadRidgePlotMerge", "Download", class = "btn-white-green", style = "width: 100%;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px; margin-top: 15px;",
            plotOutput("Ridge_plot_merge", height = "500px")
          )
        ),
        
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("table"), " Gene Expression Analysis", class = "box-title"),
          fluidRow(
            column(4,
                   textInput("gene_list_genes_expression_merge", "Selected genes (separated by , ; or spaces):", value = "")
            ),
            column(4,
                   numericInput("logfc_threshold_genes_expression_merge", "Expression Threshold:", value = 1, min = 0, max = 10)
            ),
            column(4,
                   actionButton("analyze_btn_genes_expression_merge", tagList(icon("calculator"), " Analyze Expression"), class = "btn-white-violet btn-full-width"),
                   tags$br(), tags$br(),
                   downloadButton("download_genes_number_expression_merge", "Download Table", class = "btn-white-violet", style = "width: 100%;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "margin-top: 15px;",
            dataTableOutput("expression_summary_merge")
          )
        )
      ),
      ############################## Multiple/Heatmaps ##############################
      tabItem(
        tabName = "heatmaps_dualexpression_merge",
        
        tags$div(
          class = "header-violet",
          tags$h2("Heatmaps - Multiple Datasets", class = "header-title"),
          tags$p("Advanced expression analysis across integrated datasets", class = "header-subtitle")
        ),
        tags$div(
          style = "background-color: #e8f5e9; padding: 20px; border-radius: 8px; border-left: 4px solid #2c5f2d; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(
              icon("info-circle"), 
              " About Heatmap Visualization", 
              tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;")
            )
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p(
              "Visualize average gene expression across multiple genes and clusters from all integrated datasets.",
              style = "margin-bottom: 20px; font-size: 1.05em;"
            ),
            tags$div(
              style = "padding: 15px; background-color: white; border-radius: 6px;",
              tags$h5(icon("th"), " Heatmap Analysis", style = "color: #2c5f2d; margin-top: 0; margin-bottom: 10px;"),
              tags$p(
                tags$strong("What it shows:"), " Each column represents the mean expression of all cells within a cluster (or cluster-dataset combination). Each row is one gene.",
                style = "margin-bottom: 8px;"
              ),
              tags$ul(
                style = "margin-left: 20px; margin-bottom: 0;",
                tags$li(tags$strong("Z-score scaling:"), " Shows relative expression (high/low) for each gene across clusters. Yellow = high, purple = low."),
                tags$li(tags$strong("Split by dataset:"), " Compare the same clusters across different datasets (e.g., IIa-WT vs IIa-mdx vs IIa-Htz). Each cluster becomes multiple columns, one per dataset."),
                tags$li(tags$strong("Cluster ordering:"), " Drag and reorder clusters in the selection field to group related cell types together.")
              )
            )
          )
        ),
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("th"), " Heatmap Analysis", class = "box-title"),
          fluidRow(
            column(4,
                   tags$h5(icon("cog"), " Analysis Settings", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   selectInput("split_by_heatmap_multi", "Split by (second axis):", choices = NULL),
                   conditionalPanel(
                     condition = "input.split_by_heatmap_multi != 'None' && input.split_by_heatmap_multi != ''",
                     radioButtons("heatmap_order_multi", "Panels:", choices = c("One panel per cluster" = "cluster", "One panel per split value" = "split"), selected = "cluster", inline = TRUE)
                   ),
                   tags$hr(style = "margin: 10px 0;"),
                   radioButtons("heatmap_scale_multi", "Values:", choices = c("Row z-score" = "zscore", "Average expression (log)" = "none"), selected = "zscore", inline = TRUE),
                   selectInput("color_palette_heatmap_multi", "Color Palette:",
                               choices = c("Viridis" = "viridis", "Magma" = "magma", "Inferno" = "inferno",
                                           "Plasma" = "plasma", "Red-Yellow-Blue" = "RdYlBu",
                                           "Red-White-Blue" = "RdBu", "Blue-White-Red" = "BlueRed",
                                           "Yellow-Red" = "YellowRed"),
                               selected = "viridis")
            ),
            column(4,
                   tags$h5(icon("layer-group"), " Selection", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   selectizeInput("clusters_heatmap_multi", 
                                  "Select Clusters (in order):",
                                  choices = NULL, multiple = TRUE,
                                  options = list(plugins = list('remove_button', 'drag_drop'), persist = FALSE)),
                   tags$hr(style = "margin: 10px 0;"),
                   radioButtons("gene_source_heatmap_multi", "Genes:",
                                choices = c("Manual list" = "manual", "Top markers per cluster" = "cluster_markers",
                                            "Top DE genes between two groups (per cluster)" = "condition_de"),
                                selected = "manual"),
                   conditionalPanel(
                     condition = "input.gene_source_heatmap_multi == 'manual'",
                     textInput("gene_select_heatmap_multi", "Select Genes (separated by , ; or spaces):")
                   ),
                   conditionalPanel(
                     condition = "input.gene_source_heatmap_multi == 'condition_de'",
                     selectizeInput("heatmap_group1_multi", "Group 1 (values of the split column):", choices = NULL, multiple = TRUE),
                     selectizeInput("heatmap_group2_multi", "Group 2:", choices = NULL, multiple = TRUE),
                     selectInput("heatmap_de_direction_multi", "Direction:", choices = c("Up and down" = "both", "Higher in Group 1" = "up", "Higher in Group 2" = "down"), selected = "both"),
                     tags$small("Genes ranked by a cell-level Wilcoxon test within each cluster: use it to choose genes to show, and the pseudobulk DE for statistics.", style = "color: #666;")
                   ),
                   conditionalPanel(
                     condition = "input.gene_source_heatmap_multi != 'manual'",
                     numericInput("n_top_genes_multi", "Number of genes per cluster:", value = 10, min = 1, max = 50, step = 1),
                     downloadButton("download_heatmap_markers_multi", "Export genes CSV", class = "btn-sm btn-info")
                   ),
                   tags$hr(style = "margin: 10px 0;"),
                   actionButton("generateCustomHeatmap_multi", tagList(icon("fire"), " Generate Heatmap"), 
                                class = "btn-white-green", style = "width: 100%;")
            ),
            column(4,
                   tags$h5(icon("download"), " Export Options", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   numericInput("dpi_heatmap_multi", "Resolution (DPI):", value = 300, min = 72, step = 72),
                   selectInput("format_heatmap_multi", "Format:", 
                               choices = c("PNG" = "png", "PDF" = "pdf", "TIFF" = "tiff", "JPEG" = "jpeg"), selected = "png"),
                   downloadButton("download_heatmap_multi", "Download Heatmap", class = "btn-white-green", style = "width: 100%;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "margin-top: 15px; overflow-y: auto; max-height: 1000px; border: 1px solid #ddd;",
            plotOutput("heatmap_plot_multi")
          )
        )
      ),
      
      ############################################## Co-expression Analysis - Multiple Datasets ##############################################################
      tabItem(
        tabName = "coexpression_merge",
        tags$div(
          class = "header-violet",
          tags$h2("Co-expression Analysis - Multiple Datasets", class = "header-title"),
          tags$p("Gene co-expression statistics and expression-based differential analysis", class = "header-subtitle")
        ),
        # --- Scatter plot ---
        tags$div(
          class = "box-gradient-violet",
          tags$h4(icon("braille"), " Gene Co-expression Scatter", class = "box-title"),
          fluidRow(
            column(3,
                   tags$h5(icon("dna"), " Gene Selection", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   textInput("feature1_select_multi", "First Gene:", ""),
                   textInput("feature2_select_multi", "Second Gene:", ""),
                   tags$hr(style = "margin: 10px 0;"),
                   selectInput("assay_select_scatter_multi", "Data Assay:", choices = c("RNA", "integrated"), selected = "RNA")
            ),
            column(3,
                   tags$h5(icon("palette"), " Display", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   selectInput("color_by_scatter_multi", "Color by:", choices = c("Cluster" = "cluster"), selected = "cluster"),
                   selectInput("split_by_scatter_multi", "Split by (optional):", choices = c("None" = "None"), selected = "None")
            ),
            column(3,
                   tags$h5(icon("filter"), " Filters", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   textInput("scatter_text_clusters_multi", "Clusters (optional):", placeholder = "e.g., IIa, IIb"),
                   tags$hr(style = "margin: 10px 0;"),
                   numericInput("threshold_gene1_multi", "Gene 1 Threshold:", value = 0, min = 0, step = 0.1),
                   numericInput("threshold_gene2_multi", "Gene 2 Threshold:", value = 0, min = 0, step = 0.1)
            ),
            column(3,
                   tags$h5(icon("download"), " Export", style = "margin-top: 0; color: #495057; font-weight: bold;"),
                   selectInput("format_scatter_multi", "Format:", choices = c("PNG" = "png", "JPEG" = "jpeg", "TIFF" = "tiff", "SVG" = "svg", "PDF" = "pdf"), selected = "png"),
                   numericInput("dpi_scatter_multi", "Resolution (DPI):", value = 300, min = 72, step = 72),
                   downloadButton("download_scatter_multi", "Download Plot", class = "btn-white-green", style = "width: 100%;")
            )
          ),
          fluidRow(
            column(12,
                   tags$hr(style = "margin: 15px 0;"),
                   actionButton("generateScatter_multi", tagList(icon("braille"), " Generate Plot"), class = "btn-white-green", style = "width: 100%; height: 45px; font-size: 16px; font-weight: bold;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px; margin-top: 15px;",
            plotOutput("scatter_plot_multi", height = "500px")
          )
        ),
        # --- Co-expression stats + DE ---
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("project-diagram"), " Gene Co-expression Analysis", class = "box-title"),
          fluidRow(
            column(width = 4,
                   textAreaInput("gene_text_coexpression_multiple", "Enter genes (separated by , ; or spaces):", value = "", placeholder = "e.g., Pax7, Myod1, Myog", height = "80px", width = "100%"),
                   selectInput("coexpr_assay_multiple", "Assay (normalized thresholds and DE):", choices = NULL),
                   selectInput("coexpr_plot_type_multiple", "Visualization:", choices = c("Stacked Bar" = "bar", "Heatmap" = "heat"), selected = "bar")
            ),
            column(width = 4,
                   radioButtons("coexpr_threshold_unit_multiple", "Thresholds on:", choices = c("Raw UMI counts (robust to ambient RNA)" = "counts", "Normalized expression (log)" = "data"), selected = "counts"),
                   uiOutput("coexpr_thresholds_ui_multiple")
            ),
            column(width = 4,
                   selectInput("coexpr_background_clusters_multiple", "Background cluster(s): cells that should not express these genes", choices = NULL, multiple = TRUE),
                   actionButton("estimate_coexpr_thresholds_multiple", tagList(icon("magic"), " Estimate thresholds from background"), class = "btn-white-violet btn-full-width"),
                   tags$small("Sets each threshold just above the 95th percentile of the gene in the background cells (ambient RNA level); raw UMIs are read from the RNA assay", style = "color: #666;"),
                   uiOutput("coexpr_background_info_multiple")
            )
          ),
          tags$hr(),
          tags$div(
            style = "background-color: white; padding: 15px; border-radius: 8px;",
            tabsetPanel(
              id = "coexpr_main_tabs_multiple", type = "pills",
              tabPanel(
                tagList(icon("chart-bar"), " Co-expression Statistics"), value = "coexpr_stats",
                tags$div(
                  style = "background-color: #f5f5f5; padding: 12px 15px; border-radius: 6px; margin: 15px 0; border-left: 3px solid #11998e;",
                  tags$p(tags$b("What this does:"), " For each cluster (or dataset), counts how many cells express each gene individually and how many co-express all genes simultaneously. Reports exact percentages per group.",
                         style = "margin: 0; font-size: 13px;")
                ),
                fluidRow(
                  column(width = 4,
                         selectInput("coexpr_group_by_multiple", "Primary grouping:", choices = c("By Dataset" = "dataset", "By Cluster" = "cluster"), selected = "dataset"),
                         conditionalPanel(condition = "input.coexpr_group_by_multiple == 'cluster'",
                                          selectInput("coexpr_secondary_split_multiple", "Show dataset breakdown:", choices = c("No" = "none", "Yes" = "dataset"), selected = "none")),
                         conditionalPanel(condition = "input.coexpr_group_by_multiple == 'dataset'",
                                          selectInput("coexpr_secondary_split_multiple_alt", "Show cluster breakdown:", choices = c("No" = "none", "Yes" = "cluster"), selected = "none"))
                  ),
                  column(width = 4, selectInput("coexpr_groupby_col_merge", "Group by column:", choices = NULL),
                         tags$small("Metadata column to restrict cells by", style = "color: #666;")),
                  column(width = 4, selectInput("coexpr_filter_groups_merge", "Restrict to group(s):", choices = NULL, multiple = TRUE),
                         tags$small("Leave empty to use all cells", style = "color: #666;"))
                ),
                fluidRow(
                  column(width = 4, actionButton("analyze_coexpression_multiple", tagList(icon("calculator"), " Analyze Co-expression"), class = "btn-white-green btn-full-width", style = "margin-bottom: 10px;")),
                  column(width = 4, downloadButton("download_coexpression_table_multiple", "Download Table", class = "btn-white-green", style = "width: 100%;")),
                  column(width = 4, downloadButton("download_coexpression_plot_multiple", "Download Plot", class = "btn-white-green", style = "width: 100%;"))
                ),
                tags$hr(),
                tabsetPanel(
                  tabPanel("Results Table", DTOutput("gene_coexpression_table_multiple")),
                  tabPanel("Visualization", plotOutput("gene_coexpression_plot_multiple", height = "700px")),
                  tabPanel("Summary", tableOutput("coexpression_summary_stats_multiple"))                )
              ),
              tabPanel(
                tagList(icon("dna"), " DE: Expressors vs Rest"), value = "coexpr_de",
                tags$div(
                  style = "background-color: #f5f5f5; padding: 12px 15px; border-radius: 6px; margin: 15px 0; border-left: 3px solid #11998e;",
                  tags$p(tags$b("What this does:"), " Within the selected cluster(s), cells are split into two groups based on your gene(s) and thresholds, then differential expression (FindMarkers) identifies what distinguishes the two populations.",
                         style = "margin: 0 0 6px 0; font-size: 13px;"),
                  tags$p(tags$b("1 gene:"), " Cells expressing the gene (> threshold) vs cells that do not.",
                         style = "margin: 0 0 3px 0; font-size: 12px; color: #555;"),
                  tags$p(tags$b("2+ genes:"), " Cells co-expressing ALL genes (> thresholds) vs the rest (partial expressors + negatives).",
                         style = "margin: 0; font-size: 12px; color: #555;")
                ),
                fluidRow(
                  column(4,
                         selectInput("coexpr_de_cluster_multiple", "Cluster(s) to analyze:", choices = NULL, multiple = TRUE),
                         tags$p("Select one or more clusters. Cells are pooled if multiple clusters are selected.",
                                style = "color: #888; font-size: 11px; margin-top: -6px;"),
                         tags$br(),
                         actionButton("run_coexpr_de_multiple", tagList(icon("dna"), " Run DE"), class = "btn-white-green btn-full-width"),
                         tags$br(), tags$br(),
                         uiOutput("coexpr_de_info_multiple"),
                         tags$br(),
                         downloadButton("download_coexpr_de_multiple", "Download DE table", class = "btn-white-green", style = "width: 100%;")
                  ),
                  column(8, DTOutput("coexpr_de_table_multiple"))
                )
              )
            )
          )
        )),
      ############################## assigning cell type identity merge ##############################
      tabItem(
        tabName = "assigning_cell_type_identity_merge",
        
        tags$div(
          class = "header-violet",
          tags$h2("Cell Type Assignment", class = "header-title"),
          tags$p("Annotate clusters with biological identities", class = "header-subtitle")
        ),
        tags$div(
          style = "background-color: #e8f5e9; padding: 20px; border-radius: 8px; border-left: 4px solid #11998e; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Cell Type Assignment", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Give each cluster a biological identity. Three routes are available and can be combined: automated label transfer from your own annotated reference, automated annotation against built-in references (SingleR), and manual renaming based on marker genes. Automated labels are suggestions to be confirmed, not a final annotation."),
            tags$h5(icon("route"), " Annotation methods", style = "color: #11998e; margin-top: 15px;"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(strong("Your reference dataset (recommended):"), " upload an annotated Seurat object (.rds) and choose the column holding the cell types. Seurat anchors (FindTransferAnchors + TransferData) project the reference labels onto each cell and give a per-cell prediction score. Use a reference of the same species, tissue and technology (single-nucleus vs single-cell)."),
              tags$li(strong("Built-in references (SingleR):"), " each cell is correlated with bulk expression profiles of sorted cell types from celldex (mouse: 2 references, human: 5). Labels are broad and mostly immune: suitable for blood/immune cells or as a first pass."),
              tags$li(strong("Manual renaming:"), " in the Manual renaming tab, select a cluster and type its name. Typing an existing name merges the two clusters.")
            ),
            tags$h5(icon("list-check"), " Workflow for automated annotation", style = "color: #11998e; margin-top: 15px;"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(strong("Run:"), " each cell receives a label; each cluster (current cluster column) receives the majority label of its cells. Nothing is renamed at this step."),
              tags$li(strong("Review:"), " the cluster summary gives the agreement (% of cells carrying the majority label), the unassigned fraction (cells below the score threshold or pruned), the median score/delta and a High/Low confidence flag. Red rows must be checked with marker genes."),
              tags$li(strong("Apply:"), " 'Apply Labels to Clusters' renames the clusters (ClusterIdents). Low-confidence clusters can keep their current name, and clusters sharing a label can be numbered (Myonuclei_1, Myonuclei_2)."),
              tags$li(strong("Correct:"), " rename any cluster manually afterwards, or click 'Undo Last Change' (Manual renaming tab) to restore the previous names.")
            ),
            tags$p(icon("database"), " Per-cell predictions are stored as metadata columns (RefTransfer_* and SingleR_*: label, pruned_label, score, delta, cluster_label). They can be used as grouping variables in the other tabs and downloaded as CSV (Summary / Per cell).",
                   style = "font-size: 13px; margin-top: 10px;"),
            tags$p(icon("layer-group"), tags$strong(" Integrated data:"), " annotation uses the RNA expression of each cell, not the integrated embedding; the integrated clusters are only used for the per-cluster vote. Check that a label holds in every dataset (Cluster Composition tab) before applying it.",
                   style = "font-size: 13px; margin-top: 6px;"),
            tags$div(
              style = "background-color: #ffebee; padding: 12px; border-radius: 6px; border-left: 3px solid #e53935; margin-top: 15px;",
              tags$p(icon("triangle-exclamation"), tags$strong(" Limitations of automated annotation"), style = "margin: 0 0 6px 0;"),
              tags$ul(
                style = "margin: 0; padding-left: 18px; line-height: 1.6;",
                tags$li("A method can only return labels present in its reference. A population absent from the reference receives the closest label, often with a high score (e.g. skeletal muscle myonuclei labelled 'Cardiomyocytes' and FAPs 'Fibroblasts' with MouseRNAseqData)."),
                tags$li("A high prediction score means that the reference neighbours agree, not that the cell type is correct."),
                tags$li("Rare, transitional or disease-specific states are usually merged into the nearest reference type.")
              )
            ),
            tags$h5(icon("palette"), " Display", style = "color: #11998e; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(strong("Plot Settings:"), " choose the colouring column, reduction, title, label size and point size of the final UMAP."),
              tags$li(strong("Cluster Colors:"), " assign a colour to each cluster; colours are stored in the object and reused across plots.")
            ),
            tags$div(
              style = "background-color: #fff9c4; padding: 12px; border-radius: 6px; border-left: 3px solid #fbc02d; margin-top: 15px;",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " check each label with known markers: display plots from the gene visualization section by selecting the plot type in the Alternative Visualizations box.", style = "margin: 0;")
            )
          )
        ),
        tags$div(
          style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; border: 2px solid #dee2e6; margin-bottom: 20px;",
          tags$h4(icon("cogs"), " Cluster Identity Assignment", style = "color: #495057; margin-top: 0; margin-bottom: 8px; font-weight: bold;"),
          tags$p("Name clusters manually or from an automated annotation. Automated predictions are first stored as new metadata columns; cluster names only change when you click 'Apply Labels'.",
                 style = "font-size: 13px; color: #6c757d; margin-bottom: 15px;"),
          tabsetPanel(
            id = "annotation_method_tabs_merge", selected = "manual",
            ############ Tab 1: manual renaming ############
            tabPanel(
              tagList(icon("pen"), " Manual renaming"), value = "manual",
              tags$div(
                style = "margin-top: 15px;",
                fluidRow(
                  column(width = 6,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #667eea;",
                           tags$h5("Rename Clusters", style = "color: #667eea; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           fluidRow(
                             column(7,
                                    selectInput("select_cluster_merge", "Select cluster:", choices = NULL),
                                    textInput("rename_single_cluster_merge", "New name:")
                             ),
                             column(5,
                                    actionButton("rename_single_cluster_merge_button", tagList(icon("check"), " Apply Name"), class = "btn-gradient-violet btn-full-width", style = "margin-top: 25px;"),
                                    actionButton("undo_cluster_merge", tagList(icon("undo"), " Undo"), class = "btn-white-violet btn-full-width", style = "margin-top: 40px;")
                             )
                           ),
                           tags$p(icon("info-circle"), " Undo also reverts the last 'Apply Labels' of the reference or SingleR tabs.", style = "font-size: 11px; color: #6c757d; margin: 0;")
                         )
                  ),
                  column(width = 6,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #11998e; font-size: 13px; color: #555;",
                           tags$h5(icon("lightbulb"), " How it works", style = "color: #11998e; margin-top: 0; margin-bottom: 10px; font-weight: bold;"),
                           tags$ul(
                             style = "padding-left: 18px; line-height: 1.7; margin: 0;",
                             tags$li("Select a cluster, type its name, click 'Apply Name'."),
                             tags$li("Typing a name that already exists merges the two clusters."),
                             tags$li("Names are written to the ClusterIdents column and used by all other tabs.")
                           )
                         )
                  )
                ))
            ),
            ############ Tab 2: label transfer from the user's reference ############
            tabPanel(
              tagList(icon("database"), " Your reference dataset"), value = "transfer",
              tags$div(
                style = "margin-top: 15px;",
                tags$div(
                  style = "background-color: #fff3cd; padding: 10px 15px; border-radius: 6px; border-left: 4px solid #f0ad4e; margin-bottom: 15px; font-size: 13px;",
                  tags$div(
                    style = "cursor: pointer;", onclick = "$(this).next().slideToggle();",
                    icon("triangle-exclamation"), tags$strong(" Labels can only come from your reference: a population missing from it gets the closest reference label. "),
                    tags$span("Show limitations", style = "text-decoration: underline; color: #8a6d3b;")
                  ),
                  tags$ul(
                    style = "display: none; margin: 8px 0 0 0; padding-left: 18px; line-height: 1.6;",
                    tags$li("A missing population is often assigned with a high prediction score: the score measures agreement between neighbouring reference cells, not how well the cell matches the reference."),
                    tags$li("Use a reference of the same species, tissue and technology (single-nucleus vs single-cell data differ, e.g. myonuclei are absent from single-cell atlases)."),
                    tags$li("Check every cluster with marker genes before applying labels.")
                  )
                ),
                fluidRow(
                  column(width = 4,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #667eea; min-height: 360px;",
                           tags$h5(tags$span("1", style = "display: inline-block; width: 24px; height: 24px; border-radius: 50%; background: #667eea; color: white; text-align: center; line-height: 24px; font-size: 13px; margin-right: 6px;"),
                                   "Load the reference", style = "color: #667eea; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           fileInput("annotation_reference_file_merge", "Annotated Seurat object (.rds):", accept = ".rds"),
                           selectInput("annotation_reference_label_merge", "Column containing the cell types:", choices = NULL),
                           tags$p(icon("circle-question"), " The metadata column whose values (cell types) will be transferred to your cells.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           uiOutput("annotation_reference_info_merge")
                         )
                  ),
                  column(width = 4,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #11998e; min-height: 360px;",
                           tags$h5(tags$span("2", style = "display: inline-block; width: 24px; height: 24px; border-radius: 50%; background: #11998e; color: white; text-align: center; line-height: 24px; font-size: 13px; margin-right: 6px;"),
                                   "Transfer the labels", style = "color: #11998e; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           numericInput("transfer_min_score_merge", "Min. prediction score per cell (0-1):", value = 0.5, min = 0, max = 1, step = 0.05),
                           tags$p(icon("circle-question"), " Each cell gets the reference label with the highest score (0-1). Cells below this value are marked 'Unassigned'. Higher = stricter; lower it (e.g. 0.3) if many cells end up unassigned with a reference that has many close subtypes. Re-run after changing it.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           numericInput("transfer_max_cells_merge", "Max. reference cells per label (0 = all):", value = 500, min = 0, max = 10000, step = 100),
                           tags$p(icon("circle-question"), " Each reference cell type is randomly reduced to this number of cells (seed 42): faster, and abundant types do not outweigh rare ones. 0 keeps every cell (slower, more memory).", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           numericInput("transfer_dims_merge", "PCA dimensions:", value = 30, min = 5, max = 100, step = 5),
                           tags$p(icon("circle-question"), " Number of principal components used to match your cells to the reference. 30 suits most data; 40-50 for references with many fine subtypes.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           actionButton("run_transfer_merge", tagList(icon("play"), " Transfer Labels"), class = "btn-gradient-violet btn-full-width"),
                           uiOutput("transfer_status_merge")
                         )
                  ),
                  column(width = 4,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #764ba2; min-height: 360px;",
                           tags$h5(tags$span("3", style = "display: inline-block; width: 24px; height: 24px; border-radius: 50%; background: #764ba2; color: white; text-align: center; line-height: 24px; font-size: 13px; margin-right: 6px;"),
                                   "Review & apply", style = "color: #764ba2; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           numericInput("transfer_min_agreement_merge", "Min. agreement for a confident cluster (%):", value = 50, min = 10, max = 100, step = 5),
                           tags$p(icon("circle-question"), " A cluster is confident when at least this % of its cells share the same label and at most 30% of its cells are unassigned. Changes the summary immediately, no re-run needed.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           checkboxInput("transfer_keep_low_merge", "Keep current name for low-confidence clusters", value = TRUE),
                           tags$p(icon("circle-question"), " Checked: only confident (green) clusters are renamed. Unchecked: every cluster takes its majority label.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           checkboxInput("transfer_unique_merge", "Number clusters sharing a label (Myonuclei_1, Myonuclei_2...)", value = TRUE),
                           tags$p(icon("circle-question"), " Unchecked: clusters with the same label are merged into one.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           actionButton("apply_transfer_labels_merge", tagList(icon("tags"), " Apply Labels to Clusters"), class = "btn-gradient-violet btn-full-width"),
                           fluidRow(
                             style = "margin-top: 10px;",
                             column(6, downloadButton("download_transfer_clusters_merge", "Summary", class = "btn-white-violet", style = "width: 100%;")),
                             column(6, downloadButton("download_transfer_cells_merge", "Per cell", class = "btn-white-violet", style = "width: 100%;"))
                           )
                         )
                  )
                ),
                tags$div(
                  style = "background-color: white; padding: 15px; border-radius: 6px; margin-top: 15px;",
                  tags$h5(icon("table"), " Cluster summary", style = "color: #495057; margin-top: 0; font-weight: bold;"),
                  tags$p("Agreement: % of cells carrying the majority label. Unassigned: % of cells below the minimum score. Score: median prediction score. Delta: gap between the best and second-best label. Red rows need checking with marker genes.",
                         style = "font-size: 12px; color: #6c757d;"),
                  DTOutput("transfer_summary_merge"),
                  fluidRow(
                    style = "margin-top: 15px;",
                    column(width = 6, plotOutput("transfer_umap_labels_merge", height = "500px")),
                    column(width = 6, plotOutput("transfer_umap_score_merge", height = "500px"))
                  )
                )
              )
            ),
            ############ Tab 3: SingleR with built-in celldex references ############
            tabPanel(
              tagList(icon("book"), " Built-in references (SingleR)"), value = "singler",
              tags$div(
                style = "margin-top: 15px;",
                tags$div(
                  style = "background-color: #fff3cd; padding: 10px 15px; border-radius: 6px; border-left: 4px solid #f0ad4e; margin-bottom: 15px; font-size: 13px;",
                  tags$div(
                    style = "cursor: pointer;", onclick = "$(this).next().slideToggle();",
                    icon("triangle-exclamation"), tags$strong(" Broad, mainly immune reference types: cell types absent from the reference get the closest label. "),
                    tags$span("Show limitations", style = "text-decoration: underline; color: #8a6d3b;")
                  ),
                  tags$ul(
                    style = "display: none; margin: 8px 0 0 0; padding-left: 18px; line-height: 1.6;",
                    tags$li("The celldex references are bulk RNA-seq or microarray profiles of sorted cells (mouse: 2 references; human: 5)."),
                    tags$li("A cell type absent from the reference is given the closest available label (e.g. skeletal muscle myonuclei are labelled 'Cardiomyocytes', FAPs 'Fibroblasts')."),
                    tags$li("Suitable for blood/immune cells or as a first pass. For other tissues, prefer a reference dataset of the same tissue (previous tab), and always confirm labels with marker genes.")
                  )
                ),
                fluidRow(
                  column(width = 4,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #667eea; min-height: 360px;",
                           tags$h5(tags$span("1", style = "display: inline-block; width: 24px; height: 24px; border-radius: 50%; background: #667eea; color: white; text-align: center; line-height: 24px; font-size: 13px; margin-right: 6px;"),
                                   "Choose the reference", style = "color: #667eea; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           selectInput("singler_species_merge", "Species:", choices = c("Mouse" = "mouse", "Human" = "human"), selected = "mouse"),
                           selectInput("singler_reference_merge", "Reference dataset:", choices = NULL),
                           tags$p(icon("circle-question"), " Bulk expression profiles of sorted cell types (celldex). Pick the one closest to your tissue; most are immune-centred.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           radioButtons("singler_label_level_merge", "Label resolution:",
                                        choices = c("Broad cell types" = "label.main", "Fine subtypes" = "label.fine"), selected = "label.main"),
                           tags$p(icon("circle-question"), " Broad: main types (e.g. 'T cells', 'Fibroblasts'), more robust. Fine: subtypes (e.g. 'T cells (CD4+, naive)'), more precise but more often wrong or ambiguous.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;")
                         )
                  ),
                  column(width = 4,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #11998e; min-height: 360px;",
                           tags$h5(tags$span("2", style = "display: inline-block; width: 24px; height: 24px; border-radius: 50%; background: #11998e; color: white; text-align: center; line-height: 24px; font-size: 13px; margin-right: 6px;"),
                                   "Run SingleR", style = "color: #11998e; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           radioButtons("singler_mode_merge", "Annotation mode:",
                                        choices = c("Per cell + cluster vote" = "cell", "Per cluster (faster)" = "cluster"), selected = "cell"),
                           tags$p(icon("circle-question"), " Per cell: every cell is annotated, then each cluster takes its majority label; gives agreement and per-cell confidence. Per cluster: the average profile of each cluster is annotated once; much faster, but agreement is always 100% so mixed clusters cannot be detected.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           actionButton("run_singler_merge", tagList(icon("play"), " Run SingleR"), class = "btn-gradient-violet btn-full-width"),
                           uiOutput("singler_status_merge")
                         )
                  ),
                  column(width = 4,
                         tags$div(
                           style = "background-color: white; padding: 15px; border-radius: 6px; border-top: 4px solid #764ba2; min-height: 360px;",
                           tags$h5(tags$span("3", style = "display: inline-block; width: 24px; height: 24px; border-radius: 50%; background: #764ba2; color: white; text-align: center; line-height: 24px; font-size: 13px; margin-right: 6px;"),
                                   "Review & apply", style = "color: #764ba2; margin-top: 0; margin-bottom: 15px; font-weight: bold;"),
                           numericInput("singler_min_agreement_merge", "Min. agreement for a confident cluster (%):", value = 50, min = 10, max = 100, step = 5),
                           tags$p(icon("circle-question"), " A cluster is confident when at least this % of its cells share the same label and at most 30% of its cells are unassigned (pruned by SingleR). Changes the summary immediately, no re-run needed.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           checkboxInput("singler_keep_low_merge", "Keep current name for low-confidence clusters", value = TRUE),
                           tags$p(icon("circle-question"), " Checked: only confident (green) clusters are renamed. Unchecked: every cluster takes its majority label.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           checkboxInput("singler_unique_merge", "Number clusters sharing a label (Neurons_1, Neurons_2...)", value = TRUE),
                           tags$p(icon("circle-question"), " Unchecked: clusters with the same label are merged into one.", style = "font-size: 11px; color: #6c757d; margin: -10px 0 12px 0; line-height: 1.4;"),
                           actionButton("apply_singler_labels_merge", tagList(icon("tags"), " Apply Labels to Clusters"), class = "btn-gradient-violet btn-full-width"),
                           fluidRow(
                             style = "margin-top: 10px;",
                             column(6, downloadButton("download_singler_clusters_merge", "Summary", class = "btn-white-violet", style = "width: 100%;")),
                             column(6, downloadButton("download_singler_cells_merge", "Per cell", class = "btn-white-violet", style = "width: 100%;"))
                           )
                         )
                  )
                ),
                tags$div(
                  style = "background-color: white; padding: 15px; border-radius: 6px; margin-top: 15px;",
                  tags$h5(icon("table"), " Cluster summary", style = "color: #495057; margin-top: 0; font-weight: bold;"),
                  tags$p("Agreement: % of cells carrying the majority label. Unassigned: % of cells whose label was pruned as unreliable. Delta: gap between the best and second-best label score. Red rows need checking with marker genes.",
                         style = "font-size: 12px; color: #6c757d;"),
                  DTOutput("singler_summary_merge"),
                  fluidRow(
                    style = "margin-top: 15px;",
                    column(width = 6, plotOutput("singler_umap_labels_merge", height = "500px")),
                    column(width = 6, plotOutput("singler_umap_delta_merge", height = "500px"))
                  )
                )
              )
            )
          ),
          # Display banner for the final UMAP, visible whatever the annotation tab
          tags$div(
            style = "background-color: white; padding: 15px; border-radius: 6px; margin-top: 20px; border-top: 4px solid #11998e;",
            tags$h5(icon("sliders-h"), " Plot Settings", style = "color: #11998e; margin-top: 0; margin-bottom: 10px; font-weight: bold;"),
            fluidRow(
              column(3, selectInput("umap_finale_col_merge", label = "Color by:", choices = NULL)),
              column(3, selectInput("umap_reduction_merge", label = "UMAP reduction:", choices = NULL)),
              column(2, textInput("plot_title_merge", "Plot title:", value = "")),
              column(2, numericInput("label_font_size_merge", "Label size:", value = 5, min = 1, max = 20, step = 0.5)),
              column(2, numericInput("pt_size_merge", "Point size:", value = 1.5, min = 0.1, max = 3, step = 0.1))
            ),
            tags$hr(style = "margin: 5px 0 10px 0;"),
            tags$h5(icon("palette"), " Cluster Colors", style = "color: #764ba2; margin-top: 0; margin-bottom: 10px; font-weight: bold;"),
            fluidRow(
              column(4, selectInput("select_color_merge", "Select cluster:", choices = NULL)),
              column(4, colourInput("select_cluster_merge_color", "Choose color:", value = "red")),
              column(4, actionButton("update_colour_merge_button", tagList(icon("palette"), " Update Color"), class = "btn-gradient-violet btn-full-width", style = "margin-top: 25px;"))
            )
          )
        ),
        tags$div(
          class = "box-gradient-green",
          tags$h4(icon("project-diagram"), " Interactive UMAP", class = "box-title"),
          tags$div(
            class = "plot-container",
            style = "min-height: 600px;",
            plotlyOutput("umap_finale_merge", height = "600px")
          )
        ),
        
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("chart-bar"), " Alternative Visualizations", class = "box-title"),
          fluidRow(
            column(width = 3,
                   selectInput("plot_type_select_merge", "Plot type:", choices = c("FeaturePlot", "VlnPlot", "DotPlot", "RidgePlot"))
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px; margin-top: 15px;",
            plotOutput("selected_plot_display_merge", height = "500px")
          )
        )
      ),
      ############################## Multiple/Calculation of differentially expressed genes ##############################
      tabItem(
        tabName = "cluster_comparison_merge",
        tags$div(
          class = "header-violet",
          tags$h2("Cluster Comparison - Multiple Datasets", class = "header-title"),
          tags$p("Compare clusters and identify biomarkers across integrated datasets", class = "header-subtitle")
        ),
        tags$div(
          style = "background-color: #e8eaf6; padding: 20px; border-radius: 8px; border-left: 4px solid #667eea; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Multi-Dataset Cluster Comparison", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Perform comprehensive comparisons across clusters and datasets to identify conserved and dataset-specific markers."),
            tags$h5(icon("star"), " Comparison Types:", style = "color: #667eea; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Cluster vs all others:"), " markers of one cluster"),
              tags$li(tags$strong("Cluster vs cluster:"), " genes distinguishing two clusters"),
              tags$li(tags$strong("Condition vs condition:"), " two groups of a metadata column (e.g. WT vs KO) within the selected clusters, with a cell-level (Wilcoxon) or sample-level (pseudobulk DESeq2) test")
            )
          )
        ),
        tags$div(
          class = "box-gradient-green",
          tags$h4(icon("project-diagram"), " UMAP Visualization", class = "box-title"),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px;",
            plotOutput("filtered_umap_plot", height = "500px")
          ),
          fluidRow(
            style = "margin-top: 15px;",
            column(width = 4,
                   uiOutput("dataset_filter_ui"),
                   checkboxInput("bold_labels_merge", "Bold labels", value = FALSE),
                   checkboxInput("show_labels_merge", "Show labels", value = TRUE),
                   checkboxInput("dark_mode_filtered_umap", "Dark mode", value = FALSE),
                   selectInput("filtered_umap_col_merge", label = "Color by:", choices = NULL)
                   
            ),
            column(width = 4,
                   selectInput("filtered_umap_reduction_merge", label = "UMAP reduction:", choices = NULL),
                   numericInput("filtered_umap_plot_dpi", "Resolution (DPI):", value = 300, min = 72, step = 72)
            ),
            column(width = 4,
                   selectInput("filtered_umap_format", "File Format:", choices = c("TIFF" = "tiff", "PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg", "SVG" = "svg"), selected = "tiff"),
                   downloadButton("download_filtered_umap_plot", "Download UMAP", class = "btn-white-green", style = "width: 100%;")
            )
          )
        ),
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("dna"), " Differential Expression", class = "box-title"),
          fluidRow(
            column(width = 4, selectInput("assay_de_merge", "Assay for DE analyses:", choices = "RNA")),
            column(width = 8, tags$p(icon("info-circle"), " log2FC > 0 means higher expression in the first group (target cluster or Group 1). Test thresholds apply when the analysis is run; the adjusted p-value filter only changes the table display.",
                                     style = "color: white; margin-top: 25px; font-size: 13px;"))
          ),
          tags$div(
            class = "de-tabs",
            tabsetPanel(
              id = "de_tabs_merge",
              type = "tabs",
              tabPanel(
                title = tagList(icon("dna"), " Cluster vs all others"),
                value = "one_vs_all",
                tags$div(
                  class = "de-tab-help",
                  tags$p(icon("lightbulb"), " Markers of one cluster: genes enriched in the selected cluster compared with all other cells.", style = "margin: 0;")
                ),
                fluidRow(
                  column(width = 4,
                         tags$h5("1. Groups", class = "de-step-title"),
                         selectInput("selected_cluster", label = "Target cluster:", choices = NULL)
                  ),
                  column(width = 4,
                         tags$h5("2. Test parameters (re-run to apply)", class = "de-step-title"),
                         numericInput("min_pct_merge", "Min fraction of cells expressing (0-1):", value = 0.01, min = 0, max = 1, step = 0.01),
                         numericInput("logfc_threshold_merge", label = "Min |log2FC|:", value = 0.25),
                         actionButton("calculate_DE", tagList(icon("play"), " Run analysis"), class = "btn-white-violet btn-full-width")
                  ),
                  column(width = 4,
                         tags$h5("3. Display & export", class = "de-step-title"),
                         numericInput("pval_adj_filter_single_cluster_merge", "Max adj. p-value (display):", value = 0.05, min = 0, max = 1, step = 0.01),
                         downloadButton("download_markers_single_cluster_merge", "Download results", class = "btn-white-violet", style = "width: 100%;")
                  )
                ),
                tags$div(class = "plot-container", style = "margin-top: 15px;", DTOutput("DE_genes_table"))
              ),
              tabPanel(
                title = tagList(icon("exchange-alt"), " Cluster vs cluster"),
                value = "pairwise",
                tags$div(
                  class = "de-tab-help",
                  tags$p(icon("lightbulb"), " Genes distinguishing two clusters, e.g. two closely related cell types or states.", style = "margin: 0;")
                ),
                fluidRow(
                  column(width = 4,
                         tags$h5("1. Groups", class = "de-step-title"),
                         uiOutput("cluster1_compare_ui"),
                         uiOutput("cluster2_compare_ui")
                  ),
                  column(width = 4,
                         tags$h5("2. Test parameters (re-run to apply)", class = "de-step-title"),
                         numericInput("min_pct_compare_merge", "Min fraction of cells expressing (0-1):", value = 0.01, min = 0, max = 1, step = 0.01),
                         numericInput("logfc_threshold_compare_merge", label = "Min |log2FC|:", value = 0.25),
                         actionButton("compare_clusters_button", tagList(icon("play"), " Run analysis"), class = "btn-white-violet btn-full-width")
                  ),
                  column(width = 4,
                         tags$h5("3. Display & export", class = "de-step-title"),
                         numericInput("pval_adj_filter_multiple_clusters_merge", "Max adj. p-value (display):", value = 0.05, min = 0, max = 1, step = 0.01),
                         downloadButton("download_markers_multiple_clusters_merge", "Download results", class = "btn-white-violet", style = "width: 100%;")
                  )
                ),
                tags$div(class = "plot-container", style = "margin-top: 15px;", DTOutput("diff_genes_table_compare"))
              ),
              tabPanel(
                title = tagList(icon("layer-group"), " Condition vs condition"),
                value = "conditions",
                tags$div(
                  class = "de-tab-help",
                  tags$p(icon("lightbulb"), " Compare two groups of a metadata column (WT vs KO, Day0 vs Day7...) within the selected clusters.", style = "margin: 0;"),
                  tags$p(icon("balance-scale"), " Wilcoxon treats each cell as a replicate and is descriptive when cells come from few animals. Pseudobulk sums the raw counts of each sample (animal) and tests the samples with DESeq2: it requires at least 2 samples per group (3 or more recommended).",
                         style = "margin: 6px 0 0 0; font-size: 13px;")
                ),
                fluidRow(
                  column(width = 4,
                         tags$h5("1. Groups", class = "de-step-title"),
                         uiOutput("metadata_column_compare_ui"),
                         uiOutput("dataset1_compare_ui"),
                         uiOutput("dataset2_compare_ui"),
                         uiOutput("cluster_column_selector_ui"),
                         uiOutput("cluster_compare_ui"),
                         checkboxInput("all_clusters", "Use all clusters", value = FALSE)
                  ),
                  column(width = 4,
                         tags$h5("2. Test parameters (re-run to apply)", class = "de-step-title"),
                         radioButtons("de_method_compare", "Statistical test:", choices = c("Wilcoxon (cell level)" = "wilcox", "Pseudobulk DESeq2 (sample level)" = "pseudobulk"), selected = "wilcox"),
                         numericInput("min_pct_compare_dataset_merge", "Min fraction of cells expressing (0-1):", value = 0.01, min = 0, max = 1, step = 0.01),
                         numericInput("logfc_threshold_datasets", label = "Min |log2FC|:", value = 0.25),
                         actionButton("compare_datasets_button", tagList(icon("play"), " Run analysis"), class = "btn-white-violet btn-full-width")
                  ),
                  column(width = 4,
                         tags$h5("3. Display & export", class = "de-step-title"),
                         numericInput("pval_adj_filter_datasets_merge", "Max adj. p-value (display):", value = 0.05, min = 0, max = 1, step = 0.01),
                         downloadButton("download_diff_dataset_cluster", "Download results", class = "btn-white-violet", style = "width: 100%;")
                  )
                ),
                conditionalPanel(
                  condition = "input.de_method_compare == 'pseudobulk'",
                  tags$div(
                    class = "plot-container",
                    style = "margin-top: 10px;",
                    tags$h5(icon("users"), " Pseudobulk design", style = "color: #764ba2; margin-top: 0; font-weight: bold;"),
                    fluidRow(
                      column(width = 6, uiOutput("sample_column_compare_ui")),
                      column(width = 6, numericInput("pb_min_cells_compare", "Minimum cells per sample (in the selected clusters):", value = 10, min = 1, step = 1))
                    ),
                    uiOutput("pseudobulk_design_status"),
                    tableOutput("pseudobulk_design_compare")
                  )
                ),
                tags$div(class = "plot-container", style = "margin-top: 15px;", DTOutput("diff_dataset_cluster"))
              )
            )
          )
        ),
        
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("table"), " Cluster Composition by Metadata", class = "box-title"),
          
          # Info box explaining functionality
          tags$div(
            style = "background-color: #e8f5e9; padding: 12px; border-radius: 6px; margin-bottom: 15px; border-left: 3px solid #11998e;",
            tags$p(
              icon("info-circle"),
              " Analyze how clusters are distributed across any metadata variable (e.g., dataset, condition, timepoint, tissue)",
              style = "margin: 0; font-size: 13px;"
            )
          ),
          
          fluidRow(
            column(4,
                   uiOutput("metadata_column_composition_ui")
            ),
            column(4,
                   actionButton("generate_cluster_table_multiple", 
                                tagList(icon("calculator"), " Generate Table"), 
                                class = "btn-white-green btn-full-width")
            ),
            column(4,
                   downloadButton("download_cluster_composition_multiple", 
                                  "Download Table", 
                                  class = "btn-white-green", 
                                  style = "width: 100%;")
            )
          ),
          tags$br(),
          tags$div(
            class = "plot-container",
            DTOutput("cluster_table_multiple")
          )
        ),
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("chart-pie"), " Cluster Composition - Pie Chart Visualization", class = "box-title"),
          
          # Info explaining what the pie chart shows
          tags$div(
            style = "background-color: #f3e5f5; padding: 12px; border-radius: 6px; margin-bottom: 15px; border-left: 3px solid #764ba2;",
            tags$p(
              icon("info-circle"),
              " Visualize the composition of a selected cluster across metadata groups. ",
              tags$strong("Normalization"), " adjusts for different group sizes to show enrichment patterns.",
              style = "margin: 0; font-size: 13px;"
            )
          ),
          
          fluidRow(
            column(3, 
                   uiOutput("metadata_column_pie_ui")
            ),
            column(3, 
                   selectInput("cluster_pie_select_merge", "Select Cluster:", choices = NULL, width = "100%")
            ),
            column(2, 
                   checkboxInput("normalize_by_dataset_size_merge", "Normalize by group size", value = FALSE)
            ),
            column(2, 
                   checkboxInput("dark_mode_pie_chart_merge", "Dark mode", value = FALSE), 
                   checkboxInput("remove_legend_pie_merge", "Remove Legend", value = FALSE)
            ),
            column(2, 
                   selectInput("pie_format_merge", "Format:", 
                               choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg", "JPEG" = "jpeg", "TIFF" = "tiff"), 
                               selected = "png", width = "100%")
            )
          ),
          fluidRow(
            column(12,
                   tags$br(),
                   downloadButton("download_pie_chart_merge", "Download Chart", style = "width: 100%;", class = "btn-white-violet")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 600px; margin-top: 15px;",
            plotOutput("cluster_composition_pie_merge", height = "600px")
          )
        )
      ),
      
      ############################# Multiple dataset analysis / Pathway enrichment #############################
      tabItem(
        tabName = "enrichment_merge",
        
        tags$div(
          class = "header-violet",
          tags$h2("Pathway Enrichment - Multiple Datasets", class = "header-title"),
          tags$p("Find the Gene Ontology terms or custom gene sets enriched between two groups of cells", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8eaf6; padding: 20px; border-radius: 8px; border-left: 4px solid #667eea; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " Enrichment Methods", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Pathway enrichment asks which biological processes or pathways are over-represented among the genes that differ between two groups of cells (for example a knockout and its control within one cluster). A gene set is a list of genes sharing a function, such as a Gene Ontology term (e.g. ", tags$i("muscle contraction, GO:0006936"), ")."),
            tags$h5(icon("list-check"), " 1. Choose the genes", style = "color: #667eea; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(strong("New comparison:"), " Cell-Hub compares group 1 with group 2 (Wilcoxon test on every gene detected in at least 5% of the cells of either group, no fold-change threshold). Both methods below can then run."),
              tags$li(strong("Stored table:"), " a differential expression table saved from another page (cluster markers, condition comparison, pseudobulk test). It contains only the genes that passed its thresholds, so only ORA can run.")
            ),
            tags$h5(icon("chart-bar"), " 2. Two complementary methods", style = "color: #667eea; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(strong("Over-representation analysis (ORA)."), " Takes the significant genes (adjusted p-value and |log2FC| thresholds) and asks, for each gene set, whether more of its genes are in this list than expected by chance among all tested genes (hypergeometric test). Up- and down-regulated genes are tested separately, so each set is reported as enriched among the genes that go up, or that go down, in group 1. ",
                      tags$br(), tags$i("Read: Fold enrichment = how many times more genes of the set than expected (1 = chance); GeneRatio = set genes found / significant genes; BgRatio = set genes / all tested genes; geneID = the significant genes of the set.")),
              tags$li(strong("Gene set enrichment analysis (GSEA)."), " Uses all tested genes ranked from most up to most down in group 1, without any threshold, and asks whether the genes of a set are concentrated at the top or the bottom of the ranking (fgsea, 42 as seed). It detects coordinated moderate changes that ORA misses. ",
                      tags$br(), tags$i("Read: NES > 0 = the set is shifted towards genes up in group 1, NES < 0 = towards genes down in group 1; core_enrichment (leading edge) = the genes that drive the result."))
            ),
            tags$p("In the plots and tables, ", tags$span("red = up in group 1", style = "color: #B2182B; font-weight: bold;"), " and ", tags$span("blue = down in group 1", style = "color: #2166AC; font-weight: bold;"), "; darker = more significant. Tables list the terms below the adjusted p-value cutoff (tick the box to see every tested term, in white when not significant) and are sorted by the column Signed -log10(adj. p) (positive = up, negative = down): the most significant up terms come first and the most significant down terms last. A set can appear both up and down in ORA when some of its genes increase and others decrease; GSEA then gives the dominant direction."),
            tags$h5(icon("database"), " 3. Gene sets", style = "color: #667eea; margin-top: 15px;"),
            tags$p("Gene Ontology (biological process, molecular function or cellular component) is built into the application for mouse, human and rat and works offline; its version is written in the Parameters tab. For another species, or for other collections (MSigDB hallmarks, KEGG, Reactome, your own lists), load a GMT file whose gene symbols match those of the dataset."),
            tags$h5(icon("exclamation-triangle"), " 4. Interpret with care", style = "color: #667eea; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8; color: #555;",
              tags$li("P-values are corrected for the number of gene sets tested (Benjamini-Hochberg); keep terms with adjusted p-value < 0.05."),
              tags$li("Gene Ontology terms overlap: several terms often reflect the same genes. Check the gene column before concluding that several processes change."),
              tags$li("A very high fold enrichment carried by 3 to 5 genes is specific but fragile."),
              tags$li("Cell-level tests treat every cell as an independent replicate, which inflates significance when conditions come from few animals. With several animals or donors per condition, run the pseudobulk test (Multiple datasets, Clusters comparison), store its table and use it here as a stored table.")
            )
          )
        ),
        
        # Settings
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("sliders-h"), " Enrichment Settings", class = "box-title"),
          fluidRow(
            column(4,
                   radioButtons("enrichment_source_merge", "Genes to analyse:",
                                choices = c("New comparison between two groups (ORA and GSEA)" = "compare", "Stored differential expression table (ORA only)" = "stored"),
                                selected = "compare"),
                   conditionalPanel(
                     condition = "input.enrichment_source_merge == 'compare'",
                     selectInput("enrichment_group_col_merge", "Group cells by:", choices = NULL),
                     selectizeInput("enrichment_group1_merge", "Group 1:", choices = NULL, multiple = TRUE),
                     selectizeInput("enrichment_group2_merge", "Group 2 (empty = all other cells):", choices = NULL, multiple = TRUE),
                     selectizeInput("enrichment_restrict_merge", "Restrict to cells of (optional, e.g. one cluster):", choices = NULL, multiple = TRUE),
                     radioButtons("enrichment_de_method_merge", "Differential expression:",
                                  choices = c("Wilcoxon, cell level (descriptive between conditions)" = "wilcox",
                                              "Pseudobulk DESeq2, sample level (conditions with biological replicates)" = "pseudobulk"),
                                  selected = "wilcox"),
                     conditionalPanel(
                       condition = "input.enrichment_de_method_merge == 'wilcox'",
                       numericInput("enrichment_min_pct_merge", "Genes tested: detected in at least this fraction of cells of either group:", value = 0.05, min = 0, max = 1, step = 0.01)
                     ),
                     conditionalPanel(
                       condition = "input.enrichment_de_method_merge == 'pseudobulk'",
                       selectInput("enrichment_sample_col_merge", "Sample column (one value = one animal / biological replicate):", choices = NULL),
                       numericInput("enrichment_pb_min_cells_merge", "Minimum cells per sample:", value = 10, min = 1, step = 1),
                       tags$p(icon("info-circle"), " Raw counts are summed per sample and tested with DESeq2 (at least 2 samples per group); GSEA ranks the genes by the Wald statistic.",
                              style = "font-size: 12px; color: #666;")
                     )
                   ),
                   conditionalPanel(
                     condition = "input.enrichment_source_merge == 'stored'",
                     selectInput("enrichment_stored_table_merge", "Stored table:", choices = NULL),
                     tags$p(icon("info-circle"), " Tables are stored from the differential expression pages (same list as the Venn diagrams).", style = "font-size: 12px; color: #666;")
                   )
            ),
            column(4,
                   radioButtons("enrichment_gene_sets_merge", "Gene sets:",
                                choices = c("Gene Ontology — mouse" = "mouse", "Gene Ontology — human" = "human", "Gene Ontology — rat" = "rat", "My gene sets (GMT file)" = "gmt"),
                                selected = "mouse"),
                   conditionalPanel(condition = "input.enrichment_gene_sets_merge != 'gmt'",
                                    selectInput("enrichment_ontology_merge", "Ontology:", choices = c("Biological process" = "BP", "Molecular function" = "MF", "Cellular component" = "CC"), selected = "BP")),
                   conditionalPanel(condition = "input.enrichment_gene_sets_merge == 'gmt'",
                                    fileInput("enrichment_gmt_file_merge", "GMT file:", accept = c(".gmt", ".txt")),
                                    tags$p(icon("info-circle"), " One set per line: name, description, then genes, separated by tabs (MSigDB format). Matching ignores case but does not convert between species.",
                                           style = "font-size: 12px; color: #666;"),
                                    verbatimTextOutput("enrichment_gmt_info_merge")),
                   fluidRow(column(6, numericInput("enrichment_min_size_merge", "Min. genes per set:", value = 10, min = 2)),
                            column(6, numericInput("enrichment_max_size_merge", "Max. genes per set:", value = 500, min = 10)))
            ),
            column(4,
                   checkboxGroupInput("enrichment_methods_merge", "Methods:", choices = c("Over-representation (ORA)" = "ORA", "Gene set enrichment (GSEA)" = "GSEA"), selected = c("ORA", "GSEA")),
                   selectInput("enrichment_direction_merge", "ORA genes:", choices = c("Up and down (tested separately)" = "both", "Up-regulated in group 1" = "up", "Down-regulated in group 1" = "down"), selected = "both"),
                   fluidRow(column(6, numericInput("enrichment_ora_max_p_adj_merge", "Adjusted p-value threshold:", value = 0.05, min = 0, max = 1, step = 0.01)),
                            column(6, numericInput("enrichment_ora_min_log2fc_merge", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05))),
                   numericInput("enrichment_show_max_p_adj_merge", "Show terms with adjusted p-value below:", value = 0.05, min = 0, max = 1, step = 0.01),
                   checkboxInput("enrichment_show_all_merge", "Tables: also show non-significant terms", value = FALSE)
            )
          ),
          fluidRow(
            column(12,
                   actionButton("enrichment_run_merge", tagList(icon("play"), " Run Enrichment"), class = "btn-white-violet btn-full-width", style = "margin-top: 10px; margin-bottom: 15px;"),
                   tags$p(icon("box-archive"), " Each run is stored with the other result tables: it is saved with the session and its parameters appear in the report. Delete unwanted runs in Save & Session.",
                          style = "font-size: 12px; color: #666; margin-top: -8px;")
            )
          ),
          tags$div(
            class = "plot-container",
            verbatimTextOutput("enrichment_summary_merge")
          )
        ),
        
        # Results
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("diagram-project"), " Enrichment Results", class = "box-title"),
          tags$div(
            class = "de-tabs",
            tabsetPanel(
              tabPanel("Over-representation (ORA)",
                       tags$div(class = "plot-container", style = "min-height: 550px; margin-top: 15px;", plotOutput("enrichment_ora_plot_merge", height = "550px")),
                       fluidRow(
                         column(4, numericInput("enrichment_ora_n_terms_merge", "Terms shown:", value = 15, min = 1, max = 50)),
                         column(4, downloadButton("download_enrichment_ora_plot_merge", "Download Plot (PNG)", class = "btn-white-green", style = "width: 100%; margin-top: 25px;")),
                         column(4, downloadButton("download_enrichment_ora_table_merge", "Download Table (Excel)", class = "btn-white-green", style = "width: 100%; margin-top: 25px;"))
                       ),
                       tags$div(class = "plot-container", style = "margin-top: 15px;", DT::DTOutput("enrichment_ora_table_merge"))),
              tabPanel("Gene set enrichment (GSEA)",
                       tags$div(class = "plot-container", style = "min-height: 550px; margin-top: 15px;", plotOutput("enrichment_gsea_plot_merge", height = "550px")),
                       fluidRow(
                         column(4, numericInput("enrichment_gsea_n_terms_merge", "Terms shown:", value = 15, min = 1, max = 50)),
                         column(4, downloadButton("download_enrichment_gsea_plot_merge", "Download Plot (PNG)", class = "btn-white-green", style = "width: 100%; margin-top: 25px;")),
                         column(4, downloadButton("download_enrichment_gsea_table_merge", "Download Table (Excel)", class = "btn-white-green", style = "width: 100%; margin-top: 25px;"))
                       ),
                       tags$div(class = "plot-container", style = "margin-top: 15px;", DT::DTOutput("enrichment_gsea_table_merge"))),
              tabPanel("Parameters",
                       tags$div(class = "plot-container", style = "margin-top: 15px;", tableOutput("enrichment_params_table_merge")),
                       downloadButton("download_enrichment_params_merge", "Download Parameters (CSV)", class = "btn-white-green", style = "margin-top: 10px;"))
            )
          )
        )
      ),
      ############################## Multiple/Exclusive Biomarkers ##############################
      tabItem(
        tabName = "exclusive_biomarkers_merge",
        
        tags$div(
          class = "header-violet",
          tags$h2("Exclusive Biomarkers Discovery - Multiple Datasets", class = "header-title"),
          tags$p("Identify genes predominantly expressed in specific clusters across integrated datasets", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8eaf6; padding: 20px; border-radius: 8px; border-left: 4px solid #667eea; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Exclusive Biomarkers in Integrated Data", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Find genes that are specifically expressed in target cluster(s) across multiple integrated datasets."),
            tags$h5(icon("flask"), " Analysis Options:", style = "color: #d97706; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(strong("By Cluster:"), " Find markers for specific cluster(s) across all datasets"),
              tags$li(strong("By Dataset:"), " Find markers for cluster(s) within a specific dataset"),
              tags$li(strong("Dataset + Cluster:"), " Combine both filters for maximum specificity")
            )
          )
        ),
        tags$div(
          class = "box-gradient-green",
          tags$h4(icon("sliders-h"), " Analysis Parameters", class = "box-title"),
          fluidRow(
            column(width = 3,
                   selectInput("exclusive_groupby_col_merge",
                               "Group by column:",
                               choices = NULL),
                   tags$small("Metadata column to split cells by (e.g. dataset, condition, genotype)", style = "color: #666;")
            ),
            column(width = 3,
                   selectInput("exclusive_target_groups_merge",
                               "Target group(s):",
                               choices = NULL,
                               multiple = TRUE),
                   tags$small("Group(s) to find exclusive markers for", style = "color: #666;")
            ),
            column(width = 3,
                   selectInput("exclusive_compare_groups_merge",
                               "Compare against group(s):",
                               choices = NULL,
                               multiple = TRUE),
                   tags$small("Leave empty to compare against all other groups", style = "color: #666;")
            ),
            column(width = 3,
                   selectInput("exclusive_target_cluster_merge",
                               "Restrict to cluster(s):",
                               choices = NULL,
                               multiple = TRUE),
                   tags$small("Optional — leave empty to use all cells", style = "color: #666;")
            )
          ),
          fluidRow(
            style = "margin-top: 10px;",
            column(width = 3,
                   numericInput("exclusive_min_pct_target_merge", 
                                "Min % in Target:", 
                                value = 50, min = 0, max = 100, step = 5),
                   tags$small("Minimum % of cells expressing gene in target. 50% (default) favours robust markers; use 30% for genes captured at low rates (transcription factors, single-nucleus data, Visium HD bins).", style = "color: #666;")            ),
            column(width = 3,
                   numericInput("exclusive_max_pct_other_merge", 
                                "Max % in Others:", 
                                value = 25, min = 0, max = 100, step = 5),
                   tags$small("Maximum % of cells expressing gene in others", style = "color: #666;")
            ),
            column(width = 3,
                   numericInput("exclusive_min_log2fc_merge", 
                                "Min Log2FC:", 
                                value = 1.5, min = 0, max = 5, step = 0.1),
                   tags$small("Minimum fold change (log2 scale)", style = "color: #666;")
            ),
            column(width = 3,
                   numericInput("exclusive_detection_threshold_merge", 
                                "Detection Threshold:", 
                                value = 0, min = 0, max = 5, step = 0.1),
                   tags$small("Minimum expression to consider gene as detected", style = "color: #666;")
            )
          ),
          fluidRow(
            style = "margin-top: 10px;",
            column(width = 3,
                   numericInput("exclusive_min_mean_expr_merge", 
                                "Min Mean Expression:", 
                                value = 0.5, min = 0, max = 10, step = 0.1),
                   tags$small("Minimum mean expression in target", style = "color: #666;")
            ),
            column(width = 3,
                   selectInput("exclusive_statistical_test_merge", 
                               "Statistical Test:", 
                               choices = c("Wilcoxon" = "wilcox", "None" = "none"),
                               selected = "wilcox"),
                   tags$small("Use Wilcoxon test for validation", style = "color: #666;")
            ),
            column(width = 3,
                   numericInput("exclusive_max_pvalue_merge", 
                                "Max P-value:", 
                                value = 0.05, min = 0.001, max = 1, step = 0.01),
                   tags$small("Maximum adjusted p-value", style = "color: #666;")
            )
          ),
          fluidRow(
            style = "margin-top: 15px;",
            column(width = 6,
                   actionButton("find_exclusive_markers_merge", 
                                tagList(icon("rocket"), " Find Exclusive Biomarkers"), 
                                class = "btn-white-green", 
                                style = "width: 100%;")                             
            ),
            column(width = 6,
                   downloadButton("download_exclusive_markers_merge", 
                                  "Download Results", 
                                  class = "btn-white-green", 
                                  style = "width: 100%;"),tags$br()
            )
          )
        ),
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("table"), " Exclusive Biomarkers Results", class = "box-title"),
          tags$div(
            class = "plot-container",
            DTOutput("table_exclusive_markers_merge")
          )
        ),
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("chart-bar"), " Expression Visualization", class = "box-title"),
          tags$div(
            style = "background-color: #d1fae5; padding: 15px; border-radius: 5px; margin-bottom: 15px;",
            tags$p(
              icon("info-circle"), 
              " Enter gene names from the results table above to visualize their expression pattern across clusters and datasets.",
              style = "margin: 0; color: #065f46;"
            )
          ),
          fluidRow(
            column(width = 3,
                   
                   selectInput("exclusive_plot_type_merge", 
                               "Plot Type:", 
                               choices = c("Dot Plot" = "dotplot", 
                                           "Violin Plot" = "violin",
                                           "Feature Plot" = "feature"),
                               selected = "dotplot"),
                   selectInput("exclusive_split_by_column_merge", 
                               "Split by column (optional):", 
                               choices = c("None" = ""),
                               selected = ""),
                   uiOutput("exclusive_split_values_filter_ui")
            ),
            column(width = 3,
                   textInput("exclusive_genes_to_plot_merge", 
                             "Gene Names (separated by , ; or spaces):", 
                             placeholder = "e.g., Gene1, Gene2, Gene3")
            ),
            column(width = 3,
                   tags$br(),
                   tags$br(),
                   actionButton("generate_exclusive_plot_merge", 
                                tagList(icon("chart-line"), " Generate Plot"), 
                                class = "btn-white-green",
                                style = "width: 100%;"),
                   tags$br(),
                   tags$br(),
                   
                   downloadButton("download_exclusive_plot_merge", 
                                  "Download Plot", 
                                  class = "btn-white-green", 
                                  style = "width: 100%;")
            ),
            column(width = 3,
                   numericInput("exclusive_plot_dpi_merge", "Resolution (DPI):", 
                                value = 300, min = 72, max = 600, step = 72),
                   selectInput("exclusive_plot_format_merge", "Format:", 
                               choices = c("TIFF" = "tiff", "PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg", "SVG" = "svg"), 
                               selected = "tiff")
                   
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px; margin-top: 20px;",
            plotOutput("plot_exclusive_markers_merge", height = "600px")
          )
        )
      ),
      
      
      ############################## Multiple/DEG Visualizations ##############################
      tabItem(
        tabName = "deg_visualizations_merge",
        
        tags$div(
          class = "header-violet",
          tags$h2("DEG Visualizations - Multiple Datasets", class = "header-title"),
          tags$p("Explore differential expression results through Venn diagrams and volcano plots", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8eaf6; padding: 20px; border-radius: 8px; border-left: 4px solid #667eea; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " DEG Visualization Tools", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Visualize and compare differential expression results from integrated dataset analyses."),
            tags$h5(icon("chart-bar"), " Available Visualizations:", style = "color: #667eea; margin-top: 15px;"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(strong("Venn Diagrams:"), "Compare gene lists across different comparisons to identify overlaps and unique markers."),
              tags$li(strong("Volcano Plots:"), "Visualize the relationship between fold change and statistical significance.")
            )
          )
        ),
        
        # Venn Diagram Section
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 20px;",
          tags$h4(icon("circle-notch"), " Venn Diagram Comparison", class = "box-title"),
          fluidRow(
            column(4,
                   selectInput("venn_table_1", "Select first gene list:", choices = c("None" = ""), selected = ""),
                   checkboxInput("significant_only_venn_1", "Include only significant genes", value = TRUE),
                   numericInput("log_fc_threshold_venn_1", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05)
            ),
            column(4,
                   selectInput("venn_table_2", "Select second gene list:", choices = c("None" = ""), selected = ""),
                   checkboxInput("significant_only_venn_2", "Include only significant genes", value = TRUE),
                   numericInput("log_fc_threshold_venn_2", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05)
            ),
            column(4,
                   selectInput("venn_table_3", "Select third gene list (optional):", choices = c("None" = ""), selected = ""),
                   checkboxInput("significant_only_venn_3", "Include only significant genes", value = TRUE),
                   numericInput("log_fc_threshold_venn_3", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05)
            )
          ),
          fluidRow(
            column(4,
                   selectInput("venn_direction", "Gene selection criteria:", 
                               choices = c("Up-regulated (avg_log2FC > threshold)" = "up", 
                                           "Down-regulated (avg_log2FC < -threshold)" = "down", 
                                           "Both directions (|avg_log2FC| > threshold)" = "both"), 
                               selected = "up")
            ),
            column(4,
                   numericInput("p_val_threshold_venn", "p-value threshold:", value = 0.05, min = 0, max = 1, step = 0.01)
            ),
            column(4,
                   checkboxInput("use_adjusted_p_venn", "Use adjusted p-values", value = TRUE),
                   actionButton("generate_venn_btn", tagList(icon("circle-notch"), " Generate Venn Diagram"), 
                                class = "btn-white-violet btn-full-width", style = "margin-top: 23px;")
            )
          ),
          fluidRow(
            column(4, colourInput("venn_color_1", "Color for first set:", value = "#56B4E9")),
            column(4, colourInput("venn_color_2", "Color for second set:", value = "#E69F00")),
            column(4, colourInput("venn_color_3", "Color for third set (if used):", value = "#009E73"))
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 400px; margin-top: 15px;",
            plotOutput("venn_plot", height = "400px")
          ),
          fluidRow(
            column(6,
                   selectInput("venn_diagram_format", "Export format:", 
                               choices = c("PNG" = "png", "TIFF" = "tiff", "PDF" = "pdf", "JPEG" = "jpeg"), 
                               selected = "png"),
                   numericInput("venn_diagram_dpi", "DPI:", value = 300, min = 72, step = 72)
            ),
            column(6,
                   tags$br(),
                   tags$br(),
                   
                   downloadButton("download_venn_diagram", "Download Venn Diagram", 
                                  class = "btn-white-violet", style = "width: 100%; margin-bottom: 10px;"),
                   downloadButton("download_venn_gene_lists", "Download Gene Lists", 
                                  class = "btn-white-violet", style = "width: 100%;")
            )
          ),
          tags$div(
            class = "plot-container",  
            style = "margin-top: 15px;",
            selectInput("selected_gene_set", "View gene list:", choices = NULL),
            DTOutput("venn_gene_table")
          )
        ),
        
        # Volcano Plot Section
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px;",
          tags$h4(icon("volcano"), " Volcano Plot", class = "box-title"),
          fluidRow(
            column(4,
                   selectInput("volcano_deg_source_merge", "Select DEG results:", 
                               choices = c("None" = ""), selected = ""),
                   numericInput("volcano_log2fc_threshold_merge", "Log2FC threshold:", 
                                value = 0.5, min = 0, max = 5, step = 0.1)
            ),
            column(4,
                   numericInput("volcano_pval_threshold_merge", "Adjusted p-value threshold:", 
                                value = 0.05, min = 0, max = 1, step = 0.01),
                   numericInput("volcano_label_genes_merge", "Number of genes to label:", 
                                value = 10, min = 0, max = 50, step = 1)
            ),
            column(4,
                   numericInput("volcano_point_size_merge", "Point size:", 
                                value = 2, min = 0.5, max = 5, step = 0.5),
                   numericInput("volcano_point_alpha_merge", "Point transparency:", 
                                value = 0.6, min = 0.1, max = 1, step = 0.1)
            )
          ),
          fluidRow(
            column(4,
                   colourInput("volcano_color_up_merge", "Color for upregulated genes:", value = "#00ba38")
            ),
            column(4,
                   colourInput("volcano_color_down_merge", "Color for downregulated genes:", value = "#f8766d")
            ),
            column(4,
                   colourInput("volcano_color_ns_merge", "Color for non-significant genes:", value = "#619cff")
            )
          ),
          fluidRow(
            column(12,
                   actionButton("generate_volcano_btn_merge", 
                                tagList(icon("volcano"), " Generate Volcano Plot"), 
                                class = "btn-white-green btn-full-width", 
                                style = "margin-top: 10px; margin-bottom: 15px;")
            )
          ),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px; margin-top: 15px;",
            plotOutput("volcano_plot_merge", height = "500px")
          ),
          fluidRow(
            column(6,
                   selectInput("volcano_plot_format_merge", "Export format:", 
                               choices = c("PNG" = "png", "TIFF" = "tiff", "PDF" = "pdf", "JPEG" = "jpeg", "SVG" = "svg"), 
                               selected = "png"),
                   numericInput("volcano_plot_dpi_merge", "DPI:", value = 300, min = 72, step = 72)
            ),
            column(6,
                   downloadButton("download_volcano_plot_merge", "Download Volcano Plot", 
                                  class = "btn-white-green", style = "width: 100%; margin-top: 25px;")
            )
          )
        )
      ),
      
      
      ############################# Multiple dataset analysis / Subset #############################
      
      tabItem(
        tabName = "subset_merge",
        
        tags$div(
          class = "header-violet",
          tags$h2("Data Subsetting - Multiple Datasets", class = "header-title"),
          tags$p("Create focused subsets from integrated datasets", class = "header-subtitle")
        ),
        
        tags$div(
          style = "background-color: #e8eaf6; padding: 20px; border-radius: 8px; border-left: 4px solid #667eea; margin-bottom: 25px;",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Multi-Dataset Subsetting", tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Select specific cell populations from your integrated datasets using multiple criteria."),
            tags$h5(icon("star"), " Subsetting Methods:", style = "color: #667eea; margin-top: 15px;"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("By Clusters:"), "Select specific clusters to retain across all datasets"),
              tags$li(tags$strong("By Gene Expression:"), "Filter cells based on expression of key marker genes"),
              tags$li(tags$strong("By Metadata:"), "Use dataset-specific metadata for targeted subsetting")
            )
          )
        ),
        
        tags$div(
          class = "box-gradient-green",
          tags$h4(icon("project-diagram"), " Original Integrated Dataset", class = "box-title"),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px;",
            plotOutput("global_umap_merge", height = "500px")
          )
        ),
        
        tags$h4(icon("filter"), " Choose Subsetting Method", style = "color: #495057; margin: 30px 0 20px 0; font-weight: bold;"),
        
        fluidRow(
          column(width = 4,
                 tags$div(
                   class = "box-gradient-violet box-thin",
                   tags$h5(icon("sitemap"), " Subset by Clusters", class = "box-title", style = "font-size: 18px;"),
                   tags$div(
                     class = "plot-container-thin",
                     selectInput("select_ident_subset_merge", "Select clusters to include:", choices = NULL, multiple = TRUE),
                     actionButton("apply_subset_merge", tagList(icon("cut"), " Apply Cluster Subset"), class = "btn-white-violet btn-full-width")
                   )
                 )
          ),
          
          column(width = 4,
                 tags$div(
                   class = "box-gradient-green box-thin",
                   tags$h5(icon("dna"), " Subset by Gene Expression", class = "box-title", style = "font-size: 18px;"),
                   tags$div(
                     class = "plot-container-thin",
                     numericInput("expression_threshold_merge", "Expression threshold:", value = 1, min = 0),
                     textInput("gene_list_merge", "Genes (separated by , ; or spaces):", value = ""),
                     numericInput("num_genes_to_express_merge", "Minimum expressed genes:", value = 1, min = 1),
                     checkboxInput("negative_gene_subset_merge", "Negative selection (exclude expressing cells)", value = FALSE),
                     actionButton("apply_gene_subset_merge", tagList(icon("cut"), " Apply Gene Subset"), class = "btn-white-green btn-full-width")
                   )
                 )
          ),
          
          column(width = 4,
                 tags$div(
                   class = "box-gradient-violet box-thin",
                   tags$h5(icon("tags"), " Subset by Metadata", class = "box-title", style = "font-size: 18px;"),
                   tags$div(
                     class = "plot-container-thin",
                     selectInput("metadata_column_subset", "Select metadata column:", choices = NULL),
                     uiOutput("metadata_values_ui"),
                     actionButton("apply_metadata_subset", tagList(icon("cut"), " Apply Metadata Subset"), class = "btn-white-violet btn-full-width")
                   )
                 )
          )
        ),
        
        tags$div(
          class = "box-gradient-green",
          style = "margin-top: 30px;",
          tags$h4(icon("eye"), " Subset Preview", class = "box-title"),
          tags$div(
            class = "plot-container",
            style = "min-height: 500px;",
            plotOutput("subset_umap_merge", height = "500px")
          )
        ),
        
        tags$div(
          class = "box-gradient-violet",
          style = "margin-top: 20px; text-align: center; padding: 30px;",
          tags$h4(icon("save"), " Save Subset", class = "box-title", style = "margin-bottom: 20px;"),
          downloadButton("download_subset_merge", "Save Subset as .RDS", class = "btn-white-violet", style = "padding: 12px 40px; font-size: 18px;")
        )
      ),
      ############################## Multiple/Session ##############################
      tabItem(
        tabName = "session_merge",
        tags$div(
          class = "header-violet",
          tags$h2("Save & Session", class = "header-title"),
          tags$p("Save the whole analysis in one file, export the analysis parameters, see what the object contains and manage stored result tables", class = "header-subtitle")
        ),
        fluidRow(
          column(5,
                 tags$div(
                   class = "box-gradient-violet",
                   tags$h4(icon("floppy-disk"), " Save analysis", class = "box-title"),
                   tags$div(
                     style = "background: white; border-radius: 6px; padding: 15px; color: #2c3e50;",
                     tags$p("One .rds file containing the Seurat object (clusters, cell type names, colors, metadata, reductions, analysis parameters) and all stored result tables. Plots are not stored: they are regenerated from the object."),
                     tags$p(icon("rotate-left"), " To resume: ", tags$b("Load datasets → pre-integrated Seurat object → select this file."), style = "font-size: 12px; color: #6c757d;"),
                     downloadButton("download_session_merge", tagList(icon("floppy-disk"), " Save analysis (.rds)"),
                                    class = "btn-gradient-violet", style = "width: 100%;")
                   )
                 ),
                 tags$div(
                   class = "box-gradient-green",
                   style = "margin-top: 20px;",
                   tags$h4(icon("file-export"), " Export report", class = "box-title"),
                   tags$div(
                     style = "background: white; border-radius: 6px; padding: 15px; color: #2c3e50;",
                     tags$p("Every analysis parameter actually applied (QC thresholds, normalization, integration, PCs, resolution, seeds, annotation), as a CSV table or a PDF report.",
                            style = "font-size: 12px; color: #6c757d;"),
                     fluidRow(
                       column(6, downloadButton("download_params_multi", tagList(icon("file-csv"), " Parameters (CSV)"), class = "btn-gradient-green btn-full-width")),
                       column(6, downloadButton("download_report_multi", tagList(icon("file-pdf"), " PDF report"), class = "btn-gradient-green btn-full-width"))
                     )
                   )
                 ),
                 tags$div(
                   class = "box-gradient-green",
                   style = "margin-top: 20px;",
                   tags$h4(icon("circle-info"), " Object contents", class = "box-title"),
                   tags$div(
                     style = "background: white; border-radius: 6px; padding: 10px; color: #2c3e50;",
                     tableOutput("session_object_summary_merge")
                   )
                 )
          ),
          column(7,
                 tags$div(
                   class = "box-gradient-green",
                   tags$h4(icon("table"), " Stored result tables", class = "box-title"),
                   tags$div(
                     style = "background: white; border-radius: 6px; padding: 15px; color: #2c3e50;",
                     tags$p("Tables created by cluster comparisons and exclusive biomarker searches (used by Venn diagrams and volcano plots). Select rows, then delete them if they are no longer needed.",
                            style = "font-size: 12px; color: #6c757d;"),
                     DT::DTOutput("session_tables_merge"),
                     tags$div(
                       style = "margin-top: 10px; display: flex; justify-content: space-between; align-items: center; gap: 10px; flex-wrap: wrap;",
                       tags$span(textOutput("session_size_info_merge", inline = TRUE), style = "font-size: 12px; color: #6c757d;"),
                       actionButton("delete_session_tables_merge", tagList(icon("trash-alt"), " Delete selected"),
                                    class = "btn-gradient-violet")
                     )
                   )
                 )
          )
        )
      ),
      ############################## Cell Chat load data ##############################
      tabItem(
        tabName = "load_data_cellchat",
        
        tags$div(
          class = "header-cellchat",
          tags$h2("CellChat Analysis"),
          tags$p("Uncover cell-cell communication through ligand-receptor interactions")
        ),
        
        tags$div(
          class = "info-box-cellchat",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About CellChat Analysis", 
                    tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("CellChat is a powerful tool for analyzing cell-cell communication through ligand-receptor interactions in single-cell transcriptomics data."),
            tags$h5(icon("star"), " Key Features:"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Ligand-Receptor Pairs:"), " Identify significant communication pathways"),
              tags$li(tags$strong("Cell-Cell Networks:"), " Visualize interaction patterns between cell types"),
              tags$li(tags$strong("Signaling Pathways:"), " Discover enriched communication pathways"),
              tags$li(tags$strong("Comparison:"), " Compare communication across conditions or datasets")
            ),
            tags$h5(icon("check-square"), " Requirements:"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li("Normalized Seurat object with cell type annotations"),
              tags$li("CellHub Database: 76,582 curated ligand-receptor interactions across 21 species from 10 source databases"),
              tags$li("Minimum of 2 cell types/populations")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Ensure your cell types are clearly annotated before starting CellChat analysis", style = "margin: 0;")
            )
          )
        ),
        
        fluidRow(
          column(6,
                 tags$div(
                   class = "card-cellchat-medium",
                   tags$h4(icon("database"), " Load Seurat Object", class = "box-title"),
                   tags$p("Upload your processed Seurat object with cell type annotations", 
                          style = "margin-bottom: 20px; opacity: 0.95;"),
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.15); padding: 12px; border-radius: 6px; margin-bottom: 15px;",
                     tags$label(icon("link"), " Use object already loaded in this session", style = "font-weight: bold;"),
                     radioButtons("cellchat_import_source", NULL,
                                  choices = c("Single dataset" = "single", "Multiple datasets" = "multiple"),
                                  selected = "single", inline = TRUE),
                     actionButton("import_session_object_cellchat", tagList(icon("arrow-right"), " Use current object"),
                                  class = "btn-white-green", style = "width: 100%;"),
                     tags$p(icon("info-circle"), " Snapshot taken at click time: re-import after renaming clusters upstream.",
                            style = "font-size: 11px; margin: 6px 0 0 0; opacity: 0.85;")
                   ),
                   fileInput("seurat_file_cellchat", "Or load from disk (.rds)",
                             accept = ".rds",
                             buttonLabel = tagList(icon("upload"), " Browse RDS"),
                             placeholder = "No object selected"),
                   tags$p(icon("info-circle"), " A Seurat object, or an analysis saved with 'Save Analysis' (CellChat objects are restored and ready to plot).",
                          style = "font-size: 11px; margin: -10px 0 10px 0; opacity: 0.85;"),
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.2); padding: 15px; 
                             border-radius: 6px; margin-top: 15px;",
                     tags$h5("Object Information:", style = "margin-top: 0; color: white;"),
                     verbatimTextOutput("seurat_info_cellchat")
                   )
                   
                   
                   
                   
                 )
          ),
          column(6,
                 tags$div(
                   class = "card-cellchat-rose",
                   tags$h4(icon("book"), " CellHub Database", class = "box-title"),
                   tags$p("A curated multi-species ligand-receptor interaction resource",
                          style = "font-style: italic; margin-bottom: 5px; opacity: 0.95;"),
                   tags$p(icon("circle"), " One database to unite them all",
                          style = "font-size: 11px; font-style: italic; margin-bottom: 10px; opacity: 0.7;"),
                   tags$p(
                     "76,582 interactions across 21 species, consolidated from 10 source databases ",
                     "(CellChatDB, CellPhoneDB, CellTalkDB, connectomeDB2025, NicheNet, ",
                     "NeuronChatDB, FlyPhoneDB2, cell2cell, PlantPhoneDB, PlantCellChatDB). ",
                     "Each interaction is curated, deduplicated, and traceable to its original source.",
                     style = "font-size: 12px; margin-bottom: 15px; opacity: 0.85;"
                   ),
                   # Species selector — all 21 CellHub Database species
                   selectInput("species_cellchat", "Select Species:",
                               choices = c("Mouse" = "mouse",
                                           "Human" = "human",
                                           "Rat" = "rat",
                                           "Zebrafish" = "zebrafish",
                                           "Drosophila" = "drosophila",
                                           "C. elegans" = "celegans",
                                           "Arabidopsis" = "arabidopsis",
                                           "Maize" = "maize",
                                           "Rice" = "rice",
                                           "Soybean" = "soybean",
                                           "Tomato" = "tomato",
                                           "Chimpanzee" = "chimp",
                                           "Macaque" = "macaque",
                                           "Marmoset" = "marmoset",
                                           "Cow" = "cow",
                                           "Dog" = "dog",
                                           "Horse" = "horse",
                                           "Pig" = "pig",
                                           "Sheep" = "sheep",
                                           "Chicken" = "chicken",
                                           "Frog" = "frog"),
                               selected = "mouse"),
                   # Database selection mode
                   radioButtons("db_mode_cellchat", "Database Selection:",
                                choices = c("All databases (merged)" = "merged",
                                            "Select individual databases" = "custom"),
                                selected = "merged",
                                inline = TRUE),
                   # Dynamic checkboxes — shown only when mode = "custom"
                   # Populated by server based on species
                   conditionalPanel(
                     condition = "input.db_mode_cellchat == 'custom'",
                     tags$div(
                       style = "background-color: rgba(255,255,255,0.15); padding: 10px; 
                                      border-radius: 6px; margin-bottom: 10px;",
                       checkboxGroupInput("db_sources_cellchat", "Available databases:",
                                          choices = NULL)
                     )
                   ),
                   # Interaction filters (confidence tier, predicted, complexes, gene list)
                   cellhubdb_filter_panel(),
                   # Load button
                   actionButton("load_db_cellchat",
                                tagList(icon("download"), " Load Database"),
                                class = "btn-white-cellchat-rose btn-full-width",
                                style = "margin-top: 10px;"),
                   # Database info display
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.2); padding: 15px; 
                             border-radius: 6px; margin-top: 15px;",
                     tags$h5("Database Information:", style = "margin-top: 0; color: white;"),
                     verbatimTextOutput("db_info_cellchat")
                   ),
                   # Genes of interest: named gene lists used by the database query and by the bubble and chord plots
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.2); padding: 15px; border-radius: 6px; margin-top: 15px;",
                     tags$h5(icon("list"), " Genes of interest:", style = "margin-top: 0; color: white;"),
                     tags$p("Paste gene symbols (for example differentially expressed genes), separated by spaces, commas or new lines, ",
                            "and give the list a name. Saved lists can be applied to the bubble and chord plots (Gene list selector) ",
                            "and searched in the loaded database. Matching ignores case and includes complex subunits; ",
                            "a list added with an existing name replaces it.",
                            style = "color: rgba(255,255,255,0.8); font-size: 11px; line-height: 1.35;"),
                     textInput("gene_list_name_cellchat", "List name:", placeholder = "DEG MuSC up", width = "100%"),
                     textAreaInput("gene_list_text_cellchat", "Genes:", rows = 4, width = "100%", placeholder = "Egfr Notch1 Dll1 Vcam1"),
                     actionButton("add_gene_list_cellchat", tagList(icon("plus"), " Add list"), class = "btn-white-cellchat-rose btn-full-width"),
                     tags$div(style = "background-color: white; color: #333; border-radius: 6px; padding: 6px; margin-top: 10px;",
                              tableOutput("gene_lists_table_cellchat")),
                     selectInput("gene_list_selected_cellchat", "Saved list:", choices = NULL, width = "100%"),
                     fluidRow(
                       column(8, actionButton("run_gene_query_cellchat", tagList(icon("search"), " Show interactions in the database"),
                                              class = "btn-white-cellchat-rose btn-full-width")),
                       column(4, actionButton("delete_gene_list_cellchat", tagList(icon("trash"), " Delete list"),
                                              class = "btn-white-cellchat-rose btn-full-width"))
                     ),
                     tags$div(style = "color: white; margin-top: 8px; font-size: 12px;", textOutput("gene_query_summary_cellchat")),
                     tags$div(style = "background-color: white; color: #333; border-radius: 6px; padding: 8px; margin-top: 8px;",
                              DT::dataTableOutput("gene_query_table_cellchat"))
                   )
                 )
          )
        ),
        
        tags$div(
          class = "logs-container-cellchat",
          tags$h4(icon("terminal"), " Processing Logs", 
                  style = "color: #6BB8B7; margin-top: 0; margin-bottom: 15px;"),
          tags$pre(
            style = "background-color: #1e272e; color: #6BB8B7; padding: 15px; 
                         border-radius: 6px; font-family: 'Courier New', monospace; 
                         font-size: 12px; max-height: 300px; overflow-y: auto; 
                         margin: 0; border-left: 3px solid #4AA3A2;",
            verbatimTextOutput("loading_logs")
          )
        ),
        
        fluidRow(
          style = "margin-top: 30px;",
          column(4,
                 tags$div(
                   class = "step-box-cellchat-turquoise",
                   tags$h6(icon("1"), strong(" Step 1: Load Data")),
                   tags$p("Upload your processed Seurat object with cell type annotations")
                 )
          ),
          column(4,
                 tags$div(
                   class = "step-box-cellchat-rose",
                   tags$h6(icon("2"), strong(" Step 2: Load Database")),
                   tags$p("Select species and load the CellHub Database for ligand-receptor interactions")
                 )
          ),
          column(4,
                 tags$div(
                   class = "step-box-cellchat-mixed",
                   tags$h6(icon("3"), strong(" Step 3: Analyze")),
                   tags$p("Proceed to 'Ligand-Receptor' tab to create and analyze CellChat objects")
                 )
          )
        ),
        # Demonstration data instructions (PBMC v3 processed in the Single module, then "Use current object")
        example_data_box("cellchat")
      ),
      
      ############################## Ligand receptor analysis ##############################
      tabItem(
        tabName = "ligand_receptor_cellchat",
        
        tags$div(
          class = "header-cellchat",
          tags$h2("Ligand-Receptor Analysis"),
          tags$p("Identify and visualize cell-cell communication patterns")
        ),
        
        tags$div(
          class = "info-box-cellchat",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Ligand-Receptor Analysis", 
                    tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Analyze cell-cell communication by identifying significant ligand-receptor interactions between cell populations."),
            tags$h5(icon("star"), " Analysis Workflow:"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Create Objects:"), " Generate CellChat objects from your Seurat data"),
              tags$li(tags$strong("Run Analysis:"), " Compute interaction probabilities and significance"),
              tags$li(tags$strong("Visualize:"), " Generate bubble plots showing communication patterns")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Use filters to focus on specific clusters and conditions", style = "margin: 0;")
            )
          )
        ),
        
        tags$div(
          class = "card-cellchat-light",
          tags$h4(icon("cogs"), " CellChat Object Creation & Analysis", class = "box-title"),
          
          # Explanation box
          tags$div(
            style = "background-color: rgba(50, 205, 50, 0.2); padding: 12px; border-radius: 6px; margin-bottom: 15px; border-left: 4px solid #32CD32;",
            tags$h6(icon("info-circle"), " How CellChat Groups Work", 
                    style = "color: white; margin-top: 0; font-weight: bold; font-size: 14px;"),
            tags$ul(
              style = "margin: 5px 0 0 15px; padding: 0; font-size: 11px; color: white; line-height: 1.5;",
              tags$li(tags$strong("Simple:"), " 1 column → Groups are values in that column (e.g., ClusterIdents → Myh7, MuscleStemCell...)"),
              tags$li(tags$strong("Combined:"), " 2-3 columns → Groups are unique combinations (e.g., ClusterIdents + condition → Myh7_WT, Myh7_Dmd...)")
            )
          ),
          
          tags$h5(icon("layer-group"), " Step 1: Define Cell Groups", 
                  style = "color: white; font-weight: bold; margin-top: 0; margin-bottom: 10px;"),
          
          fluidRow(
            column(4,
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.15); padding: 10px; border-radius: 6px;",
                     tags$h6("Column 1 (required)", style = "color: white; margin-top: 0; font-size: 13px;"),
                     selectInput("group_by_cellchat", 
                                 label = NULL,
                                 choices = NULL)
                   )
            ),
            column(4,
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.15); padding: 10px; border-radius: 6px;",
                     tags$h6("Column 2 (optional)", style = "color: white; margin-top: 0; font-size: 13px;"),
                     selectizeInput("combine_with_column", label = NULL, choices = c("None" = ""), selected = "",
                                    options = list(allowEmptyOption = TRUE))
                   )
            ),
            column(4,
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.15); padding: 10px; border-radius: 6px;",
                     tags$h6("Column 3 (optional)", style = "color: white; margin-top: 0; font-size: 13px;"),
                     selectizeInput("combine_with_column_2", label = NULL, choices = c("None" = ""), selected = "",
                                    options = list(allowEmptyOption = TRUE))
                   )
            )
          ),
          
          # Preview of created groups
          tags$div(
            style = "background-color: rgba(255,255,255,0.15); padding: 12px; border-radius: 6px; margin-top: 10px;",
            tags$h6(icon("eye"), " Preview of Groups", style = "color: white; margin-top: 0; font-size: 13px;"),
            uiOutput("cellchat_groups_preview")
          ),
          
          # Advanced Filtering section
          tags$div(
            style = "background-color: rgba(255,255,255,0.15); padding: 15px; border-radius: 6px; margin-top: 15px;",
            tags$div(
              style = "cursor: pointer;",
              onclick = "$(this).next().slideToggle();",
              tags$h6(icon("filter"), " Step 2: Advanced Filtering (Optional) ", 
                      tags$span(icon("chevron-down"), style = "float: right; font-size: 12px;"),
                      style = "color: white; font-weight: bold; margin: 0;")
            ),
            tags$div(
              style = "display: none; margin-top: 12px;",
              fluidRow(
                column(6,
                       tags$div(
                         style = "background-color: rgba(255,255,255,0.1); padding: 10px; border-radius: 4px; height: 100%;",
                         tags$h6("Keep Only Specific Groups", style = "color: white; margin-top: 0; font-size: 13px;"),
                         tags$p("Select which groups to include in analysis", 
                                style = "font-size: 11px; color: rgba(255,255,255,0.7); margin-bottom: 8px;"),
                         uiOutput("cellchat_groups_filter_ui")
                       )
                ),
                column(6,
                       tags$div(
                         style = "background-color: rgba(255,255,255,0.1); padding: 10px; border-radius: 4px; height: 100%;",
                         tags$h6("Filter by Additional Metadata", style = "color: white; margin-top: 0; font-size: 13px;"),
                         tags$p("Apply filters on other columns", 
                                style = "font-size: 11px; color: rgba(255,255,255,0.7); margin-bottom: 8px;"),
                         selectInput("condition_filter_column", 
                                     "Filter column:", 
                                     choices = NULL),
                         uiOutput("condition_values_ui")
                       )
                )
              )
            )
          ),
          
          tags$div(
            style = "margin-top: 15px;",
            fluidRow(
              column(6, 
                     textInput("cellchat_object_name", 
                               "Object name:", 
                               value = "my_object",
                               placeholder = "e.g., WT_vs_Dmd_7d")
              ),
              column(3, 
                     tags$br(),
                     actionButton("create_and_analyze_cellchat", 
                                  tagList(icon("play-circle"), " Create & Analyze"), 
                                  class = "btn-white-cellchat btn-full-width",
                                  style = "font-size: 16px; padding: 12px;")
              ),
              column(3,
                     tags$br(),
                     downloadButton("download_analyzed_objects", 
                                    tagList(icon("save"), " Save Analysis"), 
                                    class = "btn-white-cellchat btn-full-width",
                                    style = "font-size: 16px; padding: 12px;")
              )
            )
            
          ),
          tags$div(
            style = "background-color: rgba(255,255,255,0.1); padding: 15px; border-radius: 6px; margin-top: 15px; max-height: 400px; overflow-y: auto;",
            tags$h6(icon("list"), " Analysis Logs:", style = "color: white; margin-top: 0; font-size: 14px;"),
            verbatimTextOutput("analysis_logs_cellchat")
          )
        ),
        tags$div(
          class = "box-gradient-cellchat-rose",
          style = "margin-top: 20px; color: white;",
          tags$h4(icon("circle"), " Bubble Plot Visualization", class = "box-title"),
          
          fluidRow(
            column(3, selectInput("subset_choose_cellchat", "Select CellChat object:", choices = NULL)),
            column(3, selectizeInput("sources_use_cellchat", "Source cell types:", choices = NULL, multiple = TRUE, 
                                     options = list(maxItems = NULL, plugins = list('remove_button'), placeholder = 'Select sources'))),
            column(2, numericInput("plot_width", "Width (px):", value = 900, min = 400, max = 3000, step = 100)),
            column(4,   tags$br(),  downloadButton("download_bubble_plot", "Download Plot", class = "btn-white-cellchat-rose", style = "width: 100%;"))
          ),
          
          fluidRow(
            column(3, numericInput("threshold_cellchat", "Threshold:", value = 0.05, min = 0, max = 1, step = 0.01)),
            column(3, selectizeInput("targets_use_cellchat", "Target cell types:", choices = NULL, multiple = TRUE, 
                                     options = list(maxItems = NULL, plugins = list('remove_button'), placeholder = 'Select targets'))),
            column(2, numericInput("plot_height", "Height (px):", value = 1500, min = 400, max = 3000, step = 100)),
            
            column(4, tags$br(),  downloadButton("download_info_table", "Download Table", class = "btn-white-cellchat-rose", style = "width: 100%;"))
          ),
          fluidRow(
            column(3, tags$br(), actionButton("generate_plot_cellchat", tagList(icon("chart-line"), " Generate Plot"), 
                                              class = "btn-white-cellchat-rose btn-full-width", style = "font-size: 16px; padding: 12px;")),
            column(3,                            tags$br(),
                   checkboxInput("flip_axes", "Flip axes", value = FALSE)),
            column(2,                            tags$br(),
                   checkboxInput("exclude_intra", "Exclude intra-group", value = TRUE)),
            column(4, selectInput("bubble_plot_format", "Format:", 
                                  choices = c("PNG" = "png", "JPEG" = "jpeg", "TIFF" = "tiff", "SVG" = "svg", "PDF" = "pdf"), 
                                  selected = "png"))
          ),
          # Gene list saved in the Load Data tab, applied to the significant interactions after inference
          fluidRow(
            column(6, selectInput("bubble_gene_list_cellchat", "Gene list (from Load Data, Genes of interest):",
                                  choices = c("None (all interactions)" = ""), width = "100%", selectize = FALSE))
          )
        ),
        
        tags$div(
          class = "plot-container",
          style = "margin-top: 20px; background-color: white; padding: 20px; border-radius: 8px;",
          fluidRow(
            plotOutput("bubble_plot_cellchat", height = "800px")
            
          )
        )
      ),
      
      ############################## Circle Plot ##############################
      tabItem(
        tabName = "circle_plot_cellchat",
        
        tags$div(
          class = "header-cellchat",
          tags$h2("Circle Plot Analysis"),
          tags$p("Visualize communication networks with chord diagrams")
        ),
        
        tags$div(
          class = "info-box-cellchat",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Chord Plots", 
                    tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Chord diagrams visualize communication networks between cell types. The width of each chord represents the strength of communication through specific ligand-receptor pairs."),
            tags$h5(icon("star"), " Key Features:"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Sender-Receiver:"), " Visualize directional communication between cell types"),
              tags$li(tags$strong("Chord Thickness:"), " Width indicates interaction strength"),
              tags$li(tags$strong("Custom Colors:"), " Personalize colors for each cell type"),
              tags$li(tags$strong("Flexible Export:"), " Save in multiple formats (PNG, PDF, SVG, TIFF)")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Adjust plot dimensions for better visibility when many interactions are present", style = "margin: 0;")
            )
          )
        ),
        
        tags$div(
          class = "box-gradient-cellchat-dark",
          style = "color: white;",
          tags$h4(icon("exchange-alt"), " Chord Plot Parameters", class = "box-title"),
          fluidRow(
            column(3,
                   tags$h5("Data Selection", style = "color: white; font-weight: bold;"),
                   selectInput("cellchat_obj_chord", "CellChat object:", choices = NULL),
                   selectizeInput("sender_groups", "Sender cell types:", 
                                  choices = NULL, 
                                  multiple = TRUE,
                                  options = list(placeholder = "Select senders")),
                   selectizeInput("receiver_groups", "Receiver cell types:", 
                                  choices = NULL, 
                                  multiple = TRUE,
                                  options = list(placeholder = "Select receivers"))
            ),
            column(3,
                   tags$h5("Plot Parameters", style = "color: white; font-weight: bold;"),
                   numericInput("prob_threshold_chord", "Max p-value:", 
                                value = 0.05, min = 0, max = 1, step = 0.01),
                   numericInput("chord_top_n", "Top interactions (by probability, 0 = all):",
                                value = 30, min = 0, step = 5),
                   sliderInput("chord_label_size", "Label size:", 
                               value = 0.8, min = 0.4, max = 4, step = 0.1),
                   # Gene list saved in the Load Data tab, applied before the top-n selection
                   selectInput("chord_gene_list_cellchat", "Gene list (from Load Data):",
                               choices = c("None (all interactions)" = ""), width = "100%", selectize = FALSE)
            ),
            column(3,
                   tags$h5("Display Settings", style = "color: white; font-weight: bold;"),
                   sliderInput("chord_plot_dims", "Plot size (px):", 
                               value = 800, min = 400, max = 2000, step = 100),
                   checkboxInput("use_custom_colors", "Custom colors", value = FALSE),
                   conditionalPanel(
                     condition = "input.use_custom_colors == true",
                     uiOutput("color_inputs_chord")
                   )
            ),
            column(3,
                   tags$h5("Export Options", style = "color: white; font-weight: bold;"),
                   selectInput("chord_plot_format", "Format:",
                               choices = c("PNG" = "png", "PDF" = "pdf", 
                                           "SVG" = "svg", "TIFF" = "tiff"),
                               selected = "png"),
                   numericInput("chord_export_dpi", "Resolution (DPI):", 
                                value = 300, min = 72, max = 1200, step = 50),
                   tags$br(),
                   actionButton("generate_chord_plot", 
                                tagList(icon("chart-line"), " Generate Plot"), 
                                class = "btn-white-cellchat btn-full-width"),
                   tags$br(), tags$br(),
                   downloadButton("download_chord_plot", "Download Plot",
                                  class = "btn-white-cellchat", style = "width: 100%;")
            )
          )     ,   
          
          tags$div(
            class = "expandable-plot-container",
            id    = "chord_plot_expand",
            uiOutput("chord_plot_ui")
          )
        ),
        tags$hr(),
        tags$hr(),
        tags$div(
          class = "box-gradient-cellchat-dark",
          style = "color: white; margin-top: 20px;",
          tags$h4(icon("chart-pie"), " Split Chord Plot — Compare by Group", class = "box-title"),
          fluidRow(
            
            column(3,
                   tags$h5("Object & Groups", style = "color: white; font-weight: bold;"),
                   selectInput("cellchat_obj_split_chord", "CellChat object:", choices = NULL),
                   radioButtons("split_chord_mode", "Comparison mode:",
                                choices  = c("Split by group" = "split", "Compare pairs" = "pairs"),
                                selected = "split", inline = TRUE),
                   
                   conditionalPanel(
                     condition = "input.split_chord_mode == 'split'",
                     selectizeInput("split_chord_senders", "Sender cell types:",
                                    choices = NULL, multiple = TRUE,
                                    options = list(placeholder = "Select senders")),
                     selectizeInput("split_chord_receivers", "Receiver cell types:",
                                    choices = NULL, multiple = TRUE,
                                    options = list(placeholder = "Select receivers")),
                     selectInput("split_chord_by", "Split by:",
                                 choices  = c("Source" = "source", "Target" = "target"),
                                 selected = "target")
                   ),
                   
                   conditionalPanel(
                     condition = "input.split_chord_mode == 'pairs'",
                     numericInput("split_chord_n_groups", "Number of groups:",
                                  value = 2, min = 2, max = 8, step = 1),
                     uiOutput("split_chord_pairs_ui")
                   ),
                   
                   tags$hr(style = "border-color: rgba(255,255,255,0.3); margin: 12px 0;"),
                   tags$h5("Appearance", style = "color: white; font-weight: bold;"),
                   sliderInput("split_chord_label_distance", "Label distance:",
                               value = 0.5, min = 0.0, max = 2.0, step = 0.1),
                   sliderInput("split_chord_dims", "Plot size (px):",
                               value = 900, min = 400, max = 1600, step = 100),
                   checkboxInput("split_chord_custom_colors", "Custom colors", value = FALSE),
                   conditionalPanel(
                     condition = "input.split_chord_custom_colors == true",
                     uiOutput("split_chord_color_inputs")
                   ),
                   checkboxInput("split_chord_show_title", "Show title", value = FALSE)
            ),
            
            # --- Column 2: LR pair selection ---
            column(3,
                   tags$h5("LR Pairs", style = "color: white; font-weight: bold;"),
                   fluidRow(
                     column(6,
                            numericInput("split_chord_threshold", "Min probability:",
                                         value = 0.05, min = 0, max = 1, step = 0.01)
                     ),
                     column(6,
                            tags$br(),
                            actionButton("load_split_chord_lr",
                                         tagList(icon("sync"), " Load LR pairs"),
                                         class = "btn-white-cellchat btn-full-width")
                     )
                   ),
                   tags$br(),
                   uiOutput("split_chord_lr_selector")
            ),
            
            # --- Column 3: Gaps + Generate + Summary table ---
            column(4,
                   tags$h5("Plot Settings", style = "color: white; font-weight: bold;"),
                   fluidRow(
                     column(6,
                            sliderInput("split_chord_gap_inner", "Gap within group:",
                                        value = 1, min = 0, max = 5, step = 0.5)
                     ),
                     column(6,
                            sliderInput("split_chord_gap_outer", "Gap between groups:",
                                        value = 8, min = 2, max = 20, step = 1)
                     )
                   ),
                   actionButton("generate_split_chord",
                                tagList(icon("chart-pie"), " Generate"),
                                class = "btn-gradient-green btn-full-width"),
                   tags$hr(style = "border-color: rgba(255,255,255,0.3); margin: 12px 0;"),
                   tags$h5("LR Pair Summary", style = "color: white; font-weight: bold;"),
                   uiOutput("split_chord_lr_table_ui")
            ),
            
            # --- Column 4: Export ---
            column(2,
                   tags$h5("Export", style = "color: white; font-weight: bold;"),
                   selectInput("split_chord_format", "Format:",
                               choices = c("PNG" = "png", "PDF" = "pdf",
                                           "SVG" = "svg", "TIFF" = "tiff"),
                               selected = "png"),
                   numericInput("split_chord_export_width", "Width/plot (px):",
                                value = 1200, min = 600, max = 3000),
                   sliderInput("split_chord_export_dpi", "Resolution (DPI):",
                               value = 300, min = 150, max = 600, step = 50),
                   tags$br(),
                   downloadButton("download_split_chord", "Download",
                                  class = "btn-white-cellchat",
                                  style = "width: 100%;")
            )
          ),
          
          # --- Plot display ---
          tags$div(
            style = "margin-top: 20px; overflow-x: auto; width: 100%;",
            imageOutput("split_chord_plot_display", height = "auto")
          )
        )
      ),
      # ======================= SPATIAL TRANSCRIPTOMICS TABS (start) =======================
      ############################## Spatial / Load dataset ##############################
      tabItem(
        tabName = "load_spatial_dataset",
        
        tags$div(
          class = "header-magenta-cyan",
          tags$h2("Spatial Transcriptomics Analysis"),
          tags$p("Map gene expression across tissue architecture")
        ),
        
        tags$div(
          class = "info-box-magenta-cyan",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Spatial Transcriptomics",
                    tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Spatial transcriptomics combines gene expression with spatial location, revealing how gene expression varies across tissue architecture."),
            tags$h5(icon("star"), " Supported Technologies:"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("10x Visium:"), " High-resolution spatial profiling"),
              tags$li(tags$strong("Visium HD:"), " Higher resolution (2, 8 or 16 µm bins)")
            ),
            tags$h5(icon("check-square"), " Data Requirements:"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li("Expression matrix (genes × spots)"),
              tags$li("Spatial coordinates for each spot"),
              tags$li("Scale factors (JSON)"),
              tags$li("Optional: Tissue histology images")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Include high-resolution tissue images for better visualization", style = "margin: 0;")
            )
          )
        ),
        
        tags$div(
          class = "spatial-parameter-box",
          tags$h4(icon("cog"), " Dataset Parameters"),
          fluidRow(
            column(6,
                   tags$div(
                     class = "param-card-magenta",
                     tags$label("Species:"),
                     selectInput("spatial_species_choice", NULL,
                                 choices = c("Mouse" = "mouse", "Human" = "human", "Rat" = "rat"),
                                 selected = "mouse"),
                     tags$p(icon("info-circle"), " Select organism for gene detection")
                   )
            ),
            column(6,
                   tags$div(
                     class = "param-card-cyan",
                     tags$label("Visium HD bin size:"),
                     selectInput("spatial_hd_bin_size", NULL,
                                 choices = c("2 µm" = "2", "8 µm (recommended)" = "8", "16 µm" = "16"),
                                 selected = "8"),
                     tags$p(icon("info-circle"), " Bin folder loaded from a Visium HD ZIP (ignored for standard Visium). The technology, spot size and tissue image are read from the data.")
                   )
            )
          )
        ),
        
        tags$h4("Load Your Data", style = "color: #495057; margin-bottom: 15px; font-weight: bold;"),
        
        fluidRow(
          column(6,
                 tags$div(
                   class = "card-magenta",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("map-marked-alt")),
                     tags$h4("Load 10x Visium Data", style = "text-shadow: 2px 2px 4px rgba(0,0,0,0.2);")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     tags$strong("Required files in ZIP:"),
                     tags$ul(
                       tags$li("filtered_feature_bc_matrix.h5"),
                       tags$li("spatial/scalefactors_json.json"),
                       tags$li("spatial/tissue_positions (CSV or Parquet)"),
                       tags$li("spatial/tissue_hires_image.png (optional)"),
                       tags$li("spatial/tissue_lowres_image.png (optional)")
                     ),
                     tags$strong("Data structure:"),
                     tags$pre("your_data.zip/\n├── filtered_feature_bc_matrix.h5\n└── spatial/\n    ├── scalefactors_json.json\n    ├── tissue_positions.csv\n    └── tissue_hires_image.png"),
                     tags$div(
                       class = "spatial-card-highlight",
                       icon("magic"), tags$strong(" Auto-detects format (CSV or Parquet)")
                     )
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     fileInput('spatial_file', NULL,
                               accept = '.zip',
                               buttonLabel = tagList(icon("upload", style = "margin-right: 8px;"), "Load Visium Data"),
                               placeholder = ""),
                     tags$div(class = "spatial-card-footer-text", icon("clock"), " Processing: 3-7 minutes")
                   )
                 )
          ),
          
          column(6,
                 tags$div(
                   class = "card-cyan",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("database")),
                     tags$h4("Load Processed Spatial Object", style = "text-shadow: 2px 2px 4px rgba(0,0,0,0.2);")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     tags$strong("Perfect for:"),
                     tags$ul(
                       tags$li("Previously analyzed spatial data"),
                       tags$li("Shared spatial Seurat objects"),
                       tags$li("Continuing from checkpoint"),
                       tags$li("Skip initial processing steps")
                     ),
                     tags$strong("Requirements:"),
                     tags$ul(
                       tags$li("Spatial Seurat object (.rds)"),
                       tags$li("Contains spatial coordinates"),
                       tags$li("Normalized expression data"),
                       tags$li("Optional: Tissue images embedded")
                     ),
                     tags$div(
                       class = "spatial-card-highlight",
                       icon("zap"), tags$strong(" Instant loading with full spatial context")
                     )
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     fileInput("load_spatial_seurat_file", NULL,
                               accept = ".rds",
                               buttonLabel = tagList(icon("upload", style = "margin-right: 8px;"), "Load Spatial Object"),
                               placeholder = ""),
                     tags$div(class = "spatial-card-footer-text", icon("bolt"), " Loads in seconds")
                   )
                 )
          )
        ),
        tags$div(
          class = "spatial-image-upload-box",
          tags$h5(icon("image"), " Tissue Image"),
          fluidRow(
            column(6,
                   radioButtons("spatial_image_resolution", "Displayed tissue image",
                                choices = c("Low resolution" = "lowres", "High resolution" = "hires"),
                                selected = character(0), inline = TRUE),
                   tags$p("Both Space Ranger images are kept with the object; high resolution uses more memory while displayed"),
                   fileInput("upload_high_res_image", NULL,
                             accept = c(".png", ".jpg", ".jpeg", ".tiff"),
                             buttonLabel = tagList(icon("image"), " Browse Image"),
                             width = "100%"),
                   tags$p("Replace default tissue image with custom H&E staining (PNG, JPEG, TIFF)")
            ),
            column(6,
                   verbatimTextOutput("current_image_info", placeholder = TRUE)
            )
          )
        ),
        
        tags$div(
          style = "margin-top: 25px;",
          tags$h4("Dataset Information", style = "color: #495057; margin-bottom: 15px; font-weight: bold;"),
          fluidRow(
            column(4, infoBoxOutput("spatial_spots_count", width = 12)),
            column(4, infoBoxOutput("spatial_genes_count", width = 12)),
            column(4, infoBoxOutput("spatial_technology_info", width = 12))
          ),
          fluidRow(column(12, DTOutput("spatial_dataset_summary")))
        ),
        example_data_box("spatial")
      ),
      ############################## Spatial / Tissue sections ##############################
      tabItem(
        tabName = "spatial_sections",
        tags$div(
          class = "header-magenta-cyan",
          tags$h2("Tissue Sections"),
          tags$p("Separate the tissue sections placed on the same capture area and assign each one to a sample and a condition")
        ),
        tags$div(
          class = "info-box-magenta-cyan",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About tissue sections",
                    tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Spots closer than the maximum gap are connected; each connected group of spots is one section. Groups smaller than the minimum size are labelled Debris."),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Section:"), " one piece of tissue (spatial unit: neighbours are never taken across sections)"),
              tags$li(tags$strong("Sample:"), " the animal or donor the section comes from. Sections of one animal are not biological replicates. Sample and condition are stored with the spots (saved with the object) for sample-level analyses outside Cell-Hub"),              tags$li(tags$strong("Condition:"), " the group compared (e.g. WT, mdx)")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " A section split by a tear appears as two sections: give both the same name to merge them.", style = "margin: 0;")
            )
          )
        ),
        tags$h4("Detection", class = "spatial-section-title"),
        fluidRow(
          column(4,
                 tags$div(
                   class = "card-magenta",
                   tags$div(class = "spatial-card-header", tags$div(class = "card-icon-circle", icon("object-ungroup")), tags$h4("Parameters")),
                   tags$div(
                     class = "spatial-card-info",
                     tags$label("Maximum gap inside a section (grid pitches):", class = "spatial-parameter-label"),
                     numericInput("section_max_gap", NULL, value = 2, min = 1, max = 20, step = 0.5),
                     tags$p(class = "spatial-input-info", icon("info-circle"), " Spots further apart than this are in different sections. Increase it if a section is split in pieces"),
                     tags$label("Minimum section size (spots):", class = "spatial-parameter-label"),
                     numericInput("section_min_spots", NULL, value = 50, min = 1, max = 100000, step = 10),
                     tags$label("Colour the map by:", class = "spatial-parameter-label"),
                     radioButtons("section_plot_by", NULL, choices = c("Section" = "section", "Sample" = "sample", "Condition" = "condition"), selected = "section", inline = TRUE),
                     tags$label("Spot size:", class = "spatial-parameter-label"),
                     sliderInput("section_pt_size", NULL, min = 0.1, max = 4, value = 0.3, step = 0.1)
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     actionButton("detect_sections", tagList(icon("play"), " Detect sections"), class = "btn-white-magenta btn-full-width"),
                     actionButton("remove_debris_spots", tagList(icon("trash"), " Remove debris spots"), class = "btn-white-magenta btn-full-width"),
                     tags$p(icon("hand-pointer"), " Draw rectangles on the map to select spots (they add up)", style = "font-size: 12px; color: #6c757d; margin: 10px 0 5px 0;"),
                     actionButton("remove_selected_spots", tagList(icon("eraser"), " Remove selected spots"), class = "btn-white-magenta btn-full-width"),
                     actionButton("clear_section_selection", tagList(icon("rotate-left"), " Clear selection"), class = "btn-white-magenta btn-full-width")
                   )
                 )
          ),
          column(8,
                 tags$div(
                   class = "spatial-plot-box",
                   plotOutput("section_plot", height = "650px", brush = brushOpts(id = "section_brush", resetOnNew = TRUE, fill = "#000000", opacity = 0.15))                 )
          )
        ),
        tags$h4("Assignment", class = "spatial-section-title"),
        fluidRow(
          column(8,
                 tags$div(
                   class = "spatial-plot-box",
                   uiOutput("section_assignment_ui"),
                   actionButton("apply_sections", tagList(icon("check"), " Apply names, samples and conditions"), class = "btn-gradient-violet", style = "margin-top: 10px;"),
                   actionButton("remove_checked_sections", tagList(icon("trash"), " Remove ticked sections"), class = "btn-gradient-violet", style = "margin-top: 10px; margin-left: 10px;")
                   )
          ),
          column(4,
                 tags$div(
                   class = "spatial-plot-box",
                   tags$h5(icon("table"), " Sections", class = "spatial-plot-title"),
                   DTOutput("section_summary_table")
                 )
          )
        )
      ),
      ############################## Spatial / QC & normalization ##############################
      tabItem(
        tabName = "spatial_normalisation",
        
        tags$div(
          class = "header-magenta-cyan",
          tags$h2("Quality Control & Normalization"),
          tags$p("Filter and normalize your spatial transcriptomics data")
        ),
        
        tags$div(
          class = "info-box-magenta-cyan",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Spatial QC & Normalization",
                    tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Quality control for spatial data involves both standard scRNA-seq metrics and spatial-specific considerations."),
            tags$h5(icon("star"), " Key QC Metrics:"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Standard QC:"), " Number of genes, UMI counts, mitochondrial content per spot"),
              tags$li(tags$strong("Spatial Coverage:"), " Distribution of spots across tissue section"),
              tags$li(tags$strong("Tissue Detection:"), " Identification of spots over tissue vs. background"),
              tags$li(tags$strong("Spatial Continuity:"), " Expression continuity across neighboring spots")
            ),
            tags$h5(icon("check-square"), " Normalization Methods:"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Standard Normalization:"), " LogNormalize method - fast and robust"),
              tags$li(tags$strong("SCTransform:"), " Variance stabilizing transformation - handles technical noise better")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Spatial spots typically detect more genes and UMIs than single cells due to multiple cells per spot", style = "margin: 0;")
            )
          )
        ),
        
        fluidRow(
          column(4,
                 tags$div(
                   class = "card-magenta",
                   style = "min-height: 850px;",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("filter")),
                     tags$h4("QC Filtering")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     tags$div(
                       tags$label("Genes per Spot:", class = "spatial-parameter-label"),
                       sliderInput("spatial_nFeature_range", NULL,
                                   min = 0, max = 8000, value = c(50, 6000)),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Detect 200-5000 genes per spot")
                     ),
                     tags$div(
                       tags$label("UMI Counts per Spot:", class = "spatial-parameter-label"),
                       sliderInput("spatial_nCount_range", NULL,
                                   min = 0, max = 50000, value = c(50, 40000)),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Higher than single cells")
                     ),
                     tags$div(
                       tags$label("Maximum Mitochondrial %:", class = "spatial-parameter-label"),
                       sliderInput("spatial_mt_max", NULL,
                                   value = 20, min = 0, max = 100),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " May be higher in spatial data")
                     ),
                     checkboxInput("filter_tissue_spots", "Keep only spots over tissue", value = TRUE)
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     actionButton("spatial_qc_plots",
                                  tagList(icon("chart-bar"), " Generate QC Plots"),
                                  class = "btn-white-magenta btn-full-width"),
                     actionButton("apply_spatial_qc",
                                  tagList(icon("check-circle"), " Apply QC Filters"),
                                  class = "btn-white-magenta btn-full-width")
                   ),
                   tags$div(
                     class = "spatial-card-footer-text",
                     icon("info-circle"), " Review plots before applying"
                   )
                 )
          ),
          column(8,
                 tabsetPanel(
                   tabPanel("Spatial QC",
                            br(),
                            fluidRow(
                              column(6,
                                     tags$div(
                                       class = "spatial-plot-box",
                                       tags$h5(icon("map"), " Tissue Overview", class = "spatial-plot-title"),
                                       plotOutput("spatial_tissue_plot", height = "350px")
                                     )
                              ),
                              column(6,
                                     tags$div(
                                       class = "spatial-plot-box",
                                       tags$h5(icon("chart-bar"), " QC Metrics", class = "spatial-plot-title"),
                                       plotOutput("spatial_qc_violin", height = "350px")
                                     )
                              )
                            ),
                            br(),
                            fluidRow(
                              column(12,
                                     tags$div(
                                       class = "spatial-plot-box",
                                       tags$h5(icon("chart-line"), " Feature Scatter", class = "spatial-plot-title"),
                                       plotOutput("spatial_feature_scatter", height = "350px")
                                     )
                              )
                            )
                   ),
                   tabPanel("Standard QC",
                            br(),
                            fluidRow(
                              column(12,
                                     tags$div(
                                       class = "spatial-plot-box",
                                       tags$h5(icon("chart-area"), " Violin Plots", class = "spatial-plot-title"),
                                       plotOutput("spatial_vlnplot", height = "400px")
                                     )
                              )
                            ),
                            br(),
                            fluidRow(
                              column(6,
                                     tags$div(
                                       class = "spatial-plot-box",
                                       tags$h5(icon("braille"), " nFeature vs nCount", class = "spatial-plot-title"),
                                       plotOutput("spatial_scatter1", height = "300px")
                                     )
                              ),
                              column(6,
                                     tags$div(
                                       class = "spatial-plot-box",
                                       tags$h5(icon("braille"), " nCount vs Mito %", class = "spatial-plot-title"),
                                       plotOutput("spatial_scatter2", height = "300px")
                                     )
                              )
                            )
                   )
                 )
          )
        ),
        
        tags$div(style = "margin-top: 30px;"),
        fluidRow(
          column(4,
                 tags$div(
                   class = "card-cyan",
                   style = "min-height: 650px;",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("chart-line")),
                     tags$h4("Normalization")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     tags$div(
                       tags$label("Variable Features:", class = "spatial-parameter-label"),
                       numericInput("spatial_var_features", NULL,
                                    value = 3000, min = 1000, max = 8000, step = 500),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Highly variable genes to identify")
                     ),
                     tags$div(
                       tags$label("Normalization Method:", class = "spatial-parameter-label"),
                       shinyWidgets::sliderTextInput(
                         inputId = "norm_method",
                         label = NULL,
                         choices = c("Normalize", "SCT"),
                         selected = "Normalize",
                         grid = TRUE
                       ),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " LogNormalize (fast) or SCT (better)")
                     )
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     actionButton("normalize_spatial_data",
                                  tagList(icon("rocket"), " Normalize Spatial Data"),
                                  class = "btn-white-cyan btn-full-width")
                   )
                 )
          ),
          column(8,
                 tags$div(
                   class = "spatial-plot-box",
                   style = "min-height: 650px;",
                   tags$h5(icon("dna"), " Highly Variable Features", class = "spatial-plot-title"),
                   plotOutput("spatial_variable_features", height = "600px")
                 )
          )
        )
      ),
      
      ############################## Spatial / PCA & clustering ##############################
      tabItem(
        tabName = "spatial_clustering",
        
        tags$div(
          class = "header-magenta-cyan",
          tags$h2("Spatial Clustering"),
          tags$p("Identify spatial domains and cell populations across tissue")
        ),
        
        tags$div(
          class = "info-box-magenta-cyan",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Spatial Clustering",
                    tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Clustering spatial transcriptomics data identifies spatially coherent regions with similar expression profiles."),
            tags$h5(icon("star"), " Clustering Approaches:"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Expression-based:"), " Traditional clustering using gene expression similarity"),
              tags$li(tags$strong("Spatial-aware:"), " Methods considering spatial relationships between spots"),
              tags$li(tags$strong("Domain Detection:"), " Identification of spatially coherent tissue regions")
            ),
            tags$h5(icon("cogs"), " Workflow Steps:"),
            tags$ol(
              style = "line-height: 1.8;",
              tags$li(tags$strong("PCA:"), " Reduce dimensionality while preserving variance"),
              tags$li(tags$strong("Neighbors:"), " Calculate spot-spot similarities"),
              tags$li(tags$strong("UMAP:"), " Visualize expression space in 2D"),
              tags$li(tags$strong("Clustering:"), " Group spots into distinct populations")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Start with resolution 0.5, increase for more clusters or decrease for fewer", style = "margin: 0;")
            )
          )
        ),
        
        tags$div(
          style = "background: #e3f2fd; border-left: 4px solid #2196F3; padding: 20px; border-radius: 6px; margin-bottom: 25px;",
          tags$h5(icon("info-circle"), " Processing & Assay Status", style = "margin-top: 0; margin-bottom: 20px; color: #1976D2;"),
          fluidRow(
            column(3,
                   infoBoxOutput("spatial_processing_status", width = 12)
            ),
            column(3,
                   infoBoxOutput("spatial_assay_status", width = 12)
            ),
            column(3,
                   infoBoxOutput("spatial_current_assay_box", width = 12)
            ),
            column(3,
                   tags$div(
                     style = "background: white; padding: 15px; border-radius: 8px; height: 100%; display: flex; flex-direction: column; justify-content: center;",
                     tags$label("Default Assay:",
                                style = "margin-bottom: 8px; display: block; font-weight: bold; color: #495057; font-size: 13px;"),
                     selectInput("spatial_default_assay", NULL,
                                 choices = c("Spatial"), selected = "Spatial"),
                     actionButton("set_default_assay",
                                  tagList(icon("check"), " Set"),
                                  class = "btn-gradient-pink",
                                  style = "width: 100%; margin-top: 5px;")
                   )
            )
          )
        )
        ,
        fluidRow(
          column(4,
                 tags$div(
                   class = "card-magenta",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("project-diagram")),
                     tags$h4("PCA & Scaling")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     tags$div(
                       tags$label("Analysis Assay:", class = "spatial-parameter-label"),
                       selectInput("spatial_assay_select", NULL,
                                   choices = c("Spatial"), selected = "Spatial"),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Assay to use for analysis")
                     ),
                     tags$div(
                       tags$label("Number of PCA Dimensions:", class = "spatial-parameter-label"),
                       numericInput("spatial_pca_dims", NULL,
                                    value = 30, min = 10, max = 50),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Principal components to compute (bounded by the numbers of spots and genes)")
                     ),
                     tags$div(
                       tags$label("Dimensions for UMAP:", class = "spatial-parameter-label"),
                       numericInput("spatial_umap_dims", NULL,
                                    value = 30, min = 10, max = 50),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " PCs to use for UMAP")
                     )
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     actionButton("spatial_scale_pca",
                                  tagList(icon("calculator"), " Scale Data & Run PCA"),
                                  class = "btn-white-magenta btn-full-width")
                   ),
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.1); padding: 10px; border-radius: 6px; margin-top: 15px;",
                     tags$h6("Status:", style = "margin: 0 0 10px 0; font-weight: bold;"),
                     verbatimTextOutput("spatial_pca_status", placeholder = TRUE),
                     verbatimTextOutput("spatial_umap_status", placeholder = TRUE)
                   )
                 )
          ),
          column(8,
                 tags$div(
                   class = "spatial-plot-box",
                   tags$h5(icon("chart-line"), " Elbow Plot - PCA Variance", class = "spatial-plot-title"),
                   plotOutput("spatial_elbow_plot", height = "500px")
                 )
          )
        ),
        
        tags$div(style = "margin-top: 30px;"),
        fluidRow(
          column(4,
                 tags$div(
                   class = "card-cyan",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("circle-nodes")),
                     tags$h4("Clustering")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     tags$div(
                       tags$label("Clustering Dimensions:", class = "spatial-parameter-label"),
                       numericInput("spatial_cluster_dims", NULL,
                                    value = 30, min = 1, max = 50),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Number of PCs for clustering")
                     ),
                     tags$div(
                       tags$label("Clustering Resolution:", class = "spatial-parameter-label"),
                       numericInput("spatial_resolution", NULL,
                                    value = 0.5, min = 0.1, max = 150, step = 0.05),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Higher = more clusters")
                     ),
                     tags$div(
                       tags$label("Clustering Algorithm:", class = "spatial-parameter-label"),
                       selectInput("spatial_cluster_algorithm", NULL,
                                   choices = list("Louvain" = 1,
                                                  "Louvain (multilevel)" = 2,
                                                  "SLM" = 3),
                                   selected = 1),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Community detection method")
                     ),
                     tags$div(
                       checkboxInput("spatial_plot_labels", "Show cluster labels", value = TRUE)
                     )
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     actionButton("spatial_run_neighbors_umap",
                                  tagList(icon("network-wired"), " Run Neighbors + UMAP"),
                                  class = "btn-white-cyan btn-full-width"),
                     actionButton("spatial_find_clusters",
                                  tagList(icon("circle-nodes"), " Find Clusters"),
                                  class = "btn-white-cyan btn-full-width"),
                     downloadButton("save_spatial_object", "Save Spatial Object",
                                    class = "btn-white-cyan btn-full-width")
                   ),
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.1); padding: 10px; border-radius: 6px; margin-top: 15px;",
                     tags$h6("Status:", style = "margin: 0 0 10px 0; font-weight: bold;"),
                     verbatimTextOutput("spatial_neighbors_status", placeholder = TRUE),
                     verbatimTextOutput("spatial_clusters_status", placeholder = TRUE)
                   )
                 )
          ),
          column(8,
                 tags$div(
                   class = "spatial-plot-box",
                   tags$h5(icon("project-diagram"), " UMAP Clustering", class = "spatial-plot-title"),
                   plotOutput("spatial_umap_clusters", height = "450px")
                 ),
                 tags$div(
                   class = "spatial-plot-box",
                   style = "margin-top: 15px;",
                   tags$h5(icon("map"), " Spatial Clustering", class = "spatial-plot-title"),
                   plotOutput("spatial_tissue_clusters", height = "450px")
                 )
          )
        )
      ),
      
      ############################## Spatial / Interactive tissue viewer ##############################
      tabItem(
        tabName = "spatial_interactive_visualization",
        
        tags$div(
          class = "header-magenta-cyan",
          tags$h2("Interactive Tissue Viewer"),
          tags$p("Zoom into the tissue and explore clusters, genes and QC metrics spot by spot")
        ),
        
        tags$div(
          class = "info-box-magenta-cyan",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About the Tissue Viewer",
                    tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("The viewer uses the tissue image stored in the object (high resolution when available) and redraws the zoomed region on the server."),
            tags$h5(icon("star"), " How to use it:"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Zoom:"), " draw a rectangle on the tissue, or use the zoom buttons"),
              tags$li(tags$strong("Whole tissue:"), " double-click on the plot, or use the reset button"),
              tags$li(tags$strong("Navigate:"), " arrow buttons move the window by 40% of its size"),
              tags$li(tags$strong("Spot details:"), " click on a spot to see its barcode, cluster, counts and position"),
              tags$li(tags$strong("Coloring:"), " clusters (with highlight and palettes), gene expression or QC metrics")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Spots are drawn at their physical size, so they grow as you zoom in; the spot size slider scales them further", style = "margin: 0;")
            )
          )
        ),
        
        fluidRow(
          column(3,
                 tags$div(
                   class = "card-magenta",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("sliders-h")),
                     tags$h4("Display")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     checkboxInput("viewer_show_spots", "Show spots (untick to see the tissue image only)", value = TRUE),
                     conditionalPanel(
                       condition = "input.viewer_show_spots",
                       tags$label("Color spots by:", class = "spatial-parameter-label"),
                       selectInput("viewer_color_by", NULL,
                                   choices = c("Clusters" = "cluster", "Gene expression" = "gene", "Co-localization (last analysis)" = "coloc",
                                               "UMI count" = "nCount_Spatial", "Genes detected" = "nFeature_Spatial", "Mitochondrial %" = "percent.mt"),
                                   selected = "cluster"),
                       conditionalPanel(
                         condition = "input.viewer_color_by == 'coloc'",
                         checkboxInput("viewer_coloc_partners_only", "Only spots with a partner within the distance", value = FALSE),
                         tags$p(icon("info-circle"), " Shows the last two-gene analysis of the Co-localization tab. Dark blue / dark red: a spot of one gene with a spot of the other within the distance; purple: both genes in the same spot.",
                                style = "font-size: 11px; color: rgba(255,255,255,0.8);")
                       ),
                       conditionalPanel(
                         condition = "input.viewer_color_by == 'gene'",
                         selectizeInput("viewer_gene", "Gene:", choices = NULL, options = list(placeholder = "Type a gene name"))
                       ),
                       conditionalPanel(
                         condition = "input.viewer_color_by == 'cluster'",
                         selectInput("viewer_cluster_highlight", "Highlight cluster:", choices = "None", selected = "None"),
                         selectInput("viewer_palette", "Color palette:",
                                     choices = c("Default (Seurat)" = "default", "Polychrome" = "polychrome", "Set3" = "set3",
                                                 "Paired" = "paired", "Dark2" = "dark2"),
                                     selected = "default"),
                         checkboxInput("viewer_show_labels", "Show cluster labels", value = FALSE)
                       ),
                       sliderInput("viewer_spot_size", "Spot size:", min = 0.25, max = 4, value = 1, step = 0.25),
                       sliderInput("viewer_spot_alpha", "Spot opacity:", min = 0.1, max = 1, value = 0.8, step = 0.05)
                     ),
                     sliderInput("viewer_image_alpha", "Tissue image opacity:", min = 0, max = 1, value = 1, step = 0.05),
                     sliderInput("viewer_plot_height", "Viewer height (px):", min = 500, max = 1400, value = 750, step = 50),
                     checkboxInput("viewer_remove_legend", "Remove legend", value = FALSE)
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     fluidRow(
                       column(4, actionButton("viewer_zoom_in", icon("search-plus"), class = "btn-white-magenta btn-full-width")),
                       column(4, actionButton("viewer_zoom_out", icon("search-minus"), class = "btn-white-magenta btn-full-width")),
                       column(4, actionButton("viewer_reset", icon("expand"), class = "btn-white-magenta btn-full-width"))
                     ),
                     fluidRow(
                       column(3, actionButton("viewer_pan_left", icon("arrow-left"), class = "btn-white-magenta btn-full-width")),
                       column(3, actionButton("viewer_pan_up", icon("arrow-up"), class = "btn-white-magenta btn-full-width")),
                       column(3, actionButton("viewer_pan_down", icon("arrow-down"), class = "btn-white-magenta btn-full-width")),
                       column(3, actionButton("viewer_pan_right", icon("arrow-right"), class = "btn-white-magenta btn-full-width"))
                     ),
                     selectInput("viewer_download_format", "Export format:", choices = c("PNG" = "png", "PDF" = "pdf"), selected = "png"),
                     downloadButton("download_viewer_plot", "Download current view", class = "btn-white-magenta btn-full-width")
                   )
                 )
          ),
          column(9,
                 tags$div(
                   style = "background-color: white; padding: 15px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
                   verbatimTextOutput("viewer_view_info", placeholder = TRUE),
                   plotOutput("viewer_plot", height = "auto",
                              brush = brushOpts(id = "viewer_brush", resetOnNew = TRUE, direction = "xy"),
                              click = "viewer_click", dblclick = "viewer_dblclick"),
                   tags$h5(icon("crosshairs"), " Selected spot", style = "margin-top: 15px; color: #495057;"),
                   verbatimTextOutput("viewer_spot_info", placeholder = TRUE)
                 )
          )
        )
      ),
      
      ############################## Spatial / Gene expression visualization ##############################
      tabItem(
        tabName = "spatial_gene_expression_visualisation",
        
        tags$div(
          class = "header-magenta-cyan",
          tags$h2("Gene Expression Visualization"),
          tags$p("Explore gene expression patterns across spatial tissue")
        ),
        
        tags$div(
          class = "info-box-magenta-cyan",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Spatial Gene Expression",
                    tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Visualize gene expression patterns directly on tissue sections to understand spatial organization and cell-cell interactions."),
            tags$h5(icon("star"), " Visualization Types:"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Spatial Feature Plots:"), " Gene expression overlaid on tissue coordinates"),
              tags$li(tags$strong("Violin Plots:"), " Expression distribution across clusters"),
              tags$li(tags$strong("Dot Plots:"), " Expression level and percentage of expressing spots"),
              tags$li(tags$strong("Co-expression:"), " Identify genes with similar spatial patterns")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Start with known marker genes to validate clustering, then explore novel patterns", style = "margin: 0;")
            )
          )
        ),
        
        tags$div(
          style = "background: #e3f2fd; border-left: 4px solid #2196F3; padding: 20px; border-radius: 6px; margin-bottom: 25px;",
          tags$h5(icon("cog"), " Global Parameters", style = "margin-top: 0; margin-bottom: 20px; color: #1976D2;"),
          fluidRow(
            column(3,
                   tags$label("Select Genes (Global):", style = "font-weight: bold; color: #495057; margin-bottom: 8px;"),
                   pickerInput("spatial_genes_select", NULL,
                               choices = NULL, multiple = TRUE,
                               options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                   tags$p(style = "font-size: 12px; color: #666; margin-top: 5px;",
                          icon("info-circle"), " Selected genes will auto-fill in plot sections below")
            ),
            column(3,
                   tags$label("Assay:", style = "font-weight: bold; color: #495057; margin-bottom: 8px;"),
                   selectInput("spatial_viz_assay", NULL,
                               choices = c("Spatial"), selected = "Spatial"),
                   tags$p(style = "font-size: 12px; color: #666; margin-top: 5px;",
                          icon("info-circle"), " Data assay to visualize")
            ),
            column(3,
                   tags$label("Export DPI:", style = "font-weight: bold; color: #495057; margin-bottom: 8px;"),
                   numericInput("spatial_export_dpi", NULL,
                                value = 300, min = 72, step = 72),
                   tags$p(style = "font-size: 12px; color: #666; margin-top: 5px;",
                          icon("info-circle"), " Resolution for downloads")
            ),
            column(3,
                   tags$label("Dataset Info:", style = "font-weight: bold; color: #495057; margin-bottom: 8px;"),
                   verbatimTextOutput("spatial_assay_info", placeholder = TRUE)
            )
          )
        ),
        
        tags$h4("Dot Plot Analysis", class = "spatial-section-title"),
        fluidRow(
          column(4,
                 tags$div(
                   class = "card-magenta",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("circle")),
                     tags$h4("Dot Plot Settings")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     tags$div(
                       tags$label("Genes (separated by , ; or spaces):", class = "spatial-parameter-label"),
                       textInput("spatial_dotplot_genes_text", NULL,
                                 placeholder = "e.g., Myh4, Myh7, Tnnt3"),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Enter gene names separated by commas")
                     ),
                     tags$div(
                       checkboxInput("spatial_dot_scale", "Scale dot size", value = TRUE),
                       selectInput("spatial_dot_scale_by", "Scale by:",
                                   choices = c("radius", "size"), selected = "radius")
                     ),
                     tags$div(
                       tags$label("Dot Size Range:", class = "spatial-parameter-label"),
                       numericInput("spatial_dot_min", "Min:", value = 0, min = 0, max = 10),
                       numericInput("spatial_dot_max", "Max:", value = 6, min = 1, max = 20)
                     ),
                     tags$div(
                       checkboxInput("spatial_dot_cluster_idents", "Cluster on y-axis", value = FALSE),
                       selectInput("spatial_dot_colors", "Color palette:",
                                   choices = c("RdYlBu", "Blues", "Reds", "Greens", "Spectral"),
                                   selected = "RdYlBu")
                     )
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     actionButton("show_spatial_dotplot",
                                  tagList(icon("braille"), " Generate Dot Plot"),
                                  class = "btn-white-magenta btn-full-width"),
                     downloadButton("download_spatial_dotplot", "Download Dot Plot",
                                    class = "btn-white-magenta btn-full-width")
                   )
                 )
          ),
          column(8,
                 tags$div(
                   class = "spatial-plot-box expandable-plot-magenta",
                   tags$h5(icon("circle"), " Dot Plot - Expression by Cluster", class = "spatial-plot-title"),
                   plotOutput("spatial_dotplot")
                 )
          )
        ),
        
        tags$h4("Spatial Feature Plots", class = "spatial-section-title"),
        fluidRow(
          column(4,
                 tags$div(
                   class = "card-cyan",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("map-marked")),
                     tags$h4("Feature Plot Settings")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     tags$div(
                       tags$label("Genes for this plot:", class = "spatial-parameter-label"),
                       textInput("spatial_genes_text", NULL,
                                 placeholder = "Auto-filled or type manually"),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Uses global gene selection")
                     ),
                     tags$div(
                       tags$label("Image Size:", class = "spatial-parameter-label"),
                       sliderInput("spatial_image_size", NULL,
                                   min = 0.5, max = 2, value = 1, step = 0.1),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Scale tissue image")
                     ),
                     tags$div(
                       tags$label("Spot Size:", class = "spatial-parameter-label"),
                       numericInput("spatial_pt_size", NULL,
                                    value = 1.5, min = 0.1, max = 5, step = 0.1),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Size of expression spots")
                     ),
                     tags$div(
                       tags$label("Brightness:", class = "spatial-parameter-label"),
                       sliderInput("spatial_brightness", NULL,
                                   min = 0.5, max = 3, value = 1, step = 0.1),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Adjust expression intensity")
                     ),
                     checkboxInput("spatial_combine_plots", "Combine multiple genes", value = TRUE)
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     actionButton("show_spatial_features",
                                  tagList(icon("images"), " Generate Plots"),
                                  class = "btn-white-cyan btn-full-width"),
                     downloadButton("download_spatial_features", "Download",
                                    class = "btn-white-cyan btn-full-width")
                   )
                 )
          ),
          column(8,
                 tags$div(
                   class = "spatial-plot-box expandable-plot-cyan",
                   tags$h5(icon("map-marked"), " Spatial Feature Plots", class = "spatial-plot-title"),
                   plotOutput("spatial_feature_plots", height = "auto")
                 )
          )
        ),
        
        tags$h4("Violin Plot Analysis", class = "spatial-section-title"),
        fluidRow(
          column(4,
                 tags$div(
                   class = "card-magenta",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("chart-area")),
                     tags$h4("Violin Plot Settings")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     tags$div(
                       tags$label("Genes (separated by , ; or spaces):", class = "spatial-parameter-label"),
                       textInput("spatial_violin_genes_text", NULL,
                                 placeholder = "e.g., Myh4, Myh7, Tnnt3"),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Enter gene names separated by commas")
                     ),
                     tags$div(
                       checkboxInput("spatial_violin_points", "Show points", value = TRUE),
                       checkboxInput("spatial_violin_log", "Log scale", value = FALSE)
                     ),
                     tags$div(
                       tags$label("Cluster Order:", class = "spatial-parameter-label"),
                       selectInput("spatial_cluster_order", NULL,
                                   choices = NULL, multiple = TRUE),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Reorder clusters")
                     ),
                     tags$div(
                       tags$label("Point Size:", class = "spatial-parameter-label"),
                       numericInput("spatial_violin_pt_size", NULL,
                                    value = 0.1, min = 0, max = 2, step = 0.1),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Size of data points")
                     )
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     actionButton("show_spatial_violin",
                                  tagList(icon("chart-area"), " Generate Violin Plot"),
                                  class = "btn-white-magenta btn-full-width"),
                     downloadButton("download_spatial_violin", "Download",
                                    class = "btn-white-magenta btn-full-width")
                   )
                 )
          ),
          column(8,
                 tags$div(
                   class = "spatial-plot-box expandable-plot-magenta",
                   tags$h5(icon("chart-area"), " Violin Plot - Expression Distribution", class = "spatial-plot-title"),
                   plotOutput("spatial_violin_plot", height = "auto")
                 )
          )
        )
      ),
      
      ############################## Spatial / Cluster annotation ##############################
      tabItem(
        tabName = "spatial_cluster_annotation",
        
        tags$div(
          class = "header-magenta-cyan",
          tags$h2("Cluster Annotation"),
          tags$p("Assign biological identities to your spatial clusters")
        ),
        
        tags$div(
          class = "info-box-magenta-cyan",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Cluster Annotation",
                    tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Annotate your spatial clusters with biological identities based on marker gene expression patterns."),
            tags$h5(icon("star"), " Annotation Methods:"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Interactive Renaming:"), " Manually rename clusters based on marker analysis"),
              tags$li(tags$strong("Batch Renaming:"), " Upload CSV file with cluster names for quick annotation"),
              tags$li(tags$strong("Real-time Preview:"), " See changes instantly on UMAP and spatial plots"),
              tags$li(tags$strong("Export:"), " Save annotated object with new cluster identities")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Use marker genes from previous analysis to guide annotation decisions", style = "margin: 0;")
            )
          )
        ),
        
        tags$h4("Cluster Annotation", class = "spatial-section-title"),
        fluidRow(
          column(4,
                 tags$div(
                   class = "card-magenta",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("tags")),
                     tags$h4("Cluster Renaming")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     tags$label("Individual Cluster Names:", class = "spatial-parameter-label"),
                     tags$div(
                       style = "max-height: 500px; overflow-y: auto; padding: 10px;",
                       uiOutput("spatial_cluster_rename_inputs")
                     )
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     actionButton("apply_all_spatial_names",
                                  tagList(icon("check-circle"), " Apply All Names"),
                                  class = "btn-white-magenta btn-full-width"),
                     actionButton("reset_spatial_names",
                                  tagList(icon("undo"), " Reset to Original"),
                                  class = "btn-white-magenta btn-full-width")
                   ),
                   tags$div(
                     style = "background-color: rgba(255,255,255,0.1); padding: 15px; border-radius: 6px; margin-top: 15px;",
                     tags$label("Export Options:", class = "spatial-parameter-label"),
                     downloadButton("save_annotated_spatial_object",
                                    "Save Annotated Object",
                                    class = "btn-white-magenta btn-full-width")
                   )
                 )
          ),
          column(8,
                 tags$div(
                   class = "card-cyan",
                   tags$div(
                     class = "spatial-card-header",
                     tags$h4("Cluster Visualization", style = "text-shadow: 2px 2px 4px rgba(0,0,0,0.2);")
                   ),
                   tags$div(
                     style = "background-color: white; padding: 15px; border-radius: 8px; margin-bottom: 15px;",
                     fluidRow(
                       column(3,
                              checkboxInput("anno_dark_mode", "Dark mode", value = FALSE)
                       ),
                       column(3,
                              checkboxInput("spatial_remove_legend", "Remove Legend", value = FALSE)
                       ),
                       column(3,
                              checkboxInput("spatial_remove_axes", "Remove Axes", value = FALSE)
                       ),
                       column(3,
                              selectInput("spatial_plot_format", "Format:",
                                          choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg", "JPEG" = "jpeg", "TIFF" = "tiff"),
                                          selected = "png")
                       )
                     ),
                     fluidRow(
                       column(6,
                              downloadButton("download_spatial_umap", "Download UMAP",
                                             class = "btn-gradient-cyan",
                                             style = "width: 100%; margin-bottom: 5px;")
                       ),
                       column(6,
                              downloadButton("download_spatial_tissue", "Download Spatial",
                                             class = "btn-gradient-cyan",
                                             style = "width: 100%;")
                       )
                     )
                   ),
                   tabsetPanel(
                     tabPanel("UMAP View",
                              br(),
                              tags$div(
                                class = "spatial-plot-box",
                                plotOutput("spatial_annotation_umap", height = "600px")
                              )
                     ),
                     tabPanel("Spatial View",
                              br(),
                              fluidRow(
                                column(6,
                                       checkboxInput("spatial_anno_show_labels", "Show cluster labels", value = TRUE)
                                ),
                                column(6,
                                       selectInput("spatial_cluster_highlight",
                                                   "Highlight specific cluster:",
                                                   choices = "None",
                                                   selected = "None")
                                )
                              ),
                              br(),
                              tags$div(
                                class = "spatial-plot-box",
                                plotOutput("spatial_annotation_tissue", height = "600px")
                              )
                     ),
                     tabPanel("Gene Expression",
                              br(),
                              tags$div(
                                style = "background: #e3f2fd; padding: 15px; border-radius: 6px; border-left: 4px solid #2196F3; margin-bottom: 15px;",
                                tags$h5(icon("info-circle"), " View Gene Expression", style = "margin: 0 0 10px 0; color: #1976D2;"),
                                tags$p("Generate plots in the 'Gene Expression Visualization' tab first, then select which one to display here.",
                                       style = "margin: 0; font-size: 14px; color: #555;")
                              ),
                              fluidRow(
                                column(12,
                                       selectInput("annotation_selected_plot", "Select Plot Type:",
                                                   choices = c("None" = "none",
                                                               "Violin Plot" = "violin",
                                                               "Dot Plot" = "dot",
                                                               "Feature Plot" = "feature"),
                                                   selected = "none")
                                )
                              ),
                              tags$hr(),
                              div(
                                id = "annotation_plots_container",
                                div(
                                  id = "no_plot_panel",
                                  style = "text-align: center; padding: 50px; color: #666;",
                                  tags$h4("No plot selected"),
                                  tags$p("Choose a plot type above to display gene expression data")
                                ),
                                shinyjs::hidden(
                                  div(
                                    id = "violin_plot_panel",
                                    plotOutput("annotation_violin_display", height = "350px")
                                  )
                                ),
                                shinyjs::hidden(
                                  div(
                                    id = "dot_plot_panel",
                                    plotOutput("annotation_dot_display", height = "350px")
                                  )
                                ),
                                shinyjs::hidden(
                                  div(
                                    id = "feature_plot_panel",
                                    plotOutput("annotation_feature_display", height = "350px")
                                  )
                                )
                              )
                     )
                   )
                 )
          )
        )
      ),
      
      ############################## Spatial / Co-localization ##############################
      tabItem(
        tabName = "spatial_colocalization",
        
        tags$div(
          class = "header-magenta-cyan",
          tags$h2("Spatial Gene Co-localization"),
          tags$p("Analyze spatial proximity and interactions between genes")
        ),
        
        tags$div(
          class = "info-box-magenta-cyan",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Co-localization Analysis",
                    tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Co-localization analysis identifies genes that are expressed in close spatial proximity, revealing potential interactions and functional relationships."),
            tags$h5(icon("star"), " Key Features:"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Distance-based Counting:"), " Count spots of Gene A near spots of Gene B"),
              tags$li(tags$strong("Bidirectional Analysis:"), " Analyze both Gene A → Gene B and Gene B → Gene A"),
              tags$li(tags$strong("Cluster Information:"), " Track co-localization within and between clusters"),
              tags$li(tags$strong("Micrometer Scale:"), " Use real spatial distances based on technology resolution")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Start with 50μm distance for standard Visium, adjust based on tissue structure", style = "margin: 0;")
            )
          )
        ),
        
        tags$h4("Analysis Parameters", class = "spatial-section-title"),
        fluidRow(
          column(6,
                 tags$div(
                   class = "card-magenta",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("dna")),
                     tags$h4("Gene Selection")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     tags$div(
                       tags$label("Genes (separated by , ; or spaces):", class = "spatial-parameter-label"),
                       textInput("coloc_genes", NULL,
                                 value = "", placeholder = "e.g., Gene1, Gene2, Gene3"),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " 2 genes: bidirectional analysis. 3+ genes: simultaneous co-localization")
                     ),
                     tags$div(
                       tags$label("Assay:", class = "spatial-parameter-label"),
                       selectInput("coloc_assay", NULL,
                                   choices = c("Spatial"), selected = "Spatial"),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Data assay to use")
                     ),
                     checkboxInput("coloc_bidirectional", "Bidirectional analysis", value = TRUE)
                   )
                 )
          ),
          column(6,
                 tags$div(
                   class = "card-cyan",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("ruler")),
                     tags$h4("Distance Parameters")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     tags$div(
                       tags$label("Spot / bin diameter (µm):", class = "spatial-parameter-label"),
                       numericInput("coloc_spot_size_um", NULL, value = 55, min = 1, max = 100, step = 1),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " 55 µm for Visium, bin size for Visium HD (filled automatically when a ZIP or .rds is loaded)")
                     ),
                     tags$div(
                       style = "background-color: rgba(255,255,255,0.15); padding: 10px; border-radius: 4px; margin-bottom: 15px;",
                       verbatimTextOutput("detected_resolution_info", placeholder = TRUE)
                     ),
                     tags$div(
                       tags$label("Maximum Distance (μm):", class = "spatial-parameter-label"),
                       numericInput("coloc_distance_um", NULL,
                                    value = 100, min = 1, max = 500, step = 1),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Distance between spot / bin centres. Set at loading to the nearest neighbours: 100 µm for Visium, 2 bins for Visium HD. Large distances saturate for widely expressed genes")
                     ),
                     tags$div(
                       tags$label("Permutations (0 = no test):", class = "spatial-parameter-label"),
                       numericInput("coloc_n_permutations", NULL, value = 200, min = 0, max = 5000, step = 100),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Two genes: tests whether the genes are closer than expected from their abundance (second gene placed at random, seed 42). Proximity percentages are reported in all cases")                     )
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     actionButton("run_coloc_analysis",
                                  tagList(icon("play"), " Run Co-localization Analysis"),
                                  class = "btn-white-cyan btn-full-width"),
                     conditionalPanel(
                       condition = "output.coloc_results_ready",
                       downloadButton("download_coloc_results", "Download Results",
                                      class = "btn-white-cyan btn-full-width")
                     )
                   )
                 )
          )
        ),
        
        conditionalPanel(
          condition = "output.coloc_results_ready",
          
          tags$div(
            style = "background: #d4edda; border-left: 4px solid #28a745; padding: 15px; border-radius: 6px; margin: 25px 0;",
            tags$h5(icon("check-circle"), " Analysis Results", style = "margin-top: 0; color: #155724;"),
            verbatimTextOutput("coloc_quick_stats")
          ),
          
          tags$h4("Detailed Analysis", class = "spatial-section-title"),
          tabsetPanel(
            tabPanel("Cluster Analysis",
                     br(),
                     fluidRow(
                       column(12,
                              tags$div(
                                class = "spatial-plot-box",
                                tags$h5(icon("table"), " Co-localization per Cluster", class = "spatial-plot-title"),
                                DTOutput("coloc_cluster_table")
                              )
                       )
                     ),
                     br(),
                     fluidRow(
                       column(12,
                              tags$div(
                                class = "spatial-plot-box",
                                tags$h5(icon("exchange-alt"), " Inter-cluster Interactions", class = "spatial-plot-title"),
                                fluidRow(
                                  column(6,
                                         selectInput("inter_cluster_pairs_select",
                                                     "Select cluster pairs:",
                                                     choices = NULL,
                                                     selected = NULL,
                                                     multiple = TRUE),
                                         tags$p(class = "spatial-input-info",
                                                icon("info-circle"), " Only inter-cluster pairs shown")
                                  ),
                                  column(3,
                                         checkboxInput("inter_cluster_show_all",
                                                       "Show all inter-cluster pairs",
                                                       value = TRUE)
                                  ),
                                  column(3,
                                         actionButton("reset_cluster_filter",
                                                      tagList(icon("undo"), " Reset"),
                                                      class = "btn-gradient-bordeaux",
                                                      style = "margin-top: 25px; width: 100%;")
                                  )
                                ),
                                DTOutput("coloc_inter_cluster_table")
                              )
                       )
                     )
            ),
            tabPanel("Spatial Visualization",
                     br(),
                     fluidRow(
                       column(8,
                              tags$div(
                                class = "spatial-plot-box",
                                tags$h5(icon("map"), " Co-localization map", class = "spatial-plot-title"),
                                fluidRow(
                                  column(6, checkboxInput("coloc_map_partners_only", "Only spots with a partner within the distance", value = FALSE)),
                                  column(6, sliderInput("coloc_map_point_size", "Spot size:", min = 0.2, max = 3, value = 0.8, step = 0.1))
                                ),
                                plotOutput("coloc_map", height = "650px"),
                                tags$p(icon("info-circle"), " Spots expressing neither gene are not drawn. To zoom on a region, use the Interactive Tissue Viewer with 'Color spots by: Co-localization (last analysis)'.",
                                       style = "font-size: 12px; color: #666;"),
                                downloadButton("download_coloc_map", "Download map (PNG)", class = "btn-gradient-bordeaux")
                              )
                       ),
                       column(4,
                              tags$div(
                                class = "spatial-plot-box",
                                tags$h5(icon("chart-bar"), " Pairs per cluster", class = "spatial-plot-title"),
                                plotOutput("coloc_cluster_bar", height = "650px"),
                                downloadButton("download_coloc_bar", "Download chart (PNG)", class = "btn-gradient-bordeaux")
                              )
                       )
                     )
            )
          )
        )
      ),
      ############################## Spatial / Spatially variable genes ##############################
      tabItem(
        tabName = "spatial_svg",
        tags$div(
          class = "header-magenta-cyan",
          tags$h2("Spatially Variable Genes"),
          tags$p("Find genes whose expression is organized in space (Moran's I)")
        ),
        tags$div(
          class = "info-box-magenta-cyan",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About spatially variable genes",
                    tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Moran's I measures whether neighbouring spots have similar expression of a gene: values near 0 mean no spatial pattern, values towards 1 mean expression forms spatial domains."),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("Neighbour graph:"), " each spot is linked to its k nearest spots (6 matches the Visium hexagonal grid)"),
              tags$li(tags$strong("Test:"), " one-sided z-test under the normality assumption, Benjamini-Hochberg correction over the genes tested"),
              tags$li(tags$strong("Data:"), " normalized expression of the selected assay, on all spots (never a sketched subset)")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " With thousands of spots, almost every patterned gene is significant: rank by Moran's I rather than by p-value.", style = "margin: 0;")
            )
          )
        ),
        tags$h4("Parameters", class = "spatial-section-title"),
        fluidRow(
          column(4,
                 tags$div(
                   class = "card-magenta",
                   tags$div(class = "spatial-card-header", tags$div(class = "card-icon-circle", icon("dna")), tags$h4("Genes and graph")),
                   tags$div(
                     class = "spatial-card-info",
                     tags$label("Assay:", class = "spatial-parameter-label"),
                     selectInput("svg_assay", NULL, choices = c("Spatial"), selected = "Spatial"),
                     tags$label("Genes tested:", class = "spatial-parameter-label"),
                     radioButtons("svg_feature_source", NULL, choices = c("Variable features" = "variable", "Most expressed genes" = "expressed"), selected = "variable"),
                     tags$label("Maximum number of genes:", class = "spatial-parameter-label"),
                     numericInput("svg_n_features", NULL, value = 2000, min = 10, max = 20000, step = 100),
                     tags$label("Minimum % of spots expressing the gene:", class = "spatial-parameter-label"),
                     numericInput("svg_min_pct", NULL, value = 1, min = 0, max = 100, step = 1),
                     tags$label("Neighbours per spot (k):", class = "spatial-parameter-label"),
                     numericInput("svg_k", NULL, value = 6, min = 1, max = 50, step = 1),
                     tags$p(class = "spatial-input-info", icon("info-circle"), " 6 for Visium (hexagonal grid), 8 for Visium HD (square bins)")
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     actionButton("run_svg", tagList(icon("play"), " Find spatially variable genes"), class = "btn-white-magenta btn-full-width"),
                     downloadButton("download_svg_table", "Download table (CSV)", class = "btn-white-magenta btn-full-width")
                   )
                 )
          ),
          column(8,
                 tags$div(
                   class = "spatial-plot-box",
                   tags$h5(icon("table"), " Moran's I per gene (click a row to map the gene)", class = "spatial-plot-title"),
                   DTOutput("svg_table")
                 )
          )
        ),
        tags$h4("Spatial map", class = "spatial-section-title"),
        fluidRow(
          column(4,
                 tags$div(
                   class = "card-cyan",
                   tags$div(class = "spatial-card-header", tags$div(class = "card-icon-circle", icon("map")), tags$h4("Display")),
                   tags$div(
                     class = "spatial-card-info",
                     tags$label("Gene:", class = "spatial-parameter-label"),
                     selectizeInput("svg_plot_gene", NULL, choices = NULL),
                     tags$label("Spot size:", class = "spatial-parameter-label"),
                     sliderInput("svg_pt_size", NULL, min = 0.5, max = 4, value = 1.6, step = 0.1),
                     tags$label("Export format:", class = "spatial-parameter-label"),
                     selectInput("svg_plot_format", NULL, choices = c("png", "pdf", "tiff", "svg"), selected = "png")
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     downloadButton("download_svg_plot", "Download map", class = "btn-white-cyan btn-full-width")
                   )
                 )
          ),
          column(8,
                 tags$div(
                   class = "spatial-plot-box",
                   plotOutput("svg_plot", height = "600px")
                 )
          )
        )
      ),
      ############################## Spatial / Marker analysis ##############################
      tabItem(
        tabName = "spatial_marker_analysis",
        
        tags$div(
          class = "header-magenta-cyan",
          tags$h2("Marker Gene Analysis"),
          tags$p("Identify genes that characterize spatial clusters")
        ),
        
        tags$div(
          class = "info-box-magenta-cyan",
          tags$div(
            style = "cursor: pointer;",
            onclick = "$(this).next().slideToggle();",
            tags$h4(icon("info-circle"), " About Marker Analysis",
                    tags$span(icon("chevron-down"), style = "float: right; font-size: 14px;"))
          ),
          tags$div(
            style = "display: none; margin-top: 15px;",
            tags$p("Find differentially expressed genes that define and characterize your spatial clusters."),
            tags$h5(icon("star"), " Analysis Types:"),
            tags$ul(
              style = "line-height: 1.8;",
              tags$li(tags$strong("One vs All:"), " Compare one cluster against all other clusters combined"),
              tags$li(tags$strong("One vs Selected:"), " Compare one cluster against specific chosen clusters"),
              tags$li(tags$strong("Statistical Testing:"), " Wilcoxon rank-sum test for differential expression"),
              tags$li(tags$strong("Export Results:"), " Download marker gene tables as CSV")
            ),
            tags$div(
              class = "tip-box",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Use logFC > 0.5 and min.pct > 0.25 for stringent markers", style = "margin: 0;")
            )
          )
        ),
        
        tags$h4("Marker Analysis", class = "spatial-section-title"),
        fluidRow(
          column(4,
                 tags$div(
                   class = "card-magenta",
                   tags$div(
                     class = "spatial-card-header",
                     tags$div(class = "card-icon-circle", icon("dna")),
                     tags$h4("Analysis Parameters")
                   ),
                   tags$div(
                     class = "spatial-card-info",
                     tags$div(
                       tags$label("Comparison Type:", class = "spatial-parameter-label"),
                       radioButtons("spatial_comparison_type", NULL,
                                    choices = c("One vs All" = "one_vs_all",
                                                "One vs Selected" = "one_vs_selected"),
                                    selected = "one_vs_all"),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Choose comparison strategy")
                     ),
                     tags$div(
                       tags$label("Target Cluster (ident.1):", class = "spatial-parameter-label"),
                       selectInput("spatial_marker_cluster", NULL,
                                   choices = NULL),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Cluster to find markers for")
                     ),
                     conditionalPanel(
                       condition = "input.spatial_comparison_type == 'one_vs_selected'",
                       tags$div(
                         tags$label("Compare Against (ident.2):", class = "spatial-parameter-label"),
                         selectInput("spatial_comparison_clusters", NULL,
                                     choices = NULL,
                                     multiple = TRUE),
                         tags$p(class = "spatial-input-info",
                                icon("info-circle"), " Select comparison clusters")
                       )
                     ),
                     tags$div(
                       tags$label("Min Log2 Fold Change:", class = "spatial-parameter-label"),
                       numericInput("spatial_marker_logfc", NULL,
                                    value = 0.25, min = 0, max = 2, step = 0.1),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Minimum expression difference")
                     ),
                     tags$div(
                       tags$label("Min % Spots Expressing:", class = "spatial-parameter-label"),
                       numericInput("spatial_marker_min_pct", NULL,
                                    value = 0.1, min = 0, max = 1, step = 0.05),
                       tags$p(class = "spatial-input-info",
                              icon("info-circle"), " Detection threshold")
                     ),
                     checkboxInput("spatial_marker_only_pos", "Only positive markers", value = TRUE)
                   ),
                   tags$div(
                     class = "spatial-card-footer",
                     actionButton("run_spatial_markers",
                                  tagList(icon("play"), " Find Markers"),
                                  class = "btn-white-magenta btn-full-width"),
                     downloadButton("download_spatial_markers", "Download Results",
                                    class = "btn-white-magenta btn-full-width")
                   )
                 )
          ),
          column(8,
                 tags$div(
                   class = "spatial-plot-box",
                   tags$h5(icon("table"), textOutput("marker_results_title", inline = TRUE), class = "spatial-plot-title"),
                   DTOutput("spatial_marker_results")
                 )
          )
        )
      ),
      # ======================== SPATIAL TRANSCRIPTOMICS TABS (end) ========================
      
      ############################# Trajectory/Monocle Conversion and trajectory ##############################
      
      
      ########################################Trajectory Tab 1 / Monocle Conversion####################################
      tabItem(
        tabName = "load_trajectory",
        
        tags$div(
          class = "header-blue-gradient",
          tags$h2("Monocle 3 Setup & Conversion"),
          tags$p("Convert Seurat object to Monocle format for trajectory analysis")
        ),
        
        tags$div(
          class = "info-box-light-blue",
          tags$div(
            class = "info-toggle",
            onclick = "$(this).next().slideToggle();",
            tags$h4(
              icon("info-circle"), " About Monocle 3 Trajectory Analysis",
              tags$span(icon("chevron-down"), class = "chevron")
            )
          ),
          tags$div(
            class = "info-content-hidden",
            tags$p("Monocle 3 is a toolkit for analyzing single-cell RNA-seq data to reconstruct cellular trajectories and understand developmental processes."),
            tags$h5(icon("list-check"), " Key Concepts:"),
            tags$ul(
              tags$li(tags$strong("Trajectory:"), " A path through gene expression space that represents cellular differentiation or state transitions"),
              tags$li(tags$strong("Pseudotime:"), " A measure of progress along the trajectory, representing developmental or temporal progression"),
              tags$li(tags$strong("Principal Graph:"), " The structure that connects cells along the trajectory, revealing branching and decision points")
            ),
            tags$h5(icon("sitemap"), " Analysis Workflow:"),
            tags$ol(
              tags$li(tags$strong("Conversion:"), " Convert your Seurat object to Monocle 3 format while preserving UMAP and cluster information"),
              tags$li(tags$strong("Graph Construction:"), " Build the trajectory graph using principal curves through your data"),
              tags$li(tags$strong("Root Selection:"), " Define the starting point (root) of your trajectory"),
              tags$li(tags$strong("Pseudotime Calculation:"), " Order cells along the trajectory from the root"),
              tags$li(tags$strong("Gene Analysis:"), " Identify genes that change along the trajectory")
            ),
            tags$div(
              class = "alert-box-yellow",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Ensure your Seurat object has UMAP coordinates and clear cluster assignments"),
              tags$p(icon("chart-line"), tags$strong(" Note:"), " Conversion preserves all metadata, reductions (PCA/UMAP), and cluster identities"),
              tags$p(icon("dna"), tags$strong(" Important:"), " Choose a root cluster that represents your starting cell population (e.g., stem cells, progenitors)")
            ),
            tags$p(
              class = "docs-link",
              "For more information, see: ",
              tags$a(
                href = "https://cole-trapnell-lab.github.io/monocle3/docs/introduction/",
                "Monocle 3 Documentation", 
                target = "_blank"
              )
            )
          )
        ),
        
        tags$div(
          class = "box-blue-gradient",
          style = "margin-top: 30px;",
          tags$div(
            class = "box-header",
            tags$h3(icon("upload"), " Load & Convert Seurat Object")
          ),
          tags$div(
            class = "box-body",
            # White inner block — Seurat upload only
            tags$div(
              style = "background: white; border-radius: 8px; padding: 25px; margin-bottom: 25px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
              tags$div(
                style = "margin-bottom: 15px;",
                tags$h4(icon("file-upload"), " Load Data",
                        style = "color: #2c3e50; margin: 0 0 10px 0;"),
                tags$p("A Seurat object (with UMAP and clusters) is converted to Monocle3 automatically; a Monocle object saved earlier is restored with its graph and pseudotime.",
                       style = "color: #7f8c8d; margin: 0;")
              ),
              tags$div(
                class = "trajectory-import-row",
                style = "display: flex; align-items: center; flex-wrap: wrap; gap: 15px; margin-bottom: 20px; padding-bottom: 15px; border-bottom: 1px solid #ecf0f1;",
                # Remove Shiny's default form-group bottom margin so the radios align with the button
                tags$style(HTML(".trajectory-import-row .form-group { margin-bottom: 0; }")),
                radioButtons("trajectory_import_source", NULL,
                             choices = c("Single dataset" = "single", "Multiple datasets" = "multiple"),
                             selected = "single", inline = TRUE),
                actionButton("import_session_object_trajectory", tagList(icon("arrow-right"), " Use current object"),
                             class = "btn-gradient-blue", style = "padding: 8px 20px; font-size: 14px;"),
                tags$span(icon("info-circle"), " Snapshot taken at click time, or upload a file below.",
                          style = "font-size: 12px; color: #7f8c8d;")
              ),
              fileInput("load_trajectory_rds", NULL,
                        accept = ".rds", buttonLabel = "Browse...",
                        placeholder = "Seurat object or saved Monocle object (.rds)",
                        width = "100%")
            ),
            # Save/Load Monocle block — outside the white div
            tags$div(
              style = "background: white; border-radius: 8px; padding: 20px; margin-bottom: 25px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);",
              tags$h4(icon("save"), " Save Monocle Object",
                      style = "color: #2c3e50; margin: 0 0 15px 0;"),
              tags$p("Save the Monocle3 object with trajectory graph and pseudotime; reload it later with the upload above.",
                     style = "color: #7f8c8d; margin: 0 0 15px 0;"),
              downloadButton(
                "download_monocle_object",
                tagList(icon("download"), " Save Monocle Object (.rds)"),
                class = "btn-gradient-blue",
                style = "width: 100%;"
              )
            ),
            
            tags$div(
              style = "margin-top: 10px;",
              tags$h4(icon("database"), " Dataset Information", 
                      style = "color: #2c3e50; margin-bottom: 20px;"),
              fluidRow(
                column(6,
                       tags$div(
                         class = "info-card",
                         style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                                   border-radius: 10px; padding: 20px; color: white; 
                                   box-shadow: 0 4px 6px rgba(0,0,0,0.1); min-height: 120px;",
                         tags$div(
                           style = "display: flex; align-items: center; margin-bottom: 10px;",
                           tags$div(
                             style = "font-size: 40px; margin-right: 15px;",
                             icon("dna")
                           ),
                           tags$div(
                             tags$h3("Cells", style = "margin: 0; font-weight: 600;"),
                             tags$h2(textOutput("monocle_cells_count", inline = TRUE), 
                                     style = "margin: 5px 0 0 0; font-weight: bold;")
                           )
                         )
                       )
                ),
                column(6,
                       tags$div(
                         class = "info-card",
                         style = "background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%); 
                                   border-radius: 10px; padding: 20px; color: white; 
                                   box-shadow: 0 4px 6px rgba(0,0,0,0.1); min-height: 120px;",
                         tags$div(
                           style = "display: flex; align-items: center; margin-bottom: 10px;",
                           tags$div(
                             style = "font-size: 40px; margin-right: 15px;",
                             icon("layer-group")
                           ),
                           tags$div(
                             tags$h3("Clusters", style = "margin: 0; font-weight: 600;"),
                             tags$h2(textOutput("monocle_clusters_names", inline = TRUE), 
                                     style = "margin: 5px 0 0 0; font-weight: bold;")
                           )
                         )
                       )
                )
              )
            )
          )
        ),
        # Demonstration data instructions (PBMC v3 processed in the Single module, then "Use current object")
        example_data_box("trajectory")
      ),
      
      ######################################## Trajectory Tab 2 / Trajectory Construction ####################################
      tabItem(
        tabName = "trajectory_construction",
        
        tags$div(
          class = "header-blue-gradient",
          tags$h2("Trajectory Construction & Root Selection"),
          tags$p("Build trajectory graph and calculate pseudotime ordering")
        ),
        
        tags$div(
          class = "info-box-light-blue",
          tags$div(
            class = "info-toggle",
            onclick = "$(this).next().slideToggle();",
            tags$h4(
              icon("project-diagram"), " About Trajectory Construction",
              tags$span(icon("chevron-down"), class = "chevron")
            )
          ),
          tags$div(
            class = "info-content-hidden",
            tags$p("Construct the trajectory graph and define pseudotime ordering by selecting a biological starting point."),
            tags$h5(icon("list-check"), " Construction Steps:"),
            tags$ul(
              tags$li(tags$strong("Graph Construction:"), " Monocle learns a principal graph through your UMAP space that connects similar cells and reveals branching paths"),
              tags$li(tags$strong("Root Selection:"), " Choose the cluster representing your starting cell population (e.g., stem cells, early progenitors)"),
              tags$li(tags$strong("Pseudotime Calculation:"), " Cells are ordered along the graph starting from the root, creating a developmental timeline")
            ),
            tags$h5(icon("lightbulb"), " Choosing the Right Root:"),
            tags$ul(
              tags$li("Select the cluster with the earliest/most undifferentiated cells"),
              tags$li("Consider biological markers (e.g., stem cell markers, low differentiation markers)"),
              tags$li("The root defines pseudotime=0; cells progress away from this point")
            ),
            tags$div(
              class = "alert-box-yellow",
              tags$p(icon("chart-line"), tags$strong(" Tip:"), " The trajectory plot shows the graph structure overlaid on clusters"),
              tags$p(icon("map"), tags$strong(" Note:"), " Pseudotime values increase as cells move away from the root along the trajectory"),
              tags$p(icon("recycle"), tags$strong(" Flexibility:"), " You can recalculate pseudotime with different roots to explore alternative trajectories")
            ),
            tags$p(
              class = "docs-link",
              "For more information, see: ",
              tags$a(
                href = "https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/",
                "Monocle3 - Learning Trajectories", 
                target = "_blank"
              )
            )
          )
        ),
        
        tags$div(
          class = "box-blue-gradient",
          style = "margin-top: 30px;",
          tags$div(
            class = "box-header",
            tags$h3(icon("project-diagram"), " Trajectory Graph Construction")
          ),
          tags$div(
            class = "box-body",
            tags$p(style = "margin-bottom: 15px; color: #666;", 
                   "Build the principal graph that connects cells along the trajectory and visualize it on your UMAP."),
            tags$div(
              class = "control-panel-gray",
              fluidRow(
                column(3, actionButton("constructGraph", tagList(icon("sitemap"), " Construct Graph"), class = "btn-gradient-blue", style = "width: 100%;")),
                column(3, selectInput("trajectory_download_format_main", "File Format:", choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg", "TIFF" = "tiff"), selected = "png")),
                column(3, numericInput("trajectory_download_dpi", "Resolution (DPI):", value = 300, min = 72, max = 1200, step = 10)),
                column(3, downloadButton("download_trajectory_umap", "Download Plot", class = "btn-gradient-blue", style = "width: 100%; margin-top: 25px;"))
              )
            ),
            tags$hr(class = "separator-line"),
            tags$div(class = "plot-output-container", plotOutput("trajectoryPlot", height = "600px"))
          )
        ),
        
        tags$div(
          class = "box-blue-gradient",
          style = "margin-top: 30px;",
          tags$div(
            class = "box-header",
            tags$h3(icon("seedling"), " Root Selection & Pseudotime Distribution")
          ),
          tags$div(
            class = "box-body",
            tags$p(style = "margin-bottom: 15px; color: #666;", 
                   "Select the starting cluster (root) and calculate pseudotime ordering. The pseudotime distribution plot shows how cells are ordered along the trajectory."),
            tags$div(
              class = "control-panel-gray",
              fluidRow(
                column(3, selectInput("root_cluster_select", "Select Root Cluster:", choices = NULL)),
                column(3, actionButton("set_root_cell", tagList(icon("play-circle"), " Calculate Pseudotime"), class = "btn-gradient-blue", style = "width: 100%; margin-top: 25px;")),
                column(2, selectInput("pseudotime_download_format", "File Format:", choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg", "TIFF" = "tiff"), selected = "png")),
                column(2, numericInput("pseudotime_download_dpi", "Resolution (DPI):", value = 300, min = 72, max = 1200, step = 10)),
                column(2, downloadButton("download_pseudotime_umap", "Download Plot", class = "btn-gradient-blue", style = "width: 100%; margin-top: 25px;"))
              ),
              fluidRow(
                column(4, checkboxInput("pseudotime_hide_grid", "Hide grid", value = FALSE)),
                column(4, checkboxInput("pseudotime_hide_axes", "Hide axes (titles and numbers)", value = FALSE)),
                column(4, checkboxInput("pseudotime_hide_title", "Hide plot title", value = FALSE))
              )
            ),
            tags$hr(class = "separator-line"),
            verbatimTextOutput("root_info", placeholder = TRUE),
            tags$hr(class = "separator-line"),
            tags$div(class = "plot-output-container", plotOutput("pseudotimePlot", height = "600px"))
          )
        )
      ),
      #############################Trajectory Tab 3/Differentialy expressed genes#######################
      tabItem(
        tabName = "differentialy_expressed_genes_trajectory",
        tags$div(
          class = "header-blue-gradient",
          tags$h2("Pseudotime Gene Expression Analysis"),
          tags$p("Identify and visualize genes that change along developmental trajectories")
        ),
        tags$div(
          class = "info-box-light-blue",
          tags$div(
            class = "info-toggle",
            onclick = "$(this).next().slideToggle();",
            tags$h4(
              icon("dna"), " About Pseudotime Gene Expression Analysis",
              tags$span(icon("chevron-down"), class = "chevron")
            )
          ),
          tags$div(
            class = "info-content-hidden",
            tags$p("This analysis identifies genes that change their expression patterns along the trajectory, revealing the molecular progression of cells through biological processes."),
            tags$h5(icon("list-check"), " Analysis Steps:"),
            tags$ul(
              tags$li(tags$strong("Differential Gene Test:"), " Identifies genes that significantly change expression over pseudotime"),
              tags$li(tags$strong("Expression Visualization:"), " Plots expression patterns of selected genes along the trajectory"),
              tags$li(tags$strong("Result Export:"), " Download differential expression results for further analysis")
            ),
            tags$div(
              class = "alert-box-yellow",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Focus on genes with low q-values for the most reliable results"),
              tags$p(icon("info-circle"), tags$strong(" Note:"), " The q-value cutoff controls the false discovery rate in your analysis"),
              tags$p(icon("chart-line"), tags$strong(" Hint:"), " Visualize multiple genes together to identify co-regulated patterns")
            ),
            tags$p(
              class = "docs-link",
              "For more information, see: ",
              tags$a(
                href = "http://cole-trapnell-lab.github.io/monocle-release/docs/#differentialgetest-details-and-options",
                "Monocle Documentation - Differential Expression Testing", 
                target = "_blank"
              )
            )
          )
        ),
        tags$div(
          class = "box-blue-gradient",
          style = "margin-top: 20px;",
          tags$div(
            class = "box-header",
            tags$h3(icon("project-diagram"), " Trajectory Structure Analysis (graph_test)")
          ),
          tags$div(
            class = "box-body",
            tags$p(style = "margin-bottom: 15px; color: #666;", 
                   "Identifies genes with spatial patterns in the trajectory. This analysis is ", 
                   tags$strong("independent of root cell selection"), 
                   " and shows genes that vary across the trajectory structure."),
            tags$div(
              class = "control-panel-gray",
              fluidRow(
                column(4, numericInput("sig_genes_cutoff", "q-value cutoff:", value = 0.05, min = 0, max = 1, step = 0.01)),
                column(4, actionButton("run_graph_test", tagList(icon("play"), " Run Graph Test"), class = "btn-gradient-blue", style = "margin-top: 25px; width: 100%;")),
                column(4, downloadButton("download_graph_test", "Download Results", class = "btn-gradient-blue", style = "margin-top: 25px; width: 100%;"))
              )
            ),
            tags$hr(class = "separator-line"),
            tags$div(class = "data-table-container", DTOutput("graphTestTable"))
          )
        ),
        tags$div(
          class = "box-green-gradient",
          style = "margin-top: 20px;",
          tags$div(
            class = "box-header",
            tags$h3(icon("chart-line"), " Pseudotime Correlation Analysis")
          ),
          tags$div(
            class = "box-body",
            
            
            tags$p(style = "margin-bottom: 15px; color: #666;",
                   "Identifies genes that increase (", 
                   tags$span("Up", style = "color: #d9534f; font-weight: bold;"),
                   ") or decrease (", 
                   tags$span("Down", style = "color: #5bc0de; font-weight: bold;"),
                   ") along the selected pseudotime. ",
                   tags$strong("Results depend on the root cell selected.")),
            selectInput("select_pseudotime", "Select Pseudotime to Use:", choices = NULL),
            
            tags$div(
              class = "control-panel-gray",
              fluidRow(
                column(4, numericInput("sig_genes_cutoff_corr", "q-value cutoff:", value = 0.05, min = 0, max = 1, step = 0.01)),
                column(4, actionButton("run_correlation_test", tagList(icon("play"), " Run Correlation Test"), class = "btn-gradient-green", style = "margin-top: 25px; width: 100%;")),
                column(4, downloadButton("download_correlation", "Download Results", class = "btn-gradient-green", style = "margin-top: 25px; width: 100%;"))
              )
            ),
            tags$hr(class = "separator-line"),
            tags$div(class = "data-table-container", DTOutput("correlationTable"))
          )
        )
      ),
      
      
      ########################################Trajectory Tab 4 / Genes visualization####################################
      tabItem(
        tabName = "genes_visualization_trajectory",
        tags$div(
          class = "header-blue-gradient",
          tags$h2("Gene Expression Visualization"),
          tags$p("Visualize gene expression patterns along pseudotime trajectory")
        ),
        
        tags$div(
          class = "info-box-light-blue",
          tags$div(
            class = "info-toggle",
            onclick = "$(this).next().slideToggle();",
            tags$h4(
              icon("dna"), " About Gene Expression Visualization",
              tags$span(icon("chevron-down"), class = "chevron")
            )
          ),
          tags$div(
            class = "info-content-hidden",
            tags$p("Visualize how genes change their expression along the pseudotime trajectory to understand molecular dynamics during cellular differentiation or development."),
            tags$h5(icon("chart-line"), " Visualization Options:"),
            tags$ul(
              tags$li(tags$strong("Single Gene on Trajectory:"), " Overlay a single gene's expression directly on the UMAP trajectory to visualize spatial expression patterns"),
              tags$li(tags$strong("Multi-Gene Expression Plot:"), " Plot up to 20 genes simultaneously to see expression trends along pseudotime with smooth trend lines")
            ),
            tags$div(
              class = "alert-box-yellow",
              tags$p(icon("lightbulb"), tags$strong(" Tip:"), " Select genes from the differential expression results in the previous tab for best results"),
              tags$p(icon("chart-area"), tags$strong(" Note:"), " Adjust point size and minimum expression threshold to optimize visualization clarity"),
              tags$p(icon("dna"), tags$strong(" Hint:"), " Use single-gene plot to see spatial organization, and multi-gene plot to identify co-expression patterns")
            ),
            tags$p(
              class = "docs-link",
              "For more information, see: ",
              tags$a(
                href = "http://cole-trapnell-lab.github.io/monocle-release/docs/#visualizing-genes-that-vary-over-pseudotime",
                "Monocle Documentation - Visualizing Gene Expression", 
                target = "_blank"
              )
            )
          )
        ),
        tags$div(
          class = "box-blue-gradient",
          style = "margin-top: 30px;",
          tags$div(
            class = "box-header",
            tags$h3(icon("map"), " Single Gene Expression on Trajectory")
          ),
          tags$div(
            class = "box-body",
            tags$p(style = "margin-bottom: 15px; color: #666;", 
                   "Overlay expression of a single gene on the UMAP trajectory to visualize spatial patterns and cell type specificity."),
            tags$div(
              class = "control-panel-gray",
              fluidRow(
                column(3, selectInput("single_gene_umap_trajectory", "Select a Gene:", choices = NULL)),
                column(3,
                       numericInput("single_gene_cell_size_trajectory", "Cell size:", value = 0.5, min = 0.1, max = 2, step = 0.1),
                       checkboxInput("show_trajectory_graph_single", "Show trajectory graph", value = TRUE)
                ),
                column(3,
                       selectInput("download_format_trajectory", "File Format:", 
                                   choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg", "TIFF" = "tiff"), selected = "pdf"),
                       numericInput("dpi_trajectory", "Resolution (DPI):", value = 300, min = 72, max = 600, step = 72)
                ),
                column(3,
                       tags$div(class = "button-group-aligned",
                                actionButton("plot_single_gene_umap_trajectory", tagList(icon("map"), " Plot on Trajectory"),
                                             class = "btn-gradient-blue", style = "width: 100%; margin-bottom: 10px;"),
                                downloadButton("download_single_gene_umap_trajectory", "Download Plot",
                                               class = "btn-gradient-blue", style = "width: 100%;"))
                )
              )
            ),
            tags$hr(class = "separator-line"),
            tags$div(class = "plot-output-container", plotOutput("singleGeneUmapTrajectoryPlot", height = "600px"))
          )
        ),
        # Box 2: Multi-Gene Expression Along Pseudotime
        tags$div(
          class = "box-blue-gradient",
          style = "margin-top: 30px;",
          tags$div(
            class = "box-header",
            tags$h3(icon("chart-line"), " Multi-Gene Expression Along Pseudotime")
          ),
          tags$div(
            class = "box-body",
            tags$p(style = "margin-bottom: 15px; color: #666;", 
                   "Plot multiple genes (maximum 20) to visualize expression dynamics along pseudotime with trend lines."),
            tags$div(
              class = "control-panel-gray",
              fluidRow(
                column(4,
                       pickerInput("multi_gene_pseudotime_trajectory", "Select Genes:", choices = NULL, multiple = TRUE,
                                   options = list(`actions-box` = TRUE, `live-search` = TRUE, 
                                                  `selected-text-format` = "count > 3", `max-options` = 20,
                                                  title = "Search and select genes"))
                ),
                column(4,
                       pickerInput("clusters_filter_trajectory", "Filter by Clusters/Datasets:", 
                                   choices = NULL, multiple = TRUE,
                                   options = list(`actions-box` = TRUE, `live-search` = TRUE,
                                                  `selected-text-format` = "count > 3",
                                                  title = "All clusters/datasets"))
                ),
                column(4,
                       numericInput("multi_gene_cell_size_trajectory", "Point size:", value = 0.75, min = 0.1, max = 3, step = 0.1),
                       numericInput("multi_gene_min_expr_trajectory", "Min expression:", value = 0.1, min = 0, max = 5, step = 0.1),
                       numericInput("trend_line_size_trajectory", "Trend line width:", value = 1, min = 0.5, max = 3, step = 0.25)
                )
              ),
              fluidRow(
                column(12,
                       tags$div(class = "button-group-aligned", style = "margin-top: 10px;",
                                actionButton("plot_multi_gene_pseudotime_trajectory", tagList(icon("chart-area"), " Generate Plot"),
                                             class = "btn-gradient-blue", style = "width: 200px; margin-right: 10px;"),
                                downloadButton("download_multi_gene_pseudotime_trajectory", "Download Plot",
                                               class = "btn-gradient-blue", style = "width: 200px;"))
                )
              )
            ),
            tags$hr(class = "separator-line"),
            tags$div(class = "plot-output-container", plotOutput("multiGenePseudotimeTrajectoryPlot", height = "600px"))
          )
          
        )),
      
      ############################## Acknowlegment and Licence ##############################
      tabItem(tabName = "acknowledgement",
              div(class = "container-fluid", style = "padding: 0;",
                  div(class = "well", style = "background: linear-gradient(rgba(0,0,0,0.7), rgba(0,0,0,0.7)), url('muscle.png');
                                   background-size: cover;
                                   background-position: center;
                                   color: white;
                                   padding: 40px;
                                   min-height: 800px;",
                      h1("Acknowledgements & License", class = "text-center", style = "font-size: 36px; margin-bottom: 40px;"),
                      div(class = "row",
                          div(class = "col-md-10 col-md-offset-1",
                              div(class = "section", style = "margin-bottom: 30px;",
                                  h3("Development", style = "color: #7ACFB0;"),
                                  p("This application was developed by Gaspard Macaux and is the property of the Neuromuscular Development, Genetics and Physiopathology laboratory directed by Dr. Pascal Maire.")
                              ),
                              div(class = "section", style = "margin-bottom: 30px;",
                                  h3("Contributors", style = "color: #7ACFB0;"),
                                  p("Special thanks to:"),
                                  tags$ul(
                                    tags$li("Edgar Jauliac, Léa Delivry and Hugues Escoffier for their expertise in transcriptomic analysis"),
                                    tags$li("Maxime Di Gallo  & Valentina Taglietti for testing the application")
                                  )
                              ),
                              div(class = "section", style = "margin-bottom: 30px;",
                                  h3("Technologies", style = "color: #7ACFB0;"),
                                  tags$ul(
                                    tags$li("Seurat - Comprehensive single-cell analysis toolkit"),
                                    tags$li("Shiny - Web application framework"),
                                    tags$li("R Studio - Development environment"),
                                    tags$li("Cell Chat - Ligand-Receptor analysis"),
                                    tags$li("Monocle - Trajectory analysis")
                                  )
                              ),
                              div(class = "section", style = "margin-top: 40px;",
                                  h3("License", style = "color: #7ACFB0;"),
                                  p("Cell-Hub is free software released under the GNU General Public License v3.0 (GPL-3.0), without any warranty."),
                                  tags$a(href = "https://www.gnu.org/licenses/gpl-3.0.html", "Read the GPL-3.0 license", style = "color: #7ACFB0;")
                              )
                          )
                      )
                  )
              )
      ) # Fin du dernier tabItem (acknowledgement)
    ), # Fermeture de tabItems
    # Script pour les popovers
    tags$script(HTML('$(function () { $("[data-toggle=\'popover\']").popover(); });'))
  ) # Fermeture de dashboardBody
)
