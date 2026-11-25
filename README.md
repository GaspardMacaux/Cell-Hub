# Cell-Hub: a graphical interface for end-to-end single-cell RNA sequencing analysis

**Cell-Hub** is a free, open-source graphical platform for single-cell and single-nucleus RNA-seq analysis, from raw data import to publication-ready figures, with no programming required.

Built with R/Shiny and distributed as a Docker image, Cell-Hub brings together Seurat v5, CellChat 2 and Monocle 3 in one interface that runs on your own computer, saves and reloads complete analyses, and records the parameters applied.

---

## Installation

Cell-Hub is distributed as a Docker image. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/) (Windows, macOS, Linux), then:

```
docker pull gaspardmacaux/cell-hub:2.0.0
docker run -p 3838:3838 gaspardmacaux/cell-hub:2.0.0
```

Open http://localhost:3838 in your browser. Large datasets need memory: raise the memory limit of Docker Desktop (Settings → Resources) if the analysis stops.

**New to Cell-Hub?** The [tutorial](TUTORIAL.md) runs every module on public demonstration data ([doi:10.5281/zenodo.24198709](https://doi.org/10.5281/zenodo.24198709)) in about an hour and a half.

Memory needed, measured on subsets of the 10x 1.3 million mouse brain cells dataset (Supplementary Data 2 of the article): about 16 GB for a sample of 10,000 cells, 32 GB for a single-sample analysis of 100,000 cells, 48 GB to integrate 100,000 cells; 200,000 cells were analysed completely in 83 min with a peak of 64 GB. DoubletFinder is the most demanding step (60 GB for 50,000 cells) and is run sample by sample.

The image used in the Cell-Hub article is `gaspardmacaux/cell-hub:2.0.0`, digest `sha256:51142acfefe3bdecc0f421c67dccdbeba826c98612529533a6f187e2475bb29b`; to be sure to get exactly this image, pull it by digest: `docker pull gaspardmacaux/cell-hub@sha256:51142acfefe3bdecc0f421c67dccdbeba826c98612529533a6f187e2475bb29b`.

The image contains the application, every R package at a fixed version and CellHub Database 2.1.2. The source code of the image is in [`Code/`](Code) and its build recipe in [`docker/`](docker).

### Package versions in `gaspardmacaux/cell-hub:2.0.0`

| Component | Version |
|---|---|
| R | 4.4.2 |
| Shiny | 1.14.0 |
| Seurat / SeuratObject | 5.5.1 / 5.4.0 |
| CellChat | 2.1.2 |
| Monocle 3 | 1.3.7 |
| Harmony | 2.0.5 |
| batchelor (fastMNN) | 1.22.0 |
| DoubletFinder | 2.0.6 |
| SingleR | 2.8.0 |
| presto | 1.1.0 |
| CellHub Database | 2.1.2 |

Every stochastic step (PCA, neighbours, clustering, UMAP, Harmony, CellChat permutations, co-localization permutations) uses seed 42, so the same data and parameters give the same results.

---

## Features

### Single-dataset analysis
- **Data import**: 10x Cell Ranger folders (ZIP), Cell Ranger HDF5 (`.h5`), AnnData (`.h5ad`), Seurat objects (`.rds`)
- **Quality control**: filtering on genes per cell and mitochondrial percentage, with interactive threshold plots
- **Normalization**: LogNormalize, SCTransform, CLR, RC
- **Dimensionality reduction**: PCA with elbow plot, UMAP in 2D and 3D
- **Clustering**: Louvain, Louvain with multilevel refinement, SLM, adjustable resolution
- **Doublet detection**: DoubletFinder with automatic pK selection
- **Cell type annotation**: manual renaming and merging of clusters, SingleR, label transfer from a reference

### Differential expression and biomarkers
- **Differential expression**: cluster vs all others, pairwise clusters, conditions within clusters, pseudobulk tests across samples
- **Exclusive biomarkers**: specificity score S = ((pct_target − pct_other) / 100) × log2FC, detection thresholds, Wilcoxon test with Benjamini-Hochberg correction
- **Pathway enrichment**: over-representation analysis and GSEA
- **Visualization**: volcano plots, Venn diagrams, dot plots, heatmaps, violin, ridge and feature plots
- **Co-expression**: cells co-expressing several genes, per cluster, with differential expression of co-expressing vs non-co-expressing cells

### Multi-dataset integration
- **Methods**: Seurat CCA anchors, Seurat reciprocal PCA (RPCA), fastMNN, Harmony, or a simple merge without batch correction
- **Metadata**: assign any experimental variable (condition, genotype, time point) to each dataset
- **Comparisons**: every visualization split by any metadata variable; cluster composition across conditions

### Cell-cell communication
- **Inference**: CellChat 2 with [CellHub Database](https://github.com/GaspardMacaux/CellHub-Database) (76,582 interactions, 21 species, 10 source databases)
- **Confidence tiers**: high (2 or more curated sources), medium (1 curated source; for human and mouse, optionally restricted to confirmed ligand / receptor roles), low (only orthology-inferred, predicted or computationally aggregated support, such as NicheNet v2); high + medium by default for human and mouse
- **Database selection**: merged database or chosen source databases, filtered by confidence tier
- **Gene lists**: named lists of genes of interest, to list their interactions in the database and to restrict bubble and chord plots to them
- **Visualization**: bubble plots, gene-level chord diagrams, split chord diagrams to compare conditions

### Trajectory inference
- **Monocle 3**: trajectory graph, pseudotime ordering from a root chosen on the UMAP
- **Gene dynamics**: genes varying along the trajectory (graph test) or correlated with pseudotime, expression along pseudotime

### Spatial transcriptomics (10x Visium and Visium HD)
- **Loading**: Space Ranger outputs (ZIP) or Seurat objects; tissue sections detected and assigned to samples and conditions
- **Analysis**: normalization, clustering, markers, spatially variable genes (Moran's I) on all spots (no subsampling)
- **Interactive tissue viewer**: zoom and pan on the tissue image, spots coloured by cluster, gene, QC metric or co-localization
- **Gene co-localization**: distance-based co-localization of two or more genes in micrometres, per cluster and between clusters, with a permutation test for two genes

### Outputs and reproducibility
- **Figures**: PNG, TIFF, PDF and SVG at chosen size and resolution
- **Tables**: CSV and Excel for every result
- **Reports**: PDF analysis report and a parameter table listing every setting applied
- **Sessions**: save a complete analysis and reload it later

---

## CellHub Database

Cell-Hub includes [CellHub Database 2.1.2](https://github.com/GaspardMacaux/CellHub-Database), a ligand-receptor resource integrating 10 curated databases into 76,582 interactions across 21 species, in CellChat format. The release, its build pipeline and the per-species files are described in the CellHub-Database repository and archived on Zenodo ([doi:10.5281/zenodo.24197570](https://doi.org/10.5281/zenodo.24197570), all versions; always resolves to the latest).

| Source | Species | Reference |
|---|---|---|
| CellChatDB v2 | human, mouse, zebrafish | Jin et al., 2024 |
| CellPhoneDB v5 | human | Troulé et al., 2025 |
| NicheNet v2 (MultiNicheNet) | human, mouse | Browaeys et al., 2023 |
| CellTalkDB | human, mouse | Shao et al., 2021 |
| connectomeDB2025 | 14 vertebrates | Liu et al., 2026 |
| NeuronChatDB | human, mouse | Zhao et al., 2023 |
| FlyPhoneDB2 | Drosophila | Qadiri et al., 2025 |
| cell2cell | *C. elegans* | Armingol et al., 2022 |
| PlantPhoneDB | Arabidopsis | Xu et al., 2022 |
| PlantCellChatDB | 5 plant species | Liu et al., 2026 |

To run the code outside Docker, place the per-species folders of the release (`app_databases/`) in `Code/databases/` (see [`Code/databases/README.md`](Code/databases/README.md)).

---

## Repository layout

```
Code/        application source (app.R, *_server.R, functions/, www/)
docker/      Dockerfile and R profile used to build the image
tests/       test suites run on every image before release (interface, functional, database)
```

In the image, the same files are installed with a `_docker` suffix (for example `functions/spatial_functions.R` becomes `functions/spatial_functions_docker.R`); their content is identical, except `main_server_docker.R`, which sets the paths used inside the container. To inspect exactly what runs in a given image, open it by digest (`docker run --rm -it --entrypoint bash gaspardmacaux/cell-hub@sha256:...`, files in `/srv/shiny-server/mon_app`).

The Dockerfile documents how the image was built (base image, system libraries, package versions, installation order). It also copies local package sources (`CellHub.tar.gz`, `monocle3`, `seurat-wrappers`, `seurat-disk`) that are not in this repository, so use the published image rather than rebuilding it.

`tests/run_cellhub_tests.sh` runs the three suites on an image with the demonstration data (doi:10.5281/zenodo.24198709): `bash tests/run_cellhub_tests.sh gaspardmacaux/cell-hub:2.0.0 /path/to/cellhub_demo_data`. It writes an HTML report with every step, screenshot and download.

---

## Citation

If you use Cell-Hub, please cite:

> Macaux G, Di Gallo M, Taglietti V, Amthor H, Maire P. *Cell-Hub: a graphical interface for end-to-end single-cell RNA sequencing analysis.* (2026).

If you use CellHub Database, please also cite:

> Macaux G, Maire P. *CellHub Database: a consolidated multi-species ligand-receptor interaction resource for cell-cell communication analysis.* bioRxiv (2026). doi:10.64898/2026.09.23.753685

and the source databases listed above.

---

## License

Cell-Hub is distributed under the [GNU General Public License v3.0](LICENSE).

## Contact

Gaspard Macaux — Institut Cochin, INSERM U1016, Université Paris Cité, Paris, France
