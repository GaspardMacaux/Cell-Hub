# Cell-Hub : a graphical interface for end-to-end single-cell RNA sequencing analysis

**Cell-Hub** is a free, open-source graphical platform for single-cell and single-nucleus RNA-seq analysis. It provides a complete end-to-end workflow — from raw data import to publication-ready figures — with no programming required.

Built with R/Shiny and distributed as a Docker image, Cell-Hub integrates Seurat v5, CellChat 2, and Monocle 3 within a unified interface accessible to researchers without bioinformatics expertise.

---

## Features

### Single-dataset analysis
- **Data import** — 10X Cell Ranger (filtered/raw), H5, H5AD (Scanpy), Seurat RDS
- **Quality control** — filtering on gene count, UMI count, and mitochondrial percentage with interactive threshold visualization
- **Normalization** — LogNormalize, SCTransform, CLR, RC
- **Dimensionality reduction** — PCA with elbow plot, UMAP 2D and 3D
- **Clustering** — Louvain, multilevel, and SLM algorithms with adjustable resolution
- **Doublet detection** — DoubletFinder integration with automatic optimal pK selection
- **Cluster annotation** — interactive renaming and merging of clusters

### Differential expression & biomarker discovery
- **DEG analysis** — cluster vs. all, pairwise cluster comparison, condition vs. condition, within-cluster cross-condition
- **Exclusive biomarker detection** — custom algorithm with specificity score S = ((pct_target − pct_other) / 100) × log2FC, dual detection thresholds, and Wilcoxon + Benjamini-Hochberg correction
- **Visualization** — volcano plots, Venn diagrams, dot plots, heatmaps, violin plots, ridge plots, feature plots
- **Co-expression analysis** — multi-gene co-expression detection with cluster breakdown and within-cluster differential expression

### Multi-dataset integration
- **Integration methods** — Seurat CCA, RPCA, Harmony, or simple merge without batch correction
- **Metadata management** — flexible group labeling to map any experimental variable (condition, genotype, time point) across datasets
- **Comparative analysis** — all visualization types stratified by any metadata variable for side-by-side comparison
- **Cell composition** — cluster proportion analysis across conditions (stacked bar charts, pie charts)

### Cell-cell communication
- **Inference engine** — CellChat 2 with [**CellHub Database**](https://github.com/GaspardMacaux/CellHub-Database) (67,584 interactions, 21 species, 10 source databases)
- **Visualization** — bubble plots, chord diagrams, split chord diagrams for multi-condition comparison
- **Network analysis** — signaling pathway-level and individual ligand-receptor-level exploration
- **Species support** — user-selectable species with merged or individual source database loading

### Trajectory inference
- **Monocle 3 integration** — pseudotime ordering, trajectory graph construction
- **Root selection** — interactive root node selection on UMAP
- **Gene dynamics** — gene expression changes along pseudotime trajectories

### Spatial transcriptomics
- **Platform support** — 10X Visium and Visium HD
- **Spatial visualization** — feature plots, spatial clustering
- **Co-localization analysis** — inter-cluster ligand-receptor co-localization with distance measurement in µm, bidirectional analysis

### Output & reporting
- **Figures** — publication-quality export in PNG, TIFF, PDF, SVG with customizable dimensions
- **Tables** — CSV and Excel export for all results (DEGs, markers, interactions)
- **Automated PDF reports** — complete analysis reports with figures, tables, and parameters
- **Workspace** — full session save and reload for reproducibility

---

## Installation

Cell-Hub is distributed as a Docker image. Pull and run with [Docker Desktop](https://www.docker.com/products/docker-desktop/) (Windows, macOS, Linux):

```
docker pull gaspardmacaux/cell-hub:latest
docker run -p 3838:3838 gaspardmacaux/cell-hub:latest
```

Then open http://localhost:3838 in your browser.

---

## CellHub Database

Cell-Hub includes [**CellHub Database**](https://github.com/GaspardMacaux/CellHub-Database), a consolidated ligand-receptor interaction resource integrating 10 independently curated databases into **67,584 curated interactions** across **21 species** in CellChat-native format.

### Source databases

| Source | Species | Reference |
|---|---|---|
| CellChatDB v2 | Human, Mouse, Zebrafish | Jin et al., *Nat. Protoc.* 2025 |
| CellPhoneDB v5 | Human | Troulé et al., *Nat. Protoc.* 2025 |
| NicheNet v2 | Human, Mouse | Browaeys et al., *Nat. Methods* 2020 |
| CellTalkDB | Human, Mouse | Shao et al., *Brief. Bioinform.* 2021 |
| connectomeDB2025 | 14 vertebrates | Liu et al., *Nucleic Acids Res.* 2026 |
| NeuronChatDB | Human, Mouse | Zhao et al., *Nat. Commun.* 2023 |
| FlyPhoneDB2 | Drosophila | Qadiri et al., *Comput. Struct. Biotechnol. J.* 2025 |
| cell2cell | C. elegans | Armingol et al., *PLoS Comput. Biol.* 2022 |
| PlantPhoneDB | Arabidopsis | Xu et al., *Plant Biotechnol. J.* 2022 |
| PlantCellChatDB | 5 plant species | Liu et al., *Plant J.* 2026 |

For standalone usage, pipeline scripts, and per-species files, see the [CellHub-Database repository](https://github.com/GaspardMacaux/CellHub-Database).

---

## Citation

If you use Cell-Hub in your research, please cite:

> Macaux G, Di Gallo M, Taglietti V, Amthor H, Maire P. *Cell-Hub: a graphical interface for end-to-end single-cell RNA sequencing analysis.* (2026).

If you use CellHub Database, please also cite:

> Macaux G, Maire P. *CellHub Database: a consolidated multi-species ligand-receptor interaction resource for cell-cell communication analysis.* (2026).

---

## License

Cell-Hub is distributed under the [MIT License](LICENSE).

---

## Contact

Gaspard Macaux — Institut Cochin, INSERM U1016, Université Paris Cité, Paris, France
