# Cell-Hub

**Cell-Hub** is a free, open-source graphical platform for single-cell and single-nucleus RNA-seq analysis. It provides a complete end-to-end workflow — from raw data import to publication-ready figures — with no programming required.

Built with R/Shiny and distributed as a Docker image, Cell-Hub integrates Seurat v5, CellChat, and Monocle 3 within a unified interface accessible to researchers without bioinformatics expertise.

> Associated publication: *Cell-Hub: the democratization of single-cell analysis* — submitted to Nucleic Acids Research (2026)

---

## Features

- **Single-cell analysis** — data loading (10X, H5, RDS), QC, normalization, PCA, clustering, UMAP (2D/3D), doublet detection
- **Multi-dataset integration** — Seurat CCA, Harmony, or simple merge; flexible metadata management
- **Differential expression** — cluster vs. all, pairwise, condition vs. condition; volcano plots, Venn diagrams, exclusive biomarker detection
- **Cell-cell communication** — CellChat-based inference powered by **GaspouDB**, a consolidated ligand-receptor database (11,563 mouse / 9,604 human interactions)
- **Trajectory inference** — Monocle 3 pseudotime, branch analysis, gene dynamics along trajectories
- **Spatial transcriptomics** — 10X Visium and Visium HD support
- **Export** — publication-quality figures (PNG, TIFF, PDF, SVG), CSV/Excel result tables, full session save/reload

---

## Installation

### Requirements

- [Docker Desktop](https://www.docker.com/products/docker-desktop/) (Windows 10/11, macOS 10.15+, or Linux)
- 8 GB RAM minimum (16+ GB recommended for datasets > 50,000 cells)
- 10 GB free disk space

### Pull the image

Open Docker Desktop, go to **Images → Search**, and enter:

```
gaspardmacaux/cell-hub:latest
```

Click **Pull**.

### Run Cell-Hub

1. Go to **Images** in Docker Desktop
2. Click **Run** next to `gaspardmacaux/cell-hub`
3. Open **Optional settings** and set:
   - Host port: `3838`
   - Container port: `3838`
4. Click **Run**
5. Open your browser and go to: [http://localhost:3838](http://localhost:3838)

### Subsequent launches

Once the container has been created, go to **Containers**, find `cell-hub`, and click **Start**. No need to reconfigure ports.

### Memory allocation

For large datasets, increase Docker's memory limit:  
**Docker Desktop → Settings → Resources → Advanced → Memory**  
Recommended: 16 GB for standard use, 32+ GB for large datasets.

---

## GaspouDB

Cell-Hub includes **GaspouDB**, a consolidated ligand-receptor interaction database built from:

- [CellChat](https://github.com/jinworks/CellChat)
- [CellPhoneDB](https://www.cellphonedb.org/)
- [CellTalkDB](http://tcm.zju.edu.cn/celltalkdb/)
- [MultiNicheNet](https://github.com/saeyslab/multinichenetr)

| Species | Interactions |
|---------|-------------|
| Mouse   | 11,563      |
| Human   | 9,604       |

---

## Citation

If you use Cell-Hub in your research, please cite:

> Macaux G, Di Gallo M, Taglietti V, Amthor H, Maire P. *Cell-Hub: the democratization of single-cell analysis.* Nucleic Acids Research (2025) — in revision.

---

## License

Cell-Hub is distributed under the [MIT License](LICENSE).

---

## Contact

Gaspard Macaux — [gaspard.macaux@gmail.com](mailto:gaspard.macaux@gmail.com)  
Institut Cochin, INSERM U1016, Paris, France
