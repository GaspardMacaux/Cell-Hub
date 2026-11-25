# Cell-Hub tutorial

This tutorial goes through every module of Cell-Hub 2.0.0 with public demonstration data, from loading the files to exporting figures and tables. Allow about one hour and a half for the whole tutorial; each section can also be followed on its own. Names in **bold** are the menus, tabs, buttons and fields exactly as they appear in the application.

Screenshots: [to add from the final image; one per step marked 📷].

---

## 0. Before you start

### 0.1 Install and start Cell-Hub

Cell-Hub runs in Docker, so nothing else needs to be installed.

1. Install Docker Desktop (https://www.docker.com) and start it.
2. In Docker Desktop, **Settings → Resources**, give Docker enough memory. As a guide: 16 GB for a sample of 10,000 cells (doublet detection included), 32 GB for a single-sample analysis of 100,000 cells, 48 GB to integrate 100,000 cells. The demonstration data run with 8 GB.
3. In a terminal:

```
docker pull gaspardmacaux/cell-hub:2.0.0
docker run -p 3838:3838 gaspardmacaux/cell-hub:2.0.0
```

4. Open http://localhost:3838 in a web browser. The data stay on your computer: nothing is uploaded to a server.

The image is built for Intel/AMD computers; on Apple silicon Macs, Docker Desktop runs it through emulation, which is slower. To stop Cell-Hub, press Ctrl+C in the terminal. If port 3838 is already used, run `docker run -p 3839:3838 ...` and open http://localhost:3839.

### 0.2 Download the demonstration data

All files are on Zenodo: [doi:10.5281/zenodo.24198709](https://doi.org/10.5281/zenodo.24198709). They are public 10x Genomics datasets, redistributed unchanged.

| File | Content | Used in |
|---|---|---|
| `pbmc_1k_v3_filtered_feature_bc_matrix.h5` | About 1,000 human blood cells, 10x v3 chemistry | Single dataset, CellChat, Trajectory, Multiple datasets |
| `pbmc_1k_v2_filtered_feature_bc_matrix.h5` | Same donor, 10x v2 chemistry | Multiple datasets |
| `visium_mouse_brain.zip` | Visium section of adult mouse brain (Space Ranger output) | Spatial |

The loading page of each module has a **Try the demonstration data** box that names the file to use.

### 0.3 How the interface is organised

The menu on the left lists the modules: **Single Dataset Analysis**, **Multiple Datasets Analysis**, **CellChat**, **Trajectory** and **Transcriptomic Spatial**. Click a module to open its steps and go through them from top to bottom.

- Each page has a short **About** box (click its title to open it) and a log panel that reports what was done.
- Errors and missing settings are shown as notifications in the bottom-right corner.
- Every result has a **Download** button next to it. It stays greyed out until the result exists. Figures are exported in PNG, TIFF, PDF, SVG or JPEG at the chosen resolution, tables in CSV or Excel.
- Every random step uses seed 42, so the same data and parameters give the same results.

### 0.4 Which module should I use?

| You have | Use | What it does |
|---|---|---|
| One sample (one Cell Ranger output) | **Single Dataset Analysis** | QC, clustering, annotation, markers, enrichment and figures for that sample alone. |
| Several samples to analyse together (conditions, time points, replicates, batches) | **Multiple Datasets Analysis** | Integrates 2 or more samples, then runs the same analyses on the combined object, with plots and tests split by sample or condition, and sample-level (pseudobulk) tests. |
| An analysed object, and a question about cell–cell communication | **CellChat** | Ligand–receptor inference with CellHub Database (21 species). |
| An analysed object, and a differentiation process | **Trajectory** | Monocle 3 trajectory and pseudotime. |
| Visium or Visium HD data | **Transcriptomic Spatial** | Spot- or bin-level analysis on the tissue image. |

**CellChat** and **Trajectory** do not load raw data. They take an object analysed in Single or Multiple, either with **Use current object** (no file needed) or from a saved `.rds` file.

### 0.5 Recommended workflow with several samples: each sample alone first, then together

1. Open each sample in **Single Dataset Analysis**. Choose its QC thresholds from its own plots, since samples rarely have the same quality. Remove doublets with **Doublet Detection**: DoubletFinder must run on each sample separately, never on merged samples. Check the clusters and marker genes.
2. Save each cleaned sample with **Save & Session**.
3. In **Multiple Datasets Analysis**, load the saved files with the type **Seurat Object**. The QC thresholds of the loading window are not applied again to these objects, so each sample keeps its own filtering. Then integrate and compare the conditions.

Loading raw files directly in **Multiple Datasets Analysis** is faster, but applies the same QC thresholds to every sample and removes no doublets. Keep it for quick looks or for samples of similar quality. Comparisons between conditions (differential expression within each cluster, pseudobulk test, cluster composition per sample) are only in **Multiple Datasets Analysis**.

📷 Main page with the menu open.

---

## 1. Single dataset analysis (PBMC v3, about 25 min)

### 1.1 Load dataset

1. Under **Dataset Parameters**, set **Species** to **Human** and **Dataset Type** to **snRNA-seq**. This setting is used for single-cell data as well; **Multiome** is only for 10x Multiome files, of which only the gene expression part is analysed.
2. In the **Load Raw 10X Data** card, click **Select file with raw data** and choose `pbmc_1k_v3_filtered_feature_bc_matrix.h5`. Accepted raw formats: Cell Ranger output folder as a ZIP file, `.h5`, `.h5ad`.
3. The log panel reports the number of cells and genes loaded, and the QC plot appears.

Objects analysed before go in the **Load Processed Data** card (**Select processed object (.rds or .h5ad)**): Seurat objects, SingleCellExperiment objects saved as `.rds` (converted to Seurat v5) and AnnData `.h5ad` files.

📷 Load page after loading.

### 1.2 Data cleanup & Variable features

1. Click **Feature Scatter Plots** to see the number of genes per cell and the mitochondrial percentage.
2. Set the thresholds with the sliders (**Unique genes detected in each cell**, **Maximum value for mitochondrial genome**). For blood cells: 200 to 6,000 genes and mitochondrial reads below 20 %. The default values are set for single-nucleus data and would remove most PBMC cells.
3. Click **Apply QC Filters**. The QC plot is updated, and the thresholds applied are recorded for the report.
4. Set the **Number of Variable Features** (2,000 by default) and the **Normalization Method**: LogNormalize (default), SCTransform (regularized negative binomial regression), CLR (centered log-ratio) or RC (relative counts). Click **Normalize Data**.

### 1.3 Dimensional reduction

Click **Scale Data & Run PCA**, then **Show Elbow Plot**. Choose the number of principal components at the bend of the curve (about 15 here).

### 1.4 Clustering

1. Click **Find Neighbors** (it also computes the UMAP).
2. Choose the algorithm (Louvain, multilevel Louvain or SLM) and the resolution, then click **Find Clusters**.
3. **Compute 3D UMAP** adds an interactive 3D view.
4. Below the UMAP, adjust the cluster labels, point size, axes and dark mode, then **Download UMAP**.

📷 UMAP coloured by cluster.

### 1.5 Doublet Detection (recommended on every sample)

Set the expected doublet rate (about 0.8 % per 1,000 cells recovered for 10x), click **Detect Doublets**, check the doublets on the UMAP, then click **Remove Doublets**. The parameters applied (pK and how it was chosen, expected number of doublets, version, seed) are kept with the object and exported in the report.

### 1.6 Plot Gene expressions

- Type genes in each box, separated by commas, semicolons or spaces, for example `CD3E, MS4A1, CD14, NKG7`, and click **Display** (**Display Plot** for the feature plot). Feature, dot, violin and ridge plots each have their own gene box and options (cutoffs, cluster order, z-score, colour palette, 3D display).
- Genes absent from the dataset are named in a warning instead of being silently ignored.
- **Analyze Expression** gives the number and percentage of cells expressing each selected gene in each cluster.

### 1.7 Assigning cell types identity

Three routes, which can be combined. Automated labels are suggestions: they change the clusters only when you apply them, and every change can be undone.

- **Label transfer** from your own annotated Seurat object (any species sharing gene symbols with your data): upload it in **Annotated Seurat object (.rds)**, choose the **Column containing the cell types**, check **Min. prediction score per cell**, **Max. reference cells per label** and **PCA dimensions**, click **Transfer Labels**, then **Apply Labels to Clusters**.
- **SingleR** with built-in references: choose **Species** and **Reference dataset** (mouse: MouseRNAseqData, ImmGenData; human: Human Primary Cell Atlas, Blueprint/ENCODE, Monaco, DICE, Novershtern), the **Label resolution** (broad or fine) and the **Annotation mode** (per cell or per cluster), click **Run SingleR**, then **Apply Labels to Clusters**.
- **Manual**: **Select cluster**, type the **New name** and click **Apply Name**. Giving two clusters the same name merges them. **Undo** reverts the last change.

Confidence: each cell gets a score (prediction score or SingleR delta; cells below the threshold are "Unassigned"), and each cluster a summary (majority label, agreement, second label, fraction unassigned). Clusters below **Min. agreement for a confident cluster (%)** are flagged Low confidence and can keep their current name. **Summary** and **Per cell** download these tables.

A reference only knows its own cell types: a population absent from the reference gets the closest label, often with a high score (with MouseRNAseqData, muscle myonuclei are called "Cardiomyocytes"). Always confirm labels with marker genes.

### 1.8 Clusters comparison

- Choose the **Assay for DE analyses**.
- **Cluster vs all others**: choose the **Target cluster**, set **Min fraction of cells expressing** and **Min |log2FC|**, click **Run analysis**.
- **Pairwise**: choose **Group 1 clusters** and **Group 2 clusters**, click **Run analysis**.
- **Generate Cluster Composition** counts the cells per cluster.

Every result table is stored in the session, so Venn diagrams, volcano plots, heatmaps and enrichment can reuse it. P-values are kept at full precision; a P-value of exactly 0 means below the smallest number the computer can represent (2.2 × 10⁻³⁰⁸). Cell-level P-values treat each cell as a replicate: with few animals or donors, use them to rank genes, and use the pseudobulk test of the Multiple module to compare conditions (section 2.6).

### 1.9 Pathway Enrichment

1. **Genes to analyse**: a comparison of two groups of cells (**Group cells by**, **Group 1**, **Group 2**, optionally restricted to one cluster) or a **Stored table** from another page.
2. **Gene sets**: Gene Ontology for human, mouse or rat (included in the image, no internet needed; choose the **Ontology**), or your own **GMT file** (for example from MSigDB, KEGG or Reactome).
3. Set the size limits of the gene sets and the thresholds, then click **Run Enrichment**.

Over-representation analysis (ORA: up- and down-regulated genes tested separately with a hypergeometric test, against all tested genes) and gene set enrichment analysis (GSEA: all tested genes ranked by log2 fold change, fgsea) are shown in separate tabs, each with its plot, table and parameters.

### 1.10 DEG Visualizations

- **Venn diagram**: choose two or three stored tables, the significance criteria, and click **Generate Venn Diagram**. **View gene list** shows the genes of each part of the diagram; **Download Gene Lists** exports them.
- **Volcano plot**: choose one stored table and the thresholds, click **Generate Volcano Plot**.

### 1.11 Exclusive Biomarkers

This step finds genes that are detected in most cells of a target group and in few other cells. It does not use a new test: it runs a one-sided Wilcoxon test on every gene detected in at least 10 % of target cells, corrects the P-values over all these genes, then applies effect-size filters and ranks the genes by a specificity score S = ((% target − % other) / 100) × log2FC.

1. Choose the **Group by column**, the **Target group(s)** and, optionally, the groups to **Compare against** (all other cells by default).
2. Check the thresholds: **Min % in Target** (50), **Max % in Others** (25), **Min Log2FC** (1.5), **Min Mean Expression** (0.5), **Statistical Test** and **Max P-value**.
3. Click **Find Exclusive Biomarkers**. The table gives, for each gene, the detection rates, log2 fold change, specificity score and adjusted P-value, with the numbers of genes tested and retained.
4. Type genes in **Gene Names** and click **Generate Plot** to draw them.

The default thresholds discard genes detected in few cells, such as many transcription factors, and genes shared with a related population. For such genes, lower **Min % in Target** to 30 %.

### 1.12 Heatmaps

Choose the genes (typed, or the top markers of each cluster), the clusters in order and, optionally, a second axis (**Split by**), then click **Generate Heatmap**. **Export genes CSV** saves the genes shown.

### 1.13 Co-expression Analysis

Counts the cells that express two or more genes together.

1. **Enter genes**, for example `CD3E, CD8A`.
2. **Thresholds on**: raw UMI counts (default, robust to ambient RNA) or normalized expression. One threshold field appears per gene; 1 UMI means any detection.
3. Optional: choose **Background cluster(s)** that should not express these genes and click **Estimate thresholds from background**. Each threshold is set just above the 95th percentile of the gene in these cells (the ambient RNA level).
4. Optional: **Filter cells by column** and **Keep only** restrict the analysis to some cells.
5. Click **Analyze Co-expression**. The results table gives, per cluster, the percentage of cells expressing both genes or only one; **Visualization** shows the same as stacked bars or a heatmap.
6. **DE: Expressors vs Rest**: choose the **Cluster(s) to analyze** and click **Run DE** to compare the co-expressing cells with the others.

### 1.14 Subset

Keep some clusters, or the cells expressing (or not expressing) given genes, with **Create Subset**. The analysis continues on the subset; **Save Subset as RDS** saves it.

### 1.15 Save & Session

- **Save analysis (.rds)** saves the whole analysis (object, stored tables, QC statistics, parameters) in one file, to reload later from the **Load Processed Data** card.
- **Parameters (CSV)** and **PDF report** list every setting applied, the package versions and the seeds.
- The stored tables are listed; select the ones you no longer need and click **Delete selected**.

📷 Save & Session page.

---

## 2. Multiple datasets: integration and comparison (PBMC v3 + v2, about 20 min)

The two demonstration files come from the same donor sequenced with two chemistries: a technical batch effect that integration should remove. For your own data, follow the recommended workflow (section 0.5). To keep this tutorial short, the steps below load the two raw files directly.

### 2.1 Load datasets

1. Click **Load Multiple Datasets**. A window opens; the integration methods are in this window, below the files.
2. In the window, set **Species** = human, the normalization and the QC thresholds (here: 200 to 6,000 genes, mitochondrial reads below 20 %).
3. Set the number of datasets to 2. For each dataset, type a name (`PBMC_v3`, `PBMC_v2`), choose the type (**snRNA-seq** for raw files, **Seurat Object** for samples saved from the Single module) and select the file. The names typed here are used in every plot and table.
4. Each file is processed as soon as it is selected; follow the progress in **Dataset Loading Logs**.
5. When every file is processed, the integration buttons become active. Click one of them:

| Button | Method | When to use it |
|---|---|---|
| **Run Standard Integration** | Seurat CCA anchors | 2–3 datasets; strong correction; memory-intensive |
| **Run RPCA Integration** | Seurat reciprocal PCA | Faster and lighter than CCA; conservative correction, good when cell types differ between samples |
| **Run Harmony Integration** | Harmony on the PCA, for the **Batch variable** you choose | Fast; scales to many datasets |
| **Run fastMNN Integration** | Mutual nearest neighbours (batchelor) | Robust to differences in cell-type composition |
| **Merge Without Integration** | No correction | To see the batch effect, or when samples have none |

An object already integrated (`.rds`) can be loaded directly on the page, without the window.

📷 Loading window with the files loaded and the integration buttons.

### 2.2 Metadata Management

Add the experimental variables of each dataset (condition, genotype, time point, animal) with **Add New Metadata Field**, fill the values, then click **Apply Metadata to Datasets**. Every plot, comparison and co-expression analysis can then be split by these variables. The pseudobulk test needs a column identifying the sample (animal or donor) and one identifying the condition.

### 2.3 Clustering

1. Click **Run Scaling, PCA and Elbow Plot**, then **Find Neighbors & UMAP**, then **Find Clusters**. The corrected embedding of the integration method is used automatically.
2. Colour the UMAP by dataset to check the integration: the two PBMC samples should overlap.
3. **Free Memory** removes intermediate data that are no longer needed, on large datasets.

📷 UMAP coloured by dataset after integration.

### 2.4 Plot Gene expressions

As in the Single module, with three extra settings: **Group By (X-axis)**, **Split By (Optional)** (for example by condition) and **Select Assay**. Feature, dot, violin and ridge plots are drawn with **Run Feature Plot**, **Generate DotPlot**, **Run Vln Plot** and **Run Ridge Plot**.

### 2.5 Assigning cell types identity

Same three routes as in section 1.7.

### 2.6 Clusters comparison, and conditions

- **Cluster vs all others** and **Pairwise**: as in section 1.8.
- **Condition vs condition**: choose the metadata column (for example condition) and its two groups, the clusters (or **Use all clusters**), then the **Statistical test**:
  - **Wilcoxon (cell level)**: each cell is a replicate. Useful to rank genes, but its P-values are inflated when the conditions come from few animals or donors.
  - **Pseudobulk DESeq2 (sample level)**: counts are summed per sample, and samples are the replicates. Choose the sample column in the **Pseudobulk design** panel and the **Minimum cells per sample**. The test is blocked, with an explanation, when a condition has fewer than two samples. Use it to compare conditions.
- **Generate Table** gives the cluster composition per dataset or condition, with **Download Chart**.

### 2.7 Pathway Enrichment, Exclusive Biomarkers, Heatmaps, DEG Visualizations

As in sections 1.9 to 1.12, on the integrated object; enrichment can use a condition comparison within a cluster.

### 2.8 Co-expression Analysis

As in section 1.13, with the comparison between samples or conditions:

1. Enter the genes and thresholds as in section 1.13.
2. **Split by column**: the metadata column compared within each cluster, for example dataset, condition or time point (dataset by default).
3. **Keep only (optional)**: values of that column to analyse; empty = all cells.
4. Click **Analyze Co-expression**. The results table lists every cluster for each value of the split column, with an "Overall" row per value; **Summary** gives the totals per value.
5. In **Visualization**, **Bars per** chooses one bar per value of the split column or one bar per cluster, and the second menu splits the plot into panels. These settings change the plot only, without a new analysis.

### 2.9 Subset and Save & Session

**Apply Cluster Subset**, **Apply Gene Subset** and **Apply Metadata Subset** keep part of the cells. **Save & Session** works as in section 1.15.

---

## 3. CellChat: cell–cell communication (about 15 min)

Analyse PBMC v3 in the Single module up to annotation (section 1) before starting.

### 3.1 Load dataset

1. **CellChat → Load dataset**: choose **Single dataset** and click **Use current object** (or **Or load from disk (.rds)** for a saved object; choose **Multiple datasets** for an integrated object).
2. **Select Species** = Human. Cell-Hub uses CellHub Database 2.1.2: 76,582 ligand–receptor interactions for 21 species, merged from 10 sources (doi:10.5281/zenodo.24283592).
3. **Database Selection**: **All databases (merged)**, or **Select individual databases** to keep only some sources (the tiers are then recomputed on the sources kept).
4. **Confidence tiers to keep**:
   - **High**: supported by 2 or more curated sources.
   - **Medium**: 1 curated source. **Medium: confirmed ligand / receptor roles only** keeps the medium pairs whose ligand is used as a ligand and receptor as a receptor by 2 or more curated sources (human and mouse).
   - **Low**: only orthology-inferred or computationally predicted. These include intracellular or same-cell pairs and should be interpreted with caution.

   The species default is set automatically: high and medium for human and mouse, all tiers for the other species, where nearly all pairs are orthology projections.
5. Click **Load Database**. **Database Information** reports the number of interactions, sources, pathways and interactions per tier. Report the release and tiers used in your methods.

### 3.2 Genes of interest (optional)

To follow a list of genes, for example 20 to 50 differentially expressed genes:

1. In **Genes of interest**, type a **List name**, paste the **Genes** (separated by spaces, commas or new lines) and click **Add list**. Several lists can be saved; a list added with an existing name replaces it.
2. Select a **Saved list** and click **Show interactions in the database**: the interactions of the loaded database whose ligand or receptor, or one of their complex subunits, is in the list are shown with pathway, tier and sources, without running CellChat. Genes absent from the database are reported. The table can be exported.
3. After inference, the bubble and chord plots can be restricted to a list (sections 3.3 and 3.4). The list is applied after inference, so probabilities are computed on the whole database.

### 3.3 Ligand-Receptor

1. **Create Objects**: choose the column that defines the groups (**Column 1**, the cell types), optionally a second and third column (for example condition), optionally **Keep Only Specific Groups** or **Filter by Additional Metadata**, give an **Object name**, then click **Create & Analyze**. One CellChat object is created per condition when a condition column is chosen.
2. **Save Analysis** saves the CellChat object.
3. **Bubble plot**: **Select CellChat object**, **Source cell types**, **Target cell types**, **Threshold**, options (**Flip axes**, **Exclude intra-group**) and, optionally, a **Gene list**. Click **Generate Plot**. **Download Table** exports the significant interactions.

📷 Bubble plot.

### 3.4 Circle Plot

- Choose the object, sender and receiver groups, **Max p-value**, **Top interactions** (by probability) and, optionally, a **Gene list**, then click **Generate Plot** for the chord diagram of ligand–receptor pairs.
- The comparison panel (**Load LR pairs**, then **Generate**) compares two conditions.

---

## 4. Trajectory (Monocle 3, about 10 min)

1. Analyse PBMC v3 in the Single module up to clustering.
2. **Trajectory → Load dataset**: choose **Single dataset** and click **Use current object**.
3. **Trajectory construction**: click **Construct Graph**, choose the root cluster (the starting state) in **Select Root Cluster**, then click **Calculate Pseudotime**.
4. **Differentialy expressed genes**: **Run Graph Test** finds genes that vary along the trajectory; **Run Correlation Test** finds genes correlated with pseudotime.
5. **Gene expression visualisation**: **Plot on Trajectory** shows one gene on the graph; **Generate Plot** shows several genes along pseudotime.
6. **Save Monocle Object (.rds)** saves the result.

Blood cells are not a differentiation series: the trajectory here only shows how the module works. Use your own data, with a known starting population, for a meaningful trajectory.

---

## 5. Spatial transcriptomics (Visium mouse brain, about 20 min)

All analyses use every spot or bin; there is no subsampling.

1. **Load Spatial dataset**: set the **Species** to **Mouse**, select `visium_mouse_brain.zip` in the **Load 10x Visium Data** card as it is (no need to unzip) and click **Load Visium Data**. Visium HD binned outputs load the same way; choose the **Visium HD bin size** first (2, 8 or 16 µm; 8 µm recommended). Spatial objects saved before (`.rds`) go in the **Load Processed Spatial Object** card. The spot or bin size and the orientation of the spots on the image are checked on loading.
2. **Tissue Sections**: **Detect sections** finds the tissue sections of the capture area. Name each section, give it a sample and a condition, then click **Apply names, samples and conditions**. **Remove debris spots** and the spot selection tools remove isolated spots; **Remove ticked sections** removes whole sections.
3. **Normalisation Step**: **Generate QC Plots**, **Apply QC Filters**, then **Normalize Spatial Data**.
4. **Clustering**: **Scale Data & Run PCA**, **Run Neighbors + UMAP**, **Find Clusters**. **Save Spatial Object** saves the result.
5. **Interactive Tissue Viewer**: zoom and pan on the full-resolution image; colour the spots by cluster, gene, QC metric or co-localization; **Download current view**.
6. **Gene expression visualisation**: dot, feature and violin plots on the tissue (**Generate Dot Plot**, **Generate Plots**, **Generate Violin Plot**).
7. **Cluster Annotation**: rename the clusters and click **Apply All Names** (**Reset to Original** undoes it); download the UMAP and the tissue map.
8. **Co-localization**: choose two or more genes (for example `Snap25, Mbp`), a distance in micrometres and the number of permutations, then click **Run Co-localization Analysis**. Results are given per cluster and between clusters, with a permutation test for two genes, a map and a bar chart.
9. **Spatially Variable Genes**: choose the number of variable genes to test and click **Find spatially variable genes** (Moran's I with a one-sided test and Benjamini–Hochberg correction).
10. **Marker Analysis**: choose a cluster and click **Find Markers**.

Not implemented: spot deconvolution and the joint analysis of several capture areas.

📷 Tissue viewer coloured by cluster.

---

## 6. Reproducibility and outputs

- The image fixes every package version (R 4.4.2, Bioconductor 3.20, Seurat 5.5.1, CellChat 2.1.2, Monocle 3 1.3.7 and the others listed in the README). Use the image digest given in the README to be sure to run the exact same image.
- Every random step uses seed 42. CellChat keeps its own default (seed.use = 1).
- **Save & Session** keeps the whole analysis. The parameter table and the PDF report list every setting applied, including those of DoubletFinder, the integration method and the CellChat database (release, sources and tiers).
- Seurat, CellChat and Monocle objects saved by Cell-Hub open in any R session. The analysis functions are in the `functions/` folder of the repository and can be used from R without the interface.
- Downloaded tables keep full numerical precision; P-values are not rounded.

## 7. Troubleshooting

| Problem | Solution |
|---|---|
| The analysis stops or the page disconnects on a large dataset | Raise the memory of Docker Desktop (Settings → Resources); use **Free Memory** after integration |
| "port is already allocated" when starting | Another Cell-Hub is running: stop it, or use `-p 3839:3838` |
| A button does nothing | Read the notification in the bottom-right corner: a setting is missing or a previous step has not been run |
| A gene is missing from a plot | It is absent from the selected assay; the warning names it. Check the spelling and the species |
| Most cells removed by QC | The default thresholds are set for single-nucleus data; adapt them to your data from the QC plot |
| Very small P-values everywhere when comparing conditions | Cell-level tests treat cells as replicates; use the pseudobulk test (section 2.6) |

Bugs and feature requests: https://github.com/GaspardMacaux/Cell-Hub/issues
