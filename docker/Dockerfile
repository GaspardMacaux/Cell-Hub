# syntax=docker/dockerfile:1
###################
# Cell-Hub Docker image
# Build (BuildKit, the GitHub token is passed as a secret and never stored in the image or its history):
#   DOCKER_BUILDKIT=1 docker build --secret id=github_pat,env=GITHUB_PAT -t gaspardmacaux/cell-hub:2.0.0 .
# AMD64 images are built on the Ubuntu workstation, ARM64 images on the Mac.
#
# Base image: satijalab/seurat:5.4.0 — R 4.4.2, Ubuntu 24.04, Bioconductor 3.20
# Same package versions as the tested image:
#   - CRAN packages come from a dated Posit Package Manager snapshot (CRAN_SNAPSHOT = build date of the tested image),
#     as Ubuntu 24.04 binaries (much faster than compiling from source);
#   - Bioconductor is fixed at release 3.20;
#   - every GitHub package is installed at a fixed commit; Seurat 5.5.1 / SeuratObject 5.4.0 are pinned in Layer 4b;
#   - remotes never upgrades an installed package (R_REMOTES_UPGRADE=never).
###################
FROM satijalab/seurat:5.4.0

# Date of the CRAN snapshot (YYYY-MM-DD): build date of the tested image, see `docker inspect -f '{{.Created}}'`
ARG CRAN_SNAPSHOT=2026-09-01
ENV R_REMOTES_UPGRADE=never


###################
# System Dependencies (Layer 1 - never changes)
# gdal-bin and proj-bin provide executables needed by sf at compile time
# libmagick++-dev required by the magick R package (H&E image handling, chord plot assembly)
# CRAN packages are installed as Posit binaries: the system libraries they link to must be present at run time
# (libglpk for igraph, libgsl, libhdf5, libsodium, libgmp/libmpfr, libfftw3, libzmq, libpoppler-cpp, libnode for V8,
#  libprotobuf for otelsdk, libtesseract/libleptonica for tesseract)
# libncurses-dev replaces libncurses5-dev (Ubuntu 24.04)
# python3-pip replaces manual pip bootstrap (Ubuntu 24.04)
###################
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get update && \
    apt-get install -y \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libcairo2-dev \
        libudunits2-dev \
        libgdal-dev \
        libgeos-dev \
        libproj-dev \
        libsqlite3-dev \
        gdal-bin \
        proj-bin \
        gfortran \
        git \
        cmake \
        make \
        pandoc \
        libharfbuzz-dev \
        libfribidi-dev \
        zlib1g-dev \
        libjpeg-dev \
        libpng-dev \
        libtiff-dev \
        libicu-dev \
        libbz2-dev \
        liblzma-dev \
        libncurses-dev \
        libgit2-dev \
        libssh2-1-dev \
        libmagick++-dev \
        libglpk-dev \
        libgsl-dev \
        libhdf5-dev \
        libsodium-dev \
        libgmp-dev \
        libmpfr-dev \
        libfftw3-dev \
        libzmq3-dev \
        libpoppler-cpp-dev \
        libnode-dev \
        libsecret-1-dev \
        libglu1-mesa-dev \
        libprotobuf-dev \
        libtesseract-dev \
        libleptonica-dev \
        python3 \
        python3-pip \
        python3-setuptools \
        wget \
        && \
    pip3 install --no-cache-dir --break-system-packages reportlab==4.4.3 "PyYAML>=6.0.1" && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

###################
# R Base Setup & Package Management (Layer 2 - rarely changes)
# All CRAN installs read the dated snapshot set in Rprofile.site (binary packages need the R user agent);
# BiocManager uses the same CRAN repository for the CRAN dependencies of Bioconductor packages.
# No GitHub token is written anywhere in the image.
###################
RUN case "${CRAN_SNAPSHOT}" in [0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]) ;; *) echo "CRAN_SNAPSHOT must be a date YYYY-MM-DD, got '${CRAN_SNAPSHOT}'"; exit 1 ;; esac && \
    mkdir -p /usr/local/lib/R/etc && \
    echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/noble/${CRAN_SNAPSHOT}'))" >> /usr/local/lib/R/etc/Rprofile.site && \
    echo 'options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])))' >> /usr/local/lib/R/etc/Rprofile.site && \
    echo 'options(BiocManager.check_repositories = FALSE)' >> /usr/local/lib/R/etc/Rprofile.site && \
    R -e 'cat("CRAN repository:", getOption("repos")[["CRAN"]], "\n")' && \
    R -e "install.packages(c('BiocManager', 'remotes', 'R.utils'), dependencies=TRUE)" && \
    R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) stop('BiocManager not installed')" && \
    R -e "if (!requireNamespace('remotes', quietly = TRUE)) stop('remotes not installed')" && \
    R -e "if (!requireNamespace('R.utils', quietly = TRUE)) stop('R.utils not installed')" && \
    R -e 'BiocManager::install(version = "3.20")' && \
    R -e 'BiocManager::install("Rhtslib", dependencies = TRUE, update = FALSE, ask = FALSE)' && \
    R -e 'BiocManager::install("biovizBase", dependencies = TRUE, update = FALSE, ask = FALSE)'

###################
# R CRAN Packages (Layer 3 - rarely changes)
# Removed: rJava, RBioFormats (unused)
# Removed: stringr, dplyr, readr, ggplot2, purrr, tibble, tidyr (covered by tidyverse)
# Removed: grr (archived, installed later from archive)
# sf NOT included here - installed in Layer 9b alongside spdep
###################
ARG BUST_CRAN_CACHE=1
RUN R -e 'cat("CRAN cache bust:", Sys.getenv("BUST_CRAN_CACHE", "1"), "\n")' && \
    R -e 'install.packages(c(\
    "plotly", "colourpicker", "patchwork", "ggrepel", \
    "DT", "cowplot", "data.table", "httr", "hdf5r", \
    "clustree", "viridis", "callr", \
    "shiny", "shinyWidgets", "shinydashboard", "shinyjs", "shinyFiles", "fontawesome", "jsonlite", \
    "pbmcapply", "pscl", "rsample", "speedglm", \
    "fastDummies", "fitdistrplus", "ggridges", "ica", "leiden", "pbapply", "RANN", \
    "RcppAnnoy", "ROCR", "Rtsne", "scattermore", "sctransform", "spatstat.explore", \
    "spatstat.geom", "uwot", "RcppProgress", \
    "future", "future.apply", "progressr", "Rcpp", "spam", "RcppEigen", \
    "tidyverse", "proxyC", "igraph", "tidygraph", "ggraph", \
    "reshape2", "Matrix", "tidytext", "rlang", "progress", "circlize", "doParallel", \
    "foreach", "caret", "glmnet", "randomForest", "viridisLite", "zoo", \
    "cluster", "generics", "graphics", "grDevices", "grid", "irlba", "KernSmooth", \
    "lifecycle", "MASS", "matrixStats", "miniUI", "scales", "stats", \
    "fields", "parallel", \
    "NMF", "ggalluvial", "svglite", "leidenbase", "FNN", "RSpectra", "shape", \
    "VennDiagram", "openxlsx", "arrow", "glmGamPoi", \
    "readODS", "tiff", "png", "jpeg", "magick", \
    "lmtest", "pheatmap", "plyr", "proxy", "RhpcBLASctl", "RcppHNSW", "slam", \
    "ggnetwork", \
    "Polychrome"\
    ), dependencies=TRUE)'

###################
# Deriv (Layer 3b)
# Required by doBy (pbkrtest -> car -> ggpubr chain) but missing from the 2026-09-01 snapshot: installed at the version of the tested image
###################
RUN R -e 'if (!requireNamespace("Deriv", quietly = TRUE)) remotes::install_version("Deriv", version = "4.1.6", repos = "https://cloud.r-project.org", upgrade = "never"); cat("Deriv", as.character(packageVersion("Deriv")), "\n")'

###################
# R Bioconductor Packages (Layer 4 - rarely changes)
###################
RUN R -e 'BiocManager::install(c(\
    "BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma", "lme4", \
    "S4Vectors", "SingleCellExperiment", "SummarizedExperiment", "batchelor", \
    "HDF5Array", "terra", "ggrastr", "Biobase", "MatrixGenerics", "matrixStats", \
    "GenomeInfoDb", "IRanges", "stats4", "GenomicRanges", "circlize", \
    "RColorBrewer", "magrittr", "forcats", "tidyselect", "EnsDb.Hsapiens.v86", \
    "biomaRt", "edgeR", "DESeq2", "AUCell", "RcisTarget", "GSEABase", "qvalue", \
    "ComplexHeatmap", "ROCR", "mixtools", "viridisLite", \
    "rhdf5", \
    "muscat", "scater", "locfdr", "sva", "scran", "UpSetR", "factoextra", \
    "harmony"\
    ), update=FALSE, ask=FALSE)' && \
    R -e 'if (!requireNamespace("biovizBase", quietly = TRUE)) { \
        BiocManager::install("biovizBase", force = TRUE); \
    }'

###################
# Load check of the CRAN / Bioconductor binaries (Layer 4a)
# Every package installed so far is loaded in its own R process; a missing system library is reported here,
# with the package name, instead of breaking a later step. The list is printed in the build log; the build stops if any fails.
# Excluded: rJava and mallet (need a Java runtime) and av (needs the FFmpeg libraries): optional dependencies pulled in by
# dependencies=TRUE, never used by Cell-Hub.
###################
RUN Rscript -e 'cat(rownames(installed.packages(lib.loc = "/usr/local/lib/R/site-library")), sep = "\n")' > /tmp/pkgs.txt && \
    grep -vxF -e rJava -e mallet -e av /tmp/pkgs.txt | xargs -P 8 -I{} sh -c 'Rscript -e "suppressMessages(requireNamespace(\"{}\", quietly = TRUE)) || quit(status = 1)" > /dev/null 2>&1 || echo "{}"' > /tmp/failed_pkgs.txt; \
    echo "Packages checked: $(wc -l < /tmp/pkgs.txt)"; \
    if [ -s /tmp/failed_pkgs.txt ]; then echo "Packages that fail to load:"; cat /tmp/failed_pkgs.txt; \
      for p in $(cat /tmp/failed_pkgs.txt); do echo "== $p"; Rscript -e "requireNamespace('$p')" 2>&1 | grep -i "cannot open shared object\|error" | head -3; done; exit 1; \
    else echo "All packages load"; fi

###################
# Seurat pinned (Layer 4b)
# The app is developed and tested with Seurat 5.5.1 / SeuratObject 5.4.0 (newer than the base image):
# both are installed at fixed commits, and the build fails if another version ends up installed
###################
RUN --mount=type=secret,id=github_pat \
    export GITHUB_PAT="$(cat /run/secrets/github_pat 2>/dev/null)" && \
    R -e 'remotes::install_github("satijalab/seurat-object@61b2d5c53ce7500c7df033b129f85f95749b8b42", dependencies=FALSE, upgrade="never")' && \
    R -e 'remotes::install_github("satijalab/seurat@4c0f2dc16fad4e8f0d7b5c98321d8bcb18caa13a", dependencies=FALSE, upgrade="never")' && \
    R -e 'v <- c(Seurat = "5.5.1", SeuratObject = "5.4.0"); for (p in names(v)) { got <- as.character(packageVersion(p)); cat(p, got, "\n"); if (got != v[[p]]) stop(p, " ", got, " installed, ", v[[p]], " expected") }'

###################
# Install ComplexHeatmap from GitHub (Layer 5)
# Must be installed before CellChat (version compatible with CellChat)
# The GitHub token is read from the BuildKit secret, only for this RUN
###################
RUN --mount=type=secret,id=github_pat \
    GITHUB_PAT="$(cat /run/secrets/github_pat 2>/dev/null)" \
    R -e 'remotes::install_github("jokergoo/ComplexHeatmap@ad11b26afe8c63dc05b20035503352ca740f602e", dependencies=TRUE, upgrade="never")'

###################
# Install CellChat and presto (Layer 6)
###################
RUN --mount=type=secret,id=github_pat \
    export GITHUB_PAT="$(cat /run/secrets/github_pat 2>/dev/null)" && \
    R -e 'remotes::install_github("jinworks/CellChat@v2.1.2", dependencies=FALSE, upgrade="never")' && \
    R -e 'remotes::install_github("immunogenomics/presto@3c97180d90d330524b5b8353a342ab2f4b990d21", dependencies=TRUE, upgrade="never")' && \
    R -e 'library(CellChat); cat("CellChat version:", as.character(packageVersion("CellChat")), "\n")'

###################
# Install schard for H5AD (AnnData) reading (Layer 6b)
# Pure R, depends on hdf5r (Layer 3) and Seurat (base image)
###################
RUN --mount=type=secret,id=github_pat \
    export GITHUB_PAT="$(cat /run/secrets/github_pat 2>/dev/null)" && \
    R -e 'remotes::install_github("cellgeni/schard@f9c30c25a192175a1e6523e32f0e2047de1e8685", dependencies=TRUE, upgrade="never")' && \
    R -e 'library(schard); cat("schard installed\n")'

###################
# Install Signac Dependencies from Bioconductor (Layer 7)
# Signac is only needed to read ChromatinAssay objects (multiome RDS); only the RNA assay is analysed
###################
RUN R -e 'BiocManager::install(c(\
    "chromVAR", \
    "TFBSTools", \
    "motifmatchr", \
    "rtracklayer", \
    "Rsamtools", \
    "GenomicAlignments", \
    "BSgenome", \
    "RcppRoll" \
    ), update=FALSE, ask=FALSE)'

###################
# Copy Local Package Sources (Layer 8)
# seurat and seurat-object provided by base image — not copied here
# Change BUST_SEURAT_CACHE to rebuild local packages without rebuilding layers above
###################
ARG BUST_SEURAT_CACHE=1
RUN echo "Cache bust: ${BUST_SEURAT_CACHE}"
COPY seurat-wrappers /tmp/seurat-wrappers
COPY seurat-disk /tmp/seurat-disk
COPY monocle3 /tmp/monocle3

###################
# Install Local Packages except monocle3 (Layer 9a)
# Each in a separate RUN to isolate failures
###################
RUN R -e 'install.packages("/tmp/seurat-wrappers", repos=NULL, type="source")'

RUN R -e 'install.packages("/tmp/seurat-disk", repos=NULL, type="source")'

# DoubletFinder pinned to commit 1b244d8; the build fails if another commit is installed
RUN --mount=type=secret,id=github_pat \
    export GITHUB_PAT="$(cat /run/secrets/github_pat 2>/dev/null)" && \
    R -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder@1b244d8f0d54b4b1cb4365639931bbb16f01e1cd", dependencies=FALSE, upgrade="never")' && \
    R -e 'library(DoubletFinder); d <- packageDescription("DoubletFinder"); cat("DoubletFinder version:", d$Version, "| commit:", d$RemoteSha, "\n"); if (!identical(d$RemoteSha, "1b244d8f0d54b4b1cb4365639931bbb16f01e1cd")) stop("DoubletFinder commit mismatch")'

###################
# Install monocle3 separately (Layer 9b)
# grr must be installed from archive before monocle3 (archived on CRAN)
# sf, s2, units installed here as monocle3/spdep dependencies
###################
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/grr/grr_0.9.5.tar.gz", repos=NULL, type="source")'

RUN R -e 'install.packages(c("units", "s2", "sf", "spdep"), dependencies=TRUE)'

RUN R CMD INSTALL /tmp/monocle3 && \
    R -e 'library(monocle3); cat("monocle3 version:", as.character(packageVersion("monocle3")), "\n")'

###################
# Install Signac version 1.14.0 (Layer 10)
###################
RUN R -e 'remotes::install_version("Signac", version = "1.14.0", repos = getOption("repos"), upgrade="never")' && \
    R -e 'library(Signac); cat("Signac version:", as.character(packageVersion("Signac")), "\n")'

###################
# Spatial helpers (Layer 10b)
# nanoparquet: reads Visium HD tissue_positions.parquet (prepare_visium_directory in spatial_functions.R)
###################
RUN R -e 'install.packages("nanoparquet"); if (!requireNamespace("nanoparquet", quietly = TRUE)) stop("nanoparquet installation failed - see the compiler output above"); cat("nanoparquet version:", as.character(packageVersion("nanoparquet")), "\n")'

###################
# Automated cell-type annotation (Layer 10c)
# SingleR + celldex from Bioconductor 3.20 (matches the base image)
# The 7 celldex references used by the app are downloaded once at build time and cached as RDS,
# so annotation works offline; CELLHUB_REF_CACHE is read by load_singler_reference()
###################
ENV CELLHUB_REF_CACHE=/opt/cellhub/singler_references
RUN R -e 'BiocManager::install(c("SingleR", "celldex"), update = FALSE, ask = FALSE)' && \
    R -e 'for (p in c("SingleR", "celldex")) if (!requireNamespace(p, quietly = TRUE)) stop(p, " installation failed"); \
          cat("SingleR", as.character(packageVersion("SingleR")), "| celldex", as.character(packageVersion("celldex")), "\n")' && \
    mkdir -p ${CELLHUB_REF_CACHE} && \
    R -e 'refs <- c("HumanPrimaryCellAtlasData", "BlueprintEncodeData", "MonacoImmuneData", "DatabaseImmuneCellExpressionData", \
                    "NovershternHematopoieticData", "MouseRNAseqData", "ImmGenData"); \
          cache <- Sys.getenv("CELLHUB_REF_CACHE"); \
          for (r in refs) { ref <- getExportedValue("celldex", r)(); saveRDS(ref, file.path(cache, paste0(r, ".rds"))); cat("cached", r, ncol(ref), "samples\n") }' && \
    rm -rf /root/.cache/R/ExperimentHub /root/.cache/R/gypsum /root/.cache/gypsum && \
    chmod -R a+rX /opt/cellhub && \
    ls -lh ${CELLHUB_REF_CACHE}

###################
# Version check (Layer 10d)
# Prints the versions reported in the manuscript; the build log is the reference for software_versions.md
###################
RUN R -e 'for (p in c("Seurat", "SeuratObject", "CellChat", "monocle3", "harmony", "batchelor", "DoubletFinder", "SingleR", "celldex", "shiny")) \
            cat(sprintf("%-14s %s\n", p, as.character(packageVersion(p)))); cat("R", R.version$major, ".", R.version$minor, "\n", sep = "")'

###################
# Pathway enrichment (Layer 10e)
# clusterProfiler + OrgDb (mouse, human, rat) from Bioconductor 3.20
# Install warnings ("non-zero exit status") are turned into a build failure and every namespace is loaded explicitly,
# so the real cause (missing dependency, version conflict) is printed instead of a generic "installation failed"
###################
RUN R -e 'options(Ncpus = max(1L, parallel::detectCores() - 1L), timeout = 600); failed <- character(0); \
          withCallingHandlers(BiocManager::install(c("GO.db", "fgsea", "org.Mm.eg.db", "org.Hs.eg.db", "org.Rn.eg.db"), update = FALSE, ask = FALSE), \
            warning = function(w) if (grepl("non-zero exit status|installation of .* failed|not available|cannot open URL", conditionMessage(w))) failed <<- c(failed, conditionMessage(w))); \
          if (length(failed) > 0) { cat("FAILED INSTALLS:\n", paste(failed, collapse = "\n"), "\n"); quit(status = 1) }; \
          for (p in c("AnnotationDbi", "GO.db", "fgsea", "org.Mm.eg.db", "org.Hs.eg.db", "org.Rn.eg.db")) \
            tryCatch(loadNamespace(p), error = function(e) { cat("LOAD ERROR", p, ":", conditionMessage(e), "\n"); quit(status = 1) }); \
          md <- AnnotationDbi::metadata(org.Mm.eg.db::org.Mm.eg.db); \
          cat("fgsea", as.character(packageVersion("fgsea")), "| GO.db", as.character(packageVersion("GO.db")), "| org.Mm.eg.db", as.character(packageVersion("org.Mm.eg.db")), \
              "| org.Hs.eg.db", as.character(packageVersion("org.Hs.eg.db")), "| org.Rn.eg.db", as.character(packageVersion("org.Rn.eg.db")), \
              "| GO source date", md$value[md$name == "GOSOURCEDATE"], "| ggplot2", as.character(packageVersion("ggplot2")), "\n")' && \
    R -e 'ann <- suppressMessages(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, keys = c("Myh7", "Tnnt1", "Ckm"), keytype = "SYMBOL", columns = c("GOALL", "ONTOLOGYALL"))); \
          ann <- ann[!is.na(ann$GOALL) & ann$ONTOLOGYALL == "BP", ]; if (nrow(ann) == 0) stop("GO smoke test: no BP annotation for Myh7/Tnnt1/Ckm"); \
          term <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = "GO:0006936", keytype = "GOID", columns = "TERM"))$TERM; \
          set.seed(42); s <- setNames(rnorm(1000), paste0("g", 1:1000)); r <- fgsea::fgsea(list(a = names(s)[1:50], b = names(s)[51:120]), s, minSize = 10, eps = 0, nproc = 1); \
          cat("GO smoke test:", length(unique(ann$GOALL)), "BP terms | GO:0006936 =", term, "| fgsea:", nrow(r), "sets\n")'
###################
# Copy PDF report generation script (Layer 11)
###################
COPY generate_report.py /app/generate_report.py

###################
# Copy VOLATILE Application Files (Layer 12)
# CellHub.tar.gz and app files change frequently — kept at the end
###################
COPY Rprofile.site /Rprofile.site
COPY CellHub.tar.gz /tmp/CellHub.tar.gz

###################
# Install Application Package (Layer 13)
###################
RUN R -e 'install.packages("/tmp/CellHub.tar.gz", repos = NULL, type = "source")'

###################
# Copy Application (Layer 14)
# The app runs from this copy; LICENSE is shipped inside mon_app
###################
WORKDIR /srv/shiny-server/
RUN mkdir -p /srv/shiny-server/mon_app
COPY ./CellHub/inst/mon_app/ /srv/shiny-server/mon_app/

# Verify content
RUN ls -la /srv/shiny-server/mon_app/ && test -f /srv/shiny-server/mon_app/LICENSE

# Startup script
RUN echo 'shiny::runApp("/srv/shiny-server/mon_app", host="0.0.0.0", port=3838)' > /start_app.R

# Image metadata (kept at the end: changing it never invalidates the build cache)
LABEL org.opencontainers.image.title="Cell-Hub" \
      org.opencontainers.image.description="Graphical interface for single-cell, single-nucleus and spatial RNA-seq analysis" \
      org.opencontainers.image.version="2.0.0" \
      org.opencontainers.image.licenses="GPL-3.0-only" \
      org.opencontainers.image.source="https://github.com/GaspardMacaux/Cell-Hub" \
      org.opencontainers.image.authors="Gaspard Macaux <gaspard.macaux@gmail.com>"

EXPOSE 3838

ENTRYPOINT ["Rscript", "/start_app.R"]
