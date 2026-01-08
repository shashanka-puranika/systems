# Install BiocManager if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# ---- CRAN packages ----
install.packages(c(
  "tidyverse",
  "data.table",
  "ggplot2",
  "pheatmap",
  "EnhancedVolcano",
  "RColorBrewer",
  "factoextra",
  "cluster",    # for K-means clustering
  "gridExtra"   # for arranging plots
))

# ---- Bioconductor packages ----
BiocManager::install(c(
  "limma",
  "Biobase",
  "edgeR",
  "genefilter",
  "marray",
  "GOstats",
  "clusterProfiler",
  "org.Hs.eg.db",
  "GOSemSim",
  "GO.db"
), ask = FALSE, update = FALSE, force = TRUE)


# ---- Load packages to check if installed ----
cran_pkgs <- c("tidyverse","data.table","ggplot2","pheatmap",
               "EnhancedVolcano","RColorBrewer","factoextra","cluster","gridExtra")
bioc_pkgs <- c("limma","Biobase","edgeR","genefilter","marray",
               "GOstats","clusterProfiler","org.Hs.eg.db","GOSemSim","GO.db")

for (pkg in c(cran_pkgs, bioc_pkgs)) {
  if (suppressWarnings(require(pkg, character.only = TRUE))) {
    message("✅ Loaded: ", pkg)
  } else {
    message("❌ Failed: ", pkg)
  }
}                                                                                                                                                                                                                                                if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("maftools")


check this for RNAseq & DNAvariant calling >> #!/bin/bash

echo "========================================"
echo " RNA-seq + DNA Variant Calling Setup"
echo "========================================"

# Step 1: Configure conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Step 2: Create conda environment
conda create -y -n ngs_practical_env \
    python=3.10 \
    wget \
    fastqc \
    trimmomatic \
    star \
    subread \
    samtools \
    bwa \
    bcftools \
    vcftools \
    freebayes \
    r-base=4.3 \
    r-essentials \
    bioconductor-maftools

# Step 3: Activate environment
echo "Activate environment using:"
echo "conda activate ngs_practical_env"

echo "========================================"
echo " Installation completed successfully"
echo "========================================"