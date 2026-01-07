###############################################################################
# Somatic Variant Analysis using maftools
#
# Input  : MAF file (VCF2MAF / TCGA / Mutect2 compatible)
# Output : Mutation summaries, plots, statistics
#
# Author : Standard Academic Genomics Pipeline
###############################################################################

############################
# STEP 0: Load libraries
############################
suppressMessages({
  library(maftools)
  library(data.table)
  library(ggplot2)
})

############################
# STEP 1: Read MAF file safely
############################
maf_file <- "maf/sample.maf"   # MODIFY PATH IF NEEDED

if (!file.exists(maf_file)) {
  stop("MAF file not found at: ", maf_file)
}

maf <- read.maf(
  maf = maf_file,
  verbose = FALSE
)

############################
# STEP 2: Basic sanity checks
############################
cat("Number of samples :", length(getSampleSummary(maf)$Tumor_Sample_Barcode), "\n")
cat("Number of genes   :", length(getGeneSummary(maf)$Hugo_Symbol), "\n")

############################
# STEP 3: MAF summary statistics
############################
pdf("results/maf_summary.pdf", width = 8, height = 6)
plotmafSummary(
  maf = maf,
  rmOutlier = TRUE,
  addStat = "median",
  dashboard = TRUE
)
dev.off()

############################
# STEP 4: Mutation classification overview
############################
pdf("results/mutation_classification.pdf", width = 8, height = 6)
plotMutationTypes(maf)
dev.off()

############################
# STEP 5: Oncoplot (top mutated genes)
############################
pdf("results/oncoplot_top20.pdf", width = 10, height = 8)
oncoplot(
  maf = maf,
  top = 20,
  removeNonMutated = TRUE,
  showTumorSampleBarcodes = TRUE
)
dev.off()

############################
# STEP 6: Transition / Transversion analysis
############################
pdf("results/tv_ti_plot.pdf", width = 7, height = 6)
titv <- titv(maf = maf, plot = TRUE)
dev.off()

############################
# STEP 7: Variant classification distribution
############################
pdf("results/variant_classification_barplot.pdf", width = 8, height = 6)
plotVariantClassification(maf)
dev.off()

############################
# STEP 8: Variant allele frequency analysis
############################
if ("t_vaf" %in% colnames(maf@data)) {
  pdf("results/vaf_distribution.pdf", width = 8, height = 6)
  vafPlot(maf, vafCol = "t_vaf")
  dev.off()
}

############################
# STEP 9: Gene-level mutation burden
############################
gene_summary <- getGeneSummary(maf)

write.table(
  gene_summary,
  file = "results/gene_mutation_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################
# STEP 10: Sample-level mutation burden
############################
sample_summary <- getSampleSummary(maf)

write.table(
  sample_summary,
  file = "results/sample_mutation_summary.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################
# STEP 11: Lollipop plot for a key gene
############################
# MODIFY gene symbol as required (e.g., TP53, KRAS)
target_gene <- "TP53"

if (target_gene %in% maf@data$Hugo_Symbol) {
  pdf(paste0("results/lollipop_", target_gene, ".pdf"), width = 10, height = 6)
  lollipopPlot(
    maf = maf,
    gene = target_gene,
    AACol = "Protein_position"
  )
  dev.off()
}

############################
# STEP 12: Co-occurrence and mutual exclusivity
############################
pdf("results/cooccurrence_plot.pdf", width = 10, height = 8)
somaticInteractions(
  maf = maf,
  top = 25,
  pvalue = c(0.05)
)
dev.off()

############################
# STEP 13: Subsetting MAF (example)
############################
# Example: keep only Missense mutations
missense_maf <- subsetMaf(
  maf = maf,
  query = "Variant_Classification == 'Missense_Mutation'"
)

write.table(
  missense_maf@data,
  file = "results/missense_only.maf",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################
# STEP 14: Save R objects
############################
saveRDS(maf, file = "results/maf_object.rds")
saveRDS(missense_maf, file = "results/maf_missense_only.rds")

###############################################################################
# END OF maftools PIPELINE
###############################################################################
