###############################################################################
# edgeR Differential Expression Analysis
# Robust to ALL featureCounts TSV formats
# NO column position assumptions
###############################################################################

############################
# STEP 0: Load libraries
############################
suppressMessages({
  library(edgeR)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

############################
# STEP 1: Read featureCounts TSV
############################
counts_raw <- fread(
  "counts/gene_counts.tsv",
  sep = "\t",
  header = TRUE,
  data.table = FALSE
)

cat("Detected columns:\n")
print(colnames(counts_raw))

############################
# STEP 2: Automatically detect count columns
############################
# Rule:
# - Count columns are numeric
# - Annotation columns are character / factor

is_numeric <- sapply(counts_raw, is.numeric)

if (sum(is_numeric) == 0) {
  stop("No numeric columns detected. This is not a valid featureCounts file.")
}

count_matrix <- as.matrix(counts_raw[, is_numeric])
rownames(count_matrix) <- counts_raw[[1]]   # GeneID always first column

############################
# STEP 3: Sanity checks
############################
stopifnot(is.matrix(count_matrix))
stopifnot(nrow(count_matrix) > 0)
stopifnot(ncol(count_matrix) > 0)

cat("Genes:", nrow(count_matrix), "\n")
cat("Samples:", ncol(count_matrix), "\n")

############################
# STEP 4: Sample metadata
############################
# IMPORTANT: must match number of samples
sample_info <- data.frame(
  sample = colnames(count_matrix),
  condition = factor(c(
    "Control", "Control",
    "Treated", "Treated"
  ))
)

if (nrow(sample_info) != ncol(count_matrix)) {
  stop("Sample metadata does not match number of samples.")
}

############################
# STEP 5: Create DGEList
############################
dge <- DGEList(
  counts = count_matrix,
  group = sample_info$condition
)

############################
# STEP 6: Filter low-expression genes
############################
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

cat("Genes after filtering:", nrow(dge), "\n")

############################
# STEP 7: Normalize (TMM)
############################
dge <- calcNormFactors(dge)

############################
# STEP 8: Design matrix
############################
design <- model.matrix(~ condition, data = sample_info)

############################
# STEP 9: Estimate dispersion
############################
dge <- estimateDisp(dge, design)
plotBCV(dge)

############################
# STEP 10: Fit GLM (QL)
############################
fit <- glmQLFit(dge, design)

############################
# STEP 11: Differential testing
############################
qlf <- glmQLFTest(fit, coef = 2)

deg_table <- topTags(qlf, n = Inf)$table
deg_table$FDR <- p.adjust(deg_table$PValue, method = "BH")

############################
# STEP 12: Save all results
############################
dir.create("results", showWarnings = FALSE)

write.table(
  deg_table,
  file = "results/edgeR_all_DEGs.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

############################
# STEP 13: Significant DEGs
############################
sig_deg <- deg_table[
  abs(deg_table$logFC) >= 1 &
    deg_table$FDR < 0.05,
]

write.table(
  sig_deg,
  file = "results/edgeR_significant_DEGs.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

cat("Significant DEGs:", nrow(sig_deg), "\n")

############################
# STEP 14: Volcano plot
############################
volcano_df <- data.frame(
  logFC = deg_table$logFC,
  negLogFDR = -log10(deg_table$FDR)
)

ggplot(volcano_df, aes(logFC, negLogFDR)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot (edgeR)",
    x = "log2 Fold Change",
    y = "-log10(FDR)"
  )

############################
# STEP 15: Heatmap
############################
if (nrow(sig_deg) >= 2) {
  top_genes <- rownames(sig_deg)[1:min(50, nrow(sig_deg))]
  logCPM <- cpm(dge, log = TRUE)
  
  pheatmap(
    logCPM[top_genes, ],
    scale = "row",
    annotation_col = sample_info["condition"],
    main = "Top Differentially Expressed Genes"
  )
}

############################
# STEP 16: MDS plot
############################
plotMDS(
  dge,
  labels = sample_info$condition,
  col = as.numeric(sample_info$condition)
)

###############################################################################
# END OF SCRIPT
###############################################################################
