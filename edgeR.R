###############################################################################
# edgeR Differential Expression Analysis (ROBUST FEATURECOUNTS TSV INPUT)
#
# Fixes:
# - Does NOT assume counts start at column 7
# - Automatically detects count columns
# - Prevents "undefined columns selected" error
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
# STEP 1: Read featureCounts TSV safely
############################
counts_raw <- fread(
  "counts/gene_counts.tsv",
  sep = "\t",
  header = TRUE,
  data.table = FALSE
)

# Inspect structure (safe for debugging)
cat("Columns in featureCounts file:\n")
print(colnames(counts_raw))

############################
# STEP 2: Identify count columns automatically
############################
# featureCounts annotation columns usually include:
# Geneid, Chr, Start, End, Strand, Length
# Count columns = everything AFTER these

annotation_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

missing_cols <- setdiff(annotation_cols, colnames(counts_raw))
if (length(missing_cols) > 0) {
  stop("Missing expected annotation columns: ",
       paste(missing_cols, collapse = ", "))
}

count_col_start <- which(colnames(counts_raw) == "Length") + 1
count_cols <- count_col_start:ncol(counts_raw)

############################
# STEP 3: Build count matrix
############################
count_matrix <- as.matrix(counts_raw[, count_cols])
rownames(count_matrix) <- counts_raw$Geneid

# Sanity checks
stopifnot(is.matrix(count_matrix))
stopifnot(nrow(count_matrix) > 0)
stopifnot(ncol(count_matrix) > 0)

cat("Genes:", nrow(count_matrix), "\n")
cat("Samples:", ncol(count_matrix), "\n")

############################
# STEP 4: Define sample metadata
############################
# MODIFY THIS to match your experiment
sample_info <- data.frame(
  sample = colnames(count_matrix),
  condition = factor(c(
    "Control", "Control",
    "Treated", "Treated"
  ))
)

stopifnot(nrow(sample_info) == ncol(count_matrix))

############################
# STEP 5: Create DGEList object
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
# STEP 10: Fit GLM (QL framework)
############################
fit <- glmQLFit(dge, design)

############################
# STEP 11: Differential testing
############################
qlf <- glmQLFTest(fit, coef = 2)

deg_table <- topTags(qlf, n = Inf)$table
deg_table$FDR <- p.adjust(deg_table$PValue, method = "BH")

############################
# STEP 12: Save results
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
# STEP 15: Heatmap of top DEGs
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
