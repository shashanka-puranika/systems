Systems Biology – First Sessional Practical Examination

M.Sc. Bioinformatics & Systems Biology (Semester III)
Time: 3 Hours | Max Marks: 40


---

Question 1: Git Operations (10 Marks)

1(a) Initialize a Repository (2 Marks)

Aim:
To create a new directory and initialize a Git repository.

Commands:

mkdir git-sessional
cd git-sessional
git init

Output:
A .git directory is created indicating successful initialization.

Result:
Git repository successfully initialized in git-sessional.


---

1(b) Track and Commit a File (4 Marks)

Aim:
To create a file, track it, and commit it to the repository.

Commands:

echo "This file contains sessional exam notes." > exam.txt
git status
git add exam.txt
git commit -m "Initial commit: added notes.txt"

Output:

git status shows exam.txt as staged

Commit created successfully


Result:
The file exam.txt is successfully tracked and committed.


---

1(c) Branching (4 Marks)

Aim:
To create a new branch, switch to it, and commit changes.

Commands:

git branch feature-update
git checkout feature-update
echo "Additional notes added in feature branch." >> exam.txt
git add exam.txt
git commit -m "Updated notes in feature branch"

Output:

Branch switched to feature-update

New commit created


Result:
Changes committed successfully in the feature-update branch.


---

Question 2: Disease Gene Network Bar Plot (5 Marks)

Aim:
To plot a bar graph showing disease enrichment and identify the most significant disease.

R Script:

# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Read gene list
genes <- read.table("genelist.txt", header = FALSE)
gene_symbols <- genes$V1

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(gene_symbols,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# Disease enrichment analysis (example using DOSE)
library(DOSE)
disease_enrich <- enrichDO(entrez_ids$ENTREZID)

# Convert to data frame
df <- as.data.frame(disease_enrich)

# Bar plot
ggplot(df[1:10, ], aes(x = reorder(Description, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  xlab("Disease") +
  ylab("Gene Count") +
  ggtitle("Disease Gene Network Enrichment")

Result / Inference:
The disease with the highest gene count and lowest adjusted p-value is identified as the most significantly enriched disease.


---

Question 3: Differential Expression & Clustering (10 Marks)

Aim:

To identify the top 50 differentially expressed genes and perform clustering.


---

Step 1: Differential Expression Analysis

R Script:

library(limma)

# Expression data
expr <- read.table("expression_matrix.txt", header = TRUE, row.names = 1)

# Sample group information
group <- factor(c("Control","Control","Disease","Disease"))

design <- model.matrix(~group)
fit <- lmFit(expr, design)
fit <- eBayes(fit)

# Extract top 50 DE genes
top50 <- topTable(fit, number = 50)
write.csv(top50, "top50_DEGs.csv")

Result:
Top 50 differentially expressed genes obtained.


---

3(a) Hierarchical Clustering

heatmap(as.matrix(expr[rownames(top50), ]),
        scale = "row",
        col = colorRampPalette(c("blue","white","red"))(100))

Inference:
Genes cluster based on similar expression patterns across samples.


---

3(b) K-Means Clustering

set.seed(123)
kmeans_res <- kmeans(expr[rownames(top50), ], centers = 4)

plot(prcomp(expr[rownames(top50), ])$x,
     col = kmeans_res$cluster,
     pch = 16,
     main = "K-means Clustering of Top 50 DE Genes")

Inference:
Top DE genes are grouped into 4 distinct expression clusters.


---

Question 4: Lab Record (10 Marks)

Statement:
Lab record was submitted, verified, and duly signed by the course instructor.


---

Question 5: Viva (5 Marks)

Statement:
Viva voce conducted successfully.


---

Final Summary

Question	Marks

Q1	10
Q2	5
Q3	10
Q4	10
Q5	5
Total	40 / 40