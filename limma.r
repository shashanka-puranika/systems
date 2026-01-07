




library(limma)

data("swirl")
group <- factor(c(rep("A",3), rep("B",3)))
design <- model.matrix(~group)

fit <- lmFit(swirl, design)
fit <- eBayes(fit)

topTable(fit)

library(maftools)
