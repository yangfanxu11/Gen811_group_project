#!/usr/bin/env Rscript
# EdgeR differential junction usage analysis
library(edgeR)

cat("Loading junction counts...\n")
counts <- read.table(gzfile("../data/leafcutter/leafcutter_perind.counts.gz"), 
                     header=TRUE, sep=" ", check.names=FALSE)

# Extract junction IDs
intron_ids <- counts[, 1]

# Extract numerator from "x/y" format
extract_numerator <- function(x) {
  as.numeric(sapply(strsplit(as.character(x), "/"), `[`, 1))
}

counts_mat <- sapply(counts[, -1], extract_numerator)
rownames(counts_mat) <- intron_ids

cat("Total junctions:", nrow(counts_mat), "\n")

# Filter low expression
keep <- rowSums(counts_mat >= 10) >= 2
counts_mat <- counts_mat[keep, ]
cat("Junctions after filtering:", nrow(counts_mat), "\n")

# edgeR analysis
groups <- factor(c(rep("WT",4), rep("KD",2)))
y <- DGEList(counts=counts_mat, group=groups)
y <- calcNormFactors(y)
design <- model.matrix(~groups)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)

# Save results
results <- topTags(qlf, n=Inf)
write.table(results$table, "../data/leafcutter/edgeR_diff_junctions.txt", 
            sep="\t", quote=FALSE, row.names=TRUE)

# Summary statistics
sig <- results$table$FDR < 0.05
cat("\n=== RESULTS ===\n")
cat("Total junctions tested:", nrow(results$table), "\n")
cat("Significant (FDR<0.05):", sum(sig), "\n")
cat("  Up in KD:", sum(sig & results$table$logFC > 0), "\n")
cat("  Down in KD:", sum(sig & results$table$logFC < 0), "\n")
cat("\nOutput: data/leafcutter/edgeR_diff_junctions.txt\n")
