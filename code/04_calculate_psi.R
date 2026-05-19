#!/usr/bin/env Rscript
# Calculate PSI values for significant junctions

library(data.table)

cat("Reading junction counts...\n")
counts <- read.table(gzfile("../data/leafcutter/leafcutter_perind.counts.gz"), 
                     header=TRUE, sep=" ", check.names=FALSE)

# Read edgeR results
results <- read.table("../data/leafcutter/edgeR_diff_junctions.txt", 
                      header=TRUE, sep="\t")
sig_junctions <- rownames(results)[results$FDR < 0.05]

# Extract counts and totals
extract_counts <- function(x) {
  parts <- strsplit(as.character(x), "/")
  data.frame(
    numerator = as.numeric(sapply(parts, `[`, 1)),
    denominator = as.numeric(sapply(parts, `[`, 2))
  )
}

# Calculate PSI for each sample
psi_data <- data.frame(junction = counts[,1])
for (col in colnames(counts)[-1]) {
  vals <- extract_counts(counts[[col]])
  psi_data[[col]] <- vals$numerator / vals$denominator * 100
}

# Filter to significant junctions
psi_sig <- psi_data[psi_data$junction %in% sig_junctions, ]

# Calculate mean PSI per group
psi_sig$WT_mean <- rowMeans(psi_sig[, c("WT1","WT2","WT3","WT4")], na.rm=TRUE)
psi_sig$KD_mean <- rowMeans(psi_sig[, c("KD1","KD2")], na.rm=TRUE)
psi_sig$deltaPSI <- psi_sig$KD_mean - psi_sig$WT_mean

write.table(psi_sig, "../data/leafcutter/psi_values.txt", 
            sep="\t", quote=FALSE, row.names=FALSE)

cat("PSI values calculated for", nrow(psi_sig), "significant junctions\n")
cat("Output: data/leafcutter/psi_values.txt\n")
