#!/usr/bin/env Rscript

# 正确读取：从edgeR结果开始
results <- read.table("data/leafcutter/edgeR_diff_junctions.txt",
                      header=TRUE, sep="\t", row.names=1)

# 显著的
sig <- results[results$FDR < 0.05, ]
sig$junction <- rownames(sig)

cat("Total significant junctions:", nrow(sig), "\n\n")

# 提取cluster
sig$cluster <- sub(".*:(clu_[0-9]+)_.*", "\\1", sig$junction)

# 每个cluster有多少显著junction
cluster_counts <- table(sig$cluster)

single_clusters <- names(cluster_counts[cluster_counts == 1])
multi_clusters <- names(cluster_counts[cluster_counts > 1])

cat("=== AS vs TSS/TES 分类 ===\n\n")
cat("【Cluster复杂度】\n")
cat("Single junction clusters:", length(single_clusters), "\n")
cat("Multi-junction clusters:", length(multi_clusters), "\n\n")

# 分析multi-junction clusters
if (length(multi_clusters) > 0) {
  cat("【Multi-junction clusters = 确定的AS】\n")
  
  competitive <- 0
  consistent_up <- 0
  consistent_down <- 0
  
  for (clu in multi_clusters) {
    junc_in_clu <- sig[sig$cluster == clu, ]
    n_up <- sum(junc_in_clu$logFC > 0)
    n_down <- sum(junc_in_clu$logFC < 0)
    
    if (n_up > 0 && n_down > 0) {
      competitive <- competitive + 1
    } else if (n_up > 0) {
      consistent_up <- consistent_up + 1
    } else {
      consistent_down <- consistent_down + 1
    }
  }
  
  cat("Competitive (有up也有down):", competitive, "\n")
  cat("  → Exon inclusion/skipping\n")
  cat("Consistent up:", consistent_up, "\n")
  cat("Consistent down:", consistent_down, "\n")
  cat("  → Coordinated AS or TSS/TES\n\n")
  
  # 展示几个竞争性AS的例子
  cat("【竞争性AS例子】\n")
  count <- 0
  for (clu in multi_clusters) {
    junc_in_clu <- sig[sig$cluster == clu, ]
    if (sum(junc_in_clu$logFC > 0) > 0 && sum(junc_in_clu$logFC < 0) > 0) {
      count <- count + 1
      if (count <= 5) {
        cat("\nCluster:", clu, "\n")
        cat("  Up:", sum(junc_in_clu$logFC > 0), 
            "Down:", sum(junc_in_clu$logFC < 0), "\n")
        print(junc_in_clu[, c("logFC", "FDR")])
      }
    }
  }
}

cat("\n【Single junction clusters】\n")
cat("Total:", length(single_clusters), "\n")
cat("可能是Alternative TSS/TES\n\n")

# 汇总
as_junctions <- sum(sig$cluster %in% multi_clusters)
ambiguous_junctions <- sum(sig$cluster %in% single_clusters)

cat("=== 最终答案 ===\n")
cat("Total significant junctions:", nrow(sig), "\n")
cat("确定是AS:", as_junctions, sprintf("(%.1f%%)\n", as_junctions/nrow(sig)*100))
cat("可能是TSS/TES:", ambiguous_junctions, sprintf("(%.1f%%)\n", ambiguous_junctions/nrow(sig)*100))

