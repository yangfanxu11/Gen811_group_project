#!/usr/bin/env Rscript

library(ggplot2)

# 读取数据
results <- read.table("data/leafcutter/edgeR_diff_junctions.txt",
                      header=TRUE, sep="\t", row.names=1)
sig <- results[results$FDR < 0.05, ]
sig$junction <- rownames(sig)
sig$cluster <- sub(".*:(clu_[0-9]+)_.*", "\\1", sig$junction)

# 分类
cluster_counts <- table(sig$cluster)
sig$type <- ifelse(sig$cluster %in% names(cluster_counts[cluster_counts > 1]), 
                   "AS", "TSS/TES")

# 进一步细分AS
competitive_clusters <- character()
for (clu in names(cluster_counts[cluster_counts > 1])) {
  junc <- sig[sig$cluster == clu, ]
  if (sum(junc$logFC > 0) > 0 && sum(junc$logFC < 0) > 0) {
    competitive_clusters <- c(competitive_clusters, clu)
  }
}
sig$subtype <- sig$type
sig$subtype[sig$cluster %in% competitive_clusters] <- "Competitive AS"
sig$subtype[sig$type == "AS" & !(sig$cluster %in% competitive_clusters)] <- "Consistent AS"

dir.create("results/figures", showWarnings=FALSE)

# ============================================
# 图1: AS vs TSS/TES分类
# ============================================

cat("Generating Figure 1: AS vs TSS/TES classification...\n")

# 统计
type_counts <- data.frame(
  Category = c("Alternative Splicing", "Alternative TSS/TES"),
  Count = c(sum(sig$type == "AS"), sum(sig$type == "TSS/TES")),
  Percent = c(
    sum(sig$type == "AS") / nrow(sig) * 100,
    sum(sig$type == "TSS/TES") / nrow(sig) * 100
  )
)

# 柱状图
pdf("results/figures/fig1_AS_vs_TSS.pdf", width=8, height=6)
ggplot(type_counts, aes(x=Category, y=Count, fill=Category)) +
  geom_bar(stat="identity", width=0.6) +
  geom_text(aes(label=paste0(Count, "\n(", sprintf("%.1f", Percent), "%)")),
            vjust=-0.3, size=5) +
  scale_fill_manual(values=c("Alternative Splicing"="#E64B35", 
                             "Alternative TSS/TES"="#4DBBD5")) +
  labs(title="Classification of 1,791 Significant Junctions",
       x="", y="Number of Junctions") +
  theme_classic(base_size=14) +
  theme(legend.position="none",
        plot.title=element_text(hjust=0.5, face="bold"))
dev.off()

# 带细分的版本
subtype_counts <- data.frame(
  Category = c("Competitive AS", "Consistent AS", "TSS/TES"),
  Count = c(
    sum(sig$subtype == "Competitive AS"),
    sum(sig$subtype == "Consistent AS"),
    sum(sig$subtype == "TSS/TES")
  )
)
subtype_counts$Percent <- subtype_counts$Count / nrow(sig) * 100

pdf("results/figures/fig1_detailed.pdf", width=10, height=6)
ggplot(subtype_counts, aes(x=reorder(Category, -Count), y=Count, fill=Category)) +
  geom_bar(stat="identity", width=0.6) +
  geom_text(aes(label=paste0(Count, "\n(", sprintf("%.1f", Percent), "%)")),
            vjust=-0.3, size=4) +
  scale_fill_manual(values=c("Competitive AS"="#E64B35",
                             "Consistent AS"="#F39B7F",
                             "TSS/TES"="#4DBBD5")) +
  labs(title="Detailed Classification of Significant Junctions",
       subtitle="Competitive AS = exon inclusion/skipping; Consistent AS = coordinated regulation",
       x="", y="Number of Junctions") +
  theme_classic(base_size=14) +
  theme(legend.position="none",
        plot.title=element_text(hjust=0.5, face="bold"),
        plot.subtitle=element_text(hjust=0.5, size=10))
dev.off()

# ============================================
# 图2: Cluster复杂度分布
# ============================================

cat("Generating Figure 2: Cluster complexity distribution...\n")

complexity <- data.frame(
  n_junctions = as.numeric(cluster_counts),
  cluster = names(cluster_counts)
)

pdf("results/figures/fig2_cluster_complexity.pdf", width=10, height=6)
ggplot(complexity, aes(x=n_junctions)) +
  geom_histogram(binwidth=1, fill="#00A087", color="white", alpha=0.8) +
  geom_vline(xintercept=1, linetype="dashed", color="red", size=1) +
  annotate("text", x=1.5, y=max(table(complexity$n_junctions))*0.9,
           label=paste0("Single junction\n(TSS/TES)\n", sum(cluster_counts==1), " clusters"),
           hjust=0, color="red", size=4) +
  annotate("text", x=2.5, y=max(table(complexity$n_junctions))*0.7,
           label=paste0("Multi-junction\n(AS)\n", sum(cluster_counts>1), " clusters"),
           hjust=0, color="#00A087", size=4) +
  labs(title="Cluster Complexity Distribution",
       subtitle=paste0("Total: ", length(cluster_counts), " clusters with significant junctions"),
       x="Number of Significant Junctions per Cluster",
       y="Number of Clusters") +
  scale_x_continuous(breaks=1:max(complexity$n_junctions)) +
  theme_classic(base_size=14) +
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        plot.subtitle=element_text(hjust=0.5))
dev.off()

# ============================================
# 图4: LogFC分布对比 (AS vs TSS/TES)
# ============================================

cat("Generating Figure 4: LogFC distribution comparison...\n")

pdf("results/figures/fig4_logFC_comparison.pdf", width=10, height=6)

# 密度图
ggplot(sig, aes(x=logFC, fill=type)) +
  geom_density(alpha=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color="black") +
  scale_fill_manual(values=c("AS"="#E64B35", "TSS/TES"="#4DBBD5"),
                    name="Category") +
  labs(title="Effect Size Distribution: AS vs TSS/TES",
       x="log2 Fold Change (KD/WT)",
       y="Density") +
  theme_classic(base_size=14) +
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        legend.position=c(0.85, 0.85))

dev.off()

# 箱线图版本
pdf("results/figures/fig4_logFC_boxplot.pdf", width=8, height=6)

sig$direction <- ifelse(sig$logFC > 0, "Up in KD", "Down in KD")

ggplot(sig, aes(x=type, y=abs(logFC), fill=type)) +
  geom_boxplot(alpha=0.7, outlier.alpha=0.3) +
  geom_jitter(width=0.2, alpha=0.1, size=0.5) +
  scale_fill_manual(values=c("AS"="#E64B35", "TSS/TES"="#4DBBD5")) +
  labs(title="Effect Size Comparison",
       x="", y="|log2 Fold Change|") +
  theme_classic(base_size=14) +
  theme(plot.title=element_text(hjust=0.5, face="bold"),
        legend.position="none") +
  stat_summary(fun=median, geom="text", aes(label=sprintf("%.2f", ..y..)),
               vjust=-1, size=4, color="black")

dev.off()

# 方向性对比
direction_summary <- aggregate(junction ~ type + direction, sig, length)
colnames(direction_summary)[3] <- "count"

pdf("results/figures/fig4_direction.pdf", width=8, height=6)

ggplot(direction_summary, aes(x=type, y=count, fill=direction)) +
  geom_bar(stat="identity", position="fill") +
  geom_text(aes(label=count), position=position_fill(vjust=0.5), 
            color="white", size=5, fontface="bold") +
  scale_fill_manual(values=c("Up in KD"="#E64B35", "Down in KD"="#4DBBD5"),
                    name="Direction") +
  scale_y_continuous(labels=scales::percent) +
  labs(title="Junction Direction: AS vs TSS/TES",
       x="", y="Percentage") +
  theme_classic(base_size=14) +
  theme(plot.title=element_text(hjust=0.5, face="bold"))

dev.off()

cat("\n=== Figures generated ===\n")
cat("results/figures/fig1_AS_vs_TSS.pdf - Main classification\n")
cat("results/figures/fig1_detailed.pdf - With AS subtypes\n")
cat("results/figures/fig2_cluster_complexity.pdf - Cluster distribution\n")
cat("results/figures/fig4_logFC_comparison.pdf - Effect size density\n")
cat("results/figures/fig4_logFC_boxplot.pdf - Effect size boxplot\n")
cat("results/figures/fig4_direction.pdf - Up/Down comparison\n")

