############################################################
# 07_significant_gene_heatmap.R
# RNA-seq Analysis Course
# Step 7: Heatmap of significant genes
############################################################

############################################################
# 1. Amaç
############################################################

# Bu scriptin amacı:
# 1. Diferansiyel ekspresyon sonuçlarını okumak
# 2. Anlamlı genleri seçmek (FDR < 0.05, |logFC| >= 1)
# 3. Bu genlere ait ekspresyon değerlerini almak
# 4. log2(CPM + 1) kullanarak heatmap çizmek
#
# Heatmap, genlerin örnekler arasında nasıl değiştiğini görmeyi sağlar.

############################################################
# 2. Paketler
############################################################

library(pheatmap)
library(dplyr)

############################################################
# 3. Çıktı klasörü
############################################################

dir.create("results/RNAseq/07_significant_gene_heatmap", recursive = TRUE, showWarnings = FALSE)

############################################################
# 4. DE sonuçlarını okuma
############################################################

de_results <- read.csv(
  "results/RNAseq/05_differential_expression/DE_edgeR_all_results.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

############################################################
# 5. Filtrelenmiş CPM tablosunu okuma
############################################################

cpm_filtered <- read.csv(
  "results/RNAseq/01_preprocessing/cpm_filtered.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

############################################################
# 6. Metadata okuma
############################################################

metadata <- read.csv(
  "results/RNAseq/01_preprocessing/metadata.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

metadata$condition <- factor(metadata$condition, levels = c("SC", "Ssh1"))

############################################################
# 7. Anlamlı genleri seçme
############################################################

sig_genes <- de_results %>%
  filter(FDR < 0.05, abs(logFC) >= 1)

cat("Number of significant genes:", nrow(sig_genes), "\n")

############################################################
# 8. Sample kolonlarını belirleme
############################################################

gene_info_cols <- c("Geneid", "GeneName", "NewName")

sample_cols <- setdiff(colnames(cpm_filtered), gene_info_cols)

############################################################
# 9. CPM matrix oluşturma
############################################################

cpm_matrix <- cpm_filtered[, sample_cols]
cpm_matrix <- as.matrix(cpm_matrix)
rownames(cpm_matrix) <- cpm_filtered$Geneid

############################################################
# 10. Anlamlı genleri filtreleme
############################################################

# Sadece anlamlı genleri alıyoruz

cpm_sig <- cpm_matrix[rownames(cpm_matrix) %in% sig_genes$Geneid, ]

############################################################
# 11. Log dönüşümü
############################################################

log_cpm_sig <- log2(cpm_sig + 1)

############################################################
# 12. Z-score normalization
############################################################

# Her gen için ortalamayı 0, std sapmayı 1 yapar
# Heatmap'te patternleri daha iyi görmemizi sağlar

z_score <- t(scale(t(log_cpm_sig)))

############################################################
# 13. Metadata sıralama
############################################################

metadata <- metadata[match(colnames(z_score), metadata$sample_id), ]

annotation_col <- data.frame(
  condition = metadata$condition
)

rownames(annotation_col) <- metadata$sample_id

############################################################
# 14. Heatmap çizme
############################################################

pheatmap(
  z_score,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Significant genes heatmap",
  filename = "results/RNAseq/07_significant_gene_heatmap/significant_genes_heatmap.pdf",
  width = 8,
  height = 6
)

############################################################
# 15. PNG olarak kaydetme
############################################################

png(
  "results/RNAseq/07_significant_gene_heatmap/significant_genes_heatmap.png",
  width = 2000,
  height = 1500,
  res = 300
)

pheatmap(
  z_score,
  annotation_col = annotation_col,
  show_rownames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Significant genes heatmap"
)

dev.off()

############################################################
# 16. Özet
############################################################

cat("Heatmap completed successfully.\n")
cat("Genes used:", nrow(z_score), "\n")
cat("Samples:", ncol(z_score), "\n")

