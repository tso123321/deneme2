############################################################
# 04_sample_correlation_heatmap.R
# RNA-seq Analysis Course
# Step 4: Sample-to-sample correlation heatmap
############################################################

############################################################
# 1. Amaç
############################################################

# Bu scriptin amacı:
# 1. Filtrelenmiş CPM tablosunu okumak
# 2. log2(CPM + 1) dönüşümü yapmak
# 3. Örnekler arası korelasyon hesaplamak
# 4. Korelasyon heatmap çizmek
# 5. Heatmap çıktısını kaydetmektir
#
# RNA-seq analizlerinde sample correlation heatmap,
# örneklerin birbirine ne kadar benzediğini göstermek için kullanılır.
#
# Aynı gruptaki biyolojik tekrarların birbirine daha yüksek korelasyon göstermesi beklenir.

############################################################
# 2. Paketler
############################################################

library(pheatmap)
library(RColorBrewer)

############################################################
# 3. Çıktı klasörü
############################################################

dir.create("results/RNAseq/04_sample_correlation_heatmap", recursive = TRUE, showWarnings = FALSE)

############################################################
# 4. Filtrelenmiş CPM tablosunu okuma
############################################################

cpm_filtered <- read.csv(
  "results/RNAseq/01_preprocessing/cpm_filtered.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

head(cpm_filtered)

############################################################
# 5. Metadata okuma
############################################################

metadata <- read.csv(
  "results/RNAseq/01_preprocessing/metadata.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata$condition <- factor(metadata$condition, levels = c("SC", "Ssh1"))

metadata

############################################################
# 6. Sample kolonlarını belirleme
############################################################

gene_info_cols <- c("Geneid", "GeneName", "NewName")

sample_cols <- setdiff(colnames(cpm_filtered), gene_info_cols)

sample_cols

# Metadata sırasını sample kolonlarına göre düzenleyelim.

metadata <- metadata[match(sample_cols, metadata$sample_id), ]

all(sample_cols == metadata$sample_id)

############################################################
# 7. CPM matrix oluşturma
############################################################

cpm_matrix <- cpm_filtered[, sample_cols]
cpm_matrix <- as.matrix(cpm_matrix)
rownames(cpm_matrix) <- cpm_filtered$Geneid

head(cpm_matrix)

############################################################
# 8. Log dönüşümü
############################################################

# Korelasyon hesaplamadan önce log2(CPM + 1) dönüşümü yapıyoruz.
# Bunun nedeni RNA-seq ekspresyon değerlerinin geniş aralıkta dağılmasıdır.

log_cpm <- log2(cpm_matrix + 1)

head(log_cpm)

############################################################
# 9. Sample correlation hesaplama
############################################################

# cor() fonksiyonu korelasyon matrisi hesaplar.
#
# Burada kolonlar sample olduğu için cor(log_cpm) doğrudan
# sample-sample correlation verir.
#
# method = "spearman":
# sıralamaya dayalı korelasyondur, uç değerlere Pearson'a göre daha dayanıklıdır.

sample_correlation <- cor(
  log_cpm,
  method = "spearman"
)

sample_correlation

############################################################
# 10. Annotation oluşturma
############################################################

# Heatmap üzerinde örneklerin hangi gruba ait olduğunu göstermek için
# annotation_col tablosu oluşturuyoruz.
#
# pheatmap için rownames(annotation_col) sample isimleriyle aynı olmalıdır.

annotation_col <- data.frame(
  condition = metadata$condition
)

rownames(annotation_col) <- metadata$sample_id

annotation_col

############################################################
# 11. Renk paleti
############################################################

# colorRampPalette iki veya daha fazla rengi kullanarak ara renkler üretir.
# Burada korelasyon değerleri için mavi tonlarında bir palette kullanıyoruz.

heatmap_colors <- colorRampPalette(
  rev(c("#08306b", "#2171b5", "#6baed6", "#deebf7"))
)(100)

############################################################
# 12. Heatmap çizme
############################################################

pheatmap(
  sample_correlation,
  color = heatmap_colors,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Sample-to-sample correlation heatmap",
  filename = "results/RNAseq/04_sample_correlation_heatmap/sample_correlation_heatmap.pdf",
  width = 7,
  height = 6
)

############################################################
# 13. Korelasyon tablosunu kaydetme
############################################################

write.csv(
  sample_correlation,
  "results/RNAseq/04_sample_correlation_heatmap/sample_correlation_matrix.csv"
)

############################################################
# 14. PNG olarak da kaydetme
############################################################

png(
  filename = "results/RNAseq/04_sample_correlation_heatmap/sample_correlation_heatmap.png",
  width = 1800,
  height = 1500,
  res = 300
)

pheatmap(
  sample_correlation,
  color = heatmap_colors,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Sample-to-sample correlation heatmap"
)

dev.off()

############################################################
# 15. Özet
############################################################

cat("Sample correlation heatmap completed successfully.\n")
cat("Number of samples:", ncol(log_cpm), "\n")
cat("Output folder: results/RNAseq/04_sample_correlation_heatmap\n")

