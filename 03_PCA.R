############################################################
# 03_PCA.R
# RNA-seq Analysis Course
# Step 3: Principal Component Analysis (PCA)
############################################################

############################################################
# 1. Amaç
############################################################

# Bu scriptin amacı:
# 1. 01_preprocessing adımından gelen filtrelenmiş CPM tablosunu okumak
# 2. Metadata bilgisini okumak
# 3. CPM değerlerine log2(CPM + 1) dönüşümü yapmak
# 4. PCA analizi yapmak
# 5. Örneklerin PCA plot üzerinde gruplanmasını incelemek
# 6. PCA sonuçlarını ve grafiği kaydetmektir
#
# PCA, örnekler arasındaki genel benzerlik/farklılık yapısını görmek için kullanılır.
# RNA-seq analizlerinde kalite kontrol amaçlı çok sık kullanılır.
#
# Eğer aynı gruptaki biyolojik tekrarlar birbirine yakın,
# farklı gruplar birbirinden ayrılmış görünüyorsa,
# bu durum genellikle iyi bir biyolojik/teknik tutarlılık göstergesi olabilir.

############################################################
# 2. Paketler
############################################################

library(ggplot2)
library(dplyr)

############################################################
# 3. Çıktı klasörü
############################################################

dir.create("results/RNAseq/03_PCA", recursive = TRUE, showWarnings = FALSE)

############################################################
# 4. Filtrelenmiş CPM tablosunu okuma
############################################################

# PCA için count yerine CPM tablosu kullanıyoruz.
# Çünkü CPM, örnekler arası library size farkını azaltır.

cpm_filtered <- read.csv(
  "results/RNAseq/01_preprocessing/cpm_filtered.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

head(cpm_filtered)
str(cpm_filtered)

############################################################
# 5. Metadata okuma
############################################################

# Metadata örneklerin hangi gruba ait olduğunu gösterir.
# 01_preprocessing scriptinde oluşturulmuştu.

metadata <- read.csv(
  "results/RNAseq/01_preprocessing/metadata.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata$condition <- factor(metadata$condition, levels = c("SC", "Ssh1"))

metadata

############################################################
# 6. Gen bilgisi ve sample kolonlarını ayırma
############################################################

# İlk 3 kolon gen bilgisi:
# Geneid
# GeneName
# NewName
#
# Geri kalan kolonlar örneklerin CPM değerleridir.

gene_info_cols <- c("Geneid", "GeneName", "NewName")

sample_cols <- setdiff(colnames(cpm_filtered), gene_info_cols)

sample_cols

# Metadata sample sırası ile CPM tablosundaki sample sırası aynı mı kontrol edelim.

all(sample_cols == metadata$sample_id)

# Eğer FALSE çıkarsa, metadata'yı CPM tablosunun kolon sırasına göre düzenleyelim.

metadata <- metadata[match(sample_cols, metadata$sample_id), ]

all(sample_cols == metadata$sample_id)

############################################################
# 7. CPM matrix oluşturma
############################################################

# PCA için sadece sayısal ekspresyon değerleri gerekir.

cpm_matrix <- cpm_filtered[, sample_cols]

cpm_matrix <- as.matrix(cpm_matrix)

rownames(cpm_matrix) <- cpm_filtered$Geneid

head(cpm_matrix)

############################################################
# 8. Log dönüşümü
############################################################

# RNA-seq ekspresyon verileri genellikle çok geniş bir aralıkta dağılır.
# Bazı genler çok yüksek, bazıları çok düşük ifade edilir.
#
# log2(CPM + 1) dönüşümü:
# - büyük değerlerin etkisini azaltır
# - veriyi PCA için daha dengeli hale getirir
# - +1, log2(0) problemini engeller

log_cpm <- log2(cpm_matrix + 1)

head(log_cpm)

############################################################
# 9. PCA analizi
############################################################

# prcomp() R'da PCA yapmak için kullanılan temel fonksiyondur.
#
# PCA örnekler üzerinde yapılacağı için matrisi transpoze ediyoruz:
# Orijinal matrix:
# satırlar = genler
# kolonlar = örnekler
#
# PCA için istenen yapı:
# satırlar = örnekler
# kolonlar = genler
#
# Bu yüzden t(log_cpm) kullanılır.

pca <- prcomp(
  t(log_cpm),
  center = TRUE,
  scale. = FALSE
)

# center = TRUE:
# Her genin ortalamasını sıfıra çeker.
#
# scale. = FALSE:
# Her geni standart sapmaya bölmez.
# RNA-seq PCA'da log-normalize değerler üzerinde çoğu zaman scale = FALSE kullanılır.

############################################################
# 10. PCA varyans yüzdeleri
############################################################

# PCA sonucunda her principal component'in ne kadar varyans açıkladığını hesaplıyoruz.

percent_variance <- (pca$sdev^2) / sum(pca$sdev^2) * 100

percent_variance

pc1_label <- paste0("PC1 (", round(percent_variance[1], 2), "%)")
pc2_label <- paste0("PC2 (", round(percent_variance[2], 2), "%)")

pc1_label
pc2_label

############################################################
# 11. PCA sonuç tablosu oluşturma
############################################################

# pca$x, her örneğin PC skorlarını içerir.
# Bunu data frame'e çeviriyoruz ve metadata ile birleştiriyoruz.

pca_df <- as.data.frame(pca$x)

pca_df$sample_id <- rownames(pca_df)

pca_df <- pca_df %>%
  left_join(metadata, by = "sample_id")

pca_df

############################################################
# 12. PCA plot
############################################################

# ggplot ile PC1 ve PC2 eksenlerinde örnekleri çiziyoruz.
#
# x = PC1
# y = PC2
# color = condition
#
# Böylece aynı gruptaki örneklerin birbirine yakın olup olmadığını görebiliriz.

p_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = condition, label = sample_id)
) +
  geom_point(size = 4) +
  geom_text(
    vjust = -0.8,
    size = 3,
    color = "black"
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "PCA plot of RNA-seq samples",
    x = pc1_label,
    y = pc2_label,
    color = "Condition"
  )

p_pca

############################################################
# 13. PCA çıktıları kaydetme
############################################################

write.csv(
  pca_df,
  "results/RNAseq/03_PCA/PCA_coordinates.csv",
  row.names = FALSE
)

write.csv(
  data.frame(
    PC = paste0("PC", seq_along(percent_variance)),
    percent_variance = percent_variance
  ),
  "results/RNAseq/03_PCA/PCA_variance_explained.csv",
  row.names = FALSE
)

ggsave(
  filename = "results/RNAseq/03_PCA/PCA_plot.pdf",
  plot = p_pca,
  width = 6,
  height = 5
)

ggsave(
  filename = "results/RNAseq/03_PCA/PCA_plot.png",
  plot = p_pca,
  width = 6,
  height = 5,
  dpi = 300
)

############################################################
# 14. Özet
############################################################

cat("PCA completed successfully.\n")
cat("Number of genes used:", nrow(log_cpm), "\n")
cat("Number of samples:", ncol(log_cpm), "\n")
cat("PC1 variance:", round(percent_variance[1], 2), "%\n")
cat("PC2 variance:", round(percent_variance[2], 2), "%\n")
