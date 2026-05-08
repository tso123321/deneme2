############################################################
# 06_volcano_MA_plots.R
# RNA-seq Analysis Course
# Step 6: Volcano plot and MA plot
############################################################

############################################################
# 1. Amaç
############################################################

# Bu scriptin amacı:
# 1. edgeR diferansiyel ekspresyon sonuçlarını okumak
# 2. Genleri anlamlılık durumuna göre sınıflandırmak
# 3. Volcano plot çizmek
# 4. MA plot çizmek
# 5. Grafikleri PDF ve PNG olarak kaydetmektir
#
# Volcano plot:
# x ekseni = log2 fold change
# y ekseni = -log10(FDR veya p-value)
#
# MA plot:
# x ekseni = ortalama ekspresyon
# y ekseni = log2 fold change
#
# Bu iki grafik diferansiyel ekspresyon sonuçlarının genel dağılımını
# hızlıca görmek için kullanılır.

############################################################
# 2. Paketler
############################################################

library(ggplot2)
library(dplyr)

############################################################
# 3. Çıktı klasörü
############################################################

dir.create("results/RNAseq/06_volcano_MA_plots", recursive = TRUE, showWarnings = FALSE)

############################################################
# 4. DE sonuçlarını okuma
############################################################

de_results <- read.csv(
  "results/RNAseq/05_differential_expression/DE_edgeR_all_results.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

head(de_results)
str(de_results)

############################################################
# 5. Gerekli kolonları kontrol etme
############################################################

# edgeR sonucu içinde beklenen temel kolonlar:
# logFC  : log2 fold change
# logCPM : ortalama ekspresyon seviyesi
# PValue : ham p-değeri
# FDR    : düzeltilmiş p-değeri

required_cols <- c("Geneid", "GeneName", "logFC", "logCPM", "PValue", "FDR")

missing_cols <- setdiff(required_cols, colnames(de_results))

if (length(missing_cols) > 0) {
  stop(
    paste(
      "Eksik kolon(lar) var:",
      paste(missing_cols, collapse = ", ")
    )
  )
}

############################################################
# 6. Anlamlılık etiketi oluşturma
############################################################

# Burada genleri 3 gruba ayırıyoruz:
#
# Up_in_Ssh1:
# FDR < 0.05 ve logFC >= 1
#
# Down_in_Ssh1:
# FDR < 0.05 ve logFC <= -1
#
# Not_significant:
# diğer tüm genler
#
# logFC pozitifse gen Ssh1 grubunda daha yüksek,
# logFC negatifse gen SC grubunda daha yüksek demektir.

de_results <- de_results %>%
  mutate(
    regulation = case_when(
      FDR < 0.05 & logFC >= 1 ~ "Up_in_Ssh1",
      FDR < 0.05 & logFC <= -1 ~ "Down_in_Ssh1",
      TRUE ~ "Not_significant"
    )
  )

table(de_results$regulation)

############################################################
# 7. Volcano plot için -log10(FDR) hesaplama
############################################################

# FDR değeri küçüldükçe gen daha anlamlı kabul edilir.
# -log10(FDR) dönüşümü küçük FDR değerlerini grafikte yukarı taşır.
#
# FDR = 0 gibi değerler varsa -log10(0) sonsuz olur.
# Bu yüzden çok küçük bir sayı ekliyoruz.

de_results <- de_results %>%
  mutate(
    neg_log10_FDR = -log10(FDR + 1e-300)
  )

head(de_results[, c("GeneName", "logFC", "FDR", "neg_log10_FDR", "regulation")])

############################################################
# 8. Volcano plot
############################################################

# Volcano plot:
# x = logFC
# y = -log10(FDR)
#
# Sağ üst: Ssh1'de yüksek ve anlamlı genler
# Sol üst: SC'de yüksek ve anlamlı genler
# Alt bölgeler: anlamlı olmayan genler

p_volcano <- ggplot(
  de_results,
  aes(x = logFC, y = neg_log10_FDR, color = regulation)
) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic(base_size = 13) +
  labs(
    title = "Volcano plot: Ssh1 vs SC",
    x = "log2 fold change (Ssh1 / SC)",
    y = "-log10(FDR)",
    color = "Regulation"
  ) +
  scale_color_manual(
    values = c(
      "Up_in_Ssh1" = "#D55E00",
      "Down_in_Ssh1" = "#0072B2",
      "Not_significant" = "grey70"
    )
  )

p_volcano

############################################################
# 9. MA plot
############################################################

# MA plot:
# M = logFC
# A = ortalama ekspresyon
#
# edgeR sonucunda logCPM, genin ortalama ekspresyon seviyesini temsil eder.
#
# MA plot, düşük ve yüksek ekspresyon seviyelerinde logFC dağılımını görmeyi sağlar.

p_ma <- ggplot(
  de_results,
  aes(x = logCPM, y = logFC, color = regulation)
) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_hline(yintercept = 0, linetype = "solid") +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
  theme_classic(base_size = 13) +
  labs(
    title = "MA plot: Ssh1 vs SC",
    x = "Average expression (logCPM)",
    y = "log2 fold change (Ssh1 / SC)",
    color = "Regulation"
  ) +
  scale_color_manual(
    values = c(
      "Up_in_Ssh1" = "#D55E00",
      "Down_in_Ssh1" = "#0072B2",
      "Not_significant" = "grey70"
    )
  )

p_ma

############################################################
# 10. En anlamlı genleri etiketleme için tablo
############################################################

# İsteğe bağlı olarak en anlamlı genleri ayrı tablo olarak kaydediyoruz.
# Burada FDR'a göre ilk 20 geni seçiyoruz.

top_genes <- de_results %>%
  arrange(FDR) %>%
  slice_head(n = 20)

write.csv(
  top_genes,
  "results/RNAseq/06_volcano_MA_plots/top20_genes_by_FDR.csv",
  row.names = FALSE
)

############################################################
# 11. Grafikleri kaydetme
############################################################

ggsave(
  filename = "results/RNAseq/06_volcano_MA_plots/volcano_plot.pdf",
  plot = p_volcano,
  width = 7,
  height = 6
)

ggsave(
  filename = "results/RNAseq/06_volcano_MA_plots/volcano_plot.png",
  plot = p_volcano,
  width = 7,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "results/RNAseq/06_volcano_MA_plots/MA_plot.pdf",
  plot = p_ma,
  width = 7,
  height = 6
)

ggsave(
  filename = "results/RNAseq/06_volcano_MA_plots/MA_plot.png",
  plot = p_ma,
  width = 7,
  height = 6,
  dpi = 300
)

############################################################
# 12. Özet
############################################################

cat("Volcano and MA plots completed successfully.\n")
cat("Total genes:", nrow(de_results), "\n")
cat("Up in Ssh1:", sum(de_results$regulation == "Up_in_Ssh1"), "\n")
cat("Down in Ssh1:", sum(de_results$regulation == "Down_in_Ssh1"), "\n")
cat("Not significant:", sum(de_results$regulation == "Not_significant"), "\n")

