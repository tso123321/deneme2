############################################################
# 08_GO_enrichment.R
# RNA-seq Analysis Course
# Step 8: GO enrichment analysis
############################################################

############################################################
# 1. Amaç
############################################################

# Bu scriptin amacı:
# 1. Diferansiyel ekspresyon sonuçlarını okumak
# 2. Anlamlı genleri seçmek
# 3. Ensembl gene ID'lerini GO analizi için kullanmak
# 4. Biological Process (BP) GO enrichment analizi yapmak
# 5. Sonuçları tablo ve grafik olarak kaydetmektir
#
# GO enrichment analizi, anlamlı genlerin hangi biyolojik süreçlerde
# zenginleştiğini anlamak için kullanılır.

############################################################
# 2. Paketler
############################################################

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

############################################################
# 3. Çıktı klasörü
############################################################

dir.create("results/RNAseq/08_GO_enrichment", recursive = TRUE, showWarnings = FALSE)

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
# 5. Anlamlı genleri seçme
############################################################

sig_genes <- de_results %>%
  filter(FDR < 0.05, abs(logFC) >= 1)

cat("Significant genes:", nrow(sig_genes), "\n")

############################################################
# 6. Ensembl ID temizleme
############################################################

# Bazı Ensembl ID'lerinde versiyon olabilir:
# ENSG000001234.5 gibi
# GO analizi için versiyon kısmını kaldırıyoruz.

ensembl_ids <- sig_genes$Geneid
ensembl_ids <- gsub("\\..*", "", ensembl_ids)
ensembl_ids <- unique(ensembl_ids)

head(ensembl_ids)

############################################################
# 7. GO enrichment analizi
############################################################

ego_bp <- enrichGO(
  gene = ensembl_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_bp

############################################################
# 8. Sonuçları tabloya çevirme
############################################################

go_results <- as.data.frame(ego_bp)

head(go_results)

write.csv(
  go_results,
  "results/RNAseq/08_GO_enrichment/GO_BP_enrichment_results.csv",
  row.names = FALSE
)

############################################################
# 9. Dotplot
############################################################

if (nrow(go_results) > 0) {
  
  p_dot <- dotplot(
    ego_bp,
    showCategory = 20
  ) +
    ggtitle("GO Biological Process enrichment")
  
  ggsave(
    "results/RNAseq/08_GO_enrichment/GO_BP_dotplot.pdf",
    plot = p_dot,
    width = 8,
    height = 6
  )
  
  ggsave(
    "results/RNAseq/08_GO_enrichment/GO_BP_dotplot.png",
    plot = p_dot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

############################################################
# 10. Barplot
############################################################

if (nrow(go_results) > 0) {
  
  p_bar <- barplot(
    ego_bp,
    showCategory = 20
  ) +
    ggtitle("GO Biological Process enrichment")
  
  ggsave(
    "results/RNAseq/08_GO_enrichment/GO_BP_barplot.pdf",
    plot = p_bar,
    width = 8,
    height = 6
  )
  
  ggsave(
    "results/RNAseq/08_GO_enrichment/GO_BP_barplot.png",
    plot = p_bar,
    width = 8,
    height = 6,
    dpi = 300
  )
}

############################################################
# 11. GO için up/down ayrı analiz opsiyonu
############################################################

# İstersen anlamlı genleri yukarı ve aşağı regüle olan genler olarak
# ayrıca analiz edebiliriz.
#
# Şimdilik tüm anlamlı genler birlikte analiz edildi.

############################################################
# 12. Özet
############################################################

cat("GO enrichment completed successfully.\n")
cat("Input significant genes:", length(ensembl_ids), "\n")
cat("Significant GO terms:", nrow(go_results), "\n")

