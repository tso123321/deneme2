############################################################
# 09_KEGG_enrichment.R
# RNA-seq Analysis Course
# Step 9: KEGG enrichment analysis
############################################################

############################################################
# 1. Amaç
############################################################

# Bu scriptin amacı:
# 1. Diferansiyel ekspresyon sonuçlarını okumak
# 2. Anlamlı genleri seçmek
# 3. Gene symbol bilgisini Entrez ID'ye çevirmek
# 4. KEGG pathway enrichment analizi yapmak
# 5. Sonuçları tablo ve grafik olarak kaydetmektir
#
# KEGG enrichment analizi, anlamlı genlerin hangi biyolojik yolaklarda
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

dir.create("results/RNAseq/09_KEGG_enrichment", recursive = TRUE, showWarnings = FALSE)

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

############################################################
# 5. Anlamlı genleri seçme
############################################################

# KEGG analizi için FDR < 0.05 ve |logFC| >= 1 olan genleri kullanıyoruz.

sig_genes <- de_results %>%
  filter(FDR < 0.05, abs(logFC) >= 1)

cat("Significant genes:", nrow(sig_genes), "\n")

############################################################
# 6. Gene symbol temizleme
############################################################

# KEGG analizi Entrez ID kullanır.
# Bu nedenle GeneName kolonundaki gene symbol değerlerini Entrez ID'ye çevireceğiz.
#
# Önce boş veya NA GeneName kayıtlarını çıkarıyoruz.

gene_symbols <- sig_genes$GeneName

gene_symbols <- gene_symbols[!is.na(gene_symbols)]
gene_symbols <- gene_symbols[gene_symbols != ""]
gene_symbols <- trimws(gene_symbols)
gene_symbols <- unique(gene_symbols)

head(gene_symbols)

cat("Unique gene symbols:", length(gene_symbols), "\n")

############################################################
# 7. Gene symbol → Entrez ID dönüşümü
############################################################

# bitr(), clusterProfiler içindeki ID dönüştürme fonksiyonudur.
#
# fromType = "SYMBOL"
# toType   = "ENTREZID"
# OrgDb    = org.Hs.eg.db
#
# İnsan genleri için org.Hs.eg.db kullanılır.

gene_conversion <- bitr(
  gene_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

head(gene_conversion)

cat("Mapped genes:", nrow(gene_conversion), "\n")

# Entrez ID listesini alıyoruz.

entrez_ids <- unique(gene_conversion$ENTREZID)

head(entrez_ids)

############################################################
# 8. KEGG enrichment analizi
############################################################

# enrichKEGG(), KEGG pathway enrichment analizi yapar.
#
# organism = "hsa" insan için kullanılır.
# pvalueCutoff = 0.05 anlamlılık eşiğidir.
#
# Not:
# KEGG internet bağlantısı gerektirebilir.
# Codespaces veya lokal ortamda internet/kaynak erişimine bağlı olarak
# sonuç alınıp alınmaması değişebilir.

kegg_result <- enrichKEGG(
  gene = entrez_ids,
  organism = "hsa",
  pvalueCutoff = 0.05
)

kegg_result

############################################################
# 9. KEGG sonuçlarını readable hale getirme
############################################################

# setReadable(), Entrez ID'leri tekrar gene symbol olarak göstermeye çalışır.
# Böylece sonuç tabloları daha okunabilir olur.

kegg_result <- setReadable(
  kegg_result,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

############################################################
# 10. Sonuçları tabloya çevirme
############################################################

kegg_table <- as.data.frame(kegg_result)

head(kegg_table)

write.csv(
  kegg_table,
  "results/RNAseq/09_KEGG_enrichment/KEGG_enrichment_results.csv",
  row.names = FALSE
)

############################################################
# 11. Dotplot
############################################################

if (nrow(kegg_table) > 0) {
  
  p_dot <- dotplot(
    kegg_result,
    showCategory = 20
  ) +
    ggtitle("KEGG pathway enrichment")
  
  ggsave(
    "results/RNAseq/09_KEGG_enrichment/KEGG_dotplot.pdf",
    plot = p_dot,
    width = 8,
    height = 6
  )
  
  ggsave(
    "results/RNAseq/09_KEGG_enrichment/KEGG_dotplot.png",
    plot = p_dot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

############################################################
# 12. Barplot
############################################################

if (nrow(kegg_table) > 0) {
  
  p_bar <- barplot(
    kegg_result,
    showCategory = 20
  ) +
    ggtitle("KEGG pathway enrichment")
  
  ggsave(
    "results/RNAseq/09_KEGG_enrichment/KEGG_barplot.pdf",
    plot = p_bar,
    width = 8,
    height = 6
  )
  
  ggsave(
    "results/RNAseq/09_KEGG_enrichment/KEGG_barplot.png",
    plot = p_bar,
    width = 8,
    height = 6,
    dpi = 300
  )
}

############################################################
# 13. ID dönüşüm tablosunu kaydetme
############################################################

write.csv(
  gene_conversion,
  "results/RNAseq/09_KEGG_enrichment/gene_symbol_to_entrez_mapping.csv",
  row.names = FALSE
)

############################################################
# 14. Özet
############################################################

cat("KEGG enrichment completed successfully.\n")
cat("Input gene symbols:", length(gene_symbols), "\n")
cat("Mapped Entrez IDs:", length(entrez_ids), "\n")
cat("Significant KEGG pathways:", nrow(kegg_table), "\n")

