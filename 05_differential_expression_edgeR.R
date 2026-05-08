############################################################
# 05_differential_expression.R
# RNA-seq Analysis Course
# Step 5: Differential expression analysis with edgeR
############################################################

library(edgeR)
library(dplyr)

dir.create(
  "results/RNAseq/05_differential_expression",
  recursive = TRUE,
  showWarnings = FALSE
)

count_annotated <- read.csv(
  "results/RNAseq/02_annotation/count_annotated.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata <- read.csv(
  "results/RNAseq/01_preprocessing/metadata.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata$condition <- factor(metadata$condition, levels = c("SC", "Ssh1"))

sample_cols <- metadata$sample_id

missing_samples <- setdiff(sample_cols, colnames(count_annotated))

if (length(missing_samples) > 0) {
  stop(
    "count_annotated tablosunda şu sample kolonları bulunamadı: ",
    paste(missing_samples, collapse = ", ")
  )
}

count_matrix <- count_annotated[, sample_cols]
count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "numeric"
rownames(count_matrix) <- count_annotated$Geneid

metadata <- metadata[match(colnames(count_matrix), metadata$sample_id), ]

if (!all(colnames(count_matrix) == metadata$sample_id)) {
  stop("Count matrix kolonları metadata sample_id ile eşleşmiyor.")
}

if (any(is.na(metadata$condition))) {
  stop("Metadata içinde condition bilgisi eksik olan sample var.")
}

group <- metadata$condition

if (any(is.na(count_matrix))) {
  stop("Count matrix içinde NA değer var. Muhtemelen numeric olmayan kolonlar matrix'e girdi.")
}

if (any(count_matrix < 0)) {
  stop("Count matrix içinde negatif değer var. Count data negatif olamaz.")
}

if (!all(count_matrix == round(count_matrix))) {
  warning("Count matrix integer değil. edgeR count data bekler.")
}

gene_info_cols <- intersect(
  c("Geneid", "GeneName", "NewName", "GeneType", "chr", "start", "end", "strand"),
  colnames(count_annotated)
)

gene_info <- count_annotated[, gene_info_cols]

dge <- DGEList(
  counts = count_matrix,
  group = group,
  genes = gene_info
)

dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge, verbose = TRUE)
dge <- estimateTagwiseDisp(dge)

test <- exactTest(
  dge,
  pair = c("SC", "Ssh1")
)

de_results <- topTags(
  test,
  n = Inf
)

de_results <- as.data.frame(de_results)

de_results <- de_results %>%
  mutate(
    significance = case_when(
      FDR < 0.05 & logFC >= 1 ~ "Up_in_Ssh1",
      FDR < 0.05 & logFC <= -1 ~ "Down_in_Ssh1",
      FDR < 0.05 ~ "FDR_significant",
      TRUE ~ "Not_significant"
    )
  )

write.csv(
  de_results,
  "results/RNAseq/05_differential_expression/DE_edgeR_all_results.csv",
  row.names = FALSE
)

de_significant <- de_results %>%
  filter(FDR < 0.05)

write.csv(
  de_significant,
  "results/RNAseq/05_differential_expression/DE_edgeR_significant_FDR005.csv",
  row.names = FALSE
)

de_significant_logfc <- de_results %>%
  filter(FDR < 0.05, abs(logFC) >= 1)

write.csv(
  de_significant_logfc,
  "results/RNAseq/05_differential_expression/DE_edgeR_significant_FDR005_logFC1.csv",
  row.names = FALSE
)

de_summary <- data.frame(
  total_genes_tested = nrow(de_results),
  significant_FDR_005 = nrow(de_significant),
  significant_FDR_005_logFC1 = nrow(de_significant_logfc),
  up_in_Ssh1 = sum(de_results$significance == "Up_in_Ssh1"),
  down_in_Ssh1 = sum(de_results$significance == "Down_in_Ssh1"),
  fdr_significant_only = sum(de_results$significance == "FDR_significant")
)

write.csv(
  de_summary,
  "results/RNAseq/05_differential_expression/DE_summary.csv",
  row.names = FALSE
)

cat("Differential expression analysis completed successfully.\n")
cat("Total genes tested:", nrow(de_results), "\n")
cat("Significant genes FDR < 0.05:", nrow(de_significant), "\n")
cat("Significant genes FDR < 0.05 and |logFC| >= 1:", nrow(de_significant_logfc), "\n")
cat("Up in Ssh1:", sum(de_results$significance == "Up_in_Ssh1"), "\n")
cat("Down in Ssh1:", sum(de_results$significance == "Down_in_Ssh1"), "\n")

