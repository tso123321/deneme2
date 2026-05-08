############################################################
# 01_preprocessing_filtering.R
# RNA-seq Analysis Course
# Step 1: Count table reading, existing CPM table reading and filtering
############################################################

############################################################
# 1. Amaç
############################################################

# Bu scriptin amacı:
# 1. Ham count tablosunu okumak
# 2. Hazır CPM tablosunu okumak
# 3. Gen bilgisi kolonlarını ve örnek kolonlarını ayırmak
# 4. Hazır CPM değerleri üzerinden düşük ekspresyonlu genleri filtrelemek
# 5. Filtrelenmiş count ve CPM tablolarını kaydetmektir

# RNA-seq analizlerinde ham count tablosunda çok sayıda düşük/ifade edilmeyen gen bulunur.
# Bu genler diferansiyel ekspresyon analizinde gürültü oluşturabilir.
# Bu nedenle analiz öncesinde düşük count/CPM değerine sahip genler filtrelenir.

############################################################
# 2. Paketler
############################################################

# dplyr veri düzenleme için kullanılır.

library(dplyr)

dir.create("results/RNAseq", recursive = TRUE, showWarnings = FALSE)
dir.create("results/RNAseq/01_preprocessing", recursive = TRUE, showWarnings = FALSE)

# =========================
# 1) Input dosyalarını oku
# =========================

counts <- read.table(
  "data/counts_table.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cpm_table <- read.csv(
  "data/cpm_table.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# CPM tablosunda R'nin eklediği X prefix'ini düzelt
# Örn: X9_SC1 -> 9_SC1
colnames(cpm_table) <- gsub("^X(?=[0-9])", "", colnames(cpm_table), perl = TRUE)

# =========================
# 2) Sample kolonlarını tanımla
# =========================

gene_info_cols <- c("Geneid", "GeneName", "NewName")

sample_cols <- setdiff(colnames(counts), gene_info_cols)

# Counts ve CPM sample isimleri uyumlu mu kontrol et
if (!all(sample_cols %in% colnames(cpm_table))) {
  missing_cols <- sample_cols[!sample_cols %in% colnames(cpm_table)]
  stop(
    "CPM tablosunda şu sample kolonları bulunamadı: ",
    paste(missing_cols, collapse = ", ")
  )
}

# =========================
# 3) Matrix oluştur
# =========================

count_matrix <- counts[, sample_cols]
count_matrix <- as.matrix(count_matrix)
rownames(count_matrix) <- counts$Geneid

cpm_matrix <- cpm_table[, sample_cols]
cpm_matrix <- as.matrix(cpm_matrix)
rownames(cpm_matrix) <- cpm_table$Geneid

# Numeric kontrol
count_matrix <- apply(count_matrix, 2, as.numeric)
rownames(count_matrix) <- counts$Geneid

cpm_matrix <- apply(cpm_matrix, 2, as.numeric)
rownames(cpm_matrix) <- cpm_table$Geneid

# =========================
# 4) Metadata oluştur
# =========================

metadata <- data.frame(
  sample_id = sample_cols,
  condition = ifelse(grepl("SC", sample_cols), "SC", "Ssh1"),
  stringsAsFactors = FALSE
)

metadata$condition <- factor(metadata$condition, levels = c("SC", "Ssh1"))

# Kontrol
stopifnot(all(colnames(count_matrix) == metadata$sample_id))
stopifnot(all(colnames(cpm_matrix) == metadata$sample_id))

# =========================
# 5) Library size hesapla
# =========================

library_sizes <- colSums(count_matrix)

library_size_table <- data.frame(
  sample_id = names(library_sizes),
  library_size = as.numeric(library_sizes),
  condition = metadata$condition
)

write.csv(
  library_size_table,
  "results/RNAseq/01_preprocessing/library_sizes.csv",
  row.names = FALSE
)

# =========================
# 6) CPM'e göre gen filtreleme
# =========================
# Kural:
# Bir gen, herhangi bir gruptaki örneklerin en az %50'sinde
# CPM >= 1 ise tutulur.

filter_by_group_cpm <- function(x, groups, threshold = 1, percentage = 0.5) {
  
  for (group in unique(groups)) {
    
    group_index <- which(groups == group)
    n_required <- ceiling(length(group_index) * percentage)
    
    if (sum(x[group_index] >= threshold, na.rm = TRUE) >= n_required) {
      return(TRUE)
    }
  }
  
  return(FALSE)
}

keep_genes <- apply(
  cpm_matrix,
  1,
  filter_by_group_cpm,
  groups = metadata$condition,
  threshold = 1,
  percentage = 0.5
)

# =========================
# 7) Filtrelenmiş tablolar
# =========================

kept_gene_ids <- rownames(cpm_matrix)[keep_genes]

count_filtered <- counts[counts$Geneid %in% kept_gene_ids, ]
cpm_filtered <- cpm_table[cpm_table$Geneid %in% kept_gene_ids, ]

# Aynı gen sırasına getir
count_filtered <- count_filtered[match(kept_gene_ids, count_filtered$Geneid), ]
cpm_filtered <- cpm_filtered[match(kept_gene_ids, cpm_filtered$Geneid), ]

# =========================
# 8) Sonuçları kaydet
# =========================

write.csv(
  count_filtered,
  "results/RNAseq/01_preprocessing/count_filtered.csv",
  row.names = FALSE
)

write.csv(
  cpm_filtered,
  "results/RNAseq/01_preprocessing/cpm_filtered.csv",
  row.names = FALSE
)

write.csv(
  metadata,
  "results/RNAseq/01_preprocessing/metadata.csv",
  row.names = FALSE
)

filter_summary <- data.frame(
  input_gene_number = nrow(counts),
  filtered_gene_number = nrow(count_filtered),
  removed_gene_number = nrow(counts) - nrow(count_filtered)
)

write.csv(
  filter_summary,
  "results/RNAseq/01_preprocessing/filter_summary.csv",
  row.names = FALSE
)

# =========================
# 9) Özet yazdır
# =========================

cat("Input gene number:", nrow(counts), "\n")
cat("Filtered gene number:", nrow(count_filtered), "\n")
cat("Removed gene number:", nrow(counts) - nrow(count_filtered), "\n")

cat("\nLibrary sizes:\n")
print(library_size_table)

cat("\nPreprocessing completed successfully.\n")
