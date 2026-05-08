############################################################
# 02_annotation.R
# RNA-seq Analysis Course
# Step 2: Gene annotation using GENCODE GTF with rtracklayer
############################################################

############################################################
# 1. Packages
############################################################

library(dplyr)
library(rtracklayer)
library(GenomicRanges)

############################################################
# 2. Output folder
############################################################

dir.create("results/RNAseq/02_annotation", recursive = TRUE, showWarnings = FALSE)

############################################################
# 3. Input files
############################################################

count_file <- "results/RNAseq/01_preprocessing/count_filtered.csv"
gtf_file <- "data/gencode.v46.primary_assembly.basic.annotation.gtf"

############################################################
# 4. Read filtered count table
############################################################

count_filtered <- read.csv(
  count_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"Geneid" %in% colnames(count_filtered)) {
  stop("count_filtered tablosunda Geneid kolonu yok.")
}

############################################################
# 5. Clean count Gene IDs
############################################################

count_filtered <- count_filtered %>%
  mutate(
    Geneid_clean = gsub("\\..*", "", Geneid)
  )

gene_info_cols <- c("Geneid", "GeneName", "NewName", "Geneid_clean")
sample_cols <- setdiff(colnames(count_filtered), gene_info_cols)

############################################################
# 6. Import GENCODE GTF with rtracklayer
############################################################

gtf <- rtracklayer::import(gtf_file)

# GTF içeriğini kontrol
cat("Imported GTF records:", length(gtf), "\n")

############################################################
# 7. Keep only gene-level records
############################################################

gtf_gene <- gtf[gtf$type == "gene"]

cat("Gene records in GTF:", length(gtf_gene), "\n")

############################################################
# 8. Build annotation table
############################################################

annotation <- data.frame(
  Geneid = mcols(gtf_gene)$gene_id,
  GeneName = mcols(gtf_gene)$gene_name,
  GeneType = mcols(gtf_gene)$gene_type,
  chr = as.character(seqnames(gtf_gene)),
  start = start(gtf_gene),
  end = end(gtf_gene),
  strand = as.character(strand(gtf_gene)),
  stringsAsFactors = FALSE
)

annotation <- annotation %>%
  mutate(
    Geneid_clean = gsub("\\..*", "", Geneid),
    NewName = paste(Geneid_clean, GeneName, sep = "_")
  ) %>%
  select(
    Geneid_clean,
    GeneName,
    NewName,
    GeneType,
    chr,
    start,
    end,
    strand
  ) %>%
  distinct(Geneid_clean, .keep_all = TRUE)

############################################################
# 9. Join count table with GTF annotation
############################################################

count_annotated <- count_filtered %>%
  left_join(annotation, by = "Geneid_clean", suffix = c("_count", "_gtf"))

############################################################
# 10. Resolve annotation fields
############################################################

count_annotated <- count_annotated %>%
  mutate(
    GeneName = case_when(
      !is.na(GeneName_gtf) & GeneName_gtf != "" ~ GeneName_gtf,
      !is.na(GeneName_count) & GeneName_count != "" ~ GeneName_count,
      TRUE ~ Geneid_clean
    ),
    NewName = paste(Geneid_clean, GeneName, sep = "_")
  )

############################################################
# 11. Unannotated genes
############################################################

unannotated_genes <- count_annotated %>%
  filter(is.na(GeneName_gtf))

write.csv(
  unannotated_genes,
  "results/RNAseq/02_annotation/unannotated_genes.csv",
  row.names = FALSE
)

############################################################
# 12. Duplicate gene control
############################################################

duplicated_genes <- count_annotated$Geneid_clean[
  duplicated(count_annotated$Geneid_clean)
]

count_annotated <- count_annotated %>%
  distinct(Geneid_clean, .keep_all = TRUE)

############################################################
# 13. Final column order
############################################################

count_annotated <- count_annotated %>%
  select(
    Geneid = Geneid_clean,
    GeneName,
    NewName,
    GeneType,
    chr,
    start,
    end,
    strand,
    all_of(sample_cols)
  )

############################################################
# 14. Count matrix
############################################################

count_matrix <- count_annotated[, sample_cols]
count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "numeric"
rownames(count_matrix) <- count_annotated$Geneid

############################################################
# 15. Save outputs
############################################################

write.csv(
  count_annotated,
  "results/RNAseq/02_annotation/count_annotated.csv",
  row.names = FALSE
)

write.table(
  count_matrix,
  "results/RNAseq/02_annotation/count_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

annotation_summary <- data.frame(
  input_gene_number = nrow(count_filtered),
  gtf_gene_number = nrow(annotation),
  annotated_gene_number = sum(!is.na(count_annotated$GeneName)),
  unannotated_gene_number = nrow(unannotated_genes),
  duplicate_genes_removed = length(duplicated_genes),
  final_gene_number = nrow(count_annotated),
  sample_number = length(sample_cols)
)

write.csv(
  annotation_summary,
  "results/RNAseq/02_annotation/annotation_summary.csv",
  row.names = FALSE
)

############################################################
# 16. Summary
############################################################

cat("Input gene number:", nrow(count_filtered), "\n")
cat("GENCODE gene annotation number:", nrow(annotation), "\n")
cat("Final annotated gene number:", nrow(count_annotated), "\n")
cat("Unannotated genes:", nrow(unannotated_genes), "\n")
cat("Duplicate genes removed:", length(duplicated_genes), "\n")
cat("Sample number:", length(sample_cols), "\n")

