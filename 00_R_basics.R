############################################################
# 01_R_basics.R
# OMICS Course - R Temel Kodlama Girişi
############################################################

############################################################
# 1. R nedir?
############################################################

# R; veri analizi, istatistik ve görselleştirme için kullanılan
# bir programlama dilidir.
#
# Biyoinformatikte R özellikle şu işler için kullanılır:
# - count matrix okuma
# - metadata düzenleme
# - normalizasyon
# - PCA
# - heatmap
# - diferansiyel ekspresyon analizi
# - GO / KEGG enrichment analizi
# - grafik üretimi
#

############################################################
# 2. Basit hesaplamalar
############################################################

# R bir hesap makinesi gibi kullanılabilir.
# Aşağıdaki komutlar doğrudan console'a sonuç döndürür.

2 + 2      # toplama
10 - 3     # çıkarma
5 * 4      # çarpma
20 / 5     # bölme
2^3        # üs alma

# R komutları yukarıdan aşağıya çalışır.
# Script içindeki her satır sırayla değerlendirilir.

############################################################
# 3. Değişken oluşturma
############################################################

# Değişken, bir değeri daha sonra tekrar kullanabilmek için
# bir isim altında saklamamızı sağlar.
#
# R'da değişken oluşturmak için genellikle <- kullanılır.
# = de kullanılabilir; ancak R kodlarında <- daha yaygındır.
#
# Sağdaki değer soldaki değişken adına atanır.

x <- 10
y <- 5

# Değişken adını yazarsak içindeki değeri görürüz.

x
y

# Değişkenlerle işlem yapılabilir.

x + y
x * y

# Metin değerleri character veri tipidir ve tırnak içinde yazılır.

sample_name <- "control_1"
gene_name <- "geneA"

sample_name
gene_name

############################################################
# 4. Temel veri tipleri
############################################################

# R'da her objenin bir veri tipi vardır.
# class() fonksiyonu objenin tipini gösterir.

# Numeric: sayısal veri
# Count, yaş, ölçüm değeri gibi sayılar numeric olabilir.

count_value <- 25
count_value
class(count_value)

# Character: metin veri
# Gen adı, örnek adı, grup adı gibi metinler character olabilir.

gene_id <- "geneA"
gene_id
class(gene_id)

# Logical: TRUE / FALSE
# Koşul sonuçları veya evet/hayır bilgileri logical olabilir.

is_expressed <- TRUE
is_expressed
class(is_expressed)

############################################################
# 5. Vector
############################################################

# Vector, aynı tipte birden fazla değeri tutan veri yapısıdır.
# R'daki en temel veri yapılarından biridir.
#
# c() fonksiyonu combine anlamına gelir.
# Birden fazla değeri bir araya getirerek vector oluşturur.

counts <- c(10, 25, 15, 3, 20)
counts

genes <- c("geneA", "geneB", "geneC", "geneD", "geneE")
genes

# R'da indeksleme 1'den başlar.
# Yani ilk eleman [1] ile seçilir.

counts[1]   # counts vectorünün ilk elemanı
counts[2]   # ikinci elemanı
genes[1]    # genes vectorünün ilk elemanı
genes[3]    # üçüncü elemanı

# Vector üzerinde temel özet istatistikler hesaplanabilir.

sum(counts)       # toplam
mean(counts)      # ortalama
median(counts)    # medyan
min(counts)       # minimum
max(counts)       # maksimum
length(counts)    # eleman sayısı

############################################################
# 6. Data frame
############################################################

# Data frame, tablo yapısıdır.
# Excel tablosu gibi düşünebilirsiniz.
#
# Her kolon bir değişkeni temsil eder.
# Her satır bir gözlemi temsil eder.
#
# RNA-seq analizinde:
# - count matrix
# - metadata
# - diferansiyel ekspresyon sonuçları
# genellikle data frame olarak tutulur.

gene_counts <- data.frame(
  gene_id = c("geneA", "geneB", "geneC", "geneD", "geneE"),
  control_1 = c(10, 25, 15, 3, 20),
  control_2 = c(11, 22, 14, 4, 21),
  treated_1 = c(30, 8, 15, 12, 20),
  treated_2 = c(28, 7, 16, 11, 19)
)

gene_counts

# Data frame hakkında bilgi alma:

str(gene_counts)       # kolonların yapısını ve veri tiplerini gösterir
dim(gene_counts)       # satır ve kolon sayısını verir
nrow(gene_counts)      # satır sayısı
ncol(gene_counts)      # kolon sayısı
colnames(gene_counts)  # kolon isimleri
rownames(gene_counts)  # satır isimleri

# head() veri setinin ilk satırlarını gösterir.
# Büyük dosyalarda veriye hızlıca bakmak için kullanılır.

head(gene_counts)

############################################################
# 7. Data frame içinden seçim yapmak
############################################################

# Bir data frame içinden kolon, satır veya belirli hücreler seçilebilir.

# $ işareti kolon seçmek için kullanılır.
# Aşağıdaki kod gene_id kolonunu getirir.

gene_counts$gene_id

# control_1 kolonunu seçelim.

gene_counts$control_1

# Köşeli parantez yapısı:
# data[satır, kolon]
#
# gene_counts[1, ]  -> 1. satır, tüm kolonlar
# gene_counts[, 1]  -> tüm satırlar, 1. kolon

gene_counts[1, ]
gene_counts[2, ]

# Satır ve kolon birlikte seçilebilir.
# Örneğin 1. satırdaki control_1 değerini alalım.

gene_counts[1, "control_1"]

# 3. satırdaki treated_1 değerini alalım.

gene_counts[3, "treated_1"]

# Birden fazla kolon seçmek için c() kullanılır.

gene_counts[, c("gene_id", "control_1", "treated_1")]

# Belirli koşula göre satır seçmek:
# Burada control_1 değeri 10'dan büyük olan genler seçilir.

gene_counts[gene_counts$control_1 > 10, ]

# Bu mantık RNA-seq analizlerinde çok kullanılır.
# Örneğin FDR < 0.05 olan genleri seçmek gibi.

############################################################
# 8. Factor
############################################################

# Factor, kategorik değişkenleri temsil eder.
# Örneğin grup bilgisi: control / treated
#
# R'da istatistiksel modeller kategorik değişkenleri factor olarak
# daha doğru yorumlar.

condition <- c("control", "control", "treated", "treated")
condition

# Character vectorünü factor'a çevirelim.

condition_factor <- factor(condition)
condition_factor

class(condition_factor)
levels(condition_factor)

# levels(), factor içindeki kategorileri gösterir.
# Burada iki kategori vardır:
# control ve treated.
#
# RNA-seq analizlerinde condition bilgisinin factor olması önemlidir.
# Çünkü DESeq2 veya edgeR gibi paketler grupları bu şekilde tanır.

############################################################
# 9. Metadata oluşturma
############################################################

# Metadata, örneklere ait bilgileri tutan tablodur.
#
# RNA-seq analizlerinde count matrix ayrı, metadata ayrı olabilir.
# Count matrix genlerin sayımlarını içerir.
# Metadata ise örneklerin grup bilgisini içerir.

metadata <- data.frame(
  sample_id = c("control_1", "control_2", "treated_1", "treated_2"),
  condition = factor(c("control", "control", "treated", "treated"))
)

metadata
str(metadata)

# Burada sample_id örnek adıdır.
# condition ise örneklerin hangi deney grubuna ait olduğunu gösterir.

############################################################
# 10. Matrix
############################################################

# Matrix, sadece aynı tipte veri içeren iki boyutlu yapıdır.
# Data frame farklı tipte kolonlar içerebilir.
# Matrix ise genellikle tamamen numeric olur.
#
# RNA-seq count matrix genellikle numeric matrix olarak analiz edilir.

count_matrix <- gene_counts[, c("control_1", "control_2", "treated_1", "treated_2")]

count_matrix

# Data frame'i matrix'e çevirelim.

count_matrix <- as.matrix(count_matrix)

# Gen isimlerini satır adı olarak ekleyelim.
# Böylece her satır bir geni temsil eder.

rownames(count_matrix) <- gene_counts$gene_id

count_matrix
class(count_matrix)

############################################################
# 11. Satır ve kolon işlemleri
############################################################

# Count matrix'te:
# satırlar genleri,
# kolonlar örnekleri temsil eder.
#
# rowSums(): her satır için toplam hesaplar.
# Burada her genin tüm örneklerdeki toplam count değerini verir.

rowSums(count_matrix)

# colSums(): her kolon için toplam hesaplar.
# Burada her örneğin toplam read count değerini verir.
# RNA-seq'te buna kaba olarak library size mantığı diyebiliriz.

colSums(count_matrix)

# rowMeans(): her gen için ortalama count hesaplar.

rowMeans(count_matrix)

# colMeans(): her örnek için ortalama count hesaplar.

colMeans(count_matrix)

############################################################
# 12. Grup ortalamaları
############################################################

# Burada iki control ve iki treated örneğimiz var.
# Her gen için control grubunun ortalamasını hesaplayalım.

control_mean <- rowMeans(count_matrix[, c("control_1", "control_2")])

# Her gen için treated grubunun ortalamasını hesaplayalım.

treated_mean <- rowMeans(count_matrix[, c("treated_1", "treated_2")])

control_mean
treated_mean

# Bu işlem RNA-seq analizinde grupları karşılaştırmanın temel mantığını gösterir.
# Gerçek analizlerde sadece ortalamaya bakılmaz;
# varyasyon ve istatistiksel modelleme de gerekir.

############################################################
# 13. Basit log2 fold change
############################################################

# log2 fold change, iki grup arasındaki değişimin yönünü ve büyüklüğünü gösterir.
#
# Formül:
# log2FC = log2(treated_mean / control_mean)
#
# Burada +1 ekliyoruz çünkü bazı genlerde count 0 olabilir.
# 0'a bölme veya log2(0) problemi olmaması için +1 kullanılır.

log2FC <- log2((treated_mean + 1) / (control_mean + 1))

log2FC

# Yorum:
# log2FC > 0  → treated grubunda daha yüksek
# log2FC < 0  → treated grubunda daha düşük
# log2FC = 0  → iki grup benzer

result_table <- data.frame(
  gene_id = rownames(count_matrix),
  control_mean = control_mean,
  treated_mean = treated_mean,
  log2FC = log2FC
)

result_table

############################################################
# 14. Koşullu filtreleme
############################################################

# Data frame içinden belirli koşullara uyan satırları seçebiliriz.
# Örneğin log2FC > 1 olan genler treated grubunda artmış kabul edilebilir.

up_genes <- result_table[result_table$log2FC > 1, ]
up_genes

# log2FC < -1 olan genler treated grubunda azalmış kabul edilebilir.

down_genes <- result_table[result_table$log2FC < -1, ]
down_genes

# Gerçek RNA-seq analizlerinde sadece log2FC'ye göre karar verilmez.
# FDR / adjusted p-value gibi istatistiksel değerler de kullanılır.

############################################################
# 15. Paket nedir?
############################################################

# R paketleri, R'ın temel fonksiyonlarına ek olarak yeni fonksiyonlar sağlar.
#
# Örneğin:
# ggplot2 → grafik çizmek için
# dplyr   → tablo düzenlemek için
# edgeR   → RNA-seq diferansiyel ekspresyon analizi için
#
# Paket yüklemek için:
# install.packages("ggplot2")
#
# Paket çağırmak için:
# library(ggplot2)
#
# Bir paketi yüklemek ve çağırmak farklı şeylerdir:
# install.packages() paketi bilgisayara indirir.
# library() paketi aktif R oturumunda kullanılabilir hale getirir.

############################################################
# 16. Paket kontrolü ve yükleme
############################################################

# requireNamespace(), bir paketin kurulu olup olmadığını kontrol eder.
# quietly = TRUE, kontrol sırasında gereksiz mesajları azaltır.
#
# Eğer ggplot2 kurulu değilse install.packages() ile kurulacaktır.

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

# library() paketi aktif hale getirir.
# Bundan sonra ggplot() gibi fonksiyonları kullanabiliriz.

library(ggplot2)

############################################################
# 17. Basit grafik: base R
############################################################

# R'ın kendi içinde gelen temel grafik sistemi vardır.
# Buna base R plotting denir.
#
# Burada geneA'nın her örnekteki count değerlerini barplot olarak çiziyoruz.

barplot(
  count_matrix["geneA", ],       # çizilecek değerler: geneA satırı
  main = "geneA read counts",    # grafik başlığı
  ylab = "Count",                # y ekseni adı
  las = 2                        # x ekseni yazılarını dik gösterir
)

# count_matrix["geneA", ] ifadesi:
# geneA satırını ve tüm kolonları seçer.
#
# Bu grafik tek bir genin örnekler arasındaki count değişimini gösterir.

############################################################
# 18. ggplot2 ile grafik
############################################################

# ggplot2, daha esnek ve düzenli grafikler üretmek için kullanılır.
#
# ggplot mantığı:
# 1. Veri tablosu verilir.
# 2. aes() ile hangi değişkenin hangi eksene gideceği tanımlanır.
# 3. geom_*() ile grafik tipi seçilir.

# Önce geneA için küçük bir data frame oluşturalım.
# Çünkü ggplot2 genellikle long-format data frame ile rahat çalışır.

geneA_df <- data.frame(
  sample_id = colnames(count_matrix),
  count = as.numeric(count_matrix["geneA", ]),
  condition = metadata$condition
)

geneA_df

# ggplot(geneA_df, ...)
# Grafik için geneA_df tablosunu kullanır.
#
# aes(x = sample_id, y = count, fill = condition)
# x eksenine sample_id,
# y eksenine count,
# bar renklerine condition bilgisini koyar.
#
# geom_col()
# Kolon/bar grafiği çizer.
#
# theme_classic()
# Daha temiz bir tema kullanır.
#
# labs()
# Grafik başlığı ve eksen isimlerini düzenler.

ggplot(geneA_df, aes(x = sample_id, y = count, fill = condition)) +
  geom_col() +
  theme_classic() +
  labs(
    title = "geneA expression",
    x = "Sample",
    y = "Read count"
  )

############################################################
# 19. Dosya yazma
############################################################

# Analiz sonuçlarını dosya olarak dışarı yazmak önemlidir.
# Böylece sonuçlar daha sonra tekrar kullanılabilir.
#
# dir.create() klasör oluşturur.
# recursive = TRUE, iç içe klasörleri de oluşturur.
# showWarnings = FALSE, klasör zaten varsa uyarı vermemesini sağlar.

dir.create("results/R_basics", recursive = TRUE, showWarnings = FALSE)

# write.table() tabloyu dosyaya yazar.
# sep = "\t" dosyanın tab-separated yani TSV olmasını sağlar.
# quote = FALSE metinlerin tırnak içinde yazılmasını engeller.
# row.names = FALSE satır isimlerini dosyaya yazmaz.

write.table(
  result_table,
  file = "results/R_basics/basic_result_table.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# write.csv() CSV dosyası yazar.
# CSV dosyaları Excel ile kolay açılabilir.

write.csv(
  result_table,
  file = "results/R_basics/basic_result_table.csv",
  row.names = FALSE
)

############################################################
# 20. Dosya okuma
############################################################

# read.table() dışarıdaki bir tabloyu R içine okur.
#
# header = TRUE, ilk satırın kolon isimleri olduğunu belirtir.
# sep = "\t", dosyanın tab ile ayrılmış olduğunu belirtir.
# stringsAsFactors = FALSE, metinlerin otomatik factor'a çevrilmesini engeller.

result_from_file <- read.table(
  "results/R_basics/basic_result_table.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

result_from_file

############################################################
# 21. Fonksiyon yazma
############################################################

# Fonksiyon, tekrar eden işlemleri tek bir isim altında toplar.
# Böylece aynı formülü tekrar tekrar yazmak yerine bir kez tanımlarız.

# function(treated, control)
# Fonksiyonun iki girdisi vardır: treated ve control.
#
# Süslü parantez içindeki satır fonksiyonun ne yapacağını belirtir.

calculate_log2fc <- function(treated, control) {
  log2((treated + 1) / (control + 1))
}

# Fonksiyonu tek değerler üzerinde deneyelim.

calculate_log2fc(treated = 30, control = 10)

# Aynı fonksiyonu result_table içindeki kolonlara uygulayabiliriz.

result_table$log2FC_function <- calculate_log2fc(
  treated = result_table$treated_mean,
  control = result_table$control_mean
)

result_table

############################################################
# 22. apply ailesi
############################################################

# apply fonksiyonu matrix veya array üzerinde satır/kolon bazlı işlem yapar.
#
# Genel yapı:
# apply(X, MARGIN, FUN)
#
# X      → kullanılacak veri
# MARGIN → işlemin satıra mı kolona mı uygulanacağını belirtir
# FUN    → uygulanacak fonksiyon
#
# MARGIN = 1 → satırlar
# MARGIN = 2 → kolonlar

# Her gen için maksimum count:
# Çünkü count_matrix'te satırlar genleri temsil eder.
# apply(count_matrix, 1, max) her satırdaki maksimum değeri bulur.

apply(count_matrix, 1, max)

# Her sample için maksimum count:
# Çünkü count_matrix'te kolonlar örnekleri temsil eder.
# apply(count_matrix, 2, max) her kolondaki maksimum değeri bulur.

apply(count_matrix, 2, max)

# Bu mantık büyük RNA-seq tablolarında çok işe yarar.
# Örneğin her gen için ortalama, maksimum veya varyans hesaplanabilir.

apply(count_matrix, 1, mean)
apply(count_matrix, 1, sd)

############################################################
# 23. Basit if / else
############################################################

# if / else koşula göre farklı işlem yapmak için kullanılır.
#
# Eğer koşul TRUE ise if bloğu çalışır.
# Eğer koşul FALSE ise else bloğu çalışır.

value <- 12

if (value > 10) {
  print("value is greater than 10")
} else {
  print("value is 10 or smaller")
}

# RNA-seq bağlantısı:
# Örneğin log2FC > 1 ise gen up-regulated olarak etiketlenebilir.

result_table$regulation <- ifelse(
  result_table$log2FC > 1,
  "up",
  "not_up"
)

result_table

############################################################
# 24. for loop
############################################################

# for loop, aynı işlemi birden fazla eleman üzerinde tekrar eder.
#
# Genel yapı:
# for (eleman in liste) {
#   yapılacak işlem
# }

# Burada rownames(count_matrix) gen isimlerini verir.
# Loop her gen ismini sırayla gene değişkenine atar.

for (gene in rownames(count_matrix)) {
  print(gene)
}

# Her gen için ortalama count yazdıralım.

for (gene in rownames(count_matrix)) {
  gene_mean <- mean(count_matrix[gene, ])
  print(paste(gene, "mean count:", gene_mean))
}

# paste() metinleri birleştirir.
# Bu nedenle çıktıyı daha okunabilir hale getirir.
