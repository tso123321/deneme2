# 🧬 RNA-seq Workflow – Script Açıklamaları

## 00_R_basics.R  
**Amaç:**  
R programlama diline giriş yapmak ve temel veri işlemlerini öğretmek.

**Ne yapılıyor:**  
- Değişken tanımlama  
- Vektör ve data.frame oluşturma  
- Temel fonksiyonlar (`mean`, `sum`, vb.)  
- Veri okuma ve yazma  
- Basit veri manipülasyonu  

Bu script, RNA-seq analizine başlamadan önce gerekli R altyapısını sağlar.

---

## 01_preprocessing_filtering.R  
**Amaç:**  
Düşük ifade edilen genleri filtreleyerek analizdeki gürültüyü azaltmak.

**Ne yapılıyor:**  
- Count ve CPM tabloları okunur  
- CPM (Counts Per Million) hesaplanır  
- CPM ≥ 1 filtresi uygulanır  
- Düşük ifade edilen genler çıkarılır  

Bu adım:
- false positive sonuçları azaltır  
- istatistiksel gücü artırır  

---

## 02_annotation.R  
**Amaç:**  
Ensembl gene ID’lerini biyolojik olarak anlamlı gene isimleri ile eşleştirmek.

**Ne yapılıyor:**  
- GENCODE GTF dosyası `rtracklayer` ile okunur  
- `gene_id`, `gene_name`, `gene_type` bilgileri çıkarılır  
- Gen koordinat bilgileri eklenir (chr, start, end, strand)  
- Count tablosu ile birleştirilir  

### 📥 GENCODE Veri İndirme

Bu adımı çalıştırmadan önce GENCODE annotation dosyasını indirmeniz gerekmektedir.

İndirme linki:  
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.basic.annotation.gtf.gz  

gunzip gencode.v46.primary_assembly.basic.annotation.gtf.gz

mv gencode.v46.primary_assembly.basic.annotation.gtf data/

---

## 03_PCA.R  
**Amaç:**  
Örnekler arasındaki genel varyasyonu ve grup ayrımını incelemek.

**Ne yapılıyor:**  
- Normalized veri kullanılır  
- Principal Component Analysis (PCA) uygulanır  
- Örneklerin 2D uzaydaki dağılımı görselleştirilir  

PCA ile:
- batch effect  
- outlier örnekler  
- grup ayrımı  

gözlemlenir.

---

## 04_sample_correlation_heatmap.R  
**Amaç:**  
Örnekler arasındaki benzerliği ölçmek.

**Ne yapılıyor:**  
- Korelasyon matrisi hesaplanır  
- Heatmap oluşturulur  

Bu analiz:
- replikaların tutarlılığını gösterir  
- beklenmeyen örnekleri ortaya çıkarır  

---

## 05_differential_expression_edgeR.R  
**Amaç:**  
İki grup (SC vs Ssh1) arasında diferansiyel ifade edilen genleri belirlemek.

**Ne yapılıyor:**  
- Count matrix oluşturulur  
- edgeR ile normalizasyon (TMM) yapılır  
- Dispersion hesaplanır  
- exactTest uygulanır  
- logFC ve FDR hesaplanır  
- hangi genlerin up/down regüle olduğunu belirler  

---

## 06_volcano_MA_plots.R  
**Amaç:**  
Diferansiyel ekspresyon sonuçlarını görselleştirmek.

**Ne yapılıyor:**  
- Volcano plot oluşturulur  
- MA plot oluşturulur  

Bu grafikler:
- anlamlı genleri hızlıca görmeyi sağlar  
- effect size vs significance ilişkisini gösterir  

---

## 07_significant_gene_heatmap.R  
**Amaç:**  
Anlamlı genlerin ekspresyon patternlerini incelemek.

**Ne yapılıyor:**  
- FDR < 0.05 genler seçilir  
- Expression matrisi normalize edilir  
- Heatmap oluşturulur  

Bu analiz:
- genlerin grup bazlı clustering’ini gösterir  

---

## 08_GO_enrichment.R  
**Amaç:**  
Anlamlı genlerin hangi biyolojik süreçlerde yer aldığını belirlemek.

**Ne yapılıyor:**  
- Gene Ontology (GO) enrichment analizi yapılır  
- Biyolojik süreçler (BP), hücresel bileşenler (CC) ve moleküler fonksiyonlar (MF) incelenir  
- biyolojik yorum üretir  

---

## 09_KEGG_enrichment.R  
**Amaç:**  
Genlerin hangi biyolojik yolaklarda (pathway) yer aldığını belirlemek.

**Ne yapılıyor:**  
- KEGG pathway enrichment analizi yapılır  
- Anlamlı yolaklar belirlenir  

Bu analiz:
- hastalık mekanizmaları  
- sinyal yolakları  
hakkında bilgi verir  

---

# Genel Özet
→ Filtering  
→ Annotation  
→ Differential Expression  
→ Visualization  
→ Functional Interpretation  


