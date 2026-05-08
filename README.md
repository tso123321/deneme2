# OMICS Course - HPC Linux Conda Training

Bu repository, biyoinformatik analizlerde kullanılan **HPC (High Performance Computing)** sistemlerinin temel kullanımını öğretmek amacıyla hazırlanmıştır. Amaç, katılımcıların teoriyi okumaktan ziyade **komutları çalıştırarak öğrenmesini sağlamaktır**.

---

## Eğitim Akışı

Linux → Dosya yönetimi → Conda → Tool kurulumu → Script → HPC → SLURM

---


Bu repository’yi çalıştırmadan önce gerekli conda environment oluşturulmalıdır.

Codespaces ve bazı HPC ortamlarında conda otomatik aktif gelmez. Bu nedenle önce conda’yı shell’e bağlayın:

source /opt/conda/etc/profile.d/conda.sh

Ardından environment oluşturun ve aktive edin:

conda create -n hpc_course -y  
conda activate hpc_course  

Gerekli biyoinformatik araçları kurun:

conda install -c conda-forge -c bioconda seqkit fastqc samtools -y  

---

### Kurulum Kontrolü

which seqkit  
which fastqc  
which samtools  

Çıktı şu formatta olmalıdır:

/opt/conda/envs/hpc_course/...

---

### ÖNEMLİ

- `bioconda` tek başına yeterli değildir, `conda-forge` ile birlikte kullanılmalıdır  
- Conda environment aktif değilse tool’lar çalışmaz  
- Script içinde conda environment tekrar aktive edilmelidir:

source /opt/conda/etc/profile.d/conda.sh  
conda activate hpc_course  
---

## Eğitimi Çalıştırma

Tüm eğitim tek komutla çalıştırılabilir:

bash scripts/run_all_local.sh  

Bu script aşağıdaki adımları sırasıyla çalıştırır:

1. Linux temel komutları  
2. Conda environment kontrolü  
3. Tool çalıştırma (seqkit)  
4. Loop ile çoklu dosya işleme  

---

## Repository Yapısı

.  
├── data/        → örnek veri (FASTA)  
├── scripts/     → eğitim scriptleri  
├── results/     → analiz çıktıları  
└── slurm/       → HPC job scriptleri  

---

# Linux Nedir?

Linux, HPC sistemlerinin temel işletim sistemidir.  
Tüm işlemler terminal üzerinden yapılır.

---

## Temel Linux Komutları

pwd → bulunduğun dizini gösterir  
ls → dosyaları listeler  
ls -lh → boyutlarıyla listeler  

cd klasor → klasöre gir  
cd .. → bir üst klasör  
cd ~ → home dizini  

---

## Dosya İşlemleri

mkdir klasor → klasör oluştur  
mkdir -p a/b/c → iç içe klasör  

touch dosya.txt → boş dosya oluştur  

cp a.txt b.txt → kopyala  
mv a.txt yeni.txt → taşı/yeniden adlandır  

rm dosya.txt → dosya sil  
rm -r klasor → klasör sil  

---

## Dosya İçeriği

cat dosya.txt → tamamını göster  
head dosya.txt → ilk satırlar  
tail dosya.txt → son satırlar  
wc -l dosya.txt → satır sayısı  

---

## Komut Zinciri (Pipe)

cat dosya.txt | wc -l  

Bir komutun çıktısını diğerine verir.

---

# Conda Nedir?

Conda, farklı yazılım ortamlarını izole şekilde yönetmeyi sağlar.

### Neden Conda?

- Her proje için ayrı environment  
- Versiyon çakışması olmaz  
- Reproducibility sağlar  

---

## Conda Komutları

conda env list → ortamları listele  
conda create -n env → yeni ortam  
conda activate env → aktif et  
conda deactivate → çık  

---

## Paket Kurulumu

conda install paket  

Biyoinformatik için:

conda install -c bioconda seqkit  

---

# Tool Kullanımı

Örnek tool: seqkit

seqkit stats data/example.fasta  

Bu komut:
- kaç sequence var  
- toplam uzunluk  
- ortalama uzunluk  
gibi bilgileri verir.

---

# Bash Script Nedir?

Tek tek komut yazmak yerine otomatik çalıştırma sağlar.

Örnek:

#!/bin/bash  
echo "Hello"  

Çalıştırma:

bash script.sh  

---

# Loop (Döngü)

Aynı işlemi birden fazla dosyada çalıştırmak için:

for file in data/*.fasta  
do  
  echo $file  
done  

---

# HPC Nedir?

HPC (High Performance Computing), büyük hesaplamaların güçlü bilgisayar kümelerinde yapılmasını sağlar.

---

## HPC Mimarisi

Login node:
- sisteme giriş yapılır  
- küçük işler yapılır  

Compute node:
- gerçek hesaplama yapılır  
- job burada çalışır  

---

# SLURM Nedir?

SLURM, HPC sistemlerinde job yönetimi yapan scheduler’dır.

Kullanıcı doğrudan compute node’da çalışmaz → job gönderir.

---

## Çalıştırma Farkı

bash script.sh → lokal çalıştırma  
sbatch script.sh → HPC çalıştırma  

---

## SLURM Script

#!/bin/bash  
#SBATCH --job-name=test  
#SBATCH --output=results/slurm_%j.out  
#SBATCH --error=results/slurm_%j.err  
#SBATCH --time=00:05:00  
#SBATCH --cpus-per-task=1  
#SBATCH --mem=1G  

Bu parametreler:
- CPU  
- RAM  
- süre  
gibi kaynakları belirler.

---

## SLURM Komutları

sbatch script.sh → job gönder  
squeue → jobları listele  
squeue -u $USER → kendi jobların  
scancel JOB_ID → job iptal  

---

## Çıktı Dosyaları

results/slurm_JOBID.out  
results/slurm_JOBID.err  
