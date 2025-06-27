![image](https://github.com/user-attachments/assets/84f2237d-9614-4473-8054-01d7ba0782d6)

# MICROviso

## Mục đích

Bộ hàm hỗ trợ tiền xử lý và trực quan hóa dữ liệu vi sinh vật từ đối tượng `phyloseq`, bao gồm:

- Import dữ liệu
- Tạo cây phát sinh chủng loài (phylogenetic tree)
- Tính toán đa dạng alpha
- Biểu đồ thành phần vi sinh vật theo cấp độ phân loại (phylum, genus)
- Kiểm tra sự biến đổi theo nhóm
- Tính toán thống kê

---

## Installation

Dependencies:

```r
# CRAN
install.packages("tidyverse")
install.packages("patchwork")
install.packages("ggplot2")
install.packages("RColorBrewer")
install.packages("devtools")
install.packages("ggpubr")

# Bioconductor
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("phyloseq")
BiocManager::install("lefser")
BiocManager::install("msa")
BiocManager::install("phangorn")
BiocManager::install("SummarizedExperiment")

install.packages("tidyverse")
install.packages("patchwork")
install.packages("ggplot2")
install.packages("BiocManager") 
BiocManager::install("DESeq2")
BiocManager::install("phyloseq")
BiocManager::install("lefser")
install.packages("RColorBrewer")
install.packages("devtools")
install.packages("ggpubr")
install.packages("SummarizedExperiment")
```

MICROviso:

```bash
devtools::install_github("anegin24/MICROviso")
```

Kích hoạt các thư viện:

```r
library(phyloseq)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(DESeq2)
library(RColorBrewer)
library(devtools)
library(ggpubr)
library(MICROviso)
library(lefser)
library(SummarizedExperiment)
```

---

## 1. Import dữ liệu từ `phyloseq`

Import đối tượng phyloseq _**ps**_ vào R và xuất ra các bảng **_table, taxonomy, metadata_**:

```r
import_phyloseq("v3-v4.phyloseq")
```

---

## 2. Import metadata

Import metadata từ file bên ngoài (hoặc để cập nhật metadata):

```r
metadata <- import_metadata("sample-metadata.tsv")
```

---

## 3. Dựng cây phát sinh chủng loài (Phylogenetic Tree)

- Sử dụng hàm `make_tree()` để căn chỉnh trình tự từ `refseq(ps)`, xây dựng cây NJ, và (tuỳ chọn) tối ưu bằng Maximum Likelihood (GTR).
- Hỗ trợ vẽ cây trực quan với `ggtree`.

```r
ps <- make_tree(ps, msa_method = "ClustalOmega", optimize_ml = TRUE, plot_tree = TRUE)
```

Tham số:
- `msa_method`: "ClustalW", "Muscle", hoặc "ClustalOmega"
- `optimize_ml`: TRUE/FALSE để bật hoặc tắt tối ưu hóa bằng Maximum Likelihood
- `plot_tree`: TRUE/FALSE để vẽ cây

---

## 4. Tính đa dạng alpha (alpha diversity)

```r
alpha <- cal_alpha(ps)
```

Tạo bảng thống kê:

```r
res3 <- cal_alpha_stat(alpha, metadata, group_col = "treatment", strata = "timeline")
```

---

## 5. Vẽ biểu đồ alpha diversity

```r
alphaplot <- plot_alpha(alpha = alpha, metadata = metadata, x = "treatment", facet = "timeline")
```

---

## 6. Vẽ biểu đồ beta diversity

```r
betadata <- plot_beta(ps, color = "treatment", facet = "timeline", distance_method = "bray", method = "PCoA")
```

---

## 7. Vẽ biểu đồ thành phần Phylum

```r
plot_phylum(ps, group_vars = c("treatment", "timeline"), facet = "timeline", x_var = "treatment")
```

---

## 8. Vẽ biểu đồ thành phần Genus

```r
plot_genus(ps, group_vars = c("Sample"), top = 20, x_var = "Sample")

plot_genus(ps, group_vars = c("treatment", "timeline"), top = 20, facet = "timeline", x_var = "treatment")
```

---

## 9. Phân tích sự khác biệt vi sinh vật với LEfSe

```r
se <- phyloseq_to_se(ps)
out <- run_lefse_pairwise(se, classCol = "class", groups = c("ABX 0.2X_Week 6", "Control_Week 6"))
```

---

## 10. Phân tích khác biệt bằng DESeq2

### Cross-sectional study

```r
DEseq2_cross(physeq, group = NULL, comparison = NULL, padj_cutoff = 0.05)
```

### Detect Genus with Group × Time Interaction (LRT)

```r
DEseq2_global(physeq, taxrank = "Genus", group = NULL, time_var = NULL, alpha = 0.05)
```

### Pairwise comparison

```r
ps@sam_data$treatment <- gsub(" ", "_", ps@sam_data$treatment)
sample_data(ps)$timeline <- gsub(" ", "_", sample_data(ps)$timeline)

DEseq2_pairwise(ps, "Week 0", group = "treatment", time_var = "timeline", 
                comparison = c("ABX 0.2X", "ABX 0.5X"), padj_cutoff = 0.05)
```
