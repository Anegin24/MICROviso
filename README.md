![image](https://github.com/user-attachments/assets/84f2237d-9614-4473-8054-01d7ba0782d6)

## Mục đích

Bộ hàm hỗ trợ tiền xử lý và trực quan hóa dữ liệu vi sinh vật từ đối tượng `phyloseq`, bao gồm:

- Import dữ liệu
- Tính toán đa dạng alpha
- Biểu đồ thành phần vi sinh vật theo cấp độ phân loại (phylum, genus)
- Kiểm tra sự biến đổi theo nhóm
- Tính toán thống kê
---
## Installation

Dependencies:

```bash
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

```bash
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

## 1. Import dữ liệu từ `phyloseq`

Chức năng này sẽ import đối tượng phyloseq _**ps**_ vào R và xuất thẳng ra các bảng **_table, taxonomy, metadata_**

```r
import_phyloseq("v3-v4.phyloseq")
```

---

## 2. Import metadata

Trong trường hợp chúng ta cần thay đổi metadata, import metadata mới thì có thể sử dụng chức năng *import_metdata*

```r
metadata <- import_metadata("sample-metadata.tsv")
```

---

## 3. Tính đa dạng alpha (alpha diversity)

- Tính toán các chỉ số đa dạng alpha như `Shannon`, `Observed`, `Chao1`, `Simpson`.
- Trả về: Data frame chứa giá trị alpha diversity.

```r
alpha <- cal_alpha(ps)
```

Chúng ta có thể tính toán tạo bảng statistic

```r
res3 <- cal_alpha_stat(alpha, metadata, group_col = "treatment", strata = "timeline")
```

---

## 4. Vẽ biểu đồ alpha diversity

- Vẽ boxplot và thực hiện thống kê theo các cặp.
- Đầu ra: Các biểu đồ riêng lẻ theo từng metrics, và biểu đồ tổng hợp.

```r
alphaplot<-plot_alpha(alpha = alpha, metadata = metadata, x = "treatment", facet = "timeline")
```

---

## 5. Vẽ biểu đồ beta diversity


- Vẽ biểu đồ thành phần vi khuẩn cấp độ phylum.

```r
betadata<-plot_beta(ps, color = "treatment",facet = "timeline",distance_method = "bray",method = "PCoA")
```

---

## 6. Vẽ biểu đồ thành phần Phylum

- Vẽ biểu đồ thành phần vi khuẩn cấp độ phylum.
  
```r
plot_phylum(ps, group_vars = c("treatment", "timeline"), facet = "timeline", x_var = "treatment")
```

---
## 7. Vẽ biểu đồ thành phần Genus

- Vẽ biểu đồ thành phần vi khuẩn cấp độ genus.
  
```r
plot_genus(ps, group_vars = c("treatment", "timeline"), top = 20, facet = "timeline", x_var = "treatment")
```

---

## 8. Vẽ biểu đồ khác biệt thống kê các chi vi sinh vật (lefse)

Chuyển dữ liệu _phyloseq object_ thành dạng _SummarizedExperiment_

```bash
se <- phyloseq_to_se(ps)
```

Tính toán LDA và vẽ biểu đồ lefse

```bash
out <- run_lefse_pairwise(se, classCol = "group", groups = c("Mus musculus_laboratory", "Mus musculus_wild"))
```
