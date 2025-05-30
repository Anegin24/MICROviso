![image](https://github.com/user-attachments/assets/c928c978-755d-45cc-99a4-4fd140fe864d)# MICROVISO R Package
![image](https://github.com/user-attachments/assets/3628b8ca-1042-447b-8bbf-88692529ddb6)

## Mục đích

Bộ hàm hỗ trợ tiền xử lý và trực quan hóa dữ liệu vi sinh vật từ đối tượng `phyloseq`, bao gồm:

- Import dữ liệu
- Tính toán đa dạng alpha
- Biểu đồ thành phần vi sinh vật theo cấp độ phân loại (phylum, genus)
- Kiểm tra sự biến đổi theo nhóm
- Tính toán thống kê
---
## Setup
Download: 
Installation

Dependencies:

```bash
install.packages("tidyverse")
install.packages("patchwork")
install.packages("ggplot2")
install.packages("BiocManager") 
BiocManager::install("DESeq2")
BiocManager::install("phyloseq")
install.packages("RColorBrewer")
install.packages("devtools")
Install.packages("ggpubr") 
```

MICROviso:

```bash
devtools::install_github("anegin24/MICROviso")
```

## 📁 1. Import dữ liệu từ `phyloseq`

### `import_physeq(path)`
- Đọc file `.rds` chứa đối tượng `phyloseq` và tách thành `otu_table`, `tax_table`, `sample_data`.
- Trả về: List gồm `ps`, `taxonomy`, `table`.

```r
import_physeq("v3-v4.phyloseq")
ps 
taxonomy 
table
metadata
```

---

## 📁 2. Import metadata

### `import_metadata(path)`
- Đọc file metadata (hỗ trợ `.csv`, `.tsv`, `.txt`, `.xlsx`).
- Trả về: Data frame metadata.

```r
metadata <- import_metadata("sample-metadata.tsv")
```

---

## 🧮 3. Tính đa dạng alpha

### `alpha_cal(data, metrics)`
- Tính toán các chỉ số đa dạng alpha như `Shannon`, `Observed`, `Chao1`, `Simpson`.
- Trả về: Data frame chứa giá trị alpha diversity.

```r
alpha <- alpha_cal(ps)
```

---

## 📊 4. Vẽ biểu đồ alpha diversity

### `plot_alpha(alpha, metadata, x, facet, metrics)`
- Vẽ boxplot và biểu đồ tổng hợp alpha diversity.
- Trả về: List gồm `plots` và `combined`.

```r
alphaplot<-plot_alpha(alpha = alpha, metadata = metadata, x = "treatment", facet = "timeline")
```

---

## 🧬 5. Vẽ biểu đồ thành phần Phylum

### `plot_phylum(data, group_vars, facet = NULL, x_var = NULL)`
- Vẽ biểu đồ thành phần vi khuẩn cấp độ phylum.

```r
plot_phylum(ps, group_vars = c("treatment", "timeline"), facet = "timeline", x_var = "treatment")
```

---

## 🧬 6. Vẽ biểu đồ thành phần Genus

### `plot_genus(data, group_vars, top = 20, facet = NULL, x_var = NULL)`
- Vẽ biểu đồ thành phần genus phổ biến nhất.

```r
plot_genus(ps, group_vars = c("treatment", "timeline"), top = 20, facet = "timeline", x_var = "treatment")
```

---

## 📌 Gợi ý

- Lưu biểu đồ: `ggsave("Observed.pdf", Observed)`
- Lọc mẫu: `dplyr::filter(metadata, treatment == "Control")`
