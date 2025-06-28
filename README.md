![image](https://github.com/user-attachments/assets/84f2237d-9614-4473-8054-01d7ba0782d6)

# MICROviso

**MICROviso** is an R package designed to support preprocessing and visualization of microbiome data using `phyloseq` objects. It provides a user-friendly workflow for:

- Data import
- Phylogenetic tree construction
- Alpha diversity calculation
- Visualization of microbial composition (Phylum and Genus levels)
- Group-based variation analysis
- Statistical testing

---

## Installation

Install dependencies from CRAN:

```r
install.packages("BiocManager")
install.packages("tidyverse")
install.packages("patchwork")
install.packages("ggplot2")
install.packages("RColorBrewer")
install.packages("devtools")
install.packages("ggpubr")
install.packages("SummarizedExperiment")
```

Install dependencies from Bioconductor:

```r
BiocManager::install("DESeq2")
BiocManager::install("phyloseq")
BiocManager::install("lefser")
BiocManager::install("msa")
BiocManager::install("phangorn")
```

Install MICROviso from GitHub:

```r
devtools::install_github("anegin24/MICROviso")
```

Load required libraries:

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

## 1. Import `phyloseq` Object

Load your `phyloseq` object (e.g., `ps`) into R and extract the abundance table, taxonomy, and metadata:

```r
import_phyloseq("v3-v4.phyloseq")
```

---

## 2. Import or Update Metadata

```r
metadata <- import_metadata("sample-metadata.tsv")
```

---

## 3. Construct Phylogenetic Tree

Use `make_tree()` to align sequences from `refseq(ps)`, construct a Neighbor-Joining tree, and optionally optimize with Maximum Likelihood (GTR). Tree visualization with `ggtree` is supported.

```r
tree <- make_tree(ps, msa_method = "ClustalOmega", optimize_ml = TRUE, plot_tree = TRUE)
```

Parameters:
- `msa_method`: one of `"ClustalW"`, `"Muscle"`, or `"ClustalOmega"`
- `optimize_ml`: `TRUE/FALSE` to enable ML optimization
- `plot_tree`: `TRUE/FALSE` to plot the tree

---

## 4. Alpha Diversity Calculation

```r
alpha <- cal_alpha(ps)
```

Generate statistical summary:

```r
res3 <- cal_alpha_stat(alpha, metadata, group_col = "treatment", strata = "timeline")
```

---

## 5. Alpha Diversity Plot

```r
alphaplot <- plot_alpha(alpha = alpha, metadata = metadata, x = "treatment", facet = "timeline")
```

---

## 6. Beta Diversity Plot

```r
betadata <- plot_beta(ps, color = "treatment", facet = "timeline", distance_method = "bray", method = "PCoA")
```

---

## 7. Phylum Composition Plot

```r
plot_phylum(ps, group_vars = c("treatment", "timeline"), facet = "timeline", x_var = "treatment")
```

---

## 8. Genus Composition Plot

```r
plot_genus(ps, group_vars = c("Sample"), top = 20, x_var = "Sample")

plot_genus(ps, group_vars = c("treatment", "timeline"), top = 20, facet = "timeline", x_var = "treatment")
```

---

## 9. Microbial Differential Abundance with LEfSe

```r
se <- phyloseq_to_se(ps)
out <- run_lefse_pairwise(se, classCol = "class", groups = c("ABX 0.2X_Week 6", "Control_Week 6"))
```

---

## 10. Differential Abundance with DESeq2

### Cross-sectional comparison

```r
DEseq2_cross(ps, group = NULL, comparison = NULL, padj_cutoff = 0.05)
```

### Detect genus-level Group × Time interaction (LRT test)

```r
DEseq2_global(ps, taxrank = "Genus", group = NULL, time_var = NULL, alpha = 0.05)
```

### Pairwise comparison

```r
DEseq2_pairwise(ps, "Week 0", group = "treatment", time_var = "timeline", 
                comparison = c("ABX 0.2X", "ABX 0.5X"), padj_cutoff = 0.05)
```

**Note:** If your group values contain spaces, run the following to replace spaces with underscores:

```r
ps@sam_data$treatment <- gsub(" ", "_", ps@sam_data$treatment)
sample_data(ps)$timeline <- gsub(" ", "_", sample_data(ps)$timeline)
```


