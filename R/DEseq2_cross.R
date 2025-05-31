#' Differential Abundance Analysis (Cross-sectional) using DESeq2
#'
#' This function performs differential abundance analysis using DESeq2
#' for a cross-sectional study design where there is only one grouping variable
#' (e.g., "Early" vs "Late"), without repeated measures or time dependency.
#'
#' @param physeq A `phyloseq` object containing taxonomic count data, taxonomy, and sample metadata.
#' @param group Character. The name of the grouping variable in `sample_data(physeq)` (e.g., "When"). Default is NULL.
#' @param comparison Character vector of length 2. The two group values to compare (e.g., `c("Early", "Late")`). Default is NULL.
#' @param padj_cutoff Numeric. Adjusted p-value threshold to filter significant results (default = 0.05).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{sig_table}{A data frame of significantly differentially abundant genera.}
#'   \item{plot}{A ggplot2 bar plot showing log2 fold changes of significant genera.}
#' }
#'
#' @importFrom phyloseq tax_glom tax_table sample_data prune_samples phyloseq_to_deseq2
#' @importFrom DESeq2 DESeq results
#' @importFrom dplyr filter mutate arrange inner_join
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip scale_fill_manual labs theme_minimal theme element_text
#'
#' @export
DEseq2_cross <- function(physeq,
                         group = NULL,
                         comparison = NULL,
                         padj_cutoff = 0.05) {

  # Kiểm tra đầu vào
  if (is.null(group)) stop("❌ Please specify the 'group' variable (e.g., 'When').")
  if (is.null(comparison) || length(comparison) != 2) {
    stop("❌ 'comparison' must be a character vector of length 2, e.g., c('Early', 'Late').")
  }

  # Tổng hợp cấp độ Genus
  data_genus <- tax_glom(physeq, taxrank = "Genus")
  taxa <- as.data.frame(tax_table(data_genus))
  taxa$ASV <- rownames(taxa)

  # Lọc giữ 2 nhóm
  meta <- as.data.frame(sample_data(data_genus))
  if (!group %in% colnames(meta)) {
    stop(paste0("❌ Group variable '", group, "' not found in sample_data."))
  }
  keep <- meta[[group]] %in% comparison
  if (sum(keep) == 0) stop("❌ No samples matched the comparison groups.")

  data_sub <- prune_samples(keep, data_genus)
  meta_sub <- sample_data(data_sub)

  # Đặt lại levels
  meta_sub[[group]] <- factor(meta_sub[[group]], levels = comparison)
  sample_data(data_sub) <- meta_sub

  # DESeq2
  dds <- phyloseq_to_deseq2(data_sub, as.formula(paste("~", group)))
  dds <- DESeq(dds)
  res <- results(dds, contrast = c(group, comparison[1], comparison[2]))
  res_df <- as.data.frame(res)
  res_df$ASV <- rownames(res_df)

  # Ghép với taxonomy
  res_df <- inner_join(res_df, taxa, by = "ASV")

  # Lọc kết quả có ý nghĩa
  res_sig <- res_df %>%
    filter(!is.na(padj), padj < padj_cutoff) %>%
    mutate(direction = if_else(log2FoldChange > 0, comparison[1], comparison[2])) %>%
    arrange(log2FoldChange)

  # Vẽ biểu đồ
  p <- ggplot(res_sig, aes(x = reorder(Genus, log2FoldChange),
                           y = log2FoldChange,
                           fill = direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("grey", "blue"), labels = comparison) +
    labs(title = paste(comparison[1], "vs", comparison[2]),
         x = "Genus", y = "log2 Fold Change", fill = NULL) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))

  return(list(sig_table = res_sig, plot = p))
}
