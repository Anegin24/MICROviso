#' Run post-hoc DESeq2 comparison at a given time point (Genus level)
#'
#' @param physeq A `phyloseq` object (raw, not yet aggregated).
#' @param time_point A string like "Week 0" or "Day 3" (must match sample_data()[[time_var]]).
#' @param group Name of grouping variable (default: "treatment").
#' @param time_var Name of time variable (default: "timeline").
#' @param comparison A character vector of length 2: c("groupA", "groupB").
#' @param padj_cutoff Adjusted p-value threshold (default: 0.05).
#'
#' @return A list with sig_table and plot.
#' @importFrom phyloseq tax_glom tax_table sample_data subset_samples phyloseq_to_deseq2
#' @importFrom DESeq2 DESeq results
#' @importFrom dplyr inner_join filter mutate arrange
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip labs theme_minimal scale_fill_manual
#' @export
DEseq2_pairwise <- function(physeq,
                            time_point,
                            group = "treatment",
                            time_var = "timeline",
                            comparison = c("groupA", "groupB"),
                            padj_cutoff = 0.05) {

  if (length(comparison) != 2) stop("❌ 'comparison' must be a vector of length 2")

  # 1. Aggregate to Genus level
  data_genus <- tax_glom(physeq, taxrank = "Genus")

  # 2. Extract taxonomy table with ASV IDs
  taxa <- as.data.frame(tax_table(data_genus))
  taxa$ASV <- rownames(taxa)

  # 3. Subset by time
  meta <- as.data.frame(sample_data(data_genus))
  idx <- meta[[time_var]] == time_point
  data_sub <- prune_samples(idx, data_genus)

  # 4. Run DESeq2
  dds <- phyloseq_to_deseq2(data_sub, as.formula(paste("~", group)))
  dds <- DESeq(dds)
  res <- results(dds, contrast = c(group, comparison[1], comparison[2]))
  res_df <- as.data.frame(res)
  res_df$ASV <- rownames(res_df)
  res_df <- inner_join(res_df, taxa, by = "ASV")

  # 5. Filter significant results
  res_sig <- res_df %>%
    filter(!is.na(padj), padj < padj_cutoff) %>%
    mutate(direction = if_else(log2FoldChange > 0, comparison[1], comparison[2])) %>%
    arrange(log2FoldChange)

  # 6. Plot
  p <- ggplot(res_sig, aes(x = reorder(Genus, log2FoldChange),
                           y = log2FoldChange, fill = direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("grey", "blue"),
                      labels = comparison) +
    labs(title = paste(time_point, ":", comparison[1], "vs", comparison[2]),
         x = "Genus", y = "log2 Fold Change", fill = NULL) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))

  return(list(sig_table = res_sig, plot = p))
}

