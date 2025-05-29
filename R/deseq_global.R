# R/deseq_time_series.R

#' Run DESeq2 time-series interaction test
#'
#' @param physeq A phyloseq object (already aggregated if needed)
#' @param time_var Name of the time variable (user-defined)
#' @param group_var Name of the treatment group variable (user-defined)
#' @param taxrank Taxonomic rank for annotation (e.g. "Genus")
#' @param padj_cutoff Adjusted p-value threshold (default = 0.05)
#' @return A list with the result table and ggplot object
#' @export
deseq_global <- function(physeq, time_var, group_var, taxrank = "Genus", padj_cutoff = 0.05) {
  if (!inherits(physeq, "phyloseq")) stop("Input must be a phyloseq object")
  if (missing(time_var) || !(time_var %in% colnames(sample_data(physeq)))) stop("Time variable not found in sample_data")
  if (missing(group_var) || !(group_var %in% colnames(sample_data(physeq)))) stop("Group variable not found in sample_data")

  # Detect and convert time variable to numeric (remove non-numeric characters)
  time_char <- as.character(sample_data(physeq)[[time_var]])
  sample_data(physeq)[["time_numeric"]] <- suppressWarnings(as.numeric(gsub("[^0-9]+", "", time_char)))

  # Extract taxonomy
  taxa <- as.data.frame(tax_table(physeq))
  taxa$ASV <- rownames(taxa)

  # Run DESeq2 interaction model
  dds <- phyloseq_to_deseq2(physeq, as.formula(paste("~", group_var, "* time_numeric")))
  dds <- DESeq(dds, test = "LRT", reduced = as.formula(paste("~", group_var, "+ time_numeric")))
  res <- results(dds)
  res_df <- as.data.frame(res) %>%
    filter(!is.na(padj), padj < padj_cutoff) %>%
    arrange(padj)
  res_df$ASV <- rownames(res_df)
  res_sig <- inner_join(res_df, taxa, by = "ASV")

  # Plot
  p <- ggplot(res_sig, aes(x = reorder(.data[[taxrank]], -padj), y = -log10(padj))) +
    geom_bar(stat = "identity", fill = "coral") +
    coord_flip() +
    labs(title = paste("Significant interaction (", group_var, "×", time_var, ")"),
         x = taxrank, y = "-log10 Adjusted p-value") +
    theme_minimal()

  return(list(result_table = res_sig, plot = p))
}
