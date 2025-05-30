#' Detect Significant Genus with Group × Time Interaction (LRT Test)
#'
#' @param physeq A `phyloseq` object.
#' @param taxrank The taxonomic rank to aggregate by (default: "Genus").
#' @param group The name of the group (e.g., treatment) variable in `sample_data()` (default: NULL).
#' @param time_var The name of the time variable in `sample_data()` (e.g., "timeline", "day") (default: NULL).
#' @param alpha Adjusted p-value cutoff (default: 0.05).
#'
#' @return A list with:
#' \describe{
#'   \item{sig_table}{A `data.frame` of significant results merged with taxonomy.}
#'   \item{plot}{A ggplot2 barplot of -log10(padj) for significant taxa.}
#' }
#' @importFrom dplyr inner_join filter arrange mutate
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip labs theme_minimal
#' @importFrom phyloseq tax_glom tax_table sample_data phyloseq_to_deseq2
#' @importFrom DESeq2 DESeq results
#' @importFrom rlang sym
#' @export
detect_interaction_genus <- function(physeq, taxrank = "Genus",
                                     group = NULL,
                                     time_var = NULL,
                                     alpha = 0.05) {
  if (is.null(group) || is.null(time_var)) {
    stop("❌ Please provide both 'group' and 'time_var' column names from sample_data().")
  }

  # Aggregate at desired taxonomic rank
  data_genus <- tax_glom(physeq, taxrank = taxrank)

  # Get taxonomy table
  taxa <- as.data.frame(tax_table(data_genus))
  taxa$ASV <- rownames(taxa)

  # Extract sample_data
  meta <- as.data.frame(sample_data(data_genus))

  # Ensure numeric time variable (extract digits from time_var)
  time_numeric <- as.numeric(gsub("[^0-9]", "", meta[[time_var]]))
  if (any(is.na(time_numeric))) stop("❌ Cannot convert time_var to numeric values.")

  meta$time_numeric <- time_numeric
  sample_data(data_genus)$time_numeric <- time_numeric

  # Dynamically build formula
  formula_full <- as.formula(paste("~", group, "* time_numeric"))
  formula_reduced <- as.formula(paste("~", group, "+ time_numeric"))

  # Run DESeq2 with LRT
  dds <- phyloseq_to_deseq2(data_genus, formula_full)
  dds <- DESeq(dds, test = "LRT", reduced = formula_reduced)
  res <- results(dds)
  res_df <- as.data.frame(res)

  # Filter significant results
  res_sig <- res_df %>%
    filter(!is.na(padj), padj < alpha) %>%
    arrange(padj)
  res_sig$ASV <- rownames(res_sig)

  # Merge with taxonomy
  res_sig <- inner_join(res_sig, taxa, by = "ASV")

  # Create plot
  p <- ggplot(res_sig, aes(x = reorder(!!sym(taxrank), -padj), y = -log10(padj))) +
    geom_bar(stat = "identity", fill = "coral") +
    coord_flip() +
    labs(
      title = paste("Significant", taxrank, "Interaction (", group, " × ", time_var, ")", sep = ""),
      x = taxrank, y = "-log10 Adjusted p-value"
    ) +
    theme_minimal()

  return(list(sig_table = res_sig, plot = p))
}
