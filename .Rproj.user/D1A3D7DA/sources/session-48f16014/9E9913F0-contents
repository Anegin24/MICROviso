#' Plot mean relative abundance at genus level
#'
#' @param data A phyloseq object
#' @param group_vars Character vector of metadata variables to group by (e.g., c("treatment", "timeline"))
#' @param top Integer number of top genera to display, remaining will be grouped as 'Others'
#' @param facet (Optional) Metadata variable to facet by (e.g., "timeline")
#' @param x_var Variable to be used as x-axis (must be one of group_vars)
#' @return A ggplot2 object
#' @export
plot_genus <- function(data, group_vars, top = 20, facet = NULL, x_var = NULL) {
  if (missing(group_vars) || length(group_vars) < 1) {
    stop("Please provide one or more grouping variables (e.g., group_vars = c('treatment'))")
  }

  if (!is.null(x_var) && !(x_var %in% group_vars)) {
    stop("x_var must be one of the group_vars")
  }

  if (is.null(x_var)) {
    x_var <- group_vars[1]
  }

  data_genus <- tax_glom(data, taxrank = "Genus")
  data_rel <- transform_sample_counts(data_genus, function(x) x / sum(x))
  melted <- psmelt(data_rel)

  # Summarise by genus and grouping vars
  grouped <- melted %>%
    group_by(Genus, !!!rlang::syms(group_vars)) %>%
    summarise(mean_abundance = mean(Abundance), .groups = "drop")

  # Identify top genera
  a <- grouped %>%
    group_by(Genus) %>%
    summarise(max = max(mean_abundance), .groups = "drop") %>%
    arrange(desc(max))
  top_genera <- head(a, top)

  # Recode 'Others'
  grouped <- grouped %>%
    mutate(Genus = if_else(Genus %in% top_genera$Genus, Genus, "Others")) %>%
    group_by(Genus, !!!rlang::syms(group_vars)) %>%
    summarise(mean_abundance = sum(mean_abundance), .groups = "drop")

  # Build plot
  p <- ggplot(grouped, aes(x = .data[[x_var]], y = mean_abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = x_var, y = "Mean Relative Abundance", fill = "Genus") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = guide_legend(ncol = 1)) +
    theme_q2r() +
    scale_fill_manual(name = NULL, values = c(
      RColorBrewer::brewer.pal(8, "Dark2"),
      "blue", "gray30", "yellow", "red", "darkmagenta", "green", "pink2", "darkgreen",
      "cyan", "orange", "purple", "lightblue", "lightgreen", "salmon", "gold", "darkred", "navy", "orchid"
    ))

  # Facet if specified
  if (!is.null(facet)) {
    if (!facet %in% colnames(grouped)) {
      stop("Facet variable not found in data: ", facet)
    }
    p <- p + facet_wrap(as.formula(paste("~", facet)))
  }

  return(p)
}
