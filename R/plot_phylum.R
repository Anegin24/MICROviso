#' Plot mean relative abundance at phylum level
#'
#' @param data A phyloseq object
#' @param group_vars Character vector of metadata variables to group by (e.g., c("treatment", "timeline"))
#' @param facet (Optional) Metadata variable to facet by (e.g., "timeline")
#' @param x_var Variable to be used as x-axis (must be one of group_vars)
#' @return A ggplot2 object
#' @export
plot_phylum <- function(data, group_vars, facet = NULL, x_var = NULL) {
  if (missing(group_vars) || length(group_vars) < 1) {
    stop("Please provide one or more grouping variables (e.g., group_vars = c('treatment'))")
  }

  if (!is.null(x_var) && !(x_var %in% group_vars)) {
    stop("x_var must be one of the group_vars")
  }

  if (is.null(x_var)) {
    x_var <- group_vars[1]
  }

  data_phylum <- tax_glom(data, taxrank = "Phylum")
  data_rel <- transform_sample_counts(data_phylum, function(x) x / sum(x))
  melted <- psmelt(data_rel)

  grouped <- melted %>%
    group_by(Phylum, !!!rlang::syms(group_vars)) %>%
    summarise(mean_abundance = mean(Abundance), .groups = "drop")

  p <- ggplot(grouped, aes(x = .data[[x_var]], y = mean_abundance, fill = Phylum)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = x_var, y = "Mean Relative Abundance", fill = "Phylum") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = guide_legend(ncol = 1)) +
    theme_q2r()

  if (!is.null(facet)) {
    if (!facet %in% colnames(grouped)) {
      stop("Facet variable not found in data: ", facet)
    }
    p <- p + facet_wrap(as.formula(paste("~", facet)))
  }

  return(p)
}
