#' Plot Mean Relative Abundance at Phylum Level
#'
#' This function plots the mean relative abundance of phyla in a `phyloseq` object
#' using stacked bar charts. It allows grouping by one or more metadata variables
#' and optional faceting.
#'
#' @param data A `phyloseq` object containing taxonomic and sample metadata.
#' @param group_vars A character vector of metadata variables to group by
#'   (e.g., \code{c("treatment", "timeline")}).
#' @param facet Optional. A metadata variable to facet by (e.g., \code{"timeline"}).
#' @param x_var Optional. A variable to use as the x-axis. Must be one of \code{group_vars}.
#'   If not provided, the first variable in \code{group_vars} will be used.
#'
#' @return A `ggplot2` object showing stacked bar plots of mean phylum relative abundance.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Aggregates data at the Phylum level using \code{tax_glom}.
#'   \item Converts counts to relative abundance.
#'   \item Groups and summarizes mean relative abundance across the provided metadata variables.
#'   \item Optionally facets the plot by an additional variable.
#' }
#'
#' @examples
#' \dontrun{
#' plot_phylum(data = ps, group_vars = c("treatment", "timeline"), x_var = "treatment", facet = "timeline")
#' }
#'
#' @seealso \code{\link{plot_genus}}, \code{\link[phyloseq]{tax_glom}}, \code{\link[ggplot2]{ggplot}}
#' @import ggplot2
#' @importFrom phyloseq tax_glom transform_sample_counts psmelt
#' @importFrom dplyr group_by summarise mutate arrange
#' @importFrom rlang .data syms
#' @importFrom RColorBrewer brewer.pal
#' @importFrom magrittr %>%
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
    guides(fill = guide_legend(ncol = 1))

  if (!is.null(facet)) {
    if (!facet %in% colnames(grouped)) {
      stop("Facet variable not found in data: ", facet)
    }
    p <- p + facet_wrap(as.formula(paste("~", facet)))
  }

  return(p)
}
