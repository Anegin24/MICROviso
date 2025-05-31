#' Plot Mean Relative Abundance at Genus Level
#'
#' This function creates a stacked bar plot of the mean relative abundance of genera
#' in a `phyloseq` object, allowing grouping by multiple metadata variables, optional
#' faceting, and customizable top genera display.
#'
#' @param data A `phyloseq` object containing OTU, taxonomy, and metadata.
#' @param group_vars Character vector of metadata variables to group by
#'   (e.g., \code{c("treatment", "timeline")}).
#' @param top Integer number of top genera to display. All other genera will be grouped as \code{"Others"}.
#' @param facet Optional. A metadata variable name used for faceting the plot (e.g., \code{"timeline"}).
#' @param x_var Optional. A variable to use for the x-axis (must be one of \code{group_vars}).
#'   If not provided, the first element of \code{group_vars} is used by default.
#'
#' @return A `ggplot2` object showing stacked bar plots of mean genus relative abundance.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Aggregates abundance data to the Genus level.
#'   \item Transforms counts to relative abundances.
#'   \item Groups by selected metadata variables and computes mean abundance.
#'   \item Displays the top N most abundant genera, grouping others as "Others".
#'   \item Allows optional faceting using another metadata variable.
#' }
#'
#' @examples
#' \dontrun{
#' plot_genus(data = ps, group_vars = c("treatment", "timeline"), x_var = "treatment", facet = "timeline")
#' }
#'
#' @seealso \code{\link{plot_phylum}}, \code{\link[phyloseq]{tax_glom}}, \code{\link[ggplot2]{ggplot}}
#'
#' @import ggplot2
#' @importFrom phyloseq tax_glom transform_sample_counts psmelt
#' @importFrom dplyr group_by summarise mutate arrange
#' @importFrom rlang .data syms
#' @importFrom RColorBrewer brewer.pal
#' @importFrom magrittr %>%
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
    group_by(Genus, !!!syms(group_vars)) %>%
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
    group_by(Genus, !!!syms(group_vars)) %>%
    summarise(mean_abundance = sum(mean_abundance), .groups = "drop")
  grouped <- grouped %>%
    mutate(
      Genus = forcats::fct_reorder(Genus, mean_abundance),
      Genus = forcats::fct_relevel(Genus, "Others", after = Inf)
    )

  # Build plot
  p <- ggplot(grouped, aes(x = .data[[x_var]], y = mean_abundance, fill = Genus)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = x_var, y = "Mean Relative Abundance", fill = "Genus") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = guide_legend(ncol = 1)) +
    scale_fill_manual(name = NULL, values = c(
      brewer.pal(8, "Dark2"),
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
