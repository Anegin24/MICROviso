#' Plot Alpha Diversity Boxplots and Export Plots to Global Environment
#'
#' This function generates boxplots for specified alpha diversity metrics using `ggplot2`
#' and optionally facets them by a grouping variable. It also assigns each individual plot
#' and the combined patchwork plot to the global environment for direct access.
#'
#' @param alpha A `data.frame` containing alpha diversity metrics (e.g., from `alpha_cal()`).
#' @param metadata A `data.frame` containing sample metadata (must share an ID column with `alpha`).
#' @param metrics A character vector of metrics to visualize (e.g., c("Observed", "Shannon")).
#' @param x A string specifying the column in metadata to use for the x-axis (e.g., "treatment").
#' @param facet Optional. A string specifying the metadata column to facet by (e.g., "timeline").
#'
#' @return A named list with:
#' \describe{
#'   \item{plots}{A named list of individual `ggplot` objects, keyed by metric name.}
#'   \item{combined}{A combined `patchwork` plot of all provided metrics.}
#' }
#'
#' @details
#' This function automatically attempts to match the sample identifier columns between
#' `alpha` and `metadata` by comparing overlapping values. It uses `patchwork` to combine plots.
#' The result includes individual plots assigned to global environment (e.g., `Observed`, `Shannon`)
#' and a combined object named `alpha_combined`.
#'
#' @seealso \code{\link{alpha_cal}}, \code{\link[ggplot2]{ggplot}}, \code{\link[patchwork]{wrap_plots}}
#'
#' @import ggplot2
#' @importFrom dplyr inner_join
#' @importFrom ggpubr stat_compare_means
#' @export
plot_alpha <- function(alpha, metadata, metrics = c("Observed", "Shannon", "Chao1", "Simpson"),
                       x, facet = NULL) {

  match_col_exact <- function(df1, df2) {
    for (col1 in colnames(df1)) {
      for (col2 in colnames(df2)) {
        if (setequal(na.omit(df1[[col1]]), na.omit(df2[[col2]]))) {
          return(c(col1, col2))
        }
      }
    }
    return(NULL)
  }

  # Auto-match ID columns
  id_cols <- match_col_exact(alpha, metadata)
  if (is.null(id_cols)) {
    stop("❌ Cannot detect matching SampleID columns by exact content")
  }

  # Merge using detected ID columns
  merged <- dplyr::inner_join(alpha, metadata, by = setNames(id_cols[2], id_cols[1]))

  group_levels <- unique(na.omit(merged[[x]]))
  my_comparison <- if (length(group_levels) >= 2) combn(group_levels, 2, simplify = FALSE) else NULL

  plots <- list()
  for (metric in metrics) {
    p <- ggplot(merged, aes(x = .data[[x]], y = .data[[metric]])) +
      geom_boxplot(alpha = 0.5) +
      labs(title = paste(metric, "Diversity"), y = metric, x = x) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    if (!is.null(my_comparison)) {
      p <- p + ggpubr::stat_compare_means(comparisons = my_comparison,
                                          label = "p.format", tip.length = 0.01, step.increase = 0.1)
    }

    if (!is.null(facet)) {
      p <- p + facet_wrap(as.formula(paste("~", facet)))
    }

    plots[[metric]] <- p
    assign(metric, p, envir = globalenv())
  }

  if (all(c("Shannon", "Observed", "Chao1", "Simpson") %in% names(plots))) {
    combined <- (plots[["Shannon"]] | plots[["Observed"]]) /
      (plots[["Chao1"]] | plots[["Simpson"]]) +
      patchwork::plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(size = 16, face = "bold"))
  } else {
    combined <- patchwork::wrap_plots(plots)
  }

  assign("alpha_combined", combined, envir = globalenv())
  return(list(plots = plots, combined = combined))
}

