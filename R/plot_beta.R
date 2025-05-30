#' Plot Beta Diversity using Ordination Methods
#'
#' This function calculates distance matrices (e.g., Bray-Curtis, Jaccard, UniFrac),
#' performs ordination (e.g., PCoA, NMDS), and plots the ordination with optional
#' coloring, faceting, and ellipse overlays.
#'
#' @param data A `phyloseq` object.
#' @param color A sample metadata variable to color the points.
#' @param facet Optional. A metadata variable to facet the plot by.
#' @param distance_method Character. Distance metric to use (e.g., "bray", "jaccard", "unifrac").
#' @param method Character. Ordination method (e.g., "PCoA", "NMDS", "DCA"). Default is "PCoA".
#' @param weighted Logical. Whether to use weighted version for UniFrac or Bray. Default is FALSE.
#'
#' @return A named list with:
#' \describe{
#'   \item{distance}{The distance matrix (of class `dist`).}
#'   \item{ordination}{The ordination object.}
#'   \item{plot}{The resulting `ggplot2` ordination plot.}
#' }
#'
#' @import ggplot2
#' @importFrom phyloseq distance ordinate plot_ordination
#' @export
plot_beta <- function(data, color, facet = NULL,
                      distance_method = "bray",
                      method = "PCoA",
                      weighted = FALSE) {
  # Calculate distance matrix
  dist <- phyloseq::distance(data, method = distance_method, weighted = weighted)

  # Perform ordination
  ordination <- ordinate(data, method = method, distance = dist)

  # Create ordination plot
  p <- plot_ordination(data, ordination, color = color) +
    stat_ellipse(type = "t", linetype = 2) +
    theme(aspect.ratio = 1)

  # Add faceting if specified
  if (!is.null(facet)) {
    p <- p + facet_wrap(as.formula(paste("~", facet)))
  }

  # Return all components
  return(list(
    distance = dist,
    ordination = ordination,
    plot = p
  ))
}
