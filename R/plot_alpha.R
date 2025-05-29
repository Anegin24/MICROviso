#' Calculate alpha diversity metrics
#'
#' @param data A phyloseq object
#' @param metrics A character vector of alpha diversity metrics to calculate (default: c("Observed", "Shannon", "Chao1", "Simpson"))
#' @return A data.frame of alpha diversity values
#' @export
alpha_cal <- function(data, metrics = c("Observed", "Shannon", "Chao1", "Simpson")) {
  if (!inherits(data, "phyloseq")) stop("Input must be a phyloseq object")
  result <- estimate_richness(data, measures = metrics)
  result$SampleID <- rownames(result)
  return(result)
}

#' Plot alpha diversity boxplots and assign individual plots to global environment
#'
#' @param alpha A data.frame of alpha diversity values
#' @param metadata A data.frame of sample metadata
#' @param metrics Character vector of metrics to plot (e.g., c("Observed", "Shannon"))
#' @param x Variable name in metadata to use for x-axis (e.g., "timeline")
#' @param facet Variable name in metadata to facet by (optional)
#' @return A list of individual ggplot objects and a combined patchwork plot
#' @export
plot_alpha <- function(alpha, metadata, metrics = c("Observed", "Shannon", "Chao1", "Simpson"), x, facet = NULL) {
  match_col <- function(df1, df2) {
    for (col1 in colnames(df1)) {
      for (col2 in colnames(df2)) {
        if (all(na.omit(df1[[col1]]) %in% df2[[col2]])) {
          return(c(col1, col2))
        }
      }
    }
    return(NULL)
  }

  id_cols <- match_col(alpha, metadata)
  if (is.null(id_cols)) stop("❌ Cannot detect matching SampleID columns by content")

  merged <- dplyr::inner_join(alpha, metadata, by = setNames(id_cols[2], id_cols[1]))

  plots <- list()
  for (metric in metrics) {
    p <- ggplot(merged, aes(x = .data[[x]], y = .data[[metric]])) +
      geom_boxplot(alpha = 0.5) +
      labs(title = paste(metric, "Diversity"), y = metric, x = x) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if (!is.null(facet)) {
      p <- p + facet_wrap(as.formula(paste("~", facet)))
    }
    plots[[metric]] <- p

    # Assign plot to Global Environment with name based on metric
    assign(metric, p, envir = globalenv())
  }

  # Combine plots using patchwork
  if (all(c("Shannon", "Observed", "Chao1", "Simpson") %in% names(plots))) {
    combined <- (plots[["Shannon"]] | plots[["Observed"]]) / (plots[["Chao1"]] | plots[["Simpson"]]) +
      patchwork::plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(size = 16, face = "bold"))
  } else {
    combined <- patchwork::wrap_plots(plots)
  }

  assign("alpha_combined", combined, envir = globalenv())
  return(list(plots = plots, combined = combined))
}
