#' Compute Pairwise Statistical Tests for Alpha Diversity Metrics
#'
#' This function calculates pairwise statistical comparisons (e.g., Wilcoxon or t-test)
#' for specified alpha diversity metrics across groups defined in sample metadata.
#' It supports optional stratification by a second variable (e.g., timeline).
#'
#' @param alpha A `data.frame` containing alpha diversity metrics (e.g., output of `alpha_cal()`).
#' @param metadata A `data.frame` with sample-level metadata, including grouping and optional stratification variables.
#' @param group_col A string specifying the column in `metadata` used for group comparison (e.g., "treatment").
#' @param metrics A character vector of alpha metrics to compare (default: `c("Observed", "Shannon", "Chao1", "Simpson")`).
#' @param method Statistical test to use: `"wilcox"` (default) or `"t.test"`.
#' @param p.adjust.method Method for p-value adjustment (default: `"BH"` for Benjamini-Hochberg).
#' @param strata Optional. A string specifying a column in `metadata` used for stratification (e.g., "timeline").
#'
#' @return A named list where each element is a tibble of pairwise test results for one alpha metric.
#'
#' @details
#' The function automatically detects and matches sample identifier columns between `alpha` and `metadata`.
#' It performs pairwise statistical comparisons for each alpha diversity metric using `rstatix::pairwise_wilcox_test()`
#' or `rstatix::pairwise_t_test()` and optionally stratifies the comparisons by a second grouping variable.
#'
#' @seealso \code{\link{plot_alpha}}, \code{\link[rstatix]{pairwise_wilcox_test}}, \code{\link[rstatix]{pairwise_t_test}}
#'
#' @importFrom dplyr inner_join group_by mutate
#' @importFrom rstatix pairwise_wilcox_test pairwise_t_test
#' @export
cal_alpha_stat <- function(alpha, metadata, group_col,
                                 metrics = c("Observed", "Shannon", "Chao1", "Simpson"),
                                 method = "wilcox", p.adjust.method = "BH",
                                 strata = NULL) {

  # Internal function to match ID columns exactly
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

  # Match ID columns
  id_cols <- match_col_exact(alpha, metadata)
  if (is.null(id_cols)) stop("❌ Cannot match sample ID columns")

  # Merge data
  merged <- dplyr::inner_join(alpha, metadata, by = setNames(id_cols[2], id_cols[1]))

  results <- list()

  # Loop over each metric
  for (metric in metrics) {
    data_metric <- merged[, c(group_col, strata, metric), drop = FALSE]
    colnames(data_metric)[which(colnames(data_metric) == metric)] <- "value"

    # Build formula object
    formula_obj <- as.formula(paste("value ~", group_col))

    # Apply pairwise test (with or without stratification)
    if (!is.null(strata)) {
      stat_result <- data_metric %>%
        dplyr::group_by(.data[[strata]]) %>%
        {
          if (method == "wilcox") {
            rstatix::pairwise_wilcox_test(., formula = formula_obj, p.adjust.method = p.adjust.method)
          } else if (method == "t.test") {
            rstatix::pairwise_t_test(., formula = formula_obj, p.adjust.method = p.adjust.method)
          } else {
            stop("Unsupported method. Use 'wilcox' or 't.test'.")
          }
        } %>%
        dplyr::mutate(metric = metric)
    } else {
      stat_result <- {
        if (method == "wilcox") {
          rstatix::pairwise_wilcox_test(data_metric, formula = formula_obj, p.adjust.method = p.adjust.method)
        } else if (method == "t.test") {
          rstatix::pairwise_t_test(data_metric, formula = formula_obj, p.adjust.method = p.adjust.method)
        } else {
          stop("Unsupported method. Use 'wilcox' or 't.test'.")
        }
      } %>%
        dplyr::mutate(metric = metric)
    }

    results[[metric]] <- stat_result
  }

  return(results)
}
