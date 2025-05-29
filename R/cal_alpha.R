#' Calculate alpha diversity metrics
#'
#' @param data A phyloseq object
#' @param metrics A character vector of alpha diversity metrics to calculate (default: c("Observed", "Shannon", "Chao1", "Simpson"))
#' @return A data.frame of alpha diversity values
#' @importFrom phyloseq estimate_richness
#' @export
cal_alpha <- function(data, metrics = c("Observed", "Shannon", "Chao1", "Simpson")) {
  if (!inherits(data, "phyloseq")) stop("Input must be a phyloseq object")
  result <- estimate_richness(data, measures = metrics)
  result$SampleID <- rownames(result)
  rownames(result) <- NULL
  return(result)
}
