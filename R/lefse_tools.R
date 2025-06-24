#' Convert phyloseq to SummarizedExperiment (auto transpose if needed)
#'
#' This function converts a phyloseq object to a SummarizedExperiment object,
#' and automatically transposes the OTU table if taxa_are_rows is FALSE.
#'
#' @param ps A phyloseq object
#'
#' @return A SummarizedExperiment object
#' @export
#'
#' @examples
#' se <- phyloseq_to_se(ps)
phyloseq_to_se <- function(ps) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required.")
  }

  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) {
    otu <- t(otu)
  }

  meta <- as.data.frame(sample_data(ps))
  meta <- meta[colnames(otu), , drop = FALSE]
  rownames(meta) <- colnames(otu)

  tax <- as.data.frame(as(tax_table(ps), "matrix"))
  tax <- tax[rownames(otu), , drop = FALSE]
  rownames(tax) <- rownames(otu)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = otu),
    colData = meta,
    rowData = tax
  )
  return(se)
}


#' Run pairwise LEfSe analysis and return both result and plot
#'
#' This function runs LEfSe pairwise between two groups from a SummarizedExperiment,
#' and returns both the result table and the LDA plot.
#'
#' @param se A SummarizedExperiment object
#' @param classCol The name of the metadata column containing the group information
#' @param groups A character vector of 2 group levels to compare
#' @param lda_threshold LDA score cutoff
#'
#' @return A list with:
#'   \item{res}{data.frame of LEfSe results}
#'   \item{plot}{a ggplot2 object for LDA visualization}
#' @export
#'
#' @examples
#' # out <- run_lefse_pairwise(se, classCol = "group", groups = c("A", "B"))
#' # print(out$plot); head(out$res)
run_lefse_pairwise <- function(se, classCol, groups, lda_threshold = 2) {
  if (!(classCol %in% colnames(colData(se)))) {
    stop(paste("Column", classCol, "not found in colData(se)."))
  }

  if (length(groups) != 2) {
    stop("You must provide exactly two group names.")
  }

  subset_se <- se[, colData(se)[[classCol]] %in% groups]
  colData(subset_se)[[classCol]] <- droplevels(factor(colData(subset_se)[[classCol]]))

  subset_se <- lefser::relativeAb(subset_se)

  res <- lefser::lefser(
    subset_se,
    classCol = classCol,
    lda.threshold = lda_threshold
  )

  res$group_A <- groups[1]
  res$group_B <- groups[2]

  # Gắn taxonomy vào kết quả
  tax<-as.data.frame(tax_table(ps))
  tax$features <- rownames(tax)

  #Đảm bảo không lỗi legend title
  class_attr <- attr(res, "classes")
  res <- merge(res, tax, by = "features")
  attr(res, "classes") <- class_attr

  res <- res %>%
    mutate(
      features = ifelse(
        !is.na(Species) & Species != "" & Species != "unclassified",
        paste0(features, "_g_", Genus, ";s_", Species),
        paste0(features, "_g_", Genus)
      )
    )

  p <- lefser::lefserPlot(res, trim.names = TRUE)

  return(list(res = res, plot = p))
}
