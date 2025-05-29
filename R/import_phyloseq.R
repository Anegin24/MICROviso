#' Import phyloseq components and assign as global variables
#'
#' @param phyloseq_path Path to the saved RDS phyloseq object
#' @return NULL (assigns ps, taxonomy, table, metadata to global environment)
#' @export
import_phyloseq <- function(phyloseq_path) {
  if (!file.exists(phyloseq_path)) {
    stop("File not found: ", phyloseq_path)
  }

  message("Reading phyloseq object from: ", phyloseq_path)
  ps <- readRDS(phyloseq_path)

  taxonomy <- as.data.frame(phyloseq::otu_table(ps))
  taxonomy <- tibble::rownames_to_column(taxonomy, var = "FeatureID")

  table <- as.data.frame(phyloseq::tax_table(ps))
  table <- tibble::rownames_to_column(table, var = "FeatureID")

  metadata <- as.data.frame(phyloseq::sample_data(ps))
  metadata$SampleID <- rownames(metadata)
  rownames(metadata) <- NULL

  assign("ps", ps, envir = .GlobalEnv)
  assign("taxonomy", taxonomy, envir = .GlobalEnv)
  assign("table", table, envir = .GlobalEnv)
  assign("metadata", metadata, envir = .GlobalEnv)

  message("Variables created in global environment: ps, taxonomy, table, metadata")
  invisible(NULL)
}
