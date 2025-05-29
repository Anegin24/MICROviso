#' Import phyloseq components and assign as global variables
#'
#' This function reads a phyloseq object from an RDS file and extracts its main components:
#' OTU table, taxonomic table, and metadata. These components are then assigned as global variables
#' (`ps`, `otu_table`, `tax_table`, `metadata`) for easier downstream analysis.
#'
#' @param phyloseq_path Character. Path to the saved RDS phyloseq object.
#'
#' @return NULL. This function is called for its side effects (assigning variables to the global environment).
#'
#' @examples
#' \dontrun{
#' import_phyloseq("data/phyloseq_obj.rds")
#' }
#'
#' @export
import_phyloseq <- function(phyloseq_path) {
  if (!file.exists(phyloseq_path)) {
    stop("File not found: ", phyloseq_path)
  }

  message(" Reading phyloseq object from: ", phyloseq_path)
  ps <- readRDS(phyloseq_path)

  # Extract OTU table
  otu_df <- as.data.frame(otu_table(ps))
  otu_df <- tibble::rownames_to_column(otu_df, var = "FeatureID")

  # Extract taxonomy table
  tax_df <- as.data.frame(tax_table(ps))
  tax_df <- tibble::rownames_to_column(tax_df, var = "FeatureID")

  # Extract metadata
  meta_df <- as.data.frame(sample_data(ps))
  meta_df$SampleID <- rownames(meta_df)
  rownames(meta_df) <- NULL

  # Assign to global environment
  assign("ps", ps, envir = .GlobalEnv)
  assign("otu_table", otu_df, envir = .GlobalEnv)
  assign("tax_table", tax_df, envir = .GlobalEnv)
  assign("metadata", meta_df, envir = .GlobalEnv)

  message("Variables created in global environment: ps, otu_table, tax_table, metadata")
  invisible(NULL)
}
