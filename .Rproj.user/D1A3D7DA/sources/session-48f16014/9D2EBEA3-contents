#' Import metadata from CSV, TSV, TXT, or Excel
#'
#' @param metadata_path Path to the metadata file
#' @return A data.frame containing metadata
#' @export
import_metadata <- function(metadata_path) {
  if (!file.exists(metadata_path)) {
    stop("File not found: ", metadata_path)
  }

  message("Importing metadata from: ", metadata_path)
  ext <- tolower(tools::file_ext(metadata_path))

  if (ext %in% c("csv", "tsv", "txt", "tabular")) {
    # Smart delimiter detection for text files
    first_line <- readLines(metadata_path, n = 1)
    if (grepl(",", first_line)) {
      metadata <- readr::read_csv(metadata_path)
    } else if (grepl("\t", first_line)) {
      metadata <- readr::read_tsv(metadata_path)
    } else {
      stop("Could not auto-detect delimiter in text file.")
    }
  } else if (ext %in% c("xls", "xlsx")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required to read Excel files. Please install it with install.packages('readxl').")
    }
    metadata <- readxl::read_excel(metadata_path)
  } else {
    stop("Unsupported file format: must be .csv, .tsv, .txt, .tabular, .xls, or .xlsx")
  }

  message("Metadata imported: ", nrow(metadata), " rows and ", ncol(metadata), " columns.")
  return(metadata)
}
