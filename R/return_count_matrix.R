#' Merge Count Matrices from Different Paths
#'
#' This function reads a CSV file which contains paths to different count matrices.
#' It then reads and merges these count matrices into a single matrix.
#'
#' @param downstream_ch_file A character string indicating the path to a CSV file.
#' The CSV file should contain a column named "count_matrix_path" with paths to count matrices.
#'
#' @return A data frame representing the merged count matrix.
#' The column names of the matrix are modified to remove any substring after (and including) "_".
#'
#' @export
#'
#' @examples
#' \dontrun{
#' downstream_ch_file <- system.file("extdata", "example_path_file.csv", package = "yourPackageName")
#' return_count_matrix(downstream_ch_file)
#' }
return_count_matrix <- function(downstream_ch_file) {
  # Read meta downstream file
  meta_downstream_file <- readr::read_csv(downstream_ch_file)

  # Merge count matrices from different paths
  merged_count_matrix <- meta_downstream_file[["count_matrix_path"]] %>%
    purrr::map(read_in_feature_counts) %>%
    purrr::reduce(dplyr::inner_join)

  # Clean column names
  colnames(merged_count_matrix) <- gsub("_.*", "", colnames(merged_count_matrix))

  merged_count_matrix
}
