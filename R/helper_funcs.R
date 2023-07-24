#' Read in Feature Count Data
#'
#' This function reads a tab-separated file produced by the 'featureCounts' function
#' and removes the unnecessary columns.
#'
#' @param file A character string specifying the path to the input file.
#'
#' @return A data frame of the input file with the 'Chr', 'Start', 'End', 'Strand', and 'Length' columns removed.
#' @importFrom dplyr select
#' @importFrom readr read_tsv
read_in_feature_counts <- function(file) {
  readr::read_tsv(file, col_names = T, comment = "#") %>%
    dplyr::select(-Chr, -Start, -End, -Strand, -Length)
}

#' Export Data Table
#'
#' This function prepares the results of the differential expression analysis for export,
#' adding the sample names and the parameters used in the DESeq2 analysis.
#'
#' @param comparison_name A DESeq2 results object.
#'
#' @return A data frame ready for export. It includes information about the compared samples
#' and the parameters used for DESeq2 analysis.
#' @importFrom dplyr filter arrange mutate
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_split
#' @importFrom purrr map_dfr
export_data <- function(comparison_name) {
  # Check for compared samples
  meta_export_table <- dplyr::filter(meta, group %in%
    c(
      unlist(stringr::str_split(comparison_name@elementMetadata$description[2], " "))[6],
      unlist(stringr::str_split(comparison_name@elementMetadata$description[2], " "))[8]
    ))

  # Prepare the data for export
  export_table <- data.frame(comparison_name) %>%
    tibble::rownames_to_column(var = "Geneid") %>%
    arrange(Geneid) %>%
    mutate(
      Compared_samples = c(
        meta_export_table[["sample"]],
        rep(NA, nrow(comparison_name) - length(meta_export_table[["sample"]]))
      ),
      DESeq2_params = c(
        comparison_name@elementMetadata$description,
        rep(NA, nrow(comparison_name) - length(comparison_name@elementMetadata$description))
      )
    )

  # add hgnc
  export_table <- left_join(export_table, biomart_mapping)

  export_table
}