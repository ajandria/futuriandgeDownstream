#' Remove miRNA-related Rows from Count Matrix
#'
#' This function filters out rows related to miRNAs in a given count matrix.
#' The filtering is based on gene identifiers that match those classified as miRNAs
#' in a provided Gene Transfer Format (GTF) file.
#'
#' @param count_matrix A data frame representing the count matrix.
#' It must contain a column named "Geneid" that includes gene identifiers.
#' @param organism A character string indicating the organism.
#' It is used to properly import a respective GTF file.
#'
#' @return A data frame. It represents the count matrix after miRNA-related rows have been removed.
#' The data frame is sorted by the "Geneid" column in ascending order.
#'
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' count_matrix <- data.frame(
#'   Geneid = c("gene1", "gene2", "gene3"),
#'   sample1 = c(10, 15, 5),
#'   sample2 = c(20, 10, 10)
#' )
#' gtf_path <- system.file("extdata", "example.gtf", package = "yourPackageName")
#' remove_mirna(count_matrix, gtf_path)
#' }
remove_mirna <- function(count_matrix, organism) {
  if (organism == "homo_sapiens") {
    gtf_path <- system.file("extdata", "Homo_sapiens.GRCh38.110.chr.gtf.gz", package = "futuriandgeDownstream")
  } else if (organism == "mus_musculus") {
    gtf_path <- system.file("extdata", "Mus_musculus.GRCm39.110.chr.gtf.gz", package = "futuriandgeDownstream")
  } else {
    stop(message("Organism not supported..."))
  }

  # Import GTF and filter for miRNA rows
  gtf_df_filtered <- rtracklayer::import(gtf_path) %>%
    as.data.frame() %>%
    dplyr::select(.data$gene_id, .data$gene_biotype) %>%
    dplyr::filter(.data$gene_biotype == "miRNA")

  # Define negation of `%in%`
  `%nin%` <- Negate(`%in%`)

  # Remove miRNA rows and sort by gene ID
  count_matrix %>%
    dplyr::filter(.data$Geneid %nin% gtf_df_filtered$gene_id) %>%
    dplyr::arrange(.data$Geneid)
}
