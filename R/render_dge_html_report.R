#' Render a Differential Gene Expression Analysis Report
#'
#' This function generates a report using R Markdown for the results of
#' a differential gene expression analysis. The report includes details about
#' the comparison performed, the results, metadata, normalized counts, raw counts,
#' biomaRt tables, and organism. The output is a .html file.
#'
#' @param comparison_name A string specifying the name of the comparison used in the analysis.
#' @param dge_results_in A data frame containing the results of the differential gene expression analysis.
#' @param metadata A data frame containing metadata associated with the analysis.
#' @param norm_counts A data frame containing the normalized count matrix.
#' @param raw_counts A data frame containing the raw count matrix.
#' @param organism A string specifying the organism used in the analysis.
#' @param output_name A string specifying the name of the output .html file.
#'
#' @importFrom rmarkdown render
#' @importFrom utils system.file
#'
#' @return No return value, called for side effects
render_dge_html_report <- function(comparison_name,
                          dge_results_in,
                          metadata,
                          norm_counts,
                          raw_counts,
                          organism,
                          output_name,
                          output_dir = ".") {

  # Define the path to the R Markdown template
  template_path <- system.file("rmarkdown", "templates", "dge-report-template.Rmd", package = "futuriandgeDownstream")

  # Render the R Markdown document
  rmarkdown::render(
    input = template_path,
    params = list(
      comparison_name = comparison_name,
      dge_results_in = dge_results_in,
      metadata = metadata,
      norm_counts = norm_counts,
      raw_counts = raw_counts
    ),
    output_file = output_name,
    output_dir = output_dir)
}
