% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/render_dge_html_report.R
\name{render_dge_html_report}
\alias{render_dge_html_report}
\title{Render a Differential Gene Expression Analysis Report}
\usage{
render_dge_html_report(
  comparison_name,
  dge_results_in,
  metadata,
  norm_counts,
  raw_counts,
  organism,
  output_name,
  output_dir = "."
)
}
\arguments{
\item{comparison_name}{A string specifying the name of the comparison used in the analysis.}

\item{dge_results_in}{A data frame containing the results of the differential gene expression analysis.}

\item{metadata}{A data frame containing metadata associated with the analysis.}

\item{norm_counts}{A data frame containing the normalized count matrix.}

\item{raw_counts}{A data frame containing the raw count matrix.}

\item{organism}{A string specifying the organism used in the analysis.}

\item{output_name}{A string specifying the name of the output .html file.}

\item{output_dir}{A string specifying the output directory.}
}
\value{
No return value, called for side effects
}
\description{
This function generates a report using R Markdown for the results of
a differential gene expression analysis. The report includes details about
the comparison performed, the results, metadata, normalized counts, raw counts,
biomaRt tables, and organism. The output is a .html file.
}
