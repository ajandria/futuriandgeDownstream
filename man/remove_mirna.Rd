% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_mirna.R
\name{remove_mirna}
\alias{remove_mirna}
\title{Remove miRNA-related Rows from Count Matrix}
\usage{
remove_mirna(count_matrix, organism)
}
\arguments{
\item{count_matrix}{A data frame representing the count matrix.
It must contain a column named "Geneid" that includes gene identifiers.}

\item{organism}{A character string indicating the organism.
It is used to properly import a respective GTF file.}
}
\value{
A data frame. It represents the count matrix after miRNA-related rows have been removed.
The data frame is sorted by the "Geneid" column in ascending order.
}
\description{
This function filters out rows related to miRNAs in a given count matrix.
The filtering is based on gene identifiers that match those classified as miRNAs
in a provided Gene Transfer Format (GTF) file.
}
\examples{
\dontrun{
count_matrix <- data.frame(
  Geneid = c("gene1", "gene2", "gene3"),
  sample1 = c(10, 15, 5),
  sample2 = c(20, 10, 10)
)
gtf_path <- system.file("extdata", "example.gtf", package = "yourPackageName")
remove_mirna(count_matrix, gtf_path)
}
}
