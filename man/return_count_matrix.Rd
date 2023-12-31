% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/return_count_matrix.R
\name{return_count_matrix}
\alias{return_count_matrix}
\title{Merge Count Matrices from Different Paths}
\usage{
return_count_matrix(downstream_ch_file)
}
\arguments{
\item{downstream_ch_file}{A character string indicating the path to a CSV file.
The CSV file should contain a column named "count_matrix_path" with paths to count matrices.}
}
\value{
A data frame representing the merged count matrix.
The column names of the matrix are modified to remove any substring after (and including) "_".
}
\description{
This function reads a CSV file which contains paths to different count matrices.
It then reads and merges these count matrices into a single matrix.
}
\examples{
\dontrun{
downstream_ch_file <- system.file("extdata", "example_path_file.csv", package = "yourPackageName")
return_count_matrix(downstream_ch_file)
}
}
