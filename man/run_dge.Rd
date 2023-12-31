% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_dge.R
\name{run_dge}
\alias{run_dge}
\title{Run Differential Gene Expression Analysis}
\usage{
run_dge(count_matrix_raw, metadata_raw, comparisons)
}
\arguments{
\item{count_matrix_raw}{A data frame representing the raw count matrix.
It must contain a column named "Geneid" that includes gene identifiers.}

\item{metadata_raw}{A data frame representing the metadata.
It must contain a column named "sample" with sample identifiers
and a column named "condition" with condition identifiers.}

\item{comparisons}{A list of vectors. Each vector includes two elements:
the identifiers of two conditions to be compared.}
}
\value{
A data frame representing the results of the differential gene expression analysis.
The data frame is sorted by the "Geneid" and "Comparison" columns in ascending order.
}
\description{
This function runs a DESeq2-based differential gene expression analysis
on a raw count matrix and associated metadata.
}
