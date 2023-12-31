% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_funcs.R
\name{read_in_feature_counts}
\alias{read_in_feature_counts}
\title{Read in Feature Count Data}
\usage{
read_in_feature_counts(file)
}
\arguments{
\item{file}{A character string specifying the path to the input file.}
}
\value{
A data frame of the input file with the 'Chr', 'Start', 'End', 'Strand', and 'Length' columns removed.
}
\description{
This function reads a tab-separated file produced by the 'featureCounts' function
and removes the unnecessary columns.
}
