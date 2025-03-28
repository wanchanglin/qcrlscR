% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qcrlscR.R
\name{mv_filter_qc}
\alias{mv_filter_qc}
\title{Data filtering based on "qc" missing values}
\usage{
mv_filter_qc(x, y, thres = 0.3)
}
\arguments{
\item{x}{a data matrix.}

\item{y}{a character string with contents of "sample", "qc" and "blank".}

\item{thres}{threshold of missing values. Features less than this
threshold will be kept.}
}
\value{
a list of with contents: \itemize{
\item dat the filtered data matrix
\item idx a logical vector of index for keeping features.
}
}
\description{
Data filtering based on "qc" missing values
}
