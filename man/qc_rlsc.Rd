% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qcrlscR.R
\name{qc_rlsc}
\alias{qc_rlsc}
\title{QC based robust LOESS signal correction (QC-RLSC)}
\usage{
qc_rlsc(x, y, method = c("subtract", "divide"), opti = TRUE, ...)
}
\arguments{
\item{x}{A data frame with samples (row) and variables (column).}

\item{y}{A vector with string of "qc" and "sample".}

\item{method}{Data scaling method.}

\item{opti}{A logical value indicating whether or not optimise 'span'}

\item{...}{Other parameter for 'loess'.}
}
\value{
A corrected data frame.
}
\description{
QC based robust LOESS (locally estimated scatterplot smoothing) signal
correction (QC-RLSC)
}
\references{
Dunn et al. Procedures for large-scale metabolic profiling of serum and
plasma using gas chromatography and liquid chromatography coupled to
mass spectrometry. Nature Protocols 6, 1060–1083 (2011)
}
