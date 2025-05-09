% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qcrlscR.R
\name{batch.shift}
\alias{batch.shift}
\title{Batch shifting}
\usage{
batch.shift(x, y, method = "mean", overall_average = TRUE)
}
\arguments{
\item{x}{a data matrix.}

\item{y}{a categorical data for batch/block information.}

\item{method}{method for shifting.}

\item{overall_average}{a logical value to indicate whether or not an
overall average will be added after shifting.}
}
\value{
a shifted data matrix.
}
\description{
Remove batch effect withing each block.
}
\examples{
names(man_qc)
data <- man_qc$data
meta <- man_qc$meta
## batch shifting
cls.bl <- factor(meta$batch)
res <- batch.shift(data, cls.bl, overall_average = TRUE)
}
\references{
Silvia Wagner, et.al, Tools in Metabonomics: An Integrated Validation
Approach for LC-MS Metabolic Profiling of Mercapturic Acids in Human
Urine Anal. Chem., 2007, 79 (7), pp 2918-2926, DOI: 10.1021/ac062153w
}
