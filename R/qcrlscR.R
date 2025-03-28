## Functions for QC-RLSC
## wl-27-03-2025, Thu: Prepare for devtools::check()

## -----------------------------------------------------------------------
#' Wrapper function for QC-RLSC
#'
#' Wrapper function for QC-RLSC
#'
#' @param dat A data frame with samples (row) and variables (column).
#' @param cls.qc A vector with string of "qc" and "sample".
#' @param cls.bl A vector with string of batch indicators.
#' @param method Data scaling method. Support "subtract" and "divide"
#' @param intra A logical value indicating whether signal correction is
#'   performed inside each batch.
#' @param opti A logical value indicating whether or not 'span' parameters
#'   are optimised.
#' @param log10 A logical value indicating whether log10 transformation for
#'   the data set or not.
#' @param outl A logical value indicating whether or not QC outlier
#'   detection is employed.
#' @param shift A logical value indicating whether or not batch shift is
#' applied after signal correction.
#' @param ... Other parameter for 'loess'.
#' @return  A corrected data frame.
#' @export
## wl-19-07-2024, Fri: wrapper function for QC-RLSC
## wl-27-03-2025, Thu: remove 'tidyverse', including
##  'res <- dplyr::bind_rows(res)'
qc_rlsc_wrap <- function(dat, cls.qc, cls.bl,
                         method = c("subtract", "divide"),
                         intra = FALSE, opti = TRUE, log10 = TRUE,
                         outl = TRUE, shift = TRUE, ...) {

  ## log transformation
  if (log10) {
    dat[dat ==0] <- NA
    dat <- log10(dat)
  }

  ## QC outlier detection
  if  (outl) {
    dat <- sapply(dat, function(x){ #' x <- dat[, 6, drop = T]
      qc_ind <- grepl("qc", cls.qc, ignore.case =  TRUE, perl = TRUE)
      ## get median of qc data
      qc_dat <- x[qc_ind]
      qc_median <- median(qc_dat, na.rm = TRUE)
      ## assign other data as NA for QC outlier detection
      tmp <- x
      tmp[!qc_ind] <- NA
      ## QC outlier detection
      out_ind <- outl_det_u(tmp)
      ## assign outlier as qc median
      x[out_ind] <- qc_median
      return(x)
    })
    dat <- as.data.frame(dat)
  }

  ## QC-RLSC
  if (!intra) {
    res <- qc_rlsc(dat, cls.qc, method = method, opti = opti)
  } else { # wl-16-08-2023, Wed: do signal correction inside each batch
    res <- lapply(levels(cls.bl), function(x){
      idx <- cls.bl %in% x
      tmp <- qc_rlsc(dat[idx,], cls.qc[idx], method = method, opti = opti)
    })
    res <- do.call(rbind, res)
  }

  ## batch shift. sensitive to missing values
  if (shift) {
    res <- batch_shift(res, cls.bl, overall_average = T)
  }

  ## inverse log10 transformation
  if (log10) {
    res <- 10^res
  }

  return(res)
}

## -----------------------------------------------------------------------
#' QC based robust LOESS signal correction (QC-RLSC)
#'
#' QC based robust LOESS (locally estimated scatterplot smoothing) signal
#' correction (QC-RLSC)
#'
#' @param x A data frame with samples (row) and variables (column).
#' @param y A vector with string of "qc" and "sample".
#' @param method Data scaling method.
#' @param opti A logical value indicating whether or not optimise 'span'
#' @param ... Other parameter for 'loess'.
#' @return  A corrected data frame.
#' @references
#'  Dunn et al. Procedures for large-scale metabolic profiling of serum and
#'  plasma using gas chromatography and liquid chromatography coupled to
#'  mass spectrometry. Nature Protocols 6, 1060–1083 (2011)
#' @export
## wl-14-08-2023, Mon: QC-RLSC
##   Note that the variables are divided by predicted values using qc-based
##   'loess'. Have nothing to do with batch information. Also see:
##     1.) "statTarget"(https://bit.ly/454KIiT)
##     2.) "Rcpm" (https://github.com/ricoderks/Rcpm/)
##   see 'staTarget_shiftcor_v2.R'  and 'Rcpm_qc_rlsc.R'
## wl-20-06-2024, Thu: Add 'method'. The default method is 'subtract'
##   which gets the data back on the original scale.
## wl-02-07-2024, Tue: 'losses' fails in vector which includes only missing
##   values. Should filter based on each batch or lower down threshold
## wl-08-07-2024, Mon: call 'loess_gcv' for optimisation span
## wl-30-07-2024, Tue: use less.control for extrapolation
qc_rlsc <- function(x, y, method = c("subtract", "divide"), opti = TRUE,
                    ...) {
  method <- match.arg(method)

  ## order for qc
  ind <- grep("qc", y, ignore.case =  TRUE, perl = TRUE)
  if (length(ind) == 0) stop("No QC samples")

  ## order for interpolation (all samples and qcs)
  ord <- 1:length(y)

  ## wl-02-07-2024, Tue: problem if x[ind, i] are all missing values.
  nc <- ncol(x)
  smo <- sapply(1:nc, function(i) { # i = 60 #' cat(i, "\n")
    ## apply loess to QCs
    if (opti) { # optimise span
      loe <- loess_gcv(ind, x[ind, i, drop = TRUE],
                       span.range = c(.05, .95), ...)
    } else {    # default span: 0.75
      loe <- loess(x[ind, i, drop = TRUE] ~ ind,
                   control = loess.control(surface = "direct"), ...)
    }
    if (T) {              # predict all (sample and qc)
      yf <- predict(loe, ord)
    } else {              # approximate all the samples (interpolation only)
      yf <- approx(x = ind, y = loe$fitted, xout = ord)$y
    }
  })

  ## get average values of QC
  mn <- colMeans(x[ind, , drop = F], na.rm = TRUE)
  mn <- matrix(mn, nrow = nrow(x), ncol = ncol(x), byrow = TRUE)

  ## normalise data set
  res <-
    switch(method,
      "subtract" = (x - smo) + mn,
      "divide"   = (x / smo) * mn
    )
  res <- as.data.frame(res)

  return(res)
}

## -----------------------------------------------------------------------
#' Optimise LOESS's span parameter using generalized cross-validation (GCV)
#'
#' Optimise LOESS's span parameter using generalized cross-validation (GCV)
#'
#' @param x, y vectors giving the coordinates of the points in the scatter
#'   plot.
#' @param span.range a vector for span optimisation range.
#' @param ... Other parameter for 'loess'.
#' @return  A optimised loess model.
#' @importFrom stats optimize resid update
#' @noRd
#' @keywords internal
## wl-08-07-2024, Mon: Optimal loess model with GCV
## Modified from https://bit.ly/3zBo3Qn
## wl-30-07-2024, Tue: use less.control for extrapolation
loess_gcv <- function(x, y, span.range = c(.05, .95), ...) {

  ## ---------------------------------------------------------------------
  ## GCV for loess.
  ## Modified from code by Michael Friendly
  ## http://tolstoy.newcastle.edu.au/R/help/05/11/15899.html
  loessGCV <- function(x) {
    if (!(inherits(x, "loess")))
      stop("Error: argument must be a loess object")
    span <- x$pars$span
    n <- x$n
    traceL <- x$trace.hat
    sigma2 <- sum(resid(x)^2) / (n - 1)
    gcv <- n * sigma2 / (n - traceL)^2
    result <- list(span = span, gcv = gcv)
    result
  }

  ## ---------------------------------------------------------------------
  ## Wrapper function to link optimize() with loessGCV()
  bestLoess <- function(model, spans = c(.05, .95)) {
    f <- function(span) {
      mod <- update(model, span = span)
      loessGCV(mod)[["gcv"]]
    }
    result <- optimize(f, spans)
    result
  }

  ## ---------------------------------------------------------------------
  ## default span = 0.75
  loe <- loess(y ~ x,
               control = loess.control(surface = "direct"), ...)

  ## optimise model span
  opti <- bestLoess(loe, spans = span.range)

  ## optimal model
  loe_gcv <- update(loe, span = opti$minimum)

  return(loe_gcv)
}

## ------------------------------------------------------------------------
#' Batch shifting
#'
#' Remove batch effect withing each block.
#'
#' @param x a data matrix.
#' @param y a categorical data for batch/block information.
#' @param method method for shifting.
#' @param overall_average a logical value to indicate whether or not an
#'   overall average will be added after shifting.
#' @return  a shifted data matrix.
#' @references
#'   Silvia Wagner, et.al, Tools in Metabonomics: An Integrated Validation
#'   Approach for LC-MS Metabolic Profiling of Mercapturic Acids in Human
#'   Urine Anal. Chem., 2007, 79 (7), pp 2918-2926, DOI: 10.1021/ac062153w
#' @export
## wl-07-07-2011, Thu: Batch shifting: remove mean within each batch/block
## wl-03-07-2024, Wed: Minor changes
##  - Very sensitive with missing values.
##  - Shift to overall average
## wl-18-07-2024, Thu: fix a bug
batch_shift <- function(x, y, method = "mean", overall_average = TRUE) {
  x <- as.data.frame(x)
  ## overall
  o.mean <- sapply(x, method, na.rm = T)
  o.mean <- matrix(o.mean, nrow=nrow(x), ncol=ncol(x), byrow=TRUE)
  ## group
  g.mean <- sapply(x, function(x) tapply(x, y, method, na.rm = T))
  g.mean <- sapply(1:ncol(x), function(i) g.mean[, i][y])
  ## overall average
  if (overall_average) {
    x <- x - g.mean + o.mean
  } else {
    x <- x - g.mean
  }
  return(x)
}

## ------------------------------------------------------------------------
#' Univariate outlier detection
#'
#' Perform outlier detection using univariate method.
#'
#' @param x a numeric vector.
#' @param method method for univariate outlier detection. Only `percentile`
#'   and `median` are supported.
#' @return  a logical vector.
#' @references
#'   Wilcox R R, Fundamentals of Modern Statistical Methods: Substantially
#'   Improving Power and Accuracy, Springer 2010 (2nd edition), pages 31-35.
#' @details
#'  - `median`: the absolute difference between the observation and the sample
#'    median is larger than 2 times of the Median Absolute Deviation divided
#'    by 0.6745.
#'  - `percentile`: either smaller than the 1st quartile minus 1.5 times of
#'    IQR, or larger than the 3rd quartile plus 1.5 times of IQR.
#' @examples
#' x <- c(2, 3, 4, 5, 6, 7, NA, 9, 50, 50)
#' outl_det_u(x, "percentile")
#' @export
## wl-19-09-2020, Sat: Univariate outlier detection.
##   Modified from R package GmAMisc.
## wl-17-07-2024, Wed: add 'na.rm' for missing values
outl_det_u <- function(x, method = c("percentile", "median")) {
  method <- match.arg(method)
  if (method == "median") {
    med <- median(x, na.rm = TRUE)
    mad <- median(abs(med - x), na.rm = TRUE)
    outlier <- abs(x - med) > 2 * (mad / 0.6745)
  }
  if (method == "percentile") {
    q1 <- quantile(x, 0.25, na.rm = TRUE)
    q3 <- quantile(x, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    outlier <- x < q1 - 1.5 * iqr | x > q3 + 1.5 * iqr
  }
  outlier[is.na(outlier)] <- FALSE
  return(outlier)
}

## ------------------------------------------------------------------------
#' Missing value percentage
#'
#' Calculate missing value percentage.
#'
#' @param x an vector, matrix or data frame.
#' @return missing value percentage.
#' @export
## wl-24-11-2021, Wed: extract from 'mv.stats' in 'mt'.
mv_perc <- function(x) {
  if (is.matrix(x)) {
    apply(x, 2, mv_perc)
  } else if (is.vector(x)) {
    round(sum(is.na(x) | is.nan(x)) / length(x), digits = 3)
  } else if (is.data.frame(x)) {
    sapply(x, mv_perc)
  } else {
    round(sum(is.na(as.vector(x)) | is.nan(as.vector(x))) /
          length(as.vector(x)), digits = 3)
  }
}

## -----------------------------------------------------------------------
#' Filtering variable based on the percentage of missing values
#'
#' This function calculates the percentage of missing values and keeps those
#' features with missing values percentage less than the designed threshold.
#'
#' @param x a data matrix. The columns are features.
#' @param thres threshold of missing values. Features less than this
#'   threshold will be kept. Value has to be between 0 and 1.
#' @return a list of with contents: \itemize{
#'  \item dat the filtered data matrix
#'  \item idx a logical vector of index for keeping features.
#' }
#' @export
## wl-14-06-2011, Tue: Filter features based on the percentage of MVs
## wl-17-06-2021, Thu: several version but this one is simple. Need to test
## wl-06-11-2018, Tue: feature filter index based on missing values
mv_filter <- function(x, thres = 0.3) {

  if (!(is.matrix(x) || is.data.frame(x))) {
    stop("\n'data' is not matrix or data frame.\n")
  }

  thres <- thres * nrow(x)

  ## get number of Na, NaN and zero in each of feature/variable
  count <- apply(x, 2, function(y) {
    tmp <- sum(is.na(y) | is.nan(y) | (y == 0))
  })

  ## index of features whose number of MVs are less than threshold
  idx <- count <= thres

  x <- x[, idx, drop = FALSE]
  return(list(dat = x, idx = idx))
}

## -----------------------------------------------------------------------
#' Data filtering based on "qc" missing values
#'
#' @param x a data matrix.
#' @param y a character string with contents of "sample", "qc" and "blank".
#' @param thres threshold of missing values. Features less than this
#'  threshold will be kept.
#' @return a list of with contents: \itemize{
#'  \item dat the filtered data matrix
#'  \item idx a logical vector of index for keeping features.
#' }
#' @export
## wl-14-06-2011, Tue: Filter features based on missing values in QC
## wl-01-07-2024, Mon: Review and minor change. different from 'qc_filter'
##   in 'mtExtra'.
mv_filter_qc <- function(x, y, thres = 0.3) {
  tmp <- grep("qc", y, ignore.case =  TRUE, perl = TRUE)
  idx <- mv_filter(x[tmp, , drop = FALSE], thres = thres)$idx
  x <- x[, idx, drop = FALSE]
  return(list(dat = x, idx = idx))
}

## ------------------------------------------------------------------------
#' @description `qcrlscR` implements quality control–based robust LOESS
#' signal correction for metabolomics data analysis.
#'
#' @section Main functions:
#' The `qcrlscR` provides functions for signal correction for metabolomics
#' data analysis. It allows users to optimise `span` for LOESS. This package
#' also provides simple functions for univariate outlier detection and
#' missing value filtering. Users need to use other R packages for missing
#' value filling before applying this package if the data set has missing
#' values.
#'
#' @section Package context:
#' This package does not use "tidyverse" and only uses basic R for
#' simplicity and easy maintenance.
#'
#' @importFrom stats median approx loess loess.control predict quantile
#' @keywords internal
"_PACKAGE"

##  1) qc_rlsc_wrap
##  2) qc_rlsc
##  3) loess_gcv
##  4) batch_shift
##  5) outl_det_u
##  6) mv_perc
##  7) mv_filter
##  8) mv_filter_qc
