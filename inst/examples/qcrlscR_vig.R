#' ---
#' title: "Quality Control–based Robust LOESS Signal Correction in R"
#' author: "Wanchang Lin"
#' date: "`r Sys.Date()`"
#' output:
#'   BiocStyle::html_document:
#'     toc_depth: 3
#'     number_section: false
#'     toc_float: false
#'   BiocStyle::pdf_document:
#'     keep_tex: true
#'     toc: true
#'     toc_depth: 3
#'     number_section: false
#'     citation_package: natbib
#'     latex_engine: xelatex
#' always_allow_html: true
#' geometry: margin=1in
#' fontsize: 11pt
#' ---

#' <!--
#' # Test code for QC-RLSC
#' > wl-28-03-2025, Fri: Rscript -e 'knitr::spin("qcrlscR_vig.R")'
#' > wl-31-03-2025, Mon: get html and pdf:
#' >   Rscript -e "rmarkdown::render('qcrlscR_vig.R', BiocStyle::html_document())"
#' >   Rscript -e "rmarkdown::render('qcrlscR_vig.R', BiocStyle::pdf_document())"
#' -->

#+ common, include=F
rm(list = ls(all = TRUE))
options(help_type = "html")
library(knitr)
opts_chunk$set(
  collapse = TRUE,
  cache = TRUE,
  comment = "#>",
  messages = FALSE,
  warning = FALSE,
  tidy = FALSE,
  fig.align = "center",
  fig.width = 10,
  fig.height = 10,
  # dev = "png",
  # dpi = 100,
  # fig.margin = TRUE,
  # fig.asp = 0.618,  # 1 / phi
  # fig.keep = "none",
  # fig.path = "figure",
  fig.show = "hold"
)

## ---- Load libraries ----
#' ## Load libraries

#' R package `mt` is used to plot PCA, PLS and LDA plots to assess
#' performance of signal correction implemented in R package `qcrlscR`. R
#' package `tityverse` fulfils some data crunch operations and `tictoc`
#' records the running time, especially for the optimisation of LOESS. All
#' of these packages are available in the CRAN package repository.  

#+ message=F
pkgs <- c("qcrlscR", "mt", "tidyverse", "tictoc")
## install.packages(pkgs)
invisible(lapply(pkgs, library, character.only = TRUE))

## ---- Read data ----
#' ## Read data

#' The data `man_qc` in package `qcrlscR` is a list of two data frames,
#' `data` and `meta`. 
names(man_qc)
t(sapply(man_qc, dim))

#' Get meta and data matrix:
meta <- man_qc$meta
data <- man_qc$data %>%
  mutate_if(is.character, as.numeric)

#' Extract group information of batch and sample types from meta:
names(meta)
cls.qc <- factor(meta$sample_type)
table(cls.qc)

cls.bl <- factor(meta$batch)
table(cls.bl)

## ---- Missing value filter ----
#' ## Missing value filter

#' Before signal correction, the data should be checked based on the missing
#' values rate and filtered if the rate is higher than the threshold, such
#' as 20%.
#'
#' Check missing value rates:
tail(sort(mv.perc(data)), 20)

#' Filter data matrix based on missing values rate:
filter_qc <- FALSE      # filter on qc missing values or all missing values
thres <- 0.15           # threshold for filtering

if (filter_qc) {      # filter using all missing values
  ret <- mv.filter(data, thres = thres)
} else {              # filter using qc missing values
  ret <- mv.filter.qc(data, cls.qc, thres = thres)
}

#' Update data matrix after filtering:
dat <- ret$dat

#' Missing value imputation is not required in `qcrlscR` but visualisation
#' of a data matrix, like PCA and LDA plots, does not allow the missing
#' values. Here the missing value filling is used to data screening before
#' and after signal correction.
#'
#' `mv.fill` in R package `mt` is used here for missing value imputation.
#' It should be noted that there are a lot of R package available for such
#' purpose in CRAN repository.
dat_fill  <- dat %>% mv.fill(method = "median", ze_ne = T) %>% as_tibble()

#' Two categories methods, unsupervised method, PCA, and supervised methods,
#' PLS and PCA-LDA, are used for data screening before and after signal
#' correction.
#'
#' PCA plot for sample types:
pcaplot(dat_fill, cls.qc, pcs = c(2, 1), ep = 1)

#' PCA plot for batches:
pcaplot(dat_fill, cls.bl, pcs = c(2, 1), ep = 1)

#' LDA plot for batches:
plot(pcalda(dat_fill, cls.bl))

#' LDA plot of batches: LD1 vs LD2 (only for batch groups larger than 2)
plot(pcalda(dat_fill, cls.bl), dimen = c(1:2), ep = 2)

#' PLS plot of batches: LC1 vs LC2
plot(plslda(dat_fill, cls.bl), dimen = c(1:2), ep = 2)

## ---- Set parameters for QC-RLSC ----
#' ## Set parameters for QC-RLSC

#' Some parameters for signal correction:
log10 <- T            # log 10 transform data or not
outl <- T             # outlier detect in qc samples or not
intra <- F            # signal correction within batch or not
method <- "subtract"  # two methods: "subtract", "divide"
opti <- T             # optimise smooth parameter  or not
shift <- T            # batch shift or not

## ---- Logarithmic transformation ----
#' ## Logarithmic transformation

#' Log transformation or not:
if (log10) {
  dat[dat == 0] <- NA
  dat <- log10(dat)
}

## ---- QC outlier detection ----
#' ## QC outlier detection

#' Outlier detection based on QC or not:
if (outl) {
  dat <- sapply(dat, function(x) { #' x <- dat[, 6, drop = T]
    qc_ind <- grepl("qc", cls.qc, ignore.case =  TRUE, perl = TRUE)
    ## get median of qc data
    qc_dat <- x[qc_ind]
    qc_median <- median(qc_dat, na.rm = TRUE)
    ## assign other data as NA for QC outlier detection
    tmp <- x
    tmp[!qc_ind] <- NA
    ## QC outlier detection
    out_ind <- outl.det.u(tmp)
    ## assign outlier as qc median
    x[out_ind] <- qc_median
    return(x)
  }) %>% as_tibble()
}
dat

#' User can outlier detection methods provided by other R packages.
#' Here `qcrlscR` only implements a simple univariate method, and detected
#' outliers are replaced with the median values of QC data points.

## ---- QC-RLSC ----
#' ## QC-RLSC

#' Perform qc-rlsc within each batch or not (intra batch or inter batch):
tic()
if (!intra) {
  res <- qc.rlsc(dat, cls.qc, method = method, opti = opti)
} else { # do signal correction inside each batch
  res <- lapply(levels(cls.bl), function(x) {
    idx <- cls.bl %in% x
    tmp <- qc.rlsc(dat[idx,], cls.qc[idx], method = method, opti = opti)
  })
  res <- bind_rows(res)
}
toc()

#' `qcrlscR` can optimise smoothing span in a range of 0.05 and 0.95,
#' controlled by a binary option `opti` in function `qc.rlsc` and
#' `qc.rlsc.wrap`.
#'
#' Data visualisation after signal correction:
res_fill  <- res %>% mv.fill(method = "median", ze_ne = T) %>% as_tibble()

#' PCA plot for sample types:
pcaplot(res_fill, cls.qc, pcs = c(2, 1), ep = 1)

#' PCA plot for batches:
pcaplot(res_fill, cls.bl, pcs = c(2, 1), ep = 1)

#' LDA plot for batches:
plot(pcalda(res_fill, cls.bl))

#' LDA plot of batches: LD1 vs LD2 (only for batch groups larger than 2)
plot(pcalda(res_fill, cls.bl), dimen = c(1:2), ep = 2)

#' PLS plot of batches: LC1 vs LC2 
plot(plslda(res_fill, cls.bl), dimen = c(1:2), ep = 2)

## ---- Batch shift ----
#' ## Batch shift

#' If the batch effects are still in the data set, a straightforward batch
#' shifting method can be applied:
if (shift) {
  res <- batch.shift(res, cls.bl, overall_average = T) %>% as_tibble()
}

#' Data visualisation after batch shift:
res_fill  <- res %>% mv.fill(method = "median", ze_ne = T) %>% as_tibble()

#' PCA plot for sample types:
pcaplot(res_fill, cls.qc, pcs = c(2, 1), ep = 1)

#' PCA plot for batches:
pcaplot(res_fill, cls.bl, pcs = c(2, 1), ep = 1)

#' LDA plot for batches:
plot(pcalda(res_fill, cls.bl))

#' LDA plot of batches: LD1 vs LD2 (only for batch groups larger than 2)
plot(pcalda(res_fill, cls.bl), dimen = c(1:2), ep = 2)

#' PLS plot of batches: LC1 vs LC2
plot(plslda(res_fill, cls.bl), dimen = c(1:2), ep = 2)

## ---- Save results ----
#' ## Save results

#' Inverse log10 transformation:
if (log10) {
  res <- 10^res %>% as_tibble()
}

#+ qcrlsc_save, eval=F, include=T
tmp <- list(data =  res, meta = meta)
write.xlsx(tmp, file = here::here("data", paste0(FILE, "_res.xlsx")),
           asTable = F, overwrite = T, rowNames = F, colNames = T)

## ---- QC-RLSC wrapper function ----
#' ## QC-RLSC wrapper function

#' User can use wrapper function `qc.rlsc.wrap` directly with options of
#' intra or inter batch signal correction, optimisation of span, log
#' transformation and batch shifting.
#+ qcrlsc_save_1, eval=F, include=T
res <- qc.rlsc.wrap(dat, cls.qc, cls.bl, method, intra, opti, log10, outl,
                    shift) tmp <- list(data =  res, meta = meta)
write.xlsx(tmp, file = here::here("data", paste0(FILE, "_res.xlsx")),
           asTable = F, overwrite = T, rowNames = F, colNames = T)
