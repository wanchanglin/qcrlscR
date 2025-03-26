#' ---
#' title: "Quality control–based robust LOESS signal correction (QC-RLSC)"
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
#' > wl-20-06-2024, Thu: commence <br/>
#' > wl-21-06-2024, Fri: test qc_filter and qc_pass <br/>
#' > wl-03-07-2024, Wed: more test.
#' > wl-17-07-2024, Wed: QC outlier detection
#' > wl-18-07-2024, Thu: batch shift
#' > wl-23-07-2024, Tue: Rscript -e 'knitr::spin("qc_rlsc_vig.R")'
#' -->

#+ common, include=F
rm(list = ls(all = TRUE))
options(help_type = "html")
source("_common.R")
pkgs <- c("mt", "mtExtra", "tidyverse", "readxl", "openxlsx", "tictoc")
invisible(lapply(pkgs, library, character.only = TRUE))

## ---- Read data ----
#' ## Read data

#' Select file for signal correction
## FILE <- "data_qcmxp_b4_tidy"
## FILE <- "data_qcmxp_b3_tidy"
## FILE <- "data_qcmxp_b2_tidy"
FILE <- "data_man_b4_tidy"

PATH <- here::here("data", paste0(FILE, ".xlsx"))

#' Load into R
xls  <- PATH %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_excel, path = PATH)

#' Check the data
names(xls)
t(sapply(xls, dim))

#' Get meta and data matrix
meta <- xls$meta
data <- xls$data %>%
  mutate_if(is.character, as.numeric)

peak <- xls$peak

#' Extract group information of batch and sample types
names(meta)
if (T) {     # meta with sample_type and batch
  (cls.qc <- factor(meta$sample_type))
  (cls.bl <- factor(meta$batch))
} else {     # meta with SampleType and Batch
  (cls.qc <- factor(meta$SampleType))
  (cls.bl <- factor(meta$Batch))
}

## ---- Missing value filter and fill ----
#' ## Missing value filter and fill

#' Let zero as NA before missing value process
data[data == 0] <- NA

#' Check missing value rates
tail(sort(mv_perc(data)), 20)

#' Filter based on missing values
filter_qc <- FALSE    # filter on qc missing values or all missing values
thres <- 0.2           # threshold for filtering

if (filter_qc) {      # filter using all missing values
  ret <- mv_filter(data, thres = thres)
} else {              # filter using qc missing values
  ret <- mv_filter_qc(data, cls.qc, thres = thres)
}

#' Update data matrix and peak
dat <- ret$dat
pek <- peak[ret$idx, ]

#' Missing values filling for visualisation
dat_fill  <- dat %>% mv.fill(method = "median", ze_ne = T) %>% as_tibble()

#' Data screening before signal correction

#' PCA plot for sample types
pcaplot(dat_fill, cls.qc, pcs = c(2, 1), ep = 1)

#' PCA plot for batches
pcaplot(dat_fill, cls.bl, pcs = c(2, 1), ep = 1)

#' LDA plot for batches
plot(pcalda(dat_fill, cls.bl))

#' LDA plot of batches: LD1 vs LD2 (only for batch groups larger than 2)
plot(pcalda(dat_fill, cls.bl), dimen = c(1:2), ep = 2)

#' PLS plot of batches: LC1 vs LC2
plot(plslda(dat_fill, cls.bl), dimen = c(1:2), ep = 2)

## ---- Set parameters for QC-RLSC ----
#' ## Set parameters for QC-RLSC

method <- "subtract"  # two methods: "subtract", "divide"
intra <- F            # signal correction within batch or not
opti <- T             # optimise smooth parameter  or not
log10 <- T            # log 10 transform data or not
outl <- T             # outlier detect in qc samples or not
shift <- T            # batch shift or not

## ---- QC outlier detection ----
#' ## QC outlier detection

#' log transformation
if (log10) {
  dat[dat == 0] <- NA
  dat <- log10(dat)
}

#' outlier detection based on QC
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
    out_ind <- outl_det_u(tmp)
    ## assign outlier as qc median
    x[out_ind] <- qc_median
    return(x)
  }) %>% as_tibble()
}
dat

## ---- QC-RLSC ----
#' ## QC-RLSC

#' perform qc-rlsc within each batch or not
tic()
if (!intra) {
  res <- qc_rlsc(dat, cls.qc, method = method, opti = opti)
} else { # do signal correction inside each batch
  res <- lapply(levels(cls.bl), function(x) {
    idx <- cls.bl %in% x
    tmp <- qc_rlsc(dat[idx,], cls.qc[idx], method = method, opti = opti)
  })
  res <- bind_rows(res)
}
toc()

#' Data visualisation after signal correction
res_fill  <- res %>% mv.fill(method = "median", ze_ne = T) %>% as_tibble()

#' PCA plot for sample types
pcaplot(res_fill, cls.qc, pcs = c(2, 1), ep = 1)

#' PCA plot for batches
pcaplot(res_fill, cls.bl, pcs = c(2, 1), ep = 1)

#' LDA plot for batches
plot(pcalda(res_fill, cls.bl))

#' LDA plot of batches: LD1 vs LD2 (only for batch groups larger than 2)
plot(pcalda(res_fill, cls.bl), dimen = c(1:2), ep = 2)

#' PLS plot of batches: LC1 vs LC2 
plot(plslda(res_fill, cls.bl), dimen = c(1:2), ep = 2)

## ---- Batch shift ----
#' ## Batch shift

if (shift) {
  res <- batch_shift(res, cls.bl, overall_average = T) %>% as_tibble()
}

#' Data visualisation after batch shift
res_fill  <- res %>% mv.fill(method = "median", ze_ne = T) %>% as_tibble()

#' PCA plot for sample types
pcaplot(res_fill, cls.qc, pcs = c(2, 1), ep = 1)

#' PCA plot for batches
pcaplot(res_fill, cls.bl, pcs = c(2, 1), ep = 1)

#' LDA plot for batches
plot(pcalda(res_fill, cls.bl))

#' LDA plot of batches: LD1 vs LD2 (only for batch groups larger than 2)
plot(pcalda(res_fill, cls.bl), dimen = c(1:2), ep = 2)

#' PLS plot of batches: LC1 vs LC2
plot(plslda(res_fill, cls.bl), dimen = c(1:2), ep = 2)

## ---- Save results ----
#' ## Save results

#' inverse log10 transformation
res <- 10^res %>% as_tibble()

## tmp <- list(data =  res, meta = meta)
tmp <- list(data =  res, meta = meta, peak = pek)

## write.xlsx(tmp, file = here::here("data", paste0(FILE, "_res.xlsx")),
##            asTable = F, overwrite = T, rowNames = F, colNames = T)
