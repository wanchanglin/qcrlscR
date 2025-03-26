#' # Test code for QC-RLSC
#' > wl-20-06-2024, Thu: commence <br/>
#' > wl-21-06-2024, Fri: test qc_filter and qc_pass <br/>
#' > wl-03-07-2024, Wed: more test.
#' > wl-17-07-2024, Wed: QC outlier detection
#' > wl-18-07-2024, Thu: batch shift
#+ common, include=F
rm(list = ls(all = TRUE))
options(help_type = "html")
#' source("_common.R")
pkgs <- c("mt", "tidyverse", "readxl", "openxlsx", "tictoc")
invisible(lapply(pkgs, library, character.only = TRUE))
source("qcrlsc.R")

## 1.) ----- Read data -----
#' FILE <- "data_qcmxp_b2_tidy"
#' FILE <- "data_qcmxp_b3_tidy"
#' FILE <- "data_qcmxp_b4_tidy"
#' FILE <- "data_man_b4_tidy"
FILE <- "MTBLS79_tidy"

DIR <- "./data/"
(PATH <- paste(DIR, FILE, ".xlsx", sep = ""))

xls  <- PATH %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_excel, path = PATH)

#' Check the data
names(xls)
t(sapply(xls, dim))
data <- xls$data
meta <- xls$meta
peak <- xls$peak

#' group information for data screening
names(meta)

#' need to select right column names for batch and sample type
if (F) {   # for qcmxp convert format
  (cls.qc <- factor(meta$SampleType))
  (cls.bl <- factor(meta$Batch))
} else {  # for qcrlsc format
  (cls.qc <- factor(meta$sample_type))
  (cls.bl <- factor(meta$batch))
}

## 2.) ----- Missing value filtering and filling -----

#' let zero as NA before missing value process
data[data == 0] <- NA

#' check missing value rates
tail(sort(mv_perc(data)), 20)

#' filter based on missing values
if (F) {      # filter using all missing values
  ret <- mv_filter(data, thres = 0.2)
} else {      # filter using qc missing values
  ret <- mv_filter_qc(data, cls.qc, thres = 0.2)
}
dat <- ret$dat
pek <- peak[ret$idx,]

dat_fill  <- dat %>% mv.fill(method = "median", ze_ne = T) %>% as_tibble()

#' Data screening before batch alignment
pcaplot(dat_fill, cls.qc, pcs = c(2, 1), ep = 1)
pcaplot(dat_fill, cls.bl, pcs = c(2, 1), ep = 1)
plot(pcalda(dat_fill, cls.bl))
plot(pcalda(dat_fill, cls.bl), dimen = c(1:2), ep = 2)
plot(plslda(dat_fill, cls.bl), dimen = c(1:2), ep = 2)

#' Parameters for QC-RLSC
method <- "divide"  # "subtract" 
intra <- F
opti <- T
log10 <- T
outl <- T
shift <- T

## 3.) ----- QC outlier detection  -----

#' log transformation?
if (log10) {
  dat[dat == 0] <- NA
  dat <- log10(dat)
}

if (outl) {
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
    ## asisgn outlier as qc median
    x[out_ind] <- qc_median
    return(x)
  }) %>% as_tibble()
}
dat

## 4.) ----- QC-RLSC  -----
tic()
if (intra) {            # signal correction inside each batch
  res <- lapply(levels(cls.bl), function(x) {
    idx <- cls.bl %in% x
    tmp <- qc_rlsc(dat[idx,], cls.qc[idx], method = method, opti = opti)
  })
  res <- bind_rows(res)
} else {
  res <- qc_rlsc(dat, cls.qc, method = method, opti = opti)
}
toc()

#' Data screening after signal correction
res_fill  <- res %>% mv.fill(method = "median", ze_ne = T) %>% as_tibble()

pcaplot(res_fill, cls.qc, pcs = c(2, 1), ep = 1)
pcaplot(res_fill, cls.bl, pcs = c(2, 1), ep = 1)
plot(pcalda(res_fill, cls.bl))
plot(pcalda(res_fill, cls.bl), dimen = c(1:2), ep = 2)
plot(plslda(res_fill, cls.bl), dimen = c(1:2), ep = 2)

## 5.) ----- Batch shift  -----
#' batch shift is sensitive to missing values
if (shift) {
  res <- batch_shift(res, cls.bl, overall_average = T) %>% as_tibble()
}

#' Data screening after batch shift
res_fill  <- res %>% mv.fill(method = "median", ze_ne = T) %>% as_tibble()

pcaplot(res_fill, cls.qc, pcs = c(2, 1), ep = 1)
pcaplot(res_fill, cls.bl, pcs = c(2, 1), ep = 1)
plot(pcalda(res_fill, cls.bl))
plot(pcalda(res_fill, cls.bl), dimen = c(1:2), ep = 2)
plot(plslda(res_fill, cls.bl), dimen = c(1:2), ep = 2)

#' inverse log10 transformation
tres <- 10^res %>% as_tibble()

## 6.) ----- Save results -----
tmp <- list(data = res, meta = meta, peak = pek)
write.xlsx(tmp, file = paste0("./res/", FILE, "_res.xlsx"),
           asTable = F, overwrite = T, rowNames = F, colNames = T)
