#' # Test code for QC-RLSC
#' > wl-19-07-2024, Fri: commence <br/>
#+ common, include=F
rm(list = ls(all = TRUE))
options(help_type = "html")
# source("_common.R")
pkgs <- c("mt", "tidyverse", "readxl", "openxlsx", "tictoc")
invisible(lapply(pkgs, library, character.only = TRUE))
source("qcrlsc.R")

## 1.) ----- Data load -----
#' FILE <- "data_qcmxp_b2_tidy"
#' FILE <- "data_qcmxp_b3_tidy"
FILE <- "data_qcmxp_b4_tidy"
#' FILE <- "data_man_b4_tidy"
#' FILE <- "MTBLS79_tidy"

#' DIR <- "./data_large/"
DIR <- "./data/"
(PATH <- paste(DIR, FILE, ".xlsx", sep = ""))

xls  <- PATH %>%
  excel_sheets() %>%
  set_names() %>%
  map(read_excel, path = PATH)
#' Check the data
names(xls)
t(sapply(xls, dim))
#' The data set include both source data and corrected data by qc-mxp
data <- xls$data
meta <- xls$meta
peak <- xls$peak

#' group information for data screening
names(meta)

#' need to select right column names for batch and sample type
if (T) {   # for qcmxp convert format
  (cls.qc <- factor(meta$SampleType))
  (cls.bl <- factor(meta$Batch))
} else {  # for qcrlsc format
  (cls.qc <- factor(meta$sample_type))
  (cls.bl <- factor(meta$batch))
}

## 2.) ----- Missing value filtering -----

#' let zero as NA before missing value process
data[data == 0] <- NA

#' check missing value rates
tail(sort(mv_perc(data)), 20)

#' filter based on missing values
if (T) {      # filter using all missing value rate
  ret <- mv_filter(data, thres = 0.2)
  dat <- ret$dat
  pek <- peak[ret$idx,]
} else {      # filter using qc missing value rate
  ret <- mv_filter_qc(data, cls.qc, thres = 0.2)
  dat <- ret$dat
  pek <- peak[ret$idx,]
}

## 3.) ----- QC-RLSC  -----
method <- "divide"     # "subtract"  
intra <- T
opti <- T
log10 <- T
outl <- T
shift <- F

tic()
res <- qc_rlsc_wrap(dat, cls.qc, cls.bl, method, intra, opti, log10, outl, shift)
toc()

## 4.) ----- Data screening before and after QC-RLSC  -----
#' before
dat_fill <- dat %>% mv.fill(method = "median", ze_ne = T) %>% as_tibble()
plot(pcalda(dat_fill, cls.bl))

pcaplot(dat_fill, cls.qc, pcs = c(2, 1), ep = 1)
pcaplot(dat_fill, cls.bl, pcs = c(2, 1), ep = 1)

#' after
res_fill <- res %>% mv.fill(method = "median", ze_ne = T) %>% as_tibble()
plot(pcalda(res_fill, cls.bl))

pcaplot(res_fill, cls.qc, pcs = c(2, 1), ep = 1)
pcaplot(res_fill, cls.bl, pcs = c(2, 1), ep = 1)

## 5.) ----- Save results -----
tmp <- list(data = res, meta = meta, peak = pek)
write.xlsx(tmp, file = paste0("./res/", FILE, "_res.xlsx"),
           asTable = F, overwrite = T, rowNames = F, colNames = T)
