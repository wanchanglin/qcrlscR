#' wl-24-07-2021, Sat: set up for knitr
library(knitr)
library(here)
# library(tidyverse)
here::i_am("_common.R")
#' rm(list = ls(all = T))

#+ setup, include=FALSE
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
