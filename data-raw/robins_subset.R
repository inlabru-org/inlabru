## code to prepare `robins_subset` dataset goes here

library(tidyverse)

# get data from tmmeeha/inlaSVCBC
download_dat <- read.csv(paste0(
  "https://raw.github.com/tmeeha/inlaSVCBC/master/code/modeling_data.csv"
))

robins_subset <- download_dat %>%
  select(
    circle, bcr, state, year, std_yr, count, log_hrs,
    lon, lat, obs
  ) %>%
  mutate(year = year + 1899) %>%
  filter(
    state %in% c(
      "TEXAS", "OKLAHOMA", "KANSAS", "MISSOURI",
      "ARKANSAS", "LOUISIANA"
    ),
    year >= 1987
  )

usethis::use_data(robins_subset, overwrite = TRUE)
