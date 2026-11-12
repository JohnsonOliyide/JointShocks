## Converts the raw bank volatility CSVs into .rda objects under data/,
## so they can be loaded directly via data(Volatility_data_pre), etc.,
## matching the convention used by bvarsv (usmacro) and ConnectednessApproach.
## Run this script once (source it, or Ctrl+Shift+Enter in RStudio) from the
## package root; it only needs to be re-run if the source CSVs change.

read_volatility_csv <- function(path) {
  raw <- read.csv(path, check.names = FALSE)
  dates <- as.Date(raw[[1]], format = "%m/%d/%Y")
  data <- as.matrix(raw[, -1])
  colnames(data) <- colnames(raw)[-1]
  rownames(data) <- as.character(dates)
  attr(data, "date") <- dates   #Store the real Date vector as an attribute, so users don't need to re-parse rownames(data) back into dates.
  data
}

Volatility_data_pre <- read_volatility_csv("Data/Volatility_data_pre.csv")
Volatility_data_pandemic <- read_volatility_csv("Data/Volatility_data_pandemic.csv")
Volatility_data_fullsample <- read_volatility_csv("Data/Volatility_data_fullsample.csv")

usethis::use_data(Volatility_data_pre, overwrite = TRUE)
usethis::use_data(Volatility_data_pandemic, overwrite = TRUE)
usethis::use_data(Volatility_data_fullsample, overwrite = TRUE)
