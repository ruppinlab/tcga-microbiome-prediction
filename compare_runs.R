suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
tests <- read_tsv(args[1], col_types = cols())
covariates <- read_tsv(args[2], col_types = cols())

do_wilcox <- function(x, y) {
  test <- NULL
  try({
    test <- wilcox.test(x, y, paired = TRUE)
  })
  ifelse(is.null(test), as.double(NA), test$p.value)
}

do_onesided_wilcox <- function(x, y) {
  test <- NULL
  try({
    test <- wilcox.test(x, y, paired = TRUE, alt = "g")
  })
  ifelse(is.null(test), as.double(NA), test$p.value)
}

joined_goodness <- inner_join(
  covariates,
  tests,
  by = c("cancer", "analysis", "versus", "features", "how", "index"),
  suffix = c("_cov", "_test"),
)

joined_goodness %>%
  group_by(cancer, analysis, versus, features, how) %>%
  summarize(
    avg_test = mean(goodness_test, na.rm = TRUE),
    sd_test = sd(goodness_test, na.rm = TRUE),
    avg_cov = mean(goodness_cov, na.rm = TRUE),
    sd_cov = sd(goodness_cov, na.rm = TRUE),
    p_value = do_wilcox(goodness_test, goodness_cov),
    p_greater = do_onesided_wilcox(goodness_test, goodness_cov),
    .groups = "drop"
  ) %>%
  arrange(analysis, features, how, cancer, desc(avg_test)) %>%
  format_tsv() %>%
  cat()
