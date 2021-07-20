# Univariate analysis for survival
suppressPackageStartupMessages({
  library(tidyverse)
  library(coin)
  library(optparse)
})

source("univariate_common.R")

long_or_short_survival <- function(s, s_time) {
  mid <- median(Surv(s_time, s), na.rm = TRUE)
  mid <- as.numeric(mid[["quantile"]])
  survival <- as.factor(ifelse(s_time >= mid, "b_long", "a_short"))
  # Remove persons lost to follow up before median time
  survival[s_time < mid & s != 1] <- NA
  survival
}

option_list <- list(
  make_option("--survival_pdata",
    action = "store",
    default = "data/survival_pdata.rds",
    help = "Drug response data"
  ),
  make_option("--kraken_meta",
    action = "store",
    default = "data/knight_kraken_meta.rds",
    help = "Metadata for the microbial features"
  ),
  make_option("--microbial_features",
    action = "store",
    type = "character",
    help = "List of nominally significant microbial features"
  ),
  make_option("--feature_type",
    action = "store",
    default = "kraken",
    help = "Perform analysis for which type of feature"
  ),
  make_option("--how",
    action = "store",
    default = "RFE",
    help = "ML method used to generate the features"
  )
)

parser <- OptionParser(
  usage = "%prog [options] file",
  option_list = option_list
)

args <- parse_args(parser, positional_arguments = FALSE)

survival_pdata <- args[["survival_pdata"]]
kraken_meta <- args[["kraken_meta"]]
microbial_features <- args[["microbial_features"]]
feature_type <- args[["feature_type"]]
how <- args[["how"]]

features <- read_tsv(microbial_features, col_types = cols()) %>%
  filter(features == feature_type & what == "surv" & how == !!how)

set.seed(98765)
# A seed for each iteration
eff_seeds <- sample(1:2^15, 10000)

surv <- readRDS("data/survival_pdata.rds") %>%
  rename(sample_barcode = case_submitter_id) %>%
  mutate(cancer = sub("TCGA-", "", cancer)) %>%
  group_by(cancer) %>%
  mutate(
    PFI_survival = long_or_short_survival(PFI, PFI_time),
    OS_survival = long_or_short_survival(OS, OS_time)
  ) %>%
  ungroup()

kraken <- readRDS("data/knight_kraken_data.rds")
aliquots <- read_kraken_aliquots(kraken_meta)
results <- vector("list", nrow(features))

for (k in seq_len(nrow(features))) {
  # For each (cancer, drug, genus) triplets
  set.seed(eff_seeds[k])

  microbial_sample <- find_microbial_sample(
    kraken, aliquots, features[k, "genera", drop = TRUE]
  )
  what <- features[k, "versus", drop = TRUE]
  cancer <- features[k, "cancer", drop = TRUE]
  what_time <- paste0(what, "_time")
  what_survival <- paste0(what, "_survival")

  data <- surv %>%
    select(sample_barcode, cancer,
      s = one_of(what),
      s_time = one_of(what_time),
      survival = one_of(what_survival)
    ) %>%
    filter(cancer == !!cancer & !is.na(survival)) %>%
    join_survival_with_microbial(microbial_sample)

  if (length(unique(data$survival)) != 2) {
    next
  }

  p_value <- NA
  direction <- NA
  try({
    wc <- wilcox_test(
      formula = abundance ~ survival, data = data,
      distribution = "approximate"
    )
    if (TRUE || pvalue(wc) < 0.05) {
      wc <- wilcox_test(
        formula = abundance ~ survival, data = data,
        distribution = "exact"
      )
    }
    if (pvalue(wc) < 0.05) {
      ranks <- rank(data$abundance)
      mean_rank_affected <- mean(ranks[data$survival == "b_long"])
      mean_rank_unaffected <- mean(ranks[data$survival == "a_short"])
      # Split
      direction <- ifelse(mean_rank_affected > mean_rank_unaffected, 1, -1)
    }
    p_value <- pvalue(wc)
  })

  results[[k]] <-
    tibble(
      features[k, , drop = FALSE],
      univariate_direction = !!direction,
      univariate_p_value = !!p_value,
      seed = eff_seeds[k]
    )
}

bind_rows(results) %>%
  group_by(cancer, versus) %>%
  mutate(univariate_fdr = p.adjust(univariate_p_value, "BH")) %>%
  format_tsv() %>%
  cat()
