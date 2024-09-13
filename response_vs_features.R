# Univariate analysis for drug response
suppressPackageStartupMessages({
  library(coin)
  library(dplyr)
  library(optparse)
  library(readr)
})

source("univariate_common.R")

option_list <- list(
  make_option("--response_pdata",
    action = "store",
    default = "data/response_pdata.rds",
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
  )
)

parser <- OptionParser(
  usage = "%prog [options] file",
  option_list = option_list
)

args <- parse_args(parser, positional_arguments = FALSE)

response_pdata <- args[["response_pdata"]]
kraken_meta <- args[["kraken_meta"]]
microbial_features <- args[["microbial_features"]]
feature_type <- args[["feature_type"]]
how <- args[["how"]]

features <- read_tsv(microbial_features, col_types = cols()) %>%
  filter(features == feature_type & what == "resp") %>%
  group_by(cancer, what, versus, features, genera) %>%
  tally(name = "count") %>%
  filter(count >= 2) %>%
  select(cancer, what, versus, features, genera)

set.seed(98765)
# A seed for each iteration
eff_seeds <- sample(1:2^15, 10000)

res <- readRDS(response_pdata) %>%
  select(
    sample_barcode = case_submitter_id, cancer,
    versus = drug_name, start_time = start_time,
    response
  ) %>%
  mutate(
    response = as.factor(
      ifelse(response %in% c("Complete Response", "Partial Response"),
        "yes", "no"
      )
    ),
    cancer = sub("TCGA-", "", cancer)
  )

kraken <- readRDS("data/knight_kraken_data.rds")
aliquots <- read_kraken_aliquots(kraken_meta)
results <- vector("list", nrow(features))

for (k in seq_len(nrow(features))) {
  # For each (cancer, drug, genus) triplets
  set.seed(eff_seeds[k])

  microbial_sample <- find_microbial_sample(
    kraken, aliquots, features[k, "genera", drop = TRUE]
  )
  data <- res %>%
    semi_join(features[k, , drop = FALSE], by = c("cancer", "versus")) %>%
    join_response_with_microbial(microbial_sample)

  if (length(unique(data$response)) != 2) {
    next # Don't test if there is only one response
  }

  p_value <- NA
  direction <- NA
  try({
    wc <- wilcox_test(
      formula = abundance ~ response, data = data,
      distribution = "approximate"
    )
    if (TRUE || pvalue(wc) < 0.05) {
      wc <- wilcox_test(
        formula = abundance ~ response, data = data,
        distribution = "exact"
      )
    }
    if (pvalue(wc) < 0.05) {
      ranks <- rank(data$abundance)
      mean_rank_affected <- mean(ranks[data$response == "yes"], na.rm = TRUE)
      mean_rank_unaffected <- mean(ranks[data$response == "no"], na.rm = TRUE)
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
