suppressMessages({
  library(dplyr)
  library(optparse)
  library(readr)
  library(tibble)
})

ROW_MARGIN <- 1
COL_MARGIN <- 2
CLINICAL_FEATURES <- c("age_at_diagnosis", "gender_male", "tumor_stage")

parse_filename <- function(filename) {
  ss <- as.list(strsplit(basename(filename), "_")[[1]])
  # ignore the first element
  names(ss) <- c("ignore", "cancer", "what", "versus", "features", "how")
  ss$cancer <- toupper(ss$cancer)
  ss$versus <- versus <- ifelse(
    ss$what == "surv",
    toupper(ss$versus),
    tools::toTitleCase(ss$versus) # Drugs are title case
  )
  ss$how <- toupper(ss$how)
  ss
}

wilcox_feature_not_zero <- function(trials) {
  apply(
    trials,
    ROW_MARGIN,
    function(x) {
      tryCatch(
        {
          wilcox.test(x)$p.value
        },
        error = function(e) as.numeric(NA)
      )
    }
  )
}

option_list <- list(
  make_option("--rank-cutoff",
    action = "store",
    default = 50,
    type = "integer",
    help = "Only consider args of lower rank"
  ),
  make_option("--seen-cutoff",
    action = "store",
    default = 0.2,
    type = "double",
    help = "Fraction that must contain the feature"
  ),
  make_option("--p-cutoff",
    action = "store",
    default = 0.01,
    type = "double",
    help = "P-value that feature is away from zero"
  )
)

parser <- OptionParser(
  usage = "%prog [options] file",
  option_list = option_list
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
rank_cutoff <- opt[["rank-cutoff"]]
seen_cutoff <- opt[["seen-cutoff"]]
p_cutoff <- opt[["p-cutoff"]]

filenames <- arguments$args

results <- vector("list", length(filenames))
for (i in seq_along(filenames)) {
  filename <- filenames[i]
  metadata <- parse_filename(filename)
  trials <- readRDS(filename) %>% as.matrix()

  features_wilcox <- p.adjust(wilcox_feature_not_zero(trials), method = "holm")

  feature_ranks <- apply(-abs(trials), COL_MARGIN, rank)
  feature_ranks_used <- ifelse(feature_ranks <= rank_cutoff, feature_ranks, NA)
  feature_seen <- apply(
    feature_ranks_used, ROW_MARGIN, function(x) sum(!is.na(x))
  )

  selected <- (feature_seen >= seen_cutoff * ncol(trials)) &
    (!is.na(features_wilcox) & features_wilcox <= p_cutoff) &
    !(names(feature_seen) %in% CLINICAL_FEATURES)

  results[[i]] <- tibble(
    cancer = metadata$cancer,
    what = metadata$versus,
    features = metadata$features,
    how = metadata$how,
    genera = rownames(trials)[selected],
    seen = feature_seen[selected],
    mean = apply(
      trials[selected, , drop = FALSE], ROW_MARGIN, mean,
      na.rm = TRUE
    ),
    median_rank = apply(
      feature_ranks_used[selected, , drop = FALSE], ROW_MARGIN, median,
      na.rm = TRUE
    ),
    p_value = features_wilcox[selected]
  )
}

do.call(rbind, results) %>%
  arrange(cancer, what, desc(mean)) %>%
  mutate(
    p_value = sprintf("%.3g", p_value),
    mean = sprintf("%.3g", mean),
    median_rank = sprintf("%.1f", median_rank)
  ) %>%
  format_tsv() %>%
  cat()
