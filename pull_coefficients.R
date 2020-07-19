suppressMessages({
  library(dplyr)
  library(readr)
  library(optparse)
  library(tibble)
})

options(error = traceback)

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

parse_filename <- function(filename) {
  base <- basename(filename)
  ss <- strsplit(base, "_")[[1]]
  list(
    cancer = toupper(ss[[2]]),
    what = ss[[3]],
    versus = ss[[4]],
    features = ss[[5]],
    how = ss[[6]]
  )
}

results <- list()

for (filename in filenames) {
  metadata <- parse_filename(filename)
  trials <- readRDS(filename)

  tm <- as.matrix(trials)
  tm_wilcox <- apply(
    tm, 1,
    function(x) {
      tryCatch(
        {
          wilcox.test(x)$p.value
        },
        error = function(e) as.numeric(NA)
      )
    }
  )
  tm_wilcox <- p.adjust(tm_wilcox, method = "holm")

  t.ranks <- apply(-abs(tm), 2, rank)
  t.ranks.saturated <- ifelse(t.ranks <= rank_cutoff, t.ranks, NA)

  t.cutoff <- apply(t.ranks, 1, function(x) sum(x <= rank_cutoff))
  tm_wilcox.cutoff <- !is.na(tm_wilcox) & tm_wilcox <= p_cutoff
  selected <- t.cutoff >= seen_cutoff * ncol(tm) & tm_wilcox.cutoff
  selected <- selected & !(names(selected) %in%
    c("age_at_diagnosis", "gender_male", "tumor_stage"))


  results[[length(results) + 1]] <- tibble(
    cancer = metadata$cancer,
    what = metadata$versus,
    genera = rownames(tm)[selected],
    seen = t.cutoff[selected],
    mean = apply(tm[selected, , drop = FALSE], 1, mean, na.rm = TRUE),
    median_rank = apply(
      t.ranks.saturated[selected, , drop = FALSE], 1, median,
      na.rm = TRUE
    ),
    p_value = tm_wilcox[selected]
  ) %>%
    mutate(
      what = ifelse(
        what %in% c("os", "pfi"),
        toupper(what),
        tools::toTitleCase(what)
      )
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
