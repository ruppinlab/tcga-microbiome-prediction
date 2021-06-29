suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(reshape2)
  library(tibble)
})

methods <- c(
  lgr = "LGR", svm = "RFE", rfe = "RFE", cnet = "CNET", cox = "CNET",
  limma = "LIMMA", edger = "EDGER"
)

method_from_filename <- function(filename) {
  methods[strsplit(basename(filename), "_")[[1]][1]]
}

parse_testset <- function(testset) {
  parts <- strsplit(as.character(testset), "_")
  # The first element of parts is the string 'tcga'
  cancer <- sapply(parts, function(x) toupper(x[2]))
  analysis <- sapply(parts, function(x) x[3])
  versus <- sapply(parts, function(x) x[4])
  features <- sapply(parts, function(x) x[5])
  versus <- ifelse(
    analysis == "surv",
    toupper(versus),
    tools::toTitleCase(versus) # Drugs are title case
  )
  tibble(cancer = cancer, analysis = analysis, versus = versus, features = features)
}

filenames <- commandArgs(trailingOnly = TRUE)
results <- list()
for (i in seq_along(filenames)) {
  file <- filenames[i]
  # Models come in 0 indexed, but R doesn't appreciate that, so shift.
  model_scores <- as_tibble(readRDS(file), rownames = "index") %>%
    mutate(index = as.numeric(index) + 1) %>%
    melt("index", variable.name = "testset", value.name = "goodness")

  results[[i]] <- cbind(
    parse_testset(model_scores$testset),
    tibble(
      how = method_from_filename(file),
      index = model_scores$index,
      # Sometimes NaNs get embedded, where NA is the correct meaning.
      # Let's fix that.
      goodness = ifelse(is.nan(model_scores$goodness), NA, model_scores$goodness)
    )
  )
}
results %>%
  bind_rows() %>%
  format_tsv() %>%
  cat()
