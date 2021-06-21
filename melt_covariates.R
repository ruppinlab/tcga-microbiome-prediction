suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(reshape2)
  library(tibble)
})

filenames <- commandArgs(trailingOnly = TRUE)
results <- list()

methods = c(lgr = 'LGR', svm = 'RFE', rfe = 'RFE', cnet = 'CNET', cox = 'CNET')
for (file in filenames) {
  how = methods[strsplit(basename(file), '_')[[1]][1]]
  model_scores <- as_tibble(readRDS(file), rownames = "index")

  # Models come in 0 indexed, but R doesn't appreciate that, so shift.
  model_scores$index <- as.numeric(model_scores$index) + 1

  # Sometimes NaNs get embedded, where NA is the correct meaning.
  # Let's fix that.

  for (j in seq(ncol(model_scores))) {
    model_scores[[j]][is.nan(model_scores[[j]])] <- NA
  }

  model_scores <- melt(
    model_scores, "index",
    variable.name = "testset", value.name = "goodness"
  )

  parts <- strsplit(as.character(model_scores$testset), "_")

  results[[length(results) + 1]] <- tibble(
    # The first element of parts is the string 'tcga'
    cancer = sapply(parts, function(x) toupper(x[2])),
    analysis = sapply(parts, function(x) x[3]),
    versus = sapply(parts, function(x) x[4]),
    features = sapply(parts, function(x) x[5]),
    how = how,
    index = model_scores$index,
    goodness = model_scores$goodness
  ) %>% mutate(
    versus = ifelse(
      analysis == "surv",
      toupper(versus),
      tools::toTitleCase(versus)
    )
  )
}
dataset <- do.call(rbind, results)
cat(format_tsv(dataset))
