suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(ggpubr)
  library(ggplot2)
  library(gridExtra)
  library(rstatix)
})

color_palette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

bar_colors <- tibble(
  htseq = color_palette[c(7, 1)],
  kraken = color_palette[c(3, 1)]
)

ylabels <- c(resp = "AUROC", OS = "C-index", PFI = "C-index")
xlabels <- c(
  OS = "Overall Survival by Cancer",
  PFI = "Progression Free Interval by Cancer",
  resp = "Cancer Drug Pair"
)

find_top_model <- function(goodness_hits, selected_hits) {
  selected <- read_tsv(goodness_hits, col_types = cols()) %>%
    semi_join(read_tsv(selected_hits, col_types = cols()),
      by = c("cancer", "analysis", "versus", "features", "how")
    )

  selected %>%
    group_by(cancer, analysis, features, versus) %>%
    arrange(p_value, desc(avg_test)) %>%
    slice(1) %>% # Take the top p-value per group
    ungroup()
}

spread_covariate_rows <- function(covariate_goodness) {
  rbind(
    covariate_goodness,
    # Limma and Edger compare to LGR.  Duplicate the covariates *but*
    # remember that LIMMA is only kraken and Edger only expression.
    covariate_goodness %>%
      filter(how == "LGR" & features == "kraken") %>%
      mutate(how = "LIMMA"),
    covariate_goodness %>%
      filter(how == "LGR" & features == "htseq") %>%
      mutate(how = "EDGER")
  )
}

read_barplot_data <- function(top_model, model_goodness_file,
                              covariate_goodness_file) {
  model_goodness <- read_tsv(
    model_goodness_file,
    col_types = cols()
  )

  covariate_goodness <- read_tsv(
    covariate_goodness_file,
    col_types = cols()
  ) %>% spread_covariate_rows()

  test_points <- model_goodness %>%
    semi_join(top_model, by = c(
      "cancer", "analysis", "features", "versus", "how"
    ))

  cov_points <- covariate_goodness %>%
    semi_join(top_model, by = c(
      "cancer", "analysis", "features", "versus", "how"
    ))

  rbind(
    test_points %>% mutate(model = "full_model"),
    cov_points %>% mutate(model = "covariates_only")
  ) %>%
    mutate(
      model = factor(model, levels = c("full_model", "covariates_only")),
      analysis = ifelse(analysis == "surv", versus, analysis)
    )
}

capped_sd <- function(x, ...) {
  mean_sd(x, ...) %>% mutate(ymax = pmin(ymax, 1))
}

generate_barplot <- function(barplot_data, colors, xlabel, ylabel) {
  if (nrow(barplot_data) == 0) {
    return(NULL)
  }

  barplot_data <- barplot_data %>%
    mutate(pair = ifelse(
      versus %in% c("OS", "PFI"),
      cancer,
      paste(cancer, versus)
    ))

  ggplot(barplot_data, aes(x = pair, y = goodness)) +
    stat_summary(
      aes(fill = model),
      fun.data = "mean_sd",
      position = position_dodge(.9),
      geom = "bar",
      show.legend = FALSE
    ) +
    geom_point(
      aes(fill = model),
      position = position_jitterdodge(
        jitter.height = 0, jitter.width = .6, dodge.width = .9
      ),
      size = .4, color = "grey20", show.legend = FALSE
    ) +
    stat_summary(
      aes(fill = model),
      fun.data = "capped_sd", geom = "errorbar",
      position = position_dodge(.9), width = .5,
      color = "gold", size = 1
    ) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                       limits = c(0, 1)) +
    xlab(xlabel) +
    ylab(ylabel) +
    theme_pubr() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 50, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      panel.grid.major.y = element_line(
        color = "grey90",
        size = 0.75,
        linetype = 1
      ),
      panel.grid.minor.y = element_line(
        color = "grey90",
        size = 0.75,
        linetype = 1
      )
    )
}

hack_widths <- function(plot, scale_factor) {
  gg <- ggplot_gtable(ggplot_build(plot))
  gg$widths[5] <- unit(scale_factor, "cm")
  gg$widths[1] <- unit(1, "null")
  gg$widths[9] <- unit(1, "null")

  gg
}

hack_stats <- function(top_model, barplot_data) {
  summary_data <- top_model %>%
    semi_join(barplot_data, by = c("cancer", "versus", "features")) %>%
    mutate(pair = paste(cancer, versus)) %>%
    arrange(pair)

  # It is difficult to create a valid rstatix object, but given a valid object,
  # easy to add manually specified p-values.  So that is what we do.
  hacked_stat <- barplot_data %>%
    mutate(pair = paste(cancer, versus)) %>%
    group_by(pair) %>%
    t_test(goodness ~ model) %>%
    arrange(pair)

  hacked_stat$p_adj <- summary_data$p_adj
  hacked_stat %>%
    add_xy_position(x = "pair", dodge = .9) %>%
    add_significance(
      p.col = "p_adj",
      cutpoints = c(0, 1e-200, 1e-04, 0.001, 0.01, 1)
    )
}

args <- commandArgs(trailingOnly = TRUE)
goodness_hits <- args[[1]]
selected_hits <- args[[2]]
model_goodness_file <- args[[3]]
covariate_goodness_file <- args[[4]]
outdir <- args[[5]]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

top_model <- find_top_model(goodness_hits, selected_hits)

raw_barplot_data <- read_barplot_data(
  top_model,
  model_goodness_file,
  covariate_goodness_file
)

barplots <- raw_barplot_data %>%
  filter(features %in% c("kraken", "htseq")) %>%
  group_by(analysis, features) %>%
  tally() %>%
  mutate(n = n / 200)


for (i in seq_len(nrow(barplots))) {
  features <- barplots[i, "features", drop = TRUE]
  analysis <- barplots[i, "analysis", drop = TRUE]
  num_bars <- barplots[i, "n", drop = TRUE]

  outfile <- ifelse(features == "htseq", "expression", "microbial")
  if (analysis == "resp") {
    outfile <- paste0(outfile, "_response")
  } else {
    outfile <- paste0(outfile, "_", tolower(analysis))
  }
  barplot_data <- raw_barplot_data %>%
    filter(analysis == !!analysis & features == !!features) %>%
    select(cancer, versus, features, index, goodness, model) %>%
    arrange(model, cancer, versus, features, index)

  p <- generate_barplot(
    barplot_data,
    colors = bar_colors[[features]],
    xlabel = xlabels[[analysis]],
    ylabel = ylabels[[analysis]]
  )
  hacked_stat <- hack_stats(top_model, barplot_data)
  pp <- p +
    stat_pvalue_manual(
      hacked_stat,
      label = "p_adj.signif",
      bracket.nudge.y = .01,
      tip.length = .005
    )

  if (analysis != "resp" && features == "htseq") {
      scale_factor <- 0.9 * num_bars
  } else {
      scale_factor <- 2 * num_bars
  }
  gtable <- hack_widths(p, scale_factor)
  pdf(file.path(outdir, paste0(outfile, ".pdf")))
  plot(gtable)
  dev.off()
  plot(gtable)

  gtable <- hack_widths(p, scale_factor)
  pdf(file.path(outdir, paste0(outfile, "_signif.pdf")))
  plot(gtable)
  dev.off()

  igoodness <- colnames(barplot_data) == "goodness"
  colnames(barplot_data)[igoodness] <- ylabels[[analysis]]
  write_tsv(barplot_data, file.path(outdir, paste0(outfile, ".tsv")))
  stats <- hacked_stat %>%
    mutate(features = ifelse(features == "htseq", "expression",
      "microbial"
    )) %>%
    select(comparision = pair, features, p_adj)

  write_tsv(stats, file.path(outdir, paste0(outfile, "_stats.tsv")))
}
