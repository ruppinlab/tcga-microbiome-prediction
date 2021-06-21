library(readr)
library(ggplot2)
library(dplyr)
library(plotrix)
library(RColorBrewer)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

args <- commandArgs(trailingOnly = TRUE)
goodness_hits <- args[1]
model_goodness <- args[2]
covariate_goodness <- args[3]
outdir <- args[4]

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## Find candidate drug responses
compare_runs <- read_tsv(goodness_hits, col_types = cols())
compare_runs <- compare_runs %>%
  filter(
    analysis != "rest" &
      how %in% c('CNET', 'RFE') &
      features == "htseq" &
      avg_test >= .6
  )

## Load and preprocess the ROC scores
tests <- read_tsv(model_goodness, col_types = cols())
covariates <- read_tsv(covariate_goodness, col_types = cols())
joined_goodness <- inner_join(
  covariates,
  tests,
  by = c("cancer", "analysis", "versus", "features", "how", "index")
) %>%
  rename(cov_goodness = goodness.x, test_goodness = goodness.y) %>%
  filter(is.finite(cov_goodness) & is.finite(test_goodness)) %>%
  select(-c(how))

runs <- compare_runs %>%
  filter(analysis == "resp") %>%
  distinct(cancer, versus) %>%
  arrange(cancer, versus)
runs

plot_runs <- function(runs, analysis) {
  plots <- list()
  for (i in seq(nrow(runs))) {
    dat <- joined_goodness %>%
      filter(
        cancer == runs$cancer[i] &
          analysis == !!analysis &
          features == "htseq" &
          versus == runs$versus[i]
      )
    p_adj <- compare_runs %>%
      filter(
        cancer == runs$cancer[i] &
          analysis == !!analysis &
          features == "htseq" &
          versus == runs$versus[i]
      ) %>%
      pull(p_adj)

    if (p_adj <= 0.0001) {
      stars <- "***"
    } else if (p_adj <= 0.001) {
      stars <- "**"
    } else if (p_adj <= 0.01) {
      stars <- "*"
    } else {
      stars <- "NS"
    }

    dat_plot <- data.frame(
      ROC = c(dat$test_goodness, dat$cov_goodness),
      Features = c(
        rep("Expression + Covariates", length(dat$test_goodness)),
        rep("Covariates", length(dat$cov_goodness))
      )
    )
    dat_plot$Features <- factor(dat_plot$Features, levels = c("Expression + Covariates", "Covariates"))
    plots[[i]] <-
      ggplot(data = dat_plot) +
      theme_minimal() +
      scale_fill_manual(values = c(cbPalette[3], cbPalette[1])) +
      xlim(0, 1) +
      xlab("AUROC") +
      ylab("Density") +
      theme(
        legend.position = c(0.25, 0.85),
        plot.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20)
      ) +
      ggtitle(paste(runs$cancer[i], runs$versus[[i]])) +
      geom_density(alpha = .2, aes(x = ROC, fill = Features), show.legend = FALSE)

    y_range <- layer_scales(plots[[i]])$y$range$range
    ymid <- mean(y_range)
    xx <- mean(dat$test_goodness)
    sigma <- std.error(dat$test_goodness, na.rm = T)

    xxx <- mean(dat$cov_goodness)
    tau <- std.error(dat$cov_goodness, na.rm = T)

    # Test mean and error bars
    plots[[i]] <- plots[[i]] +
      geom_segment(
        aes(
          x = !!(xx - sigma),
          xend = !!(xx + sigma),
          y = !!(1.1 * ymid),
          yend = !!(1.1 * ymid)
        ),
        color = cbPalette[3], size = 1
      ) +
      geom_segment(
        aes(
          x = !!(xx - sigma),
          xend = !!(xx - sigma),
          y = !!(1.1 * ymid - 0.05),
          yend = !!(1.1 * ymid + 0.05)
        ),
        color = cbPalette[3], size = 1
      ) +
      geom_segment(
        aes(
          x = !!(xx + sigma),
          xend = !!(xx + sigma),
          y = !!(1.1 * ymid - 0.05),
          yend = !!(1.1 * ymid + 0.05)
        ),
        color = cbPalette[3], size = 1
      ) +
      geom_segment(aes(x = !!xx, xend = !!xx, y = 0, yend = !!(1.8 * ymid)),
        color = cbPalette[3], linetype = "dashed", size = 1
      )

    # Covariate mean and error bars
    plots[[i]] <- plots[[i]] +
      geom_segment(
        aes(
          x = !!(xxx - tau),
          xend = !!(xxx + tau),
          y = !!(.9 * ymid),
          yend = !!(.9 * ymid)
        ),
        color = cbPalette[1], size = 1
      ) +
      geom_segment(
        aes(
          x = !!(xxx - tau),
          xend = !!(xxx - tau),
          y = !!(.9 * ymid - 0.05),
          yend = !!(.9 * ymid + 0.05)
        ),
        color = cbPalette[1], size = 1
      ) +
      geom_segment(
        aes(
          x = !!(xxx + tau),
          xend = !!(xxx + tau),
          y = !!(.9 * ymid - 0.05),
          yend = !!(.9 * ymid + 0.05)
        ),
        color = cbPalette[1], size = 1
      ) +
      geom_segment(aes(x = !!xxx, xend = !!xxx, y = 0, yend = !!(1.8 * ymid)),
        color = cbPalette[1], linetype = "dashed", size = 1
      )

    # Significance marker
    plots[[i]] <- plots[[i]] +
      geom_segment(aes(x = !!xxx, xend = !!xx, y = !!(1.95 * ymid), yend = !!(1.95 * ymid)),
        size = .5
      ) +
      geom_segment(aes(x = !!xxx, xend = !!xxx, y = !!(1.95 * ymid), yend = !!(1.9 * ymid)),
        size = .5
      ) +
      geom_segment(aes(x = !!xx, xend = !!xx, y = !!(1.95 * ymid), yend = !!(1.9 * ymid)),
        size = .5
      ) +
      annotate("text", x = (xx + xxx) / 2, y = 2 * ymid, label = stars, size = 6)
  }
  plots
}

plots <- plot_runs(runs, "resp")
cat("length plots", length(plots), "\n")

runs <- compare_runs %>%
  filter(analysis == "surv") %>%
  distinct(cancer, versus) %>%
  arrange(cancer, versus)
runs

surv_plots <- plot_runs(runs, "surv")

cat("length surv_plots", length(surv_plots), "\n")

for (i in seq_along(plots)) {
  pdf(file.path(outdir, paste0("drug_response_expression_density", i, ".pdf")))
  print(plots[[i]])
  dev.off()
}

for (i in seq_along(surv_plots)) {
  pdf(file.path(outdir, paste0("survival_expression_density", i, ".pdf")))
  print(surv_plots[[i]])
  dev.off()
}
