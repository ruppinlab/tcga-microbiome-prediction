suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

model_goodness <- args[1]
outdir <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

runs <- tibble(
  cancer = c("BLCA", "CESC", "HNSC", "SARC"),
  versus = c("Cisplatin", "Cisplatin", "Carboplatin", "Docetaxel")
)

tests <- read_tsv(model_goodness, col_types = cols()) %>%
  filter((analysis == "resp" & how == "RFE") |
    (analysis == "surv" & how == "CNET"))

micro <- tests %>% filter(features == "kraken" & analysis == "resp")
expr <- tests %>% filter(features == "htseq" & analysis == "resp")
combo <- tests %>% filter(features == "combo" & analysis == "resp")
micro_expr <- inner_join(
  micro,
  expr,
  by = c("cancer", "analysis", "versus", "how", "index")
) %>%
  rename(micro_goodness = goodness.x, expr_goodness = goodness.y) %>%
  filter(is.finite(micro_goodness) & is.finite(expr_goodness)) %>%
  select(-c(how, features.x, features.y))
joined_goodness <- inner_join(
  micro_expr,
  combo,
  by = c("cancer", "analysis", "versus", "index")
) %>%
  rename(combo_goodness = goodness) %>%
  select(-c(how, features))

plots <- list()
for (i in seq(nrow(runs))) {
  dat <- joined_goodness %>%
    filter(
      cancer == runs$cancer[i] &
        versus == runs$versus[i]
    )

  dat_plot <- data.frame(
    ROC = c(dat$micro_goodness, dat$expr_goodness, dat$combo_goodness),
    Features = c(
      rep("Microbiome", length(dat$micro_goodness)),
      rep("Expression", length(dat$expr_goodness)),
      rep("Combination", length(dat$combo_goodness))
    )
  )
  dat_plot$Features <- factor(
    dat_plot$Features,
    levels = c("Microbiome", "Expression", "Combination")
  )
  plots[[i]] <-
    ggplot(dat_plot, aes(x = ROC, fill = Features)) +
    theme_minimal() +
    scale_fill_manual(values = c(cbPalette[3], cbPalette[7], cbPalette[8])) +
    geom_density(alpha = .2) +
    xlim(0, 1) +
    xlab("AUROC") +
    ylab("Density") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20)
    ) +
    ggtitle(paste(runs$cancer[i], runs$versus[[i]]))
  print(plots[[i]])
}

runs <- tibble(
  cancer = c("THYM"),
  versus = c("OS")
)

tests <- read_tsv(model_goodness, col_types = cols()) %>%
  filter((analysis == "resp" & how == "RFE") |
    (analysis == "surv" & how == "CNET"))

micro <- tests %>% filter(features == "kraken" & analysis == "surv")
expr <- tests %>% filter(features == "htseq" & analysis == "surv")
combo <- tests %>% filter(features == "combo" & analysis == "surv")
micro_expr <- inner_join(
  micro,
  expr,
  by = c("cancer", "analysis", "versus", "how", "index")
) %>%
  rename(micro_goodness = goodness.x, expr_goodness = goodness.y) %>%
  filter(is.finite(micro_goodness) & is.finite(expr_goodness)) %>%
  select(-c(how, features.x, features.y))
joined_goodness <- inner_join(
  micro_expr,
  combo,
  by = c("cancer", "analysis", "versus", "index")
) %>%
  rename(combo_goodness = goodness) %>%
  select(-c(how, features))

surv_plots <- list()
for (i in seq(nrow(runs))) {
  dat <- joined_goodness %>%
    filter(
      cancer == runs$cancer[i] &
        versus == runs$versus[i]
    )

  dat_plot <- data.frame(
    ROC = c(dat$micro_goodness, dat$expr_goodness, dat$combo_goodness),
    Features = c(
      rep("Microbiome", length(dat$micro_goodness)),
      rep("Expression", length(dat$expr_goodness)),
      rep("Combination", length(dat$combo_goodness))
    )
  )
  dat_plot$Features <- factor(dat_plot$Features,
    levels = c("Microbiome", "Expression", "Combination")
  )
  surv_plots[[i]] <-
    ggplot(dat_plot, aes(x = ROC, fill = Features)) +
    theme_minimal() +
    scale_fill_manual(values = c(cbPalette[3], cbPalette[7], cbPalette[8])) +
    geom_density(alpha = .2) +
    xlim(0, 1) +
    xlab("AUROC") +
    ylab("Density") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20)
    ) +
    ggtitle(paste(runs$cancer[i], runs$versus[[i]]))
  print(surv_plots[[i]])
}

for (i in seq_along(plots)) {
  pdf(file.path(outdir, paste0("drug_response_combo_density", i, ".pdf")))
  print(plots[[i]])
  dev.off()
}

for (i in seq_along(surv_plots)) {
  pdf(file.path(outdir, paste0("survival_combo_density", i, ".pdf")))
  print(surv_plots[[i]])
  dev.off()
}
