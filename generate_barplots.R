library(readr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)

outdir <- 'figures/barplots'
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

goodness <- read_tsv("analysis/goodness_hits.txt", col_types = cols())

response <- goodness %>%
  filter(analysis == "resp") %>%
  mutate(pair = paste(cancer, versus, sep = " "))

test_response <- response %>%
  select(pair, features, ROC = avg_test, sd_ROC = sd_test, p_adj) %>%
  mutate(Features = "aavg_test")
cov_response <- response %>%
  select(pair, features, ROC = avg_cov, sd_ROC = sd_cov) %>%
  mutate(p_adj = 1, Features = "avg_cov")

test_response

plots <- vector("list", 5)

values <- rbind(test_response, cov_response) %>% filter(features == "kraken")
values <- values %>% arrange(pair, Features)
values$Features[values$Features == "aavg_test"] <- "Microbiome + Covariates"
values$Features[values$Features == "avg_cov"] <- "Covariates"
values$Features <- factor(
  values$Features,
  levels = c("Expression + Covariates", "Microbiome + Covariates", "Covariates")
)

levels <- values %>%
  filter(Features == "Microbiome + Covariates") %>%
  transmute(x = seq(n()), y = pmin(ROC + sd_ROC, 1), p_adj = p_adj)

marks <- list()
stars <- rep("NS", nrow(levels))
for (i in seq(nrow(levels))) {
  mid_x <- levels$x[[i]]
  bar_y <- levels$y[[i]] + 0.05

  marks[[length(marks) + 1]] <-
    tibble(
      x = c(mid_x - 0.22, mid_x - 0.22, mid_x + 0.22, mid_x + 0.22),
      y = c(bar_y - 0.01, bar_y, bar_y, bar_y - 0.01),
      group = i
    )
  stars[levels$p_adj <= 0.01] <- "*"
  stars[levels$p_adj <= 0.001] <- "**"
  stars[levels$p_adj <= 0.0001] <- "***"
}
marks <- do.call(rbind, marks)


plot <-
  ggplot(values, aes(y = ROC, x = pair)) +
  geom_bar(position = "dodge", stat = "identity", aes(fill = Features)) +
  scale_fill_manual(values = c(cbPalette[3], cbPalette[1])) +
  xlab("Cancer Drug Pair") +
  ylab("AUROC") +
  geom_line(data = marks, aes(x = x, y = y, group = group)) +
  theme_minimal() +
  scale_y_continuous(breaks=c(0,.25,.5,.75,1)) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

for (i in seq_along(stars)) {
  plot <- plot +
    annotate("text", x = i, y = levels$y[[i]] + 0.06, label = stars[i], size = 5)
}
for (i in seq(nrow(values) / 2)) {
  sdbar <- tibble(
    x = c(i - .22, i - .22),
    y = c(
      pmax(values$ROC[[2 * i - 1]] - values$sd_ROC[[2 * i - 1]], 0),
      pmin(values$ROC[[2 * i - 1]] + values$sd_ROC[[2 * i - 1]], 1)
    )
  )
  plot <- plot + geom_line(data = sdbar, aes(x = x, y = y))
  sdbar <- tibble(
    x = c(i + .22, i + .22),
    y = c(
      pmax(values$ROC[[2 * i]] - values$sd_ROC[[2 * i]], 0),
      pmin(values$ROC[[2 * i]] + values$sd_ROC[[2 * i]], 1)
    )
  )
  plot <- plot + geom_line(data = sdbar, aes(x = x, y = y))
}
plots[[1]] <- plot
print(plots[[1]])


values <- rbind(test_response, cov_response) %>% filter(features == "htseq")
values <- values %>% arrange(pair, Features)
values$Features[values$Features == "aavg_test"] <- "Expression + Covariates"
values$Features[values$Features == "avg_cov"] <- "Covariates"
values$Features <- factor(
  values$Features,
  levels = c("Expression + Covariates", "Microbiome + Covariates", "Covariates")
)

levels <- values %>%
  filter(Features == "Expression + Covariates") %>%
  transmute(x = seq(n()), y = pmin(ROC + sd_ROC, 1), p_adj = p_adj)

marks <- list()
stars <- rep("NS", nrow(levels))
for (i in seq(nrow(levels))) {
  mid_x <- levels$x[[i]]
  bar_y <- levels$y[[i]] + 0.05

  marks[[length(marks) + 1]] <-
    tibble(
      x = c(mid_x - 0.22, mid_x - 0.22, mid_x + 0.22, mid_x + 0.22),
      y = c(bar_y - 0.01, bar_y, bar_y, bar_y - 0.01),
      group = i
    )
  stars[levels$p_adj <= 0.01] <- "*"
  stars[levels$p_adj <= 0.001] <- "**"
  stars[levels$p_adj <= 0.0001] <- "***"
}
marks <- do.call(rbind, marks)

plot <-
  ggplot(values, aes(y = ROC, x = pair)) +
  geom_bar(position = "dodge", stat = "identity", aes(fill = Features)) +
  scale_fill_manual(values = c(cbPalette[7], cbPalette[1])) +
  xlab("Cancer-Drug Pair") +
  ylab("AUROC") +
  geom_line(data = marks, aes(x = x, y = y, group = group)) +
  theme_minimal() +
  scale_y_continuous(breaks=c(0,.25,.5,.75,1)) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

for (i in seq_along(stars)) {
  plot <- plot +
    annotate("text", x = i, y = levels$y[[i]] + 0.06, label = stars[i], size = 5)
}
for (i in seq(nrow(values) / 2)) {
  sdbar <- tibble(
    x = c(i - .22, i - .22),
    y = c(
      pmax(values$ROC[[2 * i - 1]] - values$sd_ROC[[2 * i - 1]], 0),
      pmin(values$ROC[[2 * i - 1]] + values$sd_ROC[[2 * i - 1]], 1)
    )
  )
  plot <- plot + geom_line(data = sdbar, aes(x = x, y = y))
  sdbar <- tibble(
    x = c(i + .22, i + .22),
    y = c(
      pmax(values$ROC[[2 * i]] - values$sd_ROC[[2 * i]], 0),
      pmin(values$ROC[[2 * i]] + values$sd_ROC[[2 * i]], 1)
    )
  )
  plot <- plot + geom_line(data = sdbar, aes(x = x, y = y))
}
plots[[2]] <- plot
print(plots[[2]])

surv <- goodness %>%
  filter(analysis == "surv") %>%
  mutate(pair = paste(cancer, versus, sep = " ")) %>%
  select(
    pair, cancer, versus, features, avg_test, sd_test, avg_cov, sd_cov, p_adj
  )
surv

test_surv <- surv %>%
  select(
    pair, cancer, versus, features,
    ROC = avg_test, sd_ROC = sd_test, p_adj
  ) %>%
  mutate(Features = "aavg_test")
cov_surv <- surv %>%
  select(pair, cancer, versus, features, ROC = avg_cov, sd_ROC = sd_cov) %>%
  mutate(p_adj = 1, Features = "avg_cov")

values <- rbind(test_surv, cov_surv) %>% filter(features == "kraken")
values <- values %>% arrange(pair, Features)
values$Features[values$Features == "aavg_test"] <- "Microbiome + Covariates"
values$Features[values$Features == "avg_cov"] <- "Covariates"
values$Features <- factor(
  values$Features,
  levels = c("Expression + Covariates", "Microbiome + Covariates", "Covariates")
)

levels <- values %>%
  filter(Features == "Microbiome + Covariates") %>%
  transmute(x = seq(n()), y = pmin(ROC + sd_ROC, 1), p_adj = p_adj)

marks <- list()
stars <- rep("NS", nrow(levels))
for (i in seq(nrow(levels))) {
  mid_x <- levels$x[[i]]
  bar_y <- levels$y[[i]] + 0.05

  marks[[length(marks) + 1]] <-
    tibble(
      x = c(mid_x - 0.22, mid_x - 0.22, mid_x + 0.22, mid_x + 0.22),
      y = c(bar_y - 0.01, bar_y, bar_y, bar_y - 0.01),
      group = i
    )
  stars[levels$p_adj <= 0.01] <- "*"
  stars[levels$p_adj <= 0.001] <- "**"
  stars[levels$p_adj <= 0.0001] <- "***"
}
marks <- do.call(rbind, marks)

plot <-
  ggplot(values, aes(y = ROC, x = pair)) +
  geom_bar(position = "dodge", stat = "identity", aes(fill = Features)) +
  scale_fill_manual(values = c(cbPalette[3], cbPalette[1])) +
  xlab("Cancer Survival Measure") +
  ylab("C-index") +
  geom_line(data = marks, aes(x = x, y = y, group = group)) +
  theme_minimal() +
  scale_y_continuous(limits=c(0,1)) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

for (i in seq_along(stars)) {
  plot <- plot +
    annotate("text", x = i, y = levels$y[[i]] + 0.06, label = stars[i], size = 5)
}
for (i in seq(nrow(values) / 2)) {
  sdbar <- tibble(
    x = c(i - .22, i - .22),
    y = c(
      pmax(values$ROC[[2 * i - 1]] - values$sd_ROC[[2 * i - 1]], 0),
      pmin(values$ROC[[2 * i - 1]] + values$sd_ROC[[2 * i - 1]], 1)
    )
  )
  plot <- plot + geom_line(data = sdbar, aes(x = x, y = y))
  sdbar <- tibble(
    x = c(i + .22, i + .22),
    y = c(
      pmax(values$ROC[[2 * i]] - values$sd_ROC[[2 * i]], 0),
      pmin(values$ROC[[2 * i]] + values$sd_ROC[[2 * i]], 1)
    )
  )
  plot <- plot + geom_line(data = sdbar, aes(x = x, y = y))
}
plots[[3]] <- plot
print(plots[[3]])

test_surv <- surv %>%
  select(
    pair, cancer, versus, features,
    ROC = avg_test, sd_ROC = sd_test, p_adj
  ) %>%
  mutate(Features = "aavg_test")
cov_surv <- surv %>%
  select(pair, cancer, versus, features, ROC = avg_cov, sd_ROC = sd_cov) %>%
  mutate(p_adj = 1, Features = "avg_cov")

values <- rbind(test_surv, cov_surv) %>%
  filter(features == "htseq" & versus == "OS")
values <- values %>% arrange(pair, Features)
values$Features[values$Features == "aavg_test"] <- "Expression + Covariates"
values$Features[values$Features == "avg_cov"] <- "Covariates"
values$Features <- factor(
  values$Features,
  levels = c("Expression + Covariates", "Microbiome + Covariates", "Covariates")
)

levels <- values %>%
  filter(Features == "Expression + Covariates") %>%
  transmute(x = seq(n()), y = pmin(ROC + sd_ROC, 1), p_adj = p_adj)

marks <- list()
stars <- rep("NS", nrow(levels))
for (i in seq(nrow(levels))) {
  mid_x <- levels$x[[i]]
  bar_y <- levels$y[[i]] + 0.05

  marks[[length(marks) + 1]] <-
    tibble(
      x = c(mid_x - 0.22, mid_x - 0.22, mid_x + 0.22, mid_x + 0.22),
      y = c(bar_y - 0.01, bar_y, bar_y, bar_y - 0.01),
      group = i
    )
  stars[levels$p_adj <= 0.01] <- "*"
  stars[levels$p_adj <= 0.001] <- "**"
  stars[levels$p_adj <= 0.0001] <- "***"
}
marks <- do.call(rbind, marks)

plot <-
  ggplot(values, aes(y = ROC, x = cancer)) +
  geom_bar(position = "dodge", stat = "identity", width=.9, aes(fill = Features)) +  
  scale_fill_manual(values = c(cbPalette[7], cbPalette[1])) +
  xlab("Overall Survival by Cancer") +
  ylab("C-index") +
  geom_line(data = marks, aes(x = x, y = y, group = group)) +
  theme_minimal() +
    scale_y_continuous(limits=c(0,1)) +
  theme(
    plot.title = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 16),
    legend.position = 'none',
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

for (i in seq_along(stars)) {
  plot <- plot +
    annotate("text", x = i, y = levels$y[[i]] + 0.06, label = stars[i])
}
for (i in seq(nrow(values) / 2)) {
  sdbar <- tibble(
    x = c(i - .22, i - .22),
    y = c(
      pmax(values$ROC[[2 * i - 1]] - values$sd_ROC[[2 * i - 1]], 0),
      pmin(values$ROC[[2 * i - 1]] + values$sd_ROC[[2 * i - 1]], 1)
    )
  )
  plot <- plot + geom_line(data = sdbar, aes(x = x, y = y))
  sdbar <- tibble(
    x = c(i + .22, i + .22),
    y = c(
      pmax(values$ROC[[2 * i]] - values$sd_ROC[[2 * i]], 0),
      pmin(values$ROC[[2 * i]] + values$sd_ROC[[2 * i]], 1)
    )
  )
  plot <- plot + geom_line(data = sdbar, aes(x = x, y = y))
}
plots[[4]] <- plot
print(plots[[4]])

values <- rbind(test_surv, cov_surv) %>%
  filter(features == "htseq" & versus == "PFI")
values <- values %>% arrange(pair, Features)
values$Features[values$Features == "aavg_test"] <- "Expression + Covariates"
values$Features[values$Features == "avg_cov"] <- "Covariates"
values$Features <- factor(
  values$Features,
  levels = c("Expression + Covariates", "Microbiome + Covariates", "Covariates")
)

levels <- values %>%
  filter(Features == "Expression + Covariates") %>%
  transmute(x = seq(n()), y = ROC + sd_ROC, p_adj = p_adj)

marks <- list()
stars <- rep("NS", nrow(levels))
for (i in seq(nrow(levels))) {
  mid_x <- levels$x[[i]]
  bar_y <- levels$y[[i]] + 0.05

  marks[[length(marks) + 1]] <-
    tibble(
      x = c(mid_x - 0.22, mid_x - 0.22, mid_x + 0.22, mid_x + 0.22),
      y = c(bar_y - 0.01, bar_y, bar_y, bar_y - 0.01),
      group = i
    )
  stars[levels$p_adj <= 0.01] <- "*"
  stars[levels$p_adj <= 0.001] <- "**"
  stars[levels$p_adj <= 0.0001] <- "***"
}
marks <- do.call(rbind, marks)

plot <-
  ggplot(values, aes(y = ROC, x = cancer)) +
  geom_bar(position = "dodge", stat = "identity", aes(fill = Features)) +
    scale_fill_manual(values = c(cbPalette[7], cbPalette[1])) +
  xlab("Progression Free Interval by Cancer") +
  ylab("C-index") +
  geom_line(data = marks, aes(x = x, y = y, group = group)) +
  theme_minimal() +
    scale_y_continuous(limits=c(0,1)) +
  theme(
    plot.title = element_text(size = 16),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.position = 'none',
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

for (i in seq_along(stars)) {
  plot <- plot +
    annotate("text", x = i, y = levels$y[[i]] + 0.06, label = stars[i])
}
for (i in seq(nrow(values) / 2)) {
  sdbar <- tibble(
    x = c(i - .22, i - .22),
    y = c(
      pmax(values$ROC[[2 * i - 1]] - values$sd_ROC[[2 * i - 1]], 0),
      pmin(values$ROC[[2 * i - 1]] + values$sd_ROC[[2 * i - 1]], 1)
    )
  )
  plot <- plot + geom_line(data = sdbar, aes(x = x, y = y))
  sdbar <- tibble(
    x = c(i + .22, i + .22),
    y = c(
      pmax(values$ROC[[2 * i]] - values$sd_ROC[[2 * i]], 0),
      pmin(values$ROC[[2 * i]] + values$sd_ROC[[2 * i]], 1)
    )
  )
  plot <- plot + geom_line(data = sdbar, aes(x = x, y = y))
}
plots[[5]] <- plot
print(plots[[5]])

for (i in seq_along(plots)) {
  pdf(file.path(outdir, paste0("barplot", i, ".pdf")))
  print(plots[[i]])
  dev.off()
}

test_surv <- surv %>%
  select(
    pair, cancer, versus, features,
    ROC = avg_test, sd_ROC = sd_test, p_adj
  ) %>%
  mutate(Features = "aavg_test")
cov_surv <- surv %>%
  select(pair, cancer, versus, features, ROC = avg_cov, sd_ROC = sd_cov) %>%
  mutate(p_adj = 1, Features = "avg_cov")

values <- rbind(test_surv, cov_surv) %>% filter(features == "htseq")
values <- values %>% arrange(pair, Features)
values$Features[values$Features == "aavg_test"] <- "Expression + Covariates"
values$Features[values$Features == "avg_cov"] <- "Covariates"
values$Features <- factor(
  values$Features,
  levels = c("Expression + Covariates", "Microbiome + Covariates", "Covariates")
)

kraken_pairs <- surv %>% filter(features == "kraken")
values = values %>% filter(pair %in% kraken_pairs$pair)

levels <- values %>%
  filter(Features == "Expression + Covariates") %>%
  transmute(x = seq(n()), y = pmin(ROC + sd_ROC, 1), p_adj = p_adj)

marks <- list()
stars <- rep("NS", nrow(levels))
for (i in seq(nrow(levels))) {
  mid_x <- levels$x[[i]]
  bar_y <- levels$y[[i]] + 0.05

  marks[[length(marks) + 1]] <-
    tibble(
      x = c(mid_x - 0.22, mid_x - 0.22, mid_x + 0.22, mid_x + 0.22),
      y = c(bar_y - 0.01, bar_y, bar_y, bar_y - 0.01),
      group = i
    )
  stars[levels$p_adj <= 0.01] <- "*"
  stars[levels$p_adj <= 0.001] <- "**"
  stars[levels$p_adj <= 0.0001] <- "***"
}
marks <- do.call(rbind, marks)

plot <-
  ggplot(values, aes(y = ROC, x = pair)) +
  geom_bar(position = "dodge", stat = "identity", aes(fill = Features)) +
    scale_fill_manual(values = c(cbPalette[7], cbPalette[1])) +
  xlab("Cancer/Survival Measure") +
  ylab("C-index") +
  geom_line(data = marks, aes(x = x, y = y, group = group)) +
  theme_minimal() +
  scale_y_continuous(limits=c(0,1)) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)
  )

for (i in seq_along(stars)) {
  plot <- plot +
    annotate("text", x = i, y = levels$y[[i]] + 0.06, label = stars[i], size = 5)
}
for (i in seq(nrow(values) / 2)) {
  sdbar <- tibble(
    x = c(i - .22, i - .22),
    y = c(
      pmax(values$ROC[[2 * i - 1]] - values$sd_ROC[[2 * i - 1]], 0),
      pmin(values$ROC[[2 * i - 1]] + values$sd_ROC[[2 * i - 1]], 1)
    )
  )
  plot <- plot + geom_line(data = sdbar, aes(x = x, y = y))
  sdbar <- tibble(
    x = c(i + .22, i + .22),
    y = c(
      pmax(values$ROC[[2 * i]] - values$sd_ROC[[2 * i]], 0),
      pmin(values$ROC[[2 * i]] + values$sd_ROC[[2 * i]], 1)
    )
  )
  plot <- plot + geom_line(data = sdbar, aes(x = x, y = y))
}
print(plot)

pdf(file.path(outdir, 'shared_pairs.pdf'))
print(plot)
dev.off()


