suppressPackageStartupMessages({
    library(argparser)
    library(dplyr)
    library(extrafont)
    library(gdtools)
    library(ggplot2)
    library(ggsignif)
    library(ggstatsplot)
    library(ggtext)
    library(pairwiseComparisons)
    library(stringr)
    library(tibble)
    library(readr)
})

set.seed(777)

argp <- arg_parser("Create violin boxplots")
argp <- add_argument(
    argp, "--results-dir", default="results", help="Results directory"
)
argp <- add_argument(
    argp, "--out-dir", default="figures/violin", help="Figures output directory"
)
argp <- add_argument(
    argp, "--filter", default="compared_runs", help="Figure output filter"
)
argp <- add_argument(
    argp, "--file-format", default="pdf", help="Save file format"
)
args <- parse_args(argp)

stopifnot(
    args$filter %in% c("goodness_hits", "potential_hits", "compared_runs")
)

title_fontsize <- 20
axis_fontsize <- 18
point_size <- 2
fig_dim <- 4
fig_dpi <- 300
line_color <- "grey"

font_family <- ifelse(
    font_family_exists("Helvetica"), "Helvetica", ifelse(
        font_family_exists("Nimbus Sans"), "Nimbus Sans", ifelse(
            font_family_exists("Arial"), "Arial", ifelse(
                font_family_exists("DejaVu Sans"), "DejaVu Sans", "sans"
            )
        )
    )
)
if (args$file_format == "pdf") font_family <- "sans"

model_results_dir <- paste(args$results_dir, "models", sep="/")
analysis_results_dir <- paste(args$results_dir, "analysis", sep="/")

signif_hits_file <- paste(
    analysis_results_dir, paste0(args$filter, ".txt"), sep="/"
)

signif_hits <- read.delim(signif_hits_file, stringsAsFactors=FALSE)
signif_hits <- signif_hits %>%
    mutate_if(is.character, str_to_lower) %>%
    arrange(desc(analysis), cancer, versus, features)

surv_model_scores <- readRDS(paste(
    model_results_dir, "surv", "cnet_model_scores.rds", sep="/"
))
surv_clinical_model_scores <- readRDS(paste(
    model_results_dir, "surv", "cox_clinical_model_scores.rds", sep="/"
))
resp_model_scores <- cbind(
    readRDS(paste(
        model_results_dir, "resp", "edger_model_scores.rds", sep="/"
    )),
    readRDS(paste(
        model_results_dir, "resp", "lgr_model_scores.rds", sep="/"
    )),
    readRDS(paste(
        model_results_dir, "resp", "limma_model_scores.rds", sep="/"
    )),
    readRDS(paste(
        model_results_dir, "resp", "rfe_model_scores.rds", sep="/"
    ))
)
resp_clinical_model_scores <- cbind(
    readRDS(paste(
        model_results_dir, "resp", "lgr_clinical_model_scores.rds", sep="/"
    )),
    readRDS(paste(
        model_results_dir, "resp", "svm_clinical_model_scores.rds", sep="/"
    ))
)

dir.create(args$out_dir, showWarnings=FALSE, recursive=TRUE)

for (row_idx in seq_len(nrow(signif_hits))) {
    cancer <- signif_hits$cancer[row_idx]
    analysis <- signif_hits$analysis[row_idx]
    target <- signif_hits$versus[row_idx]
    data_type <- signif_hits$features[row_idx]
    model_code <- signif_hits$how[row_idx]
    if (data_type == "htseq") {
        model_name <- paste(
            "tcga", cancer, analysis, target, data_type, "counts", model_code,
            sep="_"
        )
    } else {
        model_name <- paste(
            "tcga", cancer, analysis, target, data_type, model_code, sep="_"
        )
    }
    cat(model_name, "\n")
    clinical_model_code <- ifelse(
        model_code %in% c("cnet"), "cox_clinical", ifelse(
            model_code %in% c("rfe"), "svm_clinical", "lgr_clinical"
        )
    )
    clinical_model_name <- str_replace(
        model_name, paste0(model_code, "$"), clinical_model_code
    )
    #
    # wihin groups paired violin plot
    #
    dtype_label <- ifelse(
        data_type == "kraken", "Microbiome",
        ifelse(data_type == "htseq", "Expression", "Combined")
    )
    abbr_dtype_label <- ifelse(
        data_type == "kraken", "Microbe",
        ifelse(data_type == "htseq", "Express", "Combo")
    )
    fig_num <- ifelse(
        data_type == "kraken", "1", ifelse(target == "os", "Ex1", "Ex2")
    )
    if (analysis == "surv") {
        model_scores <- surv_model_scores[[model_name]]
        clinical_model_scores <-
            surv_clinical_model_scores[[clinical_model_name]]
    } else {
        model_scores <- resp_model_scores[[model_name]]
        clinical_model_scores <-
            resp_clinical_model_scores[[clinical_model_name]]
    }
    data <- data.frame(
        Model=c(
            rep("Clinical", length(clinical_model_scores)),
            rep(paste0(dtype_label, "+", "Clinical"), length(model_scores))
        ),
        Score=c(clinical_model_scores, model_scores)
    )
    data$Model <- relevel(data$Model, "Clinical")
    colors <- c(
        "#6f828a",
        ifelse(
            data_type == "kraken", "#448ee4",
            ifelse(data_type == "htseq", "#c04e01", "#601ef9")
        )
    )
    y_label <- ifelse(analysis == "surv", "C-index", "AUROC")
    p_adj <- signif_hits$p_adj[row_idx]
    p_greater <- signif_hits$p_greater[row_idx]
    title <- str_to_upper(cancer)
    if (analysis == "surv") {
        title <- paste(title, str_to_upper(target))
    } else {
        title <- paste(title, target)
        title <- paste0(title, " (", str_to_upper(model_code), ")")
    }
    p_adj_md <- ifelse(
        !is.na(p_adj), paste0(
            "p<sub>adj</sub> = ", ifelse(
                p_greater <= 0.05, sprintf("%.2e", p_adj), paste0(
                    "<span style='color:red'>", sprintf("%.2e", p_adj),
                    "</span>"
                )
            )
        ), ""
    )
    if (analysis == "surv") {
        title <- paste(title, p_adj_md)
        subtitle <- NULL
        plot_title <- element_markdown(
            size=title_fontsize, family=font_family, face="plain",
            margin=margin(5, 0, 1, 0), padding=margin(0, 0, 0, 0)
        )
        plot_subtitle <- element_blank()
    } else {
        subtitle <- p_adj_md
        plot_title <- element_text(
           color="black", face="plain", size=title_fontsize,
           margin=margin(5, 0, 1, 0)
        )
        plot_subtitle <- element_markdown(
           color="black", face="plain", size=title_fontsize,
           margin=margin(2, 0, 1, 0), padding=margin(0, 0, 0, 0)
        )
    }
    if (analysis == "surv" && !(cancer %in% c(
        "chol", "dlbc", "pcpg", "prad", "tgct", "thym", "uvm"
    ))) {
        break_start <- 0.2
        ylim_min <- 0.19
        ylim_max <- 1.2
        labels <- c("0.2", "0.4", "0.6", "0.8", "1")
    } else {
        break_start <- 0
        ylim_min <- -0.01
        ylim_max <- 1.3
        labels <- c("0", "0.2", "0.4", "0.6", "0.8", "1")
    }
    p <- ggwithinstats(
        data=data, x="Model", y="Score", xlab="Model", ylab=y_label,
        type="np", centrality.plotting=TRUE, centrality.type="p",
        centrality.label.args=list(label.padding=0.15, size=point_size),
        centrality.path.args=list(size=1.25, color="red", alpha=0.5),
        centrality.point.args=list(color="darkred", size=point_size),
        point.path.args=list(size=0.6, alpha=0.8, color=line_color),
        p.adjust.method="BH", results.subtitle=FALSE, title=title,
        subtitle=subtitle, pairwise.comparisons=TRUE, pairwise.display="all"
    ) +
    theme(
        aspect.ratio=1,
        axis.title.x=element_blank(),
        axis.text.x=element_text(
            color="black", face="plain", size=axis_fontsize
        ),
        axis.title.y=element_text(
            color="black", face="plain", size=axis_fontsize
        ),
        axis.text.y=element_text(
            color="black", face="plain", size=axis_fontsize
        ),
        axis.line = element_line(color="black", size=0.75),
        # panel.border=element_rect(color="black", fill=NA, size=0.5),
        panel.border=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        plot.margin=unit(c(0, 10, 0, 2), "pt"),
        plot.title=plot_title,
        plot.subtitle=plot_subtitle,
        text=element_text(size=axis_fontsize, family=font_family)
    ) + scale_y_continuous(
        breaks=seq(break_start, 1, 0.2), expand=c(0, 0),
        limits=c(ylim_min, 1.01), labels=labels
    )
    suppressMessages({
        p <- p + scale_color_manual(values=colors) + scale_x_discrete(
            labels=c("Clinical", paste0(abbr_dtype_label, "+", "Clinical"))
        )
    })
    p$layers <- p$layers[c(4, 1:3, 5:6)]
    # alter point size
    p$layers[[2]]$aes_params$size <- point_size
    # alter boxplot line width
    p$layers[[3]]$aes_params$size <- 0.6
    # alter violin line width
    p$layers[[4]]$aes_params$size <- 0.6
    file <- paste(args$out_dir, paste0(
        str_c(c(model_name, "violin"), collapse="_"), ".", args$file_format
    ), sep="/")
    ggsave(
        file=file, plot=p, device=args$file_format, width=fig_dim,
        height=fig_dim, units="in", dpi=fig_dpi
    )
    if (args$file_format == "pdf") embed_fonts(file)
    # include the p_adj score because it is read, not computed
    data_tsv <- tibble(
        cancer = cancer, target = target,
        data_type = tolower(dtype_label),
        model_code = model_code,
        index = seq_along(clinical_model_scores),
        clinical_model_scores = clinical_model_scores,
        model_scores = model_scores,
        p_adj = ifelse(p_greater <= 0.05, p_adj, 1))
    tfile <- paste(
        args$out_dir,
        paste0(str_c(c(model_name, "violin"), collapse="_"), ".tsv"),
        sep="/"
    )
    write_tsv(data_tsv, tfile)
    #
    # between groups combo comparison violin plot
    #
    if (data_type != "combo") next
    model_name_parts <- str_split(model_name, "_")[[1]]
    kraken_model_name <- str_c(c(
        head(model_name_parts, -2), "kraken", tail(model_name_parts, 1)
    ), collapse="_")
    htseq_model_name <- str_c(c(
        head(model_name_parts, -2), "htseq_counts",
        ifelse(model_code == "limma", "edger", tail(model_name_parts, 1))
    ), collapse="_")
    if (analysis == "surv") {
        kraken_model_scores <- surv_model_scores[[kraken_model_name]]
        htseq_model_scores <- surv_model_scores[[htseq_model_name]]
        combo_model_scores <- surv_model_scores[[model_name]]
    } else {
        kraken_model_scores <- resp_model_scores[[kraken_model_name]]
        htseq_model_scores <- resp_model_scores[[htseq_model_name]]
        combo_model_scores <- resp_model_scores[[model_name]]
    }
    data <- data.frame(
        Model=c(
            rep("Microbiome", length(kraken_model_scores)),
            rep("Expression", length(htseq_model_scores)),
            rep("Combined", length(combo_model_scores))
        ),
        Score=c(kraken_model_scores, htseq_model_scores, combo_model_scores)
    )
    data$Model <- relevel(data$Model, "Expression")
    data$Model <- relevel(data$Model, "Microbiome")
    pw_cmps <- pairwise_comparisons(
        data=data, x="Model", y="Score", type="nonparametric", paired=FALSE,
        p.adjust.method="BH"
    ) %>%
    mutate(
        .data = ., groups = purrr::pmap(.l = list(group1, group2), .f = c)
    ) %>%
    arrange(.data = ., group1)
    pw_cmps$label <- str_replace(
        pw_cmps$label, fixed("italic(p)[FDR-corrected]"), "p[adj]"
    )
    title <- str_to_upper(cancer)
    if (analysis == "surv") {
        title <- paste(title, str_to_upper(target))
    } else {
        title <- paste(title, target)
        title <- paste0(title, " (", str_to_upper(model_code), ")")
    }
    colors <- c("#448ee4", "#c04e01", "#94568c")
    y_label <- ifelse(analysis == "surv", "C-index", "AUROC")
    p <- ggbetweenstats(
        data=data, x="Model", y="Score", xlab="Model", ylab=y_label,
        type="np", point.args=list(size=point_size),
        centrality.plotting=TRUE, centrality.type="p",
        centrality.label.args=list(label.padding=0.15, size=point_size),
        centrality.point.args=list(color="darkred", size=point_size),
        p.adjust.method="none", results.subtitle=FALSE,
        title=title, pairwise.comparisons=FALSE
    ) +
    geom_signif(
        comparisons=pw_cmps$groups[c(1, 2)], annotations=pw_cmps$label[c(1, 2)],
        na.rm=TRUE, parse=TRUE, test=NULL, step_increase=0.07, textsize=5,
        tip_length=0.02, vjust=0.07, y_position=c(1.005, 1.080)
    ) +
    theme(
        aspect.ratio=1,
        axis.title.x=element_blank(),
        axis.text.x=element_text(
            color="black", face="plain", size=axis_fontsize
        ),
        axis.title.y=element_text(
            color="black", face="plain", size=axis_fontsize
        ),
        axis.text.y=element_text(
            color="black", face="plain", size=axis_fontsize
        ),
        axis.line = element_line(color="black", size=0.75),
        # panel.border=element_rect(color="black", fill=NA, size=0.5),
        panel.border=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        plot.margin=unit(c(0, 10, 0, 2), "pt"),
        plot.title=element_text(
           color="black", face="plain", size=title_fontsize,
           margin=margin(5, 0, 1, 0)
       ),
        plot.subtitle=element_blank(),
        text=element_text(size=axis_fontsize, family=font_family)
    ) + scale_y_continuous(
        breaks=seq(break_start, 1, 0.2), expand=c(0, 0),
        limits=c(ylim_min, ylim_max), labels=labels
    )
    suppressMessages({
        p <- p + scale_color_manual(values=colors) + scale_x_discrete(
            labels=c("Microbe", "Express", "Combo")
        )
    })
    p$layers <- p$layers[c(1:5, 7)]
    # alter boxplot line width
    p$layers[[3]]$aes_params$size <- 0.6
    # alter violin line width
    p$layers[[4]]$aes_params$size <- 0.6
    file <- paste(args$out_dir, paste0(
        str_c(c(model_name, "violin", "comp"), collapse="_"), ".",
        args$file_format
    ), sep="/")
    ggsave(
        file=file, plot=p, device=args$file_format, width=fig_dim,
        height=fig_dim, units="in", dpi=fig_dpi
    )
    if (args$file_format == "pdf") embed_fonts(file)
    data_tsv <- tibble(
       cancer = cancer, target = target,
       data_type = tolower(dtype_label),
       index = seq_along(kraken_model_scores),
       kraken_model_scores = kraken_model_scores,
       htseq_model_scores = htseq_model_scores,
       combo_model_scores = combo_model_scores)
    tfile <- paste(
        args$out_dir,
        paste0(str_c(c(model_name, "violin", "comp"), collapse="_"), ".tsv"),
        sep="/"
    )
    write_tsv(data_tsv, tfile)
}
