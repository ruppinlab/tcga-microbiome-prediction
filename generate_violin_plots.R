options(warn=1)
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggsignif"))
suppressPackageStartupMessages(library("ggstatsplot"))
suppressPackageStartupMessages(library("pairwiseComparisons"))
suppressPackageStartupMessages(library("stringr"))

set.seed(777)

argp <- arg_parser("Create paired boxplots")
argp <- add_argument(
    argp, "--file-format", default="png", help="Save file format"
)
args <- parse_args(argp)

axis_fontsize <- 12
fig_dim <- 4
fig_dpi <- 300
line_color <- "grey"
font_family <- "Nimbus Sans"

results_dir <- "results"
out_dir <- "figures"
signif_hits_file <- "analysis/goodness_hits.txt"

cnet_model_scores <- readRDS(
    paste(results_dir, "surv", "cnet_model_scores.rds", sep="/")
)
cox_clinical_model_scores <- readRDS(
    paste(results_dir, "surv", "cox_clinical_model_scores.rds", sep="/")
)
rfe_model_scores <- readRDS(
    paste(results_dir, "resp", "rfe_model_scores.rds", sep="/")
)
svm_clinical_model_scores <- readRDS(
    paste(results_dir, "resp", "svm_clinical_model_scores.rds", sep="/")
)

signif_hits <- read.delim(signif_hits_file, stringsAsFactors=FALSE)
signif_hits <- signif_hits %>%
    mutate_if(is.character, str_to_lower) %>%
    arrange(desc(analysis), cancer, versus, features)

dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

for (row_idx in seq_len(nrow(signif_hits))) {
    cancer <- signif_hits$cancer[row_idx]
    analysis <- signif_hits$analysis[row_idx]
    target <- signif_hits$versus[row_idx]
    data_type <- signif_hits$features[row_idx]
    dataset_name <- paste("tcga", cancer, analysis, target, data_type, sep="_")
    if (data_type == "htseq")
        dataset_name <- paste(dataset_name, "counts", sep="_")
    cat(dataset_name, "\n")
    # wihin groups paired violin plot
    dataset_name_parts <- str_split(dataset_name, "_")[[1]]
    cancer <- dataset_name_parts[2]
    analysis <- dataset_name_parts[3]
    target <- dataset_name_parts[4]
    data_type <- dataset_name_parts[5]
    data_type_label <- ifelse(
        data_type == "kraken", "Microbiome",
        ifelse(data_type == "htseq", "Expression", "Combo")
    )
    model_label <- paste(data_type_label, "+", "Clinical")
    fig_num <- ifelse(
        data_type == "kraken", "1", ifelse(target == "os", "Ex1", "Ex2")
    )
    if (analysis == "surv") {
        model_scores <- cnet_model_scores[[dataset_name]]
        clinical_model_scores <- cox_clinical_model_scores[[dataset_name]]
    } else {
        model_scores <- rfe_model_scores[[dataset_name]]
        clinical_model_scores <- svm_clinical_model_scores[[dataset_name]]
    }
    data <- data.frame(
        Model=c(
            rep("Clinical", length(clinical_model_scores)),
            rep(model_label, length(model_scores))
        ),
        Score=c(clinical_model_scores, model_scores)
    )
    data$Model <- relevel(data$Model, "Clinical")
    title <- paste(str_to_upper(cancer), ifelse(
        analysis == "surv", str_to_upper(target), str_to_title(target)
    ))
    colors <- c(
        "#6f828a",
        ifelse(
            data_type == "kraken", "#448ee4",
            ifelse(data_type == "htseq", "#c04e01", "#601ef9")
        )
    )
    y_label <- ifelse(analysis == "surv", "C-index", "AUROC")
    p_adj <- signif_hits$p_adj[row_idx]
    p <- ggwithinstats(
        data=data, x="Model", y="Score", xlab="Model", ylab=y_label,
        type="nonparametric",
        centrality.plotting=TRUE, centrality.type="parametric",
        centrality.label.args=list(label.padding=0.15, size=2),
        centrality.point.args=list(color="darkred", size=2),
        point.path.args=list(alpha=0.8, color=line_color),
        p.adjust.method="BH", results.subtitle=FALSE, sample.size.label=FALSE,
        title=bquote(bold(.(title)) ~ " " ~ italic(p)[adj] == .(p_adj)),
        pairwise.comparisons=TRUE, pairwise.display="all"
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
        panel.border=element_rect(color="black", fill=NA, size=0.5),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        plot.margin=unit(c(0, 2, 0, 2), "pt"),
        plot.subtitle=element_blank(),
        plot.title=element_text(
           size=axis_fontsize, family=font_family, vjust=-1.5
        ),
        text=element_text(size=axis_fontsize, family=font_family)
    )
    if (
        analysis == "surv" && (
            (cancer == "acc" && target %in% c("os", "pfi")) ||
            (cancer == "cesc" && target == "os") ||
            (cancer == "lgg" && target == "pfi") ||
            (cancer == "skcm" && target == "os")
        )
    ) {
        break_start <- 0.2
        lim_min <- 0.18
        labels <- c("0.2", "0.4", "0.6", "0.8", "1")
    } else {
        break_start <- 0
        lim_min <- -0.01
        labels <- c("0", "0.2", "0.4", "0.6", "0.8", "1")
    }
    p <- p + scale_y_continuous(
        breaks=seq(break_start, 1, 0.2), expand=c(0, 0),
        limits=c(lim_min, 1.01), labels=labels
    )
    suppressMessages(p <- p + scale_color_manual(values=colors))
    p$layers <- p$layers[c(4, 1:3, 5:6)]
    p$layers[[2]]$aes_params$size <- 2
    file <- paste(out_dir, paste0(dataset_name, ".", args$file_format), sep="/")
    ggsave(
        file=file, plot=p, device=args$file_format, width=fig_dim,
        height=fig_dim, units="in", dpi=fig_dpi
    )
    # between groups combo comparison violin plot
    if (data_type != "combo") next
    dataset_name_parts <- str_split(dataset_name, "_")[[1]]
    kraken_dataset_name <- str_c(
        c(head(dataset_name_parts, -1), "kraken"), collapse="_"
    )
    htseq_dataset_name <- str_c(
        c(head(dataset_name_parts, -1), "htseq_counts"), collapse="_"
    )
    if (analysis == "surv") {
        kraken_model_scores <- cnet_model_scores[[kraken_dataset_name]]
        htseq_model_scores <- cnet_model_scores[[htseq_dataset_name]]
        combo_model_scores <- cnet_model_scores[[dataset_name]]
    } else {
        kraken_model_scores <- rfe_model_scores[[kraken_dataset_name]]
        htseq_model_scores <- rfe_model_scores[[htseq_dataset_name]]
        combo_model_scores <- rfe_model_scores[[dataset_name]]
    }
    data <- data.frame(
        Model=c(
            rep("Microbiome", length(kraken_model_scores)),
            rep("Expression", length(htseq_model_scores)),
            rep("Combo", length(combo_model_scores))
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
    title <- paste(str_to_upper(cancer), ifelse(
        analysis == "surv", str_to_upper(target), str_to_title(target)
    ))
    colors <- c("#448ee4", "#c04e01", "#94568c")
    y_label <- ifelse(analysis == "surv", "C-index", "AUROC")
    p <- ggbetweenstats(
        data=data, x="Model", y="Score", xlab="Model", ylab=y_label,
        type="nonparametric", point.args=list(size=2),
        centrality.plotting=TRUE, centrality.type="parametric",
        centrality.label.args=list(label.padding=0.15, size=2),
        centrality.point.args=list(color="darkred", size=2),
        p.adjust.method="none", results.subtitle=TRUE, sample.size.label=FALSE,
        title=bquote(bold(.(title))), pairwise.comparisons=FALSE
    ) +
    geom_signif(
        comparisons=pw_cmps$groups[c(1, 2)], annotations=pw_cmps$label[c(1, 2)],
        na.rm=TRUE, parse=TRUE, test=NULL, tip_length=0.02,
        y_position=c(1.005, 1.075)
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
        panel.border=element_rect(color="black", fill=NA, size=0.5),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.y=element_blank(),
        plot.margin=unit(c(0, 2, 0, 2), "pt"),
        plot.title=element_text(
           size=axis_fontsize, family=font_family, vjust=-1
        ),
        plot.subtitle=element_text(size=6, family=font_family, vjust=-1),
        text=element_text(size=axis_fontsize, family=font_family)
    ) + scale_y_continuous(
        breaks=seq(break_start, 1.2, 0.2), expand=c(0, 0),
        limits=c(lim_min, 1.2), labels=c(labels, "1.2")
    )
    suppressMessages(p <- p + scale_color_manual(values=colors))
    p$layers <- p$layers[c(1:5, 7)]
    file <- paste(out_dir, paste0(str_c(
        c(head(dataset_name_parts, -1), "comp"), collapse="_"
    ), ".", args$file_format), sep="/")
    ggsave(
        file=file, plot=p, device=args$file_format, width=5,
        height=5, units="in", dpi=fig_dpi
    )
}
