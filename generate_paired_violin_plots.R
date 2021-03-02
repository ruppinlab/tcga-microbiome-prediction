options(warn=1)
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggstatsplot"))
suppressPackageStartupMessages(library("stringr"))

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
    if (analysis == "surv") {
        model_scores <- cnet_model_scores[[dataset_name]]
        clinical_model_scores <- cox_clinical_model_scores[[dataset_name]]
    } else {
        model_scores <- rfe_model_scores[[dataset_name]]
        clinical_model_scores <- svm_clinical_model_scores[[dataset_name]]
    }
    p_adj <- signif_hits$p_adj[row_idx]
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
    data <- data.frame(
        Model=c(
            rep(model_label, length(model_scores)),
            rep("Clinical", length(clinical_model_scores))
        ),
        Score=c(model_scores, clinical_model_scores)
    )
    data$Model <- relevel(data$Model, model_label)
    file <- paste(out_dir, paste0(dataset_name, ".", args$file_format), sep="/")
    title <- paste(str_to_upper(cancer), ifelse(
        analysis == "surv", str_to_upper(target), str_to_title(target)
    ))
    colors <- c(ifelse(
        data_type == "kraken", "#448ee4",
        ifelse(data_type == "htseq", "#c04e01", "#601ef9")
    ), "#6f828a")
    y_label <- ifelse(analysis == "surv", "C-index", "AUROC")
    p <- ggwithinstats(
        data=data, x="Model", y="Score", type="np", xlab="Model",
        ylab=y_label, mean.label.args=list(size=2, label.padding=0.15),
        mean.point.args=list(size=2, color="darkred"),
        point.path.args=list(alpha=0.8, color=line_color),
        results.subtitle=FALSE, sample.size.label=FALSE,
        title=bquote(bold(.(title)) ~ " " ~ italic(p)[adj] == .(p_adj))
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
        limits=c(lim_min, 1.01), labels=labels,
    )
    suppressMessages(p <- p + scale_color_manual(values=colors))
    p$layers <- p$layers[c(4, 1:3, 5:6)]
    p$layers[[2]]$aes_params$size <- 2
    ggsave(
        file=file, plot=p, device=args$file_format, width=fig_dim,
        height=fig_dim, units="in", dpi=fig_dpi
    )
}
