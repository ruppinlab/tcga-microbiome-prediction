options(warn=1)
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggstatsplot"))
suppressPackageStartupMessages(library("stringr"))

argp <- arg_parser("Create paired boxplots")
argp <- add_argument(
    argp, "--results-dir", default="results", help="Results directory"
)
argp <- add_argument(
    argp, "--out-dir", default="figures", help="Output directory"
)
argp <- add_argument(
    argp, "--file-format", default="png", help="Save file format"
)
args <- parse_args(argp)

axis_fontsize <- 12
fig_dim <- 4
fig_dpi <- 300
line_color <- "grey"
font_family <- "Nimbus Sans"

signif_hits <- read.delim("analysis/goodness_hits.txt")
signif_hits$cancer <- str_to_lower(signif_hits$cancer)
signif_hits$analysis <- str_to_lower(signif_hits$analysis)
signif_hits$versus <- str_to_lower(signif_hits$versus)
signif_hits$features <- str_to_lower(signif_hits$features)

make_plot <- function(file, data, p_adj, title, colors, y_label) {
    p <- ggwithinstats(
        data=data, x="Model", y="Score", type="np", xlab="Model", ylab=y_label,
        mean.label.args=list(size=2, label.padding=0.15),
        mean.point.args=list(size=3, color="darkred"),
        point.path.args=list(alpha=0.8, color=line_color),
        results.subtitle=FALSE, sample.size.label=FALSE,
        title=bquote(bold(.(title)) ~ italic(p)[adj] == .(p_adj))
    ) +
    scale_y_continuous(
        breaks=seq(0, 1, 0.2), expand=c(0, 0), limits=c(-0.01, 1.01),
        labels=c("0", "0.2", "0.4", "0.6", "0.8", "1"),
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
    p$layers <- p$layers[c(4, 1:3, 5:6)]
    p$layers[[2]]$aes_params$size <- 2
    ggsave(
        file=file, plot=p, device=args$file_format, width=fig_dim,
        height=fig_dim, units="in", dpi=fig_dpi
    )
}

do_analysis <- function(model_scores, clinical_model_scores) {
    for (dataset_name in colnames(model_scores)) {
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
                rep(model_label, nrow(model_scores)),
                rep("Clinical", nrow(clinical_model_scores))
            ),
            Score=c(
                model_scores[[dataset_name]],
                clinical_model_scores[[dataset_name]]
            )
        )
        data$Model <- relevel(data$Model, model_label)
        if (any(
            signif_hits$cancer == cancer
            & signif_hits$analysis == analysis
            & signif_hits$versus == target
            & signif_hits$features == data_type
        )) {
            cat(dataset_name, "\n")
            p_adj <- signif_hits$p_adj[
                signif_hits$cancer == cancer
                & signif_hits$analysis == analysis
                & signif_hits$versus == target
                & signif_hits$features == data_type
            ]
            file <- paste(
                args$out_dir, paste0(dataset_name, ".", args$file_format),
                sep="/"
            )
            title <- paste(str_to_upper(cancer), ifelse(
                analysis == "surv", str_to_upper(target), str_to_title(target)
            ))
            colors <- c(
                ifelse(data_type == "kraken", "#448ee4", "#c04e01"), "#6f828a"
            )
            y_label <- ifelse(analysis == "surv", "C-index", "AUROC")
            make_plot(file, data, p_adj, title, colors, y_label)
        }
    }
}

dir.create(args$out_dir, showWarnings=FALSE, recursive=TRUE)

cnet_model_scores <- readRDS(
    paste(args$results_dir, "surv", "cnet_model_scores.rds", sep="/")
)
cox_clinical_model_scores <- readRDS(
    paste(args$results_dir, "surv", "cox_clinical_model_scores.rds", sep="/")
)
do_analysis(cnet_model_scores, cox_clinical_model_scores)

rfe_model_scores <- readRDS(
    paste(args$results_dir, "resp", "rfe_model_scores.rds", sep="/")
)
svm_clinical_model_scores <- readRDS(
    paste(args$results_dir, "resp", "svm_clinical_model_scores.rds", sep="/")
)
do_analysis(rfe_model_scores, svm_clinical_model_scores)
