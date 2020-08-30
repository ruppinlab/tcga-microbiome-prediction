options(warn=1)
suppressPackageStartupMessages(library("argparser"))
suppressPackageStartupMessages(library("ggpubr"))
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
line_size <- 0.3
line_color <- "grey"
font_family <- "Nimbus Sans"

make_plot <- function(file, data, p_adj, title, colors, y_label) {
    stat_test <- tibble::tribble(
        ~group1, ~group2, ~p.adj, ~p.adj.signif,
        colnames(data)[1], colnames(data)[2], p_adj,
        ifelse(
            p_adj <= 0.0001, "****", ifelse(
                p_adj <= 0.001, "***", ifelse(
                    p_adj <= 0.01, "**", ifelse(
                        p_adj <= 0.05, "*", "ns"
                    )
                )
            )
        )
    )
    p <- ggpaired(
        data, colnames(data)[1], colnames(data)[2], color="condition",
        legend="None", line.color=line_color, line.size=line_size,
        palette=colors, subtitle=title, xlab=FALSE, ylab=y_label
    ) +
    scale_y_continuous(
        breaks=seq(0, 1, 0.2), expand=c(0, 0),
        labels=c("0", "0.2", "0.4", "0.6", "0.8", "1"), limits=c(-0.05, 1.08)
    ) +
    stat_pvalue_manual(stat_test, y.position=1.05, label="p.adj.signif") +
    theme(
        aspect.ratio=1,
        text=element_text(size=axis_fontsize, family=font_family)
    )
    p$layers <- p$layers[c(2, 1, 3, 4)]
    ggsave(
        file=file, plot=p, device=args$file_format, width=fig_dim,
        height=fig_dim, units="in", dpi=fig_dpi
    )
}

dir.create(args$out_dir, showWarnings=FALSE, recursive=TRUE)

potential_hits <- read.delim("analysis/potential_hits.txt")
potential_hits$cancer <- str_to_lower(potential_hits$cancer)
potential_hits$analysis <- str_to_lower(potential_hits$analysis)
potential_hits$versus <- str_to_lower(potential_hits$versus)
potential_hits$features <- str_to_lower(potential_hits$features)

cnet_model_scores <- readRDS(
    paste(args$results_dir, "surv", "cnet_model_scores.rds", sep="/")
)
cox_clinical_model_scores <- readRDS(
    paste(args$results_dir, "surv", "cox_clinical_model_scores.rds", sep="/")
)
for (dataset_name in colnames(cnet_model_scores)) {
    cat(dataset_name, "\n")
    dataset_name_parts <- str_split(dataset_name, "_")[[1]]
    cancer <- dataset_name_parts[2]
    analysis <- dataset_name_parts[3]
    target <- dataset_name_parts[4]
    data_type <- dataset_name_parts[5]
    data_type_label <- ifelse(
        data_type == "kraken", "Microbiome",
        ifelse(data_type == "htseq", "Expression", "Combo")
    )
    fig_num <- ifelse(
        data_type == "kraken", "1", ifelse(target == "os", "Ex1", "Ex2")
    )
    data <- data.frame(
        cnet_model_scores[[dataset_name]],
        cox_clinical_model_scores[[dataset_name]]
    )
    colnames(data) <- c(
        paste(data_type_label, "Clinical", sep=" + "), "Clinical"
    )
    # exclude rows with all NAs
    data <- data[rowSums(is.na(data)) != ncol(data), , drop=FALSE]
    p_adj <- potential_hits$p_adj[
        potential_hits$cancer == cancer
        & potential_hits$analysis == analysis
        & potential_hits$versus == target
        & potential_hits$features == data_type
    ]
    cat(p_adj, "\n")
    file <- paste(
        args$out_dir, paste0(dataset_name, ".", args$file_format), sep="/"
    )
    title <- paste(str_to_upper(cancer), str_to_upper(target))
    colors <- c(ifelse(data_type == "kraken", "#448ee4", "#c04e01"), "#6f828a")
    y_label <- "C-index"
    make_plot(file, data, p_adj, title, colors, y_label)
}

rfe_model_scores <- readRDS(
    paste(args$results_dir, "resp", "rfe_model_scores.rds", sep="/")
)
svm_clinical_model_scores <- readRDS(
    paste(args$results_dir, "resp", "svm_clinical_model_scores.rds", sep="/")
)
for (dataset_name in colnames(rfe_model_scores)) {
    cat(dataset_name, "\n")
    dataset_name_parts <- str_split(dataset_name, "_")[[1]]
    cancer <- dataset_name_parts[2]
    analysis <- dataset_name_parts[3]
    target <- dataset_name_parts[4]
    data_type <- dataset_name_parts[5]
    data_type_label <- ifelse(
        data_type == "kraken", "Microbiome",
        ifelse(data_type == "htseq", "Expression", "Combo")
    )
    fig_num <- ifelse(
        data_type == "kraken", "1", ifelse(target == "os", "Ex1", "Ex2")
    )
    data <- data.frame(
        rfe_model_scores[[dataset_name]],
        svm_clinical_model_scores[[dataset_name]]
    )
    colnames(data) <- c(
        paste(data_type_label, "Clinical", sep=" + "), "Clinical"
    )
    # exclude rows with all NAs
    data <- data[rowSums(is.na(data)) != ncol(data), , drop=FALSE]
    p_adj <- potential_hits$p_adj[
        potential_hits$cancer == cancer
        & potential_hits$analysis == analysis
        & potential_hits$versus == target
        & potential_hits$features == data_type
    ]
    cat(p_adj, "\n")
    file <- paste(
        args$out_dir, paste0(dataset_name, ".", args$file_format), sep="/"
    )
    title <- paste(str_to_upper(cancer), str_to_title(target))
    colors <- c(ifelse(data_type == "kraken", "#448ee4", "#c04e01"), "#6f828a")
    y_label <- "AUROC"
    make_plot(file, data, p_adj, title, colors, y_label)
}
