options(warn = 1)
suppressPackageStartupMessages({
    library(argparser)
    library(Biobase)
    library(data.table)
    library(GenomicDataCommons)
    library(stringr)
})

stopifnot(GenomicDataCommons::status()$status == "OK")

argp <- arg_parser("Create ExpressionSets")
argp <- add_argument(
    argp, "--data-dir",
    default = "data", help = "Data directory"
)
argp <- add_argument(argp, "--cancers", nargs = Inf, help = "TCGA cancer codes")
argp <- add_argument(argp, "--surv-types", nargs = Inf, help = "Survival types")
argp <- add_argument(argp, "--drug-names", nargs = Inf, help = "Drug names")
argp <- add_argument(
    argp, "--gdc-workflow-types",
    default = c("HTSeq - Counts"),
    help = "GDC workflow types"
)
argp <- add_argument(
    argp, "--save-data-matrix",
    flag = TRUE,
    help = "Save data matrix rds files"
)
args <- parse_args(argp)

cat("Loading k2b_kraken_data.rds\n")
kraken_data <- readRDS(
    paste(args$data_dir, "k2b_kraken_data.rds", sep = "/")
)
cat("Loading k2b_kraken_sample_meta.rds\n")
kraken_sample_meta <- readRDS(
    paste(args$data_dir, "k2b_kraken_sample_meta.rds", sep = "/")
)
cat("Loading k2b_kraken_feature_meta.rds\n")
kraken_feature_meta <- readRDS(
    paste(args$data_dir, "k2b_kraken_feature_meta.rds", sep = "/")
)
if (!identical(row.names(kraken_feature_meta), row.names(kraken_data))) {
    stop("Kraken data matrix and feature annots row names not identical")
}
cat("Loading response_pdata.rds\n")
response_pdata <- readRDS(paste(args$data_dir, "response_pdata.rds", sep = "/"))
cat("Loading survival_pdata.rds\n")
survival_pdata <- readRDS(paste(args$data_dir, "survival_pdata.rds", sep = "/"))

# generate datasets
min_uniq_cases <- 16
min_uniq_cases_per_class <- 4
min_uniq_case_exceptions <- c()

create_surv_eset <- function(
    adata, pdata, fdata, cancer, surv_type, data_type, msg_prefix,
    print_sample_msgs = TRUE) {
    eset <- ExpressionSet(
        assayData = adata, phenoData = AnnotatedDataFrame(pdata)
    )
    if (!is.null(fdata)) featureData(eset) <- AnnotatedDataFrame(fdata)
    status_col <- toupper(surv_type)
    time_col <- paste(status_col, "time", sep = "_")
    eset$Status <- eset[[status_col]]
    eset$Survival_in_days <- eset[[time_col]]
    eset <- eset[
        , !is.na(eset[[status_col]]) & !is.na(eset[[time_col]]) &
            eset[[time_col]] > 0
    ]
    eset <- eset[rowSums(exprs(eset)) != 0, ]
    eset <- eset[, colSums(exprs(eset)) != 0]
    pData(eset) <- gdata::drop.levels(pData(eset))
    if (anyDuplicated(eset$case_submitter_id)) {
        eset$Group <- match(
            eset$case_submitter_id, unique(eset$case_submitter_id)
        )
        if ("experimental_strategy" %in% colnames(pData(eset))) {
            case_strategy <- as.data.frame(table(data.frame(
                case_submitter_id = factor(eset$case_submitter_id),
                experimental_strategy = factor(eset$experimental_strategy)
            )))
            case_strategy <- case_strategy[case_strategy$Freq != 0, ]
            max_freq <- aggregate(Freq ~ case_submitter_id, case_strategy, max)
            colnames(max_freq)[colnames(max_freq) == "Freq"] <- "MaxFreq"
            case_strategy <-
                merge(case_strategy, max_freq, by = "case_submitter_id")
            case_strategy$GroupWeight <-
                case_strategy$MaxFreq / case_strategy$Freq
            eset$GroupWeight <- case_strategy$GroupWeight[match(
                paste(eset$case_submitter_id, eset$experimental_strategy),
                paste(
                    case_strategy$case_submitter_id,
                    case_strategy$experimental_strategy
                )
            )]
        } else if ("GroupWeight" %in% colnames(pData(eset))) {
            eset$GroupWeight <- NULL
        }
    } else {
        if ("Group" %in% colnames(pData(eset))) eset$Group <- NULL
        if ("GroupWeight" %in% colnames(pData(eset))) eset$GroupWeight <- NULL
    }
    cancer_target <- gsub("^tcga ", "", paste(
        tolower(gsub("(-|\\s)+", " ", c(cancer, surv_type))),
        collapse = " "
    ))
    if (
        length(unique(eset$case_submitter_id)) >= min_uniq_cases ||
            cancer_target %in% min_uniq_case_exceptions
    ) {
        if (print_sample_msgs) {
            cat(
                msg_prefix, toupper(surv_type),
                length(unique(eset$case_submitter_id)), "cases",
                ncol(eset), "samples", "\n"
            )
        }
        eset_filename_parts <- c(cancer, "surv", surv_type, data_type, "eset")
        eset_filename <- paste0(paste(
            tolower(gsub("(-|\\s)+", "_", eset_filename_parts)),
            collapse = "_"
        ), ".rds")
        cat(msg_prefix, "Creating", eset_filename, "\n")
        saveRDS(eset, paste(args$data_dir, eset_filename, sep = "/"))
    } else {
        cat(msg_prefix, paste("Skipping", toupper(surv_type)), "\n")
    }
}

create_drug_eset <- function(
    adata, pdata, fdata, cancer, resp_type, drug_name, data_type, msg_prefix,
    print_sample_msgs = TRUE) {
    eset <- ExpressionSet(
        assayData = adata, phenoData = AnnotatedDataFrame(pdata)
    )
    if (!is.null(fdata)) featureData(eset) <- AnnotatedDataFrame(fdata)
    eset$Class <- ifelse(
        eset$response %in% c(
            "Complete Response", "Partial Response"
        ), 1,
        ifelse(eset$response %in% c(
            "Clinical Progressive Disease", "Stable Disease"
        ), 0, NA)
    )
    eset$Class <- as.factor(eset$Class)
    eset <- eset[, !is.na(eset$Class)]
    eset <- eset[rowSums(exprs(eset)) != 0, ]
    eset <- eset[, colSums(exprs(eset)) != 0]
    pData(eset) <- gdata::drop.levels(pData(eset))
    if (anyDuplicated(eset$case_submitter_id)) {
        eset$Group <- match(
            eset$case_submitter_id, unique(eset$case_submitter_id)
        )
        num_uniq_cases_per_class <- aggregate(
            case_submitter_id ~ Class, pData(eset),
            function(x) length(unique(x))
        )$case_submitter_id
    } else {
        if ("Group" %in% colnames(pData(eset))) eset$Group <- NULL
        num_uniq_cases_per_class <- table(eset$Class)
    }
    cancer_target <- gsub("^tcga ", "", paste(
        tolower(gsub("(-|\\s)+", " ", c(cancer, drug_name))),
        collapse = " "
    ))
    if (
        (
            length(unique(eset$case_submitter_id)) >= min_uniq_cases ||
                cancer_target %in% min_uniq_case_exceptions
        ) &&
            length(unique(eset$Class)) > 1 &&
            min(num_uniq_cases_per_class) >= min_uniq_cases_per_class
    ) {
        if (print_sample_msgs) {
            cat(
                msg_prefix, drug_name, length(unique(eset$case_submitter_id)),
                "cases", ncol(eset), "samples\n"
            )
        }
        eset_filename_parts <-
            c(cancer, resp_type, drug_name, data_type, "eset")
        eset_filename <- paste0(paste(
            tolower(gsub("(-|\\s)+", "_", eset_filename_parts)),
            collapse = "_"
        ), ".rds")
        cat(msg_prefix, "Creating", eset_filename, "\n")
        saveRDS(eset, paste(args$data_dir, eset_filename, sep = "/"))
    }
}

get_gdc_data <- function(project_id, workflow_types, msg_prefix) {
    file_query <-
        files() %>%
        filter(
            cases.project.project_id %in% project_id &
                analysis.workflow_type %in% workflow_types
        ) %>%
        GenomicDataCommons::select(c(
            "file_name",
            "analysis.workflow_type",
            "cases.project.project_id",
            "cases.case_id",
            "cases.submitter_id",
            "cases.samples.sample_id",
            "cases.samples.submitter_id",
            "cases.samples.sample_type",
            "cases.samples.is_ffpe",
            "cases.samples.portions.is_ffpe",
            "cases.samples.portions.analytes.aliquots.aliquot_id",
            "cases.samples.portions.analytes.aliquots.submitter_id"
        ))
    file_results <- results_all(file_query)
    file_meta <- data.frame(
        file_id = file_results$file_id,
        file_name = file_results$file_name,
        workflow_type = file_results$analysis$workflow_type,
        project_id = vapply(
            sapply(file_results$cases, `[[`, "project"), `[`, "project_id"
        ),
        case_id = sapply(file_results$cases, `[[`, "case_id"),
        case_submitter_id = sapply(file_results$cases, `[[`, "submitter_id"),
        sample_id = sapply(
            sapply(file_results$cases, `[[`, "samples"), `[[`, "sample_id"
        ),
        sample_submitter_id = sapply(
            sapply(file_results$cases, `[[`, "samples"), `[[`, "submitter_id"
        ),
        sample_type = sapply(
            sapply(file_results$cases, `[[`, "samples"), `[[`, "sample_type"
        ),
        sample_is_ffpe = sapply(
            sapply(file_results$cases, `[[`, "samples"), `[[`, "is_ffpe"
        ),
        portion_is_ffpe = sapply(
            sapply(
                sapply(file_results$cases, `[[`, "samples"), `[[`, "portions"
            ), `[[`, "is_ffpe"
        ),
        aliquot_id = sapply(
            sapply(
                sapply(
                    sapply(
                        sapply(
                            file_results$cases, `[[`, "samples"
                        ), `[[`, "portions"
                    ), `[[`, "analytes"
                ), `[[`, "aliquots"
            ), `[[`, "aliquot_id"
        ),
        aliquot_submitter_id = sapply(
            sapply(
                sapply(
                    sapply(
                        sapply(
                            file_results$cases, `[[`, "samples"
                        ), `[[`, "portions"
                    ), `[[`, "analytes"
                ), `[[`, "aliquots"
            ), `[[`, "submitter_id"
        ),
        row.names = file_results$file_id,
        stringsAsFactors = FALSE
    )
    file_meta <- file_meta[order(
        file_meta$project_id, file_meta$case_submitter_id,
        file_meta$sample_submitter_id, file_meta$aliquot_submitter_id
    ), ]
    cat(msg_prefix, "Downloading", length(file_results$id), "GDC files\n")
    files <- gdcdata(ids(file_query), progress = FALSE, use_cached = TRUE)
    return(list(meta = file_meta, files = files))
}

gtf_annots_filename <- "gencode_v36_ensg_v102_annots.tsv"
gtf_annots_file <- paste(args$data_dir, gtf_annots_filename, sep = "/")
cat("Loading", gtf_annots_file, "\n")
rna_annots <- read.delim(
    gtf_annots_file,
    row.names = 1, stringsAsFactors = FALSE
)
rna_annots$Symbol[grepl(".", rna_annots$Symbol, fixed = TRUE)] <- ""
combo_annots <- rbind(data.frame(
    Symbol = rep("", ncol(kraken_data)), row.names = colnames(kraken_data)
), rna_annots)

gdc_meta_pdata_cols <- c(
    "project_id", "case_id", "case_submitter_id", "sample_id",
    "sample_submitter_id", "sample_is_ffpe", "portion_is_ffpe", "aliquot_id",
    "aliquot_submitter_id", "file_id"
)
if (any(is.na(args$cancers))) {
    cancers <- sort(union(response_pdata$cancer, survival_pdata$cancer))
} else {
    cancers <- sort(paste("TCGA", toupper(args$cancers), sep = "-"))
}
cancer_msg_pad <- max(str_length(cancers))
if (!any(is.na(args$surv_types))) {
    surv_types <- sort(tolower(args$surv_types))
} else if (any(is.na(args$drug_names))) {
    surv_types <- c("os", "pfi")
} else {
    surv_types <- c()
}
resp_types <- c("resp")
type_msg_pad <- max(str_length(c("kraken", "combo")))

uniq_rna_case_ids <- c()
uniq_rna_sample_ids <- c()

cat("Generating datasets\n")
for (cancer in cancers) {
    # kraken
    data_type <- "kraken"
    msg_prefix <- paste(
        "[",
        str_pad(cancer, cancer_msg_pad, side = "right"),
        str_pad("Survival", type_msg_pad, side = "right"),
        "]"
    )
    # survival
    kraken_surv_meta <- merge(
        survival_pdata[
            survival_pdata$cancer == cancer,
            !(colnames(survival_pdata) %in% c("project_id", "case_id")),
            drop = FALSE
        ],
        kraken_sample_meta,
        by = "case_submitter_id"
    )
    if (nrow(kraken_surv_meta) == 0) {
        cat(msg_prefix, paste("No data"), "\n")
        next
    }
    kraken_surv_meta$num_uniq_read_groups <- NULL
    kraken_surv_meta$read_length <- NULL
    kraken_surv_meta$is_paired_end <- NULL
    row.names(kraken_surv_meta) <- kraken_surv_meta$file_id
    kraken_surv_data <- kraken_data[, row.names(kraken_surv_meta), drop = FALSE]
    kraken_surv_meta <-
        kraken_surv_meta[order(row.names(kraken_surv_meta)), , drop = FALSE]
    kraken_surv_data <-
        kraken_surv_data[, order(colnames(kraken_surv_data)), drop = FALSE]
    for (surv_type in surv_types) {
        create_surv_eset(
            kraken_surv_data, kraken_surv_meta, kraken_feature_meta, cancer,
            surv_type, data_type, msg_prefix
        )
    }
}
for (cancer in cancers) {
    # response
    msg_prefix <- paste(
        "[",
        str_pad(cancer, cancer_msg_pad, side = "right"),
        str_pad("Response", type_msg_pad, side = "right"),
        "]"
    )
    if (!any(is.na(args$drug_names))) {
        drug_names <- sort(str_to_title(args$drug_names))
    } else if (any(is.na(args$surv_types))) {
        drug_names <- sort(
            unique(response_pdata$drug_name[response_pdata$cancer == cancer])
        )
    } else {
        drug_names <- c()
    }
    for (drug_name in drug_names) {
        kraken_drug_meta <- merge(
            response_pdata[
                response_pdata$cancer == cancer &
                    response_pdata$drug_name == str_to_upper(drug_name),
                !(colnames(response_pdata) %in% c("project_id", "case_id")),
                drop = FALSE
            ],
            kraken_sample_meta,
            by = "case_submitter_id"
        )
        if (nrow(kraken_drug_meta) == 0) next
        kraken_drug_meta$num_uniq_read_groups <- NULL
        kraken_drug_meta$read_length <- NULL
        kraken_drug_meta$is_paired_end <- NULL
        row.names(kraken_drug_meta) <- kraken_drug_meta$file_id
        kraken_drug_data <-
            kraken_data[, row.names(kraken_drug_meta), drop = FALSE]
        kraken_drug_meta <-
            kraken_drug_meta[order(row.names(kraken_drug_meta)), , drop = FALSE]
        kraken_drug_data <-
            kraken_drug_data[, order(colnames(kraken_drug_data)), drop = FALSE]
        for (idx in seq_along(resp_types)) {
            create_drug_eset(
                kraken_drug_data, kraken_drug_meta, kraken_feature_meta, cancer,
                resp_types[idx], drug_name, data_type, msg_prefix, idx == 1
            )
        }
    }
}
