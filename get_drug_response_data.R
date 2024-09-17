options(warn = 1)
suppressPackageStartupMessages({
    library(dplyr)
    library(GenomicDataCommons)
    library(readr)
    library(stringr)
    library(TCGAbiolinks)
})

stopifnot(GenomicDataCommons::status()$status == "OK")

project_ids <-
    projects() %>%
    GenomicDataCommons::filter(program.name == "TCGA") %>%
    ids() %>%
    sort()

gdc_query <- GDCquery(
    project = project_ids,
    data.category = "Clinical",
    data.type = "Clinical Supplement",
    data.format = "BCR Biotab"
)
GDCdownload(gdc_query)
suppressMessages(clinical_data <- GDCprepare(gdc_query))

drug_names <- read.delim("data/tcga_drug_names.tsv", sep = "\t")
drug_names$tcga_name <- toupper(drug_names$tcga_name)
stopifnot(!any(duplicated(drug_names$tcga_name)))

drug_response_dfs <- list()
for (project_id in project_ids) {
    cat(paste("Processsing", project_id), "\n")
    # DLBC and LAML have no drug response clinical data
    if (project_id %in% c("TCGA-DLBC", "TCGA-LAML")) next
    cancer <- str_split_fixed(project_id, "-", n = 2)[2]
    drug_df_name <- paste("clinical_drug", tolower(cancer), sep = "_")
    stopifnot(drug_df_name %in% names(clinical_data))
    gdc_drug_data <- as.data.frame(clinical_data[[drug_df_name]])
    stopifnot(
        gdc_drug_data[1, 1] == "bcr_patient_uuid" &&
            gdc_drug_data[2, 1] == "CDE_ID:"
    )
    gdc_drug_data <- gdc_drug_data[3:nrow(gdc_drug_data), , drop = FALSE]
    gdc_drug_data$pharmaceutical_therapy_drug_name <- toupper(
        gdc_drug_data$pharmaceutical_therapy_drug_name
    )
    missing_drug_names <- setdiff(
        gdc_drug_data$pharmaceutical_therapy_drug_name, drug_names$tcga_name
    )
    if (length(missing_drug_names) > 0) {
        cat(paste("Missing drug names:", missing_drug_names), "\n")
    }
    drug_response_dfs[[project_id]] <-
        inner_join(
            gdc_drug_data, drug_names,
            by = join_by(pharmaceutical_therapy_drug_name == tcga_name)
        ) %>%
        mutate(cancer = cancer, .before = everything()) %>%
        dplyr::select(
            cancer,
            bcr_patient_barcode,
            standard_name,
            pharmaceutical_tx_started_days_to,
            pharmaceutical_tx_ended_days_to,
            treatment_best_response
        ) %>%
        rename(
            case_submitter_id = bcr_patient_barcode,
            drug_name = standard_name,
            start_time = pharmaceutical_tx_started_days_to,
            end_time = pharmaceutical_tx_ended_days_to,
            response = treatment_best_response
        )
}

drug_response_data <-
    bind_rows(drug_response_dfs) %>%
    dplyr::filter(!(response %in% c(
        "[Unknown]", "[Not Applicable]", "[Not Available]", "[Discrepancy]"
    ))) %>%
    write_tsv(
        "data/tcga_drug_response.tsv",
        col_names = TRUE, progress = FALSE, quote = "none"
    )
