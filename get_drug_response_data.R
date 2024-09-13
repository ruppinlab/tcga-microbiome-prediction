options(warn = 1)
suppressPackageStartupMessages({
    library(GenomicDataCommons)
    library(stringr)
    library(TCGAbiolinks)
})

stopifnot(GenomicDataCommons::status()$status == "OK")

project_ids <- sort(
    projects() %>%
    GenomicDataCommons::filter(program.name == "TCGA") %>%
    ids()
)
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
stopifnot(length(unique(drug_names$tcga_name)), nrow(drug_names))

drug_response_data <- data.frame()
for (project_id in project_ids) {
    cat(paste("Processsing", project_id), "\n")
    drug_df_name <- paste(
        "clinical_drug",
        tolower(str_split_fixed(project_id, "-", n = 2)[2]),
        sep = "_"
    )
    # DLBC and LAML have no drug response clinical data
    if (project_id %in% c("TCGA-DLBC", "TCGA-LAML")) next
    stopifnot(drug_df_name %in% names(clinical_data))
    gdc_drug_data <- as.data.frame(clinical_data[[drug_df_name]])
    stopifnot(
        gdc_drug_data[1, 1] == "bcr_patient_uuid" &&
        gdc_drug_data[2, 1] == "CDE_ID:"
    )
    gdc_drug_data <- gdc_drug_data[3:nrow(gdc_drug_data), ]
    gdc_drug_data$pharmaceutical_therapy_drug_name <- toupper(
        gdc_drug_data$pharmaceutical_therapy_drug_name
    )

    missing_drug_names <- setdiff(
        gdc_drug_data$pharmaceutical_therapy_drug_name, drug_names$tcga_name
    )
    if (length(missing_drug_names) > 0) {
        cat(paste(project_id, "missing drug names:", missing_drug_names), "\n")
    }

    drug_response_data <- rbind(
        drug_response_data,
        data.frame(
            project_id = rep(project_id, nrow(gdc_drug_data)),
            case_submitter_id = gdc_drug_data$bcr_patient_barcode,
            drug_name = drug_names$standard_name[match(
                gdc_drug_data$pharmaceutical_therapy_drug_name,
                drug_names$tcga_name
            )],
            start_time = gdc_drug_data$pharmaceutical_tx_started_days_to,
            end_time = gdc_drug_data$pharmaceutical_tx_ended_days_to,
            response = gdc_drug_data$treatment_best_response
        )
    )
}

drug_response_data <- drug_response_data[!(
    drug_response_data$response %in%
    c("[Unknown]", "[Not Applicable]", "[Not Available]", "[Discrepancy]"
)), ]
write.table(
    drug_response_data, "data/tcga_drug_response.tsv",
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
)
