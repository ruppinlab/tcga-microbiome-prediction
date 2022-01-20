# Authors: Leandro Cruz Hermida <hermidal@cs.umd.edu>
# License: BSD 3 clause

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
})

edger_filterbyexpr_mask <- function(
    X, y, sample_meta=NULL, model_batch=FALSE, is_classif=TRUE
) {
    dge <- DGEList(counts=t(X))
    if (
        model_batch && !is.null(sample_meta) &&
        length(unique(sample_meta$Batch)) > 1
    ) {
        sample_meta$Batch <- factor(sample_meta$Batch)
        if (is_classif) {
            sample_meta$Class <- factor(sample_meta$Class)
            design <- model.matrix(~Batch + Class, data=sample_meta)
        } else {
            design <- model.matrix(~Batch, data=sample_meta)
        }
    } else if (is_classif) {
        design <- model.matrix(~factor(y))
    } else {
        design <- NULL
    }
    return(filterByExpr(dge, design))
}

edger_feature_score <- function(
    X, y, sample_meta=NULL, lfc=0, scoring_meth="pv", robust=TRUE,
    model_batch=FALSE
) {
    counts <- t(X)
    if (
        model_batch && !is.null(sample_meta) &&
        length(unique(sample_meta$Batch)) > 1
    ) {
        sample_meta$Batch <- factor(sample_meta$Batch)
        sample_meta$Class <- factor(sample_meta$Class)
        design <- model.matrix(~Batch + Class, data=sample_meta)
    } else {
        design <- model.matrix(~factor(y))
    }
    dge <- DGEList(counts=counts)
    dge <- calcNormFactors(dge, method="TMM")
    dge <- estimateDisp(dge, design, robust=robust)
    fit <- glmQLFit(dge, design, robust=robust)
    if (lfc == 0) {
        glt <- glmQLFTest(fit, coef=ncol(design))
    } else {
        glt <- glmTreat(fit, coef=ncol(design), lfc=lfc)
    }
    results <- as.data.frame(topTags(
        glt, n=Inf, adjust.method="BH", sort.by="none"
    ))
    results <- results[order(as.integer(row.names(results))), , drop=FALSE]
    if (scoring_meth == "lfc_pv") {
        scores <- abs(results$logFC) * -log10(results$PValue)
    } else {
        scores <- results$PValue
    }
    return(list(scores, results$FDR))
}

limma_feature_score <- function(
    X, y, sample_meta=NULL, lfc=0, scoring_meth="pv", robust=FALSE, trend=FALSE,
    model_batch=FALSE
) {
    if (
        model_batch && !is.null(sample_meta) &&
        length(unique(sample_meta$Batch)) > 1
    ) {
        sample_meta$Batch <- factor(sample_meta$Batch)
        sample_meta$Class <- factor(sample_meta$Class)
        design <- model.matrix(~Batch + Class, data=sample_meta)
    } else {
        design <- model.matrix(~factor(y))
    }
    fit <- lmFit(t(X), design)
    fit <- treat(fit, lfc=lfc, robust=robust, trend=trend)
    results <- topTreat(
        fit, coef=ncol(design), number=Inf, adjust.method="BH", sort.by="none"
    )
    results <- results[order(as.integer(row.names(results))), , drop=FALSE]
    if (scoring_meth == "lfc_pv") {
        scores <- abs(results$logFC) * -log10(results$P.Value)
    } else {
        scores <- results$P.Value
    }
    return(list(scores, results$adj.P.Val))
}
