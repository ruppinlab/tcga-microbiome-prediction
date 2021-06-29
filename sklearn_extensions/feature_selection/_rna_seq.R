# Authors: Leandro Cruz Hermida <hermidal@cs.umd.edu>
# License: BSD 3 clause

suppressPackageStartupMessages(library("edgeR"))


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

# adapted from edgeR::calcNormFactors source code
edger_tmm_ref_column <- function(counts, lib.size=colSums(counts), p=0.75) {
    y <- t(t(counts) / lib.size)
    f <- apply(y, 2, function(x) quantile(x, p=p))
    ref_column <- which.min(abs(f - mean(f)))
}

edger_feature_score <- function(
    X, y, sample_meta=NULL, lfc=0, robust=TRUE, prior_count=1, model_batch=FALSE
) {
    suppressPackageStartupMessages(library("edgeR"))
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
    log_cpm <- cpm(dge, log=TRUE, prior.count=prior_count)
    ref_sample <- counts[, edger_tmm_ref_column(counts)]
    return(list(results$PValue, results$FDR, t(log_cpm), ref_sample))
}

limma_feature_score <- function(
    X, y, sample_meta=NULL, lfc=0, robust=FALSE, trend=FALSE, model_batch=FALSE
) {
    suppressPackageStartupMessages(library("limma"))
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
    return(list(results$P.Value, results$adj.P.Val))
}
