# Author: Leandro Cruz Hermida <hermidal@cs.umd.edu>
# License: BSD 3 clause

suppressPackageStartupMessages(library("edgeR"))


# adapted from edgeR::calcNormFactors source code
edger_tmm_ref_column <- function(counts, lib.size=colSums(counts), p=0.75) {
    y <- t(t(counts) / lib.size)
    f <- apply(y, 2, function(x) quantile(x, p=p))
    ref_column <- which.min(abs(f - mean(f)))
}

edger_tmm_logcpm_fit <- function(X, prior_count=1) {
    counts <- t(X)
    dge <- DGEList(counts=counts)
    dge <- calcNormFactors(dge, method="TMM")
    log_cpm <- cpm(dge, log=TRUE, prior.count=prior_count)
    ref_sample <- counts[, edger_tmm_ref_column(counts)]
    return(list(t(log_cpm), ref_sample))
}

edger_tmm_logcpm_transform <- function(X, ref_sample, prior_count=1) {
    counts <- t(X)
    counts <- cbind(counts, ref_sample)
    colnames(counts) <- NULL
    dge <- DGEList(counts=counts)
    dge <- calcNormFactors(dge, method="TMM", refColumn=ncol(dge))
    log_cpm <- cpm(dge, log=TRUE, prior.count=prior_count)
    log_cpm <- log_cpm[, -ncol(log_cpm)]
    return(t(log_cpm))
}
