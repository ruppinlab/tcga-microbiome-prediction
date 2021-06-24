# Author: Leandro Cruz Hermida <hermidal@cs.umd.edu>
# License: BSD 3 clause

import os
import numpy as np

import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr
from sklearn.base import BaseEstimator
from sklearn.utils import check_array, check_X_y
from sklearn.utils.validation import check_is_fitted, check_memory

from ._base import ExtendedSelectorMixin

numpy2ri.deactivate()
pandas2ri.deactivate()
numpy2ri.activate()
pandas2ri.activate()

r_base = importr('base')
if 'edger_filterbyexpr_mask' not in robjects.globalenv:
    r_base.source(os.path.dirname(__file__) + '/_rna_seq.R')
r_edger_filterbyexpr_mask = robjects.globalenv['edger_filterbyexpr_mask']
r_edger_feature_score = robjects.globalenv['edger_feature_score']
r_limma_feature_score = robjects.globalenv['limma_feature_score']
if 'edger_tmm_logcpm_transform' not in robjects.globalenv:
    r_base.source(os.path.dirname(__file__) + '/../preprocessing/_rna_seq.R')
r_edger_tmm_logcpm_transform = robjects.globalenv['edger_tmm_logcpm_transform']


def edger_feature_score(X, y, sample_meta, lfc, robust, prior_count,
                        model_batch):
    pv, pa, xt, rs = r_edger_feature_score(
        X, y, sample_meta=sample_meta, lfc=lfc, robust=robust,
        prior_count=prior_count, model_batch=model_batch)
    return (np.array(pv, dtype=float), np.array(pa, dtype=float),
            np.array(xt, dtype=float), np.array(rs, dtype=float))


def edger_tmm_logcpm_transform(X, ref_sample, prior_count):
    return np.array(r_edger_tmm_logcpm_transform(
        X, ref_sample=ref_sample, prior_count=prior_count), dtype=float)


def limma_feature_score(X, y, sample_meta, lfc, robust, trend, model_batch):
    pv, pa = r_limma_feature_score(
        X, y, sample_meta=sample_meta, lfc=lfc, robust=robust,
        trend=trend, model_batch=model_batch)
    return (np.array(pv, dtype=float), np.array(pa, dtype=float))


class EdgeRFilterByExpr(ExtendedSelectorMixin, BaseEstimator):
    """edgeR filterByExpr feature selector for count data

    Parameters
    ----------
    model_batch : bool (default = False)
        Model batch effect if sample_meta passed to fit and Batch column
        exists.

    is_classif : bool (default = True)
        Whether this is a classification design.
    """

    def __init__(self, model_batch=False, is_classif=True):
        self.model_batch = model_batch
        self.is_classif = is_classif

    def fit(self, X, y, sample_meta=None):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Training counts data matrix.

        y : array-like, shape = (n_samples,)
            Training sample class labels.

        sample_meta : pandas.DataFrame, pandas.Series (default = None), \
            shape = (n_samples, n_metadata)
            Training sample metadata.

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y, dtype=int)
        if sample_meta is None:
            sample_meta = robjects.NULL
        self._mask = np.array(r_edger_filterbyexpr_mask(
            X, y, sample_meta=sample_meta, model_batch=self.model_batch,
            is_classif=self.is_classif), dtype=bool)
        return self

    def transform(self, X, sample_meta=None):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Input counts data matrix.

        sample_meta : Ignored.

        Returns
        -------
        Xr : array of shape (n_samples, n_selected_features)
            edgeR filterByExpr counts data matrix with only the selected
            features.
        """
        check_is_fitted(self, '_mask')
        X = check_array(X, dtype=int)
        return super().transform(X)

    def inverse_transform(self, X, sample_meta=None):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Input transformed data matrix.

        sample_meta : Ignored.

        Returns
        -------
        Xr : array of shape (n_samples, n_original_features)
            `X` with columns of zeros inserted where features would have
            been removed by :meth:`transform`.
        """
        raise NotImplementedError('inverse_transform not implemented.')

    def _more_tags(self):
        return {'requires_positive_X': True}

    def _get_support_mask(self):
        check_is_fitted(self, '_mask')
        return self._mask


class EdgeR(ExtendedSelectorMixin, BaseEstimator):
    """edgeR differential expression feature selector and
    normalizer/transformer for RNA-seq count data

    Parameters
    ----------
    k : int or "all" (default = "all")
        Number of top features to select. Specifying k = "all" and pv = 1.0
        bypasses selection, for use in a parameter search. When pv is also
        specified then returns the intersection of both parameter results.

    pv : float (default = 1.0)
        Select top features below an adjusted p-value threshold. Specifying
        k = "all" and pv = 1.0 bypasses selection, for use in a parameter
        search. When k is also specified returns the intersection of both
        parameter results.

    fc : float (default = 1.0)
        glmTreat absolute fold change minimum threshold. Default value of 1.0
        gives glmQLFTest results.

    robust : bool (default = True)
        estimateDisp and glmQLFit robust option.

    prior_count : int (default = 1)
        Average count to add to each observation to avoid taking log of zero.
        Larger values for produce stronger moderation of the values for low
        counts and more shrinkage of the corresponding log fold changes.

    model_batch : bool (default = False)
        Model batch effect if sample_meta passed to fit and Batch column
        exists.

    memory : None, str or object with the joblib.Memory interface \
        (default = None)
        Used for internal caching. By default, no caching is done.
        If a string is given, it is the path to the caching directory.

    Attributes
    ----------
    pvals_ : array, shape (n_features,)
        Feature raw p-values.

    padjs_ : array, shape (n_features,)
        Feature adjusted p-values.

    ref_sample_ : array, shape (n_features,)
        TMM normalization reference sample feature vector.
    """

    def __init__(self, k='all', pv=1, fc=1, robust=True, prior_count=1,
                 model_batch=False, memory=None):
        self.k = k
        self.pv = pv
        self.fc = fc
        self.robust = robust
        self.prior_count = prior_count
        self.model_batch = model_batch
        self.memory = memory

    def fit(self, X, y, sample_meta=None):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Training counts data matrix.

        y : array-like, shape = (n_samples,)
            Training sample class labels.

        sample_meta : pandas.DataFrame, pandas.Series (default = None), \
            shape = (n_samples, n_metadata)
            Training sample metadata.

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y, dtype=int)
        self._check_params(X, y)
        memory = check_memory(self.memory)
        if sample_meta is None:
            sample_meta = robjects.NULL
        self.pvals_, self.padjs_, self._log_cpms, self.ref_sample_ = (
            memory.cache(edger_feature_score)(
                X, y, sample_meta=sample_meta, lfc=np.log2(self.fc),
                robust=self.robust, prior_count=self.prior_count,
                model_batch=self.model_batch))
        return self

    def transform(self, X, sample_meta=None):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Input counts data matrix.

        sample_meta : Ignored.

        Returns
        -------
        Xr : array of shape (n_samples, n_selected_features)
            edgeR TMM normalized log-CPM transformed data matrix with only the
            selected features.
        """
        check_is_fitted(self, '_log_cpms')
        X = check_array(X, dtype=int)
        if hasattr(self, '_train_done'):
            memory = check_memory(self.memory)
            X = memory.cache(edger_tmm_logcpm_transform)(
                X, ref_sample=self.ref_sample_, prior_count=self.prior_count)
        else:
            X = self._log_cpms
            self._train_done = True
        return super().transform(X)

    def inverse_transform(self, X, sample_meta=None):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Input transformed data matrix.

        sample_meta : Ignored.

        Returns
        -------
        Xr : array of shape (n_samples, n_original_features)
            `X` with columns of zeros inserted where features would have
            been removed by :meth:`transform`.
        """
        raise NotImplementedError('inverse_transform not implemented.')

    def _more_tags(self):
        return {'requires_positive_X': True}

    def _check_params(self, X, y):
        if not (self.k == 'all' or 0 <= self.k <= X.shape[1]):
            raise ValueError(
                "k should be 0 <= k <= n_features; got %r."
                "Use k='all' to return all features." % self.k)
        if not 0 <= self.pv <= 1:
            raise ValueError('pv should be 0 <= pv <= 1; got %r.' % self.pv)
        if self.fc < 1:
            raise ValueError(
                'fold change threshold should be >= 1; got %r.' % self.fc)

    def _get_support_mask(self):
        check_is_fitted(self, 'pvals_')
        mask = np.zeros_like(self.pvals_, dtype=bool)
        if self.pv > 0:
            if self.k == 'all':
                mask = np.ones_like(self.pvals_, dtype=bool)
                if self.pv < 1:
                    mask[self.padjs_ > self.pv] = False
            elif self.k > 0:
                mask[np.argsort(self.pvals_, kind='mergesort')[:self.k]] = True
                if self.pv < 1:
                    mask[self.padjs_ > self.pv] = False
        return mask


class Limma(ExtendedSelectorMixin, BaseEstimator):
    """limma differential expression feature selector for gene expression data

    Parameters
    ----------
    k : int or "all" (default = "all")
        Number of top features to select. Specifying k = "all" and pv = 1.0
        bypasses selection, for use in a parameter search. When pv is also
        specified then returns the intersection of both parameter results.

    pv : float (default = 1.0)
        Select top features below an adjusted p-value threshold. Specifying
        k = "all" and pv = 1.0 bypasses selection, for use in a parameter
        search. When k is also specified returns the intersection of both
        parameter results.

    fc : float (default = 1.0)
        treat absolute fold change minimum threshold. Default value of 1.0
        gives eBayes results.

    robust : bool (default = False)
        limma treat/eBayes robust option.

    trend : bool (default = False)
        limma treat/eBayes trend option.

    model_batch : bool (default = False)
        Model batch effect if sample_meta passed to fit and Batch column
        exists.

    memory : None, str or object with the joblib.Memory interface \
        (default = None)
        Used for internal caching. By default, no caching is done.
        If a string is given, it is the path to the caching directory.

    Attributes
    ----------
    pvals_ : array, shape (n_features,)
        Feature raw p-values.

    padjs_ : array, shape (n_features,)
        Feature adjusted p-values.
    """

    def __init__(self, k='all', pv=1, fc=1, robust=False, trend=False,
                 model_batch=False, memory=None):
        self.k = k
        self.pv = pv
        self.fc = fc
        self.robust = robust
        self.trend = trend
        self.model_batch = model_batch
        self.memory = memory

    def fit(self, X, y, sample_meta=None):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Training gene expression data matrix.

        y : array-like, shape = (n_samples,)
            Training sample class labels.

        sample_meta : pandas.DataFrame, pandas.Series (default = None), \
            shape = (n_samples, n_metadata)
            Training sample metadata.

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y)
        self._check_params(X, y)
        memory = check_memory(self.memory)
        if sample_meta is None:
            sample_meta = robjects.NULL
        self.pvals_, self.padjs_ = memory.cache(limma_feature_score)(
            X, y, sample_meta=sample_meta, lfc=np.log2(self.fc),
            robust=self.robust, trend=self.trend, model_batch=self.model_batch)
        return self

    def transform(self, X, sample_meta=None):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Gene expression data matrix.

        sample_meta : Ignored.

        Returns
        -------
        Xr : array of shape (n_samples, n_selected_features)
            Gene expression data matrix with only the selected features.
        """
        check_is_fitted(self, 'pvals_')
        return super().transform(X)

    def inverse_transform(self, X, sample_meta=None):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Input transformed data matrix.

        sample_meta : Ignored.

        Returns
        -------
        Xr : array of shape (n_samples, n_original_features)
            `X` with columns of zeros inserted where features would have
            been removed by :meth:`transform`.
        """
        raise NotImplementedError('inverse_transform not implemented.')

    def _check_params(self, X, y):
        if not (self.k == 'all' or 0 <= self.k <= X.shape[1]):
            raise ValueError(
                "k should be 0 <= k <= n_features; got %r."
                "Use k='all' to return all features." % self.k)
        if not 0 <= self.pv <= 1:
            raise ValueError('pv should be 0 <= pv <= 1; got %r.' % self.pv)
        if self.fc < 1:
            raise ValueError(
                'fold change threshold should be >= 1; got %r.' % self.fc)

    def _get_support_mask(self):
        check_is_fitted(self, 'pvals_')
        mask = np.zeros_like(self.pvals_, dtype=bool)
        if self.pv > 0:
            if self.k == 'all':
                mask = np.ones_like(self.pvals_, dtype=bool)
                if self.pv < 1:
                    mask[self.padjs_ > self.pv] = False
            elif self.k > 0:
                mask[np.argsort(self.pvals_, kind='mergesort')[:self.k]] = True
                if self.pv < 1:
                    mask[self.padjs_ > self.pv] = False
        return mask
