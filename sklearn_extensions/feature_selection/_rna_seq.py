# Author: Leandro Cruz Hermida <hermidal@cs.umd.edu>

import os
import numpy as np

import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr
from sklearn.base import BaseEstimator
from sklearn.utils import check_array, check_X_y
from sklearn.utils.validation import check_is_fitted

from ._base import ExtendedSelectorMixin

numpy2ri.deactivate()
pandas2ri.deactivate()
numpy2ri.activate()
pandas2ri.activate()

r_base = importr('base')
if 'edger_filterbyexpr_mask' not in robjects.globalenv:
    r_base.source(os.path.dirname(__file__) + '/_rna_seq.R')
r_edger_filterbyexpr_mask = robjects.globalenv['edger_filterbyexpr_mask']


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
