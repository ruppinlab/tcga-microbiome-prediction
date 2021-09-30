# Author: Leandro Cruz Hermida <hermidal@cs.umd.edu>
# License: BSD 3 clause

import os

import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr
from sklearn.base import BaseEstimator
from sklearn.utils import check_array
from sklearn.utils.validation import check_is_fitted

from ..base import ExtendedTransformerMixin

numpy2ri.deactivate()
pandas2ri.deactivate()
numpy2ri.activate()
pandas2ri.activate()

if 'edger_tmm_fit' not in robjects.globalenv:
    r_base = importr('base')
    r_base.source(os.path.dirname(__file__) + '/_rna_seq.R')
r_edger_tmm_fit = robjects.globalenv['edger_tmm_fit']
r_edger_tmm_logcpm_transform = robjects.globalenv['edger_tmm_logcpm_transform']


class EdgeRTMMLogCPM(ExtendedTransformerMixin, BaseEstimator):
    """edgeR TMM normalization and log-CPM transformation for count data

    Parameters
    ----------
    prior_count : int (default = 1)
        Average count to add to each observation to avoid taking log of zero.
        Larger values for produce stronger moderation of the values for low
        counts and more shrinkage of the corresponding log fold changes.

    Attributes
    ----------
    ref_sample_ : array, shape (n_features,)
        TMM normalization reference sample feature vector.
    """

    def __init__(self, prior_count=1):
        self.prior_count = prior_count

    def fit(self, X, y=None, sample_meta=None):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Input counts data matrix.

        y : ignored

        sample_meta: ignored
        """
        X = check_array(X, dtype=int)
        self.ref_sample_ = np.array(r_edger_tmm_fit(X), dtype=int)
        return self

    def transform(self, X, sample_meta=None):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Input counts data matrix.

        sample_meta : ignored

        Returns
        -------
        Xt : array of shape (n_samples, n_features)
            edgeR TMM normalized log-CPM transformed data matrix.
        """
        check_is_fitted(self, 'ref_sample_')
        X = check_array(X, dtype=int)
        X = np.array(r_edger_tmm_logcpm_transform(
            X, ref_sample=self.ref_sample_, prior_count=self.prior_count),
                     dtype=float)
        return X

    def inverse_transform(self, X, sample_meta=None):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Input transformed data matrix.

        sample_meta : ignored

        Returns
        -------
        Xr : array of shape (n_samples, n_original_features)
        """
        raise NotImplementedError('inverse_transform not implemented.')

    def _more_tags(self):
        return {'requires_positive_X': True}
