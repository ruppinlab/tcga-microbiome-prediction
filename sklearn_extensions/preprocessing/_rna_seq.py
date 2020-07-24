# Author: Leandro Cruz Hermida <hermidal@cs.umd.edu>

import os

import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr
from sklearn.base import BaseEstimator
from sklearn.utils import check_array
from sklearn.utils.validation import check_is_fitted, check_memory

from ..base import ExtendedTransformerMixin

numpy2ri.deactivate()
pandas2ri.deactivate()
numpy2ri.activate()
pandas2ri.activate()

if 'edger_tmm_logcpm_fit' not in robjects.globalenv:
    r_base = importr('base')
    r_base.source(os.path.dirname(__file__) + '/_rna_seq.R')
r_edger_tmm_logcpm_fit = robjects.globalenv['edger_tmm_logcpm_fit']
r_edger_tmm_logcpm_transform = robjects.globalenv['edger_tmm_logcpm_transform']


def edger_tmm_logcpm_fit(X, prior_count):
    xt, rs = r_edger_tmm_logcpm_fit(X, prior_count=prior_count)
    return np.array(xt, dtype=float), np.array(rs, dtype=float)


def edger_tmm_logcpm_transform(X, ref_sample, prior_count):
    return np.array(r_edger_tmm_logcpm_transform(
        X, ref_sample=ref_sample, prior_count=prior_count), dtype=float)


class EdgeRTMMLogCPM(ExtendedTransformerMixin, BaseEstimator):
    """edgeR TMM normalization and log-CPM transformation for count data

    Parameters
    ----------
    prior_count : int (default = 1)
        Average count to add to each observation to avoid taking log of zero.
        Larger values for produce stronger moderation of the values for low
        counts and more shrinkage of the corresponding log fold changes.

    memory : None, str or object with the joblib.Memory interface \
        (default = None)
        Used for internal caching. By default, no caching is done.
        If a string is given, it is the path to the caching directory.

    Attributes
    ----------
    ref_sample_ : array, shape (n_features,)
        TMM normalization reference sample feature vector.
    """

    def __init__(self, prior_count=1, memory=None):
        self.prior_count = prior_count
        self.memory = memory

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
        memory = check_memory(self.memory)
        self._log_cpms, self.ref_sample_ = memory.cache(edger_tmm_logcpm_fit)(
            X, prior_count=self.prior_count)
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
        check_is_fitted(self, '_log_cpms')
        X = check_array(X, dtype=int)
        if hasattr(self, '_train_done'):
            memory = check_memory(self.memory)
            X = memory.cache(edger_tmm_logcpm_transform)(
                X, ref_sample=self.ref_sample_, prior_count=self.prior_count)
        else:
            X = self._log_cpms
            self._train_done = True
        return X

    def inverse_transform(self, X, sample_meta=None):
        """
        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Input transformed data matrix.

        sample_meta: ignored

        Returns
        -------
        Xr : array of shape (n_samples, n_original_features)
        """
        raise NotImplementedError('inverse_transform not implemented.')

    def _more_tags(self):
        return {'requires_positive_X': True}
