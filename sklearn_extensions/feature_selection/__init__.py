"""
The :mod:`sklearn_extensions.feature_selection` module implements feature
selection algorithms. It currently includes univariate filter selection
methods and the recursive feature elimination algorithm.
"""

from ._rna_seq import EdgeRFilterByExpr
from ._rfe import ExtendedRFE


__all__ = ['EdgeRFilterByExpr',
           'ExtendedRFE']
