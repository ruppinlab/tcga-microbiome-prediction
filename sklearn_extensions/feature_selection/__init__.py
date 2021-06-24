"""
The :mod:`sklearn_extensions.feature_selection` module implements feature
selection algorithms. It currently includes univariate filter selection
methods and the recursive feature elimination algorithm.
"""

from ._from_model import SelectFromModel
from ._rna_seq import EdgeR, EdgeRFilterByExpr, Limma
from ._rfe import ExtendedRFE


__all__ = ['EdgeR',
           'EdgeRFilterByExpr',
           'ExtendedRFE',
           'Limma',
           'SelectFromModel']
