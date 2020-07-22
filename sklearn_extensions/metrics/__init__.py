"""
The :mod:`sklearn_extensions.metrics` module includes score functions,
performance metrics and pairwise metrics and distance computations.
"""

from ._scorer import check_scoring
from ._scorer import make_scorer
from ._scorer import SCORERS
from ._scorer import get_scorer


__all__ = ['check_scoring',
           'get_scorer',
           'make_scorer',
           'SCORERS']
