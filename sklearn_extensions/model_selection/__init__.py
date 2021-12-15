"""
The :mod:`sklearn_extensions.model_selection` module.
"""

from ._split import (StratifiedGroupKFold, RepeatedStratifiedGroupKFold,
                     StratifiedSampleFromGroupShuffleSplit)
from ._search import ExtendedGridSearchCV
from ._validation import shuffle_y


__all__ = ['ExtendedGridSearchCV',
           'StratifiedGroupKFold',
           'StratifiedSampleFromGroupShuffleSplit',
           'RepeatedStratifiedGroupKFold',
           'shuffle_y']
