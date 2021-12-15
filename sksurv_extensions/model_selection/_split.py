# Authors: Leandro Hermida <hermidal@cs.umd.edu>
#
# License: BSD 3 clause

from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit
from sklearn.model_selection._split import _RepeatedSplits

from sklearn_extensions.model_selection import (
    StratifiedSampleFromGroupShuffleSplit)


class SurvivalStratifiedKFold(StratifiedKFold):
    """Survival Stratified K-Folds cross-validator

    Provides train/test indices to split data in train/test sets.

    This cross-validation object is a variation of StratifiedKFold that returns
    stratified folds for a survival structured array. The folds are made by
    preserving the percentage of samples for each survival status.

    Parameters
    ----------
    n_splits : int, default=5
        Number of folds. Must be at least 2.

    shuffle : boolean, optional
        Whether to shuffle each status's samples before splitting into batches.

    random_state : int, RandomState instance or None, optional, default=None
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`. Only used when ``shuffle`` is True. This should be left
        to None if ``shuffle`` is False.

    Notes
    -----
    The implementation is designed to:

    * Generate test sets such that all contain the same distribution of
      ssurvival status, or as close as possible.
    * Be invariant to status label: relabelling ``y = ["Happy", "Sad"]`` to
      ``y = [1, 0]`` should not change the indices generated.
    * Preserve order dependencies in the dataset ordering, when
      ``shuffle=False``: all samples from status k in some test set were
      contiguous in y, or separated in y by samples from statuses other than k.
    * Generate test sets where the smallest and largest differ by at most one
      sample.
    """

    def split(self, X, y, groups=None):
        y = y[y.dtype.names[0]].astype(int)
        return super().split(X, y, groups)


class RepeatedSurvivalStratifiedKFold(_RepeatedSplits):
    """Repeated Survival Stratified K-Fold cross validator.

    Repeats Stratified K-Fold with non-overlapping groups n times with
    different randomization in each repetition.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n_splits : int, default=5
        Number of folds. Must be at least 2.

    n_repeats : int, default=10
        Number of times cross-validator needs to be repeated.

    random_state : int or RandomState instance, default=None
        Controls the generation of the random states for each repetition.
        Pass an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Notes
    -----
    Randomized CV splitters may return different results for each call of
    split. You can make the results identical by setting `random_state`
    to an integer.
    """

    def __init__(self, n_splits=5, n_repeats=10, random_state=None):
        super().__init__(SurvivalStratifiedKFold, n_splits=n_splits,
                         n_repeats=n_repeats, random_state=random_state)


class SurvivalStratifiedShuffleSplit(StratifiedShuffleSplit):

    def split(self, X, y, groups=None):
        y = y[y.dtype.names[0]].astype(int)
        return super().split(X, y, groups)


class SurvivalStratifiedSampleFromGroupShuffleSplit(
        StratifiedSampleFromGroupShuffleSplit):

    def split(self, X, y, groups, weights=None):
        y = y[y.dtype.names[0]].astype(int)
        return super().split(X, y, groups, weights)
