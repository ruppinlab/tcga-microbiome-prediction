# Authors: Leandro Cruz Hermida <hermidal@cs.umd.edu>

from collections import Counter, defaultdict

import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection._split import (_BaseKFold, _RepeatedSplits,
                                            _validate_shuffle_split)
from sklearn.utils import _approximate_mode
from sklearn.utils.validation import (check_array, check_random_state,
                                      indexable, _num_samples)


class StratifiedGroupKFold(_BaseKFold):
    """Stratified K-Folds iterator variant with non-overlapping groups.

    This cross-validation object is a variation of StratifiedKFold that returns
    stratified folds with non-overlapping groups. The folds are made by
    preserving the percentage of samples for each class.

    The same group will not appear in two different folds (the number of
    distinct groups has to be at least equal to the number of folds).

    The difference between GroupKFold and StratifiedGroupKFold is that
    the former attempts to create balanced folds such that the number of
    distinct groups is approximately the same in each fold, whereas
    StratifiedGroupKFold attempts to create folds which preserve the
    percentage of samples for each class.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n_splits : int, default=5
        Number of folds. Must be at least 2.

    shuffle : bool, default=False
        Whether to shuffle each class's samples before splitting into batches.
        Note that the samples within each split will not be shuffled.

    random_state : int or RandomState instance, default=None
        When `shuffle` is True, `random_state` affects the ordering of the
        indices, which controls the randomness of each fold for each class.
        Otherwise, leave `random_state` as `None`.
        Pass an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import StratifiedGroupKFold
    >>> X = np.ones((17, 2))
    >>> y = np.array([0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    >>> groups = np.array([1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5, 5, 6, 6, 7, 8, 8])
    >>> cv = StratifiedGroupKFold(n_splits=3)
    >>> for train_idxs, test_idxs in cv.split(X, y, groups):
    ...     print("TRAIN:", groups[train_idxs])
    ...     print("      ", y[train_idxs])
    ...     print(" TEST:", groups[test_idxs])
    ...     print("      ", y[test_idxs])
    TRAIN: [2 2 4 5 5 5 5 6 6 7]
           [1 1 1 0 0 0 0 0 0 0]
     TEST: [1 1 3 3 3 8 8]
           [0 0 1 1 1 0 0]
    TRAIN: [1 1 3 3 3 4 5 5 5 5 8 8]
           [0 0 1 1 1 1 0 0 0 0 0 0]
     TEST: [2 2 6 6 7]
           [1 1 0 0 0]
    TRAIN: [1 1 2 2 3 3 3 6 6 7 8 8]
           [0 0 1 1 1 1 1 0 0 0 0 0]
     TEST: [4 5 5 5 5]
           [1 0 0 0 0]

    See also
    --------
    StratifiedKFold: Takes class information into account to build folds which
        retain class distributions (for binary or multiclass classification
        tasks).

    GroupKFold: K-fold iterator variant with non-overlapping groups.
    """

    def __init__(self, n_splits=5, shuffle=False, random_state=None):
        super().__init__(n_splits=n_splits, shuffle=shuffle,
                         random_state=random_state)

    # Implementation based on this kaggle kernel:
    # https://www.kaggle.com/jakubwasikowski/stratified-group-k-fold-cross-validation
    def _iter_test_indices(self, X, y, groups):
        labels_num = np.max(y) + 1
        y_counts_per_group = defaultdict(lambda: np.zeros(labels_num))
        y_distr = Counter()
        for label, group in zip(y, groups):
            y_counts_per_group[group][label] += 1
            y_distr[label] += 1

        y_counts_per_fold = defaultdict(lambda: np.zeros(labels_num))
        groups_per_fold = defaultdict(set)

        groups_and_y_counts = list(y_counts_per_group.items())
        rng = check_random_state(self.random_state)
        if self.shuffle:
            rng.shuffle(groups_and_y_counts)

        for group, y_counts in sorted(groups_and_y_counts,
                                      key=lambda x: -np.std(x[1])):
            best_fold = None
            min_eval = None
            for i in range(self.n_splits):
                y_counts_per_fold[i] += y_counts
                std_per_label = []
                for label in range(labels_num):
                    std_per_label.append(np.std(
                        [y_counts_per_fold[j][label] / y_distr[label]
                         for j in range(self.n_splits)]))
                y_counts_per_fold[i] -= y_counts
                fold_eval = np.mean(std_per_label)
                if min_eval is None or fold_eval < min_eval:
                    min_eval = fold_eval
                    best_fold = i
            y_counts_per_fold[best_fold] += y_counts
            groups_per_fold[best_fold].add(group)

        for i in range(self.n_splits):
            test_indices = [idx for idx, group in enumerate(groups)
                            if group in groups_per_fold[i]]
            yield test_indices


class RepeatedStratifiedGroupKFold(_RepeatedSplits):
    """Repeated Stratified Group K-Fold cross validator.

    Repeats Stratified Group K-Fold with non-overlapping groups n times with
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

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import RepeatedStratifiedGroupKFold
    >>> X = np.ones((17, 2))
    >>> y = np.array([0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    >>> groups = np.array([1, 1, 2, 2, 3, 3, 3, 4, 5, 5, 5, 5, 6, 6, 7, 8, 8])
    >>> cv = RepeatedStratifiedGroupKFold(n_splits=2, n_repeats=2,
    ...                                   random_state=36851234)
    >>> for train_idxs, test_idxs in cv.split(X, y, groups):
    ...     print("TRAIN:", groups[train_idxs])
    ...     print("      ", y[train_idxs])
    ...     print(" TEST:", groups[test_idxs])
    ...     print("      ", y[test_idxs])
    TRAIN: [2 2 4 5 5 5 5 8 8]
           [1 1 1 0 0 0 0 0 0]
     TEST: [1 1 3 3 3 6 6 7]
           [0 0 1 1 1 0 0 0]
    TRAIN: [1 1 3 3 3 6 6 7]
           [0 0 1 1 1 0 0 0]
     TEST: [2 2 4 5 5 5 5 8 8]
           [1 1 1 0 0 0 0 0 0]
    TRAIN: [3 3 3 4 7 8 8]
           [1 1 1 1 0 0 0]
     TEST: [1 1 2 2 5 5 5 5 6 6]
           [0 0 1 1 0 0 0 0 0 0]
    TRAIN: [1 1 2 2 5 5 5 5 6 6]
           [0 0 1 1 0 0 0 0 0 0]
     TEST: [3 3 3 4 7 8 8]
           [1 1 1 1 0 0 0]

    Notes
    -----
    Randomized CV splitters may return different results for each call of
    split. You can make the results identical by setting `random_state`
    to an integer.

    See also
    --------
    RepeatedStratifiedKFold: Repeats Stratified K-Fold n times.
    """

    def __init__(self, n_splits=5, n_repeats=10, random_state=None):
        super().__init__(StratifiedGroupKFold, n_splits=n_splits,
                         n_repeats=n_repeats, random_state=random_state)


class StratifiedSampleFromGroupShuffleSplit(StratifiedShuffleSplit):

    def split(self, X, y, groups, weights=None):
        X, y, groups, weights = indexable(X, y, groups, weights)
        y = check_array(y, ensure_2d=False, dtype=None)
        if groups is None:
            raise ValueError('The groups parameter should not be None')
        groups = check_array(groups, ensure_2d=False, dtype=None)
        if weights is not None:
            weights = check_array(weights, ensure_2d=False, dtype=None)

        rng = check_random_state(self.random_state)
        for _ in range(self.n_splits):
            if weights is None:
                indices = (pd.DataFrame({'group': groups}).groupby('group')
                           .apply(pd.DataFrame.sample, n=1, random_state=rng)
                           .reset_index(level=0, drop=True).sort_index()
                           .index.values)
            else:
                indices = (pd.DataFrame({'group': groups, 'weight': weights})
                           .groupby('group')
                           .apply(pd.DataFrame.sample, n=1, random_state=rng,
                                  weights='weight')
                           .reset_index(level=0, drop=True).sort_index()
                           .index.values)

            Xr = X.iloc[indices] if hasattr(X, 'iloc') else X[indices]
            yr = y[indices]

            n_samples = _num_samples(Xr)
            n_train, n_test = _validate_shuffle_split(
                n_samples, self.test_size, self.train_size,
                default_test_size=self._default_test_size)

            if yr.ndim == 2:
                # for multi-label y, map each distinct row to a string repr
                # using join because str(row) uses an ellipsis if
                # len(row) > 1000
                yr = np.array([' '.join(row.astype('str')) for row in yr])

            classes, y_indices = np.unique(yr, return_inverse=True)
            n_classes = classes.shape[0]

            class_counts = np.bincount(y_indices)
            if np.min(class_counts) < 2:
                raise ValueError("The least populated class in y has only 1"
                                 " member, which is too few. The minimum"
                                 " number of groups for any class cannot"
                                 " be less than 2.")

            if n_train < n_classes:
                raise ValueError('The train_size = %d should be greater or '
                                 'equal to the number of classes = %d' %
                                 (n_train, n_classes))
            if n_test < n_classes:
                raise ValueError('The test_size = %d should be greater or '
                                 'equal to the number of classes = %d' %
                                 (n_test, n_classes))

            # Find the sorted list of instances for each class:
            # (np.unique above performs a sort, so code is O(n logn) already)
            class_indices = np.split(np.argsort(y_indices, kind='mergesort'),
                                     np.cumsum(class_counts)[:-1])

            # if there are ties in the class-counts, we want
            # to make sure to break them anew in each iteration
            n_i = _approximate_mode(class_counts, n_train, rng)
            class_counts_remaining = class_counts - n_i
            t_i = _approximate_mode(class_counts_remaining, n_test, rng)

            train = []
            test = []

            for i in range(n_classes):
                permutation = rng.permutation(class_counts[i])
                perm_indices_class_i = class_indices[i].take(permutation,
                                                             mode='clip')

                train.extend(perm_indices_class_i[:n_i[i]])
                test.extend(perm_indices_class_i[n_i[i]:n_i[i] + t_i[i]])

            train = rng.permutation(train)
            test = rng.permutation(test)

            train_index = indices[train]
            test_index = indices[test]

            yield train_index, test_index

    def get_n_splits(self, X=None, y=None, groups=None, weights=None):
        return self.n_splits
