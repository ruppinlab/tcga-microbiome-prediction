# Authors: Gilles Louppe, Mathieu Blondel, Maheshakya Wijewardena,
#          Leandro Hermida
# License: BSD 3 clause

import numbers
import numpy as np

from sklearn.base import BaseEstimator, clone, MetaEstimatorMixin
from sklearn.exceptions import NotFittedError
from sklearn.utils import check_X_y
from sklearn.utils.metaestimators import if_delegate_has_method

from ._base import ExtendedSelectorMixin


def _get_feature_importances(estimator, norm_order=1):
    """Retrieve or aggregate feature importances from estimator"""
    importances = getattr(estimator, 'feature_importances_', None)
    if importances is None:
        if hasattr(estimator, 'coef_'):
            if estimator.coef_.ndim == 1:
                importances = np.abs(estimator.coef_)
            else:
                importances = np.linalg.norm(estimator.coef_, axis=0,
                                             ord=norm_order)
        else:
            raise ValueError(
                'The underlying estimator %s has no `coef_` or '
                '`feature_importances_` attribute. Either pass a fitted '
                'estimator to SelectFromModel or call fit before calling '
                'transform.' % estimator.__class__.__name__)
    return importances


def _calculate_threshold(estimator, importances, threshold):
    """Interpret the threshold value"""
    if threshold is None:
        # determine default from estimator
        est_name = estimator.__class__.__name__
        if ((hasattr(estimator, 'penalty') and estimator.penalty == 'l1')
                or 'Lasso' in est_name):
            # the natural default threshold is 0 when l1 penalty was used
            threshold = 1e-5
        else:
            threshold = 'mean'
    if isinstance(threshold, str):
        if '*' in threshold:
            scale, reference = threshold.split('*')
            scale = float(scale.strip())
            reference = reference.strip()
            if reference == 'median':
                reference = np.median(importances)
            elif reference == 'mean':
                reference = np.mean(importances)
            else:
                raise ValueError('Unknown reference: ' + reference)
            threshold = scale * reference
        elif threshold == 'median':
            threshold = np.median(importances)
        elif threshold == 'mean':
            threshold = np.mean(importances)
        else:
            raise ValueError("Expected threshold='mean' or threshold='median' "
                             "got %s" % threshold)
    else:
        threshold = float(threshold)
    return threshold


class SelectFromModel(MetaEstimatorMixin, ExtendedSelectorMixin,
                      BaseEstimator):
    """Meta-transformer for selecting features based on importance weights.

    .. versionadded:: 0.17

    Parameters
    ----------
    estimator : object
        The base estimator from which the transformer is built.
        This can be both a fitted (if ``prefit`` is set to True)
        or a non-fitted estimator. The estimator must have either a
        ``feature_importances_`` or ``coef_`` attribute after fitting.

    threshold : string, float, optional default None
        The threshold value to use for feature selection. Features whose
        importance is greater or equal are kept while the others are
        discarded. If "median" (resp. "mean"), then the ``threshold`` value is
        the median (resp. the mean) of the feature importances. A scaling
        factor (e.g., "1.25*mean") may also be used. If None and if the
        estimator has a parameter penalty set to l1, either explicitly
        or implicitly (e.g, Lasso), the threshold used is 1e-5.
        Otherwise, "mean" is used by default.

    prefit : bool, default False
        Whether a prefit model is expected to be passed into the constructor
        directly or not. If True, ``transform`` must be called directly
        and SelectFromModel cannot be used with ``cross_val_score``,
        ``GridSearchCV`` and similar utilities that clone the estimator.
        Otherwise train the model using ``fit`` and then ``transform`` to do
        feature selection.

    norm_order : non-zero int, inf, -inf, default 1
        Order of the norm used to filter the vectors of coefficients below
        ``threshold`` in the case where the ``coef_`` attribute of the
        estimator is of dimension 2.

    max_features : int or None, optional
        The maximum number of features selected scoring above ``threshold``.
        To disable ``threshold`` and only select based on ``max_features``,
        set ``threshold=-np.inf``.

        .. versionadded:: 0.20

    Attributes
    ----------
    estimator_ : an estimator
        The base estimator from which the transformer is built.
        This is stored only when a non-fitted estimator is passed to the
        ``SelectFromModel``, i.e when prefit is False.

    scores_ : array, shape = (n_features,)
        Feature importance for each feature.

    threshold_ : float
        The threshold value used for feature selection.

    Notes
    -----
    Allows NaN/Inf in the input if the underlying estimator does as well.

    Examples
    --------
    >>> from sklearn.feature_selection import SelectFromModel
    >>> from sklearn.linear_model import LogisticRegression
    >>> X = [[ 0.87, -1.34,  0.31 ],
    ...      [-2.79, -0.02, -0.85 ],
    ...      [-1.34, -0.48, -2.55 ],
    ...      [ 1.92,  1.48,  0.65 ]]
    >>> y = [0, 1, 0, 1]
    >>> selector = SelectFromModel(estimator=LogisticRegression()).fit(X, y)
    >>> selector.estimator_.coef_
    array([[-0.3252302 ,  0.83462377,  0.49750423]])
    >>> selector.threshold_
    0.55245...
    >>> selector.get_support()
    array([False,  True, False])
    >>> selector.transform(X)
    array([[-1.34],
           [-0.02],
           [-0.48],
           [ 1.48]])
    """

    def __init__(self, estimator, threshold=None, prefit=False, norm_order=1,
                 max_features=None):
        self.estimator = estimator
        self.threshold = threshold
        self.prefit = prefit
        self.norm_order = norm_order
        self.max_features = max_features

    @property
    def scores_(self):
        return _get_feature_importances(self._fitted_estimator(),
                                        self.norm_order)

    @property
    def threshold_(self):
        return _calculate_threshold(self.estimator, self.scores_,
                                    self.threshold)

    def fit(self, X, y=None, **fit_params):
        """Fit the SelectFromModel meta-transformer.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The training input samples.

        y : array-like, shape (n_samples,)
            The target values (integers that correspond to classes in
            classification, real numbers in regression).

        **fit_params : dict of string -> object
            Parameters passed to the ``fit`` method of the estimator

        Returns
        -------
        self : object
        """
        X, y = check_X_y(X, y, dtype=None)
        self._check_params(X, y)
        if self.prefit:
            raise NotFittedError("Since prefit=True, call transform directly")
        self.estimator_ = clone(self.estimator)
        self.estimator_.fit(X, y, **fit_params)
        return self

    @if_delegate_has_method('estimator')
    def partial_fit(self, X, y=None, **fit_params):
        """Fit the SelectFromModel meta-transformer only once.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The training input samples.

        y : array-like, shape (n_samples,)
            The target values (integers that correspond to classes in
            classification, real numbers in regression).

        **fit_params : Other estimator specific parameters

        Returns
        -------
        self : object
        """
        if self.prefit:
            raise NotFittedError('Since prefit=True, call transform directly')
        if not hasattr(self, 'estimator_'):
            self.estimator_ = clone(self.estimator)
        self.estimator_.partial_fit(X, y, **fit_params)
        return self

    def _more_tags(self):
        estimator_tags = self.estimator._get_tags()
        return {'allow_nan': estimator_tags.get('allow_nan', True)}

    def _check_params(self, X, y):
        if self.max_features is not None:
            if not isinstance(self.max_features, numbers.Integral):
                raise TypeError("'max_features' should be an integer between"
                                " 0 and {} features. Got {!r} instead."
                                .format(X.shape[1], self.max_features))
            if self.max_features < 0 or self.max_features > X.shape[1]:
                raise ValueError("'max_features' should be 0 and {} features."
                                 "Got {} instead."
                                 .format(X.shape[1], self.max_features))

    def _get_support_mask(self):
        estimator = self._fitted_estimator()
        scores = _get_feature_importances(estimator, self.norm_order)
        threshold = _calculate_threshold(estimator, scores, self.threshold)
        if self.max_features is not None:
            mask = np.zeros_like(scores, dtype=bool)
            mask[np.argsort(scores, kind='mergesort')
                 [-self.max_features:]] = True
        else:
            mask = np.ones_like(scores, dtype=bool)
        mask[scores < threshold] = False
        return mask

    def _fitted_estimator(self):
        if self.prefit:
            estimator = self.estimator
        elif hasattr(self, 'estimator_'):
            estimator = self.estimator_
        else:
            raise ValueError('Either fit the model before transform or set '
                             '"prefit=True" while passing the fitted '
                             'estimator to the constructor.')
        return estimator
