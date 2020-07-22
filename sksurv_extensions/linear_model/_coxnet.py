# Author: Leandro Cruz Hermida <hermidal@cs.umd.edu>

"""Coxnet extensions."""

from sklearn.base import BaseEstimator, clone, MetaEstimatorMixin
from sklearn.utils import check_X_y
from sklearn.utils.metaestimators import if_delegate_has_method
from sklearn.utils.validation import check_is_fitted
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored


class ExtendedCoxnetSurvivalAnalysis(CoxnetSurvivalAnalysis):
    """Cox proportional hazards model with elastic net penalty and support for
    setting feature penalty_factor in a Pipeline context.

    See [1]_ for further description.

    Parameters
    ----------
    n_alphas : int (default=100)
        Number of alphas along the regularization path.

    alphas : array-like or None (default=None)
        List of alphas where to compute the models. If ``None`` alphas are set
        automatically.

    alpha_min_ratio : float (default=0.0001)
        Determines minimum alpha of the regularization path if ``alphas`` is
        ``None``. The smallest value for alpha is computed as the fraction of
        the data derived maximum alpha (i.e. the smallest value for which all
        coefficients are zero).

    l1_ratio : float (default=0.5)
        The ElasticNet mixing parameter, with ``0 < l1_ratio <= 1``.
        For ``l1_ratio = 0`` the penalty is an L2 penalty.
        For ``l1_ratio = 1`` it is an L1 penalty.
        For ``0 < l1_ratio < 1``, the penalty is a combination of L1 and L2.

    penalty_factor : array-like or None (default=None)
        Separate penalty factors can be applied to each coefficient. This is a
        number that multiplies alpha to allow differential shrinkage. Can be 0
        for some variables, which implies no shrinkage, and that variable is
        always included in the model. Default is 1 for all variables. Note: the
        penalty factors are internally rescaled to sum to n_features, and the
        alphas sequence will reflect this change.

    penalty_factor_meta_col : str (default=None)
        Feature metadata column name to use for ``penalty_factor``. This is
        ignored if ``penalty_factor`` is not None.

    normalize : boolean (default=False)
        If True, the features X will be normalized before optimization by
        subtracting the mean and dividing by the l2-norm. If you wish to
        standardize, please use :class:`sklearn.preprocessing.StandardScaler`
        before calling ``fit`` on an estimator with ``normalize=False``.

    copy_X : boolean (default=True)
        If ``True``, X will be copied; else, it may be overwritten.

    tol : float (default=1e-7)
        The tolerance for the optimization: optimization continues until all
        updates are smaller than ``tol``.

    max_iter : int (default=100000)
        The maximum number of iterations.

    verbose : bool (default=False)
        Whether to print additional information during optimization.

    fit_baseline_model : bool (default=False)
        Whether to estimate baseline survival function and baseline cumulative
        hazard function for each alpha. If enabled,
        :meth:`predict_cumulative_hazard_function` and
        :meth:`predict_survival_function` can be used to obtain predicted
        cumulative hazard function and survival function.

    Attributes
    ----------
    alphas_ : ndarray, shape=(n_alphas,)
        The actual sequence of alpha values used.

    penalty_factor_ : ndarray, shape=(n_features,)
        The actual penalty factors used.

    coef_ : ndarray, shape=(n_features, n_alphas)
        Matrix of coefficients.

    deviance_ratio_ : ndarray, shape=(n_alphas,)
        The fraction of (null) deviance explained.

    References
    ----------
    .. [1] Simon N, Friedman J, Hastie T, Tibshirani R.
           Regularization paths for Cox’s proportional hazards model via \
           coordinate descent.
           Journal of statistical software. 2011 Mar;39(5):1.
    """

    def __init__(self, n_alphas=100, alphas=None, alpha_min_ratio=0.0001,
                 l1_ratio=0.5, penalty_factor=None,
                 penalty_factor_meta_col=None, normalize=False, copy_X=True,
                 tol=1e-7, max_iter=100000, verbose=False,
                 fit_baseline_model=False):
        super().__init__(
            n_alphas=n_alphas, alphas=alphas, alpha_min_ratio=alpha_min_ratio,
            l1_ratio=l1_ratio, penalty_factor=penalty_factor,
            normalize=normalize, copy_X=copy_X, tol=tol, max_iter=max_iter,
            verbose=verbose, fit_baseline_model=fit_baseline_model)
        self.penalty_factor_meta_col = penalty_factor_meta_col

    def fit(self, X, y, feature_meta=None):
        """Fits the estimator.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            The training input data matrix.

        y : structured array, shape = (n_samples,)
            A structured array containing the binary event indicator as the
            first field, and time of event or time of censoring as the second
            field.

        feature_meta : pandas.DataFrame, pandas.Series (default = None), \
            shape = (n_features, n_metadata)
            Feature metadata.

        Returns
        -------
        self : object
        """
        X, y = check_X_y(X, y)
        self.__check_params(X, y, feature_meta)
        if (self.penalty_factor is None
                and self.penalty_factor_meta_col is not None):
            self.penalty_factor = (feature_meta[self.penalty_factor_meta_col]
                                   .to_numpy(dtype=float))
        return super().fit(X, y)

    def predict(self, X, alpha=None, feature_meta=None):
        """The linear predictor of the model.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Test data of which to calculate log-likelihood from

        alpha : float (default=None)
            Constant that multiplies the penalty terms. If the same alpha was
            used during training, exact coefficients are used, otherwise
            coefficients are interpolated from the closest alpha values that
            were used during training. If set to ``None``, the last alpha in
            the solution path is used.

        feature_meta : Ignored.

        Returns
        -------
        T : array, shape = (n_samples,)
            The predicted decision function
        """
        return super().predict(X, alpha=alpha)

    def score(self, X, y, alpha=None):
        """Returns the concordance index of the prediction.

        Parameters
        ----------
        X : array-like, shape = (n_samples, n_features)
            Test samples.

        y : structured array, shape = (n_samples,)
            A structured array containing the binary event indicator
            as first field, and time of event or time of censoring as
            second field.

        Returns
        -------
        cindex : float
            Estimated concordance index.
        """
        name_event, name_time = y.dtype.names
        result = concordance_index_censored(y[name_event], y[name_time],
                                            self.predict(X, alpha=alpha))
        return result[0]

    # double underscore to not override parent
    def __check_params(self, X, y, feature_meta):
        if self.penalty_factor_meta_col is not None:
            if feature_meta is None:
                raise ValueError('penalty_factor_meta_col specified but '
                                 'feature_meta not passed.')
            if self.penalty_factor_meta_col not in feature_meta.columns:
                raise ValueError('%s feature_meta column does not exist.'
                                 % self.penalty_factor_meta_col)
            if X.shape[1] != feature_meta.shape[0]:
                raise ValueError('X ({:d}) and feature_meta ({:d}) have '
                                 'different feature dimensions'
                                 .format(X.shape[1], feature_meta.shape[0]))


class MetaCoxnetSurvivalAnalysis(MetaEstimatorMixin, BaseEstimator):
    """Wrapper and redesign to CoxnetSurvivalAnalysis which allows for
    high-performance fitting of full alpha path and proper alpha hyperparmeter
    search within a scikit-learn SearchCV context.
    """

    def __init__(self, estimator, alpha=None):
        self.estimator = estimator
        self.alpha = alpha

    @property
    def coef_(self):
        check_is_fitted(self)
        return self.estimator_._get_coef(self.alpha)

    def fit(self, X, y, **fit_params):
        X, y = check_X_y(X, y, dtype=None)
        self._check_params(X, y)
        self.estimator_ = clone(self.estimator)
        self.estimator_.fit(X, y, **fit_params)
        return self

    @if_delegate_has_method(delegate='estimator')
    def predict(self, X, **predict_params):
        check_is_fitted(self)
        return self.estimator_.predict(X, alpha=self.alpha)

    @if_delegate_has_method(delegate='estimator')
    def predict_cumulative_hazard_function(self, X, **predict_params):
        check_is_fitted(self)
        return self.estimator_.predict_cumulative_hazard_function(
            X, alpha=self.alpha)

    @if_delegate_has_method(delegate='estimator')
    def predict_survival_function(self, X, **predict_params):
        check_is_fitted(self)
        return self.estimator_.predict_survival_function(X, alpha=self.alpha)

    @if_delegate_has_method(delegate='estimator')
    def score(self, X, y, sample_weight=None):
        check_is_fitted(self)
        score_params = {}
        if sample_weight is not None:
            score_params['sample_weight'] = sample_weight
        return self.estimator_.score(X, y, alpha=self.alpha, **score_params)

    def _more_tags(self):
        estimator_tags = self.estimator._get_tags()
        return {'allow_nan': estimator_tags.get('allow_nan', True)}

    def _check_params(self, X, y):
        if not isinstance(self.estimator, CoxnetSurvivalAnalysis):
            raise TypeError('{} estimator should be an instance of '
                            'CoxnetSurvivalAnalysis.'
                            .format(self.estimator.__class__.__name__))
        if self.estimator.alphas is not None:
            raise TypeError('{} estimator alphas parameter should be set to '
                            'None, got {} instead.'
                            .format(self.estimator.__class__.__name__,
                                    self.estimator.alphas))
