"""
The :mod:`sklearn_extensions.pipeline` module implements utilities to build a
composite estimator, as a chain of transforms and estimators.
"""
# Author: Edouard Duchesnay
#         Gael Varoquaux
#         Virgile Fritsch
#         Alexandre Gramfort
#         Lars Buitinck
#         Leandro Cruz Hermida <hermidal@cs.umd.edu>
# License: BSD 3 Clause

from collections import defaultdict
from itertools import islice

import numpy as np
import pandas as pd

from sklearn.base import BaseEstimator, clone
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.utils import Bunch, _print_elapsed_time
from sklearn.utils.metaestimators import if_delegate_has_method
from sklearn.utils.validation import check_memory
from .utils.metaestimators import check_routing


__all__ = ['ExtendedPipeline']


class ExtendedPipeline(Pipeline):
    """Pipeline of transforms with a final estimator.

    Sequentially apply a list of transforms and a final estimator.
    Intermediate steps of the pipeline must be 'transforms', that is, they
    must implement fit and transform methods.
    The final estimator only needs to implement fit.
    The transformers in the pipeline can be cached using ``memory`` argument.

    The purpose of the pipeline is to assemble several steps that can be
    cross-validated together while setting different parameters.
    For this, it enables setting parameters of the various steps using their
    names and the parameter name separated by a '__', as in the example below.
    A step's estimator may be replaced entirely by setting the parameter
    with its name to another estimator, or a transformer removed by setting
    it to 'passthrough' or ``None``.

    Read more in the :ref:`User Guide <pipeline>`.

    .. versionadded:: 0.5

    Parameters
    ----------
    steps : list
        List of (name, transform) tuples (implementing fit/transform) that are
        chained, in the order in which they are chained, with the last object
        an estimator.

    memory : None, str or object with the joblib.Memory interface, optional
        Used to cache the fitted transformers of the pipeline. By default,
        no caching is performed. If a string is given, it is the path to
        the caching directory. Enabling caching triggers a clone of
        the transformers before fitting. Therefore, the transformer
        instance given to the pipeline cannot be inspected
        directly. Use the attribute ``named_steps`` or ``steps`` to
        inspect estimators within the pipeline. Caching the
        transformers is advantageous when fitting is time consuming.

    verbose : bool, default=False
        If True, the time elapsed while fitting each step will be printed as it
        is completed.

    Attributes
    ----------
    named_steps : bunch object, a dictionary with attribute access
        Read-only attribute to access any step parameter by user given name.
        Keys are step names and values are steps parameters.

    See Also
    --------
    sklearn.pipeline.make_pipeline : Convenience function for simplified
        pipeline construction.

    Examples
    --------
    >>> from sklearn import svm
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.feature_selection import SelectKBest
    >>> from sklearn.feature_selection import f_regression
    >>> from sklearn.pipeline import Pipeline
    >>> # generate some data to play with
    >>> X, y = make_classification(
    ...     n_informative=5, n_redundant=0, random_state=42)
    >>> # ANOVA SVM-C
    >>> anova_filter = SelectKBest(f_regression, k=5)
    >>> clf = svm.SVC(kernel='linear')
    >>> anova_svm = Pipeline([('anova', anova_filter), ('svc', clf)])
    >>> # You can set the parameters using the names issued
    >>> # For instance, fit using a k of 10 in the SelectKBest
    >>> # and a parameter 'C' of the svm
    >>> anova_svm.set_params(anova__k=10, svc__C=.1).fit(X, y)
    Pipeline(steps=[('anova', SelectKBest(...)), ('svc', SVC(...))])
    >>> prediction = anova_svm.predict(X)
    >>> anova_svm.score(X, y)
    0.83
    >>> # getting the selected features chosen by anova_filter
    >>> anova_svm['anova'].get_support()
    array([False, False,  True,  True, False, False,  True,  True, False,
           True, False,  True,  True, False,  True, False,  True,  True,
           False, False])
    >>> # Another way to get selected features chosen by anova_filter
    >>> anova_svm.named_steps.anova.get_support()
    array([False, False,  True,  True, False, False,  True,  True, False,
           True, False,  True,  True, False,  True, False,  True,  True,
           False, False])
    >>> # Indexing can also be used to extract a sub-pipeline.
    >>> sub_pipeline = anova_svm[:1]
    >>> sub_pipeline
    Pipeline(steps=[('anova', SelectKBest(...))])
    >>> coef = anova_svm[-1].coef_
    >>> anova_svm['svc'] is anova_svm[-1]
    True
    >>> coef.shape
    (1, 10)
    >>> sub_pipeline.inverse_transform(coef).shape
    (1, 20)
    """

    # BaseEstimator interface
    _required_parameters = ['steps']

    def __init__(self, steps, memory=None, verbose=False, param_routing=None):
        self.steps = steps
        self.memory = memory
        self.verbose = verbose
        self.param_routing = param_routing
        self._validate_steps()
        self.router = check_routing(self.param_routing,
                                    [[name, '*'] for name, _ in self.steps],
                                    self._default_routing)

    def get_params(self, deep=True):
        """Get parameters for this estimator.

        Parameters
        ----------
        deep : boolean, optional
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : mapping of string to any
            Parameter names mapped to their values.
        """
        return self._get_params('steps', deep=deep)

    def set_params(self, **kwargs):
        """Set the parameters of this estimator.

        Valid parameter keys can be listed with ``get_params()``.

        Returns
        -------
        self
        """
        self._set_params('steps', **kwargs)
        if 'param_routing' in kwargs:
            self.router = check_routing(
                self.param_routing, [[name, '*'] for name, _ in self.steps],
                self._default_routing)
        return self

    def _validate_steps(self):
        names, estimators = zip(*self.steps)

        # validate names
        self._validate_names(names)

        # validate estimators
        transformers = estimators[:-1]
        estimator = estimators[-1]

        for t in transformers:
            if t is None or t == 'passthrough':
                continue
            if (not (hasattr(t, "fit") or hasattr(t, "fit_transform")) or not
                    hasattr(t, "transform")):
                raise TypeError("All intermediate steps should be "
                                "transformers and implement fit and transform "
                                "or be the string 'passthrough' "
                                "'%s' (type %s) doesn't" % (t, type(t)))

        # We allow last estimator to be None as an identity transformation
        if (estimator is not None and estimator != 'passthrough'
                and not hasattr(estimator, "fit")):
            raise TypeError(
                "Last step of Pipeline should implement fit "
                "or be the string 'passthrough'. "
                "'%s' (type %s) doesn't" % (estimator, type(estimator)))

    def _iter(self, with_final=True, filter_passthrough=True):
        """
        Generate (idx, (name, trans)) tuples from self.steps

        When filter_passthrough is True, 'passthrough' and None transformers
        are filtered out.
        """
        stop = len(self.steps)
        if not with_final:
            stop -= 1

        for idx, (name, trans) in enumerate(islice(self.steps, 0, stop)):
            if not filter_passthrough:
                yield idx, name, trans
            elif trans is not None and trans != 'passthrough':
                yield idx, name, trans

    def __len__(self):
        """
        Returns the length of the Pipeline
        """
        return len(self.steps)

    def __getitem__(self, ind):
        """Returns a sub-pipeline or a single esimtator in the pipeline

        Indexing with an integer will return an estimator; using a slice
        returns another Pipeline instance which copies a slice of this
        Pipeline. This copy is shallow: modifying (or fitting) estimators in
        the sub-pipeline will affect the larger pipeline and vice-versa.
        However, replacing a value in `step` will not affect a copy.
        """
        if isinstance(ind, slice):
            if ind.step not in (1, None):
                raise ValueError('Pipeline slicing only supports a step of 1')
            return self.__class__(self.steps[ind])
        try:
            name, est = self.steps[ind]
        except TypeError:
            # Not an int, try get step by name
            return self.named_steps[ind]
        return est

    @property
    def _estimator_type(self):
        return self.steps[-1][1]._estimator_type

    @property
    def named_steps(self):
        # Use Bunch object to improve autocomplete
        return Bunch(**dict(self.steps))

    @property
    def _final_estimator(self):
        estimator = self.steps[-1][1]
        return 'passthrough' if estimator is None else estimator

    def _log_message(self, step_idx):
        if not self.verbose:
            return None
        name, step = self.steps[step_idx]

        return '(step %d of %d) Processing %s' % (step_idx + 1,
                                                  len(self.steps),
                                                  name)

    def _default_routing(self, params):
        names = [name for name, _ in self.steps]
        step_params = {name: {} for name in names}
        remainder = set()
        for k, v in params.items():
            if v is not None:
                # XXX: not sure if we need to remove Nones
                name, prop = k.split('__', 1)
                try:
                    step_params[name][prop] = v
                except KeyError:
                    remainder.add(k)
        return [step_params[name] for name in names], remainder

    def _transform_feature_meta(self, estimator, feature_meta):
        transformed_feature_meta = None
        if isinstance(estimator, ColumnTransformer):
            for _, trf_transformer, trf_columns in estimator.transformers_:
                if (isinstance(trf_transformer, str)
                        and trf_transformer == 'drop'):
                    trf_feature_meta = feature_meta.iloc[
                        ~feature_meta.index.isin(trf_columns)]
                elif ((isinstance(trf_columns, slice)
                       and (isinstance(trf_columns.start, str)
                            or isinstance(trf_columns.stop, str)))
                      or isinstance(trf_columns[0], str)):
                    trf_feature_meta = feature_meta.loc[trf_columns]
                else:
                    trf_feature_meta = feature_meta.iloc[trf_columns]
                if isinstance(trf_transformer, BaseEstimator):
                    for transformer in trf_transformer:
                        if hasattr(transformer, 'get_support'):
                            trf_feature_meta = trf_feature_meta.loc[
                                transformer.get_support()]
                        elif hasattr(transformer, 'get_feature_names'):
                            new_trf_feature_names = (
                                transformer.get_feature_names(
                                    input_features=(trf_feature_meta.index
                                                    .values)).astype(str))
                            new_trf_feature_meta = None
                            for feature_name in trf_feature_meta.index:
                                f_feature_meta = pd.concat(
                                    [trf_feature_meta.loc[[feature_name]]]
                                    * np.sum(np.char.startswith(
                                        new_trf_feature_names,
                                        '{}_'.format(feature_name))),
                                    axis=0, ignore_index=True)
                                if new_trf_feature_meta is None:
                                    new_trf_feature_meta = f_feature_meta
                                else:
                                    new_trf_feature_meta = pd.concat(
                                        [new_trf_feature_meta, f_feature_meta],
                                        axis=0, ignore_index=True)
                            trf_feature_meta = new_trf_feature_meta.set_index(
                                new_trf_feature_names)
                if transformed_feature_meta is None:
                    transformed_feature_meta = trf_feature_meta
                else:
                    transformed_feature_meta = pd.concat(
                        [transformed_feature_meta, trf_feature_meta], axis=0)
        else:
            if transformed_feature_meta is None:
                transformed_feature_meta = feature_meta
            if hasattr(estimator, 'get_support'):
                transformed_feature_meta = (
                    transformed_feature_meta.loc[estimator.get_support()])
            elif hasattr(estimator, 'get_feature_names'):
                new_feature_names = estimator.get_feature_names(
                    input_features=transformed_feature_meta.index.values
                ).astype(str)
                new_transformed_feature_meta = None
                for feature_name in transformed_feature_meta.index:
                    f_feature_meta = pd.concat(
                        [transformed_feature_meta.loc[[feature_name]]]
                        * np.sum(np.char.startswith(
                            new_feature_names, '{}_'.format(feature_name))),
                        axis=0, ignore_index=True)
                    if new_transformed_feature_meta is None:
                        new_transformed_feature_meta = f_feature_meta
                    else:
                        new_transformed_feature_meta = pd.concat(
                            [new_transformed_feature_meta, f_feature_meta],
                            axis=0, ignore_index=True)
                transformed_feature_meta = (new_transformed_feature_meta
                                            .set_index(new_feature_names))
        return transformed_feature_meta

    def _transform_pipeline(self, caller_name, X, params):
        step_params, remainder = self.router(params)
        if remainder:
            raise TypeError('%s() got unexpected keyword arguments %r'
                            % (caller_name, sorted(remainder)))
        feature_meta = params.get('feature_meta', None)
        if caller_name == 'transform':
            step_iter = self._iter()
        elif caller_name == 'inverse_transform':
            step_iter = reversed(list(self._iter()))
        else:
            step_iter = self._iter(with_final=False)
        for step_idx, _, transformer in step_iter:
            if (step_idx > 0 and feature_meta is not None
                    and 'feature_meta' in step_params[step_idx]):
                step_params[step_idx]['feature_meta'] = feature_meta
            X = transformer.transform(X, **step_params[step_idx])
            if feature_meta is not None:
                feature_meta = self._transform_feature_meta(transformer,
                                                            feature_meta)
        if caller_name in ['transform', 'inverse_transform']:
            return X
        if 'feature_meta' in step_params[-1]:
            step_params[-1]['feature_meta'] = feature_meta
        return X, step_params[-1]

    # Estimator interface

    def _fit(self, X, y=None, **fit_params):
        # shallow copy of steps - this should really be steps_
        self.steps = list(self.steps)
        self._validate_steps()
        # Setup the memory
        memory = check_memory(self.memory)

        fit_transform_one_cached = memory.cache(_fit_transform_one)

        step_fit_params, remainder = self.router(fit_params)
        if remainder:
            raise TypeError('Got unexpected keyword arguments %r'
                            % sorted(remainder))
        feature_meta = fit_params.get('feature_meta', None)
        for (step_idx,
             name,
             transformer) in self._iter(with_final=False,
                                        filter_passthrough=False):
            if (transformer is None or transformer == 'passthrough'):
                with _print_elapsed_time('Pipeline',
                                         self._log_message(step_idx)):
                    continue

            if hasattr(memory, 'location'):
                # joblib >= 0.12
                if memory.location is None:
                    # we do not clone when caching is disabled to
                    # preserve backward compatibility
                    cloned_transformer = transformer
                else:
                    cloned_transformer = clone(transformer)
            elif hasattr(memory, 'cachedir'):
                # joblib < 0.11
                if memory.cachedir is None:
                    # we do not clone when caching is disabled to
                    # preserve backward compatibility
                    cloned_transformer = transformer
                else:
                    cloned_transformer = clone(transformer)
            else:
                cloned_transformer = clone(transformer)

            if step_idx > 0 and 'feature_meta' in step_fit_params[step_idx]:
                step_fit_params[step_idx]['feature_meta'] = feature_meta

            # Fit or load from cache the current transformer
            X, fitted_transformer = fit_transform_one_cached(
                cloned_transformer, X, y, None,
                message_clsname='Pipeline',
                message=self._log_message(step_idx),
                **step_fit_params[step_idx])

            if feature_meta is not None:
                feature_meta = self._transform_feature_meta(fitted_transformer,
                                                            feature_meta)
                if X.shape[1] != feature_meta.shape[0]:
                    raise ValueError(('X ({:d}) and feature_meta ({:d}) have '
                                      'different feature dimensions').format(
                                          X.shape[1], feature_meta.shape[0]))

            # Replace the transformer of the step with the fitted
            # transformer. This is necessary when loading the transformer
            # from the cache.
            self.steps[step_idx] = (name, fitted_transformer)
        if self._final_estimator == 'passthrough':
            return X, {}
        if 'feature_meta' in step_fit_params[-1]:
            step_fit_params[-1]['feature_meta'] = feature_meta
        return X, step_fit_params[-1]

    def fit(self, X, y=None, **fit_params):
        """Fit the model

        Fit all the transforms one after the other and transform the
        data, then fit the transformed data using the final estimator.

        Parameters
        ----------
        X : iterable
            Training data. Must fulfill input requirements of first step of the
            pipeline.

        y : iterable, default=None
            Training targets. Must fulfill label requirements for all steps of
            the pipeline.

        **fit_params : dict of string -> object
            Parameters passed to the ``fit`` method of each step, where
            each parameter name is prefixed such that parameter ``p`` for step
            ``s`` has key ``s__p``.

        Returns
        -------
        self : Pipeline
            This estimator
        """
        Xt, fit_params = self._fit(X, y, **fit_params)
        with _print_elapsed_time('Pipeline',
                                 self._log_message(len(self.steps) - 1)):
            if self._final_estimator != 'passthrough':
                self._final_estimator.fit(Xt, y, **fit_params)
        return self

    def fit_transform(self, X, y=None, **fit_params):
        """Fit the model and transform with the final estimator

        Fits all the transforms one after the other and transforms the
        data, then uses fit_transform on transformed data with the final
        estimator.

        Parameters
        ----------
        X : iterable
            Training data. Must fulfill input requirements of first step of the
            pipeline.

        y : iterable, default=None
            Training targets. Must fulfill label requirements for all steps of
            the pipeline.

        **fit_params : dict of string -> object
            Parameters passed to the ``fit`` method of each step, where
            each parameter name is prefixed such that parameter ``p`` for step
            ``s`` has key ``s__p``.

        Returns
        -------
        Xt : array-like of shape  (n_samples, n_transformed_features)
            Transformed samples
        """
        last_step = self._final_estimator
        Xt, fit_params = self._fit(X, y, **fit_params)
        with _print_elapsed_time('Pipeline',
                                 self._log_message(len(self.steps) - 1)):
            if last_step == 'passthrough':
                return Xt
            if hasattr(last_step, 'fit_transform'):
                return last_step.fit_transform(Xt, y, **fit_params)
            else:
                return last_step.fit(Xt, y, **fit_params).transform(
                    Xt, **fit_params)

    @if_delegate_has_method(delegate='_final_estimator')
    def predict(self, X, **predict_params):
        """Apply transforms to the data, and predict with the final estimator

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step
            of the pipeline.

        **predict_params : dict of string -> object
            Parameters to the ``predict`` called at the end of all
            transformations in the pipeline. Note that while this may be
            used to return uncertainties from some models with return_std
            or return_cov, uncertainties that are generated by the
            transformations in the pipeline are not propagated to the
            final estimator.

        Returns
        -------
        y_pred : array-like
        """
        Xt, predict_params = self._transform_pipeline('predict',
                                                      X, predict_params)
        return self.steps[-1][-1].predict(Xt, **predict_params)

    @if_delegate_has_method(delegate='_final_estimator')
    def fit_predict(self, X, y=None, **fit_params):
        """Applies fit_predict of last step in pipeline after transforms.

        Applies fit_transforms of a pipeline to the data, followed by the
        fit_predict method of the final estimator in the pipeline. Valid
        only if the final estimator implements fit_predict.

        Parameters
        ----------
        X : iterable
            Training data. Must fulfill input requirements of first step of
            the pipeline.

        y : iterable, default=None
            Training targets. Must fulfill label requirements for all steps
            of the pipeline.

        **fit_params : dict of string -> object
            Parameters passed to the ``fit`` method of each step, where
            each parameter name is prefixed such that parameter ``p`` for step
            ``s`` has key ``s__p``.

        Returns
        -------
        y_pred : array-like
        """
        Xt, fit_params = self._fit(X, y, **fit_params)
        with _print_elapsed_time('Pipeline',
                                 self._log_message(len(self.steps) - 1)):
            y_pred = self.steps[-1][-1].fit_predict(Xt, y, **fit_params)
        return y_pred

    @if_delegate_has_method(delegate='_final_estimator')
    def predict_proba(self, X, **predict_params):
        """Apply transforms, and predict_proba of the final estimator

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step
            of the pipeline.

        Returns
        -------
        y_proba : array-like of shape (n_samples, n_classes)
        """
        Xt, predict_params = self._transform_pipeline('predict_proba',
                                                      X, predict_params)
        return self.steps[-1][-1].predict_proba(Xt)

    @if_delegate_has_method(delegate='_final_estimator')
    def decision_function(self, X, **predict_params):
        """Apply transforms, and decision_function of the final estimator

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step
            of the pipeline.

        Returns
        -------
        y_score : array-like of shape (n_samples, n_classes)
        """
        Xt, predict_params = self._transform_pipeline('decision_function',
                                                      X, predict_params)
        return self.steps[-1][-1].decision_function(Xt)

    @if_delegate_has_method(delegate='_final_estimator')
    def score_samples(self, X, **predict_params):
        """Apply transforms, and score_samples of the final estimator.

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step
            of the pipeline.

        Returns
        -------
        y_score : ndarray, shape (n_samples,)
        """
        Xt, predict_params = self._transform_pipeline('score_samples',
                                                      X, predict_params)
        return self.steps[-1][-1].score_samples(Xt)

    @if_delegate_has_method(delegate='_final_estimator')
    def predict_log_proba(self, X, **predict_params):
        """Apply transforms, and predict_log_proba of the final estimator

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step
            of the pipeline.

        Returns
        -------
        y_score : array-like of shape (n_samples, n_classes)
        """
        Xt, predict_params = self._transform_pipeline('predict_log_proba',
                                                      X, predict_params)
        return self.steps[-1][-1].predict_log_proba(Xt)

    @property
    def transform(self):
        """Apply transforms, and transform with the final estimator

        This also works where final estimator is ``None``: all prior
        transformations are applied.

        Parameters
        ----------
        X : iterable
            Data to transform. Must fulfill input requirements of first step
            of the pipeline.

        Returns
        -------
        Xt : array-like of shape  (n_samples, n_transformed_features)
        """
        # _final_estimator is None or has transform, otherwise attribute error
        # XXX: Handling the None case means we can't use if_delegate_has_method
        if self._final_estimator != 'passthrough':
            self._final_estimator.transform
        return self._transform

    def _transform(self, X, **transform_params):
        Xt = self._transform_pipeline('transform', X, transform_params)
        return Xt

    @property
    def inverse_transform(self):
        """Apply inverse transformations in reverse order

        All estimators in the pipeline must support ``inverse_transform``.

        Parameters
        ----------
        Xt : array-like of shape  (n_samples, n_transformed_features)
            Data samples, where ``n_samples`` is the number of samples and
            ``n_features`` is the number of features. Must fulfill
            input requirements of last step of pipeline's
            ``inverse_transform`` method.

        Returns
        -------
        Xt : array-like of shape (n_samples, n_features)
        """
        # raise AttributeError if necessary for hasattr behaviour
        # XXX: Handling the None case means we can't use if_delegate_has_method
        for _, _, transform in self._iter():
            transform.inverse_transform
        return self._inverse_transform

    def _inverse_transform(self, X, **transform_params):
        Xt = self._transform_pipeline('inverse_transform', X, transform_params)
        return Xt

    @if_delegate_has_method(delegate='_final_estimator')
    def score(self, X, y=None, sample_weight=None, **transform_params):
        """Apply transforms, and score with the final estimator

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step
            of the pipeline.

        y : iterable, default=None
            Targets used for scoring. Must fulfill label requirements for all
            steps of the pipeline.

        sample_weight : array-like, default=None
            If not None, this argument is passed as ``sample_weight`` keyword
            argument to the ``score`` method of the final estimator.

        Returns
        -------
        score : float
        """
        Xt, _ = self._transform_pipeline('score', X, transform_params)
        score_params = {}
        if sample_weight is not None:
            score_params['sample_weight'] = sample_weight
        return self.steps[-1][-1].score(Xt, y, **score_params)

    @property
    def classes_(self):
        return self.steps[-1][-1].classes_

    @property
    def _pairwise(self):
        # check if first estimator expects pairwise input
        return getattr(self.steps[0][1], '_pairwise', False)


def _name_estimators(estimators):
    """Generate names for estimators."""

    names = [
        estimator
        if isinstance(estimator, str) else type(estimator).__name__.lower()
        for estimator in estimators
    ]
    namecount = defaultdict(int)
    for est, name in zip(estimators, names):
        namecount[name] += 1

    for k, v in list(namecount.items()):
        if v == 1:
            del namecount[k]

    for i in reversed(range(len(estimators))):
        name = names[i]
        if name in namecount:
            names[i] += "-%d" % namecount[name]
            namecount[name] -= 1

    return list(zip(names, estimators))


def _transform_one(transformer, X, y, weight, message_clsname='',
                   message=None, **transform_params):
    res = transformer.transform(X, **transform_params)
    # if we have a weight for this transformer, multiply output
    if weight is None:
        return res
    return res * weight


def _fit_transform_one(transformer,
                       X,
                       y,
                       weight,
                       message_clsname='',
                       message=None,
                       **fit_params):
    """
    Fits ``transformer`` to ``X`` and ``y``. The transformed result is returned
    with the fitted transformer. If ``weight`` is not ``None``, the result will
    be multiplied by ``weight``.
    """
    with _print_elapsed_time(message_clsname, message):
        if hasattr(transformer, 'fit_transform'):
            res = transformer.fit_transform(X, y, **fit_params)
        else:
            res = transformer.fit(X, y, **fit_params).transform(
                X, **fit_params)

    if weight is None:
        return res, transformer
    return res * weight, transformer


def _fit_one(transformer,
             X,
             y,
             weight,
             message_clsname='',
             message=None,
             **fit_params):
    """
    Fits ``transformer`` to ``X`` and ``y``.
    """
    with _print_elapsed_time(message_clsname, message):
        return transformer.fit(X, y, **fit_params)
