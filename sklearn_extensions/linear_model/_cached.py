# Author: Leandro Cruz Hermida <hermidal@cs.umd.edu>
# License: BSD 3 clause

from sklearn.linear_model import LogisticRegression
from ..cached import CachedFitMixin


class CachedLogisticRegression(CachedFitMixin, LogisticRegression):

    def __init__(self, memory, penalty='l2', dual=False, tol=1e-4, C=1.0,
                 fit_intercept=True, intercept_scaling=1, class_weight=None,
                 random_state=None, solver='lbfgs', max_iter=100,
                 multi_class='auto', verbose=0, warm_start=False, n_jobs=None,
                 l1_ratio=None):
        super().__init__(
            penalty=penalty, dual=dual, tol=tol, C=C,
            fit_intercept=fit_intercept, intercept_scaling=intercept_scaling,
            class_weight=class_weight, random_state=random_state,
            solver=solver, max_iter=max_iter, multi_class=multi_class,
            verbose=verbose, warm_start=warm_start, n_jobs=n_jobs,
            l1_ratio=l1_ratio)
        self.memory = memory
