# Author: Leandro Cruz Hermida <hermidal@cs.umd.edu>
# License: BSD 3 Clause

from sklearn_extensions.cached import CachedFitMixin

from ._coxnet import ExtendedCoxnetSurvivalAnalysis


class CachedExtendedCoxnetSurvivalAnalysis(CachedFitMixin,
                                           ExtendedCoxnetSurvivalAnalysis):

    def __init__(self, memory, n_alphas=100, alphas=None,
                 alpha_min_ratio=0.0001, l1_ratio=0.5, penalty_factor=None,
                 penalty_factor_meta_col=None, normalize=False, copy_X=True,
                 tol=1e-7, max_iter=100000, verbose=False,
                 fit_baseline_model=False):
        super().__init__(
            n_alphas=n_alphas, alphas=alphas, alpha_min_ratio=alpha_min_ratio,
            l1_ratio=l1_ratio, penalty_factor=penalty_factor,
            penalty_factor_meta_col=penalty_factor_meta_col,
            normalize=normalize, copy_X=copy_X, tol=tol, max_iter=max_iter,
            verbose=verbose, fit_baseline_model=fit_baseline_model)
        self.memory = memory
