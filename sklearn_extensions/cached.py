# Author: Leandro Cruz Hermida <hermidal@cs.umd.edu>

from sklearn.utils.validation import check_memory


class CachedFitMixin:
    """Mixin for caching pipeline nested estimator fits"""

    def fit(self, *args, **kwargs):
        memory = check_memory(self.memory)
        cached_self = memory.cache(super().fit)(*args, **kwargs)
        vars(self).update(vars(cached_self))
        return self
