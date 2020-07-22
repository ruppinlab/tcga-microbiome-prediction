"""Utilities for meta-estimators"""
# Author: Joel Nothman
# License: BSD 3 clause

from collections import defaultdict
import copy
import six


class _NonePolicy:
    def apply(self, props, unused):
        return {}

    def __repr__(self):
        return 'None'


class _AllPolicy:
    def apply(self, props, unused):
        unused.clear()
        return props

    def update(self, other):
        pass

    def __repr__(self):
        return repr('*')


class _ExcludePolicy:
    def __init__(self, exclusions):
        # exclusions is a set of strings
        self.exclusions = exclusions

    def apply(self, props, unused):
        out = {k: v for k, v in props.items()
               if k not in self.exclusions}
        unused.intersection_update(self.exclusions)
        return out

    def update(self, other):
        self.exclusions.update(other)

    def __repr__(self):
        return repr(['-' + k for k in sorted(self.exclusions)])


class _IncludePolicy:
    def __init__(self, inclusions):
        # inclusions maps {tgt: src}
        self.inclusions = inclusions

    def apply(self, props, unused):
        unused.difference_update(self.inclusions.values())
        return {tgt: props[src]
                for tgt, src in self.inclusions.items()
                if src in props}

    def update(self, other):
        intersection = set(self.inclusions) & set(other.inclusions)
        if intersection:
            raise ValueError('Target property names {!r} are used multiple '
                             'times in the routing policy for the same '
                             'destination'.format(sorted(intersection)))
        self.inclusions.update(other.inclusions)

    def __repr__(self):
        return '{%s}' % ', '.join('{!r}: {!r}'.format(tgt, src)
                                  for tgt, src
                                  in self.inclusions.items())


class _Router:
    """Matches sample props to destinations according to a routing policy
    Parameters
    ----------
    routing : dict {dest: props} for dest str, props {str, list, dict}
        User-defined routing policy.
        Maps each destination string to the properties that should be provided
        to that destination. Props may be:
        - '*': provide all properties
        - a list of property names to include
        - a list of property names, each prefixed by '-', to exclude; all
          others will be provided
        - a single name to include or exclude
        - a dict mapping each target property name to its source property name
    dests : list of {str, iterable of str}
        The ordered destinations for this router. If a set of strings
        is provided for each entry, any of these strings will route
        parameters to the destination.
        Usually this is fixed for a metaestimator.
    Notes
    -----
    Abstracting away destination names/aliases in this way allows for providing
    syntactic sugar to users (e.g. Pipeline can declare '*' or 'steps' as an
    alias for "provide these props to ``fit`` of all steps).  While this may be
    instead facilitated by some string-based pattern matching, the present
    approach is more explicit and ensures backwards compatibility can be
    maintained.
    """

    def __init__(self, routing, dests):
        # can immediately:
        #     * check that all routes have valid dests
        #     * consolidate routing across aliases, such that each must be
        #       either a blacklist of length 0 or more or a mapping
        alias_to_idx = defaultdict(list)
        for i, aliases in enumerate(dests):
            if isinstance(aliases, six.string_types):
                aliases = [aliases]
            for alias in aliases:
                alias_to_idx[alias].append(i)
        alias_to_idx.default_factory = None

        policies = [None] * len(dests)

        # Sorted so error messages are deterministic
        for dest, props in sorted(routing.items()):
            if props == '*':
                policy = _AllPolicy()
            else:
                if isinstance(props, six.string_types):
                    props = [props]

                if isinstance(props, dict):
                    policy = _IncludePolicy(props)
                else:
                    minuses = [prop[:1] == '-' for prop in props]
                    if all(minuses):
                        policy = _ExcludePolicy({prop[1:] for prop in props})
                    elif any(minuses):
                        raise ValueError('Routing props should either all '
                                         'start with "-" or none should start '
                                         'with "-". Got a mix for %r' % dest)
                    else:
                        policy = _IncludePolicy({prop: prop for prop in props})

            # raises KeyError if unknown dest
            for idx in alias_to_idx[dest]:
                if policies[idx] is None:
                    policies[idx] = copy.deepcopy(policy)
                else:
                    if type(policies[idx]) is not type(policy):
                        raise ValueError('When handling routing for '
                                         'destination {!r}, found a mix of '
                                         'inclusion, exclusion and pass all '
                                         'policies.'.format(dest))
                    policies[idx].update(policy)

        self.policies = [_NonePolicy() if policy is None else policy
                         for policy in policies]

    def __call__(self, props):
        """Apply the routing policy to the given sample props
        Parameters
        ----------
        props : dict
        Returns
        -------
        dest_props : list of dicts
            Props to be passed to each destination declared in the constructor.
        unused : set of str
            Names of props that were not routed anywhere.
        """
        unused = set(props)
        out = [policy.apply(props, unused)
               for policy in self.policies]
        return out, unused


def check_routing(routing, dests, default=None):
    """Validates a prop_routing parameter and returns a router
    A router is a function which takes a dict of sample properties and returns
    a tuple ``(dest_props, unused)``, defined by:
        dest_props : list of dicts
            Props to be passed to each destination declared by ``dests``.
        unused : set of str
            Names of props that were not routed anywhere.
    Parameters
    ----------
    routing : dict of {dest: props}, callable or None if default is given
        User-defined routing policy.  A callable is returned unchanged.
        Maps each destination string to the properties that should be provided
        to that destination. Props may be:
        - ``'*'``: provide all properties
        - a list of property names to include
        - a list of property names, each prefixed by ``-``, to exclude; all
          others will be provided
        - a single name to include or exclude (prefixed by ``-``)
        - a dict mapping each target property name to its source property name
    dests : list of {str, iterable of str}
        The ordered destination names for the router.  If a set of strings is
        provided for each entry, any of these aliases will route parameters to
        the destination.  The same string may appear in multiple dests entries.
        This should be fixed for each metaestimator.
    default : dict or callable, optional
        This replaces ``routing`` where routing is None.
        This should be fixed for each metaestimator.
    Returns
    -------
    router : callable with signature ``props -> (dest_props, unused)``
    Examples
    --------
    >>> from sklearn_extensions.utils.metaestimators import check_routing
    >>> props = {'foo': [1, 2], 'bar': [3, 4]}
    >>> dests = [['d1', 'all'], ['d2', 'all']]
    >>> # route 'foo' to d1 and 'bar' to d2
    >>> router = check_routing({'d1': 'foo', 'd2': 'bar'}, dests)
    >>> (d1_props, d2_props), unused = router(props)
    >>> d1_props
    {'foo': [1, 2]}
    >>> d2_props
    {'bar': [3, 4]}
    >>> list(unused)
    []
    >>> # rename 'bar' to 'baz'
    >>> router = check_routing({'d1': 'foo', 'd2': {'baz': 'bar'}}, dests)
    >>> (d1_props, d2_props), unused = router(props)
    >>> d1_props
    {'foo': [1, 2]}
    >>> d2_props
    {'baz': [3, 4]}
    >>> list(unused)
    []
    >>> # d2 takes all but foo
    >>> router = check_routing({'d1': 'foo', 'd2': '-foo'}, dests)
    >>> router(props)[0]
    [{'foo': [1, 2]}, {'bar': [3, 4]}]
    >>> # d1 takes all; d2 takes none
    >>> router = check_routing({'d1': '*'}, dests)
    >>> (d1_props, d2_props), unused = router(props)
    >>> sorted(d1_props.items())
    [('bar', [3, 4]), ('foo', [1, 2])]
    >>> d2_props
    {}
    >>> # the 'all' alias distributes to both dests
    >>> router = check_routing({'all': 'foo', 'd1': 'bar'}, dests)
    >>> (d1_props, d2_props), unused = router(props)
    >>> sorted(d1_props.items())
    [('bar', [3, 4]), ('foo', [1, 2])]
    >>> d2_props
    {'foo': [1, 2]}
    >>> # an unused prop
    >>> router = check_routing({'all': 'foo'}, dests)
    >>> (d1_props, d2_props), unused = router(props)
    >>> d1_props
    {'foo': [1, 2]}
    >>> d2_props
    {'foo': [1, 2]}
    >>> list(unused)
    ['bar']
    >>> # a default: both get foo
    >>> router = check_routing(None, dests, default={'all': 'foo'})
    >>> (d1_props, d2_props), unused = router(props)
    >>> d1_props
    {'foo': [1, 2]}
    >>> d2_props
    {'foo': [1, 2]}
    >>> # an overridden default: only d1 gets foo
    >>> router = check_routing({'d1': 'foo'}, dests, default={'all': 'foo'})
    >>> (d1_props, d2_props), unused = router(props)
    >>> d1_props
    {'foo': [1, 2]}
    >>> d2_props
    {}
    """
    if routing is None:
        if default is None:
            raise ValueError('Routing must be specified')
        routing = default
    if not callable(routing):
        return _Router(routing, dests)
    return routing
