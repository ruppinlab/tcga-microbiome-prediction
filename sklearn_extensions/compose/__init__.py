"""
The :mod:`sklearn_extensions.compose` module implements meta-estimators for
building composite models with transformers
"""

from ._column_transformer import ExtendedColumnTransformer


__all__ = ['ExtendedColumnTransformer']
