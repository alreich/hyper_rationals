"""
hyprat
======

Exact rational hypercomplex numbers (Cayley-Dickson construction:
reals -> complex -> quaternions -> octonions -> ...).

    >>> from hyprat import Hy
    >>> z = Hy('5/2', '-16/5')
    >>> str(z)
    '(5/2-16/5j)'

See :mod:`hyprat.hypercomplex` for the full implementation and API
documentation.
"""

from .hypercomplex import Hy

__all__ = ["Hy"]
__version__ = "0.1.0"
