"""
hyprat
======

Exact rational hypercomplex numbers (Cayley-Dickson construction:
reals -> complex -> quaternions -> octonions -> ...).

    >>> from hyprat import Hy
    >>> z = Hy('5/2', '-16/5')
    >>> str(z)
    '(5/2-16/5j)'

``Hu`` holds the Hurwitz integers (quaternions with all-integer or
all-half-odd-integer coordinates) as exact Python ints.

See :mod:`hyprat.hypercomplex` and :mod:`hyprat.hurwitz` for the full
implementation and API documentation.
"""

from .hypercomplex import Hy
from .hurwitz import Hu

__all__ = ["Hy", "Hu"]
__version__ = "0.1.0"
