"""
hypercomplex.py
================

An immutable class, ``Hy``, representing rational-valued hypercomplex
numbers of arbitrary rank, built recursively via the Cayley-Dickson
construction.

    rank 0  ->  a plain ``fractions.Fraction``               (a "real")
    rank 1  ->  Hy(real, imag) with real, imag : Fraction     (a "complex")
    rank 2  ->  Hy(h1, h2)     with h1, h2     : rank-1 Hy    (a "quaternion")
    rank 3  ->  Hy(h3, h4)     with h3, h4     : rank-2 Hy    (an "octonion")
    rank n  ->  Hy(x, y)       with x, y       : rank-(n-1) Hy

Every concrete ``Hy`` instance therefore stores exactly two components,
``real`` and ``imag``, which are *always the same "shape"*: either both
plain ``Fraction``s, or both ``Hy`` instances of the same rank.  The
constructor normalizes (embeds) whatever it is given so that this
invariant always holds -- e.g. ``Hy(3, Hy('1/2', '1/3'))`` will silently
promote the bare ``3`` into a rank-1 ``Hy('3', '0')`` before pairing it
with the rank-1 argument.

This is exactly the classical Cayley-Dickson doubling construction
(the one that turns R -> C -> H -> O -> sedenions -> ...), specialized
to exact rational coefficients.  Multiplication, conjugation, norms,
inverses and division are all defined by the standard recursive
formulas, so the algebra automatically "does the right thing" for
complex, quaternion and octonion values, and continues to make sense
(algebraically, if no longer as a division algebra) at higher ranks.

Quick tour
----------

    >>> from hyprat import Hy
    >>> z = Hy('5/2', '-16/5')          # a rational complex number
    >>> str(z)
    '(5/2-16/5j)'
    >>> q = Hy(Hy(0, 0), Hy(1, 0))      # the quaternion j
    >>> str(q)
    '(j)'
    >>> str(q * q)
    '(-1)'

Because ``Hy`` always stores exactly two components, a "plain real
number" is represented as a rank-1 ``Hy`` with a zero imaginary part
(exactly the way Python's own ``complex`` behaves: ``str(complex(3))``
is ``'(3+0j)'``).

String forms
------------

``str()`` renders a value the way Python renders ``complex`` numbers
for rank 1 (using ``j``), the customary ``a+bi+cj+dk`` notation for
rank 2, ``i, j, k, L, iL, jL, kL`` basis labels for rank 3 ("octonions"),
and ``e1 .. e(2**rank - 1)`` (etc.) imaginary units for rank >= 4::

    str(Hy('5/2', '16/5'))                 -> '(5/2+16/5j)'
    str(Hy(Hy(1, 2), Hy(3, 4)))            -> '(1+2i+3j+4k)'
    str(Hy(Hy(Hy(1,0),Hy(0,0)), Hy(Hy(0,1),Hy(0,0))))
                                            -> '(1+iL)'   (etc.)

``Hy.parse(s)`` is the inverse of ``str`` (and is also used
automatically by the constructor -- see below).

``repr()`` returns Python source that reconstructs an equal value,
e.g. ``Hy('5/2', '-16/5')``.

Random values and flat-array conversion
----------------------------------------

``Hy.random(rank)`` produces a random rank-``rank`` value; see
``Hy.seed()``/the ``rng=``/``seed=`` keywords on ``Hy.random`` for how to
make the sequence reproducible.

``Hy.from_array([...])`` builds a Hy from a flat list of ``2**rank``
coefficients (ints, floats, ``Fraction``s and/or fraction strings like
``'5/2'``, freely mixed); ``some_hy.to_array()`` is the inverse::

    >>> Hy.from_array([1, '2/3', 3.5, '-1/4'])
    Hy(Hy('1', '2/3'), Hy('7/2', '-1/4'))
    >>> Hy(Hy(1, 2), Hy(3, 4)).to_array(as_str=True)
    ['1', '2', '3', '4']

Signatures: split-complex numbers, split-quaternions, split-octonions, ...
---------------------------------------------------------------------------

Every doubling step of the Cayley-Dickson construction carries a nonzero
rational parameter ``mu``::

    (a, b)(c, d) = (a c + mu conj(d) b,  d a + b conj(c))

``mu = -1`` -- the default, so everything above is unchanged -- gives the
classical tower.  ``mu = +1`` at a step gives the *split* version of that
step, and any other nonzero rational gives the corresponding generalized
algebra (for instance the quaternion algebra ``(-1, 2)``)::

    split-complex         signs (1,)            j*j = +1
    split-quaternions     signs (-1, 1)         i*i = -1,  j*j = k*k = +1
    split-octonions       signs (-1, -1, 1)     i, j, k -> -1;  L, iL, jL, kL -> +1

``Hy(real, imag, mu=...)`` sets the ``mu`` of the *top* doubling level
(default ``-1``); each component already carries its own lower levels, so
the full signature -- one ``mu`` per level, lowest first -- is
``some_hy.signs`` and the top one is ``some_hy.mu``.  The factory methods
(``Hy.units``, ``Hy.random``, ``Hy.from_array``, ``Hy.parse``,
``Hy.from_matrix``) take a ``signs=`` tuple, or the name of a preset
(``"split-complex"``, ``"split-quaternion"``, ``"split-octonion"``,
``"complex"``, ``"quaternion"``, ``"octonion"``, ``"sedenion"``; also
available as ``Hy.SPLIT_QUATERNION`` and so on)::

    >>> j = Hy(0, 1, mu=1)               # the split-complex j
    >>> str(j * j)
    '(1)'
    >>> p, m = Hy(1, 1, mu=1), Hy(1, -1, mu=1)
    >>> str(p * m)                       # zero divisors: (1+j)(1-j) == 0
    '(0)'
    >>> p.is_null()
    True
    >>> Hy.SPLIT_QUATERNION
    (-1, 1)
    >>> u = Hy.units(signs="split-quaternion")
    >>> [str(u[name] * u[name]) for name in ("i", "j", "k")]
    ['(-1)', '(1)', '(1)']
    >>> repr(j)                          # repr shows mu only when it isn't -1
    "Hy('0', '1', mu=1)"

The unit *labels* are positional -- they name a place in the Cayley-Dickson
tower (``i*j == k`` and ``iL == i*L`` hold in every signature) -- so they
do not change with ``mu``; what changes is what each unit squares to (see
``Hy.unit_square``).  ``str()`` is therefore the same in every algebra, and
text read back with ``Hy.parse`` needs ``signs=`` to say which algebra it
is in.

With a positive ``mu`` somewhere, the quadratic form ``x * conj(x)``
(``some_hy.norm_squared()``, also still available as ``norm()``) is
*indefinite*: ``N(a, b) = N(a) - mu N(b)``.  It can be negative or zero for
a nonzero value, so ``abs(x)`` raises ``ValueError`` when it is negative,
and a nonzero *null* element (``some_hy.is_null()``) has no inverse -- it
raises ``ZeroDivisionError``.  The form stays multiplicative through rank
3, as in every composition algebra, split or not.

Values of different signatures never mix: ``+ - * /`` raise ``ValueError``
(a lower-rank value embeds into a higher-rank one only when its signature
is a prefix of the other's), ``==`` is simply ``False``, and plain
numbers (and text such as ``'1+j'``, read in the other operand's algebra)
mix freely.  The conversions to other quaternion packages, and
``complex()``, accept only the classical algebras.

Interoperability with other quaternion packages
------------------------------------------------

Rank-2 ``Hy`` values (and rank-1 "complex" values, which are embedded
as ``a + b*i``) can be converted to and from the quaternion types of
three third-party packages.  None of these packages is a dependency of
``hyprat``; each is imported only when a conversion method that needs
it is called::

    some_hy.to_sympy()                     Hy.from_sympy(q)
    some_hy.to_numpy_quaternion()          Hy.from_numpy_quaternion(q)
    some_hy.to_quaternionic()              Hy.from_quaternionic(q)

SymPy, like ``Hy``, represents rational numbers exactly, so
``to_sympy()``/``from_sympy()`` round-trip exactly with no rounding at
all (provided the SymPy quaternion's own components are exact -- see
:meth:`Hy.from_sympy`).  ``numpy-quaternion`` and ``quaternionic``, by
contrast, both store ``float64`` coordinates, so a ``Hy -> other ->
Hy`` round trip through *those* two is only exact if you ask for it
(see the ``exact=`` and ``max_denominator=`` keywords of their
``from_*`` methods).

Units and LaTeX rendering
--------------------------

``Hy.units(rank)`` returns a dict of every unit element (``+-1``,
``+-j``, ``+-i``, ...) of the rank-``rank`` algebra, keyed by their
string form; ``some_hy.is_unit()`` reports whether a value is one of
the units for its own rank::

    >>> Hy.units(1)
    {'1': Hy('1', '0'), '-1': Hy('-1', '0'), 'j': Hy('0', '1'), '-j': Hy('0', '-1')}
    >>> Hy(0, 1).is_unit()
    True

``some_hy.latex()`` renders a value as a LaTeX math expression, for
display in a Jupyter notebook (e.g. via ``IPython.display.Math``).
Basis units at rank 1-3 (``j``; ``i, j, k``; ``i, j, k, L, iL, jL,
kL``) render unchanged; the ``e1, e2, ...`` labels used from rank 4
up are subscripted (``e_{1}``, ``e_{2}``, ...); non-integer
coefficients are rendered as ``\\frac{num}{den}``
(the default, ``vinculum='horizontal'``) or as a plain slash,
``num/den`` (``vinculum='diagonal'``)::

    >>> Hy('5/2', '-16/5').latex()
    '\\\\frac{5}{2}-\\\\frac{16}{5}j'
    >>> Hy('5/2', '-16/5').latex(vinculum='diagonal')
    '5/2-16/5j'

Matrix (regular) representation
--------------------------------

``some_hy.to_matrix()`` returns the *regular representation* of a
value as a ``2**rank x 2**rank`` ``numpy`` array of exact
``Fraction``s -- the classical "complex numbers as 2x2 real matrices"
/ "quaternions as 4x4 real matrices" construction, generalized to any
rank via the value's own multiplication (column ``i`` is ``self *
units[i]``). ``Hy.from_matrix()`` is the inverse. This requires
``numpy`` (not a runtime dependency of ``hyprat``; imported lazily).

For rank 0-2 (real, complex, quaternion) this is a genuine algebra
isomorphism: ``M(x) @ M(y) == M(x * y)``. Rank >= 3 (octonions and
beyond) is **not** associative, so the same equality fails there --
``to_matrix()``/``from_matrix()`` raise ``ValueError`` at rank >= 3
unless ``allow_nonassociative=True`` is passed::

    >>> Hy('2', '3').to_matrix()  # doctest: +SKIP
    array([[Fraction(2, 1), Fraction(-3, 1)],
           [Fraction(3, 1), Fraction(2, 1)]], dtype=object)
"""

from __future__ import annotations

import importlib
import math
import numbers
import random
import re
from fractions import Fraction


__all__ = ["Hy"]


# --------------------------------------------------------------------------
# A private sentinel used to distinguish "no imag argument was given" from
# "imag was explicitly given as 0".
# --------------------------------------------------------------------------
class _Missing:
    def __repr__(self):
        return "<missing>"


_MISSING = _Missing()

_ScalarLike = (int, float, Fraction, str)

# --------------------------------------------------------------------------
# Cayley-Dickson parameters ("signs").  Every doubling step of the tower
# carries a nonzero rational parameter mu:
#
#     (a, b)(c, d) = (a c + mu * conj(d) b,  d a + b conj(c))
#
# mu = -1 (the default) gives the classical tower R, C, H, O, ...; mu = +1
# at a step gives the *split* version of that step (split-complex numbers,
# split-quaternions, split-octonions, ...); any other nonzero rational gives
# the corresponding generalized algebra (e.g. quaternion algebras (-1, 2)).
# --------------------------------------------------------------------------
_DEFAULT_MU = Fraction(-1)

# Named signature presets, as tuples ordered from the lowest doubling level
# to the highest.  Names are matched case-insensitively, with '-', '_' and
# ' ' interchangeable.
_SIGN_PRESETS = {
    "complex": (-1,),
    "quaternion": (-1, -1),
    "octonion": (-1, -1, -1),
    "sedenion": (-1, -1, -1, -1),
    "split-complex": (1,),
    "split-quaternion": (-1, 1),
    "split-octonion": (-1, -1, 1),
}

# --------------------------------------------------------------------------
# A module-wide default RNG used by Hy.random() whenever the caller doesn't
# supply its own random.Random instance or an explicit one-off seed. Call
# Hy.seed(value) to make subsequent Hy.random(...) calls reproducible.
# --------------------------------------------------------------------------
_default_rng = random.Random()


class Hy:
    """An immutable rational hypercomplex number (Cayley-Dickson tower).

    ``Hy`` is deliberately *not* a ``@dataclass``: construction has to
    normalize/embed its two arguments so that ``real`` and ``imag`` end
    up as the same "shape" (see module docstring), which is more than a
    dataclass's generated ``__init__`` can do on its own.  Instances are
    immutable (``__slots__`` + a locked-down ``__setattr__``) and
    hashable.
    """

    __slots__ = ("_real", "_imag", "_mu")

    # Named signature presets (tuples of mu values, lowest level first),
    # usable anywhere a ``signs=`` argument is accepted, either as these
    # tuples or by name, e.g. ``signs="split-quaternion"``.
    COMPLEX = (-1,)
    QUATERNION = (-1, -1)
    OCTONION = (-1, -1, -1)
    SEDENION = (-1, -1, -1, -1)
    SPLIT_COMPLEX = (1,)
    SPLIT_QUATERNION = (-1, 1)
    SPLIT_OCTONION = (-1, -1, 1)

    # ---------------------------------------------------------------- #
    # Construction
    # ---------------------------------------------------------------- #
    def __init__(self, real, imag=_MISSING, *, mu=None):
        mu_c = None if mu is None else _coerce_mu(mu)
        real_c = _coerce_component(real)

        if imag is _MISSING:
            if isinstance(real_c, Hy):
                # Hy(some_hy) is a copy/identity construction: it takes
                # on the value of `some_hy` as-is, at whatever rank that
                # already is (rather than promoting it one rank higher),
                # and with its own signature.
                if mu_c is not None and mu_c != real_c._mu:
                    raise ValueError(
                        f"cannot change the top-level mu of an existing Hy "
                        f"(it is {real_c._mu}, but mu={mu_c} was requested); "
                        "build a new value from its components instead, or "
                        "use Hy.parse(..., signs=...) / Hy.from_array(..., "
                        "signs=...)"
                    )
                object.__setattr__(self, "_real", real_c._real)
                object.__setattr__(self, "_imag", real_c._imag)
                object.__setattr__(self, "_mu", real_c._mu)
                return
            imag_c = Fraction(0)
        else:
            imag_c = _coerce_component(imag)

        # The two components must live in the same algebra: the one of
        # lower rank is embedded into the other, which is only possible if
        # its signature is a prefix of the other's.
        sr, si = _signs(real_c), _signs(imag_c)
        lower = sr if len(sr) >= len(si) else si
        object.__setattr__(self, "_real", _embed_signs(real_c, lower))
        object.__setattr__(self, "_imag", _embed_signs(imag_c, lower))
        object.__setattr__(
            self, "_mu", _DEFAULT_MU if mu_c is None else mu_c
        )

    @classmethod
    def _make(cls, real, imag, mu=_DEFAULT_MU) -> "Hy":
        """Internal fast constructor: assumes `real`/`imag` already have
        matching rank and signature and skips all normalization. Not for
        public use."""
        obj = object.__new__(cls)
        object.__setattr__(obj, "_real", real)
        object.__setattr__(obj, "_imag", imag)
        object.__setattr__(obj, "_mu", mu)
        return obj

    def __setattr__(self, name, value):
        raise AttributeError(
            f"Hy instances are immutable; cannot set attribute {name!r}"
        )

    def __delattr__(self, name):
        raise AttributeError("Hy instances are immutable")

    # ---------------------------------------------------------------- #
    # Accessors
    # ---------------------------------------------------------------- #
    @property
    def real(self):
        """The first component (a Fraction, or a Hy of the same rank as imag)."""
        return self._real

    @property
    def imag(self):
        """The second component (a Fraction, or a Hy of the same rank as real)."""
        return self._imag

    @property
    def rank(self) -> int:
        """0 for a bare Fraction; 1 for complex; 2 for quaternion; 3 for
        octonion; etc. (dimension of the algebra is 2**rank)."""
        return 1 + max(_rank(self._real), _rank(self._imag))

    @property
    def mu(self) -> Fraction:
        """The Cayley-Dickson parameter of this value's *top* doubling
        level, a nonzero ``Fraction`` (``-1``, the default, for the
        classical algebras; ``+1`` for a split step). The squares of the
        two halves' units are related by ``mu`` through the product
        ``(a, b)(c, d) = (a c + mu * conj(d) b, d a + b conj(c))``."""
        return self._mu

    @property
    def signs(self) -> tuple:
        """The full signature of this value's algebra: a tuple of ``rank``
        nonzero ``Fraction``s, one per doubling level, ordered from the
        lowest level to the highest (so ``signs[-1] == mu``). For the
        classical algebras every entry is ``-1``; for the split-quaternions
        it is ``(-1, 1)``."""
        return _signs(self)

    @property
    def dimension(self) -> int:
        """Number of real (Fraction) coordinates: 2**rank."""
        return 2 ** self.rank

    def components(self) -> tuple:
        """All of the real (Fraction) coordinates, in the canonical
        Cayley-Dickson order (e.g. for a quaternion: 1, i, j, k)."""
        return tuple(_flatten(self))

    def conjugate(self) -> "Hy":
        return conj(self)

    def inverse(self) -> "Hy":
        return inverse(self)

    def norm(self) -> Fraction:
        """The *squared* norm, computed exactly as a Fraction.

        This is the quadratic form ``x * conj(x)`` of the algebra. In the
        classical algebras (every mu equal to -1) it is the sum of the
        squares of the real coordinates, and so is never negative; with a
        positive mu somewhere it is *indefinite* and can be negative or
        zero for a nonzero value. Identical to :meth:`norm_squared`
        (which is the clearer name); kept under this name for backward
        compatibility. For the (float) square root see ``abs(x)``.
        """
        return abs2(self)

    def norm_squared(self) -> Fraction:
        """The exact quadratic form ``x * conj(x)``, a ``Fraction``.

        Recursively, ``N(a, b) = N(a) - mu * N(b)`` with ``N(r) = r*r`` for
        a rational ``r``; when every mu is ``-1`` this is the usual sum of
        squares. In a split algebra (or any algebra with a positive mu) it
        is indefinite: its sign tells you whether a value is *timelike*
        (positive), *spacelike* (negative) or *null* (zero). Multiplicative
        (``N(x*y) == N(x)*N(y)``) for ranks 0-3, i.e. for the composition
        algebras, split or not.
        """
        return abs2(self)

    def is_null(self) -> bool:
        """True iff this value is nonzero but has ``norm_squared() == 0``
        -- a *null* element, which exists only when the norm is
        indefinite (e.g. ``1 + j`` in the split-complex numbers).

        Null elements are not invertible. For ranks 1-3 (the composition
        algebras) they are exactly the zero divisors; from rank 4 up there
        are also zero divisors with nonzero norm, which this does not
        detect.
        """
        return (not _is_zero_val(self)) and abs2(self) == 0

    def is_zero(self) -> bool:
        return _is_zero_val(self)

    def is_unit(self) -> bool:
        """True iff this value is one of the *units* of its own rank --
        i.e. exactly one of its ``2**rank`` real coordinates is ``+-1``
        and every other coordinate is 0. Equivalent to (but cheaper
        than) checking membership in ``Hy.units(self.rank).values()``.

        "Unit" here means *basis unit*, not "invertible element": in a
        split algebra there are many invertible values that are not
        units in this sense, and a unit can square to ``+1`` rather than
        ``-1``.

            >>> Hy(0, 1).is_unit()
            True
            >>> Hy(1, 1).is_unit()
            False
            >>> Hy(0, 0).is_unit()
            False
        """
        nonzero = [c for c in _flatten(self) if c != 0]
        return len(nonzero) == 1 and abs(nonzero[0]) == 1

    # ---------------------------------------------------------------- #
    # Arithmetic
    # ---------------------------------------------------------------- #
    def _coerce_other(self, other):
        if isinstance(other, Hy):
            return other
        if isinstance(other, complex):
            # A Python complex is an ordinary (mu = -1) complex number.
            return Hy(Fraction(str(other.real)), Fraction(str(other.imag)))
        if isinstance(other, _ScalarLike):
            try:
                if isinstance(other, str):
                    try:
                        return Fraction(other)
                    except ValueError:
                        # Text such as '1+2j' carries no signature of its
                        # own, so it is read in *this* value's algebra.
                        return _parse(other, signs=self.signs, pad=True)
                return _coerce_component(other)
            except (TypeError, ValueError):
                return NotImplemented
        return NotImplemented

    def __add__(self, other):
        other_h = self._coerce_other(other)
        if other_h is NotImplemented:
            return NotImplemented
        return add(self, other_h)

    __radd__ = __add__

    def __sub__(self, other):
        other_h = self._coerce_other(other)
        if other_h is NotImplemented:
            return NotImplemented
        return sub(self, other_h)

    def __rsub__(self, other):
        other_h = self._coerce_other(other)
        if other_h is NotImplemented:
            return NotImplemented
        return sub(other_h, self)

    def __mul__(self, other):
        other_h = self._coerce_other(other)
        if other_h is NotImplemented:
            return NotImplemented
        return mul(self, other_h)

    def __rmul__(self, other):
        other_h = self._coerce_other(other)
        if other_h is NotImplemented:
            return NotImplemented
        return mul(other_h, self)

    def __truediv__(self, other):
        other_h = self._coerce_other(other)
        if other_h is NotImplemented:
            return NotImplemented
        return div(self, other_h)

    def __rtruediv__(self, other):
        other_h = self._coerce_other(other)
        if other_h is NotImplemented:
            return NotImplemented
        return div(other_h, self)

    def __neg__(self) -> "Hy":
        return neg(self)

    def __pos__(self) -> "Hy":
        return self

    def __abs__(self) -> float:
        """The (float) square root of :meth:`norm_squared`.

        Raises ``ValueError`` when ``norm_squared()`` is negative (only
        possible when the norm is indefinite, i.e. with a positive mu
        somewhere in the signature), since there is then no real square
        root; use :meth:`norm_squared` instead. Returns ``0.0`` for a
        null element.
        """
        n = abs2(self)
        if n < 0:
            raise ValueError(
                f"norm_squared() is negative ({n}), so abs() is undefined; "
                "use norm_squared() for the indefinite quadratic form"
            )
        return math.sqrt(float(n))

    def __pow__(self, n):
        if not isinstance(n, int):
            return NotImplemented
        if n == 0:
            return _embed_signs(Fraction(1), _signs(self))
        base = self if n > 0 else self.inverse()
        result = _embed_signs(Fraction(1), _signs(base))
        for _ in range(abs(n)):
            result = mul(result, base)
        return result

    def __bool__(self) -> bool:
        return not self.is_zero()

    # ---------------------------------------------------------------- #
    # Comparisons, hashing, iteration
    # ---------------------------------------------------------------- #
    def __eq__(self, other):
        other_h = self._coerce_other(other)
        if other_h is NotImplemented:
            return NotImplemented
        return _values_equal(self, other_h)

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __hash__(self):
        return hash(_to_nested_tuple(_canonical_trim(self)))

    def __iter__(self):
        return iter((self._real, self._imag))

    def __getitem__(self, idx):
        return (self._real, self._imag)[idx]

    def __len__(self):
        return 2

    def __complex__(self):
        if self.signs[0] != _DEFAULT_MU:
            raise ValueError(
                f"{self!r} is not an ordinary complex number "
                f"(its first doubling level has mu={self.signs[0]}, not -1)"
            )
        flat = _flatten(self)
        if any(c != 0 for c in flat[2:]):
            raise ValueError(f"{self!r} is not complex-valued (rank > 1)")
        re_, im_ = (flat + [Fraction(0), Fraction(0)])[:2]
        return complex(float(re_), float(im_))

    # ---------------------------------------------------------------- #
    # String forms
    # ---------------------------------------------------------------- #
    def __str__(self) -> str:
        coeffs = _flatten(self)
        labels = _basis_labels(self.rank)
        return f"({_format_terms(coeffs, labels)})"

    def __repr__(self) -> str:
        def comp_repr(x):
            if isinstance(x, Fraction):
                return repr(str(x))
            return repr(x)

        mu_part = ""
        if self._mu != _DEFAULT_MU:
            m = self._mu
            mu_part = f", mu={m.numerator}" if m.denominator == 1 else f", mu={str(m)!r}"
        return f"Hy({comp_repr(self._real)}, {comp_repr(self._imag)}{mu_part})"

    # ---------------------------------------------------------------- #
    # Parsing
    # ---------------------------------------------------------------- #
    @classmethod
    def parse(cls, s: str, *, signs=None) -> "Hy":
        """Parse the kind of string produced by ``str(some_hy)`` -- e.g.
        ``'(5/2-16/5j)'``, ``'1+2i+3j+4k'``, ``'1-k+2kL'`` -- into a Hy.

        The text itself says nothing about the algebra, so by default
        the result is read in the classical one (every mu equal to -1).
        Pass ``signs=`` (a tuple of nonzero rationals ordered from the
        lowest doubling level to the highest, or a preset name such as
        ``"split-quaternion"``) to read it in another algebra; its length
        must equal the rank implied by the text::

            >>> Hy.parse('1+2j+3k', signs="split-quaternion").signs
            (Fraction(-1, 1), Fraction(1, 1))
        """
        return _parse(s, signs=signs)

    from_string = parse  # convenient alias

    # ---------------------------------------------------------------- #
    # Randomness
    # ---------------------------------------------------------------- #
    @classmethod
    def seed(cls, seed_value=None) -> None:
        """Seed the module-wide default RNG used by :meth:`random`
        whenever *that* call isn't given its own ``rng=`` or ``seed=``.

        Call this once (e.g. at the top of a script or a test module) to
        make an entire sequence of subsequent, argument-free
        ``Hy.random(rank)`` calls reproducible. Pass ``seed_value=None``
        (the default) to reseed unpredictably from OS entropy, exactly
        like ``random.seed(None)``.
        """
        _default_rng.seed(seed_value)

    @classmethod
    def random(
        cls,
        rank: "int | None" = None,
        *,
        signs=None,
        lo: int = -9,
        hi: int = 9,
        dmax: int = 6,
        rng: "random.Random | None" = None,
        seed=None,
    ) -> "Hy":
        """A random rank-``rank`` Hy.

        Each of the ``2**rank`` real coefficients is an independent random
        ``Fraction(n, d)``, with the numerator ``n`` drawn uniformly from
        the inclusive range ``[lo, hi]`` and the denominator ``d`` drawn
        uniformly from ``[1, dmax]``.

        There are three, mutually exclusive ways to control reproducibility:

        * do nothing -- draws from the module-wide default RNG, which is
          shared and unseeded (i.e. non-reproducible) unless...
        * ...you've called ``Hy.seed(value)`` beforehand, which reseeds
          that shared default RNG once, making every subsequent
          argument-free ``Hy.random(rank)`` call reproducible; or
        * pass ``seed=value`` to this call for a one-off, freshly created
          ``random.Random(value)`` used just for this single call, without
          touching any shared state; or
        * pass ``rng=some_random.Random_instance`` to fully control (and
          optionally share across several calls) the random stream
          yourself.

        ``signs=`` selects the algebra (a tuple of nonzero rationals,
        lowest doubling level first, or a preset name such as
        ``"split-octonion"``); it defaults to the classical one (all
        ``-1``). If ``rank`` is omitted it is taken from ``len(signs)``.

        Raises ``ValueError`` if ``rank`` isn't a positive int (or can't
        be determined), if it disagrees with ``signs``, or if both
        ``rng`` and ``seed`` are given.
        """
        rank, signs = _resolve_rank_and_signs(rank, signs, min_rank=1)
        if rng is not None and seed is not None:
            raise ValueError("pass either `rng` or `seed`, not both")
        if seed is not None:
            rng = random.Random(seed)
        elif rng is None:
            rng = _default_rng

        def rand_coeff() -> Fraction:
            return Fraction(rng.randint(lo, hi), rng.randint(1, dmax))

        def build(r: int):
            if r == 0:
                return rand_coeff()
            return Hy._make(build(r - 1), build(r - 1), signs[r - 1])

        return build(rank)

    # ---------------------------------------------------------------- #
    # Flat-array conversion
    # ---------------------------------------------------------------- #
    @classmethod
    def from_array(cls, coeffs, *, signs=None) -> "Hy":
        """Build a Hy from a flat sequence of its ``2**rank`` real
        coefficients, in the same Cayley-Dickson order used by
        :meth:`components` / :meth:`to_array` (e.g., for a quaternion:
        1, i, j, k).

        Each element may be an ``int``, ``float``, ``Fraction``, or a
        string holding a plain fraction such as ``'5/2'`` or a decimal
        such as ``'3.2'`` -- these types may be freely mixed within a
        single call, e.g.::

            Hy.from_array([1, '2/3', 3.5, '-1/4'])   # a quaternion

        ``len(coeffs)`` must be a power of 2 that is >= 2 (2 -> rank 1
        "complex", 4 -> rank 2 "quaternion", 8 -> rank 3 "octonion", etc.)
        since every Hy has rank >= 1.

        ``signs=`` selects the algebra the coefficients live in (a tuple
        of nonzero rationals, lowest doubling level first, or a preset
        name such as ``"split-quaternion"``); its length must equal the
        rank. The default is the classical algebra (all ``-1``).
        """
        values = list(coeffs)
        n = len(values)
        if n < 2 or (n & (n - 1)) != 0:
            raise ValueError(
                "from_array() requires a length that is a power of 2 "
                f"and at least 2; got {n}"
            )
        rank = n.bit_length() - 1
        _, signs = _resolve_rank_and_signs(rank, signs, min_rank=1)
        fracs = [_coerce_flat_element(v) for v in values]
        return _unflatten(fracs, rank, signs)

    def to_array(self, as_str: bool = False) -> list:
        """This Hy's ``2**rank`` real coefficients, flattened into a plain
        Python list in Cayley-Dickson order -- the inverse of
        :meth:`from_array`.

        By default the entries are ``Fraction`` objects; pass
        ``as_str=True`` to get their string form instead (e.g. ``'5/2'``),
        which is convenient for JSON or other text-based serialization,
        and which :meth:`from_array` will happily read back in.
        """
        coeffs = _flatten(self)
        return [str(c) for c in coeffs] if as_str else list(coeffs)

    # ---------------------------------------------------------------- #
    # Units
    # ---------------------------------------------------------------- #
    @classmethod
    def units(cls, rank: "int | None" = None, *, signs=None) -> dict:
        """The unit elements of the rank-``rank`` algebra: the ``2 *
        2**rank`` values with exactly one real coordinate equal to
        ``+-1`` and every other coordinate 0 (e.g. for rank 1: ``+-1``
        and ``+-j``; for rank 2: ``+-1, +-i, +-j, +-k``).

        Returns a dict mapping each unit's string form (matching
        :func:`str`) to the value itself, in ``+1, -1, +unit, -unit,
        ...`` order::

            >>> Hy.units(1)
            {'1': Hy('1', '0'), '-1': Hy('-1', '0'), 'j': Hy('0', '1'), '-j': Hy('0', '-1')}

        For ``rank == 0`` -- the base case where a hypercomplex value is
        just a plain ``Fraction`` rather than a ``Hy`` -- the two units
        ``Fraction(1)`` and ``Fraction(-1)`` are returned directly
        (every other rank returns ``Hy`` instances).

        ``signs=`` selects the algebra (a tuple of nonzero rationals,
        lowest doubling level first, or a preset name such as
        ``"split-quaternion"``); if ``rank`` is omitted it is taken from
        ``len(signs)``. The *labels* are positional (they name a place in
        the Cayley-Dickson tower), so they do not change with ``signs``;
        what changes is what each unit squares to -- see
        :meth:`unit_square`.

        See also :meth:`is_unit`. Raises ``ValueError`` if ``rank``
        isn't a non-negative int (or can't be determined), or disagrees
        with ``signs``.
        """
        rank, signs = _resolve_rank_and_signs(rank, signs, min_rank=0)
        if rank == 0:
            return {"1": Fraction(1), "-1": Fraction(-1)}
        n = 2 ** rank
        labels = _basis_labels(rank)
        result = {}
        for idx, lbl in enumerate(labels):
            pos_key = lbl if lbl else "1"
            neg_key = f"-{lbl}" if lbl else "-1"
            coeffs = [Fraction(0)] * n
            coeffs[idx] = Fraction(1)
            result[pos_key] = _unflatten(coeffs, rank, signs)
            coeffs[idx] = Fraction(-1)
            result[neg_key] = _unflatten(coeffs, rank, signs)
        return result

    @staticmethod
    def unit_square(label_index: int, signs) -> Fraction:
        """What the basis unit with coordinate index ``label_index``
        (``0`` is the real unit ``1``, ``1`` is ``j``/``i``, and so on, in
        the order of :meth:`components`) squares to, in the algebra with
        the given ``signs`` (a tuple or preset name), as an exact
        ``Fraction``.

        The index's binary digits say which doubling levels the unit
        involves; if ``S`` is that set of levels, the unit squares to
        ``-prod(-mu_l for l in S)`` (and ``1`` squares to ``1``). With all
        ``mu = -1`` every imaginary unit squares to ``-1``::

            >>> Hy.unit_square(3, "split-quaternion")     # k, with signs (-1, 1)
            Fraction(1, 1)
            >>> Hy.unit_square(1, "split-quaternion")     # i
            Fraction(-1, 1)
        """
        _, sg = _resolve_rank_and_signs(None, signs, min_rank=1)
        if (not isinstance(label_index, int) or isinstance(label_index, bool)
                or not 0 <= label_index < 2 ** len(sg)):
            raise ValueError(
                f"label_index must be an int in [0, {2 ** len(sg) - 1}], "
                f"got {label_index!r}"
            )
        if label_index == 0:
            return Fraction(1)
        prod = Fraction(1)
        for level in range(len(sg)):
            if label_index >> level & 1:
                prod *= -sg[level]
        return -prod

    # ---------------------------------------------------------------- #
    # Matrix (regular) representation
    # ---------------------------------------------------------------- #
    def to_matrix(self, kind: str = "left", *, as_float: bool = False,
                  allow_nonassociative: bool = False):
        """The *regular representation* of this value as a ``2**rank x
        2**rank`` ``numpy`` array: the matrix of the linear map ``y ->
        self * y`` (``kind="left"``, the default) or ``y -> y * self``
        (``kind="right"``), acting on the ``2**rank``-dimensional real
        vector space of rank-``rank`` values.

        Column ``i`` of the matrix is simply ``(self * e_i).to_array()``
        (or ``(e_i * self).to_array()`` for ``kind="right"``), where
        ``e_i`` is the ``i``-th positive unit from :meth:`units` -- so
        this is literally "the images of the basis vectors", the
        textbook definition of the matrix of a linear map, computed
        directly from ``Hy``'s own multiplication rather than from any
        derived closed-form formula.

        Entries are exact ``fractions.Fraction`` objects (``dtype=object``)
        by default; pass ``as_float=True`` to get an ordinary ``float64``
        array instead (convenient for feeding to ``numpy.linalg``, but no
        longer exact). Requires ``numpy`` (``pip install numpy``), which
        is *not* a runtime dependency of ``hyprat`` and is imported lazily.

        For an associative algebra (rank 0-2: real, complex, quaternion)
        this matrix representation is a genuine algebra homomorphism::

            M(x) @ M(y) == M(x * y)

        and in particular ``det(M(x)) == x.norm() ** (2 ** (rank - 1))``.
        This is exactly the classical "complex numbers as 2x2 real
        matrices" / "quaternions as 4x4 real matrices" construction.

        **Rank >= 3 (octonions and beyond) are not associative**, so
        ``M(x) @ M(y) != M(x * y)`` in general there: the map ``x ->
        M(x)`` is still an injective *linear* embedding (a faithful
        vector-space representation), but it is **not** a ring
        homomorphism, because ordinary matrix multiplication is
        associative and ``Hy`` multiplication isn't past rank 2. (From
        rank 4, sedenions, up, ``Hy`` also has zero divisors, so ``M(x)``
        can even be *singular* -- i.e. ``det(M(x)) == 0`` -- for some
        nonzero ``x``.) Because this is a common trap, ``to_matrix()``
        raises ``ValueError`` at rank >= 3 unless you pass
        ``allow_nonassociative=True`` to build the (non-homomorphic)
        matrix anyway. See the "Matrix representation" section of the
        docs for the literature on faithful, non-matrix-multiplication
        representations of octonions (e.g. Zorn vector matrices).

        Examples
        --------
            >>> import numpy as np  # doctest: +SKIP
            >>> Hy('2', '3').to_matrix()  # doctest: +SKIP
            array([[Fraction(2, 1), Fraction(-3, 1)],
                   [Fraction(3, 1), Fraction(2, 1)]], dtype=object)
            >>> x, y = Hy.random(2), Hy.random(2)  # doctest: +SKIP
            >>> np.array_equal(  # doctest: +SKIP
            ...     np.dot(x.to_matrix(), y.to_matrix()), (x * y).to_matrix()
            ... )
            True

        See also :meth:`from_matrix`, the inverse of this method.
        """
        if kind not in ("left", "right"):
            raise ValueError(f"kind must be 'left' or 'right', got {kind!r}")
        rank = self.rank
        if rank >= 3 and not allow_nonassociative:
            raise ValueError(
                f"rank-{rank} values (octonions and beyond) are not "
                "associative, so their regular-representation matrix is a "
                "faithful *linear* embedding but NOT an algebra "
                "homomorphism: M(x) @ M(y) != M(x * y) in general, and "
                "(from rank 4 up) M(x) can even be singular for nonzero x. "
                "Pass allow_nonassociative=True to build the matrix anyway, "
                "or see the 'Matrix representation' section of the docs "
                "for faithful alternatives (e.g. Zorn vector matrices)."
            )
        np = _import_optional("numpy", "numpy")
        basis = _positive_units(rank, self.signs)
        if kind == "left":
            columns = [(self * e).to_array() for e in basis]
        else:
            columns = [(e * self).to_array() for e in basis]
        matrix = np.array(columns, dtype=object).T
        if as_float:
            matrix = matrix.astype(float)
        return matrix

    @classmethod
    def from_matrix(cls, matrix, *, kind: str = "left", rank: "int | None" = None,
                     allow_nonassociative: bool = False, validate: bool = False,
                     signs=None) -> "Hy":
        """Build a ``Hy`` from its regular-representation matrix, the
        inverse of :meth:`to_matrix`.

        Column 0 of a regular-representation matrix is always ``self *
        1 == self`` (or, for ``kind="right"``, ``1 * self == self``), so
        reconstruction is just ``Hy.from_array(matrix[:, 0])``; ``kind``
        only matters when ``validate=True`` (see below).

        ``rank`` is inferred from ``matrix``'s shape (which must be
        ``2**r x 2**r`` for some ``r >= 1``) if omitted; passing an
        explicit ``rank`` that disagrees with the shape raises
        ``ValueError``. As with :meth:`to_matrix`, rank >= 3 raises
        ``ValueError`` unless ``allow_nonassociative=True``, since a
        rank >= 3 matrix isn't the image of a true algebra homomorphism
        in the first place.

        ``signs=`` selects the algebra the matrix represents (default: the
        classical one); a regular-representation matrix does not record
        it, so a split-complex or split-quaternion matrix must be read
        back with the matching ``signs=``.

        ``validate=True`` rebuilds every other column from the
        recovered value and checks it against ``matrix`` exactly
        (``Fraction`` equality) -- a cheap way to catch a matrix that
        isn't actually anyone's regular representation. Use this with
        the exact (``as_float=False``) matrix produced by
        :meth:`to_matrix`; validating a float matrix can spuriously fail
        or pass due to rounding.

        Examples
        --------
            >>> x = Hy(Hy(1, 2), Hy(3, 4))  # doctest: +SKIP
            >>> Hy.from_matrix(x.to_matrix()) == x  # doctest: +SKIP
            True
        """
        np = _import_optional("numpy", "numpy")
        arr = np.asarray(matrix)
        if arr.ndim != 2 or arr.shape[0] != arr.shape[1]:
            raise ValueError(
                f"from_matrix() expects a square 2D matrix, got shape {arr.shape}"
            )
        n = arr.shape[0]
        if n < 2 or (n & (n - 1)) != 0:
            raise ValueError(
                "from_matrix() expects a size that is a power of 2 and "
                f"at least 2; got {n}"
            )
        inferred_rank = n.bit_length() - 1
        if rank is not None and rank != inferred_rank:
            raise ValueError(
                f"matrix of shape {arr.shape} implies rank {inferred_rank}, "
                f"not the given rank={rank}"
            )
        rank = inferred_rank
        if rank >= 3 and not allow_nonassociative:
            raise ValueError(
                f"rank-{rank} matrices (octonions and beyond) are not the "
                "image of a true algebra homomorphism (see to_matrix()); "
                "pass allow_nonassociative=True to reconstruct anyway."
            )
        if kind not in ("left", "right"):
            raise ValueError(f"kind must be 'left' or 'right', got {kind!r}")

        result = cls.from_array(list(arr[:, 0]), signs=signs)

        if validate:
            rebuilt = result.to_matrix(kind=kind, allow_nonassociative=True)
            if not np.array_equal(rebuilt, arr):
                raise ValueError(
                    "from_matrix(): the given matrix is not the regular "
                    f"representation of any rank-{rank} Hy value (its "
                    "columns are inconsistent with column 0)"
                )
        return result

    # ---------------------------------------------------------------- #
    # LaTeX rendering
    # ---------------------------------------------------------------- #
    def latex(self, *, vinculum: str = "horizontal", mode: str = "plain") -> str:
        """Render this value as a LaTeX math expression -- handy for
        displaying a ``Hy`` in a Jupyter notebook, e.g.::

            from IPython.display import Math
            Math(some_hy.latex())

        Basis-unit labels are rendered the same way :func:`str` renders
        them (``j`` at rank 1; ``i``, ``j``, ``k`` at rank 2; ``i, j,
        k, L, iL, jL, kL`` at rank 3), except that the ``e1, e2, ...``
        labels used from rank 4 up are subscripted: ``e5`` becomes
        ``e_{5}``.

        Parameters
        ----------
        vinculum : {"horizontal", "diagonal"}, default "horizontal"
            How to typeset a non-integer coefficient's fraction bar
            (its *vinculum*). ``"horizontal"`` renders it as
            ``\\frac{num}{den}``; ``"diagonal"`` renders it as a plain
            slash, ``num/den``.
        mode : {"plain", "inline", "display"}, default "plain"
            Whether to wrap the expression in LaTeX math delimiters.
            ``"plain"`` returns the bare expression (for embedding in a
            larger expression); ``"inline"`` wraps it in ``$...$``;
            ``"display"`` wraps it in ``\\[...\\]``.

        Examples
        --------
            >>> Hy('5/2', '-16/5').latex()
            '\\\\frac{5}{2}-\\\\frac{16}{5}j'
            >>> Hy('5/2', '-16/5').latex(vinculum='diagonal')
            '5/2-16/5j'
            >>> Hy(Hy(Hy(1, 0), Hy(0, 0)), Hy(Hy(0, 1), Hy(0, 0))).latex()
            '1+iL'
        """
        if vinculum not in ("horizontal", "diagonal"):
            raise ValueError(
                f"vinculum must be 'horizontal' or 'diagonal', got {vinculum!r}"
            )
        if mode not in ("plain", "inline", "display"):
            raise ValueError(
                f"mode must be 'plain', 'inline', or 'display', got {mode!r}"
            )
        coeffs = _flatten(self)
        labels = [_latex_unit_label(lbl) for lbl in _basis_labels(self.rank)]
        body = _format_terms_latex(coeffs, labels, vinculum)
        if mode == "inline":
            return f"${body}$"
        if mode == "display":
            return f"\\[{body}\\]"
        return body


    # ---------------------------------------------------------------- #
    # Interoperability with other quaternion implementations
    #
    # Only rank-2 values (and rank-1 values, embedded as a + b*i) can be
    # converted.  The other packages are imported lazily, so none of them
    # is a dependency of hyprat.  Coordinates are always ordered
    # (w, x, y, z) == (1, i, j, k) here, matching Hy's own
    # Cayley-Dickson order; the ``from_*`` methods take care of any
    # reordering the other package needs.
    # ---------------------------------------------------------------- #
    def to_sympy(self):
        """Convert this quaternion to a ``sympy.algebras.quaternion.Quaternion``
        with exact ``sympy.Rational`` (or ``Integer``) components -- unlike
        the other two packages ``Hy`` interoperates with, SymPy represents
        rational numbers exactly, so no rounding happens in this direction.

        Rank-1 values are embedded (``a + b*i``); rank >= 3 raises
        ``ValueError``.  Requires SymPy (``pip install sympy``).

            >>> Hy(Hy(1, 2), Hy(3, 4)).to_sympy()  # doctest: +SKIP
            1 + 2*i + 3*j + 4*k
        """
        sympy = _import_optional("sympy", "sympy")
        w, x, y, z = (
            sympy.Rational(f.numerator, f.denominator)
            for f in _quaternion_fractions(self, "a sympy Quaternion")
        )
        return sympy.Quaternion(w, x, y, z)

    @classmethod
    def from_sympy(cls, q, *, exact: bool = False, max_denominator=None) -> "Hy":
        """Build a rank-2 ``Hy`` from a single
        ``sympy.algebras.quaternion.Quaternion``.

        A component that is already an exact SymPy number -- a
        ``Rational`` or an ``Integer``, which is how a ``Quaternion``
        built from ints/``Fraction``/``Rational`` naturally stores its
        components -- converts to the *same* ``Fraction`` exactly, with
        no rounding at all.  A component that is a ``Float`` instead
        converts via the same policy as :meth:`from_numpy_quaternion`
        (see *Float conversion* there).  A non-real, non-finite (e.g.
        ``nan``, ``oo``), or symbolic component raises ``ValueError`` or
        ``TypeError``.

            >>> from sympy import Quaternion, Rational  # doctest: +SKIP
            >>> Hy.from_sympy(Quaternion(1, Rational(2, 3), -4, 0))  # doctest: +SKIP
            Hy(Hy('1', '2/3'), Hy('-4', '0'))
        """
        if exact and max_denominator is not None:
            raise ValueError("exact=True and max_denominator are mutually exclusive")
        try:
            components = (q.a, q.b, q.c, q.d)
        except AttributeError:
            raise TypeError(
                "from_sympy() expects a sympy.algebras.quaternion.Quaternion, "
                f"not {type(q).__name__}"
            ) from None
        fracs = [_sympy_to_fraction(v, exact, max_denominator) for v in components]
        return _unflatten(fracs, 2)

    def to_numpy_quaternion(self):
        """Convert this quaternion to a ``numpy-quaternion`` scalar,
        ``numpy.quaternion(w, x, y, z)``, with ``float64`` coordinates.

        Rank-1 values are embedded (``a + b*i``); rank >= 3 raises
        ``ValueError``.  Requires the ``numpy-quaternion`` package
        (``pip install numpy-quaternion``; imported as ``quaternion``).

            >>> Hy(Hy(1, 2), Hy(3, 4)).to_numpy_quaternion()  # doctest: +SKIP
            quaternion(1, 2, 3, 4)
        """
        quaternion = _import_optional("quaternion", "numpy-quaternion")
        if not hasattr(quaternion, "quaternion"):
            raise ImportError(
                "found a module named 'quaternion' that is not the "
                "numpy-quaternion package (pip install numpy-quaternion)"
            )
        w, x, y, z = _quaternion_floats(self, "a numpy quaternion")
        return quaternion.quaternion(w, x, y, z)

    @classmethod
    def from_numpy_quaternion(cls, q, *, exact: bool = False,
                              max_denominator=None) -> "Hy":
        """Build a rank-2 ``Hy`` from a ``numpy-quaternion`` scalar
        (``numpy.quaternion``).  ``q`` may be any object with numeric
        ``w, x, y, z`` attributes, but not an array of quaternions.

        **Float conversion.**  By default each float becomes the rational
        with the same shortest decimal representation (``0.1`` -> ``1/10``),
        just as :meth:`from_array` does.  ``max_denominator=N`` then
        replaces it with the closest fraction having denominator <= N (so
        a float that came from ``1/3`` becomes ``1/3`` again), and
        ``exact=True`` uses the exact binary value of the float instead.
        The last two are mutually exclusive.  ``nan`` and ``inf`` raise
        ``ValueError``.

            >>> import numpy as np, quaternion  # doctest: +SKIP
            >>> Hy.from_numpy_quaternion(np.quaternion(1.5, -0.25, 3, 0.1))  # doctest: +SKIP
            Hy(Hy('3/2', '-1/4'), Hy('3', '1/10'))
        """
        if getattr(q, "ndim", 0) != 0:
            raise ValueError(
                "from_numpy_quaternion() expects a single quaternion, not an "
                f"array of shape {getattr(q, 'shape', None)}; convert them one "
                "at a time"
            )
        try:
            w, x, y, z = q.w, q.x, q.y, q.z
        except AttributeError:
            raise TypeError(
                "from_numpy_quaternion() expects a numpy.quaternion, "
                f"not {type(q).__name__}"
            ) from None
        return cls._from_wxyz((w, x, y, z), exact, max_denominator)

    def to_quaternionic(self):
        """Convert this quaternion to a ``quaternionic.array`` of shape
        ``(4,)`` holding ``[w, x, y, z]`` as ``float64``.

        Rank-1 values are embedded (``a + b*i``); rank >= 3 raises
        ``ValueError``.  Requires the ``quaternionic`` package
        (``pip install quaternionic``).

            >>> Hy(Hy(1, 2), Hy(3, 4)).to_quaternionic()  # doctest: +SKIP
            quaternionic.array([1., 2., 3., 4.])
        """
        quaternionic = _import_optional("quaternionic", "quaternionic")
        return quaternionic.array(list(_quaternion_floats(self, "a quaternionic array")))

    @classmethod
    def from_quaternionic(cls, q, *, exact: bool = False,
                          max_denominator=None) -> "Hy":
        """Build a rank-2 ``Hy`` from a single ``quaternionic`` quaternion:
        a ``quaternionic.array`` (or any array-like of 4 real numbers)
        of shape ``(4,)``, ordered ``[w, x, y, z]``.  An array holding
        several quaternions (shape ``(..., 4)``) raises ``ValueError``.

        See :meth:`from_numpy_quaternion` for the meaning of ``exact`` and
        ``max_denominator`` ("Float conversion" in its docstring).

            >>> import quaternionic  # doctest: +SKIP
            >>> Hy.from_quaternionic(quaternionic.array([1, 0.5, 0, -2]))  # doctest: +SKIP
            Hy(Hy('1', '1/2'), Hy('0', '-2'))
        """
        data = getattr(q, "ndarray", q)        # unwrap a quaternionic array
        try:
            if isinstance(data, (str, bytes)):
                raise TypeError
            items = list(data)
        except TypeError:
            raise TypeError(
                "from_quaternionic() expects a quaternionic array or a "
                f"sequence of 4 numbers, not {type(q).__name__}"
            ) from None
        shape = getattr(data, "shape", (len(items),))
        if shape != (4,):
            raise ValueError(
                "from_quaternionic() expects a single quaternion, i.e. an "
                f"array of shape (4,), not shape {tuple(shape)}; convert "
                "multiple quaternions one at a time"
            )
        if not all(isinstance(v, numbers.Real) for v in items):
            raise TypeError(
                "from_quaternionic() expects 4 real numbers [w, x, y, z]"
            )
        return cls._from_wxyz(items, exact, max_denominator)

    @classmethod
    def _from_wxyz(cls, wxyz, exact, max_denominator) -> "Hy":
        """Internal: build a rank-2 Hy from four floats ordered (w, x, y, z)."""
        if exact and max_denominator is not None:
            raise ValueError("exact=True and max_denominator are mutually exclusive")
        fracs = [_float_to_fraction(v, exact, max_denominator) for v in wxyz]
        return _unflatten(fracs, 2)


# ============================================================================
# Module-level recursive algebra (Cayley-Dickson construction)
#
# These functions all operate on "raw" values, each of which is either a
# Fraction (rank 0) or a Hy (rank >= 1), and freely mix ranks by promoting
# ("embedding") the lower-rank operand up to match the higher-rank one.
# ============================================================================

def _rank(x) -> int:
    return 0 if isinstance(x, Fraction) else x.rank


# --------------------------------------------------------------------------
# Signatures (the per-level Cayley-Dickson parameters)
# --------------------------------------------------------------------------

def _coerce_mu(m) -> Fraction:
    """A single Cayley-Dickson parameter: any nonzero rational."""
    if isinstance(m, bool):
        raise TypeError(f"mu must be a nonzero rational, not {m!r}")
    if isinstance(m, Fraction):
        mu = m
    elif isinstance(m, int):
        mu = Fraction(m)
    elif isinstance(m, float):
        mu = Fraction(str(m))
    elif isinstance(m, str):
        try:
            mu = Fraction(m)
        except ValueError as e:
            raise ValueError(f"cannot parse {m!r} as a rational mu") from e
    else:
        raise TypeError(
            f"cannot use {m!r} (type {type(m).__name__}) as a mu; "
            "expected a nonzero int, Fraction, float or fraction string"
        )
    if mu == 0:
        raise ValueError("mu must be nonzero")
    return mu


def _preset_key(name: str) -> str:
    return name.strip().lower().replace("_", "-").replace(" ", "-")


def _coerce_signs(signs) -> tuple:
    """A signature given as a preset name or an iterable of mus, as a
    tuple of nonzero Fractions (lowest doubling level first)."""
    if isinstance(signs, str):
        key = _preset_key(signs)
        if key not in _SIGN_PRESETS:
            raise ValueError(
                f"unknown signature preset {signs!r}; known presets: "
                + ", ".join(sorted(_SIGN_PRESETS))
            )
        signs = _SIGN_PRESETS[key]
    try:
        return tuple(_coerce_mu(m) for m in signs)
    except TypeError:
        if isinstance(signs, (int, float, Fraction)):
            raise TypeError(
                "signs must be a preset name or a sequence of mu values, "
                f"not {signs!r}"
            ) from None
        raise


def _resolve_rank_and_signs(rank, signs, *, min_rank: int):
    """Validate a (rank, signs) pair: either may be omitted (but not
    both); ``signs`` defaults to all -1. Returns ``(rank, signs_tuple)``."""
    sg = None if signs is None else _coerce_signs(signs)
    if rank is None:
        if sg is None:
            kind = "positive" if min_rank >= 1 else "non-negative"
            raise ValueError(f"rank must be a {kind} int, got {rank!r}")
        rank = len(sg)
    if not isinstance(rank, int) or isinstance(rank, bool) or rank < min_rank:
        kind = "positive" if min_rank >= 1 else "non-negative"
        raise ValueError(f"rank must be a {kind} int, got {rank!r}")
    if sg is None:
        sg = (_DEFAULT_MU,) * rank
    elif len(sg) != rank:
        raise ValueError(
            f"signs has length {len(sg)} but the rank is {rank}; "
            "one mu per doubling level is required"
        )
    return rank, sg


def _signs(x) -> tuple:
    """The signature of a raw value: () for a Fraction, otherwise one mu
    per doubling level, lowest first (so the last entry is x's own mu)."""
    if isinstance(x, Fraction):
        return ()
    return _signs(x._real) + (x._mu,)


def _common_signs(x, y) -> tuple:
    """The signature of the algebra in which a binary operation on `x`
    and `y` takes place: the longer of their two signatures, provided the
    shorter is a prefix of it (so the lower-rank operand embeds)."""
    sx, sy = _signs(x), _signs(y)
    big, small = (sx, sy) if len(sx) >= len(sy) else (sy, sx)
    if big[: len(small)] != small:
        raise ValueError(
            "incompatible signatures: cannot combine a value with signs "
            f"{tuple(map(str, sx))} and one with signs {tuple(map(str, sy))}"
        )
    return big


def _embed_signs(x, target: tuple):
    """Promote `x` into the algebra with signature `target`, by pairing it
    with zeros at each additional doubling step. `x`'s own signature must
    be a prefix of `target`, and its rank no larger."""
    r, n = _rank(x), len(target)
    if r > n:
        raise ValueError(f"cannot embed a rank-{r} value into rank {n}")
    if r == n:
        if _signs(x) != target:
            raise ValueError(
                f"signature mismatch: value has signs "
                f"{tuple(map(str, _signs(x)))}, expected "
                f"{tuple(map(str, target))}"
            )
        return x
    return Hy._make(
        _embed_signs(x, target[:-1]),
        _embed_signs(Fraction(0), target[:-1]),
        target[-1],
    )


def _embed(x, target):
    """Promote `x` (Fraction or Hy) up to exactly `target`, by pairing it
    with zeros at each doubling step. `target` is either a rank (an int;
    any *new* doubling levels get the default mu = -1) or a full signature
    tuple. Raises if `x` is already of *higher* rank than `target` (you
    cannot un-embed)."""
    if isinstance(target, int):
        r = _rank(x)
        if r == target:
            return x
        if r > target:
            raise ValueError(
                f"cannot embed a rank-{r} value into rank {target}"
            )
        target = _signs(x) + (_DEFAULT_MU,) * (target - r)
    return _embed_signs(x, target)


# --------------------------------------------------------------------------
# The algebra itself.  The public functions (add, mul, ...) accept raw
# values of any (compatible) ranks and embed them into a common algebra;
# the underscore versions then recurse on operands already known to have
# identical rank and signature, so no further checking is needed.
# --------------------------------------------------------------------------

def add(x, y):
    T = _common_signs(x, y)
    return _add(_embed_signs(x, T), _embed_signs(y, T))


def _add(x, y):
    if isinstance(x, Fraction):
        return x + y
    return Hy._make(_add(x._real, y._real), _add(x._imag, y._imag), x._mu)


def neg(x):
    if isinstance(x, Fraction):
        return -x
    return Hy._make(neg(x._real), neg(x._imag), x._mu)


def sub(x, y):
    return add(x, neg(y))


def _sub(x, y):
    return _add(x, neg(y))


def conj(x):
    """Cayley-Dickson conjugate: conj(a, b) = (conj(a), -b). It does not
    depend on mu."""
    if isinstance(x, Fraction):
        return x
    return Hy._make(conj(x._real), neg(x._imag), x._mu)


def _scale(x, f: Fraction):
    """x * f for a rational scalar f, componentwise."""
    if isinstance(x, Fraction):
        return x * f
    return Hy._make(_scale(x._real, f), _scale(x._imag, f), x._mu)


def mul(x, y):
    """Cayley-Dickson product with parameter mu at each doubling level:
    (a,b)(c,d) = (ac + mu*conj(d)b, da + b*conj(c)). With the default
    mu = -1 this is the classical (ac - conj(d)b, da + b*conj(c))."""
    T = _common_signs(x, y)
    if not T:
        return x * y
    return _mul(_embed_signs(x, T), _embed_signs(y, T))


def _mul(x, y):
    if isinstance(x, Fraction):
        return x * y
    a, b = x._real, x._imag
    c, d = y._real, y._imag
    mu = x._mu
    cross = _mul(conj(d), b)
    if mu == -1:
        real_part = _sub(_mul(a, c), cross)
    elif mu == 1:
        real_part = _add(_mul(a, c), cross)
    else:
        real_part = _add(_mul(a, c), _scale(cross, mu))
    imag_part = _add(_mul(d, a), _mul(b, conj(c)))
    return Hy._make(real_part, imag_part, mu)


def abs2(x) -> Fraction:
    """The quadratic form x*conj(x), as an exact Fraction: N(a, b) = N(a)
    - mu*N(b), N(r) = r*r. For the classical algebras (every mu = -1)
    this is the squared Euclidean norm, the sum of squares of every real
    coordinate; otherwise it is indefinite."""
    if isinstance(x, Fraction):
        return x * x
    return abs2(x._real) - x._mu * abs2(x._imag)


def _scalar_div(x, f: Fraction):
    if isinstance(x, Fraction):
        return x / f
    return Hy._make(_scalar_div(x._real, f), _scalar_div(x._imag, f), x._mu)


def inverse(x):
    n = abs2(x)
    if n == 0:
        if _is_zero_val(x):
            raise ZeroDivisionError(
                "hypercomplex value has zero norm; not invertible"
            )
        raise ZeroDivisionError(
            "hypercomplex value has zero norm (it is a nonzero null "
            "element of a split algebra); not invertible"
        )
    return _scalar_div(conj(x), n)


def div(x, y):
    return mul(x, inverse(y))


def _is_zero_val(x) -> bool:
    if isinstance(x, Fraction):
        return x == 0
    return _is_zero_val(x._real) and _is_zero_val(x._imag)


def _values_equal(x, y) -> bool:
    """Equality of two raw values: the same element once zero-padded to a
    common rank, *and* living in the same algebra (same mu at every level
    that is actually used)."""
    return _to_nested_tuple(_canonical_trim(x)) == _to_nested_tuple(
        _canonical_trim(y)
    )


def _canonical_trim(x):
    """Strip away outer (Fraction-zero-imag) layers, for hashing purposes."""
    while isinstance(x, Hy) and _is_zero_val(x._imag):
        x = x._real
    return x


def _to_nested_tuple(x):
    if isinstance(x, Fraction):
        return x
    return (_to_nested_tuple(x._real), _to_nested_tuple(x._imag), x._mu)


def _flatten(x) -> list:
    if isinstance(x, Fraction):
        return [x]
    return _flatten(x._real) + _flatten(x._imag)


def _unflatten(coeffs, rank, signs=None):
    if signs is None:
        signs = (_DEFAULT_MU,) * rank
    if rank == 0:
        return coeffs[0]
    half = len(coeffs) // 2
    return Hy._make(
        _unflatten(coeffs[:half], rank - 1, signs[:-1]),
        _unflatten(coeffs[half:], rank - 1, signs[:-1]),
        signs[-1],
    )


# ============================================================================
# Scalar coercion (Fraction/int/float/str/Hy -> a component)
# ============================================================================

def _coerce_component(v):
    if isinstance(v, Hy):
        return v
    if isinstance(v, Fraction):
        return v
    if isinstance(v, bool):  # bool is a subclass of int; treat explicitly first
        return Fraction(int(v))
    if isinstance(v, int):
        return Fraction(v)
    if isinstance(v, float):
        # Fraction(str(v)) reproduces the "obvious" decimal value (e.g. 3.2
        # -> 16/5) rather than the exact (ugly) binary value Fraction(v)
        # would give.
        return Fraction(str(v))
    if isinstance(v, complex):
        return Hy(Fraction(str(v.real)), Fraction(str(v.imag)))
    if isinstance(v, str):
        try:
            return Fraction(v)
        except ValueError:
            return _parse(v)
    raise TypeError(f"cannot use {v!r} (type {type(v).__name__}) as a Hy component")


def _coerce_flat_element(v) -> Fraction:
    """Coerce a single ``Hy.from_array()`` element to a plain Fraction.

    This is deliberately narrower than ``_coerce_component``: a string
    element here must be a *plain* fraction/decimal like ``'5/2'`` or
    ``'3.2'``, not a composite expression like ``'1+2j'`` -- from_array()
    always supplies coefficients one real coordinate at a time.
    """
    if isinstance(v, Fraction):
        return v
    if isinstance(v, bool):  # bool is a subclass of int; handle it first
        return Fraction(int(v))
    if isinstance(v, int):
        return Fraction(v)
    if isinstance(v, float):
        return Fraction(str(v))
    if isinstance(v, str):
        try:
            return Fraction(v)
        except ValueError as e:
            raise ValueError(
                f"cannot parse {v!r} as a plain fraction for from_array()"
            ) from e
    raise TypeError(
        f"cannot use {v!r} (type {type(v).__name__}) as a from_array() element"
    )


# ============================================================================
# Helper for to_matrix()/from_matrix()
# ============================================================================

def _positive_units(rank: int, signs=None) -> list:
    """The 2**rank positive units of the rank-`rank` algebra (1, then each
    imaginary unit), in the same order as `components()`/`to_array()`.
    Just Hy.units(rank, signs=signs) filtered down to its non-negated
    entries."""
    return [u for name, u in Hy.units(rank, signs=signs).items()
            if not name.startswith("-")]


# ============================================================================
# Helpers for interoperability with other quaternion packages
# ============================================================================

def _import_optional(module_name: str, pip_name: str):
    """Import an optional third-party module, or raise a helpful ImportError."""
    try:
        return importlib.import_module(module_name)
    except ImportError as e:
        raise ImportError(
            f"this conversion requires the optional package {pip_name!r} "
            f"(pip install {pip_name})"
        ) from e


def _quaternion_fractions(h: "Hy", target: str) -> tuple:
    """The four coordinates (w, x, y, z) of ``h`` as exact Fractions.

    Rank-1 values are embedded as a + b*i; rank >= 3 cannot be
    represented as a quaternion at all.
    """
    if h.rank > 2:
        raise ValueError(
            f"cannot convert a rank-{h.rank} Hy ({h.dimension} coordinates) to "
            f"{target}; only rank 1 (embedded as a + b*i) and rank 2 "
            "(quaternions) can be converted"
        )
    if h.signs != (_DEFAULT_MU,) * h.rank:
        raise ValueError(
            f"cannot convert a Hy with signs {tuple(map(str, h.signs))} to "
            f"{target}; only ordinary complex/quaternion values (every "
            "mu equal to -1) can be converted"
        )
    return tuple(_flatten(_embed(h, 2)))


def _quaternion_floats(h: "Hy", target: str) -> tuple:
    """The four coordinates (w, x, y, z) of ``h`` as Python floats."""
    return tuple(float(c) for c in _quaternion_fractions(h, target))


def _sympy_to_fraction(v, exact: bool, max_denominator) -> Fraction:
    """Convert a single sympy quaternion component to a Fraction.

    An exact sympy number (Rational/Integer) converts exactly, with no
    rounding; a Float goes through the same policy as _float_to_fraction.
    A non-real, non-finite, or symbolic value raises ValueError/TypeError.
    """
    if not getattr(v, "is_number", False):
        raise TypeError(
            f"cannot use {v!r} (type {type(v).__name__}) as a sympy quaternion "
            "component; it must be a concrete number, not a symbolic expression"
        )
    if v.is_real is not True or v.is_finite is not True:
        raise ValueError(
            f"cannot convert the non-real or non-finite sympy value {v!r} to a Fraction"
        )
    if v.is_rational:
        return Fraction(int(v.p), int(v.q))
    return _float_to_fraction(float(v), exact, max_denominator)


def _float_to_fraction(v, exact: bool, max_denominator) -> Fraction:
    """Convert a float-like ``v`` to a Fraction (see Hy.from_numpy_quaternion)."""
    try:
        if isinstance(v, (str, bytes)):
            raise TypeError
        f = float(v)
    except (TypeError, ValueError):
        raise TypeError(
            f"cannot use {v!r} (type {type(v).__name__}) as a real coordinate"
        ) from None
    if not math.isfinite(f):
        raise ValueError(f"cannot convert the non-finite value {f!r} to a Fraction")
    if exact:
        return Fraction(f)
    frac = Fraction(str(f))     # shortest decimal repr, as Hy.from_array does
    if max_denominator is not None:
        frac = frac.limit_denominator(max_denominator)
    return frac


# ============================================================================
# String formatting
# ============================================================================

# _OCTONION_LABELS = ["", "i", "j", "k", "L", "iL", "jL", "kL"]


def _basis_labels(rank: int):
    n = 2 ** rank
    if rank <= 0:
        return [""]
    if rank == 1:
        return ["", "j"]
    if rank == 2:
        return ["", "i", "j", "k"]
    if rank == 3:
        # return list(_OCTONION_LABELS)
        return ["", "i", "j", "k", "L", "iL", "jL", "kL"]
    return [""] + [f"e{i}" for i in range(1, n)]


def _format_terms(coeffs, labels) -> str:
    parts = []
    for c, u in zip(coeffs, labels):
        if c == 0:
            continue
        negative = c < 0
        mag = -c if negative else c
        if u == "":
            body = str(mag)
        elif mag == 1:
            body = u
        else:
            body = f"{mag}{u}"
        if not parts:
            parts.append(f"-{body}" if negative else body)
        else:
            parts.append(f"{'-' if negative else '+'}{body}")
    return "".join(parts) if parts else "0"


# ============================================================================
# LaTeX formatting (used by Hy.latex())
# ============================================================================

def _latex_unit_label(lbl: str) -> str:
    """'' / 'i' / 'j' / 'k' pass through unchanged; 'e12' -> 'e_{12}'."""
    if lbl.startswith("e") and lbl[1:].isdigit():
        return f"e_{{{lbl[1:]}}}"
    return lbl


def _format_fraction_latex(c: Fraction, vinculum: str) -> str:
    if c.denominator == 1:
        return str(c.numerator)
    if vinculum == "diagonal":
        return f"{c.numerator}/{c.denominator}"
    return f"\\frac{{{c.numerator}}}{{{c.denominator}}}"


def _format_terms_latex(coeffs, labels, vinculum: str) -> str:
    parts = []
    for c, u in zip(coeffs, labels):
        if c == 0:
            continue
        negative = c < 0
        mag = -c if negative else c
        mag_str = _format_fraction_latex(mag, vinculum)
        if u == "":
            body = mag_str
        elif mag == 1:
            body = u
        else:
            body = f"{mag_str}{u}"
        if not parts:
            parts.append(f"-{body}" if negative else body)
        else:
            parts.append(f"{'-' if negative else '+'}{body}")
    return "".join(parts) if parts else "0"


# ============================================================================
# Parsing
# ============================================================================

_TERM_RE = re.compile(
    r"^([+-]?)(\d+/\d+|\d+\.\d+|\d+)?(iL|jL|kL|i|j|k|L|e\d+)?$"
)


def _parse(s: str, signs=None, pad: bool = False) -> Hy:
    s = s.strip()
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1]
    s = s.replace(" ", "")
    if s == "":
        raise ValueError("cannot parse an empty hypercomplex expression")

    tokens = re.findall(r"[+-]?[^+-]+", s)
    coeff_map: dict = {}
    max_e_index = 0
    has_i = has_k = False
    has_octonion_only = False  # 'L', 'iL', 'jL', 'kL' -- unique to rank 3

    for tok in tokens:
        m = _TERM_RE.match(tok)
        if not m:
            raise ValueError(f"cannot parse term {tok!r} in {s!r}")
        sign, num, unit = m.groups()
        if num is None and unit is None:
            raise ValueError(f"empty term in {s!r}")
        value = Fraction(num) if num is not None else Fraction(1)
        if sign == "-":
            value = -value
        label = unit or ""
        if label == "i":
            has_i = True
        elif label == "k":
            has_k = True
        elif label in ("L", "iL", "jL", "kL"):
            has_octonion_only = True
        elif label.startswith("e"):
            max_e_index = max(max_e_index, int(label[1:]))
        coeff_map[label] = coeff_map.get(label, Fraction(0)) + value

    if max_e_index > 0:
        # 'e1, e2, ...' labels are only used from rank 4 up (octonions,
        # rank 3, use the 'i/j/k/L/iL/jL/kL' labels instead).
        n = 16
        while n - 1 < max_e_index:
            n *= 2
        rank = n.bit_length() - 1
        labels = _basis_labels(rank)
    elif has_octonion_only:
        rank = 3
        labels = _basis_labels(rank)
    elif has_i or has_k:
        rank = 2
        labels = _basis_labels(rank)
    else:
        rank = 1
        labels = _basis_labels(rank)

    known = set(labels)
    unknown = set(coeff_map) - known
    if unknown:
        raise ValueError(
            f"unit(s) {sorted(unknown)} inconsistent with the rest of {s!r}"
        )

    if signs is None:
        sg = (_DEFAULT_MU,) * rank
    else:
        sg = _coerce_signs(signs)
        if pad:
            # Internal use (operand coercion): read the text in the algebra
            # of the other operand, truncating or extending with the
            # default mu as the text's own rank requires.
            sg = sg[:rank] + (_DEFAULT_MU,) * (rank - len(sg[:rank]))
        elif len(sg) != rank:
            raise ValueError(
                f"signs has length {len(sg)}, but {s!r} has rank {rank}; "
                "one mu per doubling level is required"
            )

    coeffs = [coeff_map.get(lbl, Fraction(0)) for lbl in labels]
    return _unflatten(coeffs, rank, sg)


if __name__ == "__main__":
    # A handful of sanity checks / usage examples.
    z = Hy("5/2", "-16/5")
    print("z       =", z, "  repr:", repr(z))
    assert str(z) == "(5/2-16/5j)"
    assert z == complex(2.5, -3.2)

    q = Hy(Hy(1, 2), Hy(3, 4))
    print("q       =", q)
    assert str(q) == "(1+2i+3j+4k)"

    i = Hy(Hy(0, 1), Hy(0, 0))
    j = Hy(Hy(0, 0), Hy(1, 0))
    k = Hy(Hy(0, 0), Hy(0, 1))
    assert i * i == -1 and j * j == -1 and k * k == -1
    assert i * j == k and j * i == -k
    print("Quaternion units check out: i*j =", i * j, " j*i =", j * i)

    o = Hy(q, Hy(Hy(0, 1), Hy(0, 0)))
    print("octonion o =", o, " rank", o.rank, " dim", o.dimension)

    # multiplicativity of the norm, |xy| = |x||y|, for quaternions:
    x = Hy(Hy("1/2", "1/3"), Hy("-2/5", "7"))
    y = Hy(Hy("3", "-1/7"), Hy("5/2", "1/4"))
    assert abs2(mul(x, y)) == abs2(x) * abs2(y)
    print("Quaternion norm is multiplicative: OK")

    # division / inverse
    assert x * x.inverse() == Hy(1)
    print("x * x.inverse() == 1: OK")

    # round trip through str/parse
    for val in (z, q, o):
        s = str(val)
        parsed = Hy.parse(s)
        assert parsed == val, (val, s, parsed)
    print("str/parse round trips: OK")

    # hash consistency
    assert hash(Hy("3")) == hash(Hy(Hy("3", "0"), Hy("0", "0")))
    print("hash consistency across ranks: OK")

    # random(), with a reproducible seed
    Hy.seed(42)
    r1 = Hy.random(2)
    Hy.seed(42)
    r2 = Hy.random(2)
    assert r1 == r2 and r1.rank == 2
    print("Hy.seed()/Hy.random() reproducibility: OK, e.g. random quaternion =", r1)

    # from_array() / to_array(), including a mix of numbers and strings
    arr = [1, "2/3", 3.5, "-1/4"]
    h_from_arr = Hy.from_array(arr)
    assert h_from_arr.rank == 2
    assert h_from_arr.to_array() == [Fraction(1), Fraction(2, 3), Fraction(7, 2), Fraction(-1, 4)]
    assert h_from_arr.to_array(as_str=True) == ["1", "2/3", "7/2", "-1/4"]
    assert Hy.from_array(h_from_arr.to_array(as_str=True)) == h_from_arr
    print("Hy.from_array()/to_array() round trip, mixed input types: OK")

    # Hy.units() / is_unit()
    units1 = Hy.units(1)
    assert units1 == {"1": Hy(1, 0), "-1": Hy(-1, 0), "j": Hy(0, 1), "-j": Hy(0, -1)}
    assert all(u.is_unit() for u in units1.values())
    assert not Hy(1, 1).is_unit()
    assert not Hy(0, 0).is_unit()
    assert Hy.units(0) == {"1": Fraction(1), "-1": Fraction(-1)}
    print("Hy.units()/is_unit(): OK, e.g. Hy.units(1) =", units1)

    # latex()
    assert Hy("5/2", "-16/5").latex() == r"\frac{5}{2}-\frac{16}{5}j"
    assert Hy("5/2", "-16/5").latex(vinculum="diagonal") == "5/2-16/5j"
    oct_e = Hy(Hy(Hy(1, 0), Hy(0, 0)), Hy(Hy(0, 1), Hy(0, 0)))
    assert oct_e.latex() == "1+iL"
    assert Hy("1/2").latex(mode="inline") == r"$\frac{1}{2}$"
    print("Hy.latex(): OK, e.g. Hy('5/2', '-16/5').latex() =", Hy("5/2", "-16/5").latex())

    # to_matrix() / from_matrix()
    try:
        import numpy as _np

        assert Hy.from_matrix(q.to_matrix()) == q
        assert _np.array_equal(
            _np.dot(x.to_matrix(), y.to_matrix()), (x * y).to_matrix()
        )
        try:
            o.to_matrix()
            raise AssertionError("expected ValueError for rank >= 3 to_matrix()")
        except ValueError:
            pass
        print("Hy.to_matrix()/Hy.from_matrix(): OK (roundtrip + homomorphism "
              "at rank 2, rank>=3 guard raises as expected)")
    except ImportError:
        print("Hy.to_matrix()/Hy.from_matrix(): skipped (numpy not installed)")

    # signatures: split-complex and split-quaternions
    sj = Hy(0, 1, mu=1)
    assert sj * sj == 1 and (Hy(1, 1, mu=1) * Hy(1, -1, mu=1)).is_zero()
    su = Hy.units(signs="split-quaternion")
    assert su["i"] * su["i"] == -1 and su["j"] * su["j"] == 1
    assert su["i"] * su["j"] == su["k"] and su["k"] * su["k"] == 1
    sx, sy = Hy.random(signs="split-quaternion", seed=1), Hy.random(signs="split-quaternion", seed=2)
    assert abs2(mul(sx, sy)) == abs2(sx) * abs2(sy)
    assert eval(repr(sx)) == sx
    print("split algebras / signatures: OK")

    print("\nAll self-tests passed.")
