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
"""

from __future__ import annotations

import math
import random
import re
from fractions import Fraction
# from numbers import Number

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

    __slots__ = ("_real", "_imag")

    # ---------------------------------------------------------------- #
    # Construction
    # ---------------------------------------------------------------- #
    def __init__(self, real, imag=_MISSING):
        real_c = _coerce_component(real)

        if imag is _MISSING:
            if isinstance(real_c, Hy):
                # Hy(some_hy) is a copy/identity construction: it takes
                # on the value of `some_hy` as-is, at whatever rank that
                # already is (rather than promoting it one rank higher).
                object.__setattr__(self, "_real", real_c._real)
                object.__setattr__(self, "_imag", real_c._imag)
                return
            imag_c = Fraction(0)
        else:
            imag_c = _coerce_component(imag)

        rank = max(_rank(real_c), _rank(imag_c))
        object.__setattr__(self, "_real", _embed(real_c, rank))
        object.__setattr__(self, "_imag", _embed(imag_c, rank))

    @classmethod
    def _make(cls, real, imag) -> "Hy":
        """Internal fast constructor: assumes `real`/`imag` already have
        matching rank and skips all normalization. Not for public use."""
        obj = object.__new__(cls)
        object.__setattr__(obj, "_real", real)
        object.__setattr__(obj, "_imag", imag)
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
        """The *squared* Euclidean norm, computed exactly as a Fraction."""
        return abs2(self)

    def is_zero(self) -> bool:
        return _is_zero_val(self)

    def is_unit(self) -> bool:
        """True iff this value is one of the *units* of its own rank --
        i.e. exactly one of its ``2**rank`` real coordinates is ``+-1``
        and every other coordinate is 0. Equivalent to (but cheaper
        than) checking membership in ``Hy.units(self.rank).values()``.

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
            return Hy(Fraction(str(other.real)), Fraction(str(other.imag)))
        if isinstance(other, _ScalarLike):
            try:
                return Hy(other)
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
        return math.sqrt(float(abs2(self)))

    def __pow__(self, n):
        if not isinstance(n, int):
            return NotImplemented
        if n == 0:
            return _embed(Fraction(1), self.rank)
        base = self if n > 0 else self.inverse()
        result = _embed(Fraction(1), base.rank)
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

        return f"Hy({comp_repr(self._real)}, {comp_repr(self._imag)})"

    # ---------------------------------------------------------------- #
    # Parsing
    # ---------------------------------------------------------------- #
    @classmethod
    def parse(cls, s: str) -> "Hy":
        """Parse the kind of string produced by ``str(some_hy)`` -- e.g.
        ``'(5/2-16/5j)'``, ``'1+2i+3j+4k'``, ``'1-k+2kL'`` -- into a Hy.
        """
        return _parse(s)

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
        rank: int,
        *,
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

        Raises ``ValueError`` if ``rank`` isn't a positive int, or if both
        ``rng`` and ``seed`` are given.
        """
        if not isinstance(rank, int) or isinstance(rank, bool) or rank < 1:
            raise ValueError(f"rank must be a positive int, got {rank!r}")
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
            return Hy._make(build(r - 1), build(r - 1))

        return build(rank)

    # ---------------------------------------------------------------- #
    # Flat-array conversion
    # ---------------------------------------------------------------- #
    @classmethod
    def from_array(cls, coeffs) -> "Hy":
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
        """
        values = list(coeffs)
        n = len(values)
        if n < 2 or (n & (n - 1)) != 0:
            raise ValueError(
                "from_array() requires a length that is a power of 2 "
                f"and at least 2; got {n}"
            )
        rank = n.bit_length() - 1
        fracs = [_coerce_flat_element(v) for v in values]
        return _unflatten(fracs, rank)

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
    def units(cls, rank: int) -> dict:
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

        See also :meth:`is_unit`. Raises ``ValueError`` if ``rank``
        isn't a non-negative int.
        """
        if not isinstance(rank, int) or isinstance(rank, bool) or rank < 0:
            raise ValueError(f"rank must be a non-negative int, got {rank!r}")
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
            result[pos_key] = _unflatten(coeffs, rank)
            coeffs[idx] = Fraction(-1)
            result[neg_key] = _unflatten(coeffs, rank)
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


# ============================================================================
# Module-level recursive algebra (Cayley-Dickson construction)
#
# These functions all operate on "raw" values, each of which is either a
# Fraction (rank 0) or a Hy (rank >= 1), and freely mix ranks by promoting
# ("embedding") the lower-rank operand up to match the higher-rank one.
# ============================================================================

def _rank(x) -> int:
    return 0 if isinstance(x, Fraction) else x.rank


def _embed(x, target_rank: int):
    """Promote `x` (Fraction or Hy) up to exactly `target_rank`, by pairing
    it with zeros at each doubling step. Raises if `x` is already of
    *higher* rank than `target_rank` (you cannot un-embed)."""
    r = _rank(x)
    if r == target_rank:
        return x
    if r > target_rank:
        raise ValueError(
            f"cannot embed a rank-{r} value into rank {target_rank}"
        )
    return Hy._make(_embed(x, target_rank - 1), _embed(Fraction(0), target_rank - 1))


def add(x, y):
    rx, ry = _rank(x), _rank(y)
    r = max(rx, ry)
    xa, ya = _embed(x, r), _embed(y, r)
    if r == 0:
        return xa + ya
    return Hy._make(add(xa._real, ya._real), add(xa._imag, ya._imag))


def neg(x):
    if isinstance(x, Fraction):
        return -x
    return Hy._make(neg(x._real), neg(x._imag))


def sub(x, y):
    return add(x, neg(y))


def conj(x):
    """Cayley-Dickson conjugate: conj(a, b) = (conj(a), -b)."""
    if isinstance(x, Fraction):
        return x
    return Hy._make(conj(x._real), neg(x._imag))


def mul(x, y):
    """Cayley-Dickson product: (a,b)(c,d) = (ac - conj(d)b, da + b*conj(c))."""
    rx, ry = _rank(x), _rank(y)
    r = max(rx, ry)
    if r == 0:
        return x * y
    xa, ya = _embed(x, r), _embed(y, r)
    a, b = xa._real, xa._imag
    c, d = ya._real, ya._imag
    real_part = sub(mul(a, c), mul(conj(d), b))
    imag_part = add(mul(d, a), mul(b, conj(c)))
    return Hy._make(real_part, imag_part)


def abs2(x) -> Fraction:
    """Squared Euclidean norm: sum of squares of every real coordinate."""
    if isinstance(x, Fraction):
        return x * x
    return abs2(x._real) + abs2(x._imag)


def _scalar_div(x, f: Fraction):
    if isinstance(x, Fraction):
        return x / f
    return Hy._make(_scalar_div(x._real, f), _scalar_div(x._imag, f))


def inverse(x):
    n = abs2(x)
    if n == 0:
        raise ZeroDivisionError("hypercomplex value has zero norm; not invertible")
    return _scalar_div(conj(x), n)


def div(x, y):
    return mul(x, inverse(y))


def _is_zero_val(x) -> bool:
    if isinstance(x, Fraction):
        return x == 0
    return _is_zero_val(x._real) and _is_zero_val(x._imag)


def _values_equal(x, y) -> bool:
    r = max(_rank(x), _rank(y))
    xa, ya = _embed(x, r), _embed(y, r)
    if r == 0:
        return xa == ya
    return _values_equal(xa._real, ya._real) and _values_equal(xa._imag, ya._imag)


def _canonical_trim(x):
    """Strip away outer (Fraction-zero-imag) layers, for hashing purposes."""
    while isinstance(x, Hy) and _is_zero_val(x._imag):
        x = x._real
    return x


def _to_nested_tuple(x):
    if isinstance(x, Fraction):
        return x
    return (_to_nested_tuple(x._real), _to_nested_tuple(x._imag))


def _flatten(x) -> list:
    if isinstance(x, Fraction):
        return [x]
    return _flatten(x._real) + _flatten(x._imag)


def _unflatten(coeffs, rank):
    if rank == 0:
        return coeffs[0]
    half = len(coeffs) // 2
    return Hy._make(
        _unflatten(coeffs[:half], rank - 1), _unflatten(coeffs[half:], rank - 1)
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
# String formatting
# ============================================================================

_OCTONION_LABELS = ["", "i", "j", "k", "L", "iL", "jL", "kL"]


def _basis_labels(rank: int):
    n = 2 ** rank
    if rank <= 0:
        return [""]
    if rank == 1:
        return ["", "j"]
    if rank == 2:
        return ["", "i", "j", "k"]
    if rank == 3:
        return list(_OCTONION_LABELS)
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


def _parse(s: str) -> Hy:
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

    coeffs = [coeff_map.get(lbl, Fraction(0)) for lbl in labels]
    return _unflatten(coeffs, rank)


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

    print("\nAll self-tests passed.")
