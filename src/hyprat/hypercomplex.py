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

    >>> from hypercomplex import Hy
    >>> z = Hy('5/2', '-16/5')          # a rational complex number
    >>> str(z)
    '(5/2-16/5j)'
    >>> q = Hy(Hy(1, 0), Hy(0, 1))      # the quaternion j
    >>> str(q)
    '(1j)'
    >>> str(q * q)
    '(-1)'

Because ``Hy`` always stores exactly two components, a "plain real
number" is represented as a rank-1 ``Hy`` with a zero imaginary part
(exactly the way Python's own ``complex`` behaves: ``str(complex(3))``
is ``'(3+0j)'``).

String forms
------------

``str()`` renders a value the way Python renders ``complex`` numbers
for rank 1 (using ``j``), and the customary ``a+bi+cj+dk`` notation for
rank 2, and ``e1 .. e7`` (etc.) imaginary units for rank >= 3::

    str(Hy('5/2', '16/5'))                 -> '(5/2+16/5j)'
    str(Hy(Hy(1, 2), Hy(3, 4)))            -> '(1+2i+3j+4k)'
    str(Hy(Hy(Hy(1,0),Hy(0,0)), Hy(Hy(0,1),Hy(0,0))))
                                            -> '(1+e2)'   (etc.)

``Hy.parse(s)`` is the inverse of ``str`` (and is also used
automatically by the constructor -- see below).

``repr()`` returns Python source that reconstructs an equal value,
e.g. ``Hy('5/2', '-16/5')``.
"""

from __future__ import annotations

import math
import re
from fractions import Fraction
from numbers import Number

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
        ``'(5/2-16/5j)'``, ``'1+2i+3j+4k'``, ``'1-e3+2e7'`` -- into a Hy.
        """
        return _parse(s)

    from_string = parse  # convenient alias


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


# ============================================================================
# String formatting
# ============================================================================

def _basis_labels(rank: int):
    n = 2 ** rank
    if rank <= 0:
        return [""]
    if rank == 1:
        return ["", "j"]
    if rank == 2:
        return ["", "i", "j", "k"]
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
# Parsing
# ============================================================================

_TERM_RE = re.compile(r"^([+-]?)(\d+/\d+|\d+\.\d+|\d+)?(i|j|k|e\d+)?$")


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
        elif label.startswith("e"):
            max_e_index = max(max_e_index, int(label[1:]))
        coeff_map[label] = coeff_map.get(label, Fraction(0)) + value

    if max_e_index > 0:
        n = 2
        while n - 1 < max_e_index:
            n *= 2
        rank = n.bit_length() - 1
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

    print("\nAll self-tests passed.")
