"""
hyprat.hurwitz
==============

The Hurwitz integers: the quaternions ``a + b*i + c*j + d*k`` whose four
coordinates are either all integers or all halves of odd integers, so
that ``(1 + i + j + k)/2`` is one of them.  They form a ring, and the
class :class:`Hu` stores them exactly, as plain Python integers.

``Hu`` is to ``Hy`` what ``Zi`` is to ``Qi`` in the ``gint`` package:
``Hy`` is the ambient algebra over the rationals (and the oracle the
tests check against), and ``Hu`` holds the integer points inside it::

    >>> from hyprat import Hu
    >>> q = Hu(1, 2, 3, 4)
    >>> str(q)
    '(1+2i+3j+4k)'
    >>> w = Hu('1/2', '1/2', '1/2', '1/2')     # (1+i+j+k)/2, not a Lipschitz integer
    >>> w * w * w                              # a unit of order 6
    Hu(-1, 0, 0, 0)
    >>> q.norm(), q.conjugate()
    (30, Hu(1, -2, -3, -4))

Storage and invariants
----------------------

A ``Hu`` keeps its four coordinates **doubled**, as ints ``(A, B, C, D)``
that all have the same parity (all even: a Lipschitz integer; all odd:
the half-integer kind).  The value is ``(A + B*i + C*j + D*k) / 2``.  All
arithmetic is integer arithmetic, with no ``Fraction`` objects, and a value
that is not a Hurwitz integer cannot be constructed.

There is no ``/`` operator: the ring is not closed under division.  Use
:meth:`Hu.to_hy` for exact rational division, or the Euclidean operations
below.

Division, gcds and associates
-----------------------------

The ring is Euclidean on both sides: for ``b != 0`` there is a ``q`` with
``a = b*q + r`` and ``2*N(r) <= N(b)``, and likewise ``a = q*b + r``.
Because the ring is not commutative, everything comes in left and right
forms.  ``g`` is a *left divisor* of ``a`` when ``a = g*x``, and a *right
divisor* when ``a = x*g``.

* :meth:`Hu.divmod_left`, :meth:`Hu.divmod_right` -- division with remainder
* :meth:`Hu.left_divides`, :meth:`Hu.right_divides`,
  :meth:`Hu.div_exact_left`, :meth:`Hu.div_exact_right`
* :meth:`Hu.gcld`, :meth:`Hu.gcrd` and the extended forms
  :meth:`Hu.xgcld` (``g = a*x + b*y``) and :meth:`Hu.xgcrd`
  (``g = x*a + y*b``)
* :meth:`Hu.is_left_associate`, :meth:`Hu.is_right_associate`,
  :meth:`Hu.canonical_associate`

A gcld is determined only up to a unit on its right, a gcrd up to a unit on
its left; the methods return the canonical associate (the greatest doubled
coordinate tuple, so units map to 1)::

    >>> a, b = Hu(1, 2, 3, 5), Hu(1, 1, 1, 0)
    >>> q, r = a.divmod_left(b)
    >>> a == b * q + r and 2 * r.norm() <= b.norm()
    True

Primes, content and factorization
---------------------------------

A nonzero Hurwitz integer is *prime* exactly when its norm is an ordinary
prime (:meth:`Hu.is_prime`).  Its *content* is the largest positive integer
dividing it (:meth:`Hu.content`), and it is *primitive* when the content
is 1.  A primitive element factors into primes in a way that is almost
unique and, remarkably, lets you choose the order of the prime norms
(Conway and Smith, *On Quaternions and Octonions*, 2003)::

    >>> f = Hu(1, 2, 3, 4).factor()              # norm 30 = 2 * 3 * 5
    >>> f.norms()
    (2, 3, 5)
    >>> f.reorder([5, 2, 3]).norms()
    (5, 2, 3)
    >>> f.product()
    Hu(1, 2, 3, 4)

For a non-primitive element, the rational integer part is reported
separately, as ``content``: a rational prime ``p`` is ``pi * conj(pi)`` for
many different primes ``pi``, so there is no single answer to put in a
list.  Factoring needs the factorization of the norm, which is done by
:mod:`hyprat.intfactor`; pure Python copes with norms whose second-largest
prime factor is up to about 1e15, and sympy (if installed) is used for
larger inputs.

Enumeration and sums of four squares
------------------------------------

Exactly ``24 * sigma_odd(n)`` Hurwitz integers have norm ``n``, where
``sigma_odd(n)`` is the sum of the odd divisors of ``n``.  They can be
listed with :meth:`Hu.of_norm` and counted, without listing, by
:meth:`Hu.count_of_norm`; :meth:`Hu.primes_of_norm` gives the ``p + 1``
classes of primes of norm ``p``.  The same machinery proves Lagrange's
theorem constructively: :meth:`Hu.with_norm` builds an element with integer
coordinates and norm ``n``, so :meth:`Hu.four_squares` writes ``n`` as a sum
of four squares, even for huge ``n``::

    >>> Hu.count_of_norm(6)
    96
    >>> Hu.four_squares(1000003)
    (699, 699, 151, 0)

Converting to and from ``Hy``
-----------------------------

Conversions are explicit, and a ``Hu`` is never ``==`` to a ``Hy``::

    >>> from hyprat import Hy
    >>> Hu(1, 2, 3, 4).to_hy()
    Hy(Hy('1', '2'), Hy('3', '4'))
    >>> Hu.from_hy(Hy(Hy(1, 2), Hy(3, 4)))
    Hu(1, 2, 3, 4)
    >>> Hu.from_hy(Hy(2, -3))                  # a rank-1 value is a + b*i
    Hu(2, -3, 0, 0)

``Hy.is_hurwitz()``, ``Hy.is_lipschitz()`` and ``Hy.embed(rank)`` are the
matching methods on the ``Hy`` side.

Text
----

``str(Hu)`` and :meth:`Hu.parse` always use the **quaternion** labels
``i, j, k``.  That differs from ``Hy.parse``, which reads text that
mentions only ``j`` as a *rank-1* value (the way ``gint``'s ``Zi`` and
``Qi`` print), so ``Hy.parse('(3j)')`` is ``3*j`` in the complex numbers
while ``Hu('(3j)')`` is the quaternion ``3j``.  To read ``gint``-style text
as a Hurwitz integer, go through ``Hy`` explicitly::

    >>> Hu('(3j)')
    Hu(0, 0, 3, 0)
    >>> Hu.from_hy(Hy.parse('(3j)'))
    Hu(0, 3, 0, 0)
"""

from fractions import Fraction
from itertools import product
import math
import numbers
import re
from typing import NamedTuple

from .hypercomplex import Hy, _flatten, _DEFAULT_MU
from .intfactor import factorint, is_probable_prime, sqrt_mod_prime

__all__ = ["Hu", "HuFactorization"]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _double(v, what="coordinate") -> int:
    """``2*v`` as an int, for ``v`` an integer or a half-integer given as an
    int, ``Fraction``, float or string; ``ValueError``/``TypeError``
    otherwise."""
    if isinstance(v, numbers.Integral):
        return 2 * int(v)
    if isinstance(v, Fraction):
        f = v
    elif isinstance(v, str):
        try:
            f = Fraction(v.strip())
        except (ValueError, ZeroDivisionError):
            raise ValueError(f"cannot read {v!r} as a number") from None
    elif isinstance(v, float):
        if not math.isfinite(v):
            raise ValueError(f"{what} {v!r} is not finite")
        f = Fraction(str(v))
    else:
        raise TypeError(
            f"{what} must be an int, Fraction, float or string, "
            f"not {type(v).__name__}"
        )
    d = f * 2
    if d.denominator != 1:
        raise ValueError(
            f"{what} {v!r} is neither an integer nor half of an odd integer"
        )
    return int(d)


def _show(x2: int) -> str:
    """A doubled coordinate as a readable number: ``3`` -> ``3/2``."""
    return str(Fraction(x2, 2))


def _why_not_quaternion(h: Hy) -> str:
    """Why a ``Hy`` is outside the classical quaternion subalgebra."""
    signs = h.signs
    if signs[0] != _DEFAULT_MU or (h.rank >= 2 and signs[1] != _DEFAULT_MU):
        return ("its two lowest doubling levels do not both have mu = -1, "
                "so it is not in the classical quaternions")
    return ("it has nonzero coordinates beyond the first four, so it is "
            "not a quaternion")


_TERM_RE = re.compile(r"^([+-]?)(\d+(?:\.\d+|/\d+)?)?([ijk])?$")


def _parse_coords(text: str):
    """Read ``'1+2i-j+3/2k'``-style quaternion text into four Fractions."""
    s = text.strip()
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1].strip()
    s = re.sub(r"\s*([+-])\s*", r"\1", s)       # spaces are fine around + and -
    if s == "":
        raise ValueError("cannot parse an empty string as a Hu")
    if re.search(r"\s", s):
        raise ValueError(f"cannot parse {text!r} as a Hu: unexpected whitespace")
    tokens = re.findall(r"[+-]?[^+-]+", s)
    if "".join(tokens) != s:
        raise ValueError(f"cannot parse {text!r} as a Hu")
    coeffs = {"": Fraction(0), "i": Fraction(0), "j": Fraction(0), "k": Fraction(0)}
    for tok in tokens:
        m = _TERM_RE.match(tok)
        if not m or (m.group(2) is None and m.group(3) is None):
            raise ValueError(f"cannot parse term {tok!r} in {text!r} as a Hu")
        sign, num, unit = m.groups()
        try:
            value = Fraction(num) if num is not None else Fraction(1)
        except (ValueError, ZeroDivisionError):
            raise ValueError(f"cannot parse term {tok!r} in {text!r} as a Hu") from None
        if sign == "-":
            value = -value
        coeffs[unit or ""] += value
    return coeffs[""], coeffs["i"], coeffs["j"], coeffs["k"]


def _rebuild(cls, a, b, c, d):
    """Unpickling helper: rebuild from already-valid doubled coordinates."""
    return cls._raw(a, b, c, d)


_UNIT_CACHE: list = []


def _unit_values():
    """The 24 units as a tuple of ``Hu`` (built once)."""
    if not _UNIT_CACHE:
        _UNIT_CACHE.append(tuple(Hu.units().values()))
    return _UNIT_CACHE[0]


# ---------------------------------------------------------------------------
# the class
# ---------------------------------------------------------------------------

class Hu:
    """An exact Hurwitz integer ``a + b*i + c*j + d*k``: the four
    coordinates are all integers, or all halves of odd integers.

    ``Hu(a, b, c, d)`` takes the four coordinates, each an ``int``,
    ``Fraction``, float or string such as ``'1/2'``; missing ones are 0.
    ``Hu(x)`` with a single argument also accepts a string in the form
    that ``str`` produces, another ``Hu`` (returned as is) or a ``Hy``
    (see :meth:`from_hy`).  A set of coordinates that is not a Hurwitz
    integer, such as ``Hu('1/2', 0, 0, 0)``, raises ``ValueError``.

    Instances are immutable and hashable.  Arithmetic supports ``+``,
    ``-``, ``*`` and ``**`` with other ``Hu`` and with plain ints;
    multiplication is not commutative.  A ``Hu`` is never equal to a ``Hy``
    (convert with :meth:`to_hy` / :meth:`from_hy`), but ``Hu(5) == 5``.

        >>> Hu(1, 2, 3, 4) * Hu(0, 1, 0, 0)
        Hu(-2, 1, 4, -3)
        >>> Hu(0, 1, 0, 0) * Hu(0, 0, 1, 0), Hu(0, 0, 1, 0) * Hu(0, 1, 0, 0)
        (Hu(0, 0, 0, 1), Hu(0, 0, 0, -1))
        >>> Hu('1/2', '1/2', '1/2', '1/2').is_unit()
        True
        >>> Hu('1/2', 0, 0, 0)
        Traceback (most recent call last):
            ...
        ValueError: coordinates (1/2, 0, 0, 0) must be all integers or all halves of odd integers
    """

    __slots__ = ("_a", "_b", "_c", "_d")

    # Declared for readers and type checkers; the values are stored with
    # object.__setattr__ because the class is immutable.
    _a: int
    _b: int
    _c: int
    _d: int

    # ---------------------------------------------------------------- #
    # construction
    # ---------------------------------------------------------------- #
    def __new__(cls, a=0, b=None, c=None, d=None):
        if b is None and c is None and d is None:
            if type(a) is cls:
                return a
            if isinstance(a, Hu):
                return cls._raw(a._a, a._b, a._c, a._d)
            if isinstance(a, str):
                return cls.parse(a)
            if isinstance(a, Hy):
                return cls.from_hy(a)
            coords = (a, 0, 0, 0)
        else:
            coords = tuple(0 if v is None else v for v in (a, b, c, d))
        doubled = tuple(_double(v) for v in coords)
        return cls.from_doubled(*doubled)

    @classmethod
    def _raw(cls, a2, b2, c2, d2):
        """Build from doubled coordinates already known to be valid."""
        self = object.__new__(cls)
        object.__setattr__(self, "_a", a2)
        object.__setattr__(self, "_b", b2)
        object.__setattr__(self, "_c", c2)
        object.__setattr__(self, "_d", d2)
        return self

    @classmethod
    def from_doubled(cls, a2, b2, c2, d2) -> "Hu":
        """Build the value ``(a2 + b2*i + c2*j + d2*k) / 2`` from four ints
        of equal parity -- the way a ``Hu`` stores itself.

            >>> Hu.from_doubled(1, 1, 1, -1)
            Hu('1/2', '1/2', '1/2', '-1/2')
            >>> Hu.from_doubled(2, 4, 0, -6)
            Hu(1, 2, 0, -3)
        """
        vals = (a2, b2, c2, d2)
        for v in vals:
            if not isinstance(v, numbers.Integral) or isinstance(v, bool):
                raise TypeError(
                    f"from_doubled() needs ints, not {type(v).__name__}"
                )
        ints = tuple(int(v) for v in vals)
        if len({v & 1 for v in ints}) != 1:
            raise ValueError(
                "coordinates (" + ", ".join(_show(v) for v in vals) + ") must "
                "be all integers or all halves of odd integers"
            )
        return cls._raw(*ints)

    @classmethod
    def parse(cls, text: str) -> "Hu":
        """Read the text that ``str`` produces -- ``'(1+2i+3j+4k)'``,
        ``'(1/2+1/2i+1/2j+1/2k)'``, ``'j'``, ``'-3'`` -- using the
        quaternion labels ``i``, ``j``, ``k`` (a term with no coefficient
        means 1).  The outer parentheses are optional, spaces are allowed
        around ``+`` and ``-``, and coefficients may be integers, fractions
        or decimals.  Raises
        ``ValueError`` for text that is malformed or whose coordinates do
        not make a Hurwitz integer.

            >>> Hu.parse('(1/2-1/2i+1/2j+1/2k)')
            Hu('1/2', '-1/2', '1/2', '1/2')
            >>> Hu.parse('-j')
            Hu(0, 0, -1, 0)
            >>> Hu.parse(str(Hu(3, -1, 0, 7))) == Hu(3, -1, 0, 7)
            True
        """
        return cls(*_parse_coords(text))

    from_string = parse

    @classmethod
    def from_hy(cls, h) -> "Hu":
        """Convert a ``Hy`` whose value is a Hurwitz integer, exactly.

        The value must lie in the classical quaternion subalgebra (see
        ``Hy.is_hurwitz``): a rank-1 value ``a + b*j`` is read as the
        quaternion ``a + b*i``, a value of higher rank must have zero
        coordinates beyond the first four, and the coordinates must be
        all integers or all halves of odd integers.  Otherwise
        ``ValueError``; a non-``Hy`` raises ``TypeError``.

            >>> Hu.from_hy(Hy(Hy('1/2', '1/2'), Hy('1/2', '-3/2')))
            Hu('1/2', '1/2', '1/2', '-3/2')
            >>> Hu.from_hy(Hy(Hy(1, 2), Hy(0, 0)).embed(3))   # an octonion in H
            Hu(1, 2, 0, 0)
            >>> Hu.from_hy(Hy(Hy('1/2', 0), Hy(0, 0)))
            Traceback (most recent call last):
                ...
            ValueError: not a Hurwitz integer: the coordinates (1/2, 0, 0, 0) must be all integers or all halves of odd integers
        """
        if not isinstance(h, Hy):
            raise TypeError(f"from_hy() expects a Hy, not {type(h).__name__}")
        coords = h._hamilton_coords()
        if coords is None:
            raise ValueError(f"not a Hurwitz integer: {_why_not_quaternion(h)}")
        dens = {c.denominator for c in coords}
        if dens != {1} and dens != {2}:
            raise ValueError(
                "not a Hurwitz integer: the coordinates ("
                + ", ".join(str(c) for c in coords)
                + ") must be all integers or all halves of odd integers"
            )
        return cls._raw(*(int(c * 2) for c in coords))

    def to_hy(self) -> Hy:
        """This value as a rank-2 ``Hy``, exactly.

            >>> Hu('1/2', '1/2', '1/2', '1/2').to_hy()
            Hy(Hy('1/2', '1/2'), Hy('1/2', '1/2'))
        """
        return Hy.from_array(
            [Fraction(self._a, 2), Fraction(self._b, 2),
             Fraction(self._c, 2), Fraction(self._d, 2)]
        )

    # ---------------------------------------------------------------- #
    # coordinates
    # ---------------------------------------------------------------- #
    @property
    def doubled(self) -> tuple:
        """The stored coordinates ``(A, B, C, D)``, all ints of the same
        parity, with ``value == (A + B*i + C*j + D*k) / 2``."""
        return (self._a, self._b, self._c, self._d)

    @property
    def coords(self) -> tuple:
        """The four coordinates ``(a, b, c, d)`` as ``Fraction`` objects."""
        return tuple(Fraction(x, 2) for x in self.doubled)

    @property
    def a(self) -> Fraction:
        """The real coordinate."""
        return Fraction(self._a, 2)

    @property
    def b(self) -> Fraction:
        """The coordinate of ``i``."""
        return Fraction(self._b, 2)

    @property
    def c(self) -> Fraction:
        """The coordinate of ``j``."""
        return Fraction(self._c, 2)

    @property
    def d(self) -> Fraction:
        """The coordinate of ``k``."""
        return Fraction(self._d, 2)

    def __iter__(self):
        return iter(self.coords)

    def __getitem__(self, index):
        return self.coords[index]

    def is_lipschitz(self) -> bool:
        """True iff all four coordinates are integers (the Lipschitz
        integers, a subring of the Hurwitz integers).

            >>> Hu(1, 2, 3, 4).is_lipschitz(), Hu('1/2', '1/2', '1/2', '1/2').is_lipschitz()
            (True, False)
        """
        return self._a & 1 == 0

    # ---------------------------------------------------------------- #
    # text
    # ---------------------------------------------------------------- #
    def __str__(self) -> str:
        return str(self.to_hy())

    def __repr__(self) -> str:
        if self.is_lipschitz():
            body = ", ".join(str(x // 2) for x in self.doubled)
        else:
            body = ", ".join(repr(_show(x)) for x in self.doubled)
        return f"{type(self).__name__}({body})"

    # ---------------------------------------------------------------- #
    # immutability, copying, hashing, comparison
    # ---------------------------------------------------------------- #
    def __setattr__(self, name, value):
        raise AttributeError(f"{type(self).__name__} is immutable")

    def __delattr__(self, name):
        raise AttributeError(f"{type(self).__name__} is immutable")

    def __reduce__(self):
        return (_rebuild, (type(self),) + self.doubled)

    def __copy__(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def __eq__(self, other):
        if isinstance(other, Hu):
            return self.doubled == other.doubled
        if isinstance(other, numbers.Integral):
            return (self._b == self._c == self._d == 0
                    and self._a == 2 * int(other))
        if isinstance(other, Fraction):
            return self._b == self._c == self._d == 0 and Fraction(self._a, 2) == other
        return NotImplemented

    def __hash__(self):
        if self._b == self._c == self._d == 0:
            return hash(self._a // 2)          # same as hash(int), as == requires
        return hash(self.doubled)

    def __bool__(self) -> bool:
        return any(self.doubled)

    # ---------------------------------------------------------------- #
    # ring arithmetic
    # ---------------------------------------------------------------- #
    @staticmethod
    def _coerce(other):
        """``other`` as a ``Hu`` if it is one or an integer, else None."""
        if isinstance(other, Hu):
            return other
        if isinstance(other, numbers.Integral):
            return Hu._raw(2 * int(other), 0, 0, 0)
        if isinstance(other, Fraction) and other.denominator == 1:
            return Hu._raw(2 * int(other), 0, 0, 0)
        return None

    def __neg__(self):
        return self._raw(-self._a, -self._b, -self._c, -self._d)

    def __pos__(self):
        return self

    def __add__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        return self._raw(self._a + o._a, self._b + o._b,
                         self._c + o._c, self._d + o._d)

    __radd__ = __add__

    def __sub__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        return self._raw(self._a - o._a, self._b - o._b,
                         self._c - o._c, self._d - o._d)

    def __rsub__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        return o.__sub__(self)

    @staticmethod
    def _product(x, y):
        """Hamilton product of two doubled-coordinate 4-tuples, giving the
        doubled coordinates of the product (always an exact integer)."""
        a0, a1, a2, a3 = x
        b0, b1, b2, b3 = y
        return (
            (a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3) // 2,
            (a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2) // 2,
            (a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1) // 2,
            (a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0) // 2,
        )

    def __mul__(self, other):
        if isinstance(other, numbers.Integral):
            n = int(other)
            return self._raw(self._a * n, self._b * n, self._c * n, self._d * n)
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        return self._raw(*self._product(self.doubled, o.doubled))

    def __rmul__(self, other):
        o = self._coerce(other)
        if o is None:
            return NotImplemented
        return self._raw(*self._product(o.doubled, self.doubled))

    def __pow__(self, exponent, modulus=None):
        if modulus is not None or not isinstance(exponent, numbers.Integral):
            return NotImplemented
        n = int(exponent)
        base = self
        if n < 0:
            if not self.is_unit():
                raise ValueError(
                    "only units have inverses in the Hurwitz integers, so a "
                    "negative power needs a unit; use to_hy() ** n for "
                    "rational values"
                )
            base, n = self.conjugate(), -n
        result = self._raw(2, 0, 0, 0)
        while n:
            if n & 1:
                result = result * base
            n >>= 1
            if n:
                base = base * base
        return result

    # ---------------------------------------------------------------- #
    # conjugate, norm, trace
    # ---------------------------------------------------------------- #
    def conjugate(self) -> "Hu":
        """``a - b*i - c*j - d*k``.

            >>> Hu(1, 2, 3, 4).conjugate()
            Hu(1, -2, -3, -4)
        """
        return self._raw(self._a, -self._b, -self._c, -self._d)

    def norm(self) -> int:
        """The number-theoretic norm ``N(q) = q * conj(q) = a^2 + b^2 + c^2
        + d^2``, which is an ``int`` for every Hurwitz integer, and
        multiplicative: ``N(p*q) == N(p) * N(q)``.  The same quantity that
        ``Hy.norm()`` and ``Hy.norm_squared()`` give, as an int.

            >>> Hu(1, 2, 3, 4).norm()
            30
            >>> Hu('1/2', '1/2', '1/2', '1/2').norm()
            1
        """
        return (self._a ** 2 + self._b ** 2 + self._c ** 2 + self._d ** 2) // 4

    norm_squared = norm

    def trace(self) -> int:
        """``q + conj(q) = 2*a``, an ``int``: the stored real coordinate.

            >>> Hu(1, 2, 3, 4).trace(), Hu('1/2', '1/2', '1/2', '1/2').trace()
            (2, 1)
        """
        return self._a

    # ---------------------------------------------------------------- #
    # units
    # ---------------------------------------------------------------- #
    def is_unit(self) -> bool:
        """True iff this is one of the 24 units, the values of norm 1.

            >>> Hu(0, 0, -1, 0).is_unit(), Hu(1, 1, 0, 0).is_unit()
            (True, False)
        """
        return self.norm() == 1

    @classmethod
    def units(cls) -> dict:
        """The 24 units, as a dict from each unit's string form to the
        value, in the order ``1, -1, i, -i, j, -j, k, -k`` followed by the
        sixteen ``(+-1 +-i +-j +-k)/2`` (first sign ``+`` before ``-``).

            >>> list(Hu.units())[:8]
            ['1', '-1', 'i', '-i', 'j', '-j', 'k', '-k']
            >>> len(Hu.units()), all(u.is_unit() for u in Hu.units().values())
            (24, True)
            >>> Hu.units()['1/2-1/2i-1/2j+1/2k']
            Hu('1/2', '-1/2', '-1/2', '1/2')
        """
        result = {}
        for idx in range(4):
            for sign in (1, -1):
                v = [0, 0, 0, 0]
                v[idx] = 2 * sign
                u = cls._raw(*v)
                result[str(u)[1:-1]] = u
        for signs in product((1, -1), repeat=4):
            u = cls._raw(*signs)
            result[str(u)[1:-1]] = u
        return result

    # ---------------------------------------------------------------- #
    # exact division and divisibility
    #
    # The Hurwitz integers are not commutative, so "divides" has a side.
    # In this module, "g is a LEFT divisor of a" means a = g*x for some
    # Hurwitz integer x (g stands on the left of the product), and "g is
    # a RIGHT divisor of a" means a = x*g.
    # ---------------------------------------------------------------- #
    def _operand(self, other, role):
        o = self._coerce(other)
        if o is None:
            raise TypeError(
                f"the {role} must be a Hu or an int (convert a Hy with "
                f"Hu.from_hy), not {type(other).__name__}"
            )
        return o

    def _exact(self, divisor, left):
        """``x`` with ``self == divisor*x`` (``left`` true) or ``self ==
        x*divisor`` (``left`` false), or None if there is no such Hurwitz
        integer.  ``divisor`` must be a nonzero ``Hu``."""
        n = divisor.conjugate() * self if left else self * divisor.conjugate()
        norm = divisor.norm()
        quotient = [v // norm for v in n.doubled]
        if any(v % norm for v in n.doubled):
            return None
        if len({v & 1 for v in quotient}) != 1:
            return None
        return self._raw(*quotient)

    def left_divides(self, other) -> bool:
        """True iff this value is a *left divisor* of ``other``: ``other ==
        self * x`` for some Hurwitz integer ``x``.  Zero divides only zero.

            >>> g, x = Hu(1, 1, 1, 0), Hu(0, 0, 1, 1)
            >>> g * x
            Hu(-1, 1, 0, 2)
            >>> g.left_divides(g * x), g.right_divides(g * x)
            (True, False)
        """
        o = self._operand(other, "dividend")
        if not self:
            return not o
        return o._exact(self, left=True) is not None

    def right_divides(self, other) -> bool:
        """True iff this value is a *right divisor* of ``other``: ``other ==
        x * self`` for some Hurwitz integer ``x``.  Zero divides only zero.
        """
        o = self._operand(other, "dividend")
        if not self:
            return not o
        return o._exact(self, left=False) is not None

    def div_exact_left(self, divisor) -> "Hu":
        """The ``x`` with ``self == divisor * x`` (``divisor`` is a left
        divisor of ``self``).  Raises ``ValueError`` if there is none and
        ``ZeroDivisionError`` for a zero divisor.

            >>> a, g = Hu(1, 2, 3, 4), Hu(1, 1, 1, 0)
            >>> (g * a).div_exact_left(g) == a
            True
            >>> a.div_exact_left(Hu(2))
            Traceback (most recent call last):
                ...
            ValueError: (2) is not a left divisor of (1+2i+3j+4k)
        """
        g = self._operand(divisor, "divisor")
        if not g:
            raise ZeroDivisionError("division by zero in the Hurwitz integers")
        x = self._exact(g, left=True)
        if x is None:
            raise ValueError(f"{g} is not a left divisor of {self}")
        return x

    def div_exact_right(self, divisor) -> "Hu":
        """The ``x`` with ``self == x * divisor`` (``divisor`` is a right
        divisor of ``self``).  Raises ``ValueError`` if there is none and
        ``ZeroDivisionError`` for a zero divisor.
        """
        g = self._operand(divisor, "divisor")
        if not g:
            raise ZeroDivisionError("division by zero in the Hurwitz integers")
        x = self._exact(g, left=False)
        if x is None:
            raise ValueError(f"{g} is not a right divisor of {self}")
        return x

    # ---------------------------------------------------------------- #
    # division with remainder
    # ---------------------------------------------------------------- #
    @staticmethod
    def _nearest_hurwitz(n2, norm):
        """Doubled coordinates of a Hurwitz integer nearest to the rational
        point ``t`` with ``t[i] = n2[i] / (2 * norm)``, using integer
        arithmetic only.

        There are two kinds of candidates: the nearest Lipschitz integer
        (round each coordinate to an integer) and the nearest all-halves
        point (round each coordinate to a half-odd-integer).  Whichever is
        closer wins, and the Lipschitz one wins a tie; within each kind a
        coordinate exactly halfway rounds up.
        """
        two_norm = 2 * norm
        lip = tuple(2 * ((x + norm) // two_norm) for x in n2)
        half = tuple(2 * (x // two_norm) + 1 for x in n2)
        d_lip = sum((x - norm * c) ** 2 for x, c in zip(n2, lip))
        d_half = sum((x - norm * c) ** 2 for x, c in zip(n2, half))
        return lip if d_lip <= d_half else half

    def divmod_left(self, divisor):
        """``(q, r)`` with ``self == divisor * q + r`` and ``2 * r.norm()
        <= divisor.norm()``: the divisor stands on the left of the
        quotient.  The quotient is the Hurwitz integer nearest to the
        exact rational quotient ``divisor**-1 * self``.

        Unlike integer division, the pair is not always unique: when the
        exact quotient is equally close to several Hurwitz integers, this
        picks one deterministically (a Lipschitz integer before a
        half-integer point, and a coordinate exactly halfway rounds up).
        Raises ``ZeroDivisionError`` for a zero divisor.

            >>> a, b = Hu(1, 2, 3, 5), Hu(1, 1, 1, 0)
            >>> q, r = a.divmod_left(b)
            >>> q, r
            (Hu(2, -1, 2, 1), Hu(0, 0, 0, 1))
            >>> b * q + r == a
            True
            >>> 2 * r.norm() <= b.norm()
            True
        """
        b = self._operand(divisor, "divisor")
        if not b:
            raise ZeroDivisionError("division by zero in the Hurwitz integers")
        n = b.conjugate() * self
        q = self._raw(*self._nearest_hurwitz(n.doubled, b.norm()))
        return q, self - b * q

    def divmod_right(self, divisor):
        """``(q, r)`` with ``self == q * divisor + r`` and ``2 * r.norm()
        <= divisor.norm()``: the divisor stands on the right of the
        quotient.  See :meth:`divmod_left` for what is chosen when the pair
        is not unique.  Raises ``ZeroDivisionError`` for a zero divisor.

            >>> a, b = Hu(1, 2, 3, 5), Hu(1, 1, 1, 0)
            >>> q, r = a.divmod_right(b)
            >>> q, r
            (Hu(2, 2, -1, 2), Hu(0, 0, 0, 0))
            >>> q * b + r == a
            True
        """
        b = self._operand(divisor, "divisor")
        if not b:
            raise ZeroDivisionError("division by zero in the Hurwitz integers")
        n = self * b.conjugate()
        q = self._raw(*self._nearest_hurwitz(n.doubled, b.norm()))
        return q, self - q * b

    # ---------------------------------------------------------------- #
    # associates
    #
    # u*a is a LEFT associate of a (unit on the left), a*u a RIGHT
    # associate (unit on the right); every nonzero value has 24 of each.
    # ---------------------------------------------------------------- #
    def is_left_associate(self, other) -> bool:
        """True iff ``other == u * self`` for a unit ``u`` (the unit stands
        on the left).  Zero is an associate only of zero.

            >>> a = Hu(1, 2, 3, 4)
            >>> a.is_left_associate(Hu(0, 1, 0, 0) * a)
            True
            >>> a.is_left_associate(a * Hu(0, 1, 0, 0))
            False
        """
        o = self._operand(other, "other value")
        if not self:
            return not o
        return (bool(o) and o.norm() == self.norm()
                and o._exact(self, left=False) is not None)

    def is_right_associate(self, other) -> bool:
        """True iff ``other == self * u`` for a unit ``u`` (the unit stands
        on the right).  Zero is an associate only of zero."""
        o = self._operand(other, "other value")
        if not self:
            return not o
        return (bool(o) and o.norm() == self.norm()
                and o._exact(self, left=True) is not None)

    def canonical_associate(self, side: str) -> "Hu":
        """The standard representative of this value's class of associates.

        ``side="left"`` chooses among the 24 values ``u * self`` (the unit
        on the left); ``side="right"`` among ``self * u``.  The
        representative is the one with the greatest doubled coordinates,
        compared as ``(a, b, c, d)`` in that order -- so the greatest real
        part first -- which makes every unit's representative ``1`` and
        every positive integer its own.  Zero is its own representative.

            >>> Hu(-3).canonical_associate("left")
            Hu(3, 0, 0, 0)
            >>> Hu(0, 0, 1, 0).canonical_associate("right")
            Hu(1, 0, 0, 0)
            >>> a = Hu(1, 2, 3, 4)
            >>> a.canonical_associate("left") == (Hu(0, 1, 0, 0) * a).canonical_associate("left")
            True
        """
        if side not in ("left", "right"):
            raise ValueError(f"side must be 'left' or 'right', not {side!r}")
        if not self:
            return self
        if side == "left":
            products = (u * self for u in _unit_values())
        else:
            products = (self * u for u in _unit_values())
        best = max(products, key=lambda x: x.doubled)
        return self._raw(*best.doubled)

    # ---------------------------------------------------------------- #
    # greatest common divisors (Euclid's algorithm on each side)
    # ---------------------------------------------------------------- #
    def xgcld(self, other):
        """``(g, x, y)`` with ``g`` the greatest common *left* divisor of
        ``self`` and ``other`` -- ``self == g * a'`` and ``other == g * b'``
        -- and Bezout coefficients on the right: ``g == self * x + other *
        y``.  Every common left divisor of the two is a left divisor of
        ``g``.

        ``g`` is determined only up to a unit on its right (``g * u`` is
        equally good), so the canonical representative is returned:
        ``g.canonical_associate("right")``.  ``g`` is 0 only when both
        values are 0.

            >>> a, b = Hu(-2, 1, 3, 3), Hu(3, -3, -1, -3)    # norms 23 and 28
            >>> g, x, y = a.xgcld(b)
            >>> g, x, y
            (Hu(1, 0, 0, 0), Hu(0, -1, -1, -1), Hu(0, -1, 0, -1))
            >>> a * x + b * y == g
            True
        """
        b = self._operand(other, "other value")
        r0, s0, t0 = self, self._raw(2, 0, 0, 0), self._raw(0, 0, 0, 0)
        r1, s1, t1 = b, t0, s0
        while r1:
            q, r2 = r0.divmod_left(r1)
            r0, r1 = r1, r2
            s0, s1 = s1, s0 - s1 * q
            t0, t1 = t1, t0 - t1 * q
        if r0:
            unit = r0.canonical_associate("right").div_exact_left(r0)
            r0, s0, t0 = r0 * unit, s0 * unit, t0 * unit
        return r0, s0, t0

    def xgcrd(self, other):
        """``(g, x, y)`` with ``g`` the greatest common *right* divisor of
        ``self`` and ``other`` -- ``self == a' * g`` and ``other == b' * g``
        -- and Bezout coefficients on the left: ``g == x * self + y *
        other``.  Every common right divisor of the two is a right divisor
        of ``g``.

        ``g`` is determined only up to a unit on its left, so the canonical
        representative ``g.canonical_associate("left")`` is returned.
        """
        b = self._operand(other, "other value")
        r0, s0, t0 = self, self._raw(2, 0, 0, 0), self._raw(0, 0, 0, 0)
        r1, s1, t1 = b, t0, s0
        while r1:
            q, r2 = r0.divmod_right(r1)
            r0, r1 = r1, r2
            s0, s1 = s1, s0 - q * s1
            t0, t1 = t1, t0 - q * t1
        if r0:
            unit = r0.canonical_associate("left").div_exact_right(r0)
            r0, s0, t0 = unit * r0, unit * s0, unit * t0
        return r0, s0, t0

    def gcld(self, other) -> "Hu":
        """The greatest common left divisor of ``self`` and ``other``, in
        its canonical form (see :meth:`xgcld`).

            >>> a, b, c = Hu(-2, 1, 3, 3), Hu(3, -3, -1, -3), Hu(1, 1, 1, 0)
            >>> a.gcld(b)                       # a and b have no common divisor
            Hu(1, 0, 0, 0)
            >>> (c * a).gcld(c * b) == c.canonical_associate("right")
            True
        """
        return self.xgcld(other)[0]

    def gcrd(self, other) -> "Hu":
        """The greatest common right divisor of ``self`` and ``other``, in
        its canonical form (see :meth:`xgcrd`).

            >>> a, b, c = Hu(-2, 1, 3, 3), Hu(3, -3, -1, -3), Hu(1, 1, 1, 0)
            >>> (a * c).gcrd(b * c) == c.canonical_associate("left")
            True
        """
        return self.xgcrd(other)[0]

    # ------------------------------------------------------------------
    # primes, content and factorization
    # ------------------------------------------------------------------

    def is_prime(self) -> bool:
        """Whether ``self`` is a Hurwitz prime.

        A Hurwitz integer is prime (irreducible) exactly when its norm is an
        ordinary prime number.  Units and zero are not prime.

            >>> Hu(1, 1, 0, 0).is_prime()            # norm 2
            True
            >>> Hu(1, 2, 3, 4).is_prime()            # norm 30
            False
            >>> Hu(2).is_prime()                     # norm 4: 2 is not prime here
            False
        """
        return is_probable_prime(self.norm())

    def content(self) -> int:
        """The largest positive integer ``m`` such that ``self / m`` is again
        a Hurwitz integer.

        Rational integers are central, so it does not matter on which side
        they divide.  Zero has no content (``ValueError``).

            >>> Hu(2, 4, 6, 8).content()
            2
            >>> Hu(1, 1, 1, 1).content()             # 2 * (1+i+j+k)/2
            2
            >>> Hu(1, 2, 3, 4).content()
            1
        """
        if not self:
            raise ValueError("zero has no content")
        g = math.gcd(math.gcd(self._a, self._b), math.gcd(self._c, self._d))
        # g divides every doubled coordinate.  The quotients have a common
        # parity unless g is even and they are mixed; then g/2 works.
        if g % 2 == 0 and any((x // g) % 2 == 0 for x in self.doubled):
            g //= 2
        return g

    def primitive_part(self) -> "Hu":
        """``self`` divided by its :meth:`content`: a primitive Hurwitz
        integer, meaning one not divisible by any integer greater than 1.

            >>> Hu(2, 4, 6, 8).primitive_part()
            Hu(1, 2, 3, 4)
        """
        m = self.content()
        return self._raw(self._a // m, self._b // m, self._c // m, self._d // m)

    def is_primitive(self) -> bool:
        """Whether ``self`` is nonzero with :meth:`content` 1."""
        return bool(self) and self.content() == 1

    def factor(self, order=None) -> "HuFactorization":
        """Factor ``self`` into Hurwitz primes (Conway and Smith).

        Write ``self = m * q`` with ``m = content`` and ``q`` primitive.
        Then for **any** ordering ``p1, ..., pk`` of the prime factors of
        ``N(q)`` (with multiplicity) there are primes ``pi_1, ..., pi_k`` with
        ``N(pi_i) = p_i`` and ``q = pi_1 * ... * pi_k``, and the factors are
        unique up to *unit migration*: replacing ``pi_i`` by
        ``u_(i-1)^-1 * pi_i * u_i`` (with ``u_0 = u_k = 1``) changes nothing.

        ``order`` is the wanted sequence of norms ``p1, ..., pk``; by default
        the primes in increasing order.  The factors are found one at a
        time as ``pi_1 = gcld(q, p1)``, then the same for the quotient; all
        but the last are canonical right associates, and the last absorbs
        whatever unit is left over.

        The rational integer ``m`` is **not** split into Hurwitz primes (each
        rational prime ``p`` is ``pi * conj(pi)`` in several ways); it is
        returned as :attr:`HuFactorization.content`.  Zero cannot be factored.

            >>> f = Hu(1, 2, 3, 4).factor()          # norm 30 = 2 * 3 * 5
            >>> [p.norm() for p in f.factors]
            [2, 3, 5]
            >>> f.product() == Hu(1, 2, 3, 4)
            True
            >>> g = Hu(1, 2, 3, 4).factor(order=[5, 3, 2])
            >>> [p.norm() for p in g.factors]
            [5, 3, 2]
            >>> g.product() == Hu(1, 2, 3, 4)
            True
        """
        if not self:
            raise ValueError("zero cannot be factored")
        content = self.content()
        q = self.primitive_part()
        primes = []
        for p, e in factorint(q.norm(), method="auto").items():
            primes.extend([p] * e)
        if order is None:
            order = primes
        else:
            order = [int(p) if isinstance(p, numbers.Integral) else p for p in order]
            if sorted(order) != primes:
                raise ValueError(
                    f"order must be a rearrangement of the prime factors "
                    f"{primes} of the norm, not {order}")
        if not order:
            return HuFactorization(content, q, ())
        factors = []
        for p in order[:-1]:
            pi = q.gcld(p)
            assert pi.norm() == p, "gcld of a primitive element and p has norm p"
            factors.append(pi)
            q = q.div_exact_left(pi)
        factors.append(q)               # norm is the last prime: itself prime
        return HuFactorization(content, Hu(1), tuple(factors))

    # ------------------------------------------------------------------
    # enumeration and sums of four squares
    # ------------------------------------------------------------------

    @classmethod
    def of_norm(cls, n: int, primitive: bool = False) -> list:
        """All Hurwitz integers of norm ``n``, as a list in increasing order
        of their doubled coordinates.

        There are ``24 * sigma_odd(n)`` of them, where ``sigma_odd(n)`` is
        the sum of the odd divisors of ``n`` (see :meth:`count_of_norm`).
        With ``primitive=True`` only those of content 1 are returned.  The
        work grows roughly like ``n``, so this is meant for modest ``n``.

            >>> len(Hu.of_norm(1)), len(Hu.of_norm(2)), len(Hu.of_norm(3))
            (24, 24, 96)
            >>> sorted(h.norm() for h in Hu.of_norm(5))[:2]
            [5, 5]
            >>> len(Hu.of_norm(4)), len(Hu.of_norm(4, primitive=True))
            (24, 0)
        """
        if isinstance(n, bool) or not isinstance(n, numbers.Integral):
            raise TypeError(f"n must be an int, not {type(n).__name__}")
        n = int(n)
        if n < 0:
            raise ValueError(f"n must be non-negative, not {n}")
        if n == 0:
            return [cls(0)] if not primitive else []
        four_n = 4 * n
        bound = math.isqrt(four_n)
        found = []
        for parity in (0, 1):
            by_sum = {}
            for x in range(-bound, bound + 1):
                if x % 2 != parity:
                    continue
                ybound = math.isqrt(four_n - x * x)
                for y in range(-ybound, ybound + 1):
                    if y % 2 == parity:
                        by_sum.setdefault(x * x + y * y, []).append((x, y))
            for total, pairs in by_sum.items():
                partners = by_sum.get(four_n - total)
                if partners:
                    for a, b in pairs:
                        for c, d in partners:
                            found.append((a, b, c, d))
        found.sort()
        out = [cls._raw(*v) for v in found]
        if primitive:
            out = [h for h in out if h.content() == 1]
        return out

    @staticmethod
    def count_of_norm(n: int) -> int:
        """How many Hurwitz integers have norm ``n``: ``24 * sigma_odd(n)``
        (Hurwitz), computed from the factorization of ``n`` without listing
        them.  The count for ``n = 0`` is 1.

            >>> [Hu.count_of_norm(n) for n in range(1, 8)]
            [24, 24, 96, 24, 144, 96, 192]
            >>> Hu.count_of_norm(2 ** 40)
            24
            >>> Hu.count_of_norm(10 ** 18)          # 24 * (5**19 - 1) / 4
            114440917968744
        """
        if isinstance(n, bool) or not isinstance(n, numbers.Integral):
            raise TypeError(f"n must be an int, not {type(n).__name__}")
        n = int(n)
        if n < 0:
            raise ValueError(f"n must be non-negative, not {n}")
        if n == 0:
            return 1
        total = 24
        for p, e in factorint(n, method="auto").items():
            if p != 2:
                total *= (p ** (e + 1) - 1) // (p - 1)
        return total

    @classmethod
    def primes_of_norm(cls, p: int, associates: bool = False) -> list:
        """The Hurwitz primes of norm ``p``, for a rational prime ``p``.

        By default one representative of each class of right associates
        (a canonical one, see :meth:`canonical_associate`): there are
        ``p + 1`` classes for odd ``p`` and a single class for ``p = 2``.
        With ``associates=True`` all ``24 * (p + 1)`` (or 24) primes of that
        norm are returned.  The work grows roughly like ``p``.

            >>> [len(Hu.primes_of_norm(p)) for p in (2, 3, 5, 7)]
            [1, 4, 6, 8]
            >>> [len(Hu.primes_of_norm(p, associates=True)) for p in (2, 3, 5)]
            [24, 96, 144]
        """
        if not is_probable_prime(p):
            raise ValueError(f"{p} is not a prime")
        every = cls.of_norm(p)
        if associates:
            return every
        reps = {h.canonical_associate("right").doubled: h.canonical_associate("right")
                for h in every}
        return [reps[key] for key in sorted(reps)]

    @classmethod
    def _lipschitz_prime(cls, p: int) -> "Hu":
        """A prime of norm ``p`` with integer coordinates."""
        if p == 2:
            return cls(1, 1, 0, 0)
        # x^2 + y^2 + 1 = 0 (mod p) has a solution; alpha = x + y*i + j is
        # primitive with p dividing its norm, so gcld(alpha, p) has norm p
        x = 0
        while True:
            t = (-1 - x * x) % p
            if t == 0 or pow(t, (p - 1) // 2, p) == 1:
                y = sqrt_mod_prime(t, p)
                break
            x += 1
        pi = cls(x, y, 1, 0).gcld(p)
        assert pi.norm() == p
        for u in _unit_values():
            if (pi * u).is_lipschitz():
                return pi * u
        raise AssertionError("no unit makes the prime Lipschitz")  # pragma: no cover

    @classmethod
    def with_norm(cls, n: int) -> "Hu":
        """A Hurwitz integer with integer coordinates and norm ``n >= 0``.

        This is Lagrange's four-square theorem: the four coordinates are
        integers whose squares add up to ``n``.  It is built from primes: a
        prime ``p`` gets ``gcld(x + y*i + j, p)`` where
        ``x**2 + y**2 + 1 = 0 (mod p)``, and the answers are multiplied
        together, since norms multiply.  The cost is that of factoring ``n``.

            >>> h = Hu.with_norm(2026)
            >>> h.norm(), h.is_lipschitz()
            (2026, True)
        """
        if isinstance(n, bool) or not isinstance(n, numbers.Integral):
            raise TypeError(f"n must be an int, not {type(n).__name__}")
        n = int(n)
        if n < 0:
            raise ValueError(f"n must be non-negative, not {n}")
        result = cls(1) if n else cls(0)
        if n:
            for p, e in factorint(n, method="auto").items():
                pi = cls._lipschitz_prime(p)
                for _ in range(e):
                    result = result * pi
        return result

    @classmethod
    def four_squares(cls, n: int) -> tuple:
        """Non-negative integers ``(a, b, c, d)``, in decreasing order, with
        ``a**2 + b**2 + c**2 + d**2 == n``.

            >>> Hu.four_squares(310)
            (13, 10, 5, 4)
            >>> a, b, c, d = Hu.four_squares(10 ** 30 + 7)
            >>> a * a + b * b + c * c + d * d == 10 ** 30 + 7
            True
        """
        h = cls.with_norm(n)
        return tuple(sorted((abs(x) // 2 for x in h.doubled), reverse=True))


class HuFactorization(NamedTuple):
    """The result of :meth:`Hu.factor`: ``value == content * unit * f1 * ... * fk``.

    ``content`` is a positive int, ``factors`` a tuple of Hurwitz primes in
    the order requested, and ``unit`` is ``1`` unless there are no prime
    factors at all (then the value is ``content`` times a unit).
    """

    content: int
    unit: "Hu"
    factors: tuple

    def product(self) -> "Hu":
        """The value that was factored, rebuilt from the pieces.

            >>> f = Hu(3, 1, 4, 1).factor()
            >>> f.product()
            Hu(3, 1, 4, 1)
        """
        acc = self.content * self.unit
        for p in self.factors:
            acc = acc * p
        return acc

    def norms(self) -> tuple:
        """The norms of the factors, in order."""
        return tuple(p.norm() for p in self.factors)

    def reorder(self, order) -> "HuFactorization":
        """The factorization of the same value with the primes in a new order.

            >>> f = Hu(3, 1, 4, 1).factor()
            >>> f.reorder(reversed(f.norms())).product() == f.product()
            True
        """
        return self.product().factor(list(order))
