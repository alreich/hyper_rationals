"""
hyprat.intfactor
================

Integer primality testing and factoring, in pure Python.

The factorization of a Hurwitz integer (see :mod:`hyprat.hurwitz`) starts
from the factorization of its norm, an ordinary integer, so this small
module provides what is needed and nothing more:

* :func:`is_probable_prime` -- exact below 3.3e24, and above that a
  Baillie-PSW test (a strong Fermat test to base 2 followed by a strong
  Lucas test).  No composite number is known that passes Baillie-PSW.
* :func:`sqrt_mod_prime` -- square roots modulo a prime (Tonelli-Shanks),
  needed to write a prime as a sum of four squares.
* :func:`factorint` -- trial division by the primes below 1000, then
  Brent's variant of Pollard's rho method, with the primality test above
  deciding when to stop.  If ``sympy`` happens to be installed it can be
  used for very large inputs (``method="sympy"``, or ``"auto"``); it is
  never required.

    >>> from hyprat.intfactor import factorint, is_probable_prime
    >>> factorint(360)
    {2: 3, 3: 2, 5: 1}
    >>> factorint(2**64 + 1)
    {274177: 1, 67280421310721: 1}
    >>> is_probable_prime(2**61 - 1)
    True
"""

import math

__all__ = ["is_probable_prime", "factorint", "primes_upto", "sqrt_mod_prime"]


def primes_upto(limit: int) -> list:
    """The primes ``<= limit``, by the sieve of Eratosthenes.

    >>> primes_upto(30)
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    """
    if limit < 2:
        return []
    sieve = bytearray([1]) * (limit + 1)
    sieve[0] = sieve[1] = 0
    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p::p] = bytearray(len(range(p * p, limit + 1, p)))
    return [i for i, flag in enumerate(sieve) if flag]


_SMALL_PRIMES = tuple(primes_upto(1000))

# Deterministic Miller-Rabin: these bases are enough for every n < 3.3e24
# (Sorenson & Webster, Math. Comp. 2017).
_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
_MR_LIMIT = 3317044064679887385961981


def _strong_probable_prime(n: int, a: int) -> bool:
    """Strong probable prime test of the odd number ``n > 2`` to base ``a``."""
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    x = pow(a, d, n)
    if x == 1 or x == n - 1:
        return True
    for _ in range(s - 1):
        x = x * x % n
        if x == n - 1:
            return True
    return False


def _jacobi(a: int, n: int) -> int:
    """The Jacobi symbol (a/n) for odd positive ``n``."""
    a %= n
    result = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0


def _half(x: int, n: int) -> int:
    """``x / 2`` modulo the odd number ``n``."""
    return (x + n) // 2 if x % 2 else x // 2


def _strong_lucas_probable_prime(n: int) -> bool:
    """Strong Lucas test (Selfridge's parameters) of an odd non-square ``n``."""
    d_param = 5
    while True:
        j = _jacobi(d_param, n)
        if j == -1:
            break
        if j == 0 and abs(d_param) != n:
            return False
        d_param = -d_param - 2 if d_param > 0 else -d_param + 2
    p, q = 1, (1 - d_param) // 4

    d, s = n + 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Lucas sequences U_d, V_d and Q^d modulo n, by binary doubling.
    u, v, qk = 1, p % n, q % n
    for bit in bin(d)[3:]:
        u, v, qk = u * v % n, (v * v - 2 * qk) % n, qk * qk % n
        if bit == "1":
            u, v = _half(p * u + v, n), _half(d_param * u + p * v, n)
            qk = qk * q % n
    if u % n == 0 or v % n == 0:
        return True
    for _ in range(s - 1):
        v = (v * v - 2 * qk) % n
        qk = qk * qk % n
        if v % n == 0:
            return True
    return False


def is_probable_prime(n: int) -> bool:
    """Whether ``n`` is prime.

    The answer is certain for ``n < 3.3e24``.  Beyond that it is the
    Baillie-PSW test, which has no known counterexample.

    >>> [n for n in range(30) if is_probable_prime(n)]
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    >>> is_probable_prime(3215031751)        # a strong pseudoprime to 2, 3, 5, 7
    False
    """
    if isinstance(n, bool) or not isinstance(n, int):
        raise TypeError(f"n must be an int, not {type(n).__name__}")
    if n < 2:
        return False
    for p in _SMALL_PRIMES[:25]:
        if n == p:
            return True
        if n % p == 0:
            return False
    if n < _MR_LIMIT:
        return all(_strong_probable_prime(n, a) for a in _MR_BASES)
    if not _strong_probable_prime(n, 2):
        return False
    r = math.isqrt(n)
    if r * r == n:
        return False
    return _strong_lucas_probable_prime(n)


def _pollard_brent(n: int) -> int:
    """A nontrivial factor of the odd composite ``n`` (Brent's variant)."""
    c = 1
    while True:
        y, r, q, g = 2, 1, 1, 1
        m = 128
        x = ys = y
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += m
            r *= 2
        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
        if g != n:
            return g
        c += 1


def _factor_rec(n: int, out: dict) -> None:
    if n == 1:
        return
    if is_probable_prime(n):
        out[n] = out.get(n, 0) + 1
        return
    r = math.isqrt(n)
    if r * r == n:
        factor = r
    else:
        factor = _pollard_brent(n)
    _factor_rec(factor, out)
    _factor_rec(n // factor, out)


def factorint(n: int, method: str = "python") -> dict:
    """The prime factorization of the positive integer ``n``, as a dict
    ``{prime: exponent}`` in increasing order of the primes.

    ``factorint(1)`` is ``{}``.  ``method`` is ``"python"`` (the default),
    ``"sympy"`` (requires sympy) or ``"auto"`` (sympy for inputs of more
    than 100 bits when it is installed, otherwise pure Python).  All three
    give the same answer; they differ only in speed on large inputs.

    >>> factorint(1)
    {}
    >>> factorint(600851475143)
    {71: 1, 839: 1, 1471: 1, 6857: 1}
    """
    if isinstance(n, bool) or not isinstance(n, int):
        raise TypeError(f"n must be an int, not {type(n).__name__}")
    if n < 1:
        raise ValueError(f"n must be positive, not {n}")
    if method not in ("python", "sympy", "auto"):
        raise ValueError(f"method must be 'python', 'sympy' or 'auto', not {method!r}")

    if method == "sympy" or (method == "auto" and n.bit_length() > 100):
        try:
            from sympy import factorint as sympy_factorint
        except ImportError:
            if method == "sympy":
                raise ImportError("method='sympy' needs sympy to be installed") from None
        else:
            return {int(p): int(e) for p, e in sorted(sympy_factorint(n).items())}

    out: dict = {}
    for p in _SMALL_PRIMES:
        if p * p > n:
            break
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
    if n > 1:
        if n < _SMALL_PRIMES[-1] ** 2:
            out[n] = out.get(n, 0) + 1       # no factor below 1000 and n < 997**2: prime
        else:
            _factor_rec(n, out)
    return dict(sorted(out.items()))


def sqrt_mod_prime(a: int, p: int) -> int:
    """A square root of ``a`` modulo the odd prime ``p`` (Tonelli-Shanks).

    Raises ``ValueError`` if ``a`` is not a square modulo ``p``.  The caller
    is responsible for ``p`` being an odd prime.

    >>> r = sqrt_mod_prime(10, 13)
    >>> r * r % 13
    10
    """
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        raise ValueError(f"{a} is not a square modulo {p}")
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, r = t * c % p, r * b % p
    return r
