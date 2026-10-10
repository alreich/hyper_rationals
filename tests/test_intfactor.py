"""
test_intfactor.py
=================

Unit tests for intfactor.py: primality testing (deterministic Miller-Rabin
below 3.3e24, Baillie-PSW above) and integer factoring (trial division and
Pollard-Brent rho), compared against naive code on small inputs, known
pseudoprimes, and sympy when it is installed.
"""

import math
import random
import sys
import unittest
from unittest import mock

from hyprat import intfactor
from hyprat.intfactor import factorint, is_probable_prime, primes_upto, sqrt_mod_prime


def naive_is_prime(n):
    if n < 2:
        return False
    return all(n % d for d in range(2, math.isqrt(n) + 1))


def naive_factorint(n):
    out, d = {}, 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def product_of(factors):
    return math.prod(p ** e for p, e in factors.items())


class TestPrimesUpto(unittest.TestCase):

    def test_small(self):
        self.assertEqual(primes_upto(-5), [])
        self.assertEqual(primes_upto(1), [])
        self.assertEqual(primes_upto(2), [2])
        self.assertEqual(primes_upto(10), [2, 3, 5, 7])
        self.assertEqual(primes_upto(11), [2, 3, 5, 7, 11])

    def test_against_naive(self):
        self.assertEqual(primes_upto(3000), [n for n in range(3001) if naive_is_prime(n)])

    def test_count(self):
        self.assertEqual(len(primes_upto(10 ** 5)), 9592)


class TestIsProbablePrime(unittest.TestCase):

    def test_small_exhaustive(self):
        for n in range(-5, 6000):
            self.assertEqual(is_probable_prime(n), naive_is_prime(n), n)

    def test_strong_pseudoprimes_to_base_2(self):
        for n in (2047, 3277, 4033, 4681, 8321, 15841, 29341, 42799, 49141,
                  52633, 65281, 74665, 80581, 85489, 88357, 90751):
            self.assertFalse(is_probable_prime(n), n)

    def test_strong_pseudoprimes_to_many_bases(self):
        for n in (1373653, 25326001, 3215031751, 2152302898747,
                  3474749660383, 341550071728321, 3825123056546413051):
            self.assertFalse(is_probable_prime(n), n)

    def test_carmichael_numbers(self):
        for n in (561, 1105, 1729, 2465, 2821, 6601, 8911, 41041, 825265,
                  321197185, 5394826801, 232250619601, 9746347772161):
            self.assertFalse(is_probable_prime(n), n)

    def test_mersenne_numbers(self):
        prime_exponents = (2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607)
        for p in range(2, 130):
            self.assertEqual(is_probable_prime(2 ** p - 1), p in prime_exponents, p)
        for p in (521, 607):
            self.assertTrue(is_probable_prime(2 ** p - 1))

    def test_large_composites(self):
        p, q = 2 ** 61 - 1, 2 ** 89 - 1
        self.assertFalse(is_probable_prime(p * q))
        self.assertFalse(is_probable_prime(q * q))          # a square, above the MR limit
        self.assertFalse(is_probable_prime(2 ** 128 + 1))
        self.assertFalse(is_probable_prime((2 ** 61 - 1) ** 3))

    def test_both_algorithm_paths_agree(self):
        # the Baillie-PSW pieces, run directly on small numbers
        for n in range(3, 100000, 2):
            if math.isqrt(n) ** 2 == n:
                continue
            bpsw = (intfactor._strong_probable_prime(n, 2)
                    and intfactor._strong_lucas_probable_prime(n))
            self.assertEqual(bpsw, naive_is_prime(n), n)

    def test_random_against_sympy(self):
        try:
            import sympy
        except ImportError:
            self.skipTest("sympy not installed")
        rng = random.Random(1)
        for bits in range(20, 200, 7):
            for _ in range(25):
                n = rng.getrandbits(bits)
                self.assertEqual(is_probable_prime(n), sympy.isprime(n), n)

    def test_types(self):
        for bad in (2.0, "7", None, True):
            with self.assertRaises(TypeError):
                is_probable_prime(bad)


class TestFactorint(unittest.TestCase):

    def test_trivial(self):
        self.assertEqual(factorint(1), {})
        self.assertEqual(factorint(2), {2: 1})
        self.assertEqual(factorint(97), {97: 1})

    def test_small_exhaustive(self):
        for n in range(1, 3000):
            self.assertEqual(factorint(n), naive_factorint(n), n)

    def test_medium_against_naive(self):
        rng = random.Random(3)
        for _ in range(150):
            n = rng.randint(10 ** 5, 10 ** 11)
            self.assertEqual(factorint(n), naive_factorint(n), n)

    def test_reconstruction_and_order(self):
        rng = random.Random(4)
        for _ in range(80):
            n = rng.getrandbits(rng.randint(2, 80)) + 1
            f = factorint(n)
            self.assertEqual(product_of(f), n)
            self.assertEqual(list(f), sorted(f))
            self.assertTrue(all(is_probable_prime(p) and e >= 1 for p, e in f.items()))

    def test_prime_powers(self):
        self.assertEqual(factorint(2 ** 100), {2: 100})
        self.assertEqual(factorint(1009 ** 5), {1009: 5})
        self.assertEqual(factorint((2 ** 31 - 1) ** 2), {2 ** 31 - 1: 2})
        self.assertEqual(factorint(1000003 ** 3), {1000003: 3})

    def test_known(self):
        self.assertEqual(factorint(600851475143), {71: 1, 839: 1, 1471: 1, 6857: 1})
        self.assertEqual(factorint(2 ** 64 + 1), {274177: 1, 67280421310721: 1})
        self.assertEqual(factorint(2 ** 67 - 1), {193707721: 1, 761838257287: 1})
        self.assertEqual(factorint(10 ** 12 - 1),
                         {3: 3, 7: 1, 11: 1, 13: 1, 37: 1, 101: 1, 9901: 1})

    def test_semiprimes(self):
        p, q = 1000000007, 998244353
        self.assertEqual(factorint(p * q), {q: 1, p: 1})
        self.assertEqual(factorint(p * p), {p: 2})
        p, q = 1000000000039, 1000000000061
        self.assertEqual(factorint(p * q), {p: 1, q: 1})

    def test_carmichael(self):
        self.assertEqual(factorint(561), {3: 1, 11: 1, 17: 1})
        self.assertEqual(factorint(9746347772161), {7: 1, 11: 1, 13: 1, 17: 1, 19: 1,
                                                      31: 1, 37: 1, 41: 1, 641: 1})

    def test_errors(self):
        for bad in (0, -1, -100):
            with self.assertRaises(ValueError):
                factorint(bad)
        for bad in (2.0, "6", None, True):
            with self.assertRaises(TypeError):
                factorint(bad)
        with self.assertRaises(ValueError):
            factorint(6, method="magic")

    def test_methods_agree(self):
        try:
            import sympy  # noqa: F401
        except ImportError:
            self.skipTest("sympy not installed")
        rng = random.Random(8)
        for _ in range(40):
            n = rng.getrandbits(rng.randint(2, 70)) + 1
            self.assertEqual(factorint(n, method="python"), factorint(n, method="sympy"))
            self.assertEqual(factorint(n, method="auto"), factorint(n))
        big = (2 ** 89 - 1) * 1000003 * 12345
        self.assertEqual(factorint(big, method="auto"), factorint(big, method="python"))

    def test_without_sympy(self):
        with mock.patch.dict(sys.modules, {"sympy": None}):
            with self.assertRaises(ImportError):
                factorint(12, method="sympy")
            big = (2 ** 89 - 1) * 7
            self.assertEqual(factorint(big, method="auto"), {7: 1, 2 ** 89 - 1: 1})


class TestSqrtModPrime(unittest.TestCase):

    def test_exhaustive_small_primes(self):
        for p in primes_upto(400)[1:]:
            squares = {x * x % p for x in range(p)}
            for a in range(p):
                if a in squares:
                    r = sqrt_mod_prime(a, p)
                    self.assertEqual(r * r % p, a, (a, p))
                else:
                    with self.assertRaises(ValueError):
                        sqrt_mod_prime(a, p)

    def test_large_primes(self):
        # 998244353 - 1 = 119 * 2**23 and 2**64 - 59 with 2**64 - 60 = 4 * odd
        # exercise the Tonelli-Shanks loop; 2**61 - 1 is 3 mod 4
        rng = random.Random(12)
        for p in (998244353, 2 ** 64 - 59, 2 ** 61 - 1, 1000000007):
            for _ in range(20):
                a = rng.randrange(1, p)
                r = sqrt_mod_prime(a * a % p, p)
                self.assertIn(r, (a, p - a))

    def test_reduces_and_zero(self):
        self.assertEqual(sqrt_mod_prime(0, 13), 0)
        self.assertEqual(sqrt_mod_prime(13, 13), 0)
        r = sqrt_mod_prime(10 + 13 * 5, 13)
        self.assertEqual(r * r % 13, 10)
        r = sqrt_mod_prime(-3, 13)                  # -3 = 10 (mod 13)
        self.assertEqual(r * r % 13, 10)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
