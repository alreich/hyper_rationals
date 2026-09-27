"""
test_hypercomplex.py
=====================

Unit tests for hypercomplex.py, covering the Hy class's constructor,
accessors, operators, and helper functions/methods at rank 0 (plain
Fraction) through rank 4 (2**4 == 16 real coordinates -- "sedenions").

Organization:

    TestConstructor        - Hy(...) normalization, embedding, copy
                              semantics, error cases, ranks 1-4
    TestAccessors           - .real .imag .rank .dimension .components()
                              .conjugate() .inverse() .norm() .is_zero()
    TestArithmeticFixed     - +,-,*,/ on hand-checked rank 1-4 examples
    TestArithmeticFuzz      - randomized algebraic-law fuzz tests,
                              ranks 0-4
    TestComparisonAndHash   - __eq__, __ne__, __hash__, cross-rank equality
    TestSequenceProtocol    - __iter__, __getitem__, __len__, __bool__
    TestConversions         - __abs__, __complex__, __pow__
    TestStringForms         - __str__, __repr__ (and repr round trip)
    TestUnits               - Hy.units(rank), .is_unit(), ranks 0-4
    TestLatex               - .latex(), vinculum=/mode= options, ranks 1-4
    TestParsing             - Hy.parse / Hy.from_string, ranks 1-4, errors
    TestModuleFunctions     - add, sub, neg, conj, mul, abs2, inverse, div
                              called directly (incl. on raw Fractions,
                              i.e. "rank 0")
    TestImmutability        - __slots__/locked-down __setattr__
    TestRandom              - Hy.random(rank), Hy.seed(), rng=/seed= kwargs
    TestArrayConversion     - Hy.from_array()/.to_array(), ranks 1-4,
                              mixed-type input, error cases
    TestFloatToFraction     - the float -> Fraction policy used by the
                              from_numpy_quaternion/from_quaternionic/
                              from_sympy (Float components) methods
    TestSympyToFraction     - the sympy-component -> Fraction policy used
                              by from_sympy (exact Rational/Integer, or
                              the same Float policy as above)
    TestInteropWithoutPackages
                            - interop methods exercised with stand-in
                              objects / a "missing package" simulation,
                              so they run even where sympy,
                              numpy-quaternion and quaternionic are not
                              installed
    TestSympyInterop        - Hy.to_sympy()/Hy.from_sympy()      (needs sympy)
    TestNumpyQuaternionInterop
                            - Hy.to_numpy_quaternion()/
                              Hy.from_numpy_quaternion()   (needs numpy-quaternion)
    TestQuaternionicInterop - Hy.to_quaternionic()/
                              Hy.from_quaternionic()       (needs quaternionic)
    TestInteropAcrossPackages
                            - all three packages agree with each other
                              (and with Hy) on the same value

The tests that need a third-party package are skipped (not failed) when
that package is not installed.

Algebraic laws are checked empirically via seeded random fuzzing
(deterministic across runs) in addition to fixed, hand-verified examples.
"""

import math
import random
import sys
import unittest
from fractions import Fraction
from unittest import mock

# Optional third-party quaternion packages, used only by the interoperability
# tests (which are skipped if the relevant package is missing).
try:
    import numpy as np
except ImportError:                                     # pragma: no cover
    np = None

try:
    import sympy
except ImportError:                                     # pragma: no cover
    sympy = None

try:
    import quaternion as npquat                         # the numpy-quaternion package
    if not hasattr(npquat, "quaternion"):               # some other 'quaternion' module
        npquat = None
except ImportError:                                     # pragma: no cover
    npquat = None

try:
    import quaternionic
except ImportError:                                     # pragma: no cover
    quaternionic = None

HAVE_SYMPY = sympy is not None
HAVE_NPQUAT = np is not None and npquat is not None
HAVE_QUATERNIONIC = np is not None and quaternionic is not None

from hyprat import Hy
from hyprat.hypercomplex import (
    add,
    sub,
    neg,
    conj,
    mul,
    abs2,
    inverse,
    div,
    _rank,
    _embed,
    _flatten,
    _unflatten,
    _is_zero_val,
    _parse,
    _latex_unit_label,
    _format_fraction_latex,
    _float_to_fraction,
    _sympy_to_fraction,
)


# ---------------------------------------------------------------------------
# Test fixtures / helpers
# ---------------------------------------------------------------------------

def rand_fraction(rng, lo=-9, hi=9, dmax=6):
    """A small random Fraction, using only the public rng interface."""
    num = rng.randint(lo, hi)
    den = rng.randint(1, dmax)
    return Fraction(num, den)


def rand_value(rank, rng):
    """A random value of the given rank: a Fraction for rank 0, else a Hy
    built purely through the public Hy(...) constructor."""
    if rank == 0:
        return rand_fraction(rng)
    return Hy(rand_value(rank - 1, rng), rand_value(rank - 1, rng))


def as_hy(x):
    """Wrap a raw value (Fraction or Hy) as a Hy for uniform assertions."""
    return x if isinstance(x, Hy) else Hy(x)


RANKS = (0, 1, 2, 3, 4)
FUZZ_TRIALS = 25


# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

class TestConstructor(unittest.TestCase):

    def test_from_fraction_int_float_str(self):
        self.assertEqual(Hy(Fraction(5, 2)).real, Fraction(5, 2))
        self.assertEqual(Hy(3).real, Fraction(3))
        self.assertEqual(Hy(3).imag, Fraction(0))
        self.assertEqual(Hy(2.5).real, Fraction(5, 2))
        self.assertEqual(Hy('5/2').real, Fraction(5, 2))
        self.assertEqual(Hy('-16/5').real, Fraction(-16, 5))
        self.assertEqual(Hy('3.2').real, Fraction(16, 5))

    def test_default_imag_is_zero_and_rank_1(self):
        h = Hy('5/2')
        self.assertEqual(h.imag, Fraction(0))
        self.assertEqual(h.rank, 1)

    def test_two_scalar_args_rank_1(self):
        h = Hy('5/2', '-16/5')
        self.assertEqual(h.real, Fraction(5, 2))
        self.assertEqual(h.imag, Fraction(-16, 5))
        self.assertEqual(h.rank, 1)

    def test_float_uses_decimal_not_binary_fraction(self):
        # 3.2 is not exact in binary; we want the "obvious" 16/5, not the
        # ugly exact-binary-value fraction Fraction(3.2) would give.
        h = Hy(2.5, 3.2)
        self.assertEqual(h, Hy('5/2', '16/5'))
        self.assertEqual(h.imag, Fraction(16, 5))

    def test_quaternion_from_two_complex_hys(self):
        h1 = Hy('1', '2')
        h2 = Hy('3', '4')
        q = Hy(h1, h2)
        self.assertEqual(q.rank, 2)
        self.assertEqual(q.real, h1)
        self.assertEqual(q.imag, h2)

    def test_octonion_from_two_quaternion_hys(self):
        h3 = Hy(Hy(1, 2), Hy(3, 4))
        h4 = Hy(Hy(5, 6), Hy(7, 8))
        o = Hy(h3, h4)
        self.assertEqual(o.rank, 3)
        self.assertEqual(o.dimension, 8)

    def test_rank_4_sedenion(self):
        h5 = Hy(Hy(Hy(1, 0), Hy(0, 0)), Hy(Hy(0, 0), Hy(0, 0)))
        h6 = Hy(Hy(Hy(0, 0), Hy(0, 0)), Hy(Hy(0, 0), Hy(1, 0)))
        s = Hy(h5, h6)
        self.assertEqual(s.rank, 4)
        self.assertEqual(s.dimension, 16)

    def test_single_hy_argument_is_a_copy_not_a_promotion(self):
        q = Hy(Hy(1, 2), Hy(3, 4))
        copy = Hy(q)
        self.assertEqual(copy.rank, q.rank)
        self.assertEqual(copy, q)
        self.assertIsNot(copy, q)

    def test_full_expression_string_single_arg_is_not_promoted(self):
        z = Hy("(5/2+16/5j)")
        self.assertEqual(z.rank, 1)
        self.assertEqual(z, Hy('5/2', '16/5'))

    def test_mismatched_component_ranks_are_promoted_to_match(self):
        # bare Fraction paired with a rank-1 Hy: the Fraction gets
        # embedded as a rank-1 "real" value before pairing.
        h = Hy(3, Hy('1/2', '1/3'))
        self.assertEqual(h.rank, 2)
        self.assertEqual(h.real, Hy('3', '0'))
        self.assertEqual(h.imag, Hy('1/2', '1/3'))

    def test_complex_input_type_is_accepted(self):
        h = Hy(complex(2.5, -3.2))
        self.assertEqual(h, Hy('5/2', '-16/5'))

    def test_bool_input_treated_as_0_or_1(self):
        self.assertEqual(Hy(True).real, Fraction(1))
        self.assertEqual(Hy(False).real, Fraction(0))

    def test_invalid_component_type_raises(self):
        with self.assertRaises(TypeError):
            Hy(object())

    def test_invalid_string_raises(self):
        with self.assertRaises(ValueError):
            Hy("not-a-number")

    def test_real_and_imag_always_same_shape(self):
        for rank in RANKS[1:]:
            rng = random.Random(rank)
            h = rand_value(rank, rng)
            self.assertEqual(_rank(h.real), _rank(h.imag))


# ---------------------------------------------------------------------------
# Accessors
# ---------------------------------------------------------------------------

class TestAccessors(unittest.TestCase):

    def test_real_imag_properties(self):
        h = Hy('5/2', '-16/5')
        self.assertEqual(h.real, Fraction(5, 2))
        self.assertEqual(h.imag, Fraction(-16, 5))

    def test_rank_and_dimension_for_all_ranks(self):
        rng = random.Random(0)
        for rank in RANKS[1:]:
            h = rand_value(rank, rng)
            self.assertEqual(h.rank, rank)
            self.assertEqual(h.dimension, 2 ** rank)

    def test_components_flattens_in_cayley_dickson_order(self):
        q = Hy(Hy(1, 2), Hy(3, 4))
        self.assertEqual(q.components(), (1, 2, 3, 4))
        o = Hy(q, Hy(Hy(5, 6), Hy(7, 8)))
        self.assertEqual(o.components(), (1, 2, 3, 4, 5, 6, 7, 8))

    def test_conjugate_negates_imag_only(self):
        z = Hy('5/2', '-16/5')
        self.assertEqual(z.conjugate(), Hy('5/2', '16/5'))
        q = Hy(Hy(1, 2), Hy(3, 4))
        self.assertEqual(q.conjugate(), Hy(Hy(1, -2), Hy(-3, -4)))

    def test_conjugate_is_involution_all_ranks(self):
        rng = random.Random(1)
        for rank in RANKS[1:]:
            h = rand_value(rank, rng)
            self.assertEqual(h.conjugate().conjugate(), h)

    def test_norm_matches_sum_of_squares(self):
        q = Hy(Hy(1, 2), Hy(3, 4))
        self.assertEqual(q.norm(), 1 * 1 + 2 * 2 + 3 * 3 + 4 * 4)

    def test_inverse_matches_conjugate_over_norm(self):
        h = Hy('1/2', '1/3')
        self.assertEqual(h.inverse(), h.conjugate() * Fraction(1, h.norm()))

    def test_inverse_of_zero_raises(self):
        with self.assertRaises(ZeroDivisionError):
            Hy(0, 0).inverse()

    def test_is_zero(self):
        self.assertTrue(Hy(0).is_zero())
        self.assertTrue(Hy(Hy(0, 0), Hy(0, 0)).is_zero())
        self.assertFalse(Hy('1/7').is_zero())


# ---------------------------------------------------------------------------
# Arithmetic on fixed, hand-checked examples
# ---------------------------------------------------------------------------

class TestArithmeticFixed(unittest.TestCase):

    def test_complex_matches_python_complex(self):
        a = Hy('1/2', '3')
        b = Hy('5', '-2')
        pa, pb = complex(0.5, 3), complex(5, -2)
        self.assertEqual(complex(a + b), pa + pb)
        self.assertEqual(complex(a - b), pa - pb)
        self.assertEqual(complex(a * b), pa * pb)
        qa_over_qb = a / b
        pc = pa / pb
        self.assertAlmostEqual(complex(qa_over_qb).real, pc.real)
        self.assertAlmostEqual(complex(qa_over_qb).imag, pc.imag)

    def test_quaternion_units_multiplication_table(self):
        one = Hy(1)
        i = Hy(Hy(0, 1), Hy(0, 0))
        j = Hy(Hy(0, 0), Hy(1, 0))
        k = Hy(Hy(0, 0), Hy(0, 1))
        self.assertEqual(i * i, -one)
        self.assertEqual(j * j, -one)
        self.assertEqual(k * k, -one)
        self.assertEqual(i * j, k)
        self.assertEqual(j * k, i)
        self.assertEqual(k * i, j)
        self.assertEqual(j * i, -k)
        self.assertEqual(k * j, -i)
        self.assertEqual(i * k, -j)

    def test_scalar_mixed_arithmetic(self):
        h = Hy('1', '2')
        self.assertEqual(h + 3, Hy('4', '2'))
        self.assertEqual(3 + h, Hy('4', '2'))
        self.assertEqual(2 * h, Hy('2', '4'))
        self.assertEqual(h * 2, Hy('2', '4'))
        self.assertEqual(h - 1, Hy('0', '2'))
        self.assertEqual(1 - h, Hy('0', '-2'))
        self.assertEqual(h / 2, Hy('1/2', '1'))

    def test_division_by_python_complex(self):
        h = Hy('4', '2')
        self.assertEqual(h / complex(2, 0), Hy('2', '1'))

    def test_neg_and_pos(self):
        h = Hy('1/2', '-3')
        self.assertEqual(-h, Hy('-1/2', '3'))
        self.assertIs(+h, h)

    def test_division_by_zero_raises(self):
        with self.assertRaises(ZeroDivisionError):
            Hy(1, 2) / Hy(0, 0)


# ---------------------------------------------------------------------------
# Arithmetic algebraic-law fuzz tests, ranks 0-4
# ---------------------------------------------------------------------------

class TestArithmeticFuzz(unittest.TestCase):

    def test_addition_commutative_and_associative(self):
        for rank in RANKS:
            rng = random.Random(100 + rank)
            for _ in range(FUZZ_TRIALS):
                a, b, c = (rand_value(rank, rng) for _ in range(3))
                self.assertEqual(as_hy(add(a, b)), as_hy(add(b, a)))
                self.assertEqual(
                    as_hy(add(add(a, b), c)), as_hy(add(a, add(b, c)))
                )

    def test_additive_inverse(self):
        for rank in RANKS:
            rng = random.Random(200 + rank)
            for _ in range(FUZZ_TRIALS):
                a = rand_value(rank, rng)
                zero = _embed(Fraction(0), rank)
                self.assertEqual(as_hy(add(a, neg(a))), as_hy(zero))

    def test_subtraction_consistent_with_add_neg(self):
        for rank in RANKS:
            rng = random.Random(300 + rank)
            for _ in range(FUZZ_TRIALS):
                a, b = rand_value(rank, rng), rand_value(rank, rng)
                self.assertEqual(as_hy(sub(a, b)), as_hy(add(a, neg(b))))

    def test_multiplicative_distributivity(self):
        # a(b+c) == ab + ac  and  (a+b)c == ac + bc  (holds at every rank
        # -- distributivity does not require associativity/commutativity)
        for rank in RANKS:
            rng = random.Random(400 + rank)
            for _ in range(FUZZ_TRIALS):
                a, b, c = (rand_value(rank, rng) for _ in range(3))
                lhs = mul(a, add(b, c))
                rhs = add(mul(a, b), mul(a, c))
                self.assertEqual(as_hy(lhs), as_hy(rhs))
                lhs2 = mul(add(a, b), c)
                rhs2 = add(mul(a, c), mul(b, c))
                self.assertEqual(as_hy(lhs2), as_hy(rhs2))

    def test_conjugate_of_product_reverses_order(self):
        # conj(xy) == conj(y) conj(x)  -- true for every Cayley-Dickson rank
        for rank in RANKS:
            rng = random.Random(500 + rank)
            for _ in range(FUZZ_TRIALS):
                x, y = rand_value(rank, rng), rand_value(rank, rng)
                lhs = conj(mul(x, y))
                rhs = mul(conj(y), conj(x))
                self.assertEqual(as_hy(lhs), as_hy(rhs))

    def test_x_times_conjugate_x_equals_norm(self):
        # x * conj(x) is always the (scalar) norm, at every rank.
        for rank in RANKS:
            rng = random.Random(600 + rank)
            for _ in range(FUZZ_TRIALS):
                x = rand_value(rank, rng)
                n = abs2(x)
                self.assertEqual(as_hy(mul(x, conj(x))), as_hy(_embed(n, rank)))

    def test_norm_multiplicative_up_to_octonions(self):
        # |xy| = |x| |y| holds for rank <= 3 (reals, complex, quaternions,
        # octonions are all composition algebras); we do NOT assert this
        # for rank 4 (sedenions), where it can fail.
        for rank in (0, 1, 2, 3):
            rng = random.Random(700 + rank)
            for _ in range(FUZZ_TRIALS):
                x, y = rand_value(rank, rng), rand_value(rank, rng)
                self.assertEqual(abs2(mul(x, y)), abs2(x) * abs2(y))

    def test_inverse_both_sides_up_to_rank_4(self):
        # x * x^-1 == x^-1 * x == 1 for any nonzero x, at every rank,
        # since it follows directly from x * conj(x) == norm(x) (a scalar).
        for rank in RANKS:
            rng = random.Random(800 + rank)
            one = _embed(Fraction(1), rank)
            trials = 0
            while trials < FUZZ_TRIALS:
                x = rand_value(rank, rng)
                if _is_zero_val(x):
                    continue
                trials += 1
                xi = inverse(x)
                self.assertEqual(as_hy(mul(x, xi)), as_hy(one))
                self.assertEqual(as_hy(mul(xi, x)), as_hy(one))

    def test_division_is_multiplication_by_inverse(self):
        for rank in RANKS:
            rng = random.Random(900 + rank)
            trials = 0
            while trials < FUZZ_TRIALS:
                x, y = rand_value(rank, rng), rand_value(rank, rng)
                if _is_zero_val(y):
                    continue
                trials += 1
                self.assertEqual(as_hy(div(x, y)), as_hy(mul(x, inverse(y))))

    def test_octonions_are_alternative(self):
        # An alternative algebra satisfies x(xy) == (xx)y and (yx)x == y(xx)
        # even where full associativity fails. True for octonions (rank 3).
        rank = 3
        rng = random.Random(1000)
        for _ in range(FUZZ_TRIALS):
            x, y = rand_value(rank, rng), rand_value(rank, rng)
            self.assertEqual(
                as_hy(mul(x, mul(x, y))), as_hy(mul(mul(x, x), y))
            )
            self.assertEqual(
                as_hy(mul(mul(y, x), x)), as_hy(mul(y, mul(x, x)))
            )

    def test_octonions_generally_non_associative(self):
        rng = random.Random(1001)
        found_non_associative = False
        for _ in range(FUZZ_TRIALS):
            a, b, c = (rand_value(3, rng) for _ in range(3))
            if as_hy(mul(mul(a, b), c)) != as_hy(mul(a, mul(b, c))):
                found_non_associative = True
                break
        self.assertTrue(found_non_associative)

    def test_sedenions_have_zero_divisors(self):
        # Rank 4 (dimension 16) is the first rank where the algebra stops
        # being a composition/division algebra: there exist nonzero x, y
        # with x*y == 0. We search a modest space of basis combinations
        # rather than hard-coding a published example, since the exact
        # sign/order convention is implementation-specific.
        n = 16
        basis = [_unflatten(
            [Fraction(1) if k == idx else Fraction(0) for k in range(n)], 4
        ) for idx in range(n)]
        zero = _embed(Fraction(0), 4)
        found = None
        for i in range(1, n):
            for j in range(i + 1, n):
                u = add(basis[i], basis[j])
                for k in range(1, n):
                    for l in range(k + 1, n):
                        if {k, l} == {i, j}:
                            continue
                        v = sub(basis[k], basis[l])
                        if mul(u, v) == zero:
                            found = (u, v)
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break
        self.assertIsNotNone(
            found, "expected to find at least one sedenion zero-divisor pair"
        )
        u, v = found
        self.assertFalse(_is_zero_val(u))
        self.assertFalse(_is_zero_val(v))


# ---------------------------------------------------------------------------
# Comparison and hashing
# ---------------------------------------------------------------------------

class TestComparisonAndHash(unittest.TestCase):

    def test_eq_and_ne(self):
        self.assertEqual(Hy('1/2', '1/3'), Hy('1/2', '1/3'))
        self.assertNotEqual(Hy('1/2', '1/3'), Hy('1/2', '1/4'))
        self.assertTrue(Hy(1) != Hy(2))
        self.assertFalse(Hy(1) != Hy(1))

    def test_eq_across_ranks_via_zero_padding(self):
        self.assertEqual(Hy('3'), Hy(Hy('3', '0'), Hy('0', '0')))
        self.assertEqual(
            Hy(Hy('3', '0'), Hy('0', '0')), Hy(Hy(Hy('3', '0'), Hy('0', '0')), Hy(0))
        )

    def test_eq_with_python_numeric_types(self):
        self.assertEqual(Hy(3), 3)
        self.assertEqual(Hy(3), 3.0)
        self.assertEqual(Hy('5/2', '-16/5'), complex(2.5, -3.2))
        self.assertEqual(Hy(3), Fraction(3))

    def test_eq_with_incompatible_type_returns_not_implemented_and_false(self):
        self.assertFalse(Hy(1) == "not a number")
        self.assertNotEqual(Hy(1), object())

    def test_hash_consistent_with_eq_across_ranks(self):
        pairs = [
            (Hy('3'), Hy(Hy('3', '0'), Hy('0', '0'))),
            (Hy(0), Hy(Hy(0, 0), Hy(0, 0))),
            (Hy('7/2'), Hy(Hy(Hy('7/2', '0'), Hy('0', '0')), Hy(Hy(0, 0), Hy(0, 0)))),
        ]
        for a, b in pairs:
            self.assertEqual(a, b)
            self.assertEqual(hash(a), hash(b))

    def test_hash_stable_and_usable_in_sets(self):
        s = {Hy('1/2', '1/3'), Hy('1/2', '1/3'), Hy(1, 1)}
        self.assertEqual(len(s), 2)


# ---------------------------------------------------------------------------
# Sequence-like protocol
# ---------------------------------------------------------------------------

class TestSequenceProtocol(unittest.TestCase):

    def test_iter_yields_real_then_imag(self):
        h = Hy('1/2', '-3')
        self.assertEqual(list(h), [Fraction(1, 2), Fraction(-3)])

    def test_getitem(self):
        h = Hy('1/2', '-3')
        self.assertEqual(h[0], Fraction(1, 2))
        self.assertEqual(h[1], Fraction(-3))
        self.assertEqual(h[:], (Fraction(1, 2), Fraction(-3)))

    def test_len_is_always_two(self):
        for rank in RANKS[1:]:
            rng = random.Random(rank)
            self.assertEqual(len(rand_value(rank, rng)), 2)

    def test_bool(self):
        self.assertFalse(bool(Hy(0, 0)))
        self.assertTrue(bool(Hy(0, 1)))
        self.assertTrue(bool(Hy(1)))


# ---------------------------------------------------------------------------
# Numeric conversions
# ---------------------------------------------------------------------------

class TestConversions(unittest.TestCase):

    def test_abs_matches_sqrt_norm(self):
        h = Hy('3', '4')
        self.assertAlmostEqual(abs(h), 5.0)
        q = Hy(Hy(1, 2), Hy(3, 4))
        self.assertAlmostEqual(abs(q), math.sqrt(1 + 4 + 9 + 16))

    def test_complex_conversion(self):
        h = Hy('5/2', '-16/5')
        self.assertEqual(complex(h), complex(2.5, -3.2))

    def test_complex_conversion_rejects_higher_rank_nonzero(self):
        q = Hy(Hy(1, 2), Hy(3, 4))
        with self.assertRaises(ValueError):
            complex(q)

    def test_complex_conversion_allows_higher_rank_if_extra_parts_are_zero(self):
        q = Hy(Hy(1, 2), Hy(0, 0))
        self.assertEqual(complex(q), complex(1, 2))

    def test_pow_zero_one_and_negative(self):
        h = Hy('2', '0')
        self.assertEqual(h ** 0, Hy(1))
        self.assertEqual(h ** 1, h)
        self.assertEqual(h ** 2, Hy('4', '0'))
        self.assertEqual(h ** -1, h.inverse())

    def test_pow_on_quaternion_unit(self):
        i = Hy(Hy(0, 1), Hy(0, 0))
        self.assertEqual(i ** 2, Hy(-1))
        self.assertEqual(i ** 4, Hy(1))


# ---------------------------------------------------------------------------
# String forms
# ---------------------------------------------------------------------------

class TestStringForms(unittest.TestCase):

    def test_str_rank_1(self):
        self.assertEqual(str(Hy('5/2', '-16/5')), '(5/2-16/5j)')
        self.assertEqual(str(Hy('5/2', '16/5')), '(5/2+16/5j)')
        self.assertEqual(str(Hy(0, 0)), '(0)')
        self.assertEqual(str(Hy(0, 1)), '(j)')
        self.assertEqual(str(Hy(0, -1)), '(-j)')

    def test_str_rank_2(self):
        self.assertEqual(str(Hy(Hy(1, 2), Hy(3, 4))), '(1+2i+3j+4k)')
        self.assertEqual(str(Hy(Hy(0, 1), Hy(0, 0))), '(i)')
        self.assertEqual(str(Hy(Hy(0, 0), Hy(0, -1))), '(-k)')

    def test_str_rank_3_uses_octonion_labels(self):
        o = Hy(Hy(Hy(1, 0), Hy(0, 0)), Hy(Hy(0, 1), Hy(0, 0)))
        self.assertEqual(str(o), '(1+iL)')

    def test_str_rank_4_uses_e_labels_up_to_e15(self):
        coeffs = [Fraction(0)] * 16
        coeffs[0] = Fraction(1)
        coeffs[15] = Fraction(-2)
        s = _unflatten(coeffs, 4)
        self.assertEqual(str(s), '(1-2e15)')

    def test_repr_round_trips(self):
        for h in (
            Hy('5/2', '-16/5'),
            Hy(Hy(1, 2), Hy(3, 4)),
            Hy(Hy(Hy(1, 2), Hy(3, 4)), Hy(Hy(5, 6), Hy(7, 8))),
        ):
            self.assertEqual(eval(repr(h), {"Hy": Hy}), h)

    def test_repr_uses_quoted_fraction_strings_at_leaves(self):
        self.assertEqual(repr(Hy('5/2', '-16/5')), "Hy('5/2', '-16/5')")


# ---------------------------------------------------------------------------
# Hy.units() / .is_unit()
# ---------------------------------------------------------------------------

class TestUnits(unittest.TestCase):

    def test_units_rank_0_returns_plain_fractions(self):
        u = Hy.units(0)
        self.assertEqual(u, {"1": Fraction(1), "-1": Fraction(-1)})
        for v in u.values():
            self.assertIsInstance(v, Fraction)
            self.assertNotIsInstance(v, Hy)

    def test_units_rank_1(self):
        self.assertEqual(
            Hy.units(1),
            {"1": Hy(1, 0), "-1": Hy(-1, 0), "j": Hy(0, 1), "-j": Hy(0, -1)},
        )

    def test_units_rank_1_key_order(self):
        # 1, -1, j, -j -- grouped +/- per basis element, in basis order.
        self.assertEqual(list(Hy.units(1).keys()), ["1", "-1", "j", "-j"])

    def test_units_rank_2(self):
        u = Hy.units(2)
        self.assertEqual(
            set(u.keys()), {"1", "-1", "i", "-i", "j", "-j", "k", "-k"}
        )
        self.assertEqual(u["i"], Hy(Hy(0, 1), Hy(0, 0)))
        self.assertEqual(u["-k"], Hy(Hy(0, 0), Hy(0, -1)))

    def test_units_rank_3_uses_octonion_labels(self):
        u = Hy.units(3)
        self.assertEqual(len(u), 16)  # 2 * 2**3
        self.assertIn("iL", u)
        self.assertIn("-kL", u)
        self.assertEqual(str(u["iL"]), "(iL)")

    def test_units_count_matches_2_times_dimension(self):
        for rank in RANKS:
            self.assertEqual(len(Hy.units(rank)), 2 * (2 ** rank))

    def test_units_values_all_satisfy_is_unit(self):
        for rank in RANKS[1:]:
            for v in Hy.units(rank).values():
                self.assertTrue(v.is_unit())

    def test_units_string_keys_match_str_of_values(self):
        for rank in RANKS[1:]:
            for key, val in Hy.units(rank).items():
                self.assertEqual(str(val), f"({key})")

    def test_units_rejects_bad_rank(self):
        for bad_rank in (-1, 2.5, "2", True, False):
            with self.assertRaises(ValueError):
                Hy.units(bad_rank)

    def test_is_unit_true_for_basis_units(self):
        self.assertTrue(Hy(1, 0).is_unit())
        self.assertTrue(Hy(-1, 0).is_unit())
        self.assertTrue(Hy(0, 1).is_unit())
        self.assertTrue(Hy(0, -1).is_unit())
        self.assertTrue(Hy(Hy(0, 1), Hy(0, 0)).is_unit())  # quaternion i

    def test_is_unit_false_for_zero(self):
        self.assertFalse(Hy(0, 0).is_unit())

    def test_is_unit_false_for_non_unit_values(self):
        self.assertFalse(Hy(1, 1).is_unit())
        self.assertFalse(Hy(2, 0).is_unit())
        self.assertFalse(Hy('1/2', 0).is_unit())

    def test_is_unit_consistent_with_units_membership(self):
        for rank in RANKS[1:]:
            rng = random.Random(6000 + rank)
            unit_values = list(Hy.units(rank).values())
            for _ in range(10):
                h = rand_value(rank, rng)
                self.assertEqual(h.is_unit(), h in unit_values)


# ---------------------------------------------------------------------------
# Hy.latex()
# ---------------------------------------------------------------------------

class TestLatex(unittest.TestCase):

    def test_latex_integer_coefficients_rank_1(self):
        self.assertEqual(Hy(3, -2).latex(), "3-2j")
        self.assertEqual(Hy(0, 1).latex(), "j")
        self.assertEqual(Hy(0, -1).latex(), "-j")
        self.assertEqual(Hy(0, 0).latex(), "0")

    def test_latex_default_vinculum_is_horizontal(self):
        self.assertEqual(Hy('5/2', '-16/5').latex(), r'\frac{5}{2}-\frac{16}{5}j')

    def test_latex_diagonal_vinculum(self):
        self.assertEqual(
            Hy('5/2', '-16/5').latex(vinculum='diagonal'), '5/2-16/5j'
        )

    def test_latex_rejects_bad_vinculum(self):
        with self.assertRaises(ValueError):
            Hy(1, 2).latex(vinculum='vertical')

    def test_latex_rank_2_uses_ijk_unchanged(self):
        self.assertEqual(Hy(Hy(1, 2), Hy(3, 4)).latex(), "1+2i+3j+4k")

    def test_latex_rank_3_uses_octonion_labels_unsubscripted(self):
        o = Hy(Hy(Hy(1, 0), Hy(0, 0)), Hy(Hy(0, 1), Hy(0, 0)))
        self.assertEqual(o.latex(), "1+iL")

    def test_latex_rank_4_subscripts_multidigit_e_labels(self):
        coeffs = [Fraction(0)] * 16
        coeffs[0] = Fraction(1)
        coeffs[15] = Fraction(-2)
        h = _unflatten(coeffs, 4)
        self.assertEqual(h.latex(), "1-2e_{15}")

    def test_latex_mode_plain_default(self):
        self.assertEqual(Hy(1, 2).latex(mode='plain'), Hy(1, 2).latex())

    def test_latex_mode_inline_wraps_in_dollar_signs(self):
        self.assertEqual(Hy('1/2', 0).latex(mode='inline'), r'$\frac{1}{2}$')

    def test_latex_mode_display_wraps_in_display_delimiters(self):
        self.assertEqual(Hy('1/2', 0).latex(mode='display'), r'\[\frac{1}{2}\]')

    def test_latex_rejects_bad_mode(self):
        with self.assertRaises(ValueError):
            Hy(1, 2).latex(mode='bogus')

    def test_latex_unit_coefficient_omits_the_1(self):
        # a coefficient of exactly 1 (or -1) drops the leading "1", just
        # like str() does -- e.g. "e_{5}", not "1e_{5}".
        self.assertEqual(Hy(0, 1).latex(), "j")
        self.assertEqual(Hy(0, -1).latex(), "-j")

    def test_latex_unit_label_helper(self):
        self.assertEqual(_latex_unit_label(""), "")
        self.assertEqual(_latex_unit_label("j"), "j")
        self.assertEqual(_latex_unit_label("i"), "i")
        self.assertEqual(_latex_unit_label("k"), "k")
        self.assertEqual(_latex_unit_label("L"), "L")
        self.assertEqual(_latex_unit_label("iL"), "iL")
        self.assertEqual(_latex_unit_label("jL"), "jL")
        self.assertEqual(_latex_unit_label("kL"), "kL")
        self.assertEqual(_latex_unit_label("e1"), "e_{1}")
        self.assertEqual(_latex_unit_label("e15"), "e_{15}")

    def test_format_fraction_latex_helper(self):
        self.assertEqual(_format_fraction_latex(Fraction(3), "horizontal"), "3")
        self.assertEqual(
            _format_fraction_latex(Fraction(5, 2), "horizontal"), r"\frac{5}{2}"
        )
        self.assertEqual(_format_fraction_latex(Fraction(5, 2), "diagonal"), "5/2")

    def test_latex_random_values_all_ranks_smoke_test(self):
        # No fixed expected string here -- just confirm it runs cleanly
        # and produces a non-empty string for every rank/vinculum combo.
        for rank in RANKS[1:]:
            h = Hy.random(rank, seed=7000 + rank)
            for vinculum in ("horizontal", "diagonal"):
                s = h.latex(vinculum=vinculum)
                self.assertIsInstance(s, str)
                self.assertTrue(s)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

class TestParsing(unittest.TestCase):

    def test_parse_rank_1(self):
        self.assertEqual(Hy.parse('(5/2+16/5j)'), Hy('5/2', '16/5'))
        self.assertEqual(Hy.parse('5/2-16/5j'), Hy('5/2', '-16/5'))
        self.assertEqual(Hy.parse('j'), Hy(0, 1))
        self.assertEqual(Hy.parse('-j'), Hy(0, -1))
        self.assertEqual(Hy.parse('3'), Hy(3))

    def test_parse_rank_2(self):
        self.assertEqual(Hy.parse('1+2i+3j+4k'), Hy(Hy(1, 2), Hy(3, 4)))
        self.assertEqual(Hy.parse('-k'), Hy(Hy(0, 0), Hy(0, -1)))
        self.assertEqual(Hy.parse('i+k'), Hy(Hy(0, 1), Hy(0, 1)))

    def test_parse_rank_3_and_4(self):
        self.assertEqual(Hy.parse('1+iL'), Hy(Hy(Hy(1, 0), Hy(0, 0)), Hy(Hy(0, 1), Hy(0, 0))))
        self.assertEqual(Hy.parse('1-2e15').rank, 4)

    def test_parse_rank_3_octonion_labels(self):
        self.assertEqual(
            Hy.parse('L'), Hy(Hy(Hy(0, 0), Hy(0, 0)), Hy(Hy(1, 0), Hy(0, 0)))
        )
        self.assertEqual(
            Hy.parse('-jL'), Hy(Hy(Hy(0, 0), Hy(0, 0)), Hy(Hy(0, 0), Hy(-1, 0)))
        )
        self.assertEqual(Hy.parse('1-k+2kL').rank, 3)

    def test_parse_from_string_alias(self):
        self.assertEqual(Hy.from_string('5/2+16/5j'), Hy('5/2', '16/5'))

    def test_parse_round_trips_with_str_for_random_values(self):
        for rank in (1, 2, 3, 4):
            rng = random.Random(2000 + rank)
            for _ in range(10):
                h = rand_value(rank, rng)
                self.assertEqual(Hy.parse(str(h)), h)

    def test_parse_invalid_term_raises(self):
        with self.assertRaises(ValueError):
            Hy.parse('1+2x')

    def test_parse_empty_raises(self):
        with self.assertRaises(ValueError):
            Hy.parse('')
        with self.assertRaises(ValueError):
            Hy.parse('()')

    def test_parse_inconsistent_units_raises(self):
        # 'i' is a rank <= 3 label; 'e5' only appears from rank 4 up
        # (labels e1..e15), so the two are inconsistent together.
        with self.assertRaises(ValueError):
            Hy.parse('1+i+e5')

    def test_parse_inconsistent_octonion_and_e_label_raises(self):
        # 'L' is an octonion-only (rank 3) label; 'e5' only appears from
        # rank 4 up, so the two are inconsistent together.
        with self.assertRaises(ValueError):
            Hy.parse('1+L+e5')

    def test_module_level_parse_function(self):
        self.assertEqual(_parse('5/2+16/5j'), Hy('5/2', '16/5'))


# ---------------------------------------------------------------------------
# Module-level functions, called directly (including on raw Fractions,
# i.e. "rank 0")
# ---------------------------------------------------------------------------

class TestModuleFunctions(unittest.TestCase):

    def test_add_sub_neg_on_raw_fractions(self):
        a, b = Fraction(1, 2), Fraction(1, 3)
        self.assertEqual(add(a, b), Fraction(5, 6))
        self.assertEqual(sub(a, b), Fraction(1, 6))
        self.assertEqual(neg(a), Fraction(-1, 2))

    def test_mul_div_on_raw_fractions(self):
        a, b = Fraction(2, 3), Fraction(3, 4)
        self.assertEqual(mul(a, b), Fraction(1, 2))
        self.assertEqual(div(a, b), a / b)

    def test_conj_of_raw_fraction_is_itself(self):
        self.assertEqual(conj(Fraction(5, 2)), Fraction(5, 2))

    def test_abs2_of_raw_fraction_is_its_square(self):
        self.assertEqual(abs2(Fraction(-3, 2)), Fraction(9, 4))

    def test_inverse_of_raw_fraction(self):
        self.assertEqual(inverse(Fraction(4)), Fraction(1, 4))

    def test_functions_mix_raw_fraction_and_hy_operands(self):
        # a Fraction (rank 0) combined with a rank-1 Hy should be promoted
        # automatically by every one of these functions.
        f = Fraction(3)
        h = Hy('1', '2')
        self.assertEqual(add(f, h), Hy('4', '2'))
        self.assertEqual(mul(f, h), Hy('3', '6'))
        self.assertEqual(sub(h, f), Hy('-2', '2'))

    def test_rank_helper(self):
        self.assertEqual(_rank(Fraction(1)), 0)
        self.assertEqual(_rank(Hy(1)), 1)
        self.assertEqual(_rank(Hy(Hy(1, 0), Hy(0, 0))), 2)

    def test_embed_raises_when_reducing_rank(self):
        with self.assertRaises(ValueError):
            _embed(Hy(1, 2), 0)

    def test_embed_is_identity_at_same_rank(self):
        h = Hy('1/2', '1/3')
        self.assertIs(_embed(h, 1), h)

    def test_flatten_unflatten_round_trip(self):
        for rank in RANKS[1:]:
            rng = random.Random(3000 + rank)
            h = rand_value(rank, rng)
            coeffs = _flatten(h)
            self.assertEqual(len(coeffs), 2 ** rank)
            self.assertEqual(_unflatten(coeffs, rank), h)


# ---------------------------------------------------------------------------
# Hy.random()
# ---------------------------------------------------------------------------

class TestRandom(unittest.TestCase):

    def test_random_returns_hy_of_requested_rank(self):
        for rank in RANKS[1:]:
            h = Hy.random(rank, seed=rank)
            self.assertIsInstance(h, Hy)
            self.assertEqual(h.rank, rank)
            self.assertEqual(h.dimension, 2 ** rank)

    def test_random_rejects_non_positive_or_non_int_rank(self):
        for bad_rank in (0, -1, 2.5, "2", True, False):
            with self.assertRaises(ValueError):
                Hy.random(bad_rank)

    def test_random_seed_kw_is_reproducible(self):
        a = Hy.random(3, seed=12345)
        b = Hy.random(3, seed=12345)
        self.assertEqual(a, b)

    def test_random_different_seeds_differ_with_overwhelming_probability(self):
        a = Hy.random(3, seed=1)
        b = Hy.random(3, seed=2)
        self.assertNotEqual(a, b)

    def test_random_rng_kw_is_reproducible_and_shareable(self):
        a = Hy.random(2, rng=random.Random(999))
        b = Hy.random(2, rng=random.Random(999))
        self.assertEqual(a, b)

    def test_random_rejects_both_rng_and_seed(self):
        with self.assertRaises(ValueError):
            Hy.random(1, rng=random.Random(0), seed=0)

    def test_hy_seed_makes_default_rng_reproducible(self):
        Hy.seed(2024)
        a = Hy.random(2)
        Hy.seed(2024)
        b = Hy.random(2)
        self.assertEqual(a, b)

    def test_random_respects_lo_hi_dmax_bounds(self):
        rng = random.Random(7)
        for _ in range(20):
            h = Hy.random(3, lo=-2, hi=2, dmax=3, rng=rng)
            for c in h.components():
                self.assertGreaterEqual(c.numerator, -2)
                self.assertLessEqual(c.numerator, 2)
                self.assertGreaterEqual(c.denominator, 1)
                self.assertLessEqual(c.denominator, 3)

    def test_random_values_are_exact_fractions(self):
        h = Hy.random(2, seed=0)
        for c in h.components():
            self.assertIsInstance(c, Fraction)


# ---------------------------------------------------------------------------
# Hy.from_array() / Hy.to_array()
# ---------------------------------------------------------------------------

class TestArrayConversion(unittest.TestCase):

    def test_to_array_matches_components_as_fractions(self):
        q = Hy(Hy(1, 2), Hy(3, 4))
        self.assertEqual(q.to_array(), [Fraction(1), Fraction(2), Fraction(3), Fraction(4)])
        self.assertEqual(q.to_array(), list(q.components()))

    def test_to_array_as_str(self):
        h = Hy('5/2', '-16/5')
        self.assertEqual(h.to_array(as_str=True), ['5/2', '-16/5'])

    def test_from_array_numbers_only(self):
        q = Hy.from_array([1, 2, 3, 4])
        self.assertEqual(q, Hy(Hy(1, 2), Hy(3, 4)))
        self.assertEqual(q.rank, 2)

    def test_from_array_strings_only(self):
        h = Hy.from_array(['5/2', '-16/5'])
        self.assertEqual(h, Hy('5/2', '-16/5'))

    def test_from_array_mixed_numbers_and_strings(self):
        h = Hy.from_array([1, '2/3', 3.5, '-1/4'])
        self.assertEqual(h, Hy(Hy('1', '2/3'), Hy('7/2', '-1/4')))

    def test_from_array_accepts_tuple_and_generator(self):
        self.assertEqual(Hy.from_array((1, 2)), Hy(1, 2))
        self.assertEqual(Hy.from_array(x for x in (1, 2, 3, 4)), Hy(Hy(1, 2), Hy(3, 4)))

    def test_from_array_length_must_be_power_of_two(self):
        for bad_len in (0, 1, 3, 5, 6, 7, 9):
            with self.assertRaises(ValueError):
                Hy.from_array(list(range(bad_len)))

    def test_from_array_rejects_composite_string_element(self):
        # from_array elements must be *plain* fractions, not composite
        # Hy expressions like '1+2j' (that's what Hy.parse is for).
        with self.assertRaises(ValueError):
            Hy.from_array(['1+2j', '3'])

    def test_from_array_rejects_bad_element_type(self):
        with self.assertRaises(TypeError):
            Hy.from_array([object(), 1])

    def test_from_array_octonion_length_8(self):
        vals = list(range(8))
        o = Hy.from_array(vals)
        self.assertEqual(o.rank, 3)
        self.assertEqual(o.components(), tuple(Fraction(v) for v in vals))

    def test_round_trip_to_array_from_array_all_ranks(self):
        for rank in RANKS[1:]:
            rng = random.Random(4000 + rank)
            h = rand_value(rank, rng)
            self.assertEqual(Hy.from_array(h.to_array()), h)
            self.assertEqual(Hy.from_array(h.to_array(as_str=True)), h)

    def test_round_trip_from_array_to_array_random_values(self):
        for rank in RANKS[1:]:
            h = Hy.random(rank, seed=5000 + rank)
            arr = h.to_array(as_str=True)
            self.assertEqual(len(arr), 2 ** rank)
            self.assertEqual(Hy.from_array(arr), h)


# ---------------------------------------------------------------------------
# Immutability
# ---------------------------------------------------------------------------

class TestImmutability(unittest.TestCase):

    def test_cannot_set_real_or_imag(self):
        h = Hy('1', '2')
        with self.assertRaises(AttributeError):
            h.real = Fraction(5)
        with self.assertRaises(AttributeError):
            h._real = Fraction(5)
        with self.assertRaises(AttributeError):
            h.some_new_attr = 1

    def test_cannot_delete_attributes(self):
        h = Hy('1', '2')
        with self.assertRaises(AttributeError):
            del h._real

    def test_no_instance_dict(self):
        h = Hy('1', '2')
        self.assertFalse(hasattr(h, '__dict__'))


# ============================================================================
# Interoperability with sympy, numpy-quaternion and quaternionic
# ============================================================================

INTEROP_TRIALS = 25


def rand_quat(seed):
    """A reproducible random rational quaternion (rank 2)."""
    return Hy.random(2, rng=random.Random(seed))


def as_floats(h):
    return [float(c) for c in h.to_array()]


class _FakeNumpyQuaternion:
    """Stands in for a numpy.quaternion scalar (only .w/.x/.y/.z, .ndim)."""
    ndim = 0

    def __init__(self, w, x, y, z):
        self.w, self.x, self.y, self.z = w, x, y, z


class TestFloatToFraction(unittest.TestCase):

    def test_default_uses_shortest_decimal_repr(self):
        self.assertEqual(_float_to_fraction(0.1, False, None), Fraction(1, 10))
        self.assertEqual(_float_to_fraction(-2.5, False, None), Fraction(-5, 2))
        self.assertEqual(_float_to_fraction(3, False, None), Fraction(3))

    def test_exact_uses_binary_value(self):
        self.assertEqual(_float_to_fraction(0.1, True, None), Fraction(0.1))
        self.assertNotEqual(_float_to_fraction(0.1, True, None), Fraction(1, 10))
        self.assertEqual(_float_to_fraction(0.375, True, None), Fraction(3, 8))

    def test_max_denominator_recovers_simple_fractions(self):
        self.assertEqual(_float_to_fraction(1 / 3, False, 1000), Fraction(1, 3))
        self.assertEqual(_float_to_fraction(-2 / 7, False, 1000), Fraction(-2, 7))
        self.assertEqual(_float_to_fraction(0.1, False, 5), Fraction(0))   # closest with d <= 5

    def test_round_trip_float_is_always_exact(self):
        # float -> Fraction -> float returns the very same float, in
        # both the default and the exact modes.
        rng = random.Random(1)
        for _ in range(200):
            f = rng.uniform(-1e3, 1e3)
            self.assertEqual(float(_float_to_fraction(f, False, None)), f)
            self.assertEqual(float(_float_to_fraction(f, True, None)), f)

    def test_non_finite_values_rejected(self):
        for bad in (float('nan'), float('inf'), float('-inf')):
            with self.assertRaises(ValueError):
                _float_to_fraction(bad, False, None)

    def test_non_numbers_rejected(self):
        for bad in ('1.5', b'1', None, [1], 1j):
            with self.assertRaises(TypeError):
                _float_to_fraction(bad, False, None)


class _FakeSympyRational:
    """Stands in for a sympy.Rational without needing sympy installed."""
    is_number = is_real = is_finite = is_rational = True

    def __init__(self, p, q=1):
        self.p, self.q = p, q


class _FakeSympyFloat:
    """Stands in for a sympy.Float (is_rational is None, like the real thing)."""
    is_number = is_real = is_finite = True
    is_rational = None

    def __init__(self, value):
        self._value = value

    def __float__(self):
        return float(self._value)


class _FakeSympyBad:
    """Stands in for sympy's nan/oo/I/Symbol -- rejected one way or another."""
    def __init__(self, is_number=True, is_real=False, is_finite=False):
        self.is_number, self.is_real, self.is_finite = is_number, is_real, is_finite


class _FakeSympyQuaternion:
    def __init__(self, a, b, c, d):
        self.a, self.b, self.c, self.d = a, b, c, d


class TestSympyToFraction(unittest.TestCase):
    """_sympy_to_fraction() is exercised via small fakes, so these run with
    no dependency on sympy actually being installed."""

    def test_exact_rational_converts_with_no_rounding(self):
        self.assertEqual(_sympy_to_fraction(_FakeSympyRational(2, 3), False, None),
                         Fraction(2, 3))
        self.assertEqual(_sympy_to_fraction(_FakeSympyRational(-5), False, None),
                         Fraction(-5))

    def test_float_component_uses_float_policy(self):
        f = _FakeSympyFloat(1 / 3)
        self.assertNotEqual(_sympy_to_fraction(f, False, None), Fraction(1, 3))
        self.assertEqual(_sympy_to_fraction(f, False, 1000), Fraction(1, 3))

    def test_non_real_or_non_finite_rejected(self):
        for kwargs in (dict(is_real=False), dict(is_finite=False), dict(is_real=None)):
            with self.assertRaises(ValueError):
                _sympy_to_fraction(_FakeSympyBad(**kwargs), False, None)

    def test_symbolic_rejected(self):
        with self.assertRaises(TypeError):
            _sympy_to_fraction(_FakeSympyBad(is_number=False), False, None)


class TestInteropWithoutPackages(unittest.TestCase):
    """The from_* methods are duck-typed and the to_* methods import their
    package lazily, so a good deal can be checked with no third-party
    package present at all."""

    # ---- from_numpy_quaternion (duck-typed) ------------------------- #
    def test_from_numpy_quaternion_uses_wxyz(self):
        h = Hy.from_numpy_quaternion(_FakeNumpyQuaternion(1, 2, 3, 4))
        self.assertEqual(h, Hy(Hy(1, 2), Hy(3, 4)))
        self.assertEqual(h.rank, 2)

    def test_from_numpy_quaternion_float_policies(self):
        q = _FakeNumpyQuaternion(0.1, 1 / 3, 0.375, -2.0)
        self.assertEqual(Hy.from_numpy_quaternion(q).to_array()[0], Fraction(1, 10))
        self.assertEqual(Hy.from_numpy_quaternion(q, max_denominator=100).to_array()[1],
                         Fraction(1, 3))
        self.assertEqual(Hy.from_numpy_quaternion(q, exact=True).to_array()[0],
                         Fraction(0.1))

    def test_from_numpy_quaternion_rejects_arrays_and_junk(self):
        arr = _FakeNumpyQuaternion(1, 2, 3, 4)
        arr.ndim, arr.shape = 1, (3,)
        with self.assertRaises(ValueError):
            Hy.from_numpy_quaternion(arr)
        with self.assertRaises(TypeError):
            Hy.from_numpy_quaternion(object())
        with self.assertRaises(TypeError):
            Hy.from_numpy_quaternion([1, 2, 3, 4])

    def test_exact_and_max_denominator_are_mutually_exclusive(self):
        q = _FakeNumpyQuaternion(1, 2, 3, 4)
        with self.assertRaises(ValueError):
            Hy.from_numpy_quaternion(q, exact=True, max_denominator=10)

    def test_non_finite_coordinates_rejected(self):
        with self.assertRaises(ValueError):
            Hy.from_numpy_quaternion(_FakeNumpyQuaternion(1, float('nan'), 0, 0))
        with self.assertRaises(ValueError):
            Hy.from_quaternionic([1, 0, float('inf'), 0])

    # ---- from_quaternionic (duck-typed) ----------------------------- #
    def test_from_quaternionic_accepts_any_four_reals(self):
        expected = Hy(Hy(1, 2), Hy(3, 4))
        self.assertEqual(Hy.from_quaternionic([1, 2, 3, 4]), expected)
        self.assertEqual(Hy.from_quaternionic((1.0, 2.0, 3.0, 4.0)), expected)
        self.assertEqual(Hy.from_quaternionic([Fraction(1), 2, 3.0, 4]), expected)

    def test_from_quaternionic_unwraps_ndarray_attribute(self):
        class Wrapper:
            ndarray = [1, 2, 3, 4]
        self.assertEqual(Hy.from_quaternionic(Wrapper()), Hy(Hy(1, 2), Hy(3, 4)))

    def test_from_quaternionic_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            Hy.from_quaternionic([1, 2, 3])
        with self.assertRaises(ValueError):
            Hy.from_quaternionic([1, 2, 3, 4, 5])
        with self.assertRaises(TypeError):
            Hy.from_quaternionic("abcd")
        with self.assertRaises(TypeError):
            Hy.from_quaternionic(5)
        with self.assertRaises(TypeError):
            Hy.from_quaternionic([1, 2, 3, 'x'])

    # ---- from_sympy (duck-typed) ------------------------------------- #
    def test_from_sympy_uses_abcd_in_order(self):
        r = _FakeSympyRational
        h = Hy.from_sympy(_FakeSympyQuaternion(r(1), r(2), r(3), r(4)))
        self.assertEqual(h, Hy(Hy(1, 2), Hy(3, 4)))
        self.assertEqual(h.rank, 2)

    def test_from_sympy_rejects_non_quaternion(self):
        with self.assertRaises(TypeError):
            Hy.from_sympy(object())
        with self.assertRaises(TypeError):
            Hy.from_sympy("not a quaternion")

    # ---- lazy import / missing package ------------------------------ #
    def test_missing_packages_give_helpful_import_errors(self):
        q = Hy.from_array([0, 0, 0, 1])
        cases = (
            ("sympy", q.to_sympy, "sympy"),
            ("quaternion", q.to_numpy_quaternion, "numpy-quaternion"),
            ("quaternionic", q.to_quaternionic, "quaternionic"),
        )
        for module_name, method, pip_name in cases:
            with self.subTest(module=module_name):
                with mock.patch.dict(sys.modules, {module_name: None}):
                    with self.assertRaises(ImportError) as cm:
                        method()
                self.assertIn(f"pip install {pip_name}", str(cm.exception))


@unittest.skipUnless(HAVE_SYMPY, "sympy is not installed")
class TestSympyInterop(unittest.TestCase):

    def test_component_order_of_the_four_units(self):
        expected = {
            (1, 0, 0, 0): (1, 0, 0, 0),
            (0, 1, 0, 0): (0, 1, 0, 0),
            (0, 0, 1, 0): (0, 0, 1, 0),
            (0, 0, 0, 1): (0, 0, 0, 1),
        }
        for coeffs, wxyz in expected.items():
            with self.subTest(coeffs=coeffs):
                q = Hy.from_array(coeffs).to_sympy()
                self.assertIsInstance(q, sympy.Quaternion)
                self.assertEqual((q.a, q.b, q.c, q.d), wxyz)

    def test_conversion_is_exact_no_rounding_at_all(self):
        h = Hy.from_array(['1/3', '-2/7', '5/9', '1/10'])
        q = h.to_sympy()
        self.assertEqual((q.a, q.b, q.c, q.d),
                         (Fraction(1, 3), Fraction(-2, 7), Fraction(5, 9), Fraction(1, 10)))
        for c in (q.a, q.b, q.c, q.d):
            self.assertTrue(c.is_Rational)

    def test_complex_values_are_embedded_as_a_plus_bi(self):
        q = Hy('5/2', '-16/5').to_sympy()
        self.assertEqual((q.a, q.b, q.c, q.d),
                         (Fraction(5, 2), Fraction(-16, 5), 0, 0))

    def test_higher_ranks_rejected(self):
        for rank in (3, 4):
            with self.assertRaises(ValueError):
                Hy.random(rank, seed=rank).to_sympy()

    def test_basis_units_multiply_like_hy(self):
        i, j, k = (Hy.from_array(c).to_sympy()
                   for c in ([0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]))
        self.assertEqual(i * j, k)
        self.assertEqual(j * i, -k)
        self.assertEqual(i * i, sympy.Quaternion(-1, 0, 0, 0))

    def test_products_agree_with_hy_exactly(self):
        for seed in range(INTEROP_TRIALS):
            a, b = rand_quat(2 * seed), rand_quat(2 * seed + 1)
            theirs = a.to_sympy() * b.to_sympy()
            self.assertEqual((theirs.a, theirs.b, theirs.c, theirs.d),
                             tuple((a * b).to_array()))

    def test_from_sympy_known_values(self):
        h = Hy.from_sympy(sympy.Quaternion(1, sympy.Rational(2, 3), -4, 0))
        self.assertEqual(h, Hy.from_array([1, '2/3', -4, 0]))
        self.assertEqual(h.rank, 2)

    def test_from_sympy_rejects_non_quaternion(self):
        with self.assertRaises(TypeError):
            Hy.from_sympy([1, 2, 3, 4])

    def test_round_trip_of_exact_rational_components_has_no_rounding_at_all(self):
        # Every rational quaternion round-trips exactly, with no keywords needed,
        # because sympy stores its Rational/Integer components exactly too.
        for _ in range(INTEROP_TRIALS):
            h = rand_quat(random.randint(0, 10 ** 6))
            self.assertEqual(Hy.from_sympy(h.to_sympy()), h)

    def test_float_components_use_the_same_policy_as_from_numpy_quaternion(self):
        q = sympy.Quaternion(0.1, sympy.Float(1 / 3), 0.375, -2.0)
        self.assertEqual(Hy.from_sympy(q).to_array()[0], Fraction(1, 10))
        self.assertEqual(Hy.from_sympy(q, max_denominator=1000).to_array()[1], Fraction(1, 3))
        self.assertEqual(Hy.from_sympy(q, exact=True).to_array()[0], Fraction(0.1))

    def test_exact_and_max_denominator_are_mutually_exclusive(self):
        q = sympy.Quaternion(sympy.Float(0.1), 0, 0, 0)
        with self.assertRaises(ValueError):
            Hy.from_sympy(q, exact=True, max_denominator=10)

    def test_non_real_non_finite_and_symbolic_components_rejected(self):
        for bad in (sympy.nan, sympy.oo, sympy.I):
            with self.subTest(value=bad):
                with self.assertRaises(ValueError):
                    Hy.from_sympy(sympy.Quaternion(bad, 0, 0, 0))
        with self.assertRaises(TypeError):
            Hy.from_sympy(sympy.Quaternion(sympy.Symbol('x'), 0, 0, 0))


@unittest.skipUnless(HAVE_NPQUAT, "numpy-quaternion is not installed")
class TestNumpyQuaternionInterop(unittest.TestCase):

    def test_to_numpy_quaternion(self):
        q = Hy(Hy(1, 2), Hy(3, 4)).to_numpy_quaternion()
        self.assertIsInstance(q, np.quaternion)
        self.assertEqual((q.w, q.x, q.y, q.z), (1.0, 2.0, 3.0, 4.0))
        self.assertEqual(q, np.quaternion(1, 2, 3, 4))

    def test_to_numpy_quaternion_rational_coordinates_become_floats(self):
        q = Hy.from_array(['1/2', '-1/4', '1/3', '0']).to_numpy_quaternion()
        np.testing.assert_allclose(q.components, [0.5, -0.25, 1 / 3, 0.0], atol=0)

    def test_complex_values_are_embedded_as_a_plus_bi(self):
        q = Hy('5/2', '-16/5').to_numpy_quaternion()
        self.assertEqual(q, np.quaternion(2.5, -3.2, 0, 0))

    def test_higher_ranks_rejected(self):
        with self.assertRaises(ValueError):
            Hy.random(3, seed=1).to_numpy_quaternion()

    def test_from_numpy_quaternion(self):
        h = Hy.from_numpy_quaternion(np.quaternion(1.5, -0.25, 3, 0.1))
        self.assertEqual(h, Hy.from_array(['3/2', '-1/4', '3', '1/10']))
        self.assertEqual(h.rank, 2)

    def test_from_numpy_quaternion_rejects_arrays(self):
        arr = npquat.as_quat_array(np.arange(12.0).reshape(3, 4))
        with self.assertRaises(ValueError):
            Hy.from_numpy_quaternion(arr)
        with self.assertRaises(TypeError):
            Hy.from_numpy_quaternion([1, 2, 3, 4])

    def test_basis_units_multiply_like_hy(self):
        i, j, k = (Hy.from_array(c).to_numpy_quaternion()
                   for c in ([0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]))
        self.assertEqual(i * j, k)
        self.assertEqual(j * i, -k)
        self.assertEqual(i * i, np.quaternion(-1, 0, 0, 0))

    def test_products_agree_with_hy(self):
        for seed in range(INTEROP_TRIALS):
            a, b = rand_quat(2 * seed), rand_quat(2 * seed + 1)
            theirs = a.to_numpy_quaternion() * b.to_numpy_quaternion()
            np.testing.assert_allclose(theirs.components, as_floats(a * b), rtol=1e-12, atol=1e-12)

    def test_float_round_trip_is_exact_in_float_space(self):
        rng = np.random.default_rng(2)
        for exact in (False, True):
            for _ in range(INTEROP_TRIALS):
                q = np.quaternion(*rng.normal(size=4))
                back = Hy.from_numpy_quaternion(q, exact=exact).to_numpy_quaternion()
                self.assertEqual(back, q)

    def test_rational_round_trip_policies(self):
        h = Hy.from_array(['1/3', '-2/7', '5/9', '1/10'])
        back = Hy.from_numpy_quaternion(h.to_numpy_quaternion())
        self.assertNotEqual(back, h)                                    # default: decimal digits
        self.assertEqual(Hy.from_numpy_quaternion(h.to_numpy_quaternion(),
                                                  max_denominator=1000), h)
        # dyadic rationals survive every policy
        d = Hy.from_array(['1/2', '-3/8', '5', '7/16'])
        for kw in ({}, {'exact': True}, {'max_denominator': 100}):
            self.assertEqual(Hy.from_numpy_quaternion(d.to_numpy_quaternion(), **kw), d)

    def test_non_finite_rejected(self):
        with self.assertRaises(ValueError):
            Hy.from_numpy_quaternion(np.quaternion(1, np.nan, 0, 0))


@unittest.skipUnless(HAVE_QUATERNIONIC, "quaternionic is not installed")
class TestQuaternionicInterop(unittest.TestCase):

    def test_to_quaternionic(self):
        q = Hy(Hy(1, 2), Hy(3, 4)).to_quaternionic()
        self.assertIsInstance(q, type(quaternionic.array([0, 0, 0, 1])))
        self.assertEqual(q.shape, (4,))
        self.assertEqual(q.dtype, np.float64)
        np.testing.assert_array_equal(q.ndarray, [1.0, 2.0, 3.0, 4.0])

    def test_complex_values_are_embedded_as_a_plus_bi(self):
        q = Hy('5/2', '-16/5').to_quaternionic()
        np.testing.assert_array_equal(q.ndarray, [2.5, -3.2, 0.0, 0.0])

    def test_higher_ranks_rejected(self):
        with self.assertRaises(ValueError):
            Hy.random(3, seed=1).to_quaternionic()

    def test_from_quaternionic(self):
        h = Hy.from_quaternionic(quaternionic.array([1, 0.5, 0, -2]))
        self.assertEqual(h, Hy.from_array([1, '1/2', 0, -2]))
        self.assertEqual(h.rank, 2)

    def test_from_quaternionic_accepts_plain_arrays_and_sequences(self):
        expected = Hy.from_array([1, 2, 3, 4])
        self.assertEqual(Hy.from_quaternionic(np.array([1.0, 2.0, 3.0, 4.0])), expected)
        self.assertEqual(Hy.from_quaternionic([1, 2, 3, 4]), expected)

    def test_from_quaternionic_rejects_arrays_of_quaternions(self):
        for shape in ((3, 4), (4, 4), (2, 3, 4)):
            with self.subTest(shape=shape):
                with self.assertRaises(ValueError):
                    Hy.from_quaternionic(quaternionic.array(np.ones(shape)))
        with self.assertRaises(ValueError):
            Hy.from_quaternionic([1, 2, 3])
        # ...but each row of an array can be converted individually
        rows = quaternionic.array(np.arange(8.0).reshape(2, 4))
        self.assertEqual([Hy.from_quaternionic(r) for r in rows],
                         [Hy.from_array([0, 1, 2, 3]), Hy.from_array([4, 5, 6, 7])])

    def test_basis_units_multiply_like_hy(self):
        i, j, k = (Hy.from_array(c).to_quaternionic()
                   for c in ([0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]))
        np.testing.assert_array_equal((i * j).ndarray, k.ndarray)
        np.testing.assert_array_equal((j * i).ndarray, (-k).ndarray)
        np.testing.assert_array_equal((i * i).ndarray, [-1.0, 0.0, 0.0, 0.0])

    def test_products_agree_with_hy(self):
        for seed in range(INTEROP_TRIALS):
            a, b = rand_quat(2 * seed), rand_quat(2 * seed + 1)
            theirs = a.to_quaternionic() * b.to_quaternionic()
            np.testing.assert_allclose(theirs.ndarray, as_floats(a * b), rtol=1e-12, atol=1e-12)

    def test_float_round_trip_is_exact_in_float_space(self):
        rng = np.random.default_rng(3)
        for exact in (False, True):
            for _ in range(INTEROP_TRIALS):
                q = quaternionic.array(rng.normal(size=4))
                back = Hy.from_quaternionic(q, exact=exact).to_quaternionic()
                np.testing.assert_array_equal(back.ndarray, q.ndarray)

    def test_rational_round_trip_policies(self):
        h = Hy.from_array(['1/3', '-2/7', '5/9', '1/10'])
        self.assertNotEqual(Hy.from_quaternionic(h.to_quaternionic()), h)
        self.assertEqual(Hy.from_quaternionic(h.to_quaternionic(), max_denominator=1000), h)
        d = Hy.from_array(['1/2', '-3/8', '5', '7/16'])
        for kw in ({}, {'exact': True}, {'max_denominator': 100}):
            self.assertEqual(Hy.from_quaternionic(d.to_quaternionic(), **kw), d)

    def test_non_finite_rejected(self):
        with self.assertRaises(ValueError):
            Hy.from_quaternionic(quaternionic.array([1, np.inf, 0, 0]))


@unittest.skipUnless(HAVE_SYMPY and HAVE_NPQUAT and HAVE_QUATERNIONIC,
                     "needs sympy, numpy-quaternion and quaternionic")
class TestInteropAcrossPackages(unittest.TestCase):

    def test_all_three_packages_hold_the_same_coordinates(self):
        for seed in range(INTEROP_TRIALS):
            h = rand_quat(seed)
            wxyz = np.array(as_floats(h))
            npq = h.to_numpy_quaternion().components
            qnc = h.to_quaternionic().ndarray
            sy = h.to_sympy()
            np.testing.assert_array_equal(npq, wxyz)
            np.testing.assert_array_equal(qnc, wxyz)
            self.assertEqual((sy.a, sy.b, sy.c, sy.d), tuple(h.to_array()))   # exact

    def test_converting_between_the_packages_via_hy(self):
        q = quaternionic.array([0.5, -0.25, 1.5, 2.0])
        via_hy = Hy.from_quaternionic(q).to_numpy_quaternion()
        self.assertEqual(via_hy, np.quaternion(0.5, -0.25, 1.5, 2.0))
        back = Hy.from_numpy_quaternion(via_hy).to_quaternionic()
        np.testing.assert_array_equal(back.ndarray, q.ndarray)
        # sympy round-trips the same value exactly (it's already rational)
        exact = Hy.from_array(['1/2', '-1/4', '3/2', '2'])
        self.assertEqual(Hy.from_sympy(exact.to_sympy()), exact)


def main():
    unittest.main()


if __name__ == '__main__':
    main()
