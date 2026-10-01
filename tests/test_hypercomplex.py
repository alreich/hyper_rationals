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
    TestMatrixRepresentation
                            - .to_matrix()/Hy.from_matrix(): roundtrips,
                              the algebra-homomorphism property at ranks
                              0-2, its deliberate failure at rank >= 3,
                              the allow_nonassociative=True guard, and a
                              known rank-4 (sedenion) zero-divisor giving
                              a singular matrix    (needs numpy)
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
    TestSignatureBasics     - Hy(..., mu=...), .mu / .signs, validation,
                              mixed-rank embedding, copy semantics, the
                              named presets and ``signs=`` arguments
    TestSplitComplex        - split-complex numbers (mu = +1)
    TestSplitQuaternions    - split-quaternions (signs (-1, +1))
    TestSplitOctonions      - split-octonions (signs (-1, -1, +1))
    TestGeneralRationalMu   - arbitrary nonzero rational mu
    TestSignatureAlgebraFuzz
                            - algebraic laws across many signatures,
                              ranks 1-4 (unit squares, norm
                              multiplicativity, (non)associativity, ...)
    TestSignatureEqualityHashRepr
                            - equality/hash across signatures; repr only
                              shows a non-default signature
    TestSignatureCoercionAndMixing
                            - scalars, strings, complex literals and
                              mixed signatures/ranks in arithmetic
    TestSignatureParsing    - Hy.parse(..., signs=...)
    TestSignatureNormAndPredicates
                            - norm()/norm_squared(), abs(), is_null(),
                              pow with null elements
    TestSignatureMatrixRepresentation
                            - to_matrix()/from_matrix() with signs (needs numpy)
    TestSignatureInterop    - interop converters reject non-classical algebras
    TestSignatureHelpers    - module-level functions, immutability, latex, units

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
HAVE_NUMPY = np is not None

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
# Hy.to_matrix() / Hy.from_matrix()
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAVE_NUMPY, "numpy not installed")
class TestMatrixRepresentation(unittest.TestCase):

    def test_rank1_hand_checked_matrix(self):
        z = Hy("2", "3")
        M = z.to_matrix()
        self.assertEqual(M.shape, (2, 2))
        self.assertEqual(M.dtype, object)
        expected = [[Fraction(2), Fraction(-3)], [Fraction(3), Fraction(2)]]
        self.assertTrue(np.array_equal(M, np.array(expected, dtype=object)))

    def test_roundtrip_ranks_1_and_2(self):
        for rank in (1, 2):
            rng = random.Random(7000 + rank)
            for _ in range(FUZZ_TRIALS):
                x = rand_value(rank, rng)
                self.assertEqual(Hy.from_matrix(x.to_matrix()), x)
                self.assertEqual(
                    Hy.from_matrix(x.to_matrix(kind="right"), kind="right"), x
                )

    def test_homomorphism_holds_ranks_1_and_2(self):
        for rank in (1, 2):
            rng = random.Random(7100 + rank)
            for _ in range(FUZZ_TRIALS):
                x = rand_value(rank, rng)
                y = rand_value(rank, rng)
                lhs = np.dot(x.to_matrix(), y.to_matrix())
                rhs = (x * y).to_matrix()
                self.assertTrue(np.array_equal(lhs, rhs), (x, y))

    @unittest.skipUnless(HAVE_SYMPY, "sympy not installed")
    def test_determinant_equals_norm_power_ranks_1_and_2(self):
        for rank in (1, 2):
            rng = random.Random(7200 + rank)
            for _ in range(10):
                x = rand_value(rank, rng)
                if x.is_zero():
                    continue
                det = sympy.Matrix(x.to_matrix().tolist()).det()
                n = x.norm()
                expected = sympy.Rational(n.numerator, n.denominator) ** (2 ** (rank - 1))
                self.assertEqual(det, expected, x)

    def test_kind_left_and_right_agree_for_commutative_ranks(self):
        # rank 1 (complex) is commutative, so left/right regular
        # representations coincide.
        rng = random.Random(7300)
        for _ in range(10):
            x = rand_value(1, rng)
            self.assertTrue(np.array_equal(x.to_matrix(), x.to_matrix(kind="right")))

    def test_kind_left_and_right_can_differ_for_rank2(self):
        # Quaternions are noncommutative, so left/right reps generally differ.
        i = Hy(Hy(0, 1), Hy(0, 0))
        j = Hy(Hy(0, 0), Hy(1, 0))
        self.assertFalse(np.array_equal(i.to_matrix(), j.to_matrix(kind="right")))

    def test_rank3_to_matrix_raises_without_flag(self):
        o = Hy.random(3, seed=1)
        with self.assertRaises(ValueError):
            o.to_matrix()

    def test_rank3_to_matrix_allowed_with_flag(self):
        o = Hy.random(3, seed=1)
        M = o.to_matrix(allow_nonassociative=True)
        self.assertEqual(M.shape, (8, 8))

    def test_rank3_homomorphism_fails_in_general(self):
        # Non-associativity means M(x) @ M(y) != M(x*y) for at least some
        # octonion pairs -- confirm this is actually true rather than
        # just documented.
        mismatch_found = False
        for seed in range(20):
            a = Hy.random(3, seed=seed)
            b = Hy.random(3, seed=seed + 500)
            lhs = np.dot(
                a.to_matrix(allow_nonassociative=True),
                b.to_matrix(allow_nonassociative=True),
            )
            rhs = (a * b).to_matrix(allow_nonassociative=True)
            if not np.array_equal(lhs, rhs):
                mismatch_found = True
                break
        self.assertTrue(mismatch_found)

    def test_rank3_from_matrix_raises_without_flag(self):
        o = Hy.random(3, seed=2)
        M = o.to_matrix(allow_nonassociative=True)
        with self.assertRaises(ValueError):
            Hy.from_matrix(M)

    def test_rank3_from_matrix_allowed_with_flag_roundtrips(self):
        o = Hy.random(3, seed=2)
        M = o.to_matrix(allow_nonassociative=True)
        self.assertEqual(Hy.from_matrix(M, allow_nonassociative=True), o)

    def test_sedenion_zero_divisor_gives_singular_matrix(self):
        # A known rank-4 zero divisor under this library's basis ordering:
        # (e1 + e10) * (e4 - e15) == 0, even though neither factor is zero.
        units = Hy.units(4)
        positive = [u for name, u in units.items() if not name.startswith("-")]
        a = positive[1] + positive[10]
        b = positive[4] - positive[15]
        self.assertFalse(a.is_zero())
        self.assertFalse(b.is_zero())
        self.assertTrue((a * b).is_zero())
        M = a.to_matrix(allow_nonassociative=True, as_float=True)
        self.assertAlmostEqual(np.linalg.det(M), 0.0, places=6)

    def test_as_float_option(self):
        x = Hy("5/2", "-1/4")
        M = x.to_matrix(as_float=True)
        self.assertEqual(M.dtype, np.float64)
        self.assertAlmostEqual(M[0, 0], 2.5)

    def test_invalid_kind_raises(self):
        x = Hy(1, 2)
        with self.assertRaises(ValueError):
            x.to_matrix(kind="sideways")
        with self.assertRaises(ValueError):
            Hy.from_matrix(x.to_matrix(), kind="sideways")

    def test_from_matrix_infers_rank_from_shape(self):
        x = Hy(Hy(1, 2), Hy(3, 4))
        result = Hy.from_matrix(x.to_matrix())
        self.assertEqual(result.rank, 2)

    def test_from_matrix_rank_mismatch_raises(self):
        x = Hy(Hy(1, 2), Hy(3, 4))  # rank 2 -> 4x4 matrix
        with self.assertRaises(ValueError):
            Hy.from_matrix(x.to_matrix(), rank=1)

    def test_from_matrix_rejects_bad_shapes(self):
        with self.assertRaises(ValueError):
            Hy.from_matrix(np.zeros((2, 3)))          # not square
        with self.assertRaises(ValueError):
            Hy.from_matrix(np.zeros((3, 3)))           # not a power of 2
        with self.assertRaises(ValueError):
            Hy.from_matrix(np.zeros((1, 1)))           # too small (rank 0)

    def test_from_matrix_validate_true_detects_inconsistent_matrix(self):
        x = Hy(Hy(1, 2), Hy(3, 4))
        M = x.to_matrix()
        bad = np.array(M, dtype=object, copy=True)
        bad[2, 3] = bad[2, 3] + 1  # corrupt one entry away from column 0
        with self.assertRaises(ValueError):
            Hy.from_matrix(bad, validate=True)

    def test_from_matrix_validate_true_passes_for_genuine_matrix(self):
        rng = random.Random(7400)
        for _ in range(10):
            x = rand_value(2, rng)
            M = x.to_matrix()
            self.assertEqual(Hy.from_matrix(M, validate=True), x)


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


# ---------------------------------------------------------------------------
# Signatures: split-hypercomplex numbers and general rational mu
# ---------------------------------------------------------------------------

# A battery of signatures (lowest doubling level first) used by the fuzz
# tests: classical, each split pattern through rank 3, and general rational
# mu, plus a few rank-4 ones.
SIGNATURES_UP_TO_RANK_3 = (
    (-1,), (1,), (2,), (Fraction(-3, 2),),
    (-1, -1), (-1, 1), (1, -1), (1, 1), (-1, 2), (3, -5),
    (-1, -1, -1), (-1, -1, 1), (1, 1, 1), (-1, -1, 2), (Fraction(1, 2), -1, 3),
)
SIGNATURES_RANK_4 = ((-1, -1, -1, -1), (-1, -1, 1, -1), (1, -1, 1, 1))


def rand_signed(signs, rng):
    """A random Hy living in the algebra with the given signature."""
    return Hy.random(signs=signs, rng=rng)


class TestSignatureBasics(unittest.TestCase):
    """The mu / signs properties, construction, validation, presets."""

    def test_default_signature_is_all_minus_one(self):
        for rank in (1, 2, 3, 4):
            h = Hy.random(rank, seed=rank)
            self.assertEqual(h.signs, (-1,) * rank)
            self.assertEqual(h.mu, -1)

    def test_mu_is_stored_as_a_fraction(self):
        for given in (1, Fraction(1), '1', 1.0):
            h = Hy(1, 2, mu=given)
            self.assertIsInstance(h.mu, Fraction)
            self.assertEqual(h.mu, 1)
        self.assertEqual(Hy(1, 2, mu='3/2').mu, Fraction(3, 2))
        self.assertEqual(Hy(1, 2, mu=0.25).mu, Fraction(1, 4))

    def test_mu_must_be_nonzero_rational(self):
        for bad in (0, 0.0, '0', Fraction(0)):
            with self.assertRaises(ValueError):
                Hy(1, 2, mu=bad)
        with self.assertRaises(ValueError):
            Hy(1, 2, mu='not a number')
        for bad in (True, False, [1], None.__class__):
            with self.assertRaises(TypeError):
                Hy(1, 2, mu=bad)

    def test_top_mu_is_independent_of_component_mus(self):
        inner = Hy(1, 2, mu=1), Hy(3, 4, mu=1)
        h = Hy(*inner)                       # top level defaults to -1
        self.assertEqual(h.signs, (1, -1))
        h2 = Hy(*inner, mu=3)
        self.assertEqual(h2.signs, (1, 3))
        self.assertEqual(h2.mu, 3)
        self.assertEqual(h2.real.mu, 1)

    def test_components_must_share_a_signature(self):
        with self.assertRaises(ValueError):
            Hy(Hy(1, 2, mu=1), Hy(3, 4))      # (1,) vs (-1,)

    def test_mixed_rank_components_embed_when_signature_is_a_prefix(self):
        # a bare rational pairs with a split-complex value: the rational
        # is embedded using the *other* component's signature
        h = Hy(3, Hy(1, 2, mu=1))
        self.assertEqual(h.signs, (1, -1))
        self.assertEqual(h.real, Hy(3, 0, mu=1))
        # rank-1 (-1) value pairs with a rank-2 (-1, 1) value
        q = Hy(Hy(1, 2), Hy(3, 4), mu=1)
        h3 = Hy(Hy(5, 6), Hy(q, 0))              # Hy(q, 0) is rank 3
        self.assertEqual(h3.signs, (-1, 1, -1, -1))
        self.assertEqual(h3.real.signs, (-1, 1, -1))
        with self.assertRaises(ValueError):
            Hy(Hy(5, 6, mu=1), q)            # (1,) is not a prefix of (-1, 1)

    def test_copy_construction_keeps_signature(self):
        q = Hy(Hy(1, 2), Hy(3, 4), mu=1)
        c = Hy(q)
        self.assertEqual(c, q)
        self.assertEqual(c.signs, q.signs)
        self.assertEqual(Hy(q, mu=1).signs, q.signs)   # same mu: fine

    def test_copy_construction_cannot_change_mu(self):
        q = Hy(Hy(1, 2), Hy(3, 4), mu=1)
        with self.assertRaises(ValueError):
            Hy(q, mu=-1)

    def test_plain_scalar_construction_with_mu(self):
        h = Hy(3, mu=1)
        self.assertEqual(h.rank, 1)
        self.assertEqual(h.signs, (1,))
        self.assertEqual(h, 3)

    def test_preset_constants(self):
        self.assertEqual(Hy.COMPLEX, (-1,))
        self.assertEqual(Hy.QUATERNION, (-1, -1))
        self.assertEqual(Hy.OCTONION, (-1, -1, -1))
        self.assertEqual(Hy.SEDENION, (-1, -1, -1, -1))
        self.assertEqual(Hy.SPLIT_COMPLEX, (1,))
        self.assertEqual(Hy.SPLIT_QUATERNION, (-1, 1))
        self.assertEqual(Hy.SPLIT_OCTONION, (-1, -1, 1))

    def test_presets_accepted_as_tuples_or_names(self):
        a = Hy.random(signs=Hy.SPLIT_QUATERNION, seed=1)
        b = Hy.random(signs="split-quaternion", seed=1)
        c = Hy.random(signs="Split_Quaternion", seed=1)
        d = Hy.random(signs="split quaternion", seed=1)
        self.assertEqual(a, b)
        self.assertEqual(a, c)
        self.assertEqual(a, d)
        self.assertEqual(a.signs, (-1, 1))

    def test_every_named_preset_matches_its_constant(self):
        pairs = {
            "complex": Hy.COMPLEX, "quaternion": Hy.QUATERNION,
            "octonion": Hy.OCTONION, "sedenion": Hy.SEDENION,
            "split-complex": Hy.SPLIT_COMPLEX,
            "split-quaternion": Hy.SPLIT_QUATERNION,
            "split-octonion": Hy.SPLIT_OCTONION,
        }
        for name, tup in pairs.items():
            self.assertEqual(Hy.random(signs=name, seed=0).signs, tup)

    def test_unknown_preset_raises(self):
        with self.assertRaises(ValueError):
            Hy.random(signs="split-sedenion-ish")
        with self.assertRaises(ValueError):
            Hy.units(signs="nope")

    def test_signs_must_be_a_name_or_sequence(self):
        with self.assertRaises(TypeError):
            Hy.random(2, signs=5)

    def test_rank_is_inferred_from_signs(self):
        self.assertEqual(Hy.random(signs="split-octonion", seed=1).rank, 3)
        self.assertEqual(len(Hy.units(signs=(1, -1))), 8)

    def test_rank_and_signs_must_agree(self):
        with self.assertRaises(ValueError):
            Hy.random(3, signs="split-quaternion")
        with self.assertRaises(ValueError):
            Hy.units(2, signs=(1,))
        with self.assertRaises(ValueError):
            Hy.from_array([1, 2, 3, 4], signs=(1,))

    def test_rank_or_signs_required(self):
        with self.assertRaises(ValueError):
            Hy.random()
        with self.assertRaises(ValueError):
            Hy.units()

    def test_from_array_and_to_array_with_signs(self):
        x = Hy.from_array([1, 2, 3, 4], signs="split-quaternion")
        self.assertEqual(x.signs, (-1, 1))
        self.assertEqual(x.to_array(as_str=True), ['1', '2', '3', '4'])
        self.assertEqual(Hy.from_array(x.to_array(), signs=x.signs), x)

    def test_random_with_signs_is_reproducible(self):
        a = Hy.random(signs="split-octonion", seed=42)
        b = Hy.random(signs="split-octonion", seed=42)
        self.assertEqual(a, b)
        self.assertEqual(a.signs, (-1, -1, 1))


class TestSplitComplex(unittest.TestCase):
    """The split-complex numbers: Cayley-Dickson with mu = +1."""

    def setUp(self):
        self.j = Hy(0, 1, mu=1)

    def test_j_squares_to_plus_one(self):
        self.assertEqual(self.j * self.j, 1)
        self.assertEqual(str(self.j * self.j), '(1)')

    def test_fixed_product(self):
        # (a + bj)(c + dj) = (ac + bd) + (ad + bc) j
        x, y = Hy(2, 3, mu=1), Hy(5, 7, mu=1)
        self.assertEqual(x * y, Hy(2 * 5 + 3 * 7, 2 * 7 + 3 * 5, mu=1))

    def test_commutative_and_associative(self):
        rng = random.Random(1)
        for _ in range(FUZZ_TRIALS):
            x, y, z = (rand_signed((1,), rng) for _ in range(3))
            self.assertEqual(x * y, y * x)
            self.assertEqual((x * y) * z, x * (y * z))

    def test_null_elements_and_zero_divisors(self):
        p, m = Hy(1, 1, mu=1), Hy(1, -1, mu=1)
        self.assertEqual(p * m, 0)
        self.assertTrue(p.is_null())
        self.assertTrue(m.is_null())
        self.assertFalse(bool(p * m))
        self.assertTrue(bool(p))

    def test_null_elements_are_not_invertible(self):
        for null in (Hy(1, 1, mu=1), Hy(-3, 3, mu=1)):
            with self.assertRaises(ZeroDivisionError) as ctx:
                null.inverse()
            self.assertIn('null', str(ctx.exception))
            with self.assertRaises(ZeroDivisionError):
                Hy(1, 0, mu=1) / null
            with self.assertRaises(ZeroDivisionError):
                null ** -1

    def test_zero_is_not_invertible_either(self):
        with self.assertRaises(ZeroDivisionError):
            Hy(0, 0, mu=1).inverse()

    def test_inverse_of_non_null(self):
        x = Hy(2, 1, mu=1)                       # N = 4 - 1 = 3
        self.assertEqual(x.norm_squared(), 3)
        self.assertEqual(x.inverse(), Hy(Fraction(2, 3), Fraction(-1, 3), mu=1))
        self.assertEqual(x * x.inverse(), 1)
        self.assertEqual(x.inverse() * x, 1)

    def test_norm_squared_is_the_hyperbolic_quadratic_form(self):
        self.assertEqual(Hy(5, 3, mu=1).norm_squared(), 16)      # a^2 - b^2
        self.assertEqual(Hy(3, 5, mu=1).norm_squared(), -16)
        self.assertEqual(Hy(4, 4, mu=1).norm_squared(), 0)

    def test_norm_is_multiplicative(self):
        rng = random.Random(2)
        for _ in range(FUZZ_TRIALS):
            x, y = rand_signed((1,), rng), rand_signed((1,), rng)
            self.assertEqual((x * y).norm_squared(),
                             x.norm_squared() * y.norm_squared())

    def test_conjugate_flips_the_sign_of_j_only(self):
        x = Hy(2, 3, mu=1)
        self.assertEqual(x.conjugate(), Hy(2, -3, mu=1))
        self.assertEqual(x * x.conjugate(), x.norm_squared())

    def test_idempotents(self):
        # e = (1 + j)/2 and 1 - e are orthogonal idempotents
        e = Hy(Fraction(1, 2), Fraction(1, 2), mu=1)
        f = 1 - e
        self.assertEqual(e * e, e)
        self.assertEqual(f * f, f)
        self.assertEqual(e * f, 0)

    def test_units_and_is_unit(self):
        u = Hy.units(signs=Hy.SPLIT_COMPLEX)
        self.assertEqual(list(u), ['1', '-1', 'j', '-j'])
        self.assertEqual(u['j'], self.j)
        self.assertEqual(u['j'].signs, (1,))
        self.assertTrue(all(v.is_unit() for v in u.values()))
        # an invertible value that is not a basis unit
        self.assertFalse(Hy(2, 1, mu=1).is_unit())

    def test_pow(self):
        self.assertEqual(self.j ** 2, 1)
        self.assertEqual(self.j ** 3, self.j)
        self.assertEqual(self.j ** 0, Hy(1, 0, mu=1))
        x = Hy(2, 1, mu=1)
        self.assertEqual(x ** -2, (x * x).inverse())

    def test_hyperbolic_rotation_preserves_the_form(self):
        # (cosh, sinh) with rational stand-ins: (5/4, 3/4) has N = 1
        b = Hy(Fraction(5, 4), Fraction(3, 4), mu=1)
        self.assertEqual(b.norm_squared(), 1)
        v = Hy(7, 2, mu=1)
        self.assertEqual((b * v).norm_squared(), v.norm_squared())

    def test_abs_raises_when_the_form_is_negative(self):
        with self.assertRaises(ValueError):
            abs(Hy(3, 5, mu=1))
        self.assertAlmostEqual(abs(Hy(5, 3, mu=1)), 4.0)
        self.assertEqual(abs(Hy(4, 4, mu=1)), 0.0)

    def test_not_an_ordinary_complex_number(self):
        with self.assertRaises(ValueError):
            complex(Hy(1, 2, mu=1))
        self.assertEqual(complex(Hy(1, 2)), 1 + 2j)


class TestSplitQuaternions(unittest.TestCase):
    """Split-quaternions (coquaternions): signs (-1, +1)."""

    def setUp(self):
        u = Hy.units(signs="split-quaternion")
        self.one, self.i, self.j, self.k = (u[n] for n in ('1', 'i', 'j', 'k'))

    def test_squares_match_the_standard_convention(self):
        # i^2 = -1, j^2 = k^2 = +1 (Cockle)
        self.assertEqual(self.i * self.i, -1)
        self.assertEqual(self.j * self.j, 1)
        self.assertEqual(self.k * self.k, 1)

    def test_ij_equals_k_and_anticommutes(self):
        self.assertEqual(self.i * self.j, self.k)
        self.assertEqual(self.j * self.i, -self.k)
        self.assertEqual(self.j * self.k, -self.i)
        self.assertEqual(self.k * self.j, self.i)
        self.assertEqual(self.k * self.i, self.j)
        self.assertEqual(self.i * self.k, -self.j)

    def test_ijk_product(self):
        self.assertEqual(self.i * self.j * self.k, 1)

    def test_str_labels_are_positional(self):
        x = Hy.from_array([1, 2, 3, 4], signs="split-quaternion")
        self.assertEqual(str(x), '(1+2i+3j+4k)')

    def test_associative_fuzz(self):
        rng = random.Random(3)
        for _ in range(FUZZ_TRIALS):
            x, y, z = (rand_signed((-1, 1), rng) for _ in range(3))
            self.assertEqual((x * y) * z, x * (y * z))

    def test_not_commutative(self):
        self.assertNotEqual(self.i * self.j, self.j * self.i)

    def test_norm_squared_is_indefinite_with_signature_two_two(self):
        x = Hy.from_array([1, 2, 3, 4], signs="split-quaternion")
        # N = a^2 + b^2 - c^2 - d^2
        self.assertEqual(x.norm_squared(), 1 + 4 - 9 - 16)
        self.assertEqual(self.i.norm_squared(), 1)
        self.assertEqual(self.j.norm_squared(), -1)
        self.assertEqual(self.k.norm_squared(), -1)

    def test_norm_is_multiplicative(self):
        rng = random.Random(4)
        for _ in range(FUZZ_TRIALS):
            x, y = (rand_signed((-1, 1), rng) for _ in range(2))
            self.assertEqual((x * y).norm_squared(),
                             x.norm_squared() * y.norm_squared())

    def test_zero_divisors(self):
        a, b = self.one + self.j, self.one - self.j
        self.assertEqual(a * b, 0)
        self.assertTrue(a.is_null())
        with self.assertRaises(ZeroDivisionError):
            a.inverse()

    def test_nonzero_null_elements_are_exactly_the_zero_divisors(self):
        # at rank 2 the form's null cone is exactly the zero-divisor set
        rng = random.Random(5)
        found_null = 0
        for _ in range(400):
            x = Hy.random(signs="split-quaternion", rng=rng, lo=-2, hi=2, dmax=1)
            y = Hy.random(signs="split-quaternion", rng=rng, lo=-2, hi=2, dmax=1)
            if x and y and not (x * y):
                self.assertTrue(x.is_null() and y.is_null())
            if x.is_null():
                found_null += 1
        self.assertGreater(found_null, 0)

    def test_inverse_both_sides(self):
        rng = random.Random(6)
        done = 0
        while done < FUZZ_TRIALS:
            x = rand_signed((-1, 1), rng)
            if x.is_zero() or x.is_null():
                continue
            done += 1
            self.assertEqual(x * x.inverse(), 1)
            self.assertEqual(x.inverse() * x, 1)

    def test_alternative_split_forms_all_give_split_quaternions(self):
        # (-1,+1), (+1,-1) and (+1,+1) all produce algebras with
        # indefinite norm, zero divisors, and the same dimension
        for sg in ((-1, 1), (1, -1), (1, 1)):
            a = Hy.units(signs=sg)
            lone = a['j'] * a['j']
            self.assertIn(lone.norm_squared(), (1,))
            x = (a['1'] + a['i']) * (a['1'] - a['i'])
            self.assertEqual(x.rank, 2)

    def test_isomorphic_to_2x2_matrices_homomorphism(self):
        if not HAVE_NUMPY:
            self.skipTest("numpy not installed")
        rng = random.Random(7)
        for _ in range(FUZZ_TRIALS):
            x, y = (rand_signed((-1, 1), rng) for _ in range(2))
            self.assertTrue(np.array_equal(
                np.dot(x.to_matrix(), y.to_matrix()), (x * y).to_matrix()))


class TestSplitOctonions(unittest.TestCase):
    """Split-octonions (Zorn): signs (-1, -1, +1)."""

    def setUp(self):
        self.u = Hy.units(signs="split-octonion")
        self.pos = {n: v for n, v in self.u.items() if not n.startswith('-')}

    def test_unit_squares(self):
        # i, j, k square to -1; L, iL, jL, kL to +1
        for name in ('i', 'j', 'k'):
            self.assertEqual(self.pos[name] * self.pos[name], -1)
        for name in ('L', 'iL', 'jL', 'kL'):
            self.assertEqual(self.pos[name] * self.pos[name], 1)

    def test_quaternion_subalgebra_is_the_usual_one(self):
        i, j, k = (self.pos[n] for n in 'ijk')
        self.assertEqual(i * j, k)
        self.assertEqual(j * k, i)
        self.assertEqual(k * i, j)

    def test_norm_has_four_plus_four_signature(self):
        norms = [v.norm_squared() for v in self.pos.values()]
        self.assertEqual(sorted(norms), [-1] * 4 + [1] * 4)

    def test_norm_is_multiplicative(self):
        rng = random.Random(8)
        for _ in range(FUZZ_TRIALS):
            x, y = (rand_signed((-1, -1, 1), rng) for _ in range(2))
            self.assertEqual((x * y).norm_squared(),
                             x.norm_squared() * y.norm_squared())

    def test_alternative_but_not_associative(self):
        rng = random.Random(9)
        nonassoc = 0
        for _ in range(FUZZ_TRIALS):
            x, y, z = (rand_signed((-1, -1, 1), rng) for _ in range(3))
            self.assertEqual((x * x) * y, x * (x * y))        # left alternative
            self.assertEqual((y * x) * x, y * (x * x))        # right alternative
            if (x * y) * z != x * (y * z):
                nonassoc += 1
        self.assertGreater(nonassoc, 0)

    def test_has_zero_divisors(self):
        L = self.pos['L']
        a, b = 1 + L, 1 - L
        self.assertEqual(a * b, 0)
        self.assertTrue(a.is_null())
        with self.assertRaises(ZeroDivisionError):
            a.inverse()

    def test_inverse_of_non_null(self):
        rng = random.Random(10)
        done = 0
        while done < FUZZ_TRIALS:
            x = rand_signed((-1, -1, 1), rng)
            if x.is_zero() or x.is_null():
                continue
            done += 1
            self.assertEqual(x * x.inverse(), 1)
            self.assertEqual(x.inverse() * x, 1)

    def test_str_uses_the_octonion_labels(self):
        x = Hy.from_array([1, 2, 3, 4, 5, 6, 7, 8], signs="split-octonion")
        self.assertEqual(str(x), '(1+2i+3j+4k+5L+6iL+7jL+8kL)')


class TestGeneralRationalMu(unittest.TestCase):
    """mu may be any nonzero rational: quaternion algebras (a, b) etc."""

    def test_j_squares_to_mu(self):
        for mu in (2, Fraction(-3, 2), 5, Fraction(1, 4)):
            j = Hy(0, 1, mu=mu)
            self.assertEqual(j * j, mu)

    def test_quaternion_algebra_with_mu_two(self):
        # (-1, 2): i^2 = -1, j^2 = 2, k^2 = -(i^2 j^2) = 2 ... via the formula
        u = Hy.units(signs=(-1, 2))
        i, j, k = u['i'], u['j'], u['k']
        self.assertEqual(i * i, -1)
        self.assertEqual(j * j, 2)
        self.assertEqual(i * j, k)
        self.assertEqual(k * k, Hy.unit_square(3, (-1, 2)))
        self.assertEqual(k * k, 2)

    def test_definite_forms_are_division_algebras_up_to_rank_two(self):
        # mu = -2: N(a, b) = a^2 + 2 b^2 > 0 for nonzero values
        rng = random.Random(11)
        for _ in range(FUZZ_TRIALS):
            x = rand_signed((-2,), rng)
            self.assertGreaterEqual(x.norm_squared(), 0)
            if x:
                self.assertGreater(x.norm_squared(), 0)
                self.assertEqual(x * x.inverse(), 1)

    def test_abs_for_a_positive_definite_general_mu(self):
        self.assertAlmostEqual(abs(Hy(1, 1, mu=-2)), math.sqrt(3))

    def test_non_square_mu_can_still_make_a_field(self):
        # mu = 2: N = a^2 - 2 b^2, which never vanishes for nonzero
        # rationals (sqrt 2 is irrational), so every nonzero value is
        # invertible even though the form is indefinite over the reals
        rng = random.Random(12)
        for _ in range(FUZZ_TRIALS):
            x = rand_signed((2,), rng)
            if x:
                self.assertNotEqual(x.norm_squared(), 0)
                self.assertFalse(x.is_null())
                self.assertEqual(x * x.inverse(), 1)

    def test_square_mu_is_split(self):
        # mu = 4 = 2^2: (2 + j)(2 - j) = 4 - 4 = 0, a null element
        x, y = Hy(2, 1, mu=4), Hy(2, -1, mu=4)
        self.assertEqual(x * y, 0)
        self.assertTrue(x.is_null())

    def test_repr_with_rational_mu(self):
        h = Hy(1, 2, mu='3/2')
        self.assertEqual(repr(h), "Hy('1', '2', mu='3/2')")
        self.assertEqual(eval(repr(h)), h)


class TestSignatureAlgebraFuzz(unittest.TestCase):
    """Algebraic laws across many signatures, ranks 1-4."""

    def test_unit_square_formula_matches_actual_squares(self):
        for sg in SIGNATURES_UP_TO_RANK_3 + SIGNATURES_RANK_4:
            pos = [v for n, v in Hy.units(signs=sg).items()
                   if not n.startswith('-')]
            for idx, u in enumerate(pos):
                self.assertEqual(u * u, Hy.unit_square(idx, sg),
                                 msg=f"signs={sg}, index={idx}")

    def test_unit_square_validation(self):
        with self.assertRaises(ValueError):
            Hy.unit_square(4, (-1, -1))
        with self.assertRaises(ValueError):
            Hy.unit_square(-1, (-1, -1))
        with self.assertRaises(ValueError):
            Hy.unit_square(True, (-1, -1))
        self.assertEqual(Hy.unit_square(0, "split-quaternion"), 1)

    def test_classical_signature_makes_every_imaginary_unit_square_to_minus_one(self):
        for rank in (1, 2, 3, 4):
            sg = (-1,) * rank
            for idx in range(1, 2 ** rank):
                self.assertEqual(Hy.unit_square(idx, sg), -1)

    def test_x_times_conjugate_is_norm_squared(self):
        for n, sg in enumerate(SIGNATURES_UP_TO_RANK_3 + SIGNATURES_RANK_4):
            rng = random.Random(1000 + n)
            for _ in range(FUZZ_TRIALS):
                x = rand_signed(sg, rng)
                self.assertEqual(x * x.conjugate(), x.norm_squared())
                self.assertEqual(x.conjugate() * x, x.norm_squared())

    def test_conjugate_of_product_reverses_order(self):
        for n, sg in enumerate(SIGNATURES_UP_TO_RANK_3 + SIGNATURES_RANK_4):
            rng = random.Random(1100 + n)
            for _ in range(FUZZ_TRIALS):
                x, y = rand_signed(sg, rng), rand_signed(sg, rng)
                self.assertEqual((x * y).conjugate(),
                                 y.conjugate() * x.conjugate())

    def test_distributivity(self):
        for n, sg in enumerate(SIGNATURES_UP_TO_RANK_3 + SIGNATURES_RANK_4):
            rng = random.Random(1200 + n)
            for _ in range(FUZZ_TRIALS):
                a, b, c = (rand_signed(sg, rng) for _ in range(3))
                self.assertEqual(a * (b + c), a * b + a * c)
                self.assertEqual((a + b) * c, a * c + b * c)

    def test_norm_is_multiplicative_up_to_rank_three(self):
        for n, sg in enumerate(SIGNATURES_UP_TO_RANK_3):
            rng = random.Random(1300 + n)
            for _ in range(FUZZ_TRIALS):
                x, y = rand_signed(sg, rng), rand_signed(sg, rng)
                self.assertEqual((x * y).norm_squared(),
                                 x.norm_squared() * y.norm_squared(),
                                 msg=f"signs={sg}")

    def test_associative_up_to_rank_two(self):
        for n, sg in enumerate(s for s in SIGNATURES_UP_TO_RANK_3 if len(s) <= 2):
            rng = random.Random(1400 + n)
            for _ in range(FUZZ_TRIALS):
                x, y, z = (rand_signed(sg, rng) for _ in range(3))
                self.assertEqual((x * y) * z, x * (y * z))

    def test_alternative_at_rank_three(self):
        for n, sg in enumerate(s for s in SIGNATURES_UP_TO_RANK_3 if len(s) == 3):
            rng = random.Random(1500 + n)
            for _ in range(FUZZ_TRIALS):
                x, y = rand_signed(sg, rng), rand_signed(sg, rng)
                self.assertEqual((x * x) * y, x * (x * y))
                self.assertEqual((y * x) * x, y * (x * x))

    def test_not_associative_at_rank_three(self):
        for n, sg in enumerate(s for s in SIGNATURES_UP_TO_RANK_3 if len(s) == 3):
            rng = random.Random(1600 + n)
            found = False
            for _ in range(60):
                x, y, z = (rand_signed(sg, rng) for _ in range(3))
                if (x * y) * z != x * (y * z):
                    found = True
                    break
            self.assertTrue(found, msg=f"signs={sg}")

    def test_inverse_both_sides_for_invertible_elements(self):
        for n, sg in enumerate(SIGNATURES_UP_TO_RANK_3 + SIGNATURES_RANK_4):
            rng = random.Random(1700 + n)
            done = 0
            while done < FUZZ_TRIALS:
                x = rand_signed(sg, rng)
                if x.norm_squared() == 0:
                    continue
                done += 1
                one = _embed_signs_for_test(1, sg)
                self.assertEqual(x * x.inverse(), one)
                self.assertEqual(x.inverse() * x, one)

    def test_division_round_trip_up_to_rank_two(self):
        # (x / y) * y == x needs associativity, so rank <= 2 only
        for n, sg in enumerate(s for s in SIGNATURES_UP_TO_RANK_3 if len(s) <= 2):
            rng = random.Random(1800 + n)
            done = 0
            while done < FUZZ_TRIALS:
                x, y = rand_signed(sg, rng), rand_signed(sg, rng)
                if y.norm_squared() == 0:
                    continue
                done += 1
                self.assertEqual((x / y) * y, x)

    def test_all_default_signature_agrees_with_the_classical_module_functions(self):
        # for signs all -1, the signature-aware product equals the
        # classical (ac - conj(d) b, da + b conj(c)) formula
        rng = random.Random(1900)
        for rank in (1, 2, 3, 4):
            for _ in range(FUZZ_TRIALS):
                x = Hy.random(rank, rng=rng)
                y = Hy.random(rank, rng=rng)
                a, b, c, d = x.real, x.imag, y.real, y.imag
                expected = Hy(sub(mul(a, c), mul(conj(d), b)),
                              add(mul(d, a), mul(b, conj(c))))
                self.assertEqual(x * y, expected)


def _embed_signs_for_test(value, signs):
    """The scalar `value` as an element of the algebra with `signs`."""
    return Hy.from_array([value] + [0] * (2 ** len(signs) - 1), signs=signs)


class TestSignatureEqualityHashRepr(unittest.TestCase):

    def test_different_signatures_are_unequal(self):
        self.assertNotEqual(Hy(1, 2, mu=1), Hy(1, 2))
        self.assertFalse(Hy(1, 2, mu=1) == Hy(1, 2))
        self.assertNotEqual(Hy(1, 2, mu=2), Hy(1, 2, mu=3))

    def test_same_signature_same_value_is_equal_and_hashes_equal(self):
        a, b = Hy(1, 2, mu=1), Hy(1, 2, mu=1)
        self.assertEqual(a, b)
        self.assertEqual(hash(a), hash(b))

    def test_unequal_signatures_usually_hash_differently(self):
        self.assertNotEqual(hash(Hy(1, 2, mu=1)), hash(Hy(1, 2)))

    def test_real_values_are_equal_across_signatures(self):
        # a real number lies in every algebra of the tower
        self.assertEqual(Hy(2, 0, mu=1), Hy(2, 0))
        self.assertEqual(hash(Hy(2, 0, mu=1)), hash(Hy(2, 0)))
        self.assertEqual(Hy(2, 0, mu=1), 2)
        self.assertEqual(Hy(0, 0, mu=1), 0)

    def test_cross_rank_equality_within_one_signature(self):
        sc = Hy(1, 2, mu=1)
        embedded = Hy(sc, 0)                       # same value, rank 2
        self.assertEqual(embedded.rank, 2)
        self.assertEqual(embedded, sc)
        self.assertEqual(hash(embedded), hash(sc))

    def test_value_equality_requires_matching_levels_that_are_used(self):
        # the unused top level's mu does not matter ...
        a = Hy(Hy(1, 2), 0, mu=1)
        b = Hy(Hy(1, 2), 0, mu=5)
        self.assertEqual(a, b)
        self.assertEqual(hash(a), hash(b))
        # ... but an inner level's does
        self.assertNotEqual(Hy(Hy(1, 2, mu=1), 0), Hy(Hy(1, 2), 0))

    def test_usable_in_sets_and_dicts(self):
        s = {Hy(1, 2, mu=1), Hy(1, 2), Hy(1, 2, mu=1)}
        self.assertEqual(len(s), 2)

    def test_str_is_unchanged_by_signature(self):
        self.assertEqual(str(Hy(1, 2, mu=1)), str(Hy(1, 2)))
        self.assertEqual(str(Hy(1, 2, mu=1)), '(1+2j)')

    def test_repr_omits_default_signature(self):
        self.assertEqual(repr(Hy(1, 2)), "Hy('1', '2')")
        self.assertEqual(repr(Hy(Hy(1, 2), Hy(3, 4))),
                         "Hy(Hy('1', '2'), Hy('3', '4'))")
        self.assertNotIn('mu', repr(Hy.random(3, seed=1)))

    def test_repr_shows_non_default_signature(self):
        self.assertEqual(repr(Hy(1, 2, mu=1)), "Hy('1', '2', mu=1)")
        self.assertEqual(repr(Hy(Hy(1, 2), Hy(3, 4), mu=1)),
                         "Hy(Hy('1', '2'), Hy('3', '4'), mu=1)")
        self.assertEqual(repr(Hy(Hy(1, 2, mu=1), Hy(3, 4, mu=1))),
                         "Hy(Hy('1', '2', mu=1), Hy('3', '4', mu=1))")

    def test_repr_round_trips_for_many_signatures(self):
        for n, sg in enumerate(SIGNATURES_UP_TO_RANK_3 + SIGNATURES_RANK_4):
            rng = random.Random(2000 + n)
            for _ in range(5):
                x = rand_signed(sg, rng)
                y = eval(repr(x))
                self.assertEqual(y, x)
                self.assertEqual(y.signs, x.signs)


class TestSignatureCoercionAndMixing(unittest.TestCase):

    def test_python_numbers_mix_with_any_signature(self):
        sc = Hy(1, 2, mu=1)
        for result in (sc + 3, 3 + sc, sc - 3, sc * 3, 3 * sc, sc / 2):
            self.assertEqual(result.signs, (1,))
        self.assertEqual(sc + 3, Hy(4, 2, mu=1))
        self.assertEqual(3 - sc, Hy(2, -2, mu=1))
        self.assertEqual(3 * sc, Hy(3, 6, mu=1))
        self.assertEqual(sc + Fraction(1, 2), Hy('3/2', 2, mu=1))
        self.assertEqual(sc + 0.5, Hy('3/2', 2, mu=1))
        self.assertEqual(sc + '1/2', Hy('3/2', 2, mu=1))

    def test_rational_division_by_a_scalar_keeps_the_signature(self):
        x = Hy(Hy(2, 4), Hy(6, 8), mu=1)
        self.assertEqual((x / 2).signs, (-1, 1))
        self.assertEqual(x / 2, Hy(Hy(1, 2), Hy(3, 4), mu=1))

    def test_scalar_divided_by_split_value(self):
        sc = Hy(2, 1, mu=1)
        self.assertEqual(1 / sc, sc.inverse())
        self.assertEqual((1 / sc).signs, (1,))

    def test_strings_are_read_in_the_other_operands_algebra(self):
        sq = Hy.from_array([1, 2, 3, 4], signs="split-quaternion")
        r = sq + '1+j'
        self.assertEqual(r.signs, (-1, 1))
        # the rank-1 text '1+j' embeds as 1 + i, exactly as in the
        # classical algebra (rank-1 'j' is the level-1 unit, 'i' at rank 2)
        self.assertEqual(r, Hy.from_array([2, 3, 3, 4], signs="split-quaternion"))
        self.assertEqual(sq, '1+2i+3j+4k')
        # a string of lower rank embeds; one of higher rank pads with -1
        sc = Hy(1, 2, mu=1)
        self.assertEqual((sc + '1+j').signs, (1,))
        r2 = sc + '1+i+j+k'       # rank-2 text read as signs (1, -1)
        self.assertEqual(r2.signs, (1, -1))

    def test_python_complex_literals_are_ordinary_complex(self):
        self.assertEqual(Hy(1, 2) + 1j, Hy(1, 3))
        with self.assertRaises(ValueError):
            Hy(1, 2, mu=1) + 1j
        self.assertNotEqual(Hy(1, 2, mu=1), 1 + 2j)
        self.assertEqual(Hy(1, 2), 1 + 2j)

    def test_incompatible_signatures_raise_on_arithmetic(self):
        a, b = Hy(1, 2, mu=1), Hy(3, 4)
        for op in (lambda: a + b, lambda: a - b, lambda: a * b, lambda: a / b,
                   lambda: add(a, b), lambda: mul(a, b)):
            with self.assertRaises(ValueError):
                op()

    def test_comparison_across_signatures_never_raises(self):
        self.assertFalse(Hy(1, 2, mu=1) == Hy(3, 4))
        self.assertTrue(Hy(1, 2, mu=1) != Hy(1, 2))

    def test_lower_rank_values_embed_into_a_prefix_compatible_algebra(self):
        sq = Hy.from_array([1, 2, 3, 4], signs="split-quaternion")
        z = Hy(5, 6)                              # complex, signs (-1,)
        total = sq + z
        self.assertEqual(total.signs, (-1, 1))
        self.assertEqual(total, Hy.from_array([6, 8, 3, 4], signs="split-quaternion"))
        prod = z * sq
        self.assertEqual(prod.signs, (-1, 1))
        # not a prefix: split-complex does not embed into (-1, 1)
        with self.assertRaises(ValueError):
            sq + Hy(5, 6, mu=1)

    def test_embedding_helper_accepts_signature_tuples(self):
        e = _embed(Hy(1, 2, mu=1), (Fraction(1), Fraction(-3)))
        self.assertEqual(e.signs, (1, -3))
        with self.assertRaises(ValueError):
            _embed(Hy(1, 2, mu=1), (Fraction(-1), Fraction(-1)))
        with self.assertRaises(ValueError):
            _embed(Hy(Hy(1, 2), Hy(3, 4)), (Fraction(-1),))
        # an int target gives any *new* levels the default mu = -1
        self.assertEqual(_embed(Hy(1, 2, mu=1), 2).signs, (1, -1))

    def test_sum_and_neg_keep_the_signature(self):
        x = Hy.random(signs="split-octonion", seed=5)
        self.assertEqual((-x).signs, x.signs)
        self.assertEqual((x + x).signs, x.signs)
        self.assertEqual((x - x), 0)
        self.assertEqual(x.conjugate().signs, x.signs)


class TestSignatureParsing(unittest.TestCase):

    def test_parse_default_is_classical(self):
        self.assertEqual(Hy.parse('1+2i+3j+4k').signs, (-1, -1))
        self.assertEqual(Hy('1+2j').signs, (-1,))

    def test_parse_with_signs(self):
        x = Hy.parse('1+2i+3j+4k', signs="split-quaternion")
        self.assertEqual(x.signs, (-1, 1))
        self.assertEqual(x, Hy.from_array([1, 2, 3, 4], signs="split-quaternion"))
        y = Hy.parse('(1+2j)', signs=(3,))
        self.assertEqual(y.signs, (3,))
        self.assertEqual(Hy.from_string('1+2j', signs=Hy.SPLIT_COMPLEX).signs, (1,))

    def test_parse_signs_length_must_match_rank_of_text(self):
        with self.assertRaises(ValueError):
            Hy.parse('1+2j', signs=(-1, 1))
        with self.assertRaises(ValueError):
            Hy.parse('1+2i+3j+4k', signs=(1,))

    def test_str_parse_round_trip_with_matching_signs(self):
        for n, sg in enumerate(SIGNATURES_UP_TO_RANK_3 + SIGNATURES_RANK_4):
            rng = random.Random(3000 + n)
            x = rand_signed(sg, rng)
            self.assertEqual(Hy.parse(str(x), signs=sg), x)

    def test_constructor_string_with_mu_disagreement_is_rejected(self):
        # Hy('1+2j') is read classically; asking for another mu is an
        # error rather than being silently ignored
        with self.assertRaises(ValueError):
            Hy('1+2j', mu=1)
        self.assertEqual(Hy('1+2j', mu=-1), Hy(1, 2))


class TestSignatureNormAndPredicates(unittest.TestCase):

    def test_norm_and_norm_squared_agree(self):
        for n, sg in enumerate(SIGNATURES_UP_TO_RANK_3):
            x = rand_signed(sg, random.Random(4000 + n))
            self.assertEqual(x.norm(), x.norm_squared())
            self.assertEqual(x.norm_squared(), abs2(x))

    def test_norm_squared_matches_sum_of_squares_when_classical(self):
        for rank in (1, 2, 3, 4):
            x = Hy.random(rank, seed=rank)
            self.assertEqual(x.norm_squared(),
                             sum(c * c for c in x.to_array()))

    def test_norm_squared_formula_by_signature(self):
        # N = sum over coordinates of c_i^2 * (-unit_square(i))
        for n, sg in enumerate(SIGNATURES_UP_TO_RANK_3 + SIGNATURES_RANK_4):
            x = rand_signed(sg, random.Random(4100 + n))
            expected = sum(c * c * -Hy.unit_square(i, sg) if i else c * c
                           for i, c in enumerate(x.to_array()))
            self.assertEqual(x.norm_squared(), expected)

    def test_abs_behaviour(self):
        self.assertAlmostEqual(abs(Hy(3, 4)), 5.0)             # unchanged
        self.assertAlmostEqual(abs(Hy(5, 3, mu=1)), 4.0)       # positive form
        self.assertEqual(abs(Hy(1, 1, mu=1)), 0.0)             # null
        with self.assertRaises(ValueError) as ctx:
            abs(Hy(1, 2, mu=1))                                # negative form
        self.assertIn('norm_squared', str(ctx.exception))

    def test_is_null(self):
        self.assertTrue(Hy(1, 1, mu=1).is_null())
        self.assertFalse(Hy(0, 0, mu=1).is_null())             # zero is not null
        self.assertFalse(Hy(2, 1, mu=1).is_null())
        for rank in (1, 2, 3):                                 # never, classically
            self.assertFalse(Hy.random(rank, seed=rank).is_null())

    def test_zero_divisors_without_null_norm_exist_from_rank_four(self):
        # the documented limitation: is_null() only sees the norm cone
        x = Hy.from_array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0])
        y = Hy.from_array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1])
        self.assertEqual(x * y, 0)
        self.assertFalse(x.is_null())
        self.assertFalse(y.is_null())

    def test_negative_powers_use_the_signature_aware_inverse(self):
        x = Hy.from_array([2, 1, 0, 1], signs="split-quaternion")
        self.assertEqual(x ** -1, x.inverse())
        self.assertEqual(x ** -3, (x * x * x).inverse())
        self.assertEqual(x ** 0, _embed_signs_for_test(1, (-1, 1)))
        null = Hy.from_array([1, 0, 1, 0], signs="split-quaternion")
        self.assertTrue(null.is_null())
        with self.assertRaises(ZeroDivisionError):
            null ** -2

    def test_pow_zero_keeps_signature(self):
        self.assertEqual((Hy(2, 3, mu=1) ** 0).signs, (1,))


@unittest.skipUnless(HAVE_NUMPY, "numpy not installed")
class TestSignatureMatrixRepresentation(unittest.TestCase):

    @staticmethod
    def _det_2x2(m):
        return m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]

    def test_split_complex_is_a_homomorphism(self):
        rng = random.Random(5000)
        for _ in range(FUZZ_TRIALS):
            x, y = rand_signed((1,), rng), rand_signed((1,), rng)
            self.assertTrue(np.array_equal(
                np.dot(x.to_matrix(), y.to_matrix()), (x * y).to_matrix()))

    def test_split_complex_matrix_is_symmetric_hyperbolic_form(self):
        m = Hy(2, 3, mu=1).to_matrix()
        self.assertEqual(m.tolist(), [[2, 3], [3, 2]])
        self.assertEqual(self._det_2x2(m), Hy(2, 3, mu=1).norm_squared())

    def test_general_mu_matrix_entries(self):
        m = Hy(2, 3, mu=5).to_matrix()
        self.assertEqual(m.tolist(), [[2, 15], [3, 2]])
        self.assertEqual(self._det_2x2(m), Hy(2, 3, mu=5).norm_squared())

    def test_split_quaternion_homomorphism_and_determinant(self):
        rng = random.Random(5001)
        for _ in range(FUZZ_TRIALS):
            x, y = (rand_signed((-1, 1), rng) for _ in range(2))
            self.assertTrue(np.array_equal(
                np.dot(x.to_matrix(), y.to_matrix()), (x * y).to_matrix()))
        if HAVE_SYMPY:
            x = rand_signed((-1, 1), rng)
            det = sympy.Matrix(x.to_matrix().tolist()).det()
            self.assertEqual(det, x.norm_squared() ** 2)

    def test_left_and_right_representations_with_signs(self):
        rng = random.Random(5002)
        for kind in ('left', 'right'):
            x, y = rand_signed((-1, 1), rng), rand_signed((-1, 1), rng)
            lhs = np.dot(x.to_matrix(kind=kind), y.to_matrix(kind=kind))
            rhs = ((x * y) if kind == 'left' else (y * x)).to_matrix(kind=kind)
            self.assertTrue(np.array_equal(lhs, rhs))

    def test_from_matrix_with_signs_round_trips(self):
        rng = random.Random(5003)
        for sg in ((1,), (-1, 1), (1, 1), (2,), (-1, 3)):
            x = rand_signed(sg, rng)
            back = Hy.from_matrix(x.to_matrix(), signs=sg, validate=True)
            self.assertEqual(back, x)
            self.assertEqual(back.signs, x.signs)

    def test_from_matrix_with_the_wrong_signs_reads_a_different_algebra(self):
        x = Hy(2, 3, mu=1)
        back = Hy.from_matrix(x.to_matrix())          # default signs
        self.assertEqual(back.to_array(), x.to_array())
        self.assertNotEqual(back, x)
        # validation catches the mismatch: the split-complex matrix is not
        # the regular representation of the *classical* complex number
        with self.assertRaises(ValueError):
            Hy.from_matrix(x.to_matrix(), validate=True)

    def test_rank_three_guard_still_applies_to_split_octonions(self):
        x = Hy.random(signs="split-octonion", seed=3)
        with self.assertRaises(ValueError):
            x.to_matrix()
        m = x.to_matrix(allow_nonassociative=True)
        self.assertEqual(m.shape, (8, 8))


class TestSignatureInterop(unittest.TestCase):

    def test_conversion_guard_rejects_non_classical_signatures(self):
        from hyprat.hypercomplex import _quaternion_fractions
        for sg in ((1,), (-1, 1), (1, 1), (-1, 2), (3,)):
            with self.assertRaises(ValueError):
                _quaternion_fractions(Hy.random(signs=sg, seed=1), "target")
        # ordinary complex / quaternion values still convert
        self.assertEqual(len(_quaternion_fractions(Hy(1, 2), "target")), 4)
        self.assertEqual(len(_quaternion_fractions(Hy.random(2, seed=1), "t")), 4)

    @unittest.skipUnless(HAVE_SYMPY, "sympy not installed")
    def test_to_sympy_rejects_split_quaternions(self):
        with self.assertRaises(ValueError):
            Hy.from_array([1, 2, 3, 4], signs="split-quaternion").to_sympy()
        with self.assertRaises(ValueError):
            Hy(1, 2, mu=1).to_sympy()

    @unittest.skipUnless(HAVE_SYMPY, "sympy not installed")
    def test_classical_conversion_is_unaffected(self):
        x = Hy.from_array(['1/2', 2, 3, 4])
        self.assertEqual(Hy.from_sympy(x.to_sympy()), x)
        self.assertEqual(Hy.from_sympy(x.to_sympy()).signs, (-1, -1))


class TestSignatureHelpers(unittest.TestCase):
    """Module-level functions and miscellany under signatures."""

    def test_module_functions_accept_signed_values(self):
        a, b = Hy(2, 3, mu=1), Hy(5, 7, mu=1)
        self.assertEqual(mul(a, b), a * b)
        self.assertEqual(add(a, b), a + b)
        self.assertEqual(sub(a, b), a - b)
        self.assertEqual(neg(a), -a)
        self.assertEqual(conj(a), a.conjugate())
        self.assertEqual(abs2(a), a.norm_squared())
        self.assertEqual(div(a, b), a / b)
        self.assertEqual(inverse(a), a.inverse())

    def test_module_functions_on_raw_fractions_ignore_signatures(self):
        self.assertEqual(mul(Fraction(2), Fraction(3)), 6)
        self.assertEqual(abs2(Fraction(-3)), 9)

    def test_raw_fraction_embeds_into_any_signature(self):
        r = mul(Fraction(3), Hy(1, 2, mu=1))
        self.assertEqual(r, Hy(3, 6, mu=1))
        self.assertEqual(r.signs, (1,))

    def test_hash_and_iter_protocols_still_work(self):
        x = Hy(1, 2, mu=1)
        self.assertEqual(list(x), [Fraction(1), Fraction(2)])
        self.assertEqual(len(x), 2)
        self.assertEqual(x[0], 1)

    def test_immutability_covers_mu(self):
        x = Hy(1, 2, mu=1)
        with self.assertRaises(AttributeError):
            x._mu = -1
        with self.assertRaises(AttributeError):
            x.mu = -1

    def test_latex_is_positional_and_unchanged(self):
        self.assertEqual(Hy(1, 2, mu=1).latex(), '1+2j')
        x = Hy.from_array([1, 2, 3, 4], signs="split-quaternion")
        self.assertEqual(x.latex(), '1+2i+3j+4k')

    def test_units_dict_keys_unchanged_by_signs(self):
        for rank in (1, 2, 3, 4):
            sg = (-1,) * (rank - 1) + (1,)
            self.assertEqual(list(Hy.units(signs=sg)), list(Hy.units(rank)))
        self.assertEqual(Hy.units(0), {'1': Fraction(1), '-1': Fraction(-1)})

    def test_units_have_the_requested_signature(self):
        for v in Hy.units(signs=(1, -1, 2)).values():
            self.assertEqual(v.signs, (1, -1, 2))



def main():
    unittest.main()


if __name__ == '__main__':
    main()
