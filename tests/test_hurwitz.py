"""
test_hurwitz.py
===============

Unit tests for hurwitz.py: the integer-backed ``Hu`` class holding the
Hurwitz integers (quaternions whose four coordinates are all integers or
all halves of odd integers).

Organization:

    TestConstruction        - Hu(...) from ints/Fractions/floats/strings,
                              single-argument forms, from_doubled(),
                              validation and error cases
    TestCoordinates         - .doubled .coords .a .b .c .d, iteration,
                              indexing, is_lipschitz()
    TestStrings             - str()/repr() (and the repr round trip),
                              parse()/from_string(), the quaternion-label
                              convention vs Hy.parse's rank-1 reading
    TestConversions         - to_hy()/from_hy() round trips, what from_hy
                              accepts and rejects, agreement with
                              Hy.is_hurwitz()
    TestArithmeticFixed     - hand-checked products and the quaternion
                              identities
    TestArithmeticFuzz      - seeded random checks of + - * against Hy used
                              as the oracle, plus the algebraic laws
    TestPower               - ** with non-negative and (unit) negative
                              exponents
    TestNormTraceConjugate  - norm(), norm_squared(), trace(), conjugate()
    TestUnits               - the 24 units: units(), is_unit(), group
                              structure, brute-force enumeration
    TestEqualityHashImmutability
                            - ==, hash, no equality with Hy, immutability,
                              copy/pickle, bool, subclassing
    TestUnsupportedOperations
                            - mixing with floats, Fractions, Hy; no division

Algebraic laws are checked empirically via seeded random fuzzing
(deterministic across runs), in addition to fixed hand-verified examples.
"""

import copy
import math
import pickle
import random
import unittest
from fractions import Fraction

from hyprat import Hu, Hy


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------

def rand_hu(rng, bound=8):
    """A random Hurwitz integer: half the time Lipschitz, half the time not."""
    if rng.random() < 0.5:
        return Hu(*[rng.randint(-bound, bound) for _ in range(4)])
    return Hu.from_doubled(*[2 * rng.randint(-bound, bound) + 1 for _ in range(4)])


def as_hy(x):
    """The oracle's version of a Hu."""
    return Hy.from_array(list(x.coords))


def same_parity(x):
    return len({v & 1 for v in x.doubled}) == 1


class _MyHu(Hu):
    """A subclass at module level, so that it can be pickled."""


I, J, K = Hu(0, 1, 0, 0), Hu(0, 0, 1, 0), Hu(0, 0, 0, 1)
ONE = Hu(1)
OMEGA = Hu("1/2", "1/2", "1/2", "1/2")


# ============================================================================
# construction
# ============================================================================

class TestConstruction(unittest.TestCase):
    def test_four_integers(self):
        x = Hu(1, 2, 3, 4)
        self.assertEqual(x.doubled, (2, 4, 6, 8))

    def test_missing_coordinates_are_zero(self):
        self.assertEqual(Hu().doubled, (0, 0, 0, 0))
        self.assertEqual(Hu(5).doubled, (10, 0, 0, 0))
        self.assertEqual(Hu(1, 2).doubled, (2, 4, 0, 0))
        self.assertEqual(Hu(1, 2, 3).doubled, (2, 4, 6, 0))
        self.assertEqual(Hu(a=1, d=-1).doubled, (2, 0, 0, -2))

    def test_half_integers_from_strings_fractions_and_floats(self):
        want = (1, -1, 3, -3)
        for coords in [("1/2", "-1/2", "3/2", "-3/2"),
                       (Fraction(1, 2), Fraction(-1, 2), Fraction(3, 2), Fraction(-3, 2)),
                       (0.5, -0.5, 1.5, -1.5),
                       ("0.5", "-.5" if False else "-0.5", "1.5", "-1.5"),
                       (" 1/2 ", "-1/2", "3/2", "-3/2")]:
            self.assertEqual(Hu(*coords).doubled, want, coords)

    def test_mixed_input_types(self):
        self.assertEqual(Hu(1, "2", Fraction(3), 4.0).doubled, (2, 4, 6, 8))
        self.assertEqual(Hu("1/2", Fraction(-1, 2), 0.5, "3/2").doubled, (1, -1, 1, 3))

    def test_integral_looking_fractions_are_integers(self):
        self.assertEqual(Hu(Fraction(4, 2), "6/3", 0, 0).doubled, (4, 4, 0, 0))

    def test_bools_and_numpy_style_integers_are_accepted(self):
        self.assertEqual(Hu(True, False).doubled, (2, 0, 0, 0))

    def test_mixed_parity_is_rejected(self):
        for coords in [("1/2", 0, 0, 0), ("1/2", "1/2", 0, 0), (1, 1, 1, "1/2"),
                       ("1/2", "1/2", "1/2", 1), (0, 0, 0, "3/2")]:
            with self.assertRaises(ValueError, msg=coords):
                Hu(*coords)

    def test_other_denominators_are_rejected(self):
        for coords in [("1/3", "1/3", "1/3", "1/3"), ("1/4", 0, 0, 0),
                       (1, 2, 3, "2/3"), (0.1, 0, 0, 0), ("0.25", 0, 0, 0)]:
            with self.assertRaises(ValueError, msg=coords):
                Hu(*coords)

    def test_error_message_names_the_coordinates(self):
        with self.assertRaises(ValueError) as cm:
            Hu("1/2", 0, 0, 0)
        self.assertIn("(1/2, 0, 0, 0)", str(cm.exception))
        self.assertIn("all integers or all halves of odd integers", str(cm.exception))

    def test_unreadable_values(self):
        for bad in ("abc", "", "1/0", "1/2/3"):
            with self.assertRaises(ValueError, msg=repr(bad)):
                Hu(1, bad, 0, 0)
        for bad in (float("nan"), float("inf"), float("-inf")):
            with self.assertRaises(ValueError, msg=repr(bad)):
                Hu(bad, 0, 0, 0)

    def test_wrong_types(self):
        for bad in ([1], (1, 2), {1: 2}, 1 + 2j, object()):
            with self.assertRaises(TypeError, msg=repr(bad)):
                Hu(1, bad, 0, 0)

    def test_single_argument_forms(self):
        x = Hu(1, 2, 3, 4)
        self.assertIs(Hu(x), x)
        self.assertEqual(Hu("(1+2i+3j+4k)"), x)
        self.assertEqual(Hu(Hy(Hy(1, 2), Hy(3, 4))), x)
        self.assertEqual(Hu(7), Hu(7, 0, 0, 0))

    def test_from_doubled(self):
        self.assertEqual(Hu.from_doubled(2, 4, 6, 8), Hu(1, 2, 3, 4))
        self.assertEqual(Hu.from_doubled(1, 1, 1, 1), OMEGA)
        self.assertEqual(Hu.from_doubled(-1, 3, -5, 7).doubled, (-1, 3, -5, 7))

    def test_from_doubled_validates(self):
        for bad in [(1, 0, 0, 0), (1, 1, 0, 0), (2, 2, 2, 1), (0, 1, 1, 1)]:
            with self.assertRaises(ValueError, msg=bad):
                Hu.from_doubled(*bad)
        for bad in [(1.0, 1, 1, 1), ("1", 1, 1, 1), (True, 1, 1, 1),
                    (Fraction(1), 1, 1, 1)]:
            with self.assertRaises(TypeError, msg=bad):
                Hu.from_doubled(*bad)

    def test_huge_values_are_exact(self):
        big = 10 ** 60
        x = Hu(big, -big, big + 1, 7)
        self.assertEqual(x.coords, (big, -big, big + 1, 7))
        y = Hu.from_doubled(2 * big + 1, 1, 1, -1)
        self.assertEqual(y.a, Fraction(2 * big + 1, 2))


# ============================================================================
# coordinates
# ============================================================================

class TestCoordinates(unittest.TestCase):
    def test_properties(self):
        x = Hu("1/2", "-3/2", "5/2", "7/2")
        self.assertEqual(x.doubled, (1, -3, 5, 7))
        self.assertEqual(x.coords, (Fraction(1, 2), Fraction(-3, 2),
                                    Fraction(5, 2), Fraction(7, 2)))
        self.assertEqual((x.a, x.b, x.c, x.d), x.coords)

    def test_coordinates_are_fractions_even_when_integral(self):
        x = Hu(1, 2, 3, 4)
        for c in x.coords:
            self.assertIsInstance(c, Fraction)
        self.assertIsInstance(x.doubled[0], int)

    def test_iteration_unpacking_and_indexing(self):
        x = Hu(1, 2, 3, 4)
        a, b, c, d = x
        self.assertEqual((a, b, c, d), (1, 2, 3, 4))
        self.assertEqual(list(x), [1, 2, 3, 4])
        self.assertEqual(x[0], 1)
        self.assertEqual(x[-1], 4)
        self.assertEqual(x[1:3], (2, 3))
        with self.assertRaises(IndexError):
            x[4]

    def test_is_lipschitz(self):
        self.assertTrue(Hu(1, 2, 3, 4).is_lipschitz())
        self.assertTrue(Hu(0).is_lipschitz())
        self.assertFalse(OMEGA.is_lipschitz())
        self.assertFalse(Hu("-3/2", "1/2", "1/2", "5/2").is_lipschitz())


# ============================================================================
# strings
# ============================================================================

class TestStrings(unittest.TestCase):
    def test_str_known_values(self):
        self.assertEqual(str(Hu(1, 2, 3, 4)), "(1+2i+3j+4k)")
        self.assertEqual(str(Hu(0)), "(0)")
        self.assertEqual(str(Hu(5)), "(5)")
        self.assertEqual(str(I), "(i)")
        self.assertEqual(str(-J), "(-j)")
        self.assertEqual(str(Hu(-1, 1, 0, -1)), "(-1+i-k)")
        self.assertEqual(str(OMEGA), "(1/2+1/2i+1/2j+1/2k)")
        self.assertEqual(str(Hu("1/2", "-1/2", "-3/2", "1/2")), "(1/2-1/2i-3/2j+1/2k)")

    def test_str_matches_the_hy_oracle(self):
        rng = random.Random(3)
        for _ in range(200):
            x = rand_hu(rng)
            self.assertEqual(str(x), str(as_hy(x)))

    def test_repr_known_values(self):
        self.assertEqual(repr(Hu(1, 2, 3, 4)), "Hu(1, 2, 3, 4)")
        self.assertEqual(repr(Hu(0)), "Hu(0, 0, 0, 0)")
        self.assertEqual(repr(Hu(-1, 0, 5, -2)), "Hu(-1, 0, 5, -2)")
        self.assertEqual(repr(OMEGA), "Hu('1/2', '1/2', '1/2', '1/2')")
        self.assertEqual(repr(Hu("-3/2", "1/2", "-1/2", "5/2")),
                         "Hu('-3/2', '1/2', '-1/2', '5/2')")

    def test_repr_round_trips(self):
        rng = random.Random(4)
        for _ in range(200):
            x = rand_hu(rng, bound=30)
            self.assertEqual(eval(repr(x), {"Hu": Hu}), x)

    def test_parse_known_forms(self):
        cases = {
            "(1+2i+3j+4k)": Hu(1, 2, 3, 4),
            "1+2i+3j+4k": Hu(1, 2, 3, 4),
            " ( 1 + 2i + 3j + 4k ) ": Hu(1, 2, 3, 4),
            "(1/2+1/2i+1/2j+1/2k)": OMEGA,
            "(1/2-1/2i+3/2j-5/2k)": Hu("1/2", "-1/2", "3/2", "-5/2"),
            "j": J, "-j": -J, "+k": K, "i": I, "-i": -I,
            "(i)": I, "(-k)": -K,
            "5": Hu(5), "-3": Hu(-3), "(0)": Hu(0),
            "0.5+0.5i+0.5j+0.5k": OMEGA,
            "1-i-j-k": Hu(1, -1, -1, -1),
            "2i+3i": Hu(0, 5, 0, 0),             # repeated units add up
            "i+i": Hu(0, 2, 0, 0),
        }
        for text, want in cases.items():
            self.assertEqual(Hu.parse(text), want, text)
            self.assertEqual(Hu(text), want, text)
            self.assertEqual(Hu.from_string(text), want, text)

    def test_parse_round_trips_str(self):
        rng = random.Random(5)
        for _ in range(300):
            x = rand_hu(rng, bound=40)
            self.assertEqual(Hu.parse(str(x)), x)

    def test_parse_rejects_malformed_text(self):
        for bad in ["", "   ", "()", "+", "-", "i+", "1++i", "1+-i", "ii", "1i2",
                    "1+2l", "1+2e1", "1+L", "abc", "1 2", "(1+2i", "1+2i)",
                    "1/i", "1.5.5+0i+0j+0k", "i j", "1+k+*", "((1))"]:
            with self.assertRaises(ValueError, msg=repr(bad)):
                Hu.parse(bad)

    def test_parse_rejects_non_hurwitz_values(self):
        for bad in ["1/2", "1/2+1/2i", "1/2+1/2i+1/2j+k", "1/3+i", "0.25+0i+0j+0k",
                    "1/2+1/2i+1/2j"]:
            with self.assertRaises(ValueError, msg=repr(bad)):
                Hu.parse(bad)

    def test_text_uses_quaternion_labels_unlike_hy_parse(self):
        # '(3j)' alone is a rank-1 value for Hy.parse, but the quaternion 3j for Hu
        self.assertEqual(Hu("(3j)"), Hu(0, 0, 3, 0))
        self.assertEqual(Hu.from_hy(Hy.parse("(3j)")), Hu(0, 3, 0, 0))
        self.assertEqual(Hu("(2-3j)"), Hu(2, 0, -3, 0))
        self.assertEqual(Hu.from_hy(Hy.parse("(2-3j)")), Hu(2, -3, 0, 0))

    def test_str_of_every_unit_round_trips_through_text(self):
        for key, u in Hu.units().items():
            self.assertEqual(Hu.parse(key), u, key)
            self.assertEqual(Hu.parse(str(u)), u, key)


# ============================================================================
# conversions to and from Hy
# ============================================================================

class TestConversions(unittest.TestCase):
    def test_to_hy_known(self):
        self.assertEqual(Hu(1, 2, 3, 4).to_hy(), Hy(Hy(1, 2), Hy(3, 4)))
        self.assertEqual(OMEGA.to_hy(), Hy(Hy("1/2", "1/2"), Hy("1/2", "1/2")))
        self.assertEqual(Hu(0).to_hy(), Hy(Hy(0, 0), Hy(0, 0)))

    def test_to_hy_shape(self):
        h = Hu(1, 2, 3, 4).to_hy()
        self.assertIsInstance(h, Hy)
        self.assertEqual(h.rank, 2)
        self.assertEqual(h.signs, (Fraction(-1), Fraction(-1)))

    def test_round_trip_hu_hy_hu(self):
        rng = random.Random(6)
        for _ in range(300):
            x = rand_hu(rng, bound=50)
            self.assertEqual(Hu.from_hy(x.to_hy()), x)
            self.assertEqual(Hu(x.to_hy()), x)

    def test_round_trip_hy_hu_hy(self):
        Hy.seed(7)
        got = 0
        for _ in range(400):
            h = Hy.random(2, lo=-4, hi=4, dmax=2)
            if h.is_hurwitz():
                got += 1
                self.assertEqual(Hu.from_hy(h).to_hy(), h)
        self.assertGreater(got, 20)             # the fuzz really exercised it

    def test_from_hy_rank_one_means_a_plus_bi(self):
        self.assertEqual(Hu.from_hy(Hy(2, -3)), Hu(2, -3, 0, 0))
        self.assertEqual(Hu.from_hy(Hy(0, 1)), I)
        self.assertEqual(Hu.from_hy(Hy(0, 1)) ** 2, Hu(-1))

    def test_from_hy_higher_rank_inside_the_quaternions(self):
        q = Hy(Hy(1, 2), Hy(3, 4))
        for rank in (2, 3, 4):
            self.assertEqual(Hu.from_hy(q.embed(rank)), Hu(1, 2, 3, 4))
        self.assertEqual(Hu.from_hy(OMEGA.to_hy().embed(4)), OMEGA)

    def test_from_hy_rejections(self):
        bad = [
            Hy(Hy("1/2", 0), Hy(0, 0)),                    # lone half
            Hy(Hy("1/2", "1/2"), Hy(0, 0)),                # two halves
            Hy(Hy("1/2", "1/2"), Hy("1/2", 1)),            # mixed
            Hy(Hy("1/3", "1/3"), Hy("1/3", "1/3")),        # thirds
            Hy("1/2", "1/2"),                              # rank 1, two halves
            Hy(1, 2, mu=1),                                # split-complex
            Hy.from_array([1, 2, 3, 4], signs=(-1, 1)),    # split-quaternion
            Hy.from_array([1, 2, 3, 4, 0, 0, 0, 1]),       # octonion with L part
            Hy.from_array([1, 2, 3, 4] + [0] * 11 + [1]),  # sedenion with e15
        ]
        for h in bad:
            with self.assertRaises(ValueError, msg=repr(h)):
                Hu.from_hy(h)

    def test_from_hy_error_messages_say_why(self):
        with self.assertRaises(ValueError) as cm:
            Hu.from_hy(Hy(1, 2, mu=1))
        self.assertIn("mu = -1", str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            Hu.from_hy(Hy.from_array([1, 2, 3, 4, 0, 0, 0, 1]))
        self.assertIn("beyond the first four", str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            Hu.from_hy(Hy(Hy("1/2", 0), Hy(0, 0)))
        self.assertIn("all integers or all halves of odd integers", str(cm.exception))

    def test_from_hy_rejects_non_hy(self):
        for bad in (1, "1+2i", Fraction(1), 1.5, None, [1, 2, 3, 4], 1 + 2j):
            with self.assertRaises(TypeError, msg=repr(bad)):
                Hu.from_hy(bad)

    def test_from_hy_succeeds_exactly_when_is_hurwitz(self):
        Hy.seed(11)
        rng = random.Random(11)
        samples = []
        for rank in (1, 2, 3):
            samples += [Hy.random(rank, lo=-3, hi=3, dmax=2) for _ in range(60)]
        # plus plenty of Hurwitz-looking values at several ranks / signatures
        for _ in range(60):
            c = rand_hu(rng, 3).coords
            samples.append(Hy.from_array(list(c)).embed(rng.choice([2, 3, 4])))
            samples.append(Hy.from_array(list(c), signs=(-1, 1)))
        samples += [Hy(1, 2, mu=1), Hy(3, 4), Hy("1/2", "1/2")]
        hits = 0
        for h in samples:
            if h.is_hurwitz():
                hits += 1
                Hu.from_hy(h)
            else:
                with self.assertRaises(ValueError, msg=repr(h)):
                    Hu.from_hy(h)
        self.assertGreater(hits, 60)

    def test_embedded_gaussian_integers(self):
        Hy.seed(12)
        for _ in range(40):
            g = Hy.random(1, lo=-9, hi=9, dmax=1)
            x = Hu.from_hy(g)
            self.assertEqual(x.doubled, (2 * int(g.real), 2 * int(g.imag), 0, 0))
            self.assertEqual(x.to_hy(), g.embed(2))


# ============================================================================
# arithmetic, hand-checked
# ============================================================================

class TestArithmeticFixed(unittest.TestCase):
    def test_quaternion_unit_identities(self):
        m1 = Hu(-1)
        self.assertEqual(I * I, m1)
        self.assertEqual(J * J, m1)
        self.assertEqual(K * K, m1)
        self.assertEqual(I * J * K, m1)
        self.assertEqual(I * J, K)
        self.assertEqual(J * K, I)
        self.assertEqual(K * I, J)
        self.assertEqual(J * I, -K)
        self.assertEqual(K * J, -I)
        self.assertEqual(I * K, -J)

    def test_known_product(self):
        self.assertEqual(Hu(1, 2, 3, 4) * Hu(5, 6, 7, 8), Hu(-60, 12, 30, 24))
        self.assertEqual(Hu(5, 6, 7, 8) * Hu(1, 2, 3, 4), Hu(-60, 20, 14, 32))

    def test_multiplication_is_not_commutative(self):
        self.assertNotEqual(Hu(1, 2, 3, 4) * Hu(5, 6, 7, 8),
                            Hu(5, 6, 7, 8) * Hu(1, 2, 3, 4))

    def test_half_integer_products(self):
        self.assertEqual(OMEGA * OMEGA, Hu("-1/2", "1/2", "1/2", "1/2"))
        self.assertEqual(OMEGA * OMEGA * OMEGA, Hu(-1))
        self.assertEqual(OMEGA ** 6, Hu(1))
        self.assertEqual(OMEGA * OMEGA.conjugate(), Hu(1))
        # (1+i+j+k)(1-i+j-k)/4, worked out by hand: (2-2i+2j+2k)/4
        x = Hu("1/2", "1/2", "1/2", "1/2") * Hu("1/2", "-1/2", "1/2", "-1/2")
        self.assertEqual(x, Hu("1/2", "-1/2", "1/2", "1/2"))
        # a half-integer value times its conjugate is a Lipschitz integer
        self.assertTrue((OMEGA * OMEGA.conjugate()).is_lipschitz())

    def test_addition_and_subtraction(self):
        self.assertEqual(Hu(1, 2, 3, 4) + Hu(5, 6, 7, 8), Hu(6, 8, 10, 12))
        self.assertEqual(Hu(1, 2, 3, 4) - Hu(5, 6, 7, 8), Hu(-4, -4, -4, -4))
        self.assertEqual(OMEGA + OMEGA, Hu(1, 1, 1, 1))
        self.assertEqual(OMEGA - OMEGA, Hu(0))
        # lipschitz + half-integer stays half-integer
        self.assertEqual(Hu(1, 0, 0, 0) + OMEGA, Hu("3/2", "1/2", "1/2", "1/2"))
        # two half-integer values sum to a Lipschitz integer
        self.assertTrue((OMEGA + Hu("1/2", "-1/2", "1/2", "-1/2")).is_lipschitz())

    def test_negation_and_identity(self):
        x = Hu(1, -2, 3, "-0")
        self.assertEqual(-x, Hu(-1, 2, -3, 0))
        self.assertEqual(+x, x)
        self.assertEqual(x + (-x), Hu(0))
        self.assertEqual(x * Hu(1), x)
        self.assertEqual(Hu(1) * x, x)
        self.assertEqual(x + Hu(0), x)
        self.assertEqual(x * Hu(0), Hu(0))

    def test_integer_operands(self):
        x = Hu(1, 2, 3, 4)
        self.assertEqual(x + 1, Hu(2, 2, 3, 4))
        self.assertEqual(1 + x, Hu(2, 2, 3, 4))
        self.assertEqual(x - 1, Hu(0, 2, 3, 4))
        self.assertEqual(1 - x, Hu(0, -2, -3, -4))
        self.assertEqual(x * 3, Hu(3, 6, 9, 12))
        self.assertEqual(3 * x, Hu(3, 6, 9, 12))
        self.assertEqual(OMEGA * 3, Hu("3/2", "3/2", "3/2", "3/2"))
        self.assertEqual(x * -1, -x)
        self.assertEqual(OMEGA + 1, Hu("3/2", "1/2", "1/2", "1/2"))

    def test_integral_fraction_operands(self):
        x = Hu(1, 2, 3, 4)
        self.assertEqual(x + Fraction(2), Hu(3, 2, 3, 4))
        self.assertEqual(Fraction(2) * x, Hu(2, 4, 6, 8))

    def test_results_are_hu(self):
        for r in (OMEGA + 1, OMEGA * OMEGA, 2 * OMEGA, OMEGA - 1, -OMEGA, OMEGA ** 2):
            self.assertIsInstance(r, Hu)
            self.assertTrue(same_parity(r))


# ============================================================================
# arithmetic, fuzzed against Hy
# ============================================================================

class TestArithmeticFuzz(unittest.TestCase):
    def setUp(self):
        self.rng = random.Random(20261009)

    def test_matches_the_hy_oracle(self):
        for _ in range(500):
            x, y = rand_hu(self.rng), rand_hu(self.rng)
            hx, hy = as_hy(x), as_hy(y)
            self.assertEqual(as_hy(x + y), hx + hy)
            self.assertEqual(as_hy(x - y), hx - hy)
            self.assertEqual(as_hy(x * y), hx * hy)
            self.assertEqual(as_hy(-x), -hx)
            self.assertEqual(as_hy(x.conjugate()), hx.conjugate())
            self.assertEqual(x.norm(), hx.norm_squared())
            self.assertEqual(x.norm(), hx.norm())
            self.assertEqual(x.trace(), 2 * hx.components()[0])

    def test_results_stay_in_the_ring(self):
        for _ in range(500):
            x, y = rand_hu(self.rng), rand_hu(self.rng)
            for r in (x + y, x - y, x * y, -x, x.conjugate()):
                self.assertTrue(same_parity(r), (x, y, r))
                self.assertTrue(as_hy(r).is_hurwitz())

    def test_associativity_and_distributivity(self):
        for _ in range(300):
            x, y, z = (rand_hu(self.rng) for _ in range(3))
            self.assertEqual((x * y) * z, x * (y * z))
            self.assertEqual(x * (y + z), x * y + x * z)
            self.assertEqual((x + y) * z, x * z + y * z)
            self.assertEqual(x + y, y + x)
            self.assertEqual((x + y) + z, x + (y + z))

    def test_conjugation_laws(self):
        for _ in range(300):
            x, y = rand_hu(self.rng), rand_hu(self.rng)
            self.assertEqual((x * y).conjugate(), y.conjugate() * x.conjugate())
            self.assertEqual((x + y).conjugate(), x.conjugate() + y.conjugate())
            self.assertEqual(x.conjugate().conjugate(), x)
            self.assertEqual(x * x.conjugate(), Hu(x.norm()))
            self.assertEqual(x.conjugate() * x, Hu(x.norm()))
            self.assertEqual(x + x.conjugate(), Hu(x.trace()))

    def test_norm_is_multiplicative(self):
        for _ in range(500):
            x, y = rand_hu(self.rng), rand_hu(self.rng)
            self.assertEqual((x * y).norm(), x.norm() * y.norm())

    def test_integer_scalars_are_central(self):
        for _ in range(200):
            x = rand_hu(self.rng)
            n = self.rng.randint(-9, 9)
            y = rand_hu(self.rng)
            self.assertEqual(Hu(n) * x, x * Hu(n))
            self.assertEqual((n * x) * y, x * (n * y))
            self.assertEqual(n * x, Hu(n) * x)
            self.assertEqual(as_hy(n * x), n * as_hy(x))

    def test_no_zero_divisors(self):
        for _ in range(300):
            x, y = rand_hu(self.rng), rand_hu(self.rng)
            if x and y:
                self.assertTrue(x * y)

    def test_large_coordinates(self):
        big = 10 ** 40
        for _ in range(60):
            x = Hu(*[self.rng.randint(-big, big) for _ in range(4)])
            y = Hu.from_doubled(*[2 * self.rng.randint(-big, big) + 1 for _ in range(4)])
            self.assertEqual(as_hy(x * y), as_hy(x) * as_hy(y))
            self.assertEqual((x * y).norm(), x.norm() * y.norm())


# ============================================================================
# powers
# ============================================================================

class TestPower(unittest.TestCase):
    def test_small_powers_match_repeated_multiplication(self):
        rng = random.Random(8)
        for _ in range(100):
            x = rand_hu(rng, 4)
            acc = Hu(1)
            for n in range(0, 7):
                self.assertEqual(x ** n, acc, (x, n))
                acc = acc * x

    def test_power_matches_the_oracle(self):
        rng = random.Random(9)
        for _ in range(60):
            x = rand_hu(rng, 3)
            for n in (0, 1, 2, 3, 5, 8):
                self.assertEqual(as_hy(x ** n), as_hy(x) ** n)

    def test_zero_exponent_and_zero_base(self):
        self.assertEqual(Hu(1, 2, 3, 4) ** 0, Hu(1))
        self.assertEqual(Hu(0) ** 0, Hu(1))
        self.assertEqual(Hu(0) ** 5, Hu(0))
        self.assertEqual(OMEGA ** 0, Hu(1))

    def test_exponent_laws(self):
        rng = random.Random(10)
        for _ in range(60):
            x = rand_hu(rng, 3)
            m, n = rng.randint(0, 6), rng.randint(0, 6)
            self.assertEqual(x ** (m + n), x ** m * x ** n)
            self.assertEqual(x ** (m * n), (x ** m) ** n)

    def test_negative_powers_of_units(self):
        for u in Hu.units().values():
            self.assertEqual(u ** -1, u.conjugate())
            self.assertEqual(u ** -1 * u, Hu(1))
            self.assertEqual(u ** -3, (u ** 3).conjugate())
            self.assertEqual(u ** -2 * u ** 2, Hu(1))

    def test_negative_power_of_a_non_unit_raises(self):
        for x in (Hu(2), Hu(1, 1, 0, 0), Hu(0), OMEGA * 2, Hu(1, 2, 3, 4)):
            with self.assertRaises(ValueError, msg=repr(x)):
                x ** -1
        with self.assertRaises(ValueError):
            Hu(1, 1, 0, 0) ** -3

    def test_huge_exponents_of_units_are_fast(self):
        self.assertEqual(OMEGA ** (6 * 10 ** 6), Hu(1))
        self.assertEqual(OMEGA ** (6 * 10 ** 6 + 3), Hu(-1))
        self.assertEqual(I ** (4 * 10 ** 30 + 1), I)
        self.assertEqual(OMEGA ** -(6 * 10 ** 6 + 1), OMEGA.conjugate())

    def test_unsupported_exponents(self):
        for bad in (0.5, 2.0, Fraction(1, 2), "2", None, Hu(2)):
            with self.assertRaises(TypeError, msg=repr(bad)):
                Hu(1, 2, 3, 4) ** bad
        with self.assertRaises(TypeError):
            pow(Hu(1, 2, 3, 4), 2, 5)


# ============================================================================
# norm, trace, conjugate
# ============================================================================

class TestNormTraceConjugate(unittest.TestCase):
    def test_norm_known(self):
        self.assertEqual(Hu(0).norm(), 0)
        self.assertEqual(Hu(1).norm(), 1)
        self.assertEqual(Hu(1, 2, 3, 4).norm(), 30)
        self.assertEqual(Hu(-1, -1, -1, -1).norm(), 4)
        self.assertEqual(OMEGA.norm(), 1)
        self.assertEqual(Hu("3/2", "1/2", "1/2", "1/2").norm(), 3)
        self.assertEqual(Hu("3/2", "3/2", "3/2", "3/2").norm(), 9)
        self.assertEqual(Hu(0, 3, 4, 0).norm(), 25)

    def test_norm_is_always_an_int(self):
        rng = random.Random(14)
        for _ in range(300):
            x = rand_hu(rng, 20)
            self.assertIsInstance(x.norm(), int)
            self.assertEqual(x.norm(), sum(c * c for c in x.coords))
            self.assertEqual(x.norm() * 4, sum(v * v for v in x.doubled))

    def test_norm_squared_is_the_same_thing(self):
        x = Hu(1, 2, 3, 4)
        self.assertEqual(x.norm_squared(), x.norm())
        self.assertEqual(OMEGA.norm_squared(), 1)

    def test_norm_is_zero_only_for_zero(self):
        self.assertEqual(Hu(0).norm(), 0)
        rng = random.Random(15)
        for _ in range(200):
            x = rand_hu(rng)
            self.assertEqual(x.norm() == 0, not x)

    def test_trace(self):
        self.assertEqual(Hu(1, 2, 3, 4).trace(), 2)
        self.assertEqual(OMEGA.trace(), 1)
        self.assertEqual(Hu("-3/2", "1/2", "1/2", "1/2").trace(), -3)
        self.assertEqual(Hu(0, 1, 2, 3).trace(), 0)
        self.assertIsInstance(OMEGA.trace(), int)

    def test_trace_and_norm_give_the_characteristic_polynomial(self):
        rng = random.Random(16)
        for _ in range(100):
            x = rand_hu(rng, 5)
            self.assertEqual(x * x - x.trace() * x + Hu(x.norm()), Hu(0))

    def test_conjugate_known(self):
        self.assertEqual(Hu(1, 2, 3, 4).conjugate(), Hu(1, -2, -3, -4))
        self.assertEqual(OMEGA.conjugate(), Hu("1/2", "-1/2", "-1/2", "-1/2"))
        self.assertEqual(Hu(5).conjugate(), Hu(5))
        self.assertEqual(Hu(0, 1, 2, 3).conjugate(), -Hu(0, 1, 2, 3))


# ============================================================================
# units
# ============================================================================

class TestUnits(unittest.TestCase):
    def test_there_are_24_distinct_units(self):
        units = Hu.units()
        self.assertEqual(len(units), 24)
        self.assertEqual(len(set(units.values())), 24)
        self.assertEqual(len(set(u.doubled for u in units.values())), 24)

    def test_order_and_keys(self):
        keys = list(Hu.units())
        self.assertEqual(keys[:8], ["1", "-1", "i", "-i", "j", "-j", "k", "-k"])
        self.assertEqual(keys[8], "1/2+1/2i+1/2j+1/2k")
        self.assertEqual(keys[-1], "-1/2-1/2i-1/2j-1/2k")
        self.assertIn("1/2-1/2i+1/2j-1/2k", keys)
        for key, u in Hu.units().items():
            self.assertEqual(f"({key})", str(u))

    def test_values_are_hu_of_norm_one(self):
        for u in Hu.units().values():
            self.assertIsInstance(u, Hu)
            self.assertEqual(u.norm(), 1)
            self.assertTrue(u.is_unit())

    def test_eight_lipschitz_units_and_sixteen_others(self):
        us = list(Hu.units().values())
        self.assertEqual(sum(u.is_lipschitz() for u in us), 8)
        self.assertEqual(sum(not u.is_lipschitz() for u in us), 16)

    def test_matches_brute_force_enumeration(self):
        found = set()
        for v in _all_doubled_vectors(range(-2, 3)):
            if len({x & 1 for x in v}) == 1 and sum(x * x for x in v) == 4:
                found.add(v)
        self.assertEqual(found, {u.doubled for u in Hu.units().values()})

    def test_closed_under_multiplication_and_a_group(self):
        units = list(Hu.units().values())
        table = set(units)
        for u in units:
            for v in units:
                self.assertIn(u * v, table)
            self.assertIn(u.conjugate(), table)
            self.assertEqual(u * u.conjugate(), ONE)
            self.assertIn(ONE * u, table)

    def test_every_unit_has_a_finite_order_dividing_the_group(self):
        orders = {}
        for u in Hu.units().values():
            n, p = 1, u
            while p != ONE:
                p, n = p * u, n + 1
            orders[n] = orders.get(n, 0) + 1
        # binary tetrahedral group: 1 of order 1, 1 of order 2, 6 of order 4,
        # 8 of order 3, 8 of order 6
        self.assertEqual(orders, {1: 1, 2: 1, 4: 6, 3: 8, 6: 8})

    def test_is_unit_false_elsewhere(self):
        for x in (Hu(0), Hu(2), Hu(1, 1, 0, 0), Hu(1, 1, 1, 1), Hu(-3),
                  Hu("3/2", "1/2", "1/2", "1/2")):
            self.assertFalse(x.is_unit(), x)

    def test_is_unit_agrees_with_norm_one_on_random_values(self):
        rng = random.Random(17)
        for _ in range(300):
            x = rand_hu(rng, 1)
            self.assertEqual(x.is_unit(), x.norm() == 1)
            if x.is_unit():
                self.assertIn(x, set(Hu.units().values()))

    def test_units_returns_a_fresh_dict(self):
        d = Hu.units()
        d.clear()
        self.assertEqual(len(Hu.units()), 24)

    def test_oracle_agrees_on_unit_multiplication(self):
        for u in Hu.units().values():
            for v in Hu.units().values():
                self.assertEqual(as_hy(u * v), as_hy(u) * as_hy(v))

    def test_hy_units_of_rank_two_are_the_lipschitz_units(self):
        hy_units = set(Hu.from_hy(h) for h in Hy.units(2).values())
        lip = {u for u in Hu.units().values() if u.is_lipschitz()}
        self.assertEqual(hy_units, lip)


def _all_doubled_vectors(values):
    values = list(values)
    for a in values:
        for b in values:
            for c in values:
                for d in values:
                    yield (a, b, c, d)


# ============================================================================
# equality, hash, immutability, copying, subclassing
# ============================================================================

class TestEqualityHashImmutability(unittest.TestCase):
    def test_equality_between_hu(self):
        self.assertEqual(Hu(1, 2, 3, 4), Hu(1, 2, 3, 4))
        self.assertNotEqual(Hu(1, 2, 3, 4), Hu(1, 2, 3, 5))
        self.assertEqual(OMEGA, Hu("0.5", "0.5", "0.5", "0.5"))
        self.assertNotEqual(OMEGA, -OMEGA)
        self.assertTrue(Hu(1, 2, 3, 4) == Hu(1, 2, 3, 4))
        self.assertFalse(Hu(1, 2, 3, 4) != Hu(1, 2, 3, 4))

    def test_equality_with_scalars(self):
        self.assertEqual(Hu(5), 5)
        self.assertEqual(5, Hu(5))
        self.assertEqual(Hu(0), 0)
        self.assertEqual(Hu(-3), Fraction(-3))
        self.assertNotEqual(Hu(5), 6)
        self.assertNotEqual(Hu(5, 1, 0, 0), 5)
        self.assertNotEqual(OMEGA, Fraction(1, 2))
        self.assertNotEqual(Hu(5), Fraction(11, 2))

    def test_never_equal_to_hy_or_other_types(self):
        for x in (Hu(1, 2, 3, 4), Hu(0), Hu(5), OMEGA):
            self.assertFalse(x == x.to_hy())
            self.assertFalse(x.to_hy() == x)
            self.assertTrue(x != x.to_hy())
            self.assertTrue(x.to_hy() != x)
        for other in ("(1+2i+3j+4k)", None, [1, 2, 3, 4], (1, 2, 3, 4), 1.0, 1 + 0j):
            self.assertFalse(Hu(1) == other, repr(other))

    def test_hash_is_consistent_with_equality(self):
        rng = random.Random(18)
        for _ in range(200):
            x = rand_hu(rng)
            y = Hu.from_doubled(*x.doubled)
            self.assertEqual(hash(x), hash(y))
        self.assertEqual(hash(Hu(5)), hash(5))
        self.assertEqual(hash(Hu(-7)), hash(-7))
        self.assertEqual(hash(Hu(0)), hash(0))
        self.assertEqual(hash(Hu(5)), hash(Fraction(5)))

    def test_usable_in_sets_and_dicts(self):
        s = {Hu(1, 2, 3, 4), Hu(1, 2, 3, 4), OMEGA, Hu("1/2", "1/2", "1/2", "1/2")}
        self.assertEqual(len(s), 2)
        d = {Hu(5): "five"}
        self.assertEqual(d[5], "five")
        self.assertEqual(len(set(Hu.units().values())), 24)

    def test_bool(self):
        self.assertFalse(Hu(0))
        self.assertFalse(Hu())
        for x in (Hu(1), Hu(0, 0, 0, 1), OMEGA, Hu(-1, 0, 0, 0)):
            self.assertTrue(x)

    def test_immutability(self):
        x = Hu(1, 2, 3, 4)
        for name in ("_a", "_b", "a", "doubled", "foo"):
            with self.assertRaises(AttributeError, msg=name):
                setattr(x, name, 5)
        with self.assertRaises(AttributeError):
            del x._a
        self.assertEqual(x, Hu(1, 2, 3, 4))
        with self.assertRaises(AttributeError):
            x.new_attribute = 1

    def test_no_instance_dict(self):
        self.assertFalse(hasattr(Hu(1), "__dict__"))

    def test_copy_deepcopy_pickle(self):
        for x in (Hu(1, 2, 3, 4), OMEGA, Hu(0), Hu(-5, 7, 0, 1)):
            self.assertEqual(copy.copy(x), x)
            self.assertEqual(copy.deepcopy(x), x)
            for proto in range(0, pickle.HIGHEST_PROTOCOL + 1):
                y = pickle.loads(pickle.dumps(x, proto))
                self.assertEqual(y, x)
                self.assertIsInstance(y, Hu)
                self.assertTrue(same_parity(y))

    def test_subclass_instances_keep_their_type(self):
        MyHu = _MyHu
        x = MyHu(1, 2, 3, 4)
        self.assertIsInstance(x, MyHu)
        self.assertIsInstance(x + x, MyHu)
        self.assertIsInstance(x * x, MyHu)
        self.assertIsInstance(-x, MyHu)
        self.assertIsInstance(x.conjugate(), MyHu)
        self.assertIsInstance(x ** 3, MyHu)
        self.assertIsInstance(MyHu.parse("1+i"), MyHu)
        self.assertEqual(repr(x), "_MyHu(1, 2, 3, 4)")
        self.assertEqual(x, Hu(1, 2, 3, 4))
        self.assertEqual(pickle.loads(pickle.dumps(x)).__class__, MyHu)


# ============================================================================
# things that are deliberately not supported
# ============================================================================

class TestUnsupportedOperations(unittest.TestCase):
    def test_no_division(self):
        x, y = Hu(1, 2, 3, 4), Hu(1, 1, 0, 0)
        for op in ("__truediv__", "__floordiv__", "__mod__", "__divmod__"):
            self.assertFalse(hasattr(Hu, op), op)
        with self.assertRaises(TypeError):
            x / y
        with self.assertRaises(TypeError):
            x // y
        with self.assertRaises(TypeError):
            x % y
        with self.assertRaises(TypeError):
            1 / x

    def test_floats_and_non_integral_fractions_are_not_operands(self):
        x = Hu(1, 2, 3, 4)
        for bad in (1.5, 2.0, Fraction(1, 2), 1 + 2j, "1", None, [1], (1, 2, 3, 4)):
            with self.assertRaises(TypeError, msg=repr(bad)):
                x + bad
            with self.assertRaises(TypeError, msg=repr(bad)):
                bad + x
            with self.assertRaises(TypeError, msg=repr(bad)):
                x - bad
            with self.assertRaises(TypeError, msg=repr(bad)):
                x * bad
            with self.assertRaises(TypeError, msg=repr(bad)):
                bad * x

    def test_hy_is_not_an_operand(self):
        x = Hu(1, 2, 3, 4)
        h = x.to_hy()
        for op in (lambda a, b: a + b, lambda a, b: a - b, lambda a, b: a * b):
            with self.assertRaises(TypeError):
                op(x, h)
            with self.assertRaises(TypeError):
                op(h, x)

    def test_no_ordering(self):
        with self.assertRaises(TypeError):
            Hu(1) < Hu(2)
        with self.assertRaises(TypeError):
            Hu(1) >= 0

    def test_no_int_float_complex_conversion(self):
        for conv in (int, float, complex):
            with self.assertRaises(TypeError):
                conv(Hu(1))


def main():
    unittest.main()


if __name__ == '__main__':
    main()
