Split-Hypercomplex Numbers and Interoperability
===============================================

This page continues ```hyprat``: Exact Rational Hypercomplex
Numbers <classic_hy.html>`__, which introduced the class ``Hy`` and the
classical hypercomplex numbers. It has two parts. The first extends
``Hy`` to the **split-hypercomplex numbers**, and to *signatures* in
general. The second shows how ``Hy`` values are converted to and from
the numbers of other packages: SymPy, numpy-quaternion, quaternionic,
and ``gint``.

The imports below are the ones this page uses. ``numpy``, ``sympy``,
``numpy-quaternion`` (imported as ``quaternion``) and ``quaternionic``
come from the optional ``interop`` extra, and ``gint`` is installed from
GitHub; see `Optional Packages and the ``interop``
Extra <classic_hy.html#optional-packages-and-the-interop-extra>`__ for
how to install them.

.. code:: ipython3

    >>> import math
    >>> from fractions import Fraction
    >>> import numpy as np
    >>> import sympy
    >>> import quaternion        # the numpy-quaternion package
    >>> import quaternionic
    >>> from gint import Zi, Qi
    >>> from hyprat import Hy

Split-Hypercomplex Numbers
--------------------------

The split-hypercomplex numbers come from the same doubling construction
as the classical ones, but with :math:`+1` in place of :math:`-1` in one
or more of the doubling steps. They are represented by the same class,
``Hy``. The difference is in what the products of the basis units are,
and with that in the quadratic form, which is no longer a sum of
squares. That has consequences: a nonzero value can have norm zero, and
zero divisors appear already at rank 1.

Signatures
~~~~~~~~~~

Each doubling step of the Cayley-Dickson construction takes an algebra
:math:`A` and builds the pairs :math:`(a, b)` with :math:`a, b \in A`.
In ``hyprat`` the product of two pairs carries a parameter :math:`\mu`,
which can be any nonzero rational number:

.. math::


   (a, b)(c, d) = (ac + \mu\,\overline{d}\,b,\; da + b\,\overline{c})

Addition and conjugation are the same as before,
:math:`(a, b) + (c, d) = (a + c,\, b + d)` and
:math:`\overline{(a, b)} = (\overline{a}, -b)`. With :math:`\mu = -1`
the product is exactly the classical one given in `What is the
Cayley-Dickson
construction? <classic_hy.html#what-is-the-cayley-dickson-construction>`__.
With :math:`\mu = +1` it is the *split* version of that doubling step.
The two products differ only in a sign, and the underlying vector space
is the same.

**One :math:`\mu` per doubling step.** A rank-:math:`n` value is the
result of :math:`n` doublings, one at each level of the tower, and each
doubling has its own :math:`\mu`. Level 0 is the step from the rationals
to rank 1, level 1 is the step from rank 1 to rank 2, and so on. The
list of all the :math:`\mu`\ ’s, from the lowest level to the highest,
is the **signature** of the algebra. In ``hyprat`` the signature of a
value is ``some_hy.signs``, and the :math:`\mu` of its top level is
``some_hy.mu``. Since the signature determines the multiplication,
values with different signatures belong to different algebras (see
`Mixing, Equality, and Parsing <#mixing-equality-and-parsing>`__).

The following table lists the signatures that have names. “Split” refers
to the *top* level only: the lower levels stay classical.

+-----------------+----------------------+-------------------+------------------------+-------------------------+
| Rank            | Signature            | Algebra           | Preset name            | Constant                |
+=================+======================+===================+========================+=========================+
| 1               | ``(-1,)``            | complex numbers   |                        |                         |
+-----------------+----------------------+-------------------+------------------------+-------------------------+
| 1               | ``(1,)``             | split-complex     | ``"split-complex"``    | ``Hy.SPLIT_COMPLEX``    |
|                 |                      | numbers           |                        |                         |
+-----------------+----------------------+-------------------+------------------------+-------------------------+
| 2               | ``(-1, -1)``         | quaternions       |                        |                         |
+-----------------+----------------------+-------------------+------------------------+-------------------------+
| 2               | ``(-1, 1)``          | split-quaternions | ``"split-quaternion"`` | ``Hy.SPLIT_QUATERNION`` |
+-----------------+----------------------+-------------------+------------------------+-------------------------+
| 3               | ``(-1, -1, -1)``     | octonions         |                        |                         |
+-----------------+----------------------+-------------------+------------------------+-------------------------+
| 3               | ``(-1, -1, 1)``      | split-octonions   | ``"split-octonion"``   | ``Hy.SPLIT_OCTONION``   |
+-----------------+----------------------+-------------------+------------------------+-------------------------+
| 4               | ``(-1, -1, -1, -1)`` | sedenions         |                        |                         |
+-----------------+----------------------+-------------------+------------------------+-------------------------+

Nothing restricts you to these. Any sequence of nonzero rational numbers
is a signature, for example ``(1, -1)``, ``(-1, 1, 1, -1)`` or
``(2, -1)``, and ``Hy`` computes exactly in the resulting algebra. Such
algebras need not have the pleasant properties of the named ones.

**Reading a signature from a value.** A value built with
``Hy(real, imag, mu=...)`` gets ``mu`` at its top level, and inherits
the lower levels from its components, which must agree with each other.
The default, ``mu=-1``, builds the classical numbers.

.. code:: ipython3

    >>> Hy(0, 1).signs, Hy(0, 1, mu=1).signs         # rank 1, classical and split




.. parsed-literal::

    ((Fraction(-1, 1),), (Fraction(1, 1),))



.. code:: ipython3

    >>> Hy(Hy(1, 2), Hy(3, 4)).signs                  # rank 2, classical




.. parsed-literal::

    (Fraction(-1, 1), Fraction(-1, 1))



.. code:: ipython3

    >>> Hy(Hy(1, 2), Hy(3, 4), mu=1).signs            # rank 2: classical lower level, split top level




.. parsed-literal::

    (Fraction(-1, 1), Fraction(1, 1))



**Choosing a signature when building values.** The factory methods
``Hy.units``, ``Hy.random``, ``Hy.from_array``, ``Hy.parse`` and
``Hy.from_matrix`` take a ``signs=`` argument, either a tuple or the
name of a preset. The common cases have named presets, usable as
``signs="split-quaternion"`` or through the class constants:

.. code:: ipython3

    >>> Hy.SPLIT_COMPLEX, Hy.SPLIT_QUATERNION, Hy.SPLIT_OCTONION




.. parsed-literal::

    ((1,), (-1, 1), (-1, -1, 1))



.. code:: ipython3

    >>> by_hand = Hy(Hy(1, 2), Hy(3, 4), mu=1)
    >>> from_list = Hy.from_array([1, 2, 3, 4], signs="split-quaternion")
    >>> by_hand == from_list, by_hand.signs




.. parsed-literal::

    (True, (Fraction(-1, 1), Fraction(1, 1)))



**What the signature changes.** Two things depend on the signature.

The first is *what each basis unit squares to*. The coordinates of a
rank-:math:`n` value are numbered :math:`0, 1, \dots, 2^n - 1`, in the
order of ``to_array()``, and coordinate :math:`0` is the real unit
:math:`1`. The binary digits of a coordinate’s number say which doubling
levels its unit involves: the unit at coordinate :math:`3`, binary
``11``, involves levels :math:`0` and :math:`1`. If :math:`S` is the set
of levels that a unit involves, then it squares to

.. math::


   -\prod_{l \in S} (-\mu_l).

In the classical case every :math:`-\mu_l` is :math:`1`, so every
imaginary unit squares to :math:`-1`. With :math:`\mu_l = \pm 1`, each
split level in :math:`S` contributes a factor of :math:`-1`, so a unit
that involves an even number of split levels squares to :math:`-1`, and
a unit that involves an odd number squares to :math:`+1`. The following
cell checks the formula against ``Hy.unit_square``, for a rank-4
signature chosen to be unusual:

.. code:: ipython3

    >>> def predicted_square(n, signs):
    ...     "-prod(-mu) over the doubling levels that the binary digits of n involve."
    ...     if n == 0:
    ...         return 1
    ...     return -math.prod(-s for level, s in enumerate(signs) if n >> level & 1)
    
    >>> unusual = (1, -1, 1, -1)
    >>> print([int(Hy.unit_square(n, unusual)) for n in range(16)])
    >>> print([predicted_square(n, unusual) for n in range(16)])
    >>> print(all(Hy.unit_square(n, unusual) == predicted_square(n, unusual) for n in range(16)))


.. parsed-literal::

    [1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1]
    [1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1]
    True


The second is the *quadratic form* :math:`N(x) = x\,\overline{x}`, which
``norm_squared()`` computes exactly (``norm()`` is the same thing under
an older name). For a pair it is

.. math::


   N(a, b) = N(a) - \mu\,N(b),

with :math:`N(r) = r^2` for a rational :math:`r`. In the classical
algebras this is a sum of squares. When some :math:`\mu` is positive the
form is *indefinite*: it can be negative, or zero for a nonzero value.
For rank 1 the formula says that :math:`N(a + bj) = a^2 - \mu b^2`, and
for higher ranks it applies recursively:

.. code:: ipython3

    >>> [Hy(3, 2, mu=mu).norm_squared() for mu in (-1, 1, 2)]       # 3**2 - mu * 2**2




.. parsed-literal::

    [Fraction(13, 1), Fraction(5, 1), Fraction(1, 1)]



.. code:: ipython3

    >>> x = Hy(Hy(1, 2), Hy(3, 4), mu=1)                             # signature (-1, 1)
    >>> x.norm_squared(), Hy(1, 2).norm_squared() - Hy(3, 4).norm_squared()    # N(a) - N(b)




.. parsed-literal::

    (Fraction(-20, 1), Fraction(-20, 1))



The sections below look at the split algebras one rank at a time, then
at norms, mixing and parsing, zero divisors and matrices.

Split-Complex Numbers
~~~~~~~~~~~~~~~~~~~~~

A single step with ``mu = +1`` gives the split-complex numbers, where
``j * j`` is ``+1`` instead of ``-1``.

.. code:: ipython3

    >>> js = Hy(0, 1, mu=1)
    >>> js * js




.. parsed-literal::

    Hy('1', '0', mu=1)



``repr`` shows ``mu`` only when it is not the default ``-1``, so the
output above can be pasted back in to rebuild the value. ``str`` is the
same in every algebra.

The quadratic form is no longer a sum of squares. It is *indefinite*, so
a nonzero value can have ``norm_squared()`` equal to zero. Such a value
is *null*, and is a zero divisor:

.. code:: ipython3

    >>> n_plus, n_minus = Hy(1, 1, mu=1), Hy(1, -1, mu=1)
    >>> n_plus * n_minus, n_plus.is_null(), n_plus.norm_squared()




.. parsed-literal::

    (Hy('0', '0', mu=1), True, Fraction(0, 1))



Null elements have no inverse, and asking for one raises
``ZeroDivisionError``:

.. code:: ipython3

    >>> try:
    ...     n_plus.inverse()
    ... except ZeroDivisionError as e:
    ...     print(e)


.. parsed-literal::

    hypercomplex value has zero norm (it is a nonzero null element of a split algebra); not invertible


Split-Quaternions
~~~~~~~~~~~~~~~~~

The split-quaternions are ``signs=(-1, 1)``: ``i`` squares to ``-1``,
while ``j`` and ``k`` square to ``+1``.

.. code:: ipython3

    >>> su = Hy.units(signs="split-quaternion")
    >>> {name: str(su[name] * su[name]) for name in ("i", "j", "k")}




.. parsed-literal::

    {'i': '(-1)', 'j': '(1)', 'k': '(1)'}



.. code:: ipython3

    >>> si, sj, sk = su["i"], su["j"], su["k"]
    >>> si * sj == sk, sj * si == -sk




.. parsed-literal::

    (True, True)



The split-quaternions are isomorphic to the :math:`2 \times 2` real
matrices; see `Matrices and Other
Packages <#matrices-and-other-packages>`__.

Split-Octonions
~~~~~~~~~~~~~~~

The unit **labels** are positional. They name a place in the
Cayley-Dickson tower (so ``i * j == k`` and ``iL == i * L`` hold in
every signature) and do not change with ``mu``. What changes is what
each unit squares to. For the split-octonions, ``signs=(-1, -1, 1)``,
``i``, ``j``, ``k`` square to ``-1`` and the four units involving the
top level, ``L``, ``iL``, ``jL``, ``kL``, square to ``+1``:

.. code:: ipython3

    >>> suo = Hy.units(signs="split-octonion")
    >>> {name: str(suo[name] * suo[name]) for name in ("i", "j", "k", "L", "iL", "jL", "kL")}




.. parsed-literal::

    {'i': '(-1)',
     'j': '(-1)',
     'k': '(-1)',
     'L': '(1)',
     'iL': '(1)',
     'jL': '(1)',
     'kL': '(1)'}



``Hy.unit_square(index, signs)`` gives the same answer without building
anything. The index is the position of the unit in ``to_array()`` order,
so ``0`` is the real unit ``1``:

.. code:: ipython3

    >>> [int(Hy.unit_square(n, "split-octonion")) for n in range(8)]




.. parsed-literal::

    [1, -1, -1, -1, 1, 1, 1, 1]



Any Nonzero Rational ``mu``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``mu`` is not limited to ``+1`` and ``-1``. With ``mu = 2``, ``j * j``
is ``2``, and because ``2`` is not a rational square the form
``a**2 - 2*b**2`` never vanishes on a nonzero value, so every nonzero
value is invertible:

.. code:: ipython3

    >>> xg = Hy(2, 1, mu=2)
    >>> xg.norm_squared(), xg.inverse(), xg * xg.inverse()




.. parsed-literal::

    (Fraction(2, 1), Hy('1', '-1/2', mu=2), Hy('1', '0', mu=2))



Norms
~~~~~

With a positive ``mu`` somewhere the quadratic form is *indefinite*:
``norm_squared()`` (also still available as ``norm()``) can be negative
or zero for a nonzero value. A nonzero value with
``norm_squared() == 0`` is *null*: a zero divisor with no inverse.
``abs()`` is the square root of ``norm_squared()``, so it raises
``ValueError`` when the form is negative.

``abs()`` is the square root of ``norm_squared()``, so it raises
``ValueError`` when the form is negative and there is no real square
root:

.. code:: ipython3

    >>> Hy(5, 3, mu=1).norm_squared(), abs(Hy(5, 3, mu=1))




.. parsed-literal::

    (Fraction(16, 1), 4.0)



.. code:: ipython3

    >>> try:
    ...     abs(Hy(3, 5, mu=1))
    ... except ValueError as e:
    ...     print(e)


.. parsed-literal::

    norm_squared() is negative (-16), so abs() is undefined; use norm_squared() for the indefinite quadratic form


The quadratic form is multiplicative, :math:`N(xy) = N(x)\,N(y)`, for
ranks 0 to 3, the composition algebras, split or not. From rank 4 up it
is not. For example, the classical sedenions ``e1 + e10`` and
``e4 - e15`` multiply to zero, although neither has norm zero:

.. code:: ipython3

    >>> sa, sb = Hy.random(signs="split-quaternion", seed=3), Hy.random(signs="split-quaternion", seed=4)
    >>> (sa * sb).norm_squared() == sa.norm_squared() * sb.norm_squared()




.. parsed-literal::

    True



.. code:: ipython3

    >>> so, sp = Hy.random(signs="split-octonion", seed=5), Hy.random(signs="split-octonion", seed=6)
    >>> (so * sp).norm_squared() == so.norm_squared() * sp.norm_squared()




.. parsed-literal::

    True



.. code:: ipython3

    >>> e1_e10 = Hy.from_array([0, 1] + [0] * 8 + [1] + [0] * 5)               # e1 + e10
    >>> e4_e15 = Hy.from_array([0] * 4 + [1] + [0] * 10 + [-1])                # e4 - e15
    >>> (e1_e10 * e4_e15).norm_squared(), e1_e10.norm_squared() * e4_e15.norm_squared()




.. parsed-literal::

    (Fraction(0, 1), Fraction(4, 1))



Mixing, Equality, and Parsing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Values of different signatures never combine: arithmetic raises
``ValueError`` and ``==`` is just ``False``, while plain numbers mix
freely. ``str`` is the same in every algebra and ``repr`` shows ``mu``
only when it is not ``-1``; ``Hy.parse`` needs ``signs=`` to say which
algebra a string lives in. The regular-representation matrices still
work (they are genuine homomorphisms at rank 1-2), and are read back
with the matching ``signs=``.

.. code:: ipython3

    >>> Hy(1, 2, mu=1) + 3




.. parsed-literal::

    Hy('4', '2', mu=1)



.. code:: ipython3

    >>> try:
    ...     Hy(1, 2, mu=1) + Hy(1, 2)
    ... except ValueError as e:
    ...     print(e)


.. parsed-literal::

    incompatible signatures: cannot combine a value with signs ('1',) and one with signs ('-1',)


.. code:: ipython3

    >>> Hy(1, 2, mu=1) == Hy(1, 2)




.. parsed-literal::

    False



Since text such as ``'1+2i+3j+4k'`` does not say which algebra it lives
in, ``Hy.parse`` (like ``Hy.from_array``, ``Hy.units``, ``Hy.random``
and ``Hy.from_matrix``) takes a ``signs=`` argument:

.. code:: ipython3

    >>> Hy.parse('1+2i+3j+4k', signs="split-quaternion")




.. parsed-literal::

    Hy(Hy('1', '2'), Hy('3', '4'), mu=1)



Zero Divisors
~~~~~~~~~~~~~

In the split algebras, zero divisors appear much earlier than in the
classical tower: the example ``(1 + j)(1 - j) = 0`` in `Split-Complex
Numbers <#split-complex-numbers>`__ is already at rank 1. A nonzero
value with ``norm_squared() == 0`` is *null*, and is a zero divisor with
no inverse. For ranks 1 to 3 the zero divisors are exactly the null
elements, which is what ``is_null()`` detects, and
``is_zero_divisor()``, the exact test described in `Zero
Divisors <classic_hy.html#zero-divisors>`__, agrees. From rank 4 up
there are also zero divisors with nonzero norm.

.. code:: ipython3

    >>> n_plus.is_null(), n_plus.is_zero_divisor()                                       # (1+j) in the split-complex numbers




.. parsed-literal::

    (True, True)



.. code:: ipython3

    >>> Hy(1, 2, mu=1).is_null(), Hy(1, 2, mu=1).is_zero_divisor()                       # 1 - 4 = -3: not null, invertible




.. parsed-literal::

    (False, False)



``annihilator()`` also works in every signature. A value is a zero
divisor on one side exactly when it is on the other, and the left and
right annihilators always have the same dimension. They coincide in the
classical sedenion example in `Zero
Divisors <classic_hy.html#zero-divisors>`__, but in a split algebra they
often do not:

.. code:: ipython3

    >>> zd_split = Hy.from_array([1, 0, 0, 0, 1] + [0] * 11, signs=(-1, -1, 1, -1))
    >>> left, right = zd_split.annihilator("left"), zd_split.annihilator("right")
    >>> len(left), len(right), left == right




.. parsed-literal::

    (8, 8, False)



Matrices and Other Packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The split-quaternions are isomorphic to the 2x2 real matrices, and
``to_matrix()`` still gives a genuine homomorphism at ranks 1 and 2 for
any signature. A matrix does not record its signature, so read it back
with the matching ``signs=``:

.. code:: ipython3

    >>> sa, sb = Hy.random(signs="split-quaternion", seed=10), Hy.random(signs="split-quaternion", seed=11)
    >>> np.array_equal(np.dot(sa.to_matrix(), sb.to_matrix()), (sa * sb).to_matrix())




.. parsed-literal::

    True



.. code:: ipython3

    >>> Hy.from_matrix(sa.to_matrix(), signs="split-quaternion") == sa




.. parsed-literal::

    True



The conversions to the other quaternion packages, and ``complex()``,
only make sense for the classical algebras, so they reject anything
else:

.. code:: ipython3

    >>> try:
    ...     sa.to_sympy()
    ... except ValueError as e:
    ...     print(e)


.. parsed-literal::

    cannot convert a Hy with signs ('-1', '1') to a sympy Quaternion; only ordinary complex/quaternion values (every mu equal to -1) can be converted


Interoperability with Other Packages
------------------------------------

``Hy`` values can be converted to and from the numbers of several other
packages. The conversions are explicit methods, never implicit
coercions: a ``Hy`` is never ``==`` to a value of another package, so
mixing them in arithmetic fails loudly rather than silently picking an
algebra.

Overview and Optional Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following table lists the packages, what ``Hy`` converts to and from
each of them, and whether the conversion is exact.

+------------------+------------------+-------------+--------------------------------+------------------+
| Package          | Import name      | Converts    | ``Hy`` methods                 | Exact?           |
+==================+==================+=============+================================+==================+
| SymPy            | ``sympy``        | rank 1 and  | ``to_sympy()``,                | yes              |
|                  |                  | 2           | ``Hy.from_sympy()``            |                  |
+------------------+------------------+-------------+--------------------------------+------------------+
| numpy-quaternion | ``quaternion``   | rank 1 and  | ``to_numpy_quaternion()``,     | no: ``float64``  |
|                  |                  | 2           | ``Hy.from_numpy_quaternion()`` |                  |
+------------------+------------------+-------------+--------------------------------+------------------+
| quaternionic     | ``quaternionic`` | rank 1 and  | ``to_quaternionic()``,         | no: ``float64``  |
|                  |                  | 2           | ``Hy.from_quaternionic()``     |                  |
+------------------+------------------+-------------+--------------------------------+------------------+
| ``gint``         | ``gint``         | rank 1, the | ``Hy.from_gint()``,            | yes              |
|                  |                  | Gaussian    | ``to_qi()``, ``to_zi()``,      |                  |
|                  |                  | integers    | ``is_gaussian()``              |                  |
|                  |                  | and         |                                |                  |
|                  |                  | rationals   |                                |                  |
+------------------+------------------+-------------+--------------------------------+------------------+
| NumPy            | ``numpy``        | any rank,   | ``to_matrix()``,               | yes (an array of |
|                  |                  | as a matrix | ``Hy.from_matrix()``           | ``Fraction``\ s) |
+------------------+------------------+-------------+--------------------------------+------------------+

The matrix conversions were described in `Matrix
Representation <classic_hy.html#matrix-representation>`__. This section
is about the other four packages.

**Optional dependencies.** None of these packages is required by
``hyprat``. Each is imported only when a conversion method that needs it
is called, so install whichever you use. The first three come from the
``interop`` extra, as explained in `Optional Packages and the
``interop``
Extra <classic_hy.html#optional-packages-and-the-interop-extra>`__:

::

   pip install sympy numpy-quaternion quaternionic

``gint`` is installed from GitHub (see `gint (Gaussian Integers and
Rationals) <#gint-gaussian-integers-and-rationals>`__).

**The quaternion packages.** All three quaternion packages use the same
Hamilton convention as ``Hy`` (``i*j == k``), and their coordinates
always correspond to ``1, i, j, k``, in that order, so the ``from_*``
methods need no reordering. A rank-1 ``Hy`` (a rational complex number)
is embedded as the quaternion ``a + b*i`` when it is converted. Ranks 3
and up (octonions, sedenions, …) have no quaternion equivalent, so
converting them raises a ``ValueError``. So does a value of a split
algebra: only the classical algebras, those in which every :math:`\mu`
is :math:`-1`, can be converted.

Each ``from_*`` method converts a *single* quaternion. To convert
several at once, such as an array of ``numpy-quaternion`` or
``quaternionic`` quaternions, convert them one at a time, for example
``[Hy.from_numpy_quaternion(q) for q in array]``.

Here are the quaternion and the rank-1 value that the following sections
use, and the two cases that cannot be converted:

.. code:: ipython3

    >>> q = Hy(Hy(1, 2), Hy(3, 4))          # the quaternion 1 + 2i + 3j + 4k
    >>> z = Hy('5/2', '-16/5')              # a rank-1 value, embedded as a + b*i
    >>> repr(z.to_numpy_quaternion())




.. parsed-literal::

    'quaternion(2.5, -3.2, 0, 0)'



.. code:: ipython3

    >>> try:
    ...     Hy.random(3, seed=1).to_numpy_quaternion()
    ... except ValueError as e:
    ...     print(e)


.. parsed-literal::

    cannot convert a rank-3 Hy (8 coordinates) to a numpy quaternion; only rank 1 (embedded as a + b*i) and rank 2 (quaternions) can be converted


SymPy Quaternions
~~~~~~~~~~~~~~~~~

SymPy’s ``sympy.algebras.quaternion.Quaternion`` is converted with
``some_hy.to_sympy()`` and ``Hy.from_sympy(sq)``. SymPy, like ``Hy``,
represents rational numbers exactly, so this pair round-trips with no
rounding at all.

.. code:: ipython3

    >>> q.to_sympy()   # exact -- no rounding at all




.. math::

    \displaystyle 1 + 2 i + 3 j + 4 k



.. code:: ipython3

    >>> Hy.from_sympy(sympy.Quaternion(1, sympy.Rational(2, 3), -4, 0))




.. parsed-literal::

    Hy(Hy('1', '2/3'), Hy('-4', '0'))



SymPy, unlike the other two packages, represents rational numbers
exactly, so ``to_sympy()``/``from_sympy()`` need no rounding decision at
all: a ``Hy -> sympy -> Hy`` round trip always returns the original
value.

.. code:: ipython3

    >>> h = Hy.from_array(['1/3', '-2/7', '5/9', '1/10'])
    
    >>> Hy.from_sympy(h.to_sympy()) == h   # exact, no keywords needed




.. parsed-literal::

    True



numpy-quaternion
~~~~~~~~~~~~~~~~

The package **numpy-quaternion**, imported as ``quaternion``, is
converted with ``some_hy.to_numpy_quaternion()`` and
``Hy.from_numpy_quaternion(nq)``. It stores
``numpy.quaternion(w, x, y, z)`` as ``float64``.

.. code:: ipython3

    >>> q.to_numpy_quaternion()




.. parsed-literal::

    quaternion(1, 2, 3, 4)



.. code:: ipython3

    >>> Hy.from_numpy_quaternion(np.quaternion(1.5, -0.25, 3, 0.1))




.. parsed-literal::

    Hy(Hy('3/2', '-1/4'), Hy('3', '1/10'))



quaternionic
~~~~~~~~~~~~

The package **quaternionic** is converted with
``some_hy.to_quaternionic()`` and ``Hy.from_quaternionic(qa)``. It
stores an array ``[w, x, y, z]`` as ``float64``.

.. code:: ipython3

    >>> q.to_quaternionic()




.. parsed-literal::

    quaternionic.array([1., 2., 3., 4.])



.. code:: ipython3

    >>> Hy.from_quaternionic(quaternionic.array([1, 0.5, 0, -2]))




.. parsed-literal::

    Hy(Hy('1', '1/2'), Hy('0', '-2'))



Floats versus Exact Rationals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``numpy-quaternion`` and ``quaternionic`` both store ``float64``
coordinates, whereas a ``Hy`` is exact. Converting *to* either package
rounds, and converting *back* requires a choice of which rational each
float stands for. Their ``from_*`` methods offer three choices, selected
with two optional keyword arguments:

- By default, each float becomes the rational with the same shortest
  decimal representation (``0.1`` becomes ``1/10``), just as
  ``Hy.from_array()`` does.
- ``max_denominator=N`` chooses the closest fraction whose denominator
  is at most ``N``. This recovers simple fractions such as ``1/3``.
- ``exact=True`` uses the exact binary value of each float. (This cannot
  be combined with ``max_denominator``.)

.. code:: ipython3

    >>> h = Hy.from_array(['1/3', '-2/7', '5/9', '1/10'])
    >>> nq = h.to_numpy_quaternion()
    
    >>> Hy.from_numpy_quaternion(nq) == h                        # default: decimal digits, not 1/3




.. parsed-literal::

    False



.. code:: ipython3

    >>> Hy.from_numpy_quaternion(nq, max_denominator=1000) == h  # recovers 1/3, -2/7, ...




.. parsed-literal::

    True



.. code:: ipython3

    >>> Hy.from_numpy_quaternion(nq, exact=True).to_array()[0]   # the exact value of the float nearest 1/3




.. parsed-literal::

    Fraction(6004799503160661, 18014398509481984)



Whichever option is used, a float survives the round trip *float to
``Hy`` to float* unchanged.

``nan`` and infinite coordinates cannot be represented as rationals, and
raise a ``ValueError``.

.. code:: ipython3

    >>> Hy.from_numpy_quaternion(nq, exact=True).to_numpy_quaternion() == nq      # floats survive the trip




.. parsed-literal::

    True



``sympy``, by contrast, stores rational numbers exactly too, so
``to_sympy()``/``from_sympy()`` need none of this: a
``Hy -> sympy -> Hy`` round trip is always exact (see `SymPy
Quaternions <#sympy-quaternions>`__).

``from_sympy`` still accepts the same ``exact``/``max_denominator``
keywords as the other two ``from_*`` methods, for the less common case
of a ``Quaternion`` that was itself built from ``sympy.Float``
components (rather than ints or ``Rational``\ s), which SymPy does *not*
store exactly.

.. code:: ipython3

    >>> approx = sympy.Quaternion(sympy.Float(1 / 3), 0, 0, 0)   # a Quaternion built from a float
    
    >>> print(Hy.from_sympy(approx))                                  # default: decimal digits
    >>> print(Hy.from_sympy(approx, max_denominator=1000))             # recovers 1/3


.. parsed-literal::

    (3333333333333333/10000000000000000)
    (1/3)


Multiplication agrees across ``Hy`` and all three packages – a useful
sanity check on the conventions:

.. code:: ipython3

    >>> a, b = Hy.random(2, seed=1), Hy.random(2, seed=2)
    
    >>> print("Hy               :", a * b)
    >>> print("sympy (exact)    :", a.to_sympy() * b.to_sympy())
    >>> print("numpy-quaternion :", a.to_numpy_quaternion() * b.to_numpy_quaternion())
    >>> print("quaternionic     :", a.to_quaternionic() * b.to_quaternionic())


.. parsed-literal::

    Hy               : (14/9+131/6i+39/4j-215/18k)
    sympy (exact)    : 14/9 + 131/6*i + 39/4*j + (-215/18)*k
    numpy-quaternion : quaternion(1.55555555555555, 21.8333333333333, 9.75, -11.9444444444444)
    quaternionic     : [  1.55555556  21.83333333   9.75       -11.94444444]


gint (Gaussian Integers and Rationals)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The separate package ``gint``
(https://github.com/alreich/gaussian-integers) implements the Gaussian
integers, ``Zi``, and the Gaussian rationals, ``Qi``. A ``Qi`` is
exactly a rank 1 ``Hy`` whose ``mu`` is :math:`-1`, and a ``Zi`` is the
same thing restricted to integer coordinates, so the two can be
converted into one another, exactly, in either direction. As with the
quaternion packages above, ``gint`` is not required by ``hyprat``; it is
imported only when a method that needs it is called. Install it from
GitHub:

::

   pip install git+https://github.com/alreich/gaussian-integers.git

(Do not ``pip install gint``: that name belongs to an unrelated package
on PyPI. ``hyprat`` notices this and says so.)

The conversions are explicit. A ``Zi`` or ``Qi`` is never ``==`` to a
``Hy``, so mixing them in arithmetic fails loudly instead of silently
picking an algebra.

- ``Hy.from_gint(x)`` builds a rank 1 ``Hy`` from a ``Zi`` or ``Qi``. It
  needs no ``gint`` import at all.
- ``some_hy.to_qi()`` returns a ``Qi``, which, as ``Qi`` itself does,
  collapses to a ``Zi`` when both coordinates are integers.
- ``some_hy.to_zi()`` returns a ``Zi``, or raises ``ValueError`` if a
  coordinate is not an integer.
- ``some_hy.is_gaussian()`` is a cheap predicate that needs no ``gint``.

.. code:: ipython3

    >>> Hy.from_gint(Qi('1/2', '-3/5'))




.. parsed-literal::

    Hy('1/2', '-3/5')



.. code:: ipython3

    >>> Hy.from_gint(Zi(2, -3))




.. parsed-literal::

    Hy('2', '-3')



.. code:: ipython3

    >>> Hy('1/2', '-3/5').to_qi()




.. parsed-literal::

    Qi('1/2', '-3/5')



.. code:: ipython3

    >>> Hy(2, 3).to_qi()   # integer coordinates: Qi collapses to a Zi




.. parsed-literal::

    Zi(2, 3)



.. code:: ipython3

    >>> Hy(2, -3).to_zi()




.. parsed-literal::

    Zi(2, -3)



.. code:: ipython3

    >>> try:
    ...     Hy('1/2', 3).to_zi()
    ... except ValueError as e:
    ...     print(e)


.. parsed-literal::

    cannot convert to a gint Zi: the coordinates are not all integers (use to_qi() for a Gaussian rational)


Which Values Are Gaussian?
^^^^^^^^^^^^^^^^^^^^^^^^^^

``is_gaussian()`` asks whether a value lies in the Gaussian rationals,
:math:`\mathbb{Q}[i]`: its lowest doubling level must have
:math:`\mu = -1` and every coordinate beyond the first two must be zero.
It is a statement about the *value*, not just its rank, so a quaternion
such as :math:`1 + 2i` qualifies (it equals the complex number
:math:`1 + 2j` in ``Hy``\ ’s rank 1 notation), while a split-complex
number or a quaternion with a nonzero :math:`j` or :math:`k` part does
not. ``to_qi()`` succeeds exactly when ``is_gaussian()`` is true.

.. code:: ipython3

    >>> print(Hy('1/2', '-3/5').is_gaussian())              # a rank 1 value
    >>> print(Hy(Hy(1, 2), 0).is_gaussian())                 # the quaternion 1+2i
    >>> print(Hy(Hy(1, 2), Hy(0, 1)).is_gaussian())          # 1+2i+0j+1k: a k part
    >>> print(Hy(1, 2, mu=1).is_gaussian())                  # split-complex


.. parsed-literal::

    True
    True
    False
    False


.. code:: ipython3

    >>> Hy(Hy(1, 2), 0).to_zi()    # a quaternion that happens to be Gaussian




.. parsed-literal::

    Zi(1, 2)



Strings
^^^^^^^

``Hy.parse`` reads everything ``Zi`` and ``Qi`` print. With the default
unit symbol, ``j``, nothing special is needed. A ``Zi`` or ``Qi`` whose
unit symbol has been switched to ``i`` prints ``(2-3i)``, which by
itself means a *quaternion* to ``Hy.parse``; pass ``unit='i'`` to say
that the text uses that one imaginary unit, so it is read as a rank 1
value.

.. code:: ipython3

    >>> Hy.parse(str(Qi('1/2', '-3/5')))




.. parsed-literal::

    Hy('1/2', '-3/5')



.. code:: ipython3

    >>> Zi.set_unit_symbol('i'); Qi.set_unit_symbol('i')
    
    >>> text = str(Zi(2, -3))
    >>> print(text)
    >>> print(repr(Hy.parse(text)))                 # a quaternion, by default
    >>> print(repr(Hy.parse(text, unit='i')))       # the rank 1 value that was meant
    
    >>> Zi.set_unit_symbol('j'); Qi.set_unit_symbol('j')


.. parsed-literal::

    (2-3i)
    Hy(Hy('2', '-3'), Hy('0', '0'))
    Hy('2', '-3')


Cross-Checking the Two Implementations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because ``Qi`` and a rank 1 ``Hy`` are two independent implementations
of the same arithmetic, either one is a check on the other. Convert,
operate, and compare:

.. code:: ipython3

    >>> a, b = Qi('1/2', '-3/5'), Qi(2, '1/3')
    >>> ha, hb = Hy.from_gint(a), Hy.from_gint(b)
    
    >>> print(Hy.from_gint(a * b) == ha * hb)
    >>> print(Hy.from_gint(a / b) == ha / hb)
    >>> print(Fraction(a.norm) == ha.norm_squared())


.. parsed-literal::

    True
    True
    True


The next page, `Hurwitz Integers and
Examples <hurwitz_and_examples.html>`__, introduces the integer-valued
class ``Hu``, which plays the same role for the quaternions that ``Zi``
plays for the complex numbers.
