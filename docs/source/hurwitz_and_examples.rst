Hurwitz Integers and Examples
=============================

This page continues `Split-Hypercomplex Numbers and
Interoperability <split_hy_and_interop.html>`__. It introduces the class
``Hu``, which holds the **Hurwitz integers**, the quaternion analogue of
the Gaussian integers, and with them the number theory of quaternions:
division with remainder, greatest common divisors, primes and
factorization. The pure-Python integer factoring module
``hyprat.intfactor``, on which the factorization depends, has a section
of its own. The page ends with six worked **Examples** and a **Summary**
of the whole package.

``Hu`` builds on ``Hy``, so the basic ideas are in ```hyprat``: Exact
Rational Hypercomplex Numbers <classic_hy.html>`__. Nothing on this page
needs the optional packages, except that the last check in the example
on the binary tetrahedral group uses ``finalg`` if it happens to be
installed.

.. code:: ipython3

    >>> import contextlib
    >>> import io
    >>> import random
    >>> from collections import Counter
    >>> from fractions import Fraction
    >>> from hyprat import Hu, Hy
    >>> from hyprat.intfactor import factorint, is_probable_prime, primes_upto, sqrt_mod_prime

Hurwitz Integers: the Class ``Hu``
----------------------------------

The *Hurwitz integers* are the quaternions
:math:`a + b\,i + c\,j + d\,k` whose four coordinates are either **all
integers** or **all halves of odd integers**, so that
:math:`(1 + i + j + k)/2` is one of them. They form a ring in which
every element has an integer norm, and which has 24 units. That makes
them the quaternion analogue of the Gaussian integers, with ``Hy``
playing the role that ``Qi`` plays for ``Zi``: ``Hy`` is the ambient
algebra over the rationals, and the class ``Hu`` holds the integer
points inside it.

A ``Hu`` stores its four coordinates **doubled**, as plain Python
integers of the same parity, so its arithmetic involves no fractions at
all. A set of coordinates that is not a Hurwitz integer cannot be
constructed.

What Are the Hurwitz Integers?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adolf Hurwitz noticed in 1896 that the quaternions with integer
coordinates (the *Lipschitz integers*) are not the best quaternion
analogue of the ordinary integers. They lack division with remainder.
Adding the half-integer points :math:`(\pm 1 \pm i \pm j \pm k)/2` gives
a ring, the *Hurwitz integers*, that has it, and in which Euclid’s
algorithm, greatest common divisors and unique factorization (in a
suitable sense) all work. The `Examples section on
division <#division-with-remainder-and-the-euclidean-algorithm>`__ shows
what goes wrong without the half-integers.

The rest of this part follows the order in which one meets the ideas:
construction and arithmetic, the units, conversion to and from ``Hy``,
divisibility and greatest common divisors, then (after a section on the
integer factoring that everything else depends on) primes,
factorization, counting, and the sum of four squares.

.. code:: ipython3

    >>> q = Hu(1, 2, 3, 4)
    >>> print(q)
    >>> q


.. parsed-literal::

    (1+2i+3j+4k)




.. parsed-literal::

    Hu(1, 2, 3, 4)



.. code:: ipython3

    >>> w = Hu('1/2', '1/2', '1/2', '1/2')   # (1+i+j+k)/2: not a Lipschitz integer
    >>> print(w)
    >>> print(w.is_lipschitz(), w.doubled)


.. parsed-literal::

    (1/2+1/2i+1/2j+1/2k)
    False (1, 1, 1, 1)


.. code:: ipython3

    >>> try:
    ...     Hu('1/2', 0, 0, 0)
    ... except ValueError as e:
    ...     print(e)


.. parsed-literal::

    coordinates (1/2, 0, 0, 0) must be all integers or all halves of odd integers


Arithmetic
~~~~~~~~~~

``+``, ``-``, ``*`` and ``**`` work between ``Hu`` values and with plain
integers. Multiplication is not commutative. There is deliberately no
``/``: the Hurwitz integers are not closed under division, so use
``to_hy()`` for exact rational division.

.. code:: ipython3

    >>> i, j, k = Hu(0, 1, 0, 0), Hu(0, 0, 1, 0), Hu(0, 0, 0, 1)
    
    >>> print(i * j, j * i)
    >>> print(q * w)
    >>> print(w * w, w ** 3, w ** 6)


.. parsed-literal::

    (k) (-k)
    (-4+i+3j+2k)
    (-1/2+1/2i+1/2j+1/2k) (-1) (1)


Norm, Trace and Conjugate
~~~~~~~~~~~~~~~~~~~~~~~~~

The norm :math:`N(q) = q\,\bar q = a^2 + b^2 + c^2 + d^2` is always an
``int`` (the same quantity that ``Hy.norm()`` gives), and it is
multiplicative. The trace :math:`q + \bar q = 2a` is an ``int`` too.

.. code:: ipython3

    >>> print(q.norm(), q.trace(), q.conjugate())
    >>> print(w.norm(), w.trace())
    >>> print(q * q.conjugate())
    >>> print((q * w).norm() == q.norm() * w.norm())


.. parsed-literal::

    30 2 (1-2i-3j-4k)
    1 1
    (30)
    True


The 24 Units
~~~~~~~~~~~~

The units are the elements of norm 1: the eight Lipschitz units
:math:`\pm 1, \pm i, \pm j, \pm k` and the sixteen
:math:`(\pm 1 \pm i \pm j \pm k)/2`. Only the first eight are
``Hy.units(2)``. A unit has an inverse in the ring, its conjugate, so
negative powers are allowed for units and only for units.

.. code:: ipython3

    >>> units = Hu.units()
    >>> print(len(units))
    >>> list(units)[:10]


.. parsed-literal::

    24




.. parsed-literal::

    ['1',
     '-1',
     'i',
     '-i',
     'j',
     '-j',
     'k',
     '-k',
     '1/2+1/2i+1/2j+1/2k',
     '1/2+1/2i+1/2j-1/2k']



.. code:: ipython3

    >>> print(w ** -1, w ** -1 == w.conjugate())
    >>> print(w.is_unit(), q.is_unit())
    
    >>> try:
    ...     q ** -1
    ... except ValueError as e:
    ...     print(e)


.. parsed-literal::

    (1/2-1/2i-1/2j-1/2k) True
    True False
    only units have inverses in the Hurwitz integers, so a negative power needs a unit; use to_hy() ** n for rational values


Converting To and From ``Hy``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The conversions are explicit, and a ``Hu`` is never ``==`` to a ``Hy``.
``Hy`` has matching predicates, ``is_hurwitz()`` and ``is_lipschitz()``,
which are about the *value*: a rank 1 value :math:`a + b\,j` counts as
the quaternion :math:`a + b\,i`, and so does a value of higher rank
whose coordinates beyond the first four are all zero. ``Hy.embed(rank)``
pads with zeros to a higher rank.

.. code:: ipython3

    >>> h = q.to_hy()
    >>> print(repr(h))
    >>> print(h == q, Hu.from_hy(h) == q)


.. parsed-literal::

    Hy(Hy('1', '2'), Hy('3', '4'))
    False True


.. code:: ipython3

    >>> print(Hy(Hy('1/2', '1/2'), Hy('1/2', '-3/2')).is_hurwitz())
    >>> print(Hy(Hy('1/2', '1/2'), Hy(0, 0)).is_hurwitz())
    >>> print(Hy(2, -3).is_hurwitz(), Hy(2, -3).is_lipschitz())    # the Gaussian integer 2-3i
    >>> Hu.from_hy(Hy(2, -3))


.. parsed-literal::

    True
    False
    True True




.. parsed-literal::

    Hu(2, -3, 0, 0)



.. code:: ipython3

    >>> Hu.from_hy(Hy(Hy(1, 2), Hy(0, 0)).embed(3))               # an octonion inside the quaternions




.. parsed-literal::

    Hu(1, 2, 0, 0)



.. code:: ipython3

    >>> try:
    ...     Hu.from_hy(Hy(1, 2, mu=1))        # split-complex: not in the classical quaternions
    ... except ValueError as e:
    ...     print(e)


.. parsed-literal::

    not a Hurwitz integer: its two lowest doubling levels do not both have mu = -1, so it is not in the classical quaternions


Text Forms and Parsing
~~~~~~~~~~~~~~~~~~~~~~

``str(Hu)`` and ``Hu('...')`` always use the quaternion labels ``i``,
``j``, ``k``. That differs from ``Hy.parse``, which reads text that
mentions only ``j`` as a *rank 1* value (the way ``gint``\ ’s ``Zi`` and
``Qi`` print). To read ``gint``-style text as a Hurwitz integer, go
through ``Hy`` explicitly.

.. code:: ipython3

    >>> print(repr(Hu('(3j)')))                        # the quaternion 3j
    >>> print(repr(Hu.from_hy(Hy.parse('(3j)'))))     # the rank 1 value 3j, read as 3i
    >>> print(repr(Hu('(1/2-1/2i+1/2j+1/2k)')))


.. parsed-literal::

    Hu(0, 0, 3, 0)
    Hu(0, 3, 0, 0)
    Hu('1/2', '-1/2', '1/2', '1/2')


Cross-Checking against ``Hy``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because ``Hy`` and ``Hu`` are independent implementations of the same
arithmetic, either is a check on the other:

.. code:: ipython3

    >>> a, b = Hu(1, -2, 3, 4), Hu('3/2', '1/2', '-5/2', '1/2')
    >>> ha, hb = a.to_hy(), b.to_hy()
    
    >>> print((a * b).to_hy() == ha * hb)
    >>> print((a + b).to_hy() == ha + hb)
    >>> print((a * b).norm() == (ha * hb).norm())


.. parsed-literal::

    True
    True
    True


Division with Remainder
~~~~~~~~~~~~~~~~~~~~~~~

The Hurwitz integers are a Euclidean domain on *both* sides, but the
ring is not commutative, so each notion comes in a left and a right
form. In this class, ``g`` is a **left divisor** of ``a`` when
``a = g*x``, and a **right divisor** when ``a = x*g``. For ``b != 0``,
``divmod_left`` finds ``q``, ``r`` with ``a = b*q + r``, and
``divmod_right`` finds ``a = q*b + r``. In both cases
``2*N(r) <= N(b)``, which is what makes Euclid’s algorithm terminate.
There is no ``/`` operator, because the ring is not closed under
division.

.. code:: ipython3

    >>> a, b = Hu(1, 2, 3, 5), Hu(1, 1, 1, 0)
    
    >>> q, r = a.divmod_left(b)
    >>> print(repr(q), repr(r))
    >>> print(a == b * q + r, 2 * r.norm() <= b.norm())
    
    >>> q, r = a.divmod_right(b)
    >>> print(repr(q), repr(r))
    >>> print(a == q * b + r)


.. parsed-literal::

    Hu(2, -1, 2, 1) Hu(0, 0, 0, 1)
    True True
    Hu(2, 2, -1, 2) Hu(0, 0, 0, 0)
    True


Divisibility
~~~~~~~~~~~~

Because multiplication does not commute, ``g*x`` and ``x*g`` are
different numbers, and a product is generally divisible on one side
only.

.. code:: ipython3

    >>> g, x = Hu(1, 1, 1, 0), Hu(0, 0, 1, 1)
    >>> a = g * x
    >>> print(repr(a))
    >>> print(g.left_divides(a), g.right_divides(a))
    >>> print(repr(a.div_exact_left(g)))
    >>> try:
    ...     a.div_exact_right(g)
    ... except ValueError as e:
    ...     print(e)


.. parsed-literal::

    Hu(-1, 1, 0, 2)
    True False
    Hu(0, 0, 1, 1)
    (1+i+j) is not a right divisor of (-1+i+2k)


Associates and the Canonical Associate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two numbers are **left associates** if ``b = u*a`` for a unit ``u``, and
**right associates** if ``b = a*u``. A nonzero number has 24 associates
on each side, and ``canonical_associate(side)`` picks one of them, the
one whose doubled coordinates are largest, so that every unit becomes
``1``.

.. code:: ipython3

    >>> a = Hu(1, 2, 3, 0)
    >>> i = Hu(0, 1, 0, 0)
    >>> print(a.is_left_associate(i * a), a.is_right_associate(i * a))
    >>> print(repr(a.canonical_associate('left')))
    >>> print(repr((i * a).canonical_associate('left')))
    >>> print(repr(Hu('1/2', '1/2', '-1/2', '1/2').canonical_associate('right')))


.. parsed-literal::

    True False
    Hu(3, 2, 0, -1)
    Hu(3, 2, 0, -1)
    Hu(1, 0, 0, 0)


Greatest Common Divisors
~~~~~~~~~~~~~~~~~~~~~~~~

``gcld`` is the greatest common *left* divisor, determined only up to a
unit on its right, so the canonical right associate is returned.
``gcrd`` is the mirror image. The extended forms also return Bezout
coefficients: ``g = a*x + b*y`` for ``xgcld`` and ``g = x*a + y*b`` for
``xgcrd``.

.. code:: ipython3

    >>> a, b = Hu(-2, 1, 3, 3), Hu(3, -3, -1, -3)
    >>> g, x, y = a.xgcld(b)
    >>> print(repr(g), repr(x), repr(y))
    >>> print(g == a * x + b * y)
    
    >>> c = Hu(1, 1, 1, 0)
    >>> print(repr((c * a).gcld(c * b)))
    >>> print(repr(c.canonical_associate('right')))


.. parsed-literal::

    Hu(1, 0, 0, 0) Hu(0, -1, -1, -1) Hu(0, -1, 0, -1)
    True
    Hu('3/2', '1/2', '-1/2', '1/2')
    Hu('3/2', '1/2', '-1/2', '1/2')


Integer Primality Testing and Factoring in Pure Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Everything that follows depends on one fact: a Hurwitz integer :math:`q`
is prime exactly when its norm :math:`N(q)` is an ordinary prime, and
the primes that make up a factorization of :math:`q` have the prime
factors of :math:`N(q)` as their norms. Factoring a Hurwitz integer
therefore starts with factoring an ordinary integer. Writing a prime
:math:`p` as a sum of four squares needs a square root modulo :math:`p`.
The module ``hyprat.intfactor`` supplies exactly these tools, in pure
Python and with no required dependencies. It exports four functions:

+-----------------------------------+-----------------------------------+
| Function                          | What it does                      |
+===================================+===================================+
| ``primes_upto(limit)``            | The list of primes up to          |
|                                   | ``limit``, by the sieve of        |
|                                   | Eratosthenes                      |
+-----------------------------------+-----------------------------------+
| ``is_probable_prime(n)``          | Whether ``n`` is prime. The       |
|                                   | answer is certain for             |
|                                   | :math:`n < 3.3 \times 10^{24}`;   |
|                                   | beyond that it is a Baillie-PSW   |
|                                   | test                              |
+-----------------------------------+-----------------------------------+
| ``factorint(n, method="python")`` | The prime factorization of ``n``  |
|                                   | as a dict ``{prime: exponent}``,  |
|                                   | in increasing order of the primes |
+-----------------------------------+-----------------------------------+
| ``sqrt_mod_prime(a, p)``          | A square root of ``a`` modulo an  |
|                                   | odd prime ``p`` (Tonelli-Shanks)  |
+-----------------------------------+-----------------------------------+

The module is small on purpose. It is not meant to compete with a number
theory library. It is meant to make ``Hu`` self-contained and fast
enough for the norms one meets in practice.

Primes and Primality Testing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``primes_upto`` sieves, so it is meant for modest limits.
``is_probable_prime`` first removes small prime factors, then, below
:math:`3.3 \times 10^{24}`, runs strong probable-prime tests to a fixed
set of bases that is known to be correct for every number in that range,
so the answer there is a proof and not merely a probability. Above that
range it uses the **Baillie-PSW** test, a strong Fermat test to base 2
followed by a strong Lucas test. No composite number is known that
passes Baillie-PSW, although it has not been proved that none exists.
The word “probable” in the name refers to that last, very large range.

.. code:: ipython3

    >>> print(primes_upto(50))
    >>> print([k for k in range(30) if is_probable_prime(k)])
    >>> print(is_probable_prime(3215031751))        # composite, but a strong pseudoprime to the bases 2, 3, 5 and 7
    >>> print(is_probable_prime(2**61 - 1), is_probable_prime(2**89 - 1), is_probable_prime(2**89 + 1))


.. parsed-literal::

    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    False
    True True False


Factoring Integers
^^^^^^^^^^^^^^^^^^

``factorint`` divides out the primes below 1000 first. What remains is
either a prime, which ``is_probable_prime`` recognizes, or it is split
by **Pollard’s rho method** in Brent’s variant, and the pieces are
treated the same way. The keyword ``method`` chooses the engine:

- ``"python"`` (the default) is pure Python and never imports anything.
- ``"sympy"`` uses ``sympy.factorint`` and raises ``ImportError`` if
  SymPy is not installed.
- ``"auto"`` uses SymPy for inputs of more than 100 bits if it is
  installed, and pure Python otherwise.

All three return the same answer. They differ only in speed on large
inputs, and only the first is guaranteed to be available. Pollard’s rho
finds a factor in about :math:`\sqrt{p}` steps, where :math:`p` is the
*smallest* prime factor, so a product of two primes of 30 digits each is
out of reach, while a number with a small prime factor is easy however
large it is.

.. code:: ipython3

    >>> print(factorint(1), factorint(360))
    >>> print(factorint(600851475143))
    >>> print(factorint(2**64 + 1))
    >>> print(factorint(10**12 + 39))


.. parsed-literal::

    {} {2: 3, 3: 2, 5: 1}
    {71: 1, 839: 1, 1471: 1, 6857: 1}
    {274177: 1, 67280421310721: 1}
    {1000000000039: 1}


.. code:: ipython3

    >>> try:
    ...     print(factorint(2**64 + 1, method="sympy") == factorint(2**64 + 1))
    ... except ImportError as e:
    ...     print(e)


.. parsed-literal::

    True


The functions check their arguments. ``factorint`` accepts only positive
integers, and none of the functions accepts a float or a ``bool`` where
an ``int`` is expected.

.. code:: ipython3

    >>> for bad in (0, -6, 6.0):
    ...     try:
    ...         factorint(bad)
    ...     except (TypeError, ValueError) as e:
    ...         print(type(e).__name__, "-", e)


.. parsed-literal::

    ValueError - n must be positive, not 0
    ValueError - n must be positive, not -6
    TypeError - n must be an int, not float


Square Roots Modulo a Prime
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Writing an odd prime :math:`p` as a sum of four squares starts with
:math:`x, y` such that :math:`x^2 + y^2 + 1 \equiv 0 \pmod p`, as
described in `Sums of Four Squares <#sums-of-four-squares>`__ below.
Finding :math:`x` needs a square root of :math:`-1 - y^2` modulo
:math:`p`. ``sqrt_mod_prime(a, p)`` returns *a* root, one of the two
:math:`\pm r`, and raises ``ValueError`` if ``a`` is not a square modulo
``p``. The caller is responsible for ``p`` being an odd prime.

.. code:: ipython3

    >>> root = sqrt_mod_prime(10, 13)
    >>> print(root, root * root % 13)
    >>> p = 10**9 + 9                                     # p % 4 == 1, so -1 is a square mod p
    >>> r = sqrt_mod_prime(-1, p)
    >>> print(r * r % p == p - 1)
    >>> try:
    ...     sqrt_mod_prime(2, 5)
    ... except ValueError as e:
    ...     print(e)


.. parsed-literal::

    7 10
    True
    2 is not a square modulo 5


Big Numbers and the Integer Factoring Helper
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As explained above, factoring ``q`` means factoring its norm, an
ordinary integer, and ``hyprat.intfactor`` does that. Norms of the size
one gets from coordinates with 9 digits, as below, are routine.

.. code:: ipython3

    >>> print(factorint(2**64 + 1))
    >>> print(is_probable_prime(2**89 - 1))
    
    >>> q = Hu(123456789, 987654321, 192837465, 918273645)
    >>> print(q.norm(), factorint(q.norm()))
    >>> f = q.factor()
    >>> print(f.content, f.norms(), f.product() == q)


.. parsed-literal::

    {274177: 1, 67280421310721: 1}
    True
    1871115411549373812 {2: 2, 3: 5, 2855261: 1, 674199611: 1}
    18 (3, 2855261, 674199611) True


Primes, Content and Primitive Parts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A nonzero Hurwitz integer is **prime** exactly when its norm is an
ordinary prime number. Notice that ``2`` is *not* prime here: it has
norm 4, and factors as a product of two primes of norm 2. The
**content** of an element is the largest positive integer that divides
it, and it is **primitive** when its content is 1.

.. code:: ipython3

    >>> print(Hu(1, 1, 0, 0).is_prime(), Hu(2).is_prime(), Hu(1, 2, 3, 4).is_prime())
    
    >>> a = Hu(2, 4, 6, 8)
    >>> print(a.content(), repr(a.primitive_part()), a.is_primitive())
    >>> print(Hu(1, 1, 1, 1).content())          # 2 times the unit (1+i+j+k)/2


.. parsed-literal::

    True False False
    2 Hu(1, 2, 3, 4) False
    2


Factoring into Primes
~~~~~~~~~~~~~~~~~~~~~

By the theorem of Conway and Smith, a primitive Hurwitz integer ``q``
can be written as a product of primes whose norms are the prime factors
of ``N(q)`` **in any order you like**. The factors are unique up to
*unit migration*, in which a unit moves from one side of a factor to the
other. ``factor`` returns the pieces as a named tuple with a
``product()`` method that rebuilds the value.

.. code:: ipython3

    >>> q = Hu(1, 2, 3, 4)                      # norm 30 = 2 * 3 * 5
    >>> f = q.factor()
    >>> print(f.norms())
    >>> for p in f.factors:
    ...     print(repr(p), p.is_prime())
    >>> print(f.product() == q)
    
    >>> g = q.factor(order=[5, 3, 2])
    >>> print(g.norms())
    >>> for p in g.factors:
    ...     print(repr(p))
    >>> print(g.product() == q)


.. parsed-literal::

    (2, 3, 5)
    Hu(1, 1, 0, 0) True
    Hu('3/2', '1/2', '1/2', '1/2') True
    Hu('3/2', '1/2', '3/2', '-1/2') True
    True
    (5, 3, 2)
    Hu(2, 0, 1, 0)
    Hu('3/2', '1/2', '-1/2', '-1/2')
    Hu(0, 0, 1, 1)
    True


``reorder`` gives the factorization of the same value with the primes in
a new order. Each factor is found as a greatest common left divisor,
``pi = gcld(q, p)``, and then divided out of ``q`` on the left.

.. code:: ipython3

    >>> print(f.reorder([3, 5, 2]).norms())
    >>> print(f.reorder([3, 5, 2]).product() == q)


.. parsed-literal::

    (3, 5, 2)
    True


When the Element Is Not Primitive
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A rational prime ``p`` is ``pi * conj(pi)`` for many different primes
``pi``, so a non-primitive element has no single factorization. Instead,
``factor`` reports the integer content separately and factors the
primitive part. If nothing is left but a unit, the unit is reported too:
the value is always ``content * unit * f1 * ... * fk``.

.. code:: ipython3

    >>> f = Hu(6, 12, 18, 24).factor()
    >>> print(f.content, f.norms())
    >>> print(f.product() == Hu(6, 12, 18, 24))
    
    >>> f = Hu(14, 0, 0, 0).factor()
    >>> print(f.content, repr(f.unit), f.factors)


.. parsed-literal::

    6 (2, 3, 5)
    True
    14 Hu(1, 0, 0, 0) ()


Enumerating Hurwitz Integers of a Given Norm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Exactly :math:`24\,\sigma_{\mathrm{odd}}(n)` Hurwitz integers have norm
:math:`n`, where :math:`\sigma_{\mathrm{odd}}(n)` is the sum of the odd
divisors of :math:`n`. ``Hu.of_norm(n)`` lists them (the work grows
roughly like :math:`n`), and ``Hu.count_of_norm(n)`` counts them from
the factorization of :math:`n` without listing anything.

.. code:: ipython3

    >>> print([len(Hu.of_norm(n)) for n in range(1, 9)])
    >>> print([Hu.count_of_norm(n) for n in range(1, 9)])
    >>> print(Hu.count_of_norm(10**18))
    
    >>> for h in Hu.of_norm(2)[:4]:
    ...     print(repr(h))


.. parsed-literal::

    [24, 24, 96, 24, 144, 96, 192, 24]
    [24, 24, 96, 24, 144, 96, 192, 24]
    114440917968744
    Hu(-1, -1, 0, 0)
    Hu(-1, 0, -1, 0)
    Hu(-1, 0, 0, -1)
    Hu(-1, 0, 0, 1)


Restricting to elements with integer coordinates recovers Jacobi’s
four-square count, and ``primitive=True`` keeps only the elements of
content 1:

.. code:: ipython3

    >>> n = 3
    >>> print(len(Hu.of_norm(n)), len([h for h in Hu.of_norm(n) if h.is_lipschitz()]))
    >>> print(len(Hu.of_norm(n, primitive=True)), len(Hu.of_norm(4, primitive=True)))


.. parsed-literal::

    96 32
    96 0


Primes of a Given Norm
~~~~~~~~~~~~~~~~~~~~~~

For an odd prime :math:`p` the :math:`24(p+1)` primes of norm :math:`p`
fall into :math:`p+1` classes of right associates, and for :math:`p=2`
there is a single class. ``primes_of_norm`` returns one canonical
representative of each class, or every prime with ``associates=True``.

.. code:: ipython3

    >>> for p in (2, 3, 5):
    ...     print(p, len(Hu.primes_of_norm(p)), len(Hu.primes_of_norm(p, associates=True)))
    
    >>> for h in Hu.primes_of_norm(3):
    ...     print(repr(h))


.. parsed-literal::

    2 1 24
    3 4 96
    5 6 144
    Hu('3/2', '1/2', '-1/2', '-1/2')
    Hu('3/2', '1/2', '-1/2', '1/2')
    Hu('3/2', '1/2', '1/2', '-1/2')
    Hu('3/2', '1/2', '1/2', '1/2')


Sums of Four Squares
~~~~~~~~~~~~~~~~~~~~

Lagrange’s four-square theorem falls out of the same machinery. For an
odd prime :math:`p` one finds :math:`x, y` with
:math:`x^2 + y^2 + 1 \equiv 0 \pmod p` and takes
:math:`\pi = \gcd_{\text{left}}(x + y\,i + j,\ p)`, a prime of norm
:math:`p`; multiplying a unit onto it gives integer coordinates. Norms
multiply, so ``Hu.with_norm(n)`` gets an element of norm :math:`n` by
multiplying such primes, and ``Hu.four_squares(n)`` reads off the four
squares. The cost is that of factoring :math:`n`.

.. code:: ipython3

    >>> print(Hu.four_squares(310))
    >>> a, b, c, d = Hu.four_squares(2**89 - 1)
    >>> print(a, b, c, d)
    >>> print(a*a + b*b + c*c + d*d == 2**89 - 1)
    
    >>> h = Hu.with_norm(1000003)
    >>> print(repr(h), h.norm())


.. parsed-literal::

    (13, 10, 5, 4)
    16505493659067 16505493659067 7162758386013 4775172257342
    True
    Hu(699, 151, 699, 0) 1000003


Examples
--------

The six examples below apply the machinery above. They are written as
ordinary scripts, so the code is plain Python and the results are
printed with ``print``. Each example can be read on its own. The
sections above give the mechanics, and the examples show what they are
good for.

Factoring a Product of Hurwitz Primes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

References
^^^^^^^^^^

- [Conway & Smith 2003] “On Quaternions and Octonions: Their Geometry,
  Arithmetic, and Symmetry”, Conway, John H. and Smith, Derek A., A K
  Peters, 2003
- [De Feo et al. 2020] “SQISign: Compact Post-Quantum Signatures from
  Quaternions and Isogenies”, De Feo, Luca; Kohel, David; Leroux,
  Antonin; Petit, Christophe; Wesolowski, Benjamin, ASIACRYPT 2020
  (`ePrint 2020/1240 <https://eprint.iacr.org/2020/1240>`__)

Introduction
^^^^^^^^^^^^

Quaternion algebras show up in post-quantum cryptography. The
isogeny-based signature scheme SQIsign [De Feo et al. 2020], for
example, relies on the *Deuring correspondence* between supersingular
elliptic curves and maximal orders in a quaternion algebra. This example
doesn’t implement anything like that. Instead, it uses ``Hu`` to look at
one basic ingredient: a product of Hurwitz primes, taken in a particular
order.

Since quaternion multiplication isn’t commutative, the *order* of the
factors matters. According to the Conway–Smith theorem [Conway & Smith
2003], a primitive Hurwitz integer can be factored into primes following
**any** order of prime norms we like, and the factors are unique up to
**unit migration**. That is, the prime :math:`\pi_i` can be replaced by
:math:`u_{i-1}^{-1}\,\pi_i\,u_i`, where each :math:`u_i` is one of the
24 units and :math:`u_0 = u_k = 1`, without changing the product.

In what follows, Alice chooses a secret sequence of primes and
multiplies them together, in order, to get a public key. Then we’ll see
what ``Hu.factor`` recovers from that product and the sequence of norms.

.. code:: ipython3

    # A seeded generator keeps the output reproducible.
    # (A real system would draw its secrets from the `secrets` module.)
    rng = random.Random(2026)
    
    # Public: the sequence of prime norms, in order.
    public_norms = [2, 3, 5, 7]
    
    # Secret: one Hurwitz prime of each norm, in that order.
    # (Any Hurwitz integer of prime norm is a Hurwitz prime.)
    secret_path = [rng.choice(Hu.primes_of_norm(p, associates=True)) for p in public_norms]
    
    # Public: the product of the secret primes, multiplied in order.
    public_key = secret_path[0]
    for prime in secret_path[1:]:
        public_key = public_key * prime
    
    print("Secret path:")
    for p, prime in zip(public_norms, secret_path):
        print(f"  norm {p}: {prime}")
    print(f"\nPublic key: {public_key}   (norm {public_key.norm()})")
    
    # Given the public key and the sequence of norms, Hu.factor finds a factorization:
    found = public_key.factor(order=public_norms)
    print("\nFactors found by Hu.factor(order=public_norms):")
    for p, prime in zip(public_norms, found.factors):
        print(f"  norm {p}: {prime}")
    
    print("\nTheir product is the public key:", found.product() == public_key)
    print("They equal the secret primes, entry by entry:", list(found.factors) == secret_path)


.. parsed-literal::

    Secret path:
      norm 2: (-1+k)
      norm 3: (-1/2+3/2i-1/2j-1/2k)
      norm 5: (3/2+1/2i-1/2j-3/2k)
      norm 7: (1-2i+j+k)
    
    Public key: (-4-7i+12j-k)   (norm 210)
    
    Factors found by Hu.factor(order=public_norms):
      norm 2: (1+i)
      norm 3: (3/2+1/2i-1/2j+1/2k)
      norm 5: (2-i)
      norm 7: (-2-i+j-k)
    
    Their product is the public key: True
    They equal the secret primes, entry by entry: False


The product is right, but the factors aren’t the ones Alice started
with. The difference is the unit migration described above. The
following function finds the units, :math:`u_0, \dots, u_k`, that relate
the two lists. Since :math:`y_i = u_{i-1}^{-1}\,x_i\,u_i`, each
:math:`u_i` can be computed from the one before it, as
:math:`u_i = x_i^{-1}\,(u_{i-1}\,y_i)`, using ``div_exact_left``.

.. code:: ipython3

    def unit_migration(xs, ys):
        "Units u_0 = 1, u_1, ..., u_k = 1 with ys[i] == conj(u[i]) * xs[i] * u[i+1], or None."
        u = Hu(1)
        units = [u]
        for x, y in zip(xs, ys):
            try:
                u = (u * y).div_exact_left(x)        # solves  x * u_next == u * y
            except ValueError:
                return None
            if not u.is_unit():
                return None
            units.append(u)
        return units if units[-1] == Hu(1) else None
    
    units = unit_migration(secret_path, found.factors)
    print("The found factors are the secret primes, migrated by these units:")
    for i, u in enumerate(units):
        print(f"  u_{i} = {u}")
    
    # Check the relationship directly (for a unit, the conjugate is the inverse):
    print(all(units[i].conjugate() * x * units[i + 1] == y
              for i, (x, y) in enumerate(zip(secret_path, found.factors))))


.. parsed-literal::

    The found factors are the secret primes, migrated by these units:
      u_0 = (1)
      u_1 = (-1/2-1/2i-1/2j-1/2k)
      u_2 = (-1/2+1/2i-1/2j+1/2k)
      u_3 = (i)
      u_4 = (1)
    True


Observations
^^^^^^^^^^^^

- The two lists of primes differ only by units, and either list
  multiplies out to the same public key.
- Non-commutativity does play a role: it determines *which* factors we
  get for a given order of norms. Here, the order is part of the public
  information, so there’s nothing left to choose.
- That makes this a toy puzzle, rather than a cryptosystem. In SQIsign
  and related schemes, the difficulty isn’t in factoring an element of a
  fixed quaternion order. It comes from the Deuring correspondence, at
  very large primes :math:`p`: given two supersingular elliptic curves,
  finding an isogeny between them, or computing an endomorphism ring, is
  believed to be hard, even for quantum computers.
- The Hurwitz integers are the maximal order of the quaternion algebra
  ramified at 2 and :math:`\infty`. The corresponding supersingular
  curve is unique (in characteristic 2), so there’s no graph of curves
  in which to take a secret walk.

The next computation gives another view of why the factorization is
determined. For an odd prime :math:`p`, the primitive Hurwitz integers
of norm :math:`p^k`, counted modulo a unit on the right, number
:math:`(p+1)\,p^{k-1}`. That’s also the number of non-backtracking walks
of length :math:`k` in the :math:`(p+1)`-regular tree. A factorization
corresponds to such a walk, and the product determines which one.

.. code:: ipython3

    print(" p  k   classes   (p+1)*p**(k-1)")
    for p, k in [(3, 1), (3, 2), (3, 3), (3, 4), (5, 1), (5, 2), (5, 3), (7, 2)]:
        classes = len(Hu.of_norm(p**k, primitive=True)) // 24      # 24 units per class
        print(f"{p:2d} {k:2d} {classes:9d} {(p + 1) * p**(k - 1):10d}")


.. parsed-literal::

     p  k   classes   (p+1)*p**(k-1)
     3  1         4          4
     3  2        12         12
     3  3        36         36
     3  4       108        108
     5  1         6          6
     5  2        30         30
     5  3       150        150
     7  2        56         56


A Diffie-Hellman-Style Key Exchange
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

References
^^^^^^^^^^

- [Castryck & Decru 2022] “An Efficient Key Recovery Attack on SIDH”,
  Castryck, Wouter and Decru, Thomas, EUROCRYPT 2023 (`ePrint
  2022/975 <https://eprint.iacr.org/2022/975>`__)
- [Maino et al. 2023] “A Direct Key Recovery Attack on SIDH”, Maino,
  Luciano; Martindale, Chloe; Panny, Lorenz; Pope, Giacomo; Wesolowski,
  Benjamin, EUROCRYPT 2023 (`ePrint
  2023/640 <https://eprint.iacr.org/2023/640>`__)
- [Robert 2023] “Breaking SIDH in Polynomial Time”, Robert, Damien,
  EUROCRYPT 2023

Introduction
^^^^^^^^^^^^

In the standard Diffie–Hellman key exchange, Alice and Bob end up with
the same secret because :math:`(g^a)^b = (g^b)^a`, which relies on
commutativity. With quaternions, we can try something similar by
multiplying on *opposite sides* of a public element, :math:`G`. Alice
publishes :math:`A\,G` and Bob publishes :math:`G\,B`. Then Alice
multiplies Bob’s public value on the left by her secret, :math:`A`, and
Bob multiplies Alice’s public value on the right by his secret,
:math:`B`. Associativity does the rest: :math:`A\,(G\,B) = (A\,G)\,B`.

First, here’s the exchange:

.. code:: ipython3

    rng = random.Random(7)
    
    # Public: a base element, G.
    G = Hu(1, 2, 3, 4)
    
    # Secret: Alice's A (norm 7) and Bob's B (norm 11).
    A = rng.choice(Hu.primes_of_norm(7, associates=True))
    B = rng.choice(Hu.primes_of_norm(11, associates=True))
    
    # Public: Alice multiplies G on the left, Bob multiplies it on the right.
    alice_public = A * G
    bob_public = G * B
    print(f"G            = {G}")
    print(f"A * G        = {alice_public}")
    print(f"G * B        = {bob_public}")
    
    # Shared: both arrive at A * G * B, by associativity.
    alice_shared = A * bob_public
    bob_shared = alice_public * B
    print(f"\nAlice computes A * (G * B) = {alice_shared}")
    print(f"Bob computes   (A * G) * B = {bob_shared}")
    print("They agree:", alice_shared == bob_shared)


.. parsed-literal::

    G            = (1+2i+3j+4k)
    A * G        = (10+5i-9j-2k)
    G * B        = (-4+7i+3j-16k)
    
    Alice computes A * (G * B) = (-40+10i-9j+23k)
    Bob computes   (A * G) * B = (-40+10i-9j+23k)
    They agree: True


What the Public Values Reveal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now let’s see what can be computed from the public values alone:
:math:`G`, :math:`A\,G`, and :math:`G\,B`. Since :math:`G` is nonzero,
it has an inverse among the quaternions over :math:`\mathbb{Q}`, which
``Hy`` computes exactly. With :math:`G^{-1}` in hand, one product gives
:math:`(A\,G)\,G^{-1}\,(G\,B) = A\,G\,B`, the shared element.

.. code:: ipython3

    # Using only the public values, G, A*G, and G*B:
    G_inv = G.to_hy().inverse()                       # exact rational inverse
    computed = alice_public.to_hy() * G_inv * bob_public.to_hy()      # A G G^-1 G B
    print("(A G) G^-1 (G B) =", Hu.from_hy(computed))
    print("The shared element:", Hu.from_hy(computed) == alice_shared)
    
    # The same idea recovers each of the secrets:
    print("A recovered from A*G:", Hu.from_hy(alice_public.to_hy() * G_inv) == A)
    print("B recovered from G*B:", Hu.from_hy(G_inv * bob_public.to_hy()) == B)
    
    # The norm of the shared element is the product of the three norms:
    print("\nN(shared) =", alice_shared.norm(), "= N(A)*N(G)*N(B) =",
          A.norm() * G.norm() * B.norm())
    print("Number of elements of norm 7:", Hu.count_of_norm(7))


.. parsed-literal::

    (A G) G^-1 (G B) = (-40+10i-9j+23k)
    The shared element: True
    A recovered from A*G: True
    B recovered from G*B: True
    
    N(shared) = 2310 = N(A)*N(G)*N(B) = 2310
    Number of elements of norm 7: 192


Observations
^^^^^^^^^^^^

- Since :math:`G` is invertible, the public values determine both
  secrets: :math:`A = (A\,G)\,G^{-1}` and :math:`B = G^{-1}\,(G\,B)`.
  That takes one multiplication each, with no factoring or searching.
- The norm of the shared element, :math:`N(A)\,N(G)\,N(B)`, is
  determined by the norms of the three factors, which are all public
  here. So a key derived from the norm alone wouldn’t depend on anything
  secret.
- The candidates for each secret are few: ``Hu.count_of_norm(7)`` is
  192, and ``Hu.count_of_norm(11)`` is 288.
- So this exchange is a toy, and the same observation applies to any
  scheme where the public data live in a ring in which the base can be
  inverted. Isogeny-based key exchange is a different kind of
  construction, and even there the details matter: SIDH was broken in
  2022 by classical attacks [Castryck & Decru 2022], [Maino et
  al. 2023], [Robert 2023].

Exact 3D Rotations
~~~~~~~~~~~~~~~~~~

A nonzero quaternion :math:`q` rotates a vector :math:`v = (x, y, z)` by
:math:`v \mapsto q\,v\,q^{-1}`, where :math:`v` is treated as the pure
quaternion :math:`xi + yj + zk`. (Scaling :math:`q` doesn’t change the
rotation, so :math:`q` doesn’t have to be a unit.) When :math:`q` and
:math:`v` are rational, the result is exactly rational, so there’s no
floating-point error to accumulate.

The Hurwitz unit :math:`\omega = \tfrac12(1 + i + j + k)` makes a nice
first example. It rotates by :math:`120^\circ` about the axis
:math:`(1, 1, 1)`.

Here’s a function that does the rotation, followed by :math:`\omega`
applied to a few vectors:

.. code:: ipython3

    def rotate(q, v):
        "Rotate the 3-vector v by the quaternion q (a Hu or Hy), exactly."
        qh = q.to_hy() if isinstance(q, Hu) else q
        r = qh * Hy.from_array([0, *v]) * qh.inverse()
        real, *vector = r.to_array()
        assert real == 0
        return tuple(vector)
    
    def show(v):
        return "(" + ", ".join(str(c) for c in v) + ")"
    
    omega = Hu('1/2', '1/2', '1/2', '1/2')
    
    print("omega takes (x, y, z) to (z, x, y):")
    for v in [(1, 1, -1), (1, -1, 1), (-1, 1, 1), (1, 2, 3)]:
        print(f"  {show(v):>12} -> {show(rotate(omega, v))}")
    
    # Three 120-degree turns make a full turn. Exactly, since omega**3 is -1:
    print("\nomega**3 =", omega**3)
    print("rotating by omega**3 fixes every vertex of the cube:",
          all(rotate(omega**3, v) == tuple(map(Fraction, v))
              for v in [(x, y, z) for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]))


.. parsed-literal::

    omega takes (x, y, z) to (z, x, y):
        (1, 1, -1) -> (-1, 1, 1)
        (1, -1, 1) -> (1, 1, -1)
        (-1, 1, 1) -> (1, -1, 1)
         (1, 2, 3) -> (3, 1, 2)
    
    omega**3 = (-1)
    rotating by omega**3 fixes every vertex of the cube: True


Next, a rotation that isn’t of finite order. The quaternion
:math:`1 + 2i` rotates about the :math:`x`-axis through an angle whose
cosine is :math:`-3/5` and whose sine is :math:`4/5`. Repeating a
rotation like that never brings a vector back to where it started, and
the denominators keep growing, but the length stays exactly the same, as
the second computation below shows. The last line checks that rotating
by :math:`q_2` and then by :math:`q_1` is the same as rotating by the
product :math:`q_1 q_2`.

.. code:: ipython3

    print("1 + 2i takes (0, 1, 0) to", show(rotate(Hu(1, 2, 0, 0), (0, 1, 0))))
    
    # Rotate 12 times by q = 1 + 2i + 2j. Its angle isn't a rational multiple of pi,
    # so the vector never returns, but its squared length, 14, never changes.
    q = Hu(1, 2, 2, 0)
    v = (1, 2, 3)
    for step in range(1, 13):
        v = rotate(q, v)
    print("\nafter 12 rotations:", show(v))
    print("squared length:", sum(c * c for c in v))
    
    # Rotating by q2 and then by q1 is rotating by q1 * q2:
    q1, q2, v = Hu(1, 2, 3, 1), Hu(2, -1, 1, 1), (1, 2, 3)
    print("\nrotate(q1*q2, v) == rotate(q1, rotate(q2, v)):",
          rotate(q1 * q2, v) == rotate(q1, rotate(q2, v)))


.. parsed-literal::

    1 + 2i takes (0, 1, 0) to (0, -3/5, 4/5)
    
    after 12 rotations: (-106330839439/282429536481, 953619448882/282429536481, -442757159597/282429536481)
    squared length: 14
    
    rotate(q1*q2, v) == rotate(q1, rotate(q2, v)): True


The Hurwitz Units and the Binary Tetrahedral Group
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 24 units of ``Hu`` are :math:`\pm 1, \pm i, \pm j, \pm k` and the
sixteen values :math:`\tfrac12(\pm 1 \pm i \pm j \pm k)`. They’re closed
under multiplication, and they form a group, the **binary tetrahedral
group**, which is isomorphic to :math:`SL(2,3)`. Modulo :math:`\pm 1`,
they act as the 12 rotations of a tetrahedron (the group :math:`A_4`),
as we saw with :math:`\omega` in the `previous
example <#exact-3d-rotations>`__.

The following computation checks closure, finds the orders of the
elements and the center, and counts the distinct rotations:

.. code:: ipython3

    units = list(Hu.units().values())
    unit_set = set(units)
    
    print("closed under multiplication:", all(a * b in unit_set for a in units for b in units))
    
    def order(u):
        m, x = 1, u
        while x != Hu(1):
            x, m = x * u, m + 1
        return m
    
    print("element orders {order: count}:", dict(sorted(Counter(order(u) for u in units).items())))
    print("center:", [str(u) for u in units if all(u * w == w * u for w in units)])
    
    # Modulo +-1, the units give 12 distinct rotations of 3-space:
    rotations = {rotate(u, (1, 2, 3)) for u in units}
    print("distinct rotations of the vector (1, 2, 3):", len(rotations))


.. parsed-literal::

    closed under multiplication: True
    element orders {order: count}: {1: 1, 2: 1, 3: 8, 4: 6, 6: 8}
    center: ['(1)', '(-1)']
    distinct rotations of the vector (1, 2, 3): 12


``finalg`` has a catalogue of small groups (all groups of order 40 or
less; see ``finalg.small_groups``), and the binary tetrahedral group is
one of them, named ``SL2_3``. In the next cell, we build the Cayley
table of the Hurwitz units, hand it to ``finalg``, and ask which small
group it is. (If ``finalg`` isn’t installed, the cell just prints a note
and moves on.)

.. code:: ipython3

    try:
        with contextlib.redirect_stdout(io.StringIO()):     # finalg prints a banner on import
            from finalg import make_finite_algebra, identify_small_group
    except ImportError:
        print("finalg is not installed (https://github.com/alreich/finalg); skipping this check.")
    else:
        index = {u: i for i, u in enumerate(units)}
        table = [[index[a * b] for b in units] for a in units]
        hurwitz_units = make_finite_algebra("Hu24", "The 24 Hurwitz units",
                                            [str(u) for u in units], table)
        small = identify_small_group(hurwitz_units)
        print(type(hurwitz_units).__name__, "of order", hurwitz_units.order)
        print("identified as:", small.name, "--", small.description)


.. parsed-literal::

    Group of order 24
    identified as: SL2_3 -- SL(2,3) = Q8 X| Z3, the binary tetrahedral group


Example: Sums of Four Squares as Hurwitz Norms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lagrange’s four-square theorem says that every non-negative integer is a
sum of four squares. ``Hu.four_squares(n)`` finds the squares by
constructing a quaternion of norm :math:`n` out of Hurwitz primes, so
even a very large :math:`n` costs about as much as factoring it.

Jacobi’s formula, :math:`r_4(n) = 8\sum_{d \mid n,\ 4 \nmid d} d`,
counts the ways to write :math:`n` as a sum of four squares, where order
and signs matter. That’s the number of Lipschitz integers of norm
:math:`n`. And, ``Hu.count_of_norm(n)`` counts *all* of the Hurwitz
integers of norm :math:`n`, which includes the ones with half-integer
coordinates. After a few examples, the cell below checks both counts for
:math:`n = 1, \dots, 30`.

.. code:: ipython3

    for n in [7, 310, 2026, 10**30 + 7]:
        a, b, c, d = Hu.four_squares(n)
        print(f"{n} = {a}^2 + {b}^2 + {c}^2 + {d}^2   ({a*a + b*b + c*c + d*d == n})")
    
    def r4(n):                                  # Jacobi's four-square formula
        return 8 * sum(d for d in range(1, n + 1) if n % d == 0 and d % 4 != 0)
    
    print("\nJacobi's count matches the Lipschitz integers of norm n (n = 1..30):",
          all(sum(h.is_lipschitz() for h in Hu.of_norm(n)) == r4(n) for n in range(1, 31)))
    print("count_of_norm(n) matches the enumeration of Hurwitz integers (n = 1..30):",
          all(len(Hu.of_norm(n)) == Hu.count_of_norm(n) for n in range(1, 31)))


.. parsed-literal::

    7 = 2^2 + 1^2 + 1^2 + 1^2   (True)
    310 = 13^2 + 10^2 + 5^2 + 4^2   (True)
    2026 = 23^2 + 23^2 + 22^2 + 22^2   (True)
    1000000000000000000000000000007 = 745714311982110^2 + 563651243969653^2 + 327381580838907^2 + 137944701253543^2   (True)
    
    Jacobi's count matches the Lipschitz integers of norm n (n = 1..30): True
    count_of_norm(n) matches the enumeration of Hurwitz integers (n = 1..30): True


Division with Remainder and the Euclidean Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Why use the Hurwitz integers, rather than just the Lipschitz integers
(the quaternions with integer coordinates)? One reason is that the
Lipschitz integers don’t have division with remainder. If we divide
:math:`1 + i + j + k` by 2, the best Lipschitz quotient leaves a
remainder whose norm is 4, which is not less than :math:`N(2) = 4`.

With the half-integer points of the Hurwitz integers included, division
with remainder always works, and so does the Euclidean algorithm,
including its extended version, which produces a Bézout identity. Since
multiplication isn’t commutative, each of these has a left version and a
right version.

First, the division:

.. code:: ipython3

    x, y = Hu(1, 1, 1, 1), Hu(2)
    
    # The best Lipschitz quotient, by brute force over a box of integer quotients:
    best = min((x - y * Hu(a, b, c, d)).norm()
               for a in range(-2, 3) for b in range(-2, 3) for c in range(-2, 3) for d in range(-2, 3))
    print(f"smallest remainder norm with a Lipschitz quotient: {best}   (N(2) = {y.norm()})")
    
    # With Hurwitz integers, the division here is exact:
    quotient, remainder = x.divmod_left(y)
    print("Hurwitz:  x = y * q + r   with q =", quotient, " r =", remainder)


.. parsed-literal::

    smallest remainder norm with a Lipschitz quotient: 4   (N(2) = 4)
    Hurwitz:  x = y * q + r   with q = (1/2+1/2i+1/2j+1/2k)  r = (0)


Next, a greatest common left divisor, computed with ``xgcld``, along
with the Bézout coefficients, which multiply on the right. We start with
two values that have a known common left factor, :math:`c`, so we know
what the answer should be, up to a unit on the right.

.. code:: ipython3

    c = Hu(1, 1, 1, 0)                                    # a common left factor, norm 3
    a0, b0 = Hu(-2, 1, 3, 3), Hu(3, -3, -1, -3)
    a, b = c * a0, c * b0
    
    g, s, t = a.xgcld(b)
    print("gcld(a, b)            =", g)
    print("c, up to a right unit =", c.canonical_associate("right"))
    print("Bezout:  a*s + b*t == g:", a * s + b * t == g)


.. parsed-literal::

    gcld(a, b)            = (3/2+1/2i-1/2j+1/2k)
    c, up to a right unit = (3/2+1/2i-1/2j+1/2k)
    Bezout:  a*s + b*t == g: True


Summary
-------

The three pages of this documentation follow the structure of the
package.

**Classical hypercomplex numbers** (```hyprat``: Exact Rational
Hypercomplex Numbers <classic_hy.html>`__)

- ``Hy`` is a single, immutable class for exact rational hypercomplex
  numbers at *any* Cayley-Dickson rank: reals (0), complex (1),
  quaternions (2), octonions (3), sedenions (4), and beyond. Every
  coordinate is a ``Fraction``, so nothing is ever rounded.
- Arithmetic (``+ - * /``), conjugation, norm and inversion follow the
  classical recursive formulas. Zero divisors appear from rank 4 on, and
  ``Hy`` reports them honestly: an element with a zero divisor has no
  inverse.
- ``str()``, ``repr()`` and ``Hy.parse()`` (also ``Hy.from_string()``)
  give readable round-trip text forms. ``from_array()`` and
  ``to_array()`` convert to and from flat lists of numbers or fraction
  strings. ``to_latex()`` and the matrix representation cover display
  and linear algebra.
- ``Hy.random(rank, ...)`` makes random values of any rank.
  ``Hy.seed(...)``, or the per-call ``seed=`` and ``rng=`` keywords,
  make them reproducible, and the randomized tests check the algebraic
  laws on thousands of values.

**Split-hypercomplex numbers and other packages** (`Split-Hypercomplex
Numbers and Interoperability <split_hy_and_interop.html>`__)

- The same class ``Hy`` also builds the *split* variants, by choosing a
  **signature**, one nonzero rational :math:`\mu` per doubling level.
  The presets ``Hy.SPLIT_COMPLEX``, ``Hy.SPLIT_QUATERNION`` and
  ``Hy.SPLIT_OCTONION`` cover the usual cases, and any nonzero rational
  is allowed. Split numbers have an indefinite norm, so zero divisors
  appear already at rank 1.
- Rank 2 values convert to and from SymPy quaternions (exactly),
  numpy-quaternion and quaternionic (floating point), and the Gaussian
  and rational Gaussian integers of ``gint`` convert to and from rank 1
  values. Floats are only as exact as the package that holds them.
  ``Hy`` stays exact. None of these packages is required: each is
  imported only when one of the conversion methods is called.

**Hurwitz integers and number theory** (this page)

- ``Hu`` holds the Hurwitz integers as exact Python integers (doubled
  coordinates, so no fractions) and adds number theory: the 24 units,
  division with remainder on the left and on the right, left and right
  greatest common divisors with Bézout coefficients, associates, primes,
  content, and factorization into primes in any order of norms, with the
  unit migration that non-commutativity requires.
- ``Hu.of_norm``, ``Hu.count_of_norm``, ``Hu.primes_of_norm``,
  ``Hu.with_norm`` and ``Hu.four_squares`` enumerate and construct
  elements of a given norm, and give Lagrange’s four-square theorem as a
  by-product, even for very large integers.
- ``hyprat.intfactor`` does the integer arithmetic that this needs, in
  pure Python: ``primes_upto``, ``is_probable_prime`` (certain below
  :math:`3.3 \times 10^{24}`, Baillie-PSW above), ``factorint`` (trial
  division and Pollard-Brent, with SymPy optional) and
  ``sqrt_mod_prime``.
- The examples used these tools to factor a product of Hurwitz primes,
  to try a Diffie-Hellman-style key exchange, to compute exact 3D
  rotations, to identify the binary tetrahedral group, to write integers
  as sums of four squares, and to run the Euclidean algorithm.

See the package’s `GitHub
repository <https://github.com/alreich/hyper_rationals>`__ and the `Read
the Docs <https://hyper-rationals.readthedocs.io/>`__ page for the full
API reference and source.
