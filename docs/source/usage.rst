Usage
=====

Constructing values
--------------------

``Hy`` always stores exactly two components, ``real`` and ``imag``,
which are either both plain ``Fraction``\ s (a rank-1, complex-like
value) or both ``Hy`` instances of the same rank (a higher-rank
value). The constructor normalizes whatever it is given so this
invariant always holds:

.. code-block:: python

    from fractions import Fraction
    from hyprat import Hy

    Hy('5/2', '-16/5')      # a rational complex number, 5/2 - 16/5 j
    Hy(2.5, -3.2)           # same value, from floats
    Hy(3)                   # 3 + 0j  (a bare number, still rank 1)

    h1 = Hy('1', '2')
    h2 = Hy('3', '4')
    q = Hy(h1, h2)          # a rational quaternion: 1 + 2i + 3j + 4k

    h3 = Hy(Hy(1, 2), Hy(3, 4))
    h4 = Hy(Hy(5, 6), Hy(7, 8))
    o = Hy(h3, h4)          # a rational octonion

Two more ways to build a value: from a flat list of coefficients, or
drawn at random.

.. code-block:: python

    # From a flat list of 2**rank coefficients -- ints, floats,
    # Fractions, and fraction strings like '5/2' may be freely mixed:
    Hy.from_array([1, '2/3', 3.5, '-1/4'])   # a quaternion

    # A random value of a given rank, with a seed for reproducibility:
    Hy.random(2, seed=0)                      # a random quaternion

Arithmetic
----------

``+``, ``-``, ``*`` and ``/`` are all defined recursively via the
Cayley-Dickson construction, so they work uniformly at every rank:

.. code-block:: python

    z1, z2 = Hy('1', '2'), Hy('3', '-1')
    z1 + z2, z1 - z2, z1 * z2, z1 / z2

Multiplication is non-commutative for quaternions and octonions, and
non-associative for octonions, exactly as it should be:

.. code-block:: python

    i = Hy(Hy(0, 1), Hy(0, 0))
    j = Hy(Hy(0, 0), Hy(1, 0))
    i * j   # ->  k
    j * i   # -> -k

Parsing and printing
---------------------

``str()`` renders a value the way Python renders ``complex`` numbers
for rank 1 (using ``j``), and the customary ``a+bi+cj+dk`` notation for
rank 2, and ``e1 .. e_{2^rank - 1}`` imaginary units for rank >= 3.
``Hy.parse`` is the inverse operation:

.. code-block:: python

    str(Hy('5/2', '16/5'))          # '(5/2+16/5j)'
    Hy.parse('5/2+16/5j')           # == Hy('5/2', '16/5')

    str(Hy(Hy(1, 2), Hy(3, 4)))     # '(1+2i+3j+4k)'
    Hy.parse('1+2i+3j+4k')          # a rational quaternion

``repr()`` returns Python source that reconstructs an equal value,
e.g. ``Hy('5/2', '-16/5')``.

A value's coefficients can also be read out as a flat list with
``to_array()`` -- the inverse of ``Hy.from_array()`` above -- either as
``Fraction``\ s, or as strings (handy for JSON or other text-based
serialization):

.. code-block:: python

    q = Hy(Hy(1, 2), Hy(3, 4))
    q.to_array()               # [Fraction(1), Fraction(2), Fraction(3), Fraction(4)]
    q.to_array(as_str=True)    # ['1', '2', '3', '4']
    Hy.from_array(q.to_array(as_str=True)) == q   # True

Random values
--------------

``Hy.random(rank)`` draws a random rank-``rank`` value: each of its
``2**rank`` coefficients is an independent ``Fraction(n, d)``, with
``n`` uniform over ``[lo, hi]`` (default ``[-9, 9]``) and ``d`` uniform
over ``[1, dmax]`` (default ``[1, 6]``). This is handy for examples,
demos, and property-based ("fuzz") testing.

There are a few ways to control reproducibility:

.. code-block:: python

    # A one-off seed, scoped to just this call:
    Hy.random(2, seed=7)

    # Hy.seed(...) fixes a shared default RNG for everything that
    # follows, so bare Hy.random(rank) calls become reproducible too:
    Hy.seed(2026)
    Hy.random(2)

    # Or bring your own random.Random for full control:
    import random
    Hy.random(2, rng=random.Random(123))

``rank`` must be a positive int (every ``Hy`` has rank >= 1 by
construction, so there's no rank-0 ``Hy``), and ``rng``/``seed`` are
mutually exclusive.

See :doc:`api` for the full reference.
