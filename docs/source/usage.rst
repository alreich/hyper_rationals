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

See :doc:`api` for the full reference.
