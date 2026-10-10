Usage
=====

.. code:: ipython3

    >>> from hyprat import Hy
    >>> from IPython.display import display, Math  # for LaTeX output

Constructing Hypercomplex Numbers, Recursively, via the Class ``Hy``
--------------------------------------------------------------------

The ``Hy`` class represents a “tower” of multi-dimensional numbers, and
every ``Hy`` has only two components: ``real`` and ``imag``.

At the lowest level of the tower, the two components are rational
numbers, represented by the ``fractions.Fraction`` class.

Every instance of a ``Hy`` has a property called ``rank``, which will be
a non-negative integer (rank = 0, 1, 2, …). Although the ``Fraction``
class is different than the ``Hy`` class, for consistency it is
considered here to have rank 0.

If a ``Hy`` has rank :math:`n`, where :math:`n=1, 2, \dots`, then its
two components, ``real`` and ``imag``, will both always have rank
:math:`n-1`. The ``Hy`` constructor will normalize whatever it is given
so this invariant holds. That is, if ``real`` has rank :math:`m` and
``imag`` has rank :math:`n`, and :math:`m \ne n`, then the constructor
will coerce the lower rank ``Hy`` into an instance of the higher rank
``Hy`` before constructing the new, higher rank ``Hy``.

So, starting at the lowest level of the tower, as noted above, a ``Hy``
made up of two rational numbers represents a **rational complex number**
and has rank 1, as shown in the following table.

============ ==== ===============
Hypercomplex Rank Dimension
============ ==== ===============
Rational     0    :math:`1 = 2^0`
Complex      1    :math:`2 = 2^1`
Quaternion   2    :math:`4 = 2^2`
Octonion     3    :math:`8 = 2^3`
in general   n    :math:`d = 2^n`
============ ==== ===============

Rational Complex Numbers
~~~~~~~~~~~~~~~~~~~~~~~~

Combine two rational numbers (``Fraction``\ s) to create a **rational
complex number**. Floats, ints, and strings can be used for the two
rational numbers, and they can be mixed:

.. code:: ipython3

    >>> z1 = Hy("2/3", 1.5)
    >>> print(f"{z1 = } has rank {z1.rank}\n")
    
    >>> print(f"{str(z1) = }")


.. parsed-literal::

    z1 = Hy('2/3', '3/2') has rank 1
    
    str(z1) = '(2/3+3/2j)'


Rational Quaternions
~~~~~~~~~~~~~~~~~~~~

Combine two rational complex numbers (``Hy``\ s of rank 1) to create a
**rational quaternion** (rank 2):

.. code:: ipython3

    >>> z2 = Hy(4, "-1/7")  # another rank 1 Hy (complex number)
    
    >>> q1 = Hy(z1, z2)  # rational quaternion
    >>> print(f"{q1 = } has rank {q1.rank}\n")
    
    >>> print(f"{str(q1) = }")


.. parsed-literal::

    q1 = Hy(Hy('2/3', '3/2'), Hy('4', '-1/7')) has rank 2
    
    str(q1) = '(2/3+3/2i+4j-1/7k)'


Rational Octonions
~~~~~~~~~~~~~~~~~~

Combine two rational quaternions (``Hy``\ s of rank 2) to create a
**rational octonion** (rank 3):

.. code:: ipython3

    >>> q2 = Hy(Hy('4/3', '-9/4'), Hy('-6/5', '2/5'))  # another rank 2 Hy (quaternion)
    
    >>> o1 = Hy(q1, q2)  # rational octonion
    >>> print(f"{o1 = } has rank {o1.rank}\n")
    
    >>> print(f"{str(o1) = }")


.. parsed-literal::

    o1 = Hy(Hy(Hy('2/3', '3/2'), Hy('4', '-1/7')), Hy(Hy('4/3', '-9/4'), Hy('-6/5', '2/5'))) has rank 3
    
    str(o1) = '(2/3+3/2i+4j-1/7k+4/3L-9/4iL-6/5jL+2/5kL)'


Other Ways to Construct Hypercomplex Numbers
--------------------------------------------

There are three more ways to build a Hy:

-  from a flat list of coefficients (``to_array`` and ``Hy.from_array``)
-  from a string representation (``Hy.parse``)
-  randomly generated (``Hy.random``)

Examples follow:

From a Flat List of Coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Hy``\ s can be converted both **to** and **from** flat lists of
coefficients.

**To a Flat List**:

.. code:: ipython3

    >>> q1_coef = q1.to_array(as_str=True)  # setting as_str to False (default) returns Fractions
    >>> q1_coef




.. parsed-literal::

    ['2/3', '3/2', '4', '-1/7']



**From a Flat List**:

.. code:: ipython3

    >>> q1_copy = Hy.from_array(q1_coef)
    >>> q1_copy == q1




.. parsed-literal::

    True



From a String Representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    >>> o1_str = str(o1)
    
    >>> o1_copy = Hy.parse(o1_str)
    >>> o1_copy == o1




.. parsed-literal::

    True



From a Random Number Generator (RNG)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    >>> print(f"   complex: {Hy.random(1)}")
    >>> print(f"quaternion: {Hy.random(2)}")
    >>> print(f"  octonion: {Hy.random(3)}")
    >>> print(f"  sedenion: {Hy.random(4)}")


.. parsed-literal::

       complex: (7-8/5j)
    quaternion: (1-3i-3/5j-8/3k)
      octonion: (4+7/5i+6/5j-k+8L+2iL-9jL+5kL)
      sedenion: (9/2-3/5e1+6/5e2-4/3e3-5e4-e5-8/5e6+e7+7e8+1/2e9-6e10+6e11+5/3e12-e13-7/5e14+1/2e15)


There’s more below on random generation of hypercomplex numbers.

Arithmetic
----------

``+``, ``-``, ``*`` and ``/`` are all defined recursively via the
Cayley-Dickson construction, so they work uniformly at every rank, with
the exception of division at ranks :math:`\ge 4` (i.e., sedenions,
pathions, …).

.. code:: ipython3

    >>> print(f"{z1 + z2 = }")
    >>> print(f"{z1 - z2 = }")
    >>> print(f"{z1 * z2 = }")
    >>> print(f"{z1 / z2 = }")


.. parsed-literal::

    z1 + z2 = Hy('14/3', '19/14')
    z1 - z2 = Hy('-10/3', '23/14')
    z1 * z2 = Hy('121/42', '124/21')
    z1 / z2 = Hy('721/4710', '896/2355')


Multiplication is non-commutative for quaternions and octonions, and
non-associative for octonions, exactly as it should be:

.. code:: ipython3

    >>> i = Hy(Hy(0, 1), Hy(0, 0))
    >>> j = Hy(Hy(0, 0), Hy(1, 0))
    >>> print(f"{str(i * j) = }")
    >>> print(f"{str(j * i) = }")


.. parsed-literal::

    str(i * j) = '(k)'
    str(j * i) = '(-k)'


More on Random Hypercomplex Numbers
-----------------------------------

``Hy.random(rank)``, where ``rank`` is a positive integer, draws a
random rank-``rank`` value: each of its ``2**rank`` coefficients is an
independent ``Fraction(n, d)``, with ``n`` uniform over ``[lo, hi]``
(default ``[-9, 9]``) and ``d`` uniform over ``[1, dmax]`` (default
``[1, 6]``). The defaults for ``lo``, ’hi\ ``, and``\ dmax\` are set to
small values to make examples, demos, and tests easy to read.

There are a few ways to control reproducibility of random output.

.. code:: ipython3

    # A one-off seed, scoped to just a single call:
    >>> Hy.random(2, seed=7)
    
    # Hy.seed(...) sets a shared default RNG for everything that
    # follows, so plain Hy.random(rank) calls become reproducible also:
    >>> Hy.seed(2026)
    >>> Hy.random(2)




.. parsed-literal::

    Hy(Hy('-2', '7/5'), Hy('-3', '2'))



.. code:: ipython3

    # Or use your own random.Random for full control:
    >>> import random
    >>> Hy.random(2, rng=random.Random(123))




.. parsed-literal::

    Hy(Hy('-8/3', '-7/4'), Hy('-1', '-2'))



Units
-----

``Hy.units(rank)`` returns a dictionary of all units for the particular,
``rank``, where the keys are the string representation of the units and
the values are the ``Hy``\ s.

.. code:: ipython3

    >>> Hy.units(1)




.. parsed-literal::

    {'1': Hy('1', '0'),
     '-1': Hy('-1', '0'),
     'j': Hy('0', '1'),
     '-j': Hy('0', '-1')}



.. code:: ipython3

    >>> Hy.units(2).keys()




.. parsed-literal::

    dict_keys(['1', '-1', 'i', '-i', 'j', '-j', 'k', '-k'])



``rank == 0`` is the one case where the “hypercomplex value” in question
is actually a ``Fraction`` rather than a ``Hy``, so ``Hy.units(0)``
returns ``{'1': Fraction(1, 1), '-1': Fraction(-1, 1)}``.

.. code:: ipython3

    >>> Hy.units(0)




.. parsed-literal::

    {'1': Fraction(1, 1), '-1': Fraction(-1, 1)}



``some_hy.is_unit()`` answers the corresponding membership question for
a single value, without building the whole dict:

.. code:: ipython3

    >>> Hy(0, 1).is_unit()   # j is a unit




.. parsed-literal::

    True



.. code:: ipython3

    >>> Hy(1, 1).is_unit()   # not a unit




.. parsed-literal::

    False



LaTeX rendering
---------------

``some_hy.latex()`` renders a value as a LaTeX math expression, which is
useful in a Jupyter notebook.

Two keyword-only options are available:

-  ``vinculum`` controls how a non-integer coefficient’s fraction bar is
   typeset: ``"horizontal"`` (the default) uses ``\frac{num}{den}``;
   ``"diagonal"`` uses a plain slash, ``num/den``:
-  ``mode`` controls whether the result is wrapped in LaTeX math
   delimiters: ``"plain"`` (the default) returns the bare expression,
   ``"inline"`` wraps it in ``$...$``, and ``"display"`` wraps it in
   ``\[...\]``.

Basis units render the same way `the API reference <api.rst>`__
describes for ``str()``: ``j`` at rank 1, ``i``/``j``/``k`` at rank 2,
and ``i``/``j``/``k``/``L``/``iL``/``jL``/``kL`` at rank 3, all
unchanged, except that the ``e1, e2, ...`` labels used from rank 4 up
are subscripted (rendered as ``e_{1}``, ``e_{2}``, …).

See `api <api.rst>`__ for the full reference.

.. code:: ipython3

    >>> Hy('5/2', '-16/5').latex()




.. parsed-literal::

    '\\frac{5}{2}-\\frac{16}{5}j'



.. code:: ipython3

    >>> Hy('5/2', '-16/5').latex(vinculum='diagonal')




.. parsed-literal::

    '5/2-16/5j'



The following shows how to display ``Hy``\ s in a Jupyter notebook:

.. code:: ipython3

    >>> from IPython.display import display, Math
    
    >>> print(f"{z1 = }")
    >>> display(Math(z1.latex()))


.. parsed-literal::

    z1 = Hy('2/3', '3/2')



.. math::

    \displaystyle \frac{2}{3}+\frac{3}{2}j


.. code:: ipython3

    >>> print(f"{q1 = }")
    >>> display(Math(q1.latex()))


.. parsed-literal::

    q1 = Hy(Hy('2/3', '3/2'), Hy('4', '-1/7'))



.. math::

    \displaystyle \frac{2}{3}+\frac{3}{2}i+4j-\frac{1}{7}k


.. code:: ipython3

    >>> print(f"{o1 = }")
    >>> display(Math(o1.latex()))


.. parsed-literal::

    o1 = Hy(Hy(Hy('2/3', '3/2'), Hy('4', '-1/7')), Hy(Hy('4/3', '-9/4'), Hy('-6/5', '2/5')))



.. math::

    \displaystyle \frac{2}{3}+\frac{3}{2}i+4j-\frac{1}{7}k+\frac{4}{3}L-\frac{9}{4}iL-\frac{6}{5}jL+\frac{2}{5}kL


.. code:: ipython3

    >>> s1 = Hy.random(4)  # a random sedenion
    >>> print(f"{s1 = }")
    >>> display(Math(s1.latex()))


.. parsed-literal::

    s1 = Hy(Hy(Hy(Hy('9/5', '6/5'), Hy('5/2', '-9/5')), Hy(Hy('-7', '0'), Hy('5', '1'))), Hy(Hy(Hy('1/2', '1'), Hy('2/3', '1/2')), Hy(Hy('7/6', '-7/6'), Hy('1', '8/5'))))



.. math::

    \displaystyle \frac{9}{5}+\frac{6}{5}e_{1}+\frac{5}{2}e_{2}-\frac{9}{5}e_{3}-7e_{4}+5e_{6}+e_{7}+\frac{1}{2}e_{8}+e_{9}+\frac{2}{3}e_{10}+\frac{1}{2}e_{11}+\frac{7}{6}e_{12}-\frac{7}{6}e_{13}+e_{14}+\frac{8}{5}e_{15}


Matrix representation
---------------------

``some_hy.to_matrix()`` returns the *regular representation* of a value
as a ``2**rank x 2**rank`` NumPy array of exact ``Fraction``\ s – this
is the classical “complex numbers as 2x2 real matrices” / “quaternions
as 4x4 real matrices” construction, computed directly from ``Hy``\ ’s
own multiplication rather than from a derived formula: column ``i`` is
``self * units[i]``. ``Hy.from_matrix()`` is the inverse. This needs
NumPy (``pip install numpy``), which is not a runtime dependency of
``hyprat`` and is imported lazily.

For rank 0-2 (real, complex, quaternion) ``to_matrix()`` is a genuine
algebra isomorphism onto a matrix subalgebra:

::

   M(x) @ M(y) == M(x * y)

.. code:: ipython3

    >>> import numpy as np
    >>> z1.to_matrix()




.. parsed-literal::

    array([[Fraction(2, 3), Fraction(-3, 2)],
           [Fraction(3, 2), Fraction(2, 3)]], dtype=object)



.. code:: ipython3

    >>> Hy.from_matrix(z1.to_matrix()) == z1




.. parsed-literal::

    True



The homomorphism property holds at rank 1 and rank 2:

.. code:: ipython3

    >>> a, b = Hy.random(2, seed=10), Hy.random(2, seed=11)
    >>> np.array_equal(np.dot(a.to_matrix(), b.to_matrix()), (a * b).to_matrix())




.. parsed-literal::

    True



It’s also a nice cross-check on ``.norm()``: the determinant of the
regular representation equals the (squared) norm raised to a power that
depends on rank (``norm(x) ** (2 ** (rank - 1))``):

.. code:: ipython3

    >>> import sympy
    >>> det_q1 = sympy.Matrix(q1.to_matrix().tolist()).det()
    >>> det_q1 == sympy.Rational(q1.norm().numerator, q1.norm().denominator) ** 2




.. parsed-literal::

    True



Rank >= 3 (octonions and beyond): not a homomorphism
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Octonions (and everything beyond) are **not associative**, so ordinary
matrix multiplication – which *is* associative – cannot faithfully
represent their multiplication. ``to_matrix()``/``from_matrix()`` raise
a ``ValueError`` at rank >= 3 by default, to flag this rather than let
it be a silent trap:

.. code:: ipython3

    >>> try:
    ...     o1.to_matrix()
    ... except ValueError as e:
    ...     print(e)


.. parsed-literal::

    rank-3 values (octonions and beyond) are not associative, so their regular-representation matrix is a faithful *linear* embedding but NOT an algebra homomorphism: M(x) @ M(y) != M(x * y) in general, and (from rank 4 up) M(x) can even be singular for nonzero x. Pass allow_nonassociative=True to build the matrix anyway, or see the 'Matrix representation' section of the docs for faithful alternatives (e.g. Zorn vector matrices).


Passing ``allow_nonassociative=True`` builds the matrix anyway – it’s
still a faithful *linear* embedding (useful for some purposes), just not
a ring homomorphism:

.. code:: ipython3

    >>> Mo = o1.to_matrix(allow_nonassociative=True)
    >>> Mo.shape




.. parsed-literal::

    (8, 8)



.. code:: ipython3

    >>> c, d = Hy.random(3, seed=20), Hy.random(3, seed=21)
    >>> lhs = np.dot(c.to_matrix(allow_nonassociative=True), d.to_matrix(allow_nonassociative=True))
    >>> rhs = (c * d).to_matrix(allow_nonassociative=True)
    >>> np.array_equal(lhs, rhs)   # False in general -- non-associativity breaks the homomorphism




.. parsed-literal::

    False



Starting at rank 4 (sedenions), the algebra even has *zero divisors* –
nonzero values whose product is zero – so ``to_matrix()`` can be
singular for a nonzero input. For example, with this library’s basis
ordering, ``(e1 + e10) * (e4 - e15) == 0``:

.. code:: ipython3

    >>> units4 = [u for name, u in Hy.units(4).items() if not name.startswith('-')]
    >>> a4 = units4[1] + units4[10]
    >>> b4 = units4[4] - units4[15]
    >>> a4 * b4




.. parsed-literal::

    Hy(Hy(Hy(Hy('0', '0'), Hy('0', '0')), Hy(Hy('0', '0'), Hy('0', '0'))), Hy(Hy(Hy('0', '0'), Hy('0', '0')), Hy(Hy('0', '0'), Hy('0', '0'))))



.. code:: ipython3

    >>> Ma4 = a4.to_matrix(allow_nonassociative=True, as_float=True)
    >>> round(np.linalg.det(Ma4), 6)   # singular, even though a4 != 0




.. parsed-literal::

    np.float64(0.0)



See the `bibliography <bibliography.rst>`__ for references on faithful,
non-matrix-multiplication representations of octonions (e.g. Zorn vector
matrices), which sidestep this associativity obstruction.

Split-hypercomplex numbers and other signatures
-----------------------------------------------

Each doubling step of the Cayley-Dickson construction has a parameter
``mu``, a nonzero rational, that enters the product as

::

   (a, b)(c, d) = (a c + mu conj(d) b,  d a + b conj(c))

With ``mu = -1`` (the default, and what every example above uses) the
tower is the classical one: complex numbers, quaternions, octonions, and
so on. A step with ``mu = +1`` is the *split* version of that step, and
any other nonzero rational gives the corresponding generalized algebra.

``Hy(real, imag, mu=...)`` sets ``mu`` for the **top** level of the new
value; its components already carry the lower levels. The full
*signature* of a value, one ``mu`` per level from the lowest up, is
``some_hy.signs``, and the top one is ``some_hy.mu``.

Split-complex numbers
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


Split-quaternions and split-octonions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Which levels are split is given by the signature, ordered from the
lowest level to the highest. The common cases have named presets, usable
as ``signs="split-quaternion"`` or through the class constants:

.. code:: ipython3

    >>> Hy.SPLIT_COMPLEX, Hy.SPLIT_QUATERNION, Hy.SPLIT_OCTONION




.. parsed-literal::

    ((1,), (-1, 1), (-1, -1, 1))



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



Any nonzero rational ``mu``
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


Mixing, equality and parsing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Plain numbers mix freely with any signature, but values from different
algebras never combine. Arithmetic raises ``ValueError``, while ``==``
is just ``False``:

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



Zero divisors
~~~~~~~~~~~~~

For ranks 1 to 3 the zero divisors are exactly the null elements, which
is what ``is_null()`` detects. From rank 4 up (the sedenions and beyond)
there are also zero divisors with *nonzero* norm. In the classical
sedenions every zero divisor is like this, so ``is_null()`` never sees
them.

``is_zero_divisor()`` is the exact test, for any rank and signature. It
checks whether the map ``y -> x * y`` is singular, using exact rational
elimination on that map’s matrix. ``0`` itself is not counted.

.. code:: ipython3

    >>> zd = Hy.from_array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0])   # e1 + e10, a sedenion
    >>> zd.is_null(), zd.is_zero_divisor()




.. parsed-literal::

    (False, True)



``annihilator()`` returns an exact basis of the values it kills. With
the default ``kind="left"`` these are the ``y`` with ``x * y == 0``, and
with ``kind="right"`` the ``y`` with ``y * x == 0``:

.. code:: ipython3

    >>> zd_ann = zd.annihilator()
    >>> len(zd_ann), all(zd * y == 0 for y in zd_ann)




.. parsed-literal::

    (4, True)



There is no left/right argument for ``is_zero_divisor()`` itself: a
value is a zero divisor on one side exactly when it is on the other, and
exactly when its conjugate is. The two annihilators always have the same
dimension, but they can be different sets of values. In the classical sedenion cases checked they coincided, while in a split algebra they often do not:

.. code:: ipython3

    >>> zd_split = Hy.from_array([1, 0, 0, 0, 1] + [0] * 11, signs=(-1, -1, 1, -1))
    >>> left, right = zd_split.annihilator("left"), zd_split.annihilator("right")
    >>> len(left), len(right), left == right




.. parsed-literal::

    (8, 8, False)

Matrices and other packages
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The split-quaternions are isomorphic to the 2x2 real matrices, and
``to_matrix()`` still gives a genuine homomorphism at ranks 1 and 2 for
any signature. A matrix does not record its signature, so read it back
with the matching ``signs=``:

.. code:: ipython3

    >>> import numpy as np
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

Interoperability with Other Quaternion Packages
-----------------------------------------------

A rational quaternion (a ``Hy`` of rank 2) can be converted to and from
the quaternion types of three third-party packages. None of them is
required by ``hyprat``; each is imported only when a conversion method
that needs it is called, so install whichever you use:

::

   pip install sympy numpy-quaternion quaternionic

-  **SymPy**, ``sympy.algebras.quaternion.Quaternion``:
   ``some_hy.to_sympy()`` and ``Hy.from_sympy(sq)``. SymPy, like ``Hy``,
   represents rational numbers exactly, so this pair round-trips with no
   rounding at all.
-  **numpy-quaternion**, imported as ``quaternion``:
   ``some_hy.to_numpy_quaternion()`` and
   ``Hy.from_numpy_quaternion(nq)``. It stores
   ``numpy.quaternion(w, x, y, z)`` as ``float64``.
-  **quaternionic**: ``some_hy.to_quaternionic()`` and
   ``Hy.from_quaternionic(qa)``. It stores an array ``[w, x, y, z]`` as
   ``float64``.

All three packages use the same Hamilton convention as ``Hy``
(``i*j == k``) and the coordinates always correspond to ``1, i, j, k``,
in that order, so the ``from_*`` methods need no reordering.

A rank 1 ``Hy`` (a rational complex number) is embedded as the
quaternion ``a + b*i`` when it is converted. Ranks 3 and up (octonions,
sedenions, …) have no quaternion equivalent, so converting them raises a
``ValueError``.

Converting a ``Hy`` to Another Package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    >>> import numpy as np
    >>> import quaternion            # the numpy-quaternion package
    >>> import quaternionic
    >>> import sympy
    
    >>> q = Hy(Hy(1, 2), Hy(3, 4))   # the quaternion 1 + 2i + 3j + 4k
    
    >>> q.to_numpy_quaternion()




.. parsed-literal::

    quaternion(1, 2, 3, 4)



.. code:: ipython3

    >>> q.to_quaternionic()




.. parsed-literal::

    quaternionic.array([1., 2., 3., 4.])



.. code:: ipython3

    >>> q.to_sympy()   # exact -- no rounding at all




.. math::

    \displaystyle 1 + 2 i + 3 j + 4 k



Converting Back to a ``Hy``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: ipython3

    >>> Hy.from_numpy_quaternion(np.quaternion(1.5, -0.25, 3, 0.1))




.. parsed-literal::

    Hy(Hy('3/2', '-1/4'), Hy('3', '1/10'))



.. code:: ipython3

    >>> Hy.from_quaternionic(quaternionic.array([1, 0.5, 0, -2]))




.. parsed-literal::

    Hy(Hy('1', '1/2'), Hy('0', '-2'))



.. code:: ipython3

    >>> Hy.from_sympy(sympy.Quaternion(1, sympy.Rational(2, 3), -4, 0))




.. parsed-literal::

    Hy(Hy('1', '2/3'), Hy('-4', '0'))



Each ``from_*`` method converts a *single* quaternion. To convert
several at once – an array of ``numpy-quaternion`` or ``quaternionic``
quaternions, say – convert them one at a time,
e.g. ``[Hy.from_numpy_quaternion(q) for q in array]``.

Floating Point versus Exact Rationals
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``numpy-quaternion`` and ``quaternionic`` both store ``float64``
coordinates, whereas a ``Hy`` is exact. Converting *to* either package
rounds, and converting *back* requires a choice of which rational each
float stands for. Their ``from_*`` methods offer three choices, selected
with two optional keyword arguments:

-  By default, each float becomes the rational with the same shortest
   decimal representation (``0.1`` becomes ``1/10``), just as
   ``Hy.from_array()`` does.
-  ``max_denominator=N`` chooses the closest fraction whose denominator
   is at most ``N``. This recovers simple fractions such as ``1/3``.
-  ``exact=True`` uses the exact binary value of each float. (This
   cannot be combined with ``max_denominator``.)

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

SymPy Round-Trips Exactly
~~~~~~~~~~~~~~~~~~~~~~~~~

SymPy, unlike the other two packages, represents rational numbers
exactly, so ``to_sympy()``/``from_sympy()`` need no rounding decision at
all: a ``Hy -> sympy -> Hy`` round trip always returns the original
value.

.. code:: ipython3

    >>> h = Hy.from_array(['1/3', '-2/7', '5/9', '1/10'])
    
    >>> Hy.from_sympy(h.to_sympy()) == h   # exact, no keywords needed




.. parsed-literal::

    True



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

Interoperability with ``gint`` (Gaussian Integers and Rationals)
----------------------------------------------------------------

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

-  ``Hy.from_gint(x)`` builds a rank 1 ``Hy`` from a ``Zi`` or ``Qi``.
   It needs no ``gint`` import at all.
-  ``some_hy.to_qi()`` returns a ``Qi``, which, as ``Qi`` itself does,
   collapses to a ``Zi`` when both coordinates are integers.
-  ``some_hy.to_zi()`` returns a ``Zi``, or raises ``ValueError`` if a
   coordinate is not an integer.
-  ``some_hy.is_gaussian()`` is a cheap predicate that needs no
   ``gint``.

.. code:: ipython3

    >>> from gint import Zi, Qi
    
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
~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because ``Qi`` and a rank 1 ``Hy`` are two independent implementations
of the same arithmetic, either one is a check on the other. Convert,
operate, and compare:

.. code:: ipython3

    >>> from fractions import Fraction
    
    >>> a, b = Qi('1/2', '-3/5'), Qi(2, '1/3')
    >>> ha, hb = Hy.from_gint(a), Hy.from_gint(b)
    
    >>> print(Hy.from_gint(a * b) == ha * hb)
    >>> print(Hy.from_gint(a / b) == ha / hb)
    >>> print(Fraction(a.norm) == ha.norm_squared())


.. parsed-literal::

    True
    True
    True

Hurwitz Integers (the Class ``Hu``)
-----------------------------------

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

.. code:: ipython3

    >>> from hyprat import Hu

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


Text
~~~~

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

Division with Remainder, GCDs and Associates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

    >>> from hyprat import Hu
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

Primes, Content and Primitive Parts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A nonzero Hurwitz integer is **prime** exactly when its norm is an
ordinary prime number. Notice that ``2`` is *not* prime here: it has
norm 4, and factors as a product of two primes of norm 2. The
**content** of an element is the largest positive integer that divides
it, and it is **primitive** when its content is 1.

.. code:: ipython3

    >>> from hyprat import Hu
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


Big Numbers and the Integer Factoring Helper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Factoring ``q`` means factoring its norm, an ordinary integer. The
module ``hyprat.intfactor`` does that in pure Python (trial division, a
Baillie-PSW primality test and Pollard’s rho method); sympy is used for
very large inputs only if it happens to be installed.

.. code:: ipython3

    >>> from hyprat.intfactor import factorint, is_probable_prime
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
