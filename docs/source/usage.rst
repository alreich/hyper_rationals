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

- from a flat list of coefficients (``to_array`` and ``Hy.from_array``)
- from a string representation (``Hy.parse``)
- randomly generated (``Hy.random``)

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

       complex: (9-8/3j)
    quaternion: (-8-3/2i+6/5j+1/2k)
      octonion: (2-3/2i-9/2j+7/5k+3/2L+5/2iL-7/3jL-5/6kL)
      sedenion: (1+1/2e1+4e2-2e3-6e4+4/3e5-e6+4e7-5/6e8-2e9+9e10+3e11+2/3e12+3/2e13+6/5e14+3e15)


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

- ``vinculum`` controls how a non-integer coefficient’s fraction bar is
  typeset: ``"horizontal"`` (the default) uses ``\frac{num}{den}``;
  ``"diagonal"`` uses a plain slash, ``num/den``:
- ``mode`` controls whether the result is wrapped in LaTeX math
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


Interoperability with Other Quaternion Packages
-----------------------------------------------

A rational quaternion (a ``Hy`` of rank 2) can be converted to and from
the quaternion types of three third-party packages. None of them is
required by ``hyprat``; each is imported only when a conversion method
that needs it is called, so install whichever you use:

::

   pip install sympy numpy-quaternion quaternionic

- **SymPy**, ``sympy.algebras.quaternion.Quaternion``:
  ``some_hy.to_sympy()`` and ``Hy.from_sympy(sq)``. SymPy, like ``Hy``,
  represents rational numbers exactly, so this pair round-trips with no
  rounding at all.
- **numpy-quaternion**, imported as ``quaternion``:
  ``some_hy.to_numpy_quaternion()`` and
  ``Hy.from_numpy_quaternion(nq)``. It stores
  ``numpy.quaternion(w, x, y, z)`` as ``float64``.
- **quaternionic**: ``some_hy.to_quaternionic()`` and
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

