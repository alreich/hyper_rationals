``hyprat``: Exact Rational Hypercomplex Numbers
===============================================

Introduction
------------

``hyprat`` is a pure-Python package for working with **exact,
rational-valued hypercomplex numbers**: the real numbers, complex
numbers, quaternions, octonions, sedenions, and beyond, all with
rational coordinates, all built recursively by the `Cayley-Dickson
construction <https://en.wikipedia.org/wiki/Cayley%E2%80%93Dickson_construction>`__.
Every coordinate is a ``fractions.Fraction``, so there is no
floating-point rounding anywhere: results are exactly right, however
deep the tower goes. The package needs nothing beyond the Python
standard library.

The package has two number classes and one small helper module:

-  ``Hy`` represents hypercomplex numbers of every rank, classical or
   split.
-  ``Hu`` represents the Hurwitz integers (the integer points of the
   quaternions) and adds number theory.
-  ``hyprat.intfactor`` provides the integer primality testing and
   factoring that the number theory needs.

This documentation is a guided, hands-on tour, divided into three pages
that follow one another:

1. **This page** introduces ``Hy``: the Cayley-Dickson construction,
   building values of every rank, arithmetic, strings, random values,
   units, LaTeX, and matrices.
2. `Split-Hypercomplex Numbers and
   Interoperability <split_hy_and_interop.html>`__ covers the
   split-hypercomplex numbers and the conversions between ``Hy`` and
   other packages (SymPy, numpy-quaternion, quaternionic, and ``gint``).
3. `Hurwitz Integers and Examples <hurwitz_and_examples.html>`__ covers
   ``Hu`` and ``hyprat.intfactor``, followed by worked examples.

The `API reference <api.html>`__ and the
`bibliography <bibliography.html>`__ complete the documentation.

Overview of the Class ``Hy``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A single immutable class, ``Hy``, represents hypercomplex numbers of
every rank. A ``Hy`` always has exactly two components, ``real`` and
``imag``, and both are ``Hy``\ s of one rank lower, or, at the bottom of
the tower, rational numbers. So a rank-1 ``Hy`` is a rational complex
number, a rank-2 ``Hy`` is a rational quaternion made of two complex
numbers, a rank-3 ``Hy`` is a rational octonion made of two quaternions,
and so on. A ``Hy`` of rank :math:`n` has :math:`2^n` rational
coordinates.

A ``Hy`` supports ``+``, ``-``, ``*``, ``/`` and ``**``, conjugation,
inverses and norms, all defined recursively by the Cayley-Dickson
formulas, and all exact. It can be:

-  built from numbers, strings, flat lists of coefficients, or a
   (seedable) random number generator,
-  printed and parsed back, and rendered as LaTeX,
-  turned into its regular-representation matrix,
-  tested for being a unit or a zero divisor,
-  converted to and from the quaternions of SymPy, numpy-quaternion and
   quaternionic, and the Gaussian integers and rationals of ``gint``.

Classical and Split Hypercomplex Numbers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each doubling step of the Cayley-Dickson construction has a parameter
:math:`\mu`, which can be any nonzero rational number. The default,
:math:`\mu = -1`, gives the **classical** tower of complex numbers,
quaternions, octonions, and so on. A step with :math:`\mu = +1` gives
the **split** version of that step, so ``Hy`` also represents
split-complex numbers, split-quaternions, split-octonions, and more. The
two families differ in important ways:

+-----------------------+-----------------------+-----------------------+
|                       | Classical             | Split                 |
|                       | (:math:`\mu = -1`)    | (:math:`\mu = +1` in  |
|                       |                       | the top step)         |
+=======================+=======================+=======================+
| Rank 1                | complex numbers, with | split-complex         |
|                       | :math:`j^2 = -1`      | numbers, with         |
|                       |                       | :math:`j^2 = +1`      |
+-----------------------+-----------------------+-----------------------+
| Rank 2                | quaternions           | split-quaternions,    |
|                       |                       | which are isomorphic  |
|                       |                       | to the                |
|                       |                       | :math:`2 \times 2`    |
|                       |                       | real matrices         |
+-----------------------+-----------------------+-----------------------+
| Rank 3                | octonions             | split-octonions       |
+-----------------------+-----------------------+-----------------------+
| Quadratic form        | positive definite: a  | indefinite: it can be |
|                       | sum of squares        | negative, or zero for |
|                       |                       | a nonzero value       |
+-----------------------+-----------------------+-----------------------+
| Zero divisors         | none up to rank 3,    | already present at    |
|                       | but they appear at    | rank 1                |
|                       | rank 4 (the           |                       |
|                       | sedenions)            |                       |
+-----------------------+-----------------------+-----------------------+

The classical tower is the subject of this page. The split algebras, and
the *signature* that says which steps are split, are explained in
`Split-Hypercomplex
Numbers <split_hy_and_interop.html#split-hypercomplex-numbers>`__.

Overview of the Class ``Hu``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The class ``Hu`` holds the **Hurwitz integers**: the quaternions
:math:`a + b\,i + c\,j + d\,k` whose four coordinates are either all
integers or all halves of odd integers. They are to the quaternions what
the Gaussian integers are to the complex numbers. They form a ring in
which every element has an integer norm, in which there are exactly 24
units, and in which division with remainder works, so that greatest
common divisors, primes and factorizations exist.

A ``Hu`` stores its four coordinates *doubled*, as plain Python
integers, so its arithmetic involves no fractions at all. ``Hy`` plays
the role of the ambient algebra over the rationals, and the two classes
convert into each other explicitly. On top of the ring operations,
``Hu`` provides left and right division with remainder, greatest common
divisors (with Bezout coefficients), associates, primes, content,
factorization into primes in any order of norms, enumeration of the
elements of a given norm, and Lagrange’s four-square theorem. See
`Hurwitz Integers and Examples <hurwitz_and_examples.html>`__.

Integer Factoring and Primality
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A Hurwitz integer is prime exactly when its norm is an ordinary prime
number, and factoring a Hurwitz integer begins with factoring its norm,
a positive integer that may have dozens of digits. So the number theory
in ``Hu`` rests on integer primality testing and integer factoring.

``hyprat`` brings its own small module for this, written in pure Python:
``hyprat.intfactor``. It provides ``is_probable_prime``, ``factorint``,
``primes_upto`` and ``sqrt_mod_prime``. They can be used on their own,
and ``sympy`` is used for very large inputs only when it is installed
and asked for. See `Integer Primality Testing and Factoring in Pure
Python <hurwitz_and_examples.html#integer-primality-testing-and-factoring-in-pure-python>`__.

Setup
-----

Requirements and Installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``hyprat`` needs Python 3.9 or later and has **no runtime
dependencies**: ``Hy``, ``Hu`` and ``hyprat.intfactor`` use only the
Python standard library. To install the package from GitHub:

.. code:: bash

   pip install git+https://github.com/alreich/hyper_rationals.git

For development, install an editable copy of a clone, together with
everything needed for the tests, the documentation and the
interoperability features:

.. code:: bash

   git clone https://github.com/alreich/hyper_rationals.git
   cd hyper_rationals
   pip install -e ".[dev]"

Optional Packages and the ``interop`` Extra
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Several features connect ``Hy`` to other software. None of that software
is a dependency of ``hyprat``. Each package is imported lazily, by the
one method that needs it, so you only need to install the packages whose
features you use. If a package is missing, the method raises an
``ImportError`` that names the package and the ``pip install`` command
for it.

The optional packages can be installed one by one, or in groups through
the *extras* of ``hyprat``:

+-----------------------+-----------------------+-----------------------+
| Extra                 | Installs              | Needed for            |
+=======================+=======================+=======================+
| ``interop``           | ``sympy``,            | the quaternion        |
|                       | ``numpy-quaternion``, | conversions, and      |
|                       | ``quaternionic``,     | ``to_matrix(          |
|                       | ``numpy``             | )``/``from_matrix()`` |
+-----------------------+-----------------------+-----------------------+
| ``test``              | ``pytest``            | running the test      |
|                       |                       | suite                 |
+-----------------------+-----------------------+-----------------------+
| ``docs``              | ``sphinx``,           | building this         |
|                       | ``sphinx-rtd-theme``, | documentation         |
|                       | ``ipython``           |                       |
+-----------------------+-----------------------+-----------------------+
| ``dev``               | everything above      | development           |
+-----------------------+-----------------------+-----------------------+

For example, to install ``hyprat`` together with the interoperability
packages:

.. code:: bash

   pip install "hyprat[interop] @ git+https://github.com/alreich/hyper_rationals.git"

The table below says which package each optional feature needs.

+-----------------------------------+-----------------------------------+
| Feature                           | Needs                             |
+===================================+===================================+
| ``some_hy.to_sympy()`` and        | ``sympy``                         |
| ``Hy.from_sympy()``               |                                   |
+-----------------------------------+-----------------------------------+
| ``some_hy.to_numpy_quaternion()`` | ``numpy-quaternion`` (imported as |
| and                               | ``quaternion``), which uses       |
| ``Hy.from_numpy_quaternion()``    | ``numpy``                         |
+-----------------------------------+-----------------------------------+
| ``some_hy.to_quaternionic()`` and | ``quaternionic``                  |
| ``Hy.from_quaternionic()``        |                                   |
+-----------------------------------+-----------------------------------+
| ``some_hy.to_matrix()`` and       | ``numpy``                         |
| ``Hy.from_matrix()``              |                                   |
+-----------------------------------+-----------------------------------+
| ``some_hy.to_qi()`` and           | ``gint``                          |
| ``some_hy.to_zi()``               |                                   |
+-----------------------------------+-----------------------------------+
| ``Hy.from_gint()`` and            | nothing extra (you only need      |
| ``some_hy.is_gaussian()``         | ``gint`` to make the ``Zi`` or    |
|                                   | ``Qi`` you pass in)               |
+-----------------------------------+-----------------------------------+
| ``some_hy.latex()``               | nothing, since it returns a       |
|                                   | string; to *typeset* it in a      |
|                                   | Jupyter notebook you use          |
|                                   | ``IPython.display.Math``          |
+-----------------------------------+-----------------------------------+
| ``hyprat.intfact                  | ``sympy`` (never required: the    |
| or.factorint(n, method="sympy")`` | default is pure Python)           |
+-----------------------------------+-----------------------------------+
| the group-theory check in `The    | ``finalg`` (optional: the cell    |
| Hurwitz Units and the Binary      | skips the check when it is        |
| Tetrahedral                       | missing)                          |
| Group <hurwitz_and_ex             |                                   |
| amples.html#the-hurwitz-units-and |                                   |
| -the-binary-tetrahedral-group>`__ |                                   |
+-----------------------------------+-----------------------------------+

The ``gint`` package (Gaussian integers and rationals) is not on PyPI
under that name: ``pip install gint`` installs an unrelated package.
Install it from GitHub instead:

.. code:: bash

   pip install git+https://github.com/alreich/gaussian-integers.git

Imports Used on This Page
~~~~~~~~~~~~~~~~~~~~~~~~~

Each page of this documentation begins with its own cell of imports, so
that the pages can be read in any order. On this page, ``numpy`` and
``sympy`` come from the optional ``interop`` extra, and are used only in
`Matrix Representation <#matrix-representation>`__. ``IPython`` is used
to typeset LaTeX.

.. code:: ipython3

    >>> import json
    >>> import random
    >>> from fractions import Fraction
    >>> from IPython.display import display, Math    # typeset LaTeX output in Jupyter
    >>> import numpy as np      # optional (interop extra): to_matrix() and from_matrix()
    >>> import sympy            # optional (interop extra): one determinant cross-check
    >>> import hyprat
    >>> from hyprat import Hy

To check which version of ``hyprat`` is installed:

.. code:: ipython3

    >>> hyprat.__version__




.. parsed-literal::

    '0.2.0'



What is the Cayley-Dickson construction?
----------------------------------------

Starting from the real numbers :math:`\mathbb{R}`, the Cayley-Dickson
construction repeatedly *doubles* an algebra :math:`A` into a new
algebra :math:`A^2 = A \times A` by defining, for pairs :math:`(a, b)`
and :math:`(c, d)` with :math:`a,b,c,d \in A`:

.. math::


   (a, b) + (c, d) = (a + c,\; b + d)

.. math::


   (a, b)(c, d) = (ac - \overline{d}\,b,\; da + b\,\overline{c})

.. math::


   \overline{(a, b)} = (\overline{a}, -b)

where :math:`\overline{a}` denotes the conjugate of :math:`a` (with
:math:`\overline{a} = a` for real :math:`a`). Applying this once to
:math:`\mathbb{R}` gives :math:`\mathbb{C}`; applying it again gives the
quaternions :math:`\mathbb{H}`; again gives the octonions
:math:`\mathbb{O}`; again gives the sedenions :math:`\mathbb{S}`; and so
on forever, with the dimension doubling (:math:`1, 2, 4, 8, 16, \dots`)
at every step. Each doubling costs the algebra an algebraic property:

+--------------+---+---+---------+---------+--------------------------+
| Algebra      | R | D | Commu   | Assoc   | Division algebra (no     |
|              | a | i | tative? | iative? | zero divisors)?          |
|              | n | m |         |         |                          |
|              | k |   |         |         |                          |
+==============+===+===+=========+=========+==========================+
| Reals        | 0 | 1 | yes     | yes     | yes                      |
| :math:       |   |   |         |         |                          |
| `\mathbb{R}` |   |   |         |         |                          |
+--------------+---+---+---------+---------+--------------------------+
| Complex      | 1 | 2 | yes     | yes     | yes                      |
| :math:       |   |   |         |         |                          |
| `\mathbb{C}` |   |   |         |         |                          |
+--------------+---+---+---------+---------+--------------------------+
| Quaternions  | 2 | 4 | **no**  | yes     | yes                      |
| :math:       |   |   |         |         |                          |
| `\mathbb{H}` |   |   |         |         |                          |
+--------------+---+---+---------+---------+--------------------------+
| Octonions    | 3 | 8 | no      | **no**  | yes                      |
| :math:       |   |   |         | (only   |                          |
| `\mathbb{O}` |   |   |         | *altern |                          |
|              |   |   |         | ative*) |                          |
+--------------+---+---+---------+---------+--------------------------+
| Sedenions    | 4 | 1 | no      | no      | **no**                   |
| :math:       |   | 6 |         |         |                          |
| `\mathbb{S}` |   |   |         |         |                          |
+--------------+---+---+---------+---------+--------------------------+

``hyprat``\ ’s ``Hy`` class implements exactly this tower, using
Python’s exact ``fractions.Fraction`` for every real coordinate – so
every result is *exact*, with no floating-point rounding, no matter how
deep the tower goes.

The formulas above describe the classical construction. ``hyprat``
generalizes the doubling step by a parameter :math:`\mu` that replaces
the :math:`-1` in the product, giving the split-hypercomplex numbers.
That generalization is explained, with the notion of a *signature*, in
`Signatures <split_hy_and_interop.html#signatures>`__.

Hypercomplex Numbers
--------------------

This section is about the class ``Hy``: how to construct its values at
each rank, how to do arithmetic with them, how to convert them to and
from other representations, and how to inspect them. Everything here is
classical (every :math:`\mu = -1`); split values are covered on the next
page.

Representation and Construction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Every ``Hy`` value, of every rank, is therefore described by the
following table, in which a rank-:math:`n` value ``Hy(x, y)`` has
components ``x`` and ``y`` of rank :math:`n-1`:

::

   rank 0  ->  a plain fractions.Fraction                (a rational number)
   rank 1  ->  Hy(real, imag), real/imag : Fraction      (a complex number)
   rank 2  ->  Hy(h1, h2),     h1/h2     : rank-1 Hy     (a quaternion)
   rank 3  ->  Hy(h3, h4),     h3/h4     : rank-2 Hy     (an octonion)
   rank n  ->  Hy(x, y),       x/y       : rank-(n-1) Hy

The constructor automatically promotes (“embeds”) whatever you pass it,
so that this invariant always holds. You will see plenty of examples
below. The dimension of a rank-:math:`n` value is :math:`d = 2^n`; the
table in `What is the Cayley-Dickson
construction? <#what-is-the-cayley-dickson-construction>`__ lists it for
each rank.

Rationals and Complex Numbers (Ranks 0 and 1)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the lowest level of the tower the components are rational numbers,
and a ``Hy`` made up of two rational numbers represents a **rational
complex number**. It has rank 1, as the table in `What is the
Cayley-Dickson
construction? <#what-is-the-cayley-dickson-construction>`__ shows.

A bare number given to ``Hy(...)`` becomes a **rank-1** value (a
“complex” number) with an implicit zero imaginary part – exactly the way
Python’s own ``complex(3)`` is ``(3+0j)``.

``Hy`` accepts ``int``, ``float``, ``Fraction``, ``str`` (including
fraction strings like ``'5/2'`` and decimal strings like ``'3.2'``), and
even Python’s built-in ``complex``, for either argument. Two more ways
to build a value: ``Hy.from_array([...])``, from a flat list of
coefficients, and ``Hy.random(rank)``, which draws one at random – both
covered in more detail as we go.

**Equivalent ways to build a rank-1 ``Hy``**

**One Argument**

``imag`` component defaults to ``0``.

.. code:: ipython3

    >>> print(f"             {Hy(3) = }    # from an int")
    >>> print(f"         {Hy('5/2') = }  # from a fraction string")
    >>> print(f"           {Hy(2.5) = }  # from a float")
    >>> print(f"{Hy(Fraction(5, 2)) = }  # from a Fraction instance")
    >>> print(f"{  Hy( (2.5 - 7j) ) = } # from a single Python complex value")


.. parsed-literal::

                 Hy(3) = Hy('3', '0')    # from an int
             Hy('5/2') = Hy('5/2', '0')  # from a fraction string
               Hy(2.5) = Hy('5/2', '0')  # from a float
    Hy(Fraction(5, 2)) = Hy('5/2', '0')  # from a Fraction instance
      Hy( (2.5 - 7j) ) = Hy('5/2', '-7') # from a single Python complex value


--------------

**Important Note on Printing and Parsing**

For readability and parsability, ``hyprat`` prints and parses according
the following simple rule:

*“Coefficient first, then unit symbol”*

That is, if :math:`n,d \in \mathbb{Z}` and :math:`\mu` is a unit symbol,
then ``hyprat`` considers the string “:math:`n/d\mu`” equivalent to the
expression :math:`(n/d)\mu`.

So, in the example below, the string ‘16/5j’ is parsed as (16/5)j.

.. code:: ipython3

    >>> hy1 = Hy('5/2', '-16/5')
    >>> print(f"{hy1 = }")
    >>> print(f"{str(hy1) = }")
    >>> print(f"{Hy.from_string('(5/2-16/5j)') = }\n")


.. parsed-literal::

    hy1 = Hy('5/2', '-16/5')
    str(hy1) = '(5/2-16/5j)'
    Hy.from_string('(5/2-16/5j)') = Hy('5/2', '-16/5')
    


--------------

**Two Arguments**

Input the two components, ``real`` and ``imag``, using a mix of ``str``,
``int``, ``float``, or ``Fraction``.

.. code:: ipython3

    >>> z = Hy('5/2', '-16/5')
    >>> print(f"{z = }")
    >>> print(f"{z.real = }")
    >>> print(f"{z.imag = }")
    >>> print(f"{z.rank = }")
    >>> print(f"{z.dimension = }")
    >>> print(f"{z == Hy(2.5, '-16/5') = }  # Mixed input (string & float)")
    >>> print(f"{z == Hy.from_array([2.5, '-16/5']) = }  # From an array")


.. parsed-literal::

    z = Hy('5/2', '-16/5')
    z.real = Fraction(5, 2)
    z.imag = Fraction(-16, 5)
    z.rank = 1
    z.dimension = 2
    z == Hy(2.5, '-16/5') = True  # Mixed input (string & float)
    z == Hy.from_array([2.5, '-16/5']) = True  # From an array


Floats are interpreted via their “obvious” decimal value, not the exact
“ugly” binary value that a naive Fraction(float) would give. See below:

.. code:: ipython3

    >>> print(f"{Fraction(3.2).real = }")  # naive Fraction
    >>> print(f"{      Hy(3.2).real = }")  # hyprat instead


.. parsed-literal::

    Fraction(3.2).real = Fraction(3602879701896397, 1125899906842624)
          Hy(3.2).real = Fraction(16, 5)


``hyprat`` plays nice with Python’s built-in complex type, for
construction and comparison:

.. code:: ipython3

    >>> z = Hy('5/2', '-16/5')
    >>> w = (2.5 - 3.2j)
    >>> print(f"{        z  = }")
    >>> print(f"{        w  = }  # a Python complex value")
    >>> print(f"{    z == w = }        # cross-type equality")
    >>> print(f"{complex(z) = }  # convert to a Python complex")


.. parsed-literal::

            z  = Hy('5/2', '-16/5')
            w  = (2.5-3.2j)  # a Python complex value
        z == w = True        # cross-type equality
    complex(z) = (2.5-3.2j)  # convert to a Python complex


Complex numbers, like everything else in ``hyprat``, are combined
exactly. Python numbers can be mixed in on either side:

.. code:: ipython3

    >>> # Ordinary arithmetic, exactly, with no floating-point error:
    >>> a = Hy('1/3', '1/7')
    >>> b = Hy('2/3', '-1/7')
    >>> print("a + b =", a + b)
    >>> print("a - b =", a - b)
    >>> print("a * b =", a * b)
    >>> print("a / b =", a / b)
    
    >>> # Hy also happily mixes with plain Python numbers on either side:
    >>> print(a + 1, 1 - a, 2 * a, a / 2)


.. parsed-literal::

    a + b = (1)
    a - b = (-1/3+2/7j)
    a * b = (107/441+1/21j)
    a / b = (89/205+63/205j)
    (4/3+1/7j) (2/3-1/7j) (2/3+2/7j) (1/6+1/14j)


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


Random values of any rank are drawn with ``Hy.random``, which is covered
below; the seed makes the example reproducible:

.. code:: ipython3

    >>> print(f"{Hy.random(1, seed=0) = }")


.. parsed-literal::

    Hy.random(1, seed=0) = Hy('3/4', '-8/3')


Quaternions (Rank 2)
^^^^^^^^^^^^^^^^^^^^

A **quaternion** is built by pairing up two rank-1 (complex) ``Hy``
values: ``Hy(h1, h2)`` where ``h1, h2`` are each rank-1. The result has
rank 2 and 4 real coordinates: :math:`a + bi + cj + dk`.

Quaternion multiplication is **associative but not commutative** – the
famous relations :math:`ij = k`, :math:`ji = -k`, etc.

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


.. code:: ipython3

    >>> one = Hy(1)
    >>> i = Hy(Hy(0, 1), Hy(0, 0))
    >>> j = Hy(Hy(0, 0), Hy(1, 0))
    >>> k = Hy(Hy(0, 0), Hy(0, 1))
    
    >>> for name, val in [("1", one), ("i", i), ("j", j), ("k", k)]:
    ...     print(f"{name} = {val}   (rank {val.rank}, dim {val.dimension})")


.. parsed-literal::

    1 = (1)   (rank 1, dim 2)
    i = (i)   (rank 2, dim 4)
    j = (j)   (rank 2, dim 4)
    k = (k)   (rank 2, dim 4)


.. code:: ipython3

    >>> # The defining quaternion multiplication table:
    >>> print("i*i =", i * i, "   j*j =", j * j, "   k*k =", k * k)
    >>> print("i*j =", i * j, "   j*i =", j * i, "   <-- non-commutative!")
    >>> print("j*k =", j * k, "   k*j =", k * j)
    >>> print("k*i =", k * i, "   i*k =", i * k)
    
    >>> assert i * j == k and j * i == -k
    >>> assert (i * j) * k == i * (j * k)  # associativity DOES hold for quaternions
    >>> print("\nAssociativity check ((ij)k == i(jk)):", (i * j) * k == i * (j * k))


.. parsed-literal::

    i*i = (-1)    j*j = (-1)    k*k = (-1)
    i*j = (k)    j*i = (-k)    <-- non-commutative!
    j*k = (i)    k*j = (-i)
    k*i = (j)    i*k = (-j)
    
    Associativity check ((ij)k == i(jk)): True


.. code:: ipython3

    >>> # A general quaternion, its conjugate, norm, and inverse -- built here
    >>> # straight from a flat array of its 4 coefficients:
    >>> q = Hy.from_array([1, 2, 3, 4])     # 1 + 2i + 3j + 4k
    >>> print("q          =", q)
    >>> print("q.to_array() ->", q.to_array(), " (the inverse of from_array())")
    >>> print("conj(q)    =", q.conjugate())
    >>> print("norm(q)    =", q.norm(), " (== 1^2+2^2+3^2+4^2 =", 1+4+9+16, ")")
    >>> print("q^-1       =", q.inverse())
    >>> print("q * q^-1   =", q * q.inverse(), " (should be 1)")


.. parsed-literal::

    q          = (1+2i+3j+4k)
    q.to_array() -> [Fraction(1, 1), Fraction(2, 1), Fraction(3, 1), Fraction(4, 1)]  (the inverse of from_array())
    conj(q)    = (1-2i-3j-4k)
    norm(q)    = 30  (== 1^2+2^2+3^2+4^2 = 30 )
    q^-1       = (1/30-1/15i-1/10j-2/15k)
    q * q^-1   = (1)  (should be 1)


Octonions (Rank 3)
^^^^^^^^^^^^^^^^^^

An **octonion** pairs two quaternions: ``Hy(h3, h4)`` with ``h3, h4``
rank-2. Octonions have 8 real coordinates and lose **associativity** –
but they remain an *alternative* algebra, meaning any sub-algebra
generated by just two elements *is* associative. Concretely:

.. math:: x(xy) = (xx)y \qquad\text{and}\qquad (yx)x = y(xx)

always hold, even though :math:`x(yz) \ne (xy)z` in general.

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


.. code:: ipython3

    >>> h3 = Hy(Hy(1, 0), Hy(0, 0))
    >>> h4 = Hy(Hy(0, 1), Hy(0, 0))
    >>> o = Hy(h3, h4)
    >>> print("o =", o, "  rank =", o.rank, "  dimension =", o.dimension)
    >>> print("o components:", o.components())
    >>> print("o as a flat array:", o.to_array(as_str=True))


.. parsed-literal::

    o = (1+iL)   rank = 3   dimension = 8
    o components: (Fraction(1, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 1), Fraction(0, 1), Fraction(0, 1))
    o as a flat array: ['1', '0', '0', '0', '0', '1', '0', '0']


.. code:: ipython3

    >>> # Octonions are generally non-associative. Rather than hand-writing three
    >>> # example values, let's just draw them at random -- Hy.random(rank, seed=...)
    >>> # gives a reproducible random value of any rank, handy for exactly this
    >>> # kind of throwaway example (and for the fuzz-testing style examples later
    >>> # in this notebook):
    >>> x = Hy.random(3, seed=1)
    >>> y = Hy.random(3, seed=2)
    >>> z = Hy.random(3, seed=3)
    
    >>> lhs, rhs = (x * y) * z, x * (y * z)
    >>> print("x       =", x)
    >>> print("y       =", y)
    >>> print("z       =", z)
    >>> print("(xy)z   =", lhs)
    >>> print("x(yz)   =", rhs)
    >>> print("equal?  ", lhs == rhs, " <-- associativity fails in general for octonions")


.. parsed-literal::

    x       = (-1-7/3i-3/2j+5/4k+3/2L-3/2iL-9/4jL+4/5kL)
    y       = (-8-7/3i-2/3j-3/5L-8/5iL-jL+1/2kL)
    z       = (-2/5+4i+2/5j+k+9L-9/4iL-1/5jL-kL)
    (xy)z   = (740327/9000-1348769/12000i-328771/1125j-86021/3000k+3341861/36000L+3632917/18000iL+199481/3600jL+9449/720kL)
    x(yz)   = (740327/9000-1608551/12000i-1255609/4500j+42033/1000k+3086419/36000L+3432473/18000iL+327359/3600jL-156191/3600kL)
    equal?   False  <-- associativity fails in general for octonions


.. code:: ipython3

    >>> # ...but they ARE alternative: x(xy) == (xx)y and (yx)x == y(xx) always hold.
    >>> print("x(xy) == (xx)y ?", x * (x * y) == (x * x) * y)
    >>> print("(yx)x == y(xx) ?", (y * x) * x == y * (x * x))


.. parsed-literal::

    x(xy) == (xx)y ? True
    (yx)x == y(xx) ? True


Sedenions and Beyond (Rank 4+)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One more doubling gives the **sedenions** (16 real coordinates, rank 4).
Here the algebra loses its *last* nice property: sedenions have **zero
divisors** – nonzero :math:`u, v` with :math:`uv = 0` – so they’re not
even a division algebra anymore. ``hyprat`` doesn’t stop you here: the
same Cayley–Dickson formulas keep working at every rank, they just stop
guaranteeing certain algebraic properties.

.. code:: ipython3

    >>> # Build the 16 sedenion basis units e_0 (=1), e_1, ..., e_15 via from_array
    >>> # -- a flat list of 15 zeros and a single 1 is an easy way to name a basis
    >>> # vector by index:
    >>> basis = [Hy.from_array([1 if k == idx else 0 for k in range(16)]) for idx in range(16)]
    
    >>> # A known example of sedenion zero divisors: (e_a+e_b)(e_c-e_d) == 0 for
    >>> # certain index combinations (the exact indices are convention-dependent,
    >>> # so we search for one rather than hard-coding a published example):
    >>> zero = Hy.from_array([0] * 16)
    >>> found = None
    >>> for a_idx in range(1, 16):
    ...     for b_idx in range(a_idx + 1, 16):
    ...         u = basis[a_idx] + basis[b_idx]
    ...         for c_idx in range(1, 16):
    ...             for d_idx in range(c_idx + 1, 16):
    ...                 if {c_idx, d_idx} == {a_idx, b_idx}:
    ...                     continue
    ...                 v = basis[c_idx] - basis[d_idx]
    ...                 if u * v == zero:
    ...                     found = (u, v)
    ...                     break
    ...             if found:
    ...                 break
    ...         if found:
    ...             break
    ...     if found:
    ...         break
    
    >>> u, v = found
    >>> print("u =", u)
    >>> print("v =", v)
    >>> print("u * v =", u * v, " <-- zero, even though both u and v are nonzero!")


.. parsed-literal::

    u = (e1+e10)
    v = (e4-e15)
    u * v = (0)  <-- zero, even though both u and v are nonzero!


.. code:: ipython3

    >>> # hyprat doesn't limit you to sedenions -- rank 5 ("pathions", dim 32),
    >>> # rank 6 ("chingons", dim 64), etc. all work identically. Hy.random(rank)
    >>> # is the easiest way to get your hands on one:
    >>> p = Hy.random(5, seed=42)
    >>> print("A random rank-5 value ('pathion'):")
    >>> print("  rank =", p.rank, " dimension =", p.dimension)
    >>> print("  as an array:", p.to_array())


.. parsed-literal::

    A random rank-5 value ('pathion'):
      rank = 5  dimension = 32
      as an array: [Fraction(-6, 1), Fraction(-1, 2), Fraction(-1, 1), Fraction(-1, 1), Fraction(8, 1), Fraction(9, 4), Fraction(-8, 1), Fraction(-7, 2), Fraction(-2, 5), Fraction(-9, 5), Fraction(-1, 2), Fraction(2, 1), Fraction(-1, 2), Fraction(3, 1), Fraction(-9, 2), Fraction(4, 3), Fraction(-1, 2), Fraction(-1, 1), Fraction(-6, 1), Fraction(3, 1), Fraction(2, 3), Fraction(-1, 1), Fraction(1, 1), Fraction(-3, 2), Fraction(-7, 5), Fraction(0, 1), Fraction(2, 5), Fraction(-1, 2), Fraction(-7, 1), Fraction(-2, 3), Fraction(-7, 2), Fraction(-3, 2)]


Zero Divisors
'''''''''''''

A **zero divisor** is a nonzero value :math:`x` for which some nonzero
:math:`y` has :math:`x\,y = 0`. The classical algebras of rank 0 to 3
(the rationals, complex numbers, quaternions and octonions) are division
algebras, and have none. From rank 4 up (the sedenions and beyond) there
are zero divisors, and in the classical sedenions every zero divisor has
*nonzero* norm, so ``is_null()``, which detects values whose norm is
zero, never sees them. (In the split algebras, which are covered on the
next page, the null values are exactly the zero divisors up to rank 3.
See `Zero Divisors <split_hy_and_interop.html#zero-divisors>`__.)

``is_zero_divisor()`` is the exact test, for any rank and signature. It
checks whether the map ``y -> x * y`` is singular, using exact rational
elimination on that map’s matrix. ``0`` itself is not counted.

.. code:: ipython3

    >>> zd = Hy.from_array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0])   # e1 + e10, a sedenion
    >>> print(f"{zd.is_null() = }    {zd.is_zero_divisor() = }")


.. parsed-literal::

    zd.is_null() = False    zd.is_zero_divisor() = True


``annihilator()`` returns an exact basis of the values it kills. With
the default ``kind="left"`` these are the ``y`` with ``x * y == 0``, and
with ``kind="right"`` the ``y`` with ``y * x == 0``:

.. code:: ipython3

    >>> zd_ann = zd.annihilator()
    >>> print(f"{len(zd_ann) = }    {all(zd * y == 0 for y in zd_ann) = }")


.. parsed-literal::

    len(zd_ann) = 4    all(zd * y == 0 for y in zd_ann) = True


There is no left/right argument for ``is_zero_divisor()`` itself: a
value is a zero divisor on one side exactly when it is on the other, and
exactly when its conjugate is. The two annihilators always have the same
dimension, but they can be different sets of values. For the sedenion
above they happen to coincide:

.. code:: ipython3

    >>> left, right = zd.annihilator("left"), zd.annihilator("right")
    >>> len(left), len(right), left == right




.. parsed-literal::

    (4, 4, True)



Arithmetic and the Core API at a Glance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Arithmetic
~~~~~~~~~~

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


Regardless of rank, every ``Hy`` supports the same core API. Here it is
for a quaternion:

.. code:: ipython3

    >>> q = Hy(Hy(1, 2), Hy(3, 4))    # 1 + 2i + 3j + 4k
    
    >>> print("q.real           ->", q.real)
    >>> print("q.imag           ->", q.imag)
    >>> print("q.rank           ->", q.rank)
    >>> print("q.dimension      ->", q.dimension)
    >>> print("q.components()   ->", q.components())
    >>> print("q.to_array()     ->", q.to_array(), " (like components(), but a plain list)")
    >>> print("q.to_array(as_str=True) ->", q.to_array(as_str=True), " (as strings, e.g. for JSON)")
    >>> print("q.conjugate()    ->", q.conjugate())
    >>> print("q.norm()         ->", q.norm())
    >>> print("q.inverse()      ->", q.inverse())
    >>> print("q.is_zero()      ->", q.is_zero())
    >>> print("abs(q)           ->", abs(q), " (a float, since the exact norm is generally irrational)")
    >>> print("bool(q)          ->", bool(q))
    >>> print("q ** 2           ->", q ** 2)
    >>> print("q ** -1          ->", q ** -1, " (== q.inverse())")
    >>> print("list(q)          ->", list(q), " (iterates real, imag)")
    >>> print("q[0], q[1]       ->", q[0], ",", q[1])


.. parsed-literal::

    q.real           -> (1+2j)
    q.imag           -> (3+4j)
    q.rank           -> 2
    q.dimension      -> 4
    q.components()   -> (Fraction(1, 1), Fraction(2, 1), Fraction(3, 1), Fraction(4, 1))
    q.to_array()     -> [Fraction(1, 1), Fraction(2, 1), Fraction(3, 1), Fraction(4, 1)]  (like components(), but a plain list)
    q.to_array(as_str=True) -> ['1', '2', '3', '4']  (as strings, e.g. for JSON)
    q.conjugate()    -> (1-2i-3j-4k)
    q.norm()         -> 30
    q.inverse()      -> (1/30-1/15i-1/10j-2/15k)
    q.is_zero()      -> False
    abs(q)           -> 5.477225575051661  (a float, since the exact norm is generally irrational)
    bool(q)          -> True
    q ** 2           -> (-28+4i+6j+8k)
    q ** -1          -> (1/30-1/15i-1/10j-2/15k)  (== q.inverse())
    list(q)          -> [Hy('1', '2'), Hy('3', '4')]  (iterates real, imag)
    q[0], q[1]       -> (1+2j) , (3+4j)


Other Ways to Construct Hypercomplex Numbers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are three more ways to build a Hy:

-  from a flat list of coefficients (``to_array`` and ``Hy.from_array``)
-  from a string representation (``Hy.parse``)
-  randomly generated (``Hy.random``)

Examples follow:

To/From a Flat List of Coefficients
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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



To recap, ``to_array()`` flattens a ``Hy`` into a plain Python list of
its ``2**rank`` coordinates (as ``Fraction``\ s, or as strings with
``as_str=True``). ``from_array()`` is the inverse, and happily accepts
``int``, ``float``, ``Fraction``, and plain-fraction strings like
``'5/2'``, freely mixed in the same call. The length must be a power of
2 that is at least 2 (every ``Hy`` has rank :math:`\ge 1`, so there is
no length-1 case). This pair is convenient for interop with things like
``numpy`` arrays, or for text-based serialization:

.. code:: ipython3

    >>> # Numbers, strings, and a mix of both -- all equivalent ways in:
    >>> r1 = Hy.from_array([1, 2, 3, 4])
    >>> r2 = Hy.from_array(['1', '2', '3', '4'])
    >>> r3 = Hy.from_array([1, '2/3', 3.5, Fraction(-1, 4)])
    >>> print(r1, "==", r2, "?", r1 == r2)
    >>> print("r3 =", r3)
    
    >>> # A round trip through JSON, via as_str=True:
    >>> o = Hy.random(3, seed=99)
    >>> payload = json.dumps(o.to_array(as_str=True))
    >>> print("Serialized octonion:", payload)
    >>> restored = Hy.from_array(json.loads(payload))
    >>> print("Restored == original:", restored == o)


.. parsed-literal::

    (1+2i+3j+4k) == (1+2i+3j+4k) ? True
    r3 = (1+2/3i+7/2j-1/4k)
    Serialized octonion: ["3/4", "-3/5", "-2", "-1", "-7/3", "3/5", "8", "3"]
    Restored == original: True


Errors are reported for lengths that are not a power of 2 (at least 2),
and for elements that are composite expressions rather than plain
fractions or decimals (use ``Hy.parse()`` or the ``Hy(str)`` constructor
for those):

.. code:: ipython3

    >>> # Error cases: the length must be a power of 2 (>= 2)...
    >>> for bad in ([], [1], [1, 2, 3]):
    ...     try:
    ...         Hy.from_array(bad)
    ...     except ValueError as e:
    ...         print(f"Hy.from_array({bad!r}) raised ValueError: {e}")
    
    >>> # ...and each element must be a *plain* fraction/decimal, not a composite
    >>> # expression like '1+2j' (use Hy.parse()/Hy(str) for those instead):
    >>> try:
    ...     Hy.from_array(['1+2j', '3'])
    ... except ValueError as e:
    ...     print("Hy.from_array(['1+2j', '3']) raised ValueError:", e)


.. parsed-literal::

    Hy.from_array([]) raised ValueError: from_array() requires a length that is a power of 2 and at least 2; got 0
    Hy.from_array([1]) raised ValueError: from_array() requires a length that is a power of 2 and at least 2; got 1
    Hy.from_array([1, 2, 3]) raised ValueError: from_array() requires a length that is a power of 2 and at least 2; got 3
    Hy.from_array(['1+2j', '3']) raised ValueError: cannot parse '1+2j' as a plain fraction for from_array()


To/From a String Representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The string form of a ``Hy`` can be turned back into the value with
``Hy.parse`` (also available as ``Hy.from_string``):

.. code:: ipython3

    >>> o1_str = str(o1)
    
    >>> o1_copy = Hy.parse(o1_str)
    >>> o1_copy == o1




.. parsed-literal::

    True



The forms that ``str``, ``repr`` and ``Hy.parse`` use are described in
`String Forms and Parsing <#string-forms-and-parsing>`__.

From a Random Number Generator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Hy.random(rank, ...)`` draws a random value of the given rank. The
``seed`` argument makes the draw reproducible, and every call in this
documentation passes one. (The details are in
`Reproducibility <#reproducibility>`__ below.)

.. code:: ipython3

    >>> print(f"   complex: {Hy.random(1, seed=1)}")
    >>> print(f"quaternion: {Hy.random(2, seed=2)}")
    >>> print(f"  octonion: {Hy.random(3, seed=3)}")
    >>> print(f"  sedenion: {Hy.random(4, seed=4)}")


.. parsed-literal::

       complex: (-1-7/3j)
    quaternion: (-8-7/3i-2/3j)
      octonion: (-2/5+4i+2/5j+k+9L-9/4iL-1/5jL-kL)
      sedenion: (-2/3-e1+3/4e2-5e3-7e4+3/5e5-2/5e7+8/3e8-1/2e9-2e10-3e11-1/3e12-3/2e13+2e15)


Reproducibility
'''''''''''''''

``Hy.random(rank)``, where ``rank`` is a positive integer, draws a
random rank-``rank`` value: each of its ``2**rank`` coefficients is an
independent ``Fraction(n, d)``, with ``n`` uniform over ``[lo, hi]``
(default ``[-9, 9]``) and ``d`` uniform over ``[1, dmax]`` (default
``[1, 6]``). The defaults for ``lo``, ``hi``, and ``dmax`` are set to
small values to make examples, demos, and tests easy to read.

There are a few ways to control the reproducibility of random output: a
one-off ``seed=`` for a single call, a shared default generator set with
``Hy.seed(...)``, and your own ``random.Random`` passed as ``rng=``.

.. code:: ipython3

    >>> # A one-off seed=, scoped to just this call:
    >>> a = Hy.random(2, seed=7)
    >>> b = Hy.random(2, seed=7)
    >>> print("Hy.random(2, seed=7) is reproducible:", a == b, " ->", a)
    
    >>> # Hy.seed(...) instead fixes a shared default RNG for everything that
    >>> # follows, so bare Hy.random(rank) calls become reproducible too:
    >>> Hy.seed(2026)
    >>> run1 = [Hy.random(2) for _ in range(3)]
    >>> Hy.seed(2026)
    >>> run2 = [Hy.random(2) for _ in range(3)]
    >>> print("Two Hy.seed(2026)-seeded runs match:", run1 == run2)
    
    >>> # Or bring your own random.Random for full control / sharing across calls:
    >>> my_rng = random.Random(123)
    >>> c1 = Hy.random(2, rng=my_rng)
    >>> c2 = Hy.random(2, rng=my_rng)   # draws the *next* values from the same stream
    >>> print("Two draws from a shared rng (generally different):", c1, c2)
    
    >>> # You can also narrow the coefficient ranges with lo=/hi=/dmax=:
    >>> small = Hy.random(2, lo=-1, hi=1, dmax=2, seed=0)
    >>> print("Small-coefficient quaternion:", small)
    
    >>> # rank must be a positive int -- rank 0 isn't a valid Hy (every Hy has
    >>> # rank >= 1 by construction):
    >>> try:
    ...     Hy.random(0, seed=0)
    ... except ValueError as e:
    ...     print("Hy.random(0, seed=0) raised ValueError:", e)


.. parsed-literal::

    Hy.random(2, seed=7) is reproducible: True  -> (1/2+1/2i-8j+8k)
    Two Hy.seed(2026)-seeded runs match: True
    Two draws from a shared rng (generally different): (-8/3-7/4i-j-2k) (8/5+1/3i-4j-5/3k)
    Small-coefficient quaternion: (-1/2i+1/2j)
    Hy.random(0, seed=0) raised ValueError: rank must be a positive int, got 0


String Forms and Parsing
~~~~~~~~~~~~~~~~~~~~~~~~

``str()`` renders a ``Hy`` the way Python renders ``complex`` for rank 1
(``a+bj``), the customary ``a+bi+cj+dk`` for rank 2,
``i, j, k, L, iL, jL, kL`` labelled basis units for rank 3
(“octonions”), and ``e1, e2, ...`` labelled basis units for rank ≥ 4.
``Hy.parse`` (aliased as ``Hy.from_string``) is the exact inverse.

.. code:: ipython3

    >>> values = [
    ...     Hy('5/2', '-16/5'),
    ...     Hy(Hy(1, 2), Hy(3, 4)),
    ...     Hy(Hy(Hy(1, 0), Hy(0, 0)), Hy(Hy(0, 1), Hy(0, 0))),
    ... ]
    >>> for v in values:
    ...     s = str(v)
    ...     r = repr(v)
    ...     parsed = Hy.parse(s)
    ...     print(f"str  = {s!r:20}  repr = {r:35}  parse(str) == v: {parsed == v}")


.. parsed-literal::

    str  = '(5/2-16/5j)'         repr = Hy('5/2', '-16/5')                   parse(str) == v: True
    str  = '(1+2i+3j+4k)'        repr = Hy(Hy('1', '2'), Hy('3', '4'))       parse(str) == v: True
    str  = '(1+iL)'              repr = Hy(Hy(Hy('1', '0'), Hy('0', '0')), Hy(Hy('0', '1'), Hy('0', '0')))  parse(str) == v: True


Units
~~~~~

Every rank has a set of “unit” values – ``+-1``, plus the ``+-1`` basis
elements like ``j``, or ``i``/``j``/``k``, or ``i, j, k, L, iL, jL, kL``
(rank 3), or ``e1, e2, ...`` (rank ≥ 4) – each with exactly one real
coordinate equal to ``+-1`` and every other coordinate 0.
``Hy.units(rank)`` returns all of them for a given rank, as a dict keyed
by their string form; ``some_hy.is_unit()`` answers the same question
for a single value:

.. code:: ipython3

    >>> print("Hy.units(1)  ->", Hy.units(1))
    >>> print("Hy.units(2)  ->", Hy.units(2))
    >>> print("Hy.units(3) has", len(Hy.units(3)), "entries, e.g.:", {k: Hy.units(3)[k] for k in ("1", "iL", "-kL")})
    
    >>> print()
    >>> print("i.is_unit()      ->", i.is_unit())
    >>> print("(i + j).is_unit()->", (i + j).is_unit())
    >>> print("Hy(0, 0).is_unit() ->", Hy(0, 0).is_unit())


.. parsed-literal::

    Hy.units(1)  -> {'1': Hy('1', '0'), '-1': Hy('-1', '0'), 'j': Hy('0', '1'), '-j': Hy('0', '-1')}
    Hy.units(2)  -> {'1': Hy(Hy('1', '0'), Hy('0', '0')), '-1': Hy(Hy('-1', '0'), Hy('0', '0')), 'i': Hy(Hy('0', '1'), Hy('0', '0')), '-i': Hy(Hy('0', '-1'), Hy('0', '0')), 'j': Hy(Hy('0', '0'), Hy('1', '0')), '-j': Hy(Hy('0', '0'), Hy('-1', '0')), 'k': Hy(Hy('0', '0'), Hy('0', '1')), '-k': Hy(Hy('0', '0'), Hy('0', '-1'))}
    Hy.units(3) has 16 entries, e.g.: {'1': Hy(Hy(Hy('1', '0'), Hy('0', '0')), Hy(Hy('0', '0'), Hy('0', '0'))), 'iL': Hy(Hy(Hy('0', '0'), Hy('0', '0')), Hy(Hy('0', '1'), Hy('0', '0'))), '-kL': Hy(Hy(Hy('0', '0'), Hy('0', '0')), Hy(Hy('0', '0'), Hy('0', '-1')))}
    
    i.is_unit()      -> True
    (i + j).is_unit()-> False
    Hy(0, 0).is_unit() -> False


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



LaTeX Rendering
~~~~~~~~~~~~~~~

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

Basis units render the same way `the API reference <api.html>`__
describes for ``str()``: ``j`` at rank 1, ``i``/``j``/``k`` at rank 2,
and ``i``/``j``/``k``/``L``/``iL``/``jL``/``kL`` at rank 3, all
unchanged, except that the ``e1, e2, ...`` labels used from rank 4 up
are subscripted (rendered as ``e_{1}``, ``e_{2}``, …).

See the `API reference <api.html>`__ for the full reference.

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

    >>> s1 = Hy.random(4, seed=5)  # a random sedenion
    >>> print(f"{s1 = }")
    >>> display(Math(s1.latex()))


.. parsed-literal::

    s1 = Hy(Hy(Hy(Hy('-1/6', '1/3'), Hy('7', '5/2')), Hy(Hy('-4', '-2'), Hy('3', '3/5'))), Hy(Hy(Hy('-6/5', '-2'), Hy('-3/4', '-1/2')), Hy(Hy('3/2', '-7/2'), Hy('5/2', '-5'))))



.. math::

    \displaystyle -\frac{1}{6}+\frac{1}{3}e_{1}+7e_{2}+\frac{5}{2}e_{3}-4e_{4}-2e_{5}+3e_{6}+\frac{3}{5}e_{7}-\frac{6}{5}e_{8}-2e_{9}-\frac{3}{4}e_{10}-\frac{1}{2}e_{11}+\frac{3}{2}e_{12}-\frac{7}{2}e_{13}+\frac{5}{2}e_{14}-5e_{15}


The ``mode`` option wraps the expression in LaTeX math delimiters, for
when you need the delimiters included rather than adding them yourself.
Basis labels such as ``iL`` render unchanged at rank 3, while the
``e1, e2, ...`` labels from rank 4 up are subscripted:

.. code:: ipython3

    >>> z = Hy('5/2', '-16/5')
    >>> print("z.latex(mode='inline')  ->", z.latex(mode='inline'))
    >>> print("z.latex(mode='display') ->", z.latex(mode='display'))
    >>> print("o.latex()               ->", o.latex(), "  (rank 3: unchanged labels)")
    >>> print("s1.latex()              ->", s1.latex()[:60], "...  (rank 4: subscripted labels)")


.. parsed-literal::

    z.latex(mode='inline')  -> $\frac{5}{2}-\frac{16}{5}j$
    z.latex(mode='display') -> \[\frac{5}{2}-\frac{16}{5}j\]
    o.latex()               -> \frac{3}{4}-\frac{3}{5}i-2j-k-\frac{7}{3}L+\frac{3}{5}iL+8jL+3kL   (rank 3: unchanged labels)
    s1.latex()              -> -\frac{1}{6}+\frac{1}{3}e_{1}+7e_{2}+\frac{5}{2}e_{3}-4e_{4} ...  (rank 4: subscripted labels)


Matrix Representation
~~~~~~~~~~~~~~~~~~~~~

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

    >>> z1.to_matrix()




.. parsed-literal::

    array([[Fraction(2, 3), Fraction(-3, 2)],
           [Fraction(3, 2), Fraction(2, 3)]], dtype=object)



.. code:: ipython3

    >>> Hy.from_matrix(z1.to_matrix()) == z1




.. parsed-literal::

    True



For a quaternion the regular representation is a :math:`4 \times 4`
matrix. Its columns are ``q * 1``, ``q * i``, ``q * j`` and ``q * k``:

.. code:: ipython3

    >>> M = q.to_matrix()
    >>> M




.. parsed-literal::

    array([[Fraction(1, 1), Fraction(-2, 1), Fraction(-3, 1),
            Fraction(-4, 1)],
           [Fraction(2, 1), Fraction(1, 1), Fraction(-4, 1), Fraction(3, 1)],
           [Fraction(3, 1), Fraction(4, 1), Fraction(1, 1), Fraction(-2, 1)],
           [Fraction(4, 1), Fraction(-3, 1), Fraction(2, 1), Fraction(1, 1)]],
          dtype=object)



.. code:: ipython3

    >>> Hy.from_matrix(M) == q




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

    >>> det_q1 = sympy.Matrix(q1.to_matrix().tolist()).det()
    >>> det_q1 == sympy.Rational(q1.norm().numerator, q1.norm().denominator) ** 2




.. parsed-literal::

    True



Rank >= 3 (Octonions and Beyond): Not a Homomorphism
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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



Starting at rank 4 (sedenions), the algebra even has *zero divisors*
(see `Zero Divisors <#zero-divisors>`__), so ``to_matrix()`` can be
singular for a nonzero input. For example, the pair ``u = e1 + e10`` and
``v = e4 - e15`` found earlier satisfies ``u * v == 0`` (with this
library’s basis ordering):

.. code:: ipython3

    >>> u * v == 0




.. parsed-literal::

    False



.. code:: ipython3

    >>> Mu = u.to_matrix(allow_nonassociative=True, as_float=True)
    >>> round(np.linalg.det(Mu), 6)   # singular, even though u != 0




.. parsed-literal::

    np.float64(0.0)



See the `bibliography <bibliography.html>`__ for references on faithful,
non-matrix-multiplication representations of octonions (e.g. Zorn vector
matrices), which sidestep this associativity obstruction.

Randomized Testing
------------------

Empirically Verifying Algebraic Laws
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A nice way to build confidence in an implementation like this is to
generate many random values and check that the properties we *expect* to
hold (or fail!) at each rank actually do – this is exactly the strategy
``hyprat``\ ’s own test suite uses, and ``Hy.random()`` makes it easy to
write.

.. code:: ipython3

    >>> def check_norm_is_multiplicative(rank, trials=200, seed=0):
    ...     "|xy| == |x||y| should hold for ranks 0-3 (reals/complex/quaternions/octonions)."
    ...     ok = True
    ...     for t in range(trials):
    ...         x = Hy.random(rank, seed=f"{seed}-{rank}-{t}-x")
    ...         y = Hy.random(rank, seed=f"{seed}-{rank}-{t}-y")
    ...         if (x * y).norm() != x.norm() * y.norm():
    ...             ok = False
    ...             break
    ...     return ok
    
    >>> for rank in (1, 2, 3):
    ...     print(f"rank {rank}: norm multiplicative over 200 random trials? ", check_norm_is_multiplicative(rank))


.. parsed-literal::

    rank 1: norm multiplicative over 200 random trials?  True
    rank 2: norm multiplicative over 200 random trials?  True
    rank 3: norm multiplicative over 200 random trials?  True


.. code:: ipython3

    >>> def find_a_sedenion_zero_divisor_via_random_basis_search(trials=2000, seed=0):
    ...     "Randomly sample sums/differences of sedenion basis units until we hit a zero divisor."
    ...     rng = random.Random(seed)
    ...     zero = Hy.from_array([0] * 16)
    ...     for t in range(trials):
    ...         i1, i2, i3, i4 = (rng.randint(1, 15) for _ in range(4))
    ...         if len({i1, i2, i3, i4}) < 4:
    ...             continue
    ...         u = basis[i1] + basis[i2]
    ...         v = basis[i3] - basis[i4]
    ...         if (u * v) == zero:
    ...             return u, v, t
    ...     return None
    
    >>> result = find_a_sedenion_zero_divisor_via_random_basis_search()
    >>> if result:
    ...     u, v, t = result
    ...     print(f"Found a zero-divisor pair after {t + 1} random trials:")
    ...     print("  u =", u)
    ...     print("  v =", v)
    ...     print("  u * v =", u * v)
    ... else:
    ...     print("No zero divisor found in this many trials (try increasing `trials`).")


.. parsed-literal::

    Found a zero-divisor pair after 17 random trials:
      u = (e1+e15)
      v = (e2-e12)
      u * v = (0)


.. code:: ipython3

    >>> def check_associativity_failure_rate(rank, trials=500, seed=0):
    ...     "How often does random (xy)z land exactly on x(yz)? Should be common for"
    ...     "quaternions (rank 2, always associative) and rare for octonions (rank 3)."
    ...     hits = 0
    ...     for t in range(trials):
    ...         x = Hy.random(rank, seed=f"{seed}-{rank}-{t}-x")
    ...         y = Hy.random(rank, seed=f"{seed}-{rank}-{t}-y")
    ...         z = Hy.random(rank, seed=f"{seed}-{rank}-{t}-z")
    ...         if (x * y) * z == x * (y * z):
    ...             hits += 1
    ...     return hits, trials
    
    >>> for rank, label in [(2, "quaternions"), (3, "octonions")]:
    ...     hits, trials = check_associativity_failure_rate(rank)
    ...     print(f"{label:12} (rank {rank}): associative on {hits}/{trials} random triples")


.. parsed-literal::

    quaternions  (rank 2): associative on 500/500 random triples
    octonions    (rank 3): associative on 0/500 random triples


The next page, `Split-Hypercomplex Numbers and
Interoperability <split_hy_and_interop.html>`__, extends the tower to
the split algebras.
