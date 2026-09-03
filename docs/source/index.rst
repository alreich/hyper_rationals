hyprat: rational hypercomplex numbers
======================================

``hyprat`` provides an immutable ``Hy`` class for exact, rational-valued
hypercomplex numbers built via the `Cayley-Dickson construction
<https://en.wikipedia.org/wiki/Cayley%E2%80%93Dickson_construction>`_:

.. code-block:: text

    rank 0  ->  a plain fractions.Fraction        (a "real")
    rank 1  ->  Hy(real, imag)                    (a "complex" number)
    rank 2  ->  Hy(h1, h2), h1/h2 rank 1           (a "quaternion")
    rank 3  ->  Hy(h3, h4), h3/h4 rank 2           (an "octonion")
    rank n  ->  Hy(x, y),   x/y  rank (n-1)

All arithmetic (``+ - * /``), conjugation, norms, and inverses are
defined recursively by the standard Cayley-Dickson formulas, so the
same class correctly represents rational complex numbers, quaternions,
octonions, and beyond -- all with exact ``fractions.Fraction``
arithmetic, no floating-point rounding.

Installation
------------

.. code-block:: bash

    pip install git+https://github.com/alreich/hyper_rationals.git

Quickstart
----------

.. code-block:: python

    >>> from hyprat import Hy
    >>> z = Hy('5/2', '-16/5')
    >>> str(z)
    '(5/2-16/5j)'

    >>> q = Hy(Hy(1, 2), Hy(3, 4))   # a rational quaternion
    >>> str(q)
    '(1+2i+3j+4k)'

    >>> i = Hy(Hy(0, 1), Hy(0, 0))
    >>> j = Hy(Hy(0, 0), Hy(1, 0))
    >>> str(i * j)
    '(k)'

    >>> Hy.from_array([1, 2, 3, 4])   # same quaternion, from a flat list
    Hy(Hy('1', '2'), Hy('3', '4'))

.. toctree::
   :maxdepth: 2
   :caption: Contents

   usage
   api

Indices and tables
-------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
