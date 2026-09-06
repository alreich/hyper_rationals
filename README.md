# hyper_rationals

[![CI](https://github.com/alreich/hyper_rationals/actions/workflows/ci.yml/badge.svg)](https://github.com/alreich/hyper_rationals/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/hyper-rationals/badge/?version=latest)](https://hyper-rationals.readthedocs.io/en/latest/?badge=latest)

Exact, rational-valued **hypercomplex numbers** -- reals, complex
numbers, quaternions, octonions, and beyond -- built via the
[Cayley-Dickson construction](https://en.wikipedia.org/wiki/Cayley%E2%80%93Dickson_construction),
implemented in the `hyprat` package.


### Rational Complex Numbers (rank 1)


```python
>>> from hyprat import Hy
>>> from IPython.display import display, Math

>>> z = Hy('5/2', '-16/5')
>>> print(f"{z = }\n")
>>> print(f"{str(z) = }\n")
>>> display(Math(z.latex()))
```

    z = Hy('5/2', '-16/5')
    
    str(z) = '(5/2-16/5j)'
    



$\displaystyle \frac{5}{2}-\frac{16}{5}j$


### Rational Quaternions (rank 2)


```python
>>> quat = Hy(Hy(1.5, '2/3'), Hy('3/7', 4))
>>> print(f"{quat = }\n")
>>> print(f"{str(quat) = }\n")
>>> display(Math(quat.latex()))
```

    quat = Hy(Hy('3/2', '2/3'), Hy('3/7', '4'))
    
    str(quat) = '(3/2+2/3i+3/7j+4k)'
    



$\displaystyle \frac{3}{2}+\frac{2}{3}i+\frac{3}{7}j+4k$



```python
>>> Hy.from_array([1.5, '2/3', '3/7', 4]) == quat
```




    True




```python
>>> Hy.parse('(3/2+2/3i+3/7j+4k)') == quat
```




    True



### Rational Octonions (rank 3)


```python
>>> oct = Hy.from_array(['2/3', 0, 3, -5, '-2/3', '2/5', '-4/5', 2])
>>> print(f"{oct = }\n")
>>> print(f"{str(oct) = }\n")
>>> display(Math(oct.latex()))
>>> print(f"\n{Hy.parse(str(oct)) == oct = }")
```

    oct = Hy(Hy(Hy('2/3', '0'), Hy('3', '-5')), Hy(Hy('-2/3', '2/5'), Hy('-4/5', '2')))
    
    str(oct) = '(2/3+3j-5k-2/3L+2/5iL-4/5jL+2kL)'
    



$\displaystyle \frac{2}{3}+3j-5k-\frac{2}{3}L+\frac{2}{5}iL-\frac{4}{5}jL+2kL$


    
    Hy.parse(str(oct)) == oct = True


### And So On ...

-----------------------

The single immutable `Hy` class represents every rank:

```text
rank 0  ->  a plain fractions.Fraction             (a "real")
rank 1  ->  Hy(real, imag)                         (a "complex")
rank 2  ->  Hy(h1, h2), where h1 & h2 are rank 1   (a "quaternion")
rank 3  ->  Hy(h3, h4), where h3 & h4 are rank 2   (an "octonion")
rank n  ->  Hy(x, y),   where x & y are rank (n-1) ("sedenion", "pathion", ...)
```

`+ - * /`, conjugation, norms, and inverses all follow the standard
recursive Cayley-Dickson formulas, using exact `fractions.Fraction`
arithmetic throughout -- no floating-point rounding.

Full documentation, including the API reference, is on
[Read the Docs](https://hyper-rationals.readthedocs.io/).

## Installation

```bash
pip install git+https://github.com/alreich/hyper_rationals.git
```

Or, for local development:

```bash
git clone https://github.com/alreich/hyper_rationals.git
cd hyper_rationals
pip install -e .[dev]
```

## Running the tests

```bash
pytest
```

(or `python -m unittest discover -s tests`)

## Building the docs locally

```bash
pip install -e .[docs]
sphinx-build -b html docs/source docs/_build/html
```

## Project layout

```text
hyper_rationals/
+-- src/hyprat/          the hyprat package  (import as `from hyprat import Hy`)
+-- tests/                unit tests (unittest, run via pytest or unittest)
+-- docs/source/          Sphinx documentation source
+-- notebooks/            Jupyter notebooks (examples, Claude dialog, ...)
+-- papers/               Papers on hypercomplex numbers (quaternions, octonions, ...)
+-- .github/workflows/    CI (tests + docs build)
+-- pyproject.toml        packaging / metadata
+-- .readthedocs.yaml     Read the Docs build config
```

## License

MIT -- see [LICENSE](LICENSE).

