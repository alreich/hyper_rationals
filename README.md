# hyper_rationals (`hyprat`)

[![CI](https://github.com/alreich/hyper_rationals/actions/workflows/ci.yml/badge.svg)](https://github.com/alreich/hyper_rationals/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/hyper-rationals/badge/?version=latest)](https://hyper-rationals.readthedocs.io/en/latest/?badge=latest)

Exact, rational-valued **hypercomplex numbers** -- reals, complex
numbers, quaternions, octonions, and beyond -- built via the
[Cayley-Dickson construction](https://en.wikipedia.org/wiki/Cayley%E2%80%93Dickson_construction),
implemented in the `hyprat` package (pure Python, 3.9+, no runtime dependencies).

![hypercomplex numbers](papers/hypercomplex_numbers.png)
<center><sub><sup>(Figure Source: "Hypercomplex Numbers - A Tool for Enhanced Efficiency and Intelligence in Digital Signal Processing" by Zlatka Valkova-Jarvis, et al.)</sup></sub></center>

### Rational Complex Numbers (rank 1)


```python
>>> from hyprat import Hy
>>> from IPython.display import display, Math

>>> z1 = Hy('5/2', '-16/5')
>>> print(f"{z1 = }\n")
>>> print(f"{str(z1) = }\n")
>>> print(f"{complex(z1) = }\n")
>>> display(Math(z1.latex()))
```

    z1 = Hy('5/2', '-16/5')
    
    str(z1) = '(5/2-16/5j)'
    
    complex(z1) = (2.5-3.2j)
    



$\displaystyle \frac{5}{2}-\frac{16}{5}j$



```python
>>> z2 = Hy('1/3', '3/5')
>>> print(f"{str(z1 * z2) = }")
>>> print(f"{str(z2.norm()) = }")
>>> print(f"{str(z2.inverse()) = }")
```

    str(z1 * z2) = '(413/150+13/30j)'
    str(z2.norm()) = '106/225'
    str(z2.inverse()) = '(75/106-135/106j)'


### Rational Quaternions (rank 2)


```python
>>> q1 = Hy(Hy(1.5, '2/3'), Hy('3/7', 4))
>>> print(f"{q1 = }\n")
>>> print(f"{str(q1) = }\n")
>>> display(Math(q1.latex()))
```

    q1 = Hy(Hy('3/2', '2/3'), Hy('3/7', '4'))
    
    str(q1) = '(3/2+2/3i+3/7j+4k)'
    



$\displaystyle \frac{3}{2}+\frac{2}{3}i+\frac{3}{7}j+4k$



```python
>>> Hy.from_array([1.5, '2/3', '3/7', 4]) == q1
```




    True




```python
>>> Hy.parse('(3/2+2/3i+3/7j+4k)') == q1
```




    True




```python
>>> q2 = Hy(z1, z2)
>>> print(f"{str(q1 * q2) = }")
>>> print(f"{str(q2.norm()) = }")
>>> print(f"{str(q2.inverse()) = }")
```

    str(q1 * q2) = '(1403/420-442/105i-407/35j+7871/630k)'
    str(q2.norm()) = '3053/180'
    str(q2.inverse()) = '(450/3053+576/3053i-60/3053j-108/3053k)'


### Rational Octonions (rank 3)


```python
>>> o1 = Hy.from_array(['2/3', 0, 3, -5, '-2/3', '2/5', '-4/5', 2])
>>> print(f"{o1 = }\n")
>>> print(f"{str(o1) = }\n")
>>> display(Math(o1.latex()))
>>> print(f"\n{Hy.parse(str(o1)) == o1 = }")
```

    o1 = Hy(Hy(Hy('2/3', '0'), Hy('3', '-5')), Hy(Hy('-2/3', '2/5'), Hy('-4/5', '2')))
    
    str(o1) = '(2/3+3j-5k-2/3L+2/5iL-4/5jL+2kL)'
    



$\displaystyle \frac{2}{3}+3j-5k-\frac{2}{3}L+\frac{2}{5}iL-\frac{4}{5}jL+2kL$


    
    Hy.parse(str(o1)) == o1 = True



```python
>>> o2 = Hy(q1, q2)
>>> print(f"{str(o1 * o2) = }")
>>> print(f"{str(o2.norm()) = }")
>>> print(f"{str(o2.inverse()) = }")
```

    str(o1 * o2) = '(11407/525+26566/1575i+31091/3150j-1471/150k+1112/105L-157/315iL-5623/630jL-703/42kL)'
    str(o2.norm()) = '158051/4410'
    str(o2.inverse()) = '(6615/158051-2940/158051i-1890/158051j-17640/158051k-11025/158051L+14112/158051iL-1470/158051jL-2646/158051kL)'


### And So On ...

---

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

Basis labels in `str()`, `Hy.parse()`, and `.latex()` are `j` (rank 1),
`i, j, k` (rank 2), `i, j, k, L, iL, jL, kL` (rank 3), and `e1, e2, ...`
from rank 4 (sedenions) up.

## More Features

### Unit elements


```python
>>> print(list(Hy.units(2)))   # the 8 unit quaternions
>>> print(f"{Hy.units(2)['j'].is_unit() = }")
>>> print(f"{q1.is_unit() = }")
```

    ['1', '-1', 'i', '-i', 'j', '-j', 'k', '-k']
    Hy.units(2)['j'].is_unit() = True
    q1.is_unit() = False


### Random values (seedable)


```python
>>> print(Hy.random(2, seed=42))   # a one-off seed for a single call
>>> Hy.seed(2026)                  # or seed the shared default RNG
>>> a = Hy.random(3)
>>> Hy.seed(2026)
>>> print(f"{a == Hy.random(3) = }")
```

    (-6-1/2i-j-k)
    a == Hy.random(3) = True


### Flat arrays and LaTeX options


```python
>>> print(q1.to_array(as_str=True))
>>> display(Math(q1.latex(vinculum="diagonal")))
```

    ['3/2', '2/3', '3/7', '4']



$\displaystyle 3/2+2/3i+3/7j+4k$


### Matrix representation

`to_matrix()` returns the exact (`Fraction`-valued) $2^n \times 2^n$ matrix of
left (or right) multiplication, and `Hy.from_matrix()` is its inverse.
It requires NumPy, which is an optional dependency (see [Installation](#installation)).
For ranks 0-2 (real, complex, quaternion) it is a genuine algebra homomorphism;
from rank 3 on (non-associative) it raises `ValueError` unless you pass
`allow_nonassociative=True`.


```python
>>> qa = Hy.from_array([1, 2, 3, 4])
>>> qb = Hy.from_array([1, '1/2', 0, -1])
>>> print(qa.to_matrix(as_float=True))
>>> print(f"{Hy.from_matrix(qa.to_matrix()) == qa = }")
>>> print(f"{bool((qa.to_matrix() @ qb.to_matrix() == (qa * qb).to_matrix()).all()) = }")
```

    [[ 1. -2. -3. -4.]
     [ 2.  1. -4.  3.]
     [ 3.  4.  1. -2.]
     [ 4. -3.  2.  1.]]
    Hy.from_matrix(qa.to_matrix()) == qa = True
    bool((qa.to_matrix() @ qb.to_matrix() == (qa * qb).to_matrix()).all()) = True


### Interoperability with other quaternion packages

Quaternions (rank 2, or rank-1 complex numbers embedded as `a + bi`) convert to and from
SymPy quaternions (exactly), `numpy-quaternion`, and `quaternionic` via `to_sympy`/`from_sympy`, `to_numpy_quaternion`/`from_numpy_quaternion`,
and `to_quaternionic`/`from_quaternionic`. All three packages are optional and imported lazily.


```python
>>> s = qa.to_sympy()
>>> print(s)
>>> print(f"{Hy.from_sympy(s) == qa = }")
```

    1 + 2*i + 3*j + 4*k
    Hy.from_sympy(s) == qa = True


## Documentation

Full documentation, including a usage guide, the API reference, and a bibliography of
papers and lecture notes on hypercomplex numbers, is on
[Read the Docs](https://hyper-rationals.readthedocs.io/).

Also see the Jupyter notebook
['hyprat_examples.ipynb'](https://github.com/alreich/hyper_rationals/blob/main/notebooks/hyprat_examples.ipynb)
in the `notebooks/` directory.

## Installation

`hyprat` requires Python 3.9 or later and has no runtime dependencies.

```bash
pip install git+https://github.com/alreich/hyper_rationals.git
```

Or, for local development:

```bash
git clone https://github.com/alreich/hyper_rationals.git
cd hyper_rationals
pip install -e ".[dev]"
```

Optional extras (quote them in zsh):

| Extra | Installs | Needed for |
| --- | --- | --- |
| `interop` | `sympy`, `numpy-quaternion`, `quaternionic`, `numpy` | the `to_*`/`from_*` interoperability methods and `to_matrix`/`from_matrix` |
| `test` | `pytest` | running the test suite |
| `docs` | `sphinx`, `sphinx-rtd-theme`, `ipython` | building the docs |
| `dev` | all of the above | development |

For example: `pip install "hyprat[interop] @ git+https://github.com/alreich/hyper_rationals.git"`

## Running the tests

```bash
pytest                                   # unit tests
pytest --doctest-modules src/hyprat      # doctests in the source
```

(or `python -m unittest discover -s tests`)

Tests for the interoperability and matrix methods are skipped unless the optional
packages are installed; use `pip install -e ".[test,interop]"` (or `.[dev]`) to run them.

## Building the docs locally

```bash
pip install -e ".[docs]"
sphinx-build -b html docs/source docs/_build/html
```

## Project layout

```text
hyper_rationals/
+-- src/hyprat/          the hyprat package  (import as `from hyprat import Hy`)
+-- tests/               unit tests (unittest, run via pytest or unittest)
+-- docs/source/         Sphinx documentation source (usage guide, API, bibliography)
+-- notebooks/           Jupyter notebooks (examples, editable sources of the docs, ...)
+-- bibliography/        Freely available papers and books cited in the bibliography
+-- papers/              Additional papers on hypercomplex numbers
+-- .github/workflows/   CI (tests, interoperability tests, docs build)
+-- pyproject.toml       packaging / metadata
+-- .readthedocs.yaml    Read the Docs build config
```

## License

MIT -- see [LICENSE](LICENSE).
