# hyper_rationals

[![CI](https://github.com/alreich/hyper_rationals/actions/workflows/ci.yml/badge.svg)](https://github.com/alreich/hyper_rationals/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/hyper-rationals/badge/?version=latest)](https://hyper-rationals.readthedocs.io/en/latest/?badge=latest)

Exact, rational-valued **hypercomplex numbers** -- reals, complex
numbers, quaternions, octonions, and beyond -- built via the
[Cayley-Dickson construction](https://en.wikipedia.org/wiki/Cayley%E2%80%93Dickson_construction),
implemented in the `hyprat` package.

```python
from hyprat import Hy

z = Hy('5/2', '-16/5')          # a rational complex number
str(z)                          # '(5/2-16/5j)'

q = Hy(Hy(1, 2), Hy(3, 4))      # a rational quaternion: 1 + 2i + 3j + 4k
i = Hy(Hy(0, 1), Hy(0, 0))
j = Hy(Hy(0, 0), Hy(1, 0))
i * j                           # -> Hy(Hy('0','0'), Hy('0','1'))   (== k)
```

The single immutable `Hy` class represents every rank:

```text
rank 0  ->  a plain fractions.Fraction        (a "real")
rank 1  ->  Hy(real, imag)                    (a "complex" number)
rank 2  ->  Hy(h1, h2), h1/h2 rank 1          (a "quaternion")
rank 3  ->  Hy(h3, h4), h3/h4 rank 2          (an "octonion")
rank n  ->  Hy(x, y),   x/y  rank (n-1)
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
+-- .github/workflows/    CI (tests + docs build)
+-- pyproject.toml        packaging / metadata
+-- .readthedocs.yaml     Read the Docs build config
```

## License

MIT -- see [LICENSE](LICENSE).
