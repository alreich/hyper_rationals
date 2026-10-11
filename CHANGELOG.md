# Changelog

All notable changes to `hyprat` are listed here, newest first.  The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.2.0] - 2026-10-11

### Added

- **Hurwitz integers.**  New module `hyprat.hurwitz` with the immutable,
  integer-backed class `Hu` (stores doubled coordinates, so every value is
  an exact Python `int` computation):
  - construction, parsing, `from_hy` / `to_hy`, arithmetic, conjugate,
    norm, trace, and the 24 units (`Hu.units()`);
  - division with remainder (`divmod_left`, `divmod_right`), exact division,
    left/right divisibility and associates, canonical associates;
  - greatest common left/right divisors (`gcld`, `gcrd`) and the extended
    versions (`xgcld`, `xgcrd`);
  - primes, content and primitive part, and Conway-Smith factorization
    (`Hu.factor()` returning a `HuFactorization`);
  - enumeration of elements of a given norm (`of_norm`, `count_of_norm`,
    `primes_of_norm`) and Lagrange's four-square theorem
    (`with_norm`, `four_squares`).
- **Integer factoring helper.**  New module `hyprat.intfactor`
  (`is_probable_prime`, `factorint`, `primes_upto`, `sqrt_mod_prime`) using
  the Baillie-PSW test and Pollard-Brent rho; `sympy` is used only when
  asked for.
- `Hy.is_lipschitz()`, `Hy.is_hurwitz()`, and `Hy.embed(rank)`.
- **Split-hypercomplex numbers.**  Each level of a `Hy` carries a sign `mu`
  (any nonzero rational, default -1): `Hy(real, imag, mu=...)`, `signs=`
  keywords on the factory methods, `norm_squared()`, and exact zero-divisor
  tests (`is_null()`, `is_zero_divisor()`, `annihilator()`).
- **`gint` interoperability.**  `Hy.is_gaussian()`, `Hy.from_gint()`,
  `Hy.to_qi()`, `Hy.to_zi()`, and `Hy.parse(..., unit="i"|"j")`.
- Regular-representation matrices: `Hy.to_matrix()` and `Hy.from_matrix()`.
- Interoperability with SymPy, `numpy-quaternion`, and `quaternionic`.
- Documentation: usage guide sections for split-hypercomplex numbers and
  Hurwitz integers, and a bibliography.

### Changed

- Full type annotations on `Hu`, `intfactor`, and the module helpers in
  `hyprat.hypercomplex`; the `__main__` demo became `_self_test()`.

## [0.1.0]

- Initial release: the immutable `Hy` class (Cayley-Dickson construction
  over the rationals), arithmetic, parsing, random values, LaTeX output.
