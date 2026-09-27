Bibliography
============

A curated reading list on hypercomplex number systems (complex numbers,
quaternions, octonions, sedenions, and beyond) and the Cayley-Dickson
construction that generates them. The list runs roughly from
gentlest to most demanding, and from most recent back to the
19th-century sources where this subject began.

Where a work is freely readable online, its title links to that copy
(usually the published version, an arXiv preprint, or the author's own
page). Works without a link are generally behind a publisher paywall
or are books without a legal free copy; the citation is still given in
full so it can be tracked down through a library.

This same bibliography is also maintained as a Jupyter notebook,
``notebooks/bibliography.ipynb``, in the source repository -- edit
that copy and regenerate this page from it, rather than editing this
file directly.

.. contents:: Sections
   :local:
   :depth: 1

Introductory Lecture Notes
---------------------------

Notes written by professors for specific university courses, giving
step-by-step pedagogical progressions and carefully paced
explanations.

* Keith Conrad, `"Quaternion Algebras"
  <https://kconrad.math.uconn.edu/blurbs/ringtheory/quaternionalg.pdf>`_,
  expository notes, University of Connecticut. A self-contained,
  exercise-free walk from :math:`\mathbb{H}` itself up through
  quaternion algebras over a general field -- an approachable route
  into the subject for anyone comfortable with linear algebra and
  ring theory.

* R. Gardner, `"Supplement: The Cayley-Dickson Construction and
  Nonassociative Algebras"
  <https://faculty.etsu.edu/gardnerr/4127/notes-Herstein/notes-Cayley-Dickson.pdf>`_,
  course supplement for Modern Algebra 2, East Tennessee State
  University. Presents the doubling construction itself
  (:math:`\mathbb{R} \to \mathbb{C} \to \mathbb{H} \to \mathbb{O}`)
  at the level of a first ring-theory course, largely following
  Baez's *The Octonions* below.

* John Voight, `*Quaternion Algebras*
  <https://jvoight.github.io/quat-book.pdf>`_, Graduate Texts in
  Mathematics 288, Springer, 2021. Grew out of graduate courses at
  McGill, Dartmouth, and ICERM/Brown; a full, free textbook that
  starts from the same elementary construction as Conrad's notes and
  develops the arithmetic theory in depth.

Expository / Review Articles
------------------------------

Papers that synthesize existing literature and map the landscape of a
field without getting bogged down in cutting-edge, niche proofs.

* John C. Baez, `"The Octonions"
  <https://arxiv.org/abs/math/0105155>`_, *Bulletin of the American
  Mathematical Society* **39** (2002), 145-205 (errata in **42**
  (2005), 213). The standard modern entry point to the whole subject:
  builds the octonions several different ways and surveys their
  connections to Clifford algebras, Bott periodicity, projective and
  Lorentzian geometry, and the exceptional Lie groups.

* John C. Baez, `"Division Algebras and Quantum Theory"
  <https://arxiv.org/abs/1101.5690>`_, *Foundations of Physics*
  **42** (2012), 819-855. Reviews Hurwitz's classification of the
  normed division algebras and argues that real, complex, and
  quaternionic quantum mechanics are best understood as three faces
  of one structure.

* Richard D. Schafer, `*An Introduction to Nonassociative Algebras*
  <https://math.mit.edu/~hrm/palestine/schafer-nonassociative-algebras.pdf>`_,
  Academic Press, 1966 (Dover reprint, 1995). A short, classic
  monograph surveying alternative and Jordan algebras -- including
  the Cayley-Dickson process and the octonions -- written for
  graduate students meeting nonassociative algebra for the first
  time.

Survey Papers
--------------

High-level overviews of a major topic or open problems, showing how
different theorems and sub-fields connect.

* John C. Baez and John Huerta, `"The Algebra of Grand Unified
  Theories" <https://arxiv.org/abs/0904.1556>`_, *Bulletin of the
  American Mathematical Society* **47** (2010), 483-552. Surveys how
  the division algebras and their associated Lie groups organize the
  grand unified theories of particle physics (SU(5), Spin(10),
  Pati-Salam), connecting representation theory to octonionic
  structure.

* John C. Baez and John Huerta, `"Division Algebras and
  Supersymmetry I" <https://arxiv.org/abs/0909.0551>`_, in
  *Superstrings, Geometry, Topology, and C\*-Algebras*, Proc. Symp.
  Pure Math. **81**, American Mathematical Society, 2010, 65-80.
  Surveys why supersymmetric Yang-Mills theory and the superstring
  work only in spacetime dimensions 3, 4, 6, and 10 -- exactly two
  more than the dimensions 1, 2, 4, 8 of the four normed division
  algebras.

* John H. Conway and Derek A. Smith, *On Quaternions and Octonions:
  Their Geometry, Arithmetic, and Symmetry*, A K Peters, 2003. A
  wide-ranging survey connecting the geometry of quaternions and
  octonions to the classification of finite symmetry groups, Hurwitz
  integral quaternions, and the octonion projective plane; no legal
  free copy is available online.

Tutorials & Pedagogical Introductions
----------------------------------------

Resources focused explicitly on building intuition for beginners or
cross-disciplinary researchers, often from conference proceedings or
community sites rather than journals.

* Cohl Furey, `"Standard Model Physics from an Algebra?"
  <https://arxiv.org/abs/1611.09182>`_, PhD thesis, University of
  Waterloo, 2015. Develops quaternions, octonions, and their tensor
  products from scratch for a physics audience with no prior exposure
  to division algebras, building up to their use in describing
  Standard Model particles.

* `"Introduction to Quaternions"
  <https://mil.ufl.edu/nechyba/www/__eel6667.f2003/course_materials/t3.quaternions/intro_quaternions.pdf>`_,
  course material for EEL6667: Kinematics, Dynamics and Control of
  Robot Manipulators, University of Florida. A short, practical
  tutorial aimed at engineers who need quaternions for representing
  3D rotations rather than for their algebraic structure per se.

Seminal / Foundational Papers
--------------------------------

Original papers where breakthrough theorems or definitions were
introduced.

* William Rowan Hamilton, `"On a New Species of Imaginary Quantities
  Connected with a Theory of Quaternions"
  <https://www.maths.tcd.ie/pub/HistMath/People/Hamilton/Quatern1/Quatern1.pdf>`_,
  *Proceedings of the Royal Irish Academy* **2** (1844), 424-434
  (read November 13, 1843). Hamilton's own first announcement of the
  quaternions, three weeks after the October 16, 1843 walk on which
  he carved :math:`i^2 = j^2 = k^2 = ijk = -1` into Broom Bridge.

* Adolf Hurwitz, "Über die Composition der quadratischen Formen von
  beliebig vielen Variabeln", *Nachrichten von der Gesellschaft der
  Wissenschaften zu Göttingen* (1898), 309-316. Proves that a normed
  composition algebra over the reals can exist only in dimension 1,
  2, 4, or 8 -- the theorem that makes :math:`\mathbb{R}`,
  :math:`\mathbb{C}`, :math:`\mathbb{H}`, and :math:`\mathbb{O}` a
  closed list.

* Max Zorn, "Theorie der alternativen Ringe", *Abhandlungen aus dem
  Mathematischen Seminar der Universität Hamburg* **8** (1930),
  123-147. Introduces alternative rings and the "vector-matrix"
  construction now called Zorn's vector-matrix algebra, giving a
  2x2-matrix model of the split-octonions.

* R. D. Schafer, "On the Algebras Formed by the Cayley-Dickson
  Process", *American Journal of Mathematics* **76** (1954),
  435-445. The paper that names and gives the first general algebraic
  treatment of the doubling process this library is built around,
  generalizing Cayley's and Dickson's original constructions to an
  arbitrary base field and scalar.

* L. E. Dickson, "On Quaternions and Their Generalization and the
  History of the Eight Square Theorem", *Annals of Mathematics*
  **20** (1919), 155-171. Dickson's own generalization of Cayley's
  doubling procedure (hence "Cayley-Dickson"), together with a
  history of the eight-square identity that the octonion norm
  encodes.

Short Research Notes & Expository Overviews
-----------------------------------------------

Concise notes, shorter than a full monograph, that often bridge the
gap between textbooks and active research for undergraduate or
graduate newcomers.

* Guillermo Moreno, `"The Zero Divisors of the Cayley-Dickson
  Algebras over the Real Numbers"
  <https://arxiv.org/abs/q-alg/9710013>`_, *Boletín de la Sociedad
  Matemática Mexicana* **4** (1998), 13-28. The first paper to give
  an algebraic description of the zero divisors that appear once the
  Cayley-Dickson process is carried one step past the octonions, to
  the 16-dimensional sedenions.

* Silvio Reggiani, `"The Geometry of Sedenion Zero Divisors"
  <https://arxiv.org/abs/2411.18881>`_, 2024. A short, recent, and
  very readable follow-up to Moreno: shows that the space of
  normalized sedenion zero-divisor pairs is isometric to the
  exceptional Lie group :math:`G_2`.

* John Cowles and Ruben Gamboa, `"The Cayley-Dickson Construction in
  ACL2" <https://arxiv.org/abs/1705.06822>`_, 2017. A short note
  formalizing the Cayley-Dickson construction (complex numbers,
  quaternions, octonions) in the ACL2 theorem prover -- a nice
  complement to a Python implementation like ``hyprat``, since it
  spells out exactly which algebraic identities each level of the
  construction does and does not satisfy.
