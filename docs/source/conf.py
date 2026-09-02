"""Sphinx configuration for the hyprat documentation."""

import os
import sys

# Make the `src/` layout package importable without installing it.
sys.path.insert(0, os.path.abspath("../../src"))

# -- Project information -----------------------------------------------

project = "hyprat"
copyright = "2026, Alfred J. Reich"
author = "Alfred J. Reich"

try:
    from hyprat import __version__ as release
except ImportError:  # pragma: no cover
    release = "0.1.0"
version = release

# -- General configuration -----------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

autodoc_member_order = "bysource"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

# -- Options for HTML output ----------------------------------------------

try:
    import sphinx_rtd_theme  # noqa: F401

    html_theme = "sphinx_rtd_theme"
except ImportError:  # pragma: no cover
    html_theme = "alabaster"

html_static_path = []
