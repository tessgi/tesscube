# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# Try to import the global configuration from sphinx-astropy
try:
    from sphinx_astropy.conf.v2 import *  # noqa: F403
except ImportError:
    print("sphinx_astropy not found, default configuration will be used.")
import os
import sys

# Add the modules source to the path
sys.path.insert(0, os.path.abspath("../../src"))
# sphinx-action still cant find things, try:
for x in os.walk("../../src"):
    sys.path.insert(0, x[0])
print(sys.path)
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../src/"))
sys.path.insert(0, os.path.abspath("../src/tesscube"))
sys.path.insert(0, os.path.abspath("../.."))
sys.path.insert(0, os.path.abspath("../../src/"))
sys.path.insert(0, os.path.abspath("../../src/tesscube"))

project = "tesscube"
copyright = "2024, TESS Science Support Center"
author = "TESS Science Support Center"
release = "0.1.dev"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "nbsphinx",
    # "sphinx_autodoc_typehints",
    "sphinx.ext.autosummary",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx_automodapi.automodapi",
    "numpydoc",
    "sphinx.ext.intersphinx",
]


templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_sidebars = {
    "**": [],
}

html_theme_options = {
    "external_links": [],
    "github_url": "https://github.com/tessgi/tesscube",
}

html_title = "tesscube"


# -- MyST-parser configuration--------------------------------------------------
myst_heading_anchors = 3

# -- Plot configuration -------------------------------------------------------
plot_rcparams = {
    "axes.labelsize": "large",
    "figure.figsize": (6, 6),
    "figure.subplot.hspace": 0.5,
    "savefig.bbox": "tight",
    "savefig.facecolor": "none",
}
plot_apply_rcparams = True
plot_html_show_source_link = False
plot_formats = ["png", "svg", "pdf"]
# Don't use the default - which includes a numpy and matplotlib import
plot_pre_code = ""

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ["**/.ipynb_checkpoints"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# Raw files we want to copy using the sphinxcontrib-rawfiles extension:
# - CNAME tells GitHub the domain name to use for hosting the docs
# - .nojekyll prevents GitHub from hiding the `_static` dir
rawfiles = ["CNAME", ".nojekyll"]

# Make sure text marked up `like this` will be interpreted as Python objects
default_role = "py:obj"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "astropy": ("https://docs.astropy.org/en/latest/", None),
    "astroquery": ("https://astroquery.readthedocs.io/en/latest/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
}
