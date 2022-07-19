# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

import sphinx_rtd_theme

sys.path.append(os.path.abspath('..'))
sys.path.append(os.path.abspath('../ase'))


# -- Project information -----------------------------------------------------

project = 'koopmans'
copyright = '2020, Edward Linscott, Riccardo De Gennaro, and Nicola Colonna'
author = 'Edward Linscott, Riccardo De Gennaro, and Nicola Colonna'
language = None
with open('../koopmans/__init__.py', 'r') as f:
    # The full version, including alpha/beta/rc tags
    release = f.readlines()[1].split('=')[-1].strip()
version = ''.join([c for c in release if not c.isalpha()])

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['recommonmark', 'sphinx.ext.mathjax', 'sphinx.ext.autosectionlabel',
              'sphinxcontrib.bibtex', 'sphinx_toolbox.collapse', 'sphinx.ext.autodoc',
              'numpydoc', 'sphinx.ext.autosummary']
bibtex_bibfiles = ['refs.bib']
autosectionlabel_prefix_document = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

html_theme = 'sphinx_rtd_theme'
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_logo = "_static/logos/k_white_on_transparent_banner.png"
html_favicon = "_static/logos/k_white_on_grey_tiny.png"
html_theme_options = {'logo_only': True, 'display_version': False}

master_doc = 'index'

# -- Autodoc options ----------------------------------------------------------
autodoc_typehints = 'none'
numpydoc_show_class_members = False
