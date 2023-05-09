# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

import os
import sys
from calendar import month_name
from datetime import date

import pybtex.plugin
import sphinx_rtd_theme
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.style.sorting import BaseSortingStyle
from sphinx_pyproject import SphinxConfig

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

sys.path.insert(0, os.path.abspath('../src'))

# -- Project information -----------------------------------------------------
config = SphinxConfig("../pyproject.toml", globalns=globals())
project = config.name
version = config.version
release = 'v' + version
author = config.author
language = "en"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['recommonmark', 'sphinx.ext.mathjax', 'sphinx.ext.autosectionlabel',
              'sphinxcontrib.bibtex', 'sphinx_toolbox.collapse', 'sphinx.ext.autodoc',
              'numpydoc', 'sphinx.ext.autosummary']
bibtex_bibfiles = ['../src/koopmans/references.bib']
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
autodoc_mock_imports = ['ase']
numpydoc_show_class_members = False

# -- Chronological bibliography style -----------------------------------------


class ChronoSortingStyle(BaseSortingStyle):
    def sort(self, entries):
        def get_date(entry):
            month_lookup = list(month_name)
            year = int(entry.fields['year'])
            month = month_lookup.index(entry.fields.get('month', 'January'))
            return date(year, month, 1)
        return sorted(entries, key=get_date)


class MyChronoStyle(UnsrtStyle):
    def __init__(self, *args, **kwargs):
        kwargs.update(sorting_style=ChronoSortingStyle, abbreviate_names=True)
        return super().__init__(*args, **kwargs)


pybtex.plugin.register_plugin('pybtex.style.formatting', 'chrono', MyChronoStyle)


class MyAbbrevUnsrtStyle(UnsrtStyle):
    def __init__(self, *args, **kwargs):
        kwargs.update(abbreviate_names=True)
        return super().__init__(*args, **kwargs)


pybtex.plugin.register_plugin('pybtex.style.formatting', 'abbrevunsrt', MyAbbrevUnsrtStyle)
