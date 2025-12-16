# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
from datetime import datetime
import os
import sys
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'interface-builder'
author = 'A.V.C.'
current_year = datetime.now().year
copyright = f'{current_year}'
release = '0.0.2'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["myst_parser",
            "sphinx.ext.autodoc",
            "sphinx.ext.mathjax"]

# enables ```python fenced blocks
myst_enable_extensions = ["colon_fence",
                        "dollarmath",
                        "amsmath"]

templates_path = ['_templates']
exclude_patterns = ['Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
