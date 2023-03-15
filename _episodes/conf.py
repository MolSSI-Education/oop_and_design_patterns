# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/stable/config

# -- Path setup --------------------------------------------------------------

import sys
import os

# Get the absolute path of the current directory
dir_path = os.path.abspath(os.path.dirname(__file__))

# Append the directory path to sys.path
sys.path.append(dir_path)

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

# Incase the project was not installed
import os
import sys
# -- Project information -----------------------------------------------------

project = 'molssi_best_practices'
author = 'The Molecular Sciences Software Institute'

# The short X.Y version
version = ''
# The full version, including alpha/beta/rc tags
release = ''

import lexerpatch
from sphinx.highlighting import lexers
from pygments.lexers.python import PythonLexer

# Associate *.py files with the Python lexer
lexers.setdefault('py', PythonLexer)

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.extlinks',
    'sphinx_design',
    'sphinx_copybutton',
    'myst_parser',
    'sphinx_togglebutton'
]

autosummary_generate = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'default'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pydata_sphinx_theme"


# Project logo option
html_logo = "_static/molssi_main_logo.png"
html_favicon = "_static/molssi_square.png"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
html_theme_options = {
	"github_url": "https://github.com/jchen0506/molssi_doc_theme",
	"twitter_url": "https://twitter.com/MolSSI_NSF",

	"logo": {
      "image_light": "_static/molssi_main_logo.png",
      "image_dark": "_static/molssi_main_logo_inverted_white.png",
      "text": "MolSSI Object Oriented Programming and Design Patterns",
      "molssi_light": "molssi_main_logo.png",
      "molssi_dark": "molssi_main_logo_inverted_white.png",
    },
	"show_toc_level": 2,
	"header_links_before_dropdown": 4,
	"external_links": [
      {"name": "MolSSI", "url": "https://molssi.org"}
  ],

	"secondary_sidebar_items": ["page-toc", "sourcelink"],
    "footer_start": [ "molssi_footer" ],
    "footer_end": [],
    "icon_links":[],
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_css_files = [
  'css/custom.css',
]
# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#

# html_sidebars = {
#   "**": ['globaltoc.html', 'search-field.html']
# }



# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'molssi_doc_themedoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'molssi_doc_theme.tex', 'molssi_doc_theme Documentation',
     'molssi_doc_theme', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'molssi_doc_theme', 'molssi_doc_theme Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'molssi_doc_theme', 'molssi_doc_theme Documentation',
     author, 'molssi_doc_theme', 'A short description of the project (less than one line).',
     'Miscellaneous'),
]


# -- Extension configuration -------------------------------------------------