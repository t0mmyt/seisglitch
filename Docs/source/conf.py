#!/usr/bin/env python

# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html




# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
import os
import sys
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('./../../')) # needed to show docstrings via `automodules`

#import mock
#MOCK_MODULES = ['numpy', 'scipy', 'matplotlib', 'mpl_toolkits', 'obspy', 'pandas', 'yaml']
#for mod_name in MOCK_MODULES:
#    sys.modules[mod_name] = mock.Mock()




# -- Project information -----------------------------------------------------

project   = 'seisglitch'
copyright = '2020, John-Robert Scholz'
author    = 'John-Robert Scholz'

# The short X.Y version.
release = version = ''  # Is set by calling `setup.py docs`

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
try:
    from seisglitch import __version__ as version
except ImportError:
    pass
else:
    release = version




# -- General configuration ---------------------------------------------------

# Read the Docs will set master doc to index instead (or whatever it is you have 
# specified in your settings). Try adding this to your conf.py:
master_doc = 'index'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.coverage', 'sphinx.ext.napoleon']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []




# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'  # pyramid, haiku, alabaster, scrolls, classic, nature


# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
    'logo'                  : 'glitch.png',
    'logo_name'             : False,

    'show_relbars'          : True,
    'fixed_sidebar'         : True,
    'sidebar_width'         : '300px',
    'page_width'            : '1200px',
    'sidebar_includehidden' : True
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']