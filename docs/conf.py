# Configuration file for the Sphinx documentation builder.
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Add the source directory to Python path for autodoc
sys.path.insert(0, os.path.abspath("../src"))

# -- Project information -----------------------------------------------------

project = 'noLZSS'
copyright = '2024, Omer Kerner'
author = 'Omer Kerner'

# The full version, including alpha/beta/rc tags
import noLZSS
release = noLZSS.__version__

version = release

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'myst_parser',
    'breathe',
    'exhale'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The suffix(es) of source filenames.
source_suffix = ['.rst', '.md']

# The master toctree document.
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Theme options
html_theme_options = {
    'canonical_url': '',
    'analytics_id': '',
    'logo_only': False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# -- Extension configuration -------------------------------------------------

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__'
}

# Napoleon settings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True

# Intersphinx configuration
# intersphinx_mapping = {
#     'python': ('https://docs.python.org/3/', None),
#     'numpy': ('https://numpy.org/doc/stable/', None),
# }

# MyST parser settings
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

# Breathe configuration
breathe_projects = {
    "noLZSS": "_build/doxygen/xml"
}
breathe_default_project = "noLZSS"

# Exhale configuration
exhale_args = {
    # These arguments are required
    "containmentFolder":     "./api",
    "rootFileName":          "cpp_api.rst",
    "rootFileTitle":         "C++ API Reference",
    "doxygenStripFromPath":  "..",
    # Suggested optional arguments
    "createTreeView":        True,
    # TIP: if using the sphinx-bootstrap-theme, you need
    # "treeViewIsBootstrap": True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin":    """
        INPUT            = ../src/cpp
        EXTRACT_ALL      = YES
        EXTRACT_PRIVATE  = YES
        EXTRACT_STATIC   = YES
        GENERATE_HTML    = NO
        GENERATE_LATEX   = NO
        GENERATE_XML     = YES
        XML_OUTPUT       = xml
        RECURSIVE        = YES
        FILE_PATTERNS    = *.h *.hpp *.cpp *.cc
        EXCLUDE_PATTERNS = */.*
        QUIET            = YES
        JAVADOC_AUTOBRIEF = YES
        OPTIMIZE_OUTPUT_FOR_C = NO
        ENABLE_PREPROCESSING = YES
        MACRO_EXPANSION  = YES
        EXPAND_ONLY_PREDEF = NO
        SEARCH_INCLUDES  = YES
        INCLUDE_PATH     = ../src/cpp
        PREDEFINED       = DOXYGEN_SHOULD_SKIP_THIS
        SKIP_FUNCTION_MACROS = YES
        WARN_IF_UNDOCUMENTED = NO
        WARN_IF_DOC_ERROR = YES
        WARN_NO_PARAMDOC = NO
    """
}

# Treat warnings as errors to catch documentation issues
nitpicky = False  # Disabled for now due to missing external dependencies
nitpick_ignore = [
    ('py:class', 'noLZSS._noLZSS.*'),  # Ignore C++ binding internal classes
    ('py:class', 'pathlib.Path'),     # Standard library class
    ('cpp:identifier', 'size_t'),     # Standard C++ types
    ('cpp:identifier', 'uint64_t'),
    ('cpp:identifier', 'std'),
]

# Enable warnings as errors for CI (can be overridden with SPHINXOPTS=-W)
# warnings_as_errors = True  # Uncomment for strict CI builds