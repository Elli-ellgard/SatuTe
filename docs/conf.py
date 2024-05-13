# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys

sys.path.insert(
    0, os.path.abspath("..")
)  # This points to the project root where your Python files are.


project = "Satute"
copyright = "2024, CIBIV"
author = "CIBIV"

version = "0.1"
release = "0.1"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add any Sphinx extension module names here, as strings.
# They can be extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    "sphinx.ext.autodoc",  # Document code automatically
    # Add other extensions you might be using like 'sphinx.ext.napoleon' if you use Google style docstrings,
    # 'sphinx.ext.todo', etc.
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

language = "en"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
