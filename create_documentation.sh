#!/bin/bash
set -e
# Ensure this script is run from the project root directory

# Create a directory for Sphinx documentation if it doesn't already exist
mkdir -p docs

# Navigate into the 'docs' directory to initialize Sphinx
cd docs || exit

# Generate .rst files from the Python code
echo running sphinx-apidoc
poetry run sphinx-apidoc -f -o source ../satute

echo running pandoc
pandoc --from=markdown+smart --standalone --to=rst --wrap=none --output=source/introduction.rst ../README.md

# Ensure README.rst is included in index.rst only once
echo Checking if readme is already included
if ! grep -q "README.rst" index.rst; then
    echo Readme is not included, appending it
    echo ".. include:: introduction.rst" >> index.rst
else
    echo Readme is already included
fi

# Build the HTML documentation within the 'docs' directory
echo building sphinx docs
poetry run sphinx-build -b html source build

# Check if build was successful
if [ $? -eq 0 ]; then
    echo "Documentation generated. Open ./docs/_build/html/index.html in your browser to view it."
else
    echo "Failed to generate documentation."
    exit 1
fi
