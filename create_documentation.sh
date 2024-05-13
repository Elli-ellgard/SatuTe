#!/bin/bash
# Ensure this script is run from the project root directory

# Create a directory for Sphinx documentation if it doesn't already exist
mkdir -p docs

# Navigate into the 'docs' directory to initialize Sphinx
cd docs || exit

# Initialize Sphinx if conf.py doesn't exist to avoid overwriting
if [ ! -f "conf.py" ]; then
    sphinx-quickstart --quiet --project "Satute" --author "CIBIV" -v 0.1 --release 0.1 --language en --suffix .rst --master index --makefile --batchfile
fi

# Navigate back to the project root
cd ..

# Generate .rst files from the Python code
sphinx-apidoc -f -o docs ./

pandoc --from=markdown+smart --to=rst --wrap=none --atx-headers --output=docs/README.rst README.md

# Ensure README.rst is included in index.rst only once
if ! grep -q "README.rst" docs/index.rst; then
    echo ".. include:: README.rst" >> docs/index.rst
fi

# Build the HTML documentation within the 'docs' directory
sphinx-build -b html docs docs/_build/html

# Check if build was successful
if [ $? -eq 0 ]; then
    echo "Documentation generated. Open ./docs/_build/html/index.html in your browser to view it."
else
    echo "Failed to generate documentation."
    exit 1
fi
