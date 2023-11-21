#!/bin/bash

# Check if 'pre-commit' package is installed in the current conda environment
if conda list | grep -q "^pre-commit[[:space:]]"; then
    echo "'pre-commit' package is already installed."
else
    # Install 'pre-commit' package using conda
    conda install -c conda-forge -y pre-commit
    echo "'pre-commit' package has been installed."
fi

# Check pre-commit version
pre-commit --version
# Install the git hook scripts
pre-commit install
# Run pre-commit against all files
pre-commit run --all-files
