#!/bin/bash

# scDGE Snakemake Workflow Setup Script

set -e # Exit immediately if a command fails

echo "-------------------------------------------------------"
echo " Starting Setup for scDGE_workflow"
echo "-------------------------------------------------------"

# Check for Mamba or Conda
if command -v mamba &> /dev/null; then
    PACKAGE_MANAGER="mamba"
    echo "Found Mamba"
elif command -v conda &> /dev/null; then
    PACKAGE_MANAGER="conda"
    echo "Mamba not found. Falling back to Conda (this may be slower)."
else
    echo "Error: Neither 'mamba' nor 'conda' was found."
    echo "Please install Miniconda or Miniforge first."
    exit 1
fi

# Create the Master Environment (snakemake)
ENV_NAME="snakemake"

echo "Creating the environment '$ENV_NAME' using $PACKAGE_MANAGER..."

# Check if environment already exists to avoid errors
if $PACKAGE_MANAGER env list | grep -q "$ENV_NAME"; then
    echo "Environment '$ENV_NAME' already exists. Updating it..."
    $PACKAGE_MANAGER env update -f environment.yaml --name $ENV_NAME
else
    $PACKAGE_MANAGER env create -f environment.yaml --name $ENV_NAME
fi

# Final Instructions
echo "-------------------------------------------------------"
echo "Setup Complete!"
echo ""
echo "To start using the pipeline, run:"
echo "  1. conda activate $ENV_NAME"
echo "  2. snakemake --use-conda --cores 4"
echo ""
echo "Note: The first time you run snakemake, it will take"
echo "longer as it builds the R and Python environments."
echo "-------------------------------------------------------"