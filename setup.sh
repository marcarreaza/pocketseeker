#!/bin/bash

# Nombre del entorno (puedes cambiarlo si quieres)
ENV_NAME="pocketseeker"

# Revisa si el archivo environment.yml existe
if [ ! -f programs/environment.yml ]; then
    echo "'environment.yml' was not found. Make sure you have it in the same directory."
    exit 1
fi

# Crea el entorno Conda
echo "Creating Conda environment called '$ENV_NAME'..."
conda env create -n $ENV_NAME -f environment.yml