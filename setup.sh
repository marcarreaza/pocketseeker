#!/bin/bash

# Nombre del entorno (puedes cambiarlo si quieres)
ENV_NAME="pocketseeker"

# Detectar sistema operativo
OS_TYPE=$(uname)

if [ "$OS_TYPE" == "Darwin" ]; then
    ENV_FILE="programs/environment_mac.yml"
elif [ "$OS_TYPE" == "Linux" ]; then
    ENV_FILE="programs/environment_linux.yml"
else
    echo "Unsupported OS: $OS_TYPE"
    exit 1
fi

# Verificar si el archivo de entorno existe
if [ ! -f "$ENV_FILE" ]; then
    echo "'$ENV_FILE' was not found. Make sure it exists."
    exit 1
fi

# Crear el entorno Conda
echo "Creating Conda environment called '$ENV_NAME' using '$ENV_FILE'..."
conda env create -n $ENV_NAME -f "$ENV_FILE"

# Descomprimir archivos de la base de datos swissprot
SWISSPROT_DIR="programs/swissprot"

if [ -d "$SWISSPROT_DIR" ]; then
    echo "Uncompressing files in '$SWISSPROT_DIR'..."
    gunzip -v "$SWISSPROT_DIR"/*.gz
else
    echo "Directory '$SWISSPROT_DIR' not found. Skipping decompression."
fi
