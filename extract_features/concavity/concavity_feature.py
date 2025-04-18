import os 
import sys
import numpy as np
from Bio import PDB
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def open_pdb(pdb_file):
    """ Abre el archivo PDB y devuelve la estructura """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    return structure

def calculate_curvature(residue1, residue2, residue3):
    """ Calcula la curvatura local entre tres residuos consecutivos """
    # Obtener las coordenadas de Cα
    coord1 = residue1['CA'].get_coord()
    coord2 = residue2['CA'].get_coord()
    coord3 = residue3['CA'].get_coord()

    # Calcular vectores
    vec1 = coord2 - coord1
    vec2 = coord3 - coord2

    # Normalizar los vectores
    vec1 = vec1 / np.linalg.norm(vec1)
    vec2 = vec2 / np.linalg.norm(vec2)

    # Calcular el ángulo entre los dos vectores
    dot_product = np.dot(vec1, vec2)
    angle = np.arccos(dot_product)

    # La curvatura puede definirse como el ángulo dividido por la distancia entre los residuos
    distance = np.linalg.norm(coord3 - coord1)
    curvature = angle / distance
    return curvature

def save_concavity_to_csv(concavity_list, output_csv):
    """ Guarda los valores de concavidad en un archivo CSV único. """
    output_path = os.path.join(os.getcwd(), output_csv)
    df = pd.DataFrame(concavity_list, columns=["File", "Res", "Concavity"])
    df.to_csv(output_path, index=False)
    print(f"Concavity values saved in {output_csv.split('/')[-1]}")

def concavity_extraction(pdb_file):
    """ Extrae la concavidad de cada residuo en la estructura """
    print('Calculating concavity')
    structure = open_pdb(pdb_file)
    model = structure[0]
    concavity = []

    for chain in model:
        residues = [res for res in chain if PDB.is_aa(res, standard=True)]  # Filtramos solo residuos estándar
        if len(residues) < 3:
            continue  # Si hay menos de 3 residuos en la cadena, no podemos calcular la curvatura
        
        # Agregar el primer residuo con "-"
        concavity.append((pdb_file, residues[0].get_id()[1], "-"))

        for i in range(1, len(residues) - 1):
            # Calculamos la curvatura entre tres residuos consecutivos
            curvature = calculate_curvature(residues[i - 1], residues[i], residues[i + 1])
            res_id = residues[i].get_id()[1]
            concavity.append((pdb_file, res_id, curvature))

        # Agregar el último residuo con "-"
        concavity.append((pdb_file, residues[-1].get_id()[1], "-"))

    return concavity

def main(pdb_file, output_csv):
    # Check if the protein.pdb file exists
    if os.path.exists(pdb_file):
        concavity_data = concavity_extraction(pdb_file)
        save_concavity_to_csv(concavity_data, output_csv)
    
    else:
        print(f"{file} does not exist")



### Function to extract input features
def concavity_feature (file):
    try:
        concavity_data = concavity_extraction(file)
        return pd.DataFrame(concavity_data, columns=["File", "Res", "Concavity"])
        
    except Exception as e:
        print(f"Error calculating concavity values: {e}")


### For executing as script
if __name__ == "__main__":
    file = sys.argv[1]
    main(file, os.path.join(BASE_DIR,"concavity.csv"))