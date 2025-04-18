import os
import sys
import numpy as np
from Bio import PDB
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def calculate_distance_to_core(residue, centroid):
    """ Calcula la distancia entre un residuo y el centroide (núcleo) de la proteína """
    coord_residue = residue['CA'].get_coord()
    distance = np.linalg.norm(coord_residue - centroid)
    return distance

def calculate_centroid(structure):
    """ Calcula el centroide (promedio) de la proteína a partir de las coordenadas de Cα """
    total_coord = np.zeros(3)
    count = 0

    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue, standard=True) and 'CA' in residue:
                    total_coord += residue['CA'].get_coord()
                    count += 1

    if count > 0:
        return total_coord / count
    else:
        raise ValueError("No se encontraron residuos con átomos Cα en la estructura.")

def save_distances_to_csv(distance_list, output_csv):
    """ Guarda los valores de distancia en un archivo CSV único. """
    output_path = os.path.join(os.getcwd(), output_csv)
    df = pd.DataFrame(distance_list, columns=["Res", "Distance_to_Core"])
    df.to_csv(output_path, index=False)
    print(f"Core distances saved in {output_csv.split('/')[-1]}")

def distance_to_core_extraction(pdb_file):
    """ Extrae la distancia de cada residuo al núcleo de la estructura """
    print(f"Calculating core distances")
    structure = PDB.PDBParser(QUIET=True).get_structure("protein", pdb_file)
    centroid = calculate_centroid(structure)
    distances = []

    for model in structure:
        for chain in model:
            residues = [res for res in chain if PDB.is_aa(res, standard=True) and 'CA' in res]

            if not residues:
                continue

            for residue in residues:
                res_id = residue.get_id()[1]
                distance = calculate_distance_to_core(residue, centroid)
                distances.append((res_id, distance))

    return distances

def main(pdb_file, output_csv):
    # Check if the protein.pdb file exists
    if os.path.exists(pdb_file):
        distance_data = distance_to_core_extraction(pdb_file)
        save_distances_to_csv(distance_data, output_csv)
    else:
        print(f"{file} does not exist")

    

### Function to extract input features
def core_distance(file):
    try:
        distance_data = distance_to_core_extraction(file)
        return pd.DataFrame(distance_data, columns=["Res", "Distance_to_Core"])
        
    except Exception as e:
        print(f"Error calculating core distances: {e}")


### For executing as script
if __name__ == "__main__":
    file = sys.argv[1]
    main(file, os.path.join(BASE_DIR, "core_distances.csv"))  