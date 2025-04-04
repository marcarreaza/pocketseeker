import os
import numpy as np
from Bio import PDB
import pandas as pd

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
    df = pd.DataFrame(distance_list, columns=["File", "Res", "Distance_to_Core"])
    df.to_csv(output_path, index=False)
    print(f"Distancias al núcleo guardadas en {output_csv}")

def distance_to_core_extraction(pdb_file):
    """ Extrae la distancia de cada residuo al núcleo de la estructura """
    print(f"Procesando archivo PDB: {pdb_file}")
    structure = PDB.PDBParser(QUIET=True).get_structure("protein", pdb_file)
    centroid = calculate_centroid(structure)
    distances = []
    if __name__ == "__main__":
        pdb = pdb_file.split("/")[-2]
    else:
        pdb = pdb_file

    for model in structure:
        for chain in model:
            residues = [res for res in chain if PDB.is_aa(res, standard=True) and 'CA' in res]

            if not residues:
                continue

            for residue in residues:
                res_id = residue.get_id()[1]
                distance = calculate_distance_to_core(residue, centroid)
                distances.append((pdb, res_id, distance))

    print(f"Distancias calculadas para {len(distances)} residuos en {pdb}.")
    return distances

def main(dir, output_csv):
    all_distances = []
    for folder in os.listdir(dir):
        folder_path = os.path.join(dir, folder)
        if os.path.isdir(folder_path):
            # Construct the path to the protein.pdb inside the folder
            pdb_file = os.path.join(folder_path, 'protein.pdb')
            # Check if the protein.pdb file exists
            if os.path.exists(pdb_file):
                distance_data = distance_to_core_extraction(pdb_file)
                all_distances.extend(distance_data)

    save_distances_to_csv(all_distances, output_csv)


### Para extraer los features del training set
if __name__ == "__main__":
    dir = "../../input/"
    output_csv = "distances_to_core.csv"

    main(dir, output_csv)


### Para extraer los features del input
def core_distance(file):
    try:
        # Intentamos parsear el archivo PDB
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", file)

        distance_data = distance_to_core_extraction(file)
        return pd.DataFrame(distance_data, columns=["File", "Res", "Distance_to_Core"])
        
    except Exception as e:
        print(f"{file} is not a pdb file: {e}")

    