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
    df = pd.DataFrame(distance_list, columns=["File", "Res", "Distance_to_Core"])
    df.to_csv(output_csv, index=False)
    print(f"Distancias al núcleo guardadas en {output_csv}")

def distance_to_core_extraction(pdb_file):
    """ Extrae la distancia de cada residuo al núcleo de la estructura """
    print(f"Procesando archivo PDB: {pdb_file}")
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
                distances.append((pdb_file, res_id, distance))

    print(f"Distancias calculadas para {len(distances)} residuos en {pdb_file}.")
    return distances

def main(pdb_files, output_csv):
    all_distances = []

    for pdb_file in pdb_files:
        distance_data = distance_to_core_extraction(pdb_file)
        all_distances.extend(distance_data)

    save_distances_to_csv(all_distances, output_csv)

if __name__ == "__main__":
    pdb_files = [
        "../../input/1a2b_1/protein.pdb",
        "../../input/1a2n_1/protein.pdb"
    ]  # Lista de archivos PDB a procesar
    output_csv = "distances_to_core.csv"

    main(pdb_files, output_csv)