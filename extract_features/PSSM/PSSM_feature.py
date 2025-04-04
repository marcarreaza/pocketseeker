import os
import subprocess
import numpy as np
import pandas as pd
from Bio import PDB, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1  # Conversión de código de tres letras a una

# Definir el orden correcto de los aminoácidos según el encabezado
correct_amino_acids = "ARNDCQEGHILKMFPSTWYV"

def extract_sequence_from_pdb(pdb_file):
    """ Extrae la secuencia de aminoácidos de un archivo PDB. """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    seq = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue, standard=True):
                    seq.append(seq1(residue.get_resname()))  # Convierte a código de una letra
    return "".join(seq)

def save_fasta(sequence, output_fasta):
    """ Guarda la secuencia en formato FASTA. """
    if sequence:
        seq_record = SeqRecord(Seq(sequence), id="Protein")
        SeqIO.write(seq_record, output_fasta, "fasta")
    else:
        raise ValueError("La secuencia extraída del PDB está vacía.")

def run_psiblast(fasta_file, output_pssm, db, iterations=3):
    """ Ejecuta PSI-BLAST para generar la matriz PSSM. """
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"El archivo FASTA {fasta_file} no existe.")

    cmd = [
        "psiblast", "-query", fasta_file, "-db", db, "-num_iterations", str(iterations),
        "-out_ascii_pssm", output_pssm
    ]
    
    try:
        subprocess.run(cmd, check=True)
        if not os.path.exists(output_pssm):
            raise FileNotFoundError(f"Error: PSI-BLAST no generó el archivo {output_pssm}")
    except subprocess.CalledProcessError as e:
        print(f"Error ejecutando PSI-BLAST: {e}")
        exit(1)

def parse_pssm_and_calculate_MI_DI(pssm_file, sequence, pdb_file):
    """ Extrae la PSSM de PSI-BLAST y calcula MI y DI, añadiendo el path completo del archivo. """
    with open(pssm_file) as f:
        lines = f.readlines()

    pssm_data = []
    residues = []

    # Procesar las líneas relevantes del archivo PSSM
    for line in lines[3:]:
        tokens = line.strip().split()
        if len(tokens) < 23:  # Verifica que la línea tenga suficientes datos
            continue

        residues.append(tokens[1])  # Extrae el residuo (letra de aminoácido)
        pssm_data.append(list(map(int, tokens[2:22])))  # Extrae las 20 primeras columnas (de la 2 a la 21)

    if len(pssm_data) != len(sequence):
        print(f"Advertencia: El número de filas en la matriz PSSM no coincide con la longitud de la secuencia ({len(sequence)}).")

    pssm_matrix = np.array(pssm_data)

    # Calcular Información Mutua (MI)
    mi_matrix = np.sum(pssm_matrix, axis=1).tolist()

    # Calcular Información Directa (DI)
    di_matrix = [0] + [abs(mi_matrix[i] - mi_matrix[i - 1]) for i in range(1, len(mi_matrix))]

    # Crear el DataFrame con las columnas en el orden correcto
    df = pd.DataFrame(pssm_matrix, columns=list(correct_amino_acids))
    if __name__ == "__main__":
        pdb_file = pdb_file.split("/")[-2]

    df.insert(0, "File", pdb_file)  # Agregar la columna con el path completo del archivo PDB
    df.insert(1, "Res", range(1, len(sequence) + 1))  # Agregar la columna de posición
    df["MI"] = mi_matrix  # Asignar MI a cada residuo
    df["DI"] = di_matrix  # Asignar DI a cada residuo

    return df

def save_df_to_csv(df, output_file, first_file=False):
    """ Guarda el DataFrame en un archivo CSV, sobrescribiéndolo si es el primer archivo procesado. """
    output_path = os.path.join(os.getcwd(), output_file)
    mode = "w" if first_file else "a"  # "w" para el primer archivo (sobrescribe), "a" para el resto (añade)
    header = first_file  # Solo escribe el encabezado en el primer archivo
    df.to_csv(output_path, sep=",", index=False, mode=mode, header=header)
    print(f"Datos guardados en {output_file}")

def process_pdb_file(pdb_file):
    """ Procesa un solo archivo PDB y devuelve su DataFrame de resultados. """
    fasta_file = "protein.fasta"
    pssm_file = "protein.pssm"

    # Extraer secuencia del PDB y guardar en FASTA
    sequence = extract_sequence_from_pdb(pdb_file)
    save_fasta(sequence, fasta_file)
    
    # Ejecutar PSI-BLAST
    print(f"Ejecutando PSI-BLAST para {pdb_file}...")
    
    if __name__ == "__main__":
        db = "../swissprot/swissprot"
    else:
        db = "../extract_features/swissprot/swissprot"

    run_psiblast(fasta_file, pssm_file, db)
    
    # Procesar PSSM y calcular MI y DI
    print(f"Procesando PSSM y calculando MI/DI para {pdb_file}...")
    df = parse_pssm_and_calculate_MI_DI(pssm_file, sequence, pdb_file)

    return df

def main(dir, output_file):
    """ Procesa múltiples archivos PDB y guarda los resultados en un único CSV, sobrescribiendo si ya existe. """
    for i, folder in enumerate(os.listdir(dir)):
        folder_path = os.path.join(dir, folder)
        if os.path.isdir(folder_path):
            # Construct the path to the protein.pdb inside the folder
            pdb_file = os.path.join(folder_path, 'protein.pdb')
            # Check if the protein.pdb file exists
            if os.path.exists(pdb_file):
                df = process_pdb_file(pdb_file)
                save_df_to_csv(df, output_file, first_file=(i == 0))  # El primer archivo sobrescribe el CSV


### Para extraer los features del training set
if __name__ == "__main__":
    dir = "../../input/"
    output_file = "pssm_and_coevolution.csv"

    main(dir, output_file)


### Para extraer los features del input
def PSSM_feature(file):
    try:
        # Intentamos parsear el archivo PDB
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", file)

        return process_pdb_file(file)
        
    except Exception as e:
        print(f"{file} is not a pdb file: {e}")
