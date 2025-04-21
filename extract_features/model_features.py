import pandas as pd
import os
import argparse
import sys
from Bio.PDB import PDBParser, PDBIO

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(BASE_DIR)


from concavity.concavity_feature import concavity_feature
from core_distance.core_distance_feature import core_distance
from physicochemical.physicochemical_features import physicochemical_feature
from PSSM.PSSM_feature import PSSM_feature
from SASA.sasa_features import sasa_feature
from secondary_structure.ss_feature import ss_feature
from solvent.solvent_exposure_features import solvent_feature
from programs.site_identifier import site_identifier
from programs.processing import processing_dataset
from programs.processing import processing
from model.random_forest import random_forest

pd.set_option('future.no_silent_downcasting', True)


 
def extract_features (protein_pdb):
    print(f"\nProcessing {protein_pdb}")
    processing(protein_pdb)
    # Cargar el archivo PDB
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", protein_pdb)

    # Filter by SASA and Rsidual Depth before extract the others features
    sasa = sasa_feature (protein_pdb)
    filtered_sasa = sasa[(sasa[sasa.columns[sasa.columns.get_loc('SASA'):]]==0).all(axis=1)]
    pos1 = list(filtered_sasa['Position'])

    solvent = solvent_feature (protein_pdb)
    filtered_solvent = solvent[(solvent["Residue_Depth"]>4.0)]
    pos2 = list(filtered_solvent['Position'])

    positions = list(set(pos1) & set(pos2))

    # Crear una nueva estructura eliminando los residuos no deseados
    for model in structure:
        for chain in model:
            residues_a_quitar = []
            for residue in chain:
                res_id = residue.get_id()[1]  # El ID del residuo (número)
                if res_id in positions:
                    residues_a_quitar.append(residue)
            # Eliminar los residuos marcados
            for residue in residues_a_quitar:
                chain.detach_child(residue.id)

    # Guardar el nuevo archivo PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save("filter_pdb.pdb")
    processing("filter_pdb.pdb")

    concavity = concavity_feature("filter_pdb.pdb")
    distance_to_core = core_distance("filter_pdb.pdb")
    physicochemical = physicochemical_feature ("filter_pdb.pdb")
    pssm = PSSM_feature ("filter_pdb.pdb")
    ss = ss_feature ("filter_pdb.pdb")
    sasa = sasa_feature("filter_pdb.pdb")
    solvent = solvent_feature ("filter_pdb.pdb")
    df_final = pd.concat([
        concavity,
        distance_to_core["Distance_to_Core"],
        physicochemical.drop(columns=["Position", "Residue"]),
        pssm.drop(columns=["File", "Res"]),
        sasa.drop(columns=["Position", "Residue"]),
        ss.drop(columns=["Res"]),
        solvent.drop(columns=["Position", "Residue"]),
    ], axis=1)

    # Delete the intermediate files 
    files = ['filter_pdb.pdb', 'protein.fasta', 'protein.pssm']
    for file in files:
        if os.path.exists(file):
            os.remove(file)

    return df_final

def display_progress_bar(current, total):
    bar_length = 40  # Length of the progress bar
    progress = (current / total) * bar_length
    bar = "█" * int(progress) + "-" * (bar_length - int(progress))
    sys.stdout.write(f"\r[{bar}] {current}/{total} folders processed\n")
    sys.stdout.flush()


def model_features(input_folder, output_folder):
    all_features = []
    processing_dataset(input_folder)
    folders = [folder for folder in os.listdir(input_folder) if os.path.isdir(os.path.join(input_folder, folder))]
    total_folders = len(folders)
    for i,folder in enumerate(folders, start = 1):
        folder_path = os.path.join(input_folder, folder)
        try:
            if os.path.isdir(folder_path):
                display_progress_bar(i, total_folders)
                # Path to the unbound (protein) PDB file
                protein_pdb = os.path.join(folder_path, 'protein.pdb')
                binding_site_pdb = os.path.join(folder_path, 'site.pdb')
                binding_site_txt = os.path.join(folder_path, 'site.txt')
            
                if os.path.exists(protein_pdb):
                    df_features = extract_features (protein_pdb)
                    if os.path.exists(binding_site_pdb):
                        binding_sites = site_identifier(protein_pdb, binding_site_pdb=binding_site_pdb)
                    elif os.path.exists(binding_site_txt):
                        binding_sites = site_identifier(protein_pdb, binding_site_txt=binding_site_txt)
                    
                    protein_output = os.path.join(output_folder, "features", folder)
                    os.makedirs(protein_output, exist_ok=True)
                    binding_sites.to_csv(os.path.join(protein_output, "binding.csv"), index=False)
                    binding_sites = binding_sites.drop(columns=["Residue"])
                    binding_sites["Binding_Site"] = (
                        binding_sites["Binding_Site"]
                        .replace({"Yes": 1, "No": 0})
                        .astype(int)
                    )


                    df_final = df_features.merge(binding_sites, how='left', left_on='Res', right_on='Position')
                    df_final = df_final.drop(columns='Position')
                    df_final = df_final.drop(columns=["File"])
                    df_final.to_csv(os.path.join(protein_output, "features.csv"), index=False)

                    all_features.append(df_final)
        except Exception as e:
            print(f"Failed to process {folder}: {e}")
    if all_features:
        final_df = pd.concat(all_features, ignore_index=True)
        final_path = os.path.join(output_folder, "features", "total_features.csv")
        final_df.to_csv(final_path, index=False)
        print(f"All features saved to {final_path}")

def check_protein_folders(input_folder):
    for protein_name in os.listdir(input_folder):
        protein_path = os.path.join(input_folder, protein_name)
        if not os.path.isdir(protein_path):
            continue
        required_file = os.path.join(protein_path, "protein.pdb")
        site_files = [
            os.path.join(protein_path, "site.pdb"),
            os.path.join(protein_path, "site.txt")
        ]
        if not os.path.exists(required_file) or not any(os.path.exists(f) for f in site_files):
            print(f"Skipping {protein_name}: missing required files.")
            continue
        yield protein_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process protein folders and extract features.")
    parser.add_argument("-i", "--input_folder", required=True, help="Directory containing protein folders.")
    parser.add_argument("-o", "--output_folder", default='results', help="Directory to store output features.")
    args = parser.parse_args()

    output_features_path = os.path.join(args.output_folder, "features")
    os.makedirs(output_features_path, exist_ok=True)

    # Run feature extraction
    valid_folders = list(check_protein_folders(args.input_folder))
    if not valid_folders:
        print("No valid protein folders found. Exiting.")
        exit(1)

    model_features(valid_folders, args.output_folder)

    # Path to final feature CSV
    final_path = os.path.join(output_features_path, "total_features.csv")

    # If the file exists, train the model
    if os.path.exists(final_path):
        print(f"Running model training using: {final_path}")
        df = pd.read_csv(final_path)
        random_forest(df, args.output_folder)
    else:
        print(f"total_features.csv not found at {final_path}. Skipping training.")

