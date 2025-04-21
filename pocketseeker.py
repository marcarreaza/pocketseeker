import sys
import os
import argparse
import pandas as pd
import joblib
import numpy as np

from Bio import PDB
import warnings
import requests
import shutil
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
print(BASE_DIR)
from extract_features.model_features import extract_features, check_protein_folders
from programs.run_chimera import run_chimera
from extract_features.model_features import model_features
from model.random_forest import random_forest

warnings.filterwarnings("ignore")
BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def predict_binding_sites(file, output_folder=None):
    df_final = extract_features(file)
    model = joblib.load(os.path.join(BASE_DIR, 'model/random_forest_binding_site_model.joblib'))
    X = df_final.drop(columns=['Res'])
    X.replace('-', np.nan, inplace=True)
    X = pd.get_dummies(X, columns=['SS'])

    for col in model.feature_names_in_:
        if col not in X.columns:
            X[col] = 0
    X = X[model.feature_names_in_]

    y_proba = model.predict_proba(X)[:, 1]
    y_pred = (y_proba >= 0.5).astype(int)
    df_predicted = pd.DataFrame({
        'Res': df_final["Res"],
        'Binding_sites': y_pred,
        'Score': y_proba
    })

    # Fill the GAPs in the binding sites
    binding_residues = df_predicted[df_predicted['Binding_sites'] == 1]["Res"].tolist()
    for i in range(len(binding_residues)):
        j = i + 1
        try:
            if binding_residues[i] + 2 == binding_residues[j]:
                binding_residues.append(binding_residues[i] + 1)
        except IndexError:
            break

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", file)
    last_residue = max([list(chain.get_residues())[-1].get_id()[1] for model in structure for chain in model])

    data = []
    for i in range(1, last_residue):
        match = df_predicted[df_predicted["Res"] == i]
        if not match.empty:
            score = float(match["Score"].values[0])
        else:
            score = 0.0
        label = "Yes" if i in binding_residues else "No"
        data.append({'Res': str(i), 'Binding_sites': label, 'Score': round(score, 4)})

    df_final = pd.DataFrame(data)

    filename = os.path.basename(file)
    file_basename = os.path.splitext(filename)[0]

    # Usar la carpeta directamente si el usuario la especificó
    output_dir = output_folder if output_folder else os.path.join(BASE_DIR, "output", file_basename)
    os.makedirs(output_dir, exist_ok=True)

    df_final.to_csv(os.path.join(output_dir, "binding_sites_predictions.csv"), index=False)

    io = PDB.PDBIO()
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_id = residue.get_id()[1]
                for atom in residue:
                    atom.set_bfactor(0)
                if residue_id in binding_residues:
                    for atom in residue:
                        atom.set_bfactor(50)

    output_pdb = os.path.join(output_dir, f"modified_{file_basename}.pdb")
    io.set_structure(structure)
    io.save(output_pdb)

    print(f"\nResults saved in: {output_dir}")
    return output_pdb


def main():
    
    parser = argparse.ArgumentParser(
        prog="pocketseeker",
        description="Predict protein binding sites from a PDB file using a Random Forest ML model.\n"
                    "You can use local files or download them from the PDB server.\n" \
                    "Also with -train_model argument you can train your own model",
        epilog="Output: A CSV file is generated with the predictions and a PDB with marked B-factors of the residues identified as binding sites.\n"
               "Requires that the structure of this folder be maintained.\n" \
               "The output for -train_model is the features retrieved and the model"
    )
    parser.add_argument("file", nargs="?", help="Local PDB file to analyse (default)")
    parser.add_argument('--d', '--directory', help='Analyse all files located in one local directory.')
    parser.add_argument('--local_many', nargs='+', help='Analyse many local pdb files.')
    parser.add_argument('--online', help='Get a pdb file from the pdb server and analyse that.')
    parser.add_argument('--online_many', nargs='+', help='Get many pdb files from the pdb server and analyse them.')
    parser.add_argument("-o", "--output", default='results', help="Directory to store output features.")
    parser.add_argument('-ch', '--chimera', action='store_true', help='Open Chimera after analysis (single file only).')
    parser.add_argument("-train_model", action="store_true", help="Train a random forest model from extracted features.")
    parser.add_argument("-i", "--input_folder", help="Input directory with protein folders for feature extraction (required if -train_model).")

    args = parser.parse_args()
    print("\n" + "#" * 60)
    print("#" + " " * 18 + "🔬 PocketSeeker! 🔬" + " " * 21 + "#")

    # --- Training Mode ---
    if args.train_model:
        print("\n" + "#" * 60)
        print("#" + " " * 15 + "🧠 Training Mode: Feature Extraction & Model Training  🧠" + " " * 10 + "#")
        print("#" * 60 + "\n")
        if not args.input_folder:
            print("Error: -i (input folder) is required when using -train_model.")
            exit(1)

        output_features_path = os.path.join(args.output, "features")
        os.makedirs(output_features_path, exist_ok=True)

        # Check valid folders
        valid_folders = list(check_protein_folders(args.input_folder))
        if not valid_folders:
            print("No valid protein folders found. Exiting.")
            exit(1)
        final_path = os.path.join(output_features_path, "total_features.csv")
        if os.path.exists(final_path):
            print(f"🚀  Starting model training using : {final_path}")
            df = pd.read_csv(final_path)
            random_forest(df, args.output)
        else:
            model_features(args.input_folder, args.output)
            print(f"🚀  Starting model training using : {final_path}")
            df = pd.read_csv(final_path)
            random_forest(df, args.output)
        
    # --- Prediction Mode ---
    else:
        print("\n" + "#" * 60)
        print("#" + " " * 9 + "Prediction Mode: Binding Site Detection" + " " * 10 + "#")
        print("#" * 60 + "\n")
        
        if args.file:
            output_pdb = predict_binding_sites(args.file, args.output)

        elif args.directory:
            for file in os.listdir(args.directory):
                if file.endswith(".pdb"):
                    predict_binding_sites(os.path.join(args.directory, file), args.output)

        elif args.local_many:
            for file in args.local_many:
                predict_binding_sites(file, args.output)

        elif args.online:
            file = download_pdb(args.online)
            output_pdb = predict_binding_sites(file, args.output)

        elif args.online_many:
            for pdb_id in args.online_many:
                file = download_pdb(pdb_id)
                predict_binding_sites(file, args.output)

        else:
            print("No valid prediction input provided. Use -h for help.")

        if args.chimera and args.file:
            run_chimera(output_pdb)


if __name__ == "__main__":
    try:
        main()
    finally:
        downloads_dir = os.path.join(BASE_DIR, "downloads")
        if os.path.exists(downloads_dir):
            shutil.rmtree(downloads_dir)