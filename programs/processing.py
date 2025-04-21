import glob
import os
import mdtraj as md
from Bio import PDB
import numpy as np
import re

def processing(pdb_file):
    try:
        with open(pdb_file, "r") as f:
            lines = f.readlines()
            if not lines:
                print(f"Warning: {pdb_file} is empty.")
                return

            # Check if already processed
            if lines[0].startswith("HEADER") and "MODEL" in lines[1]:
                print(f"{pdb_file} is already processed. Skipping.")
                return

            # Filter out lines that start with "ATOM"
            atom_lines = [line for line in lines if line.startswith("ATOM")]

            if not atom_lines:
                print(f"No ATOM lines found in {pdb_file}.")
                return

            # Optional: Remove first line (header), add new header
            new_header = ["HEADER    Sample PDB file\n", "MODEL        1\n"]
            
            # Write the filtered ATOM lines and new header to the file
            with open(pdb_file, "w") as file:
                file.writelines(new_header + atom_lines + ["ENDMDL\n"])
                print(f"Processed {pdb_file}")

    except Exception as e:
        print(f"Failed to process {pdb_file}: {e}")



def ConvertMol2toPDB(mol_file):
    pdb_file = f"{os.path.splitext(mol_file)[0]}.pdb"
    try:
        traj = md.load_mol2(mol_file)
        traj.save_pdb(pdb_file)
        print(f"Converted {mol_file} → {pdb_file}")
        processing(pdb_file)
    except Exception as e:
        print(f"Failed to convert {mol_file}: {e}")


def validate_site_txt(site_txt_path):
    try:
        with open(site_txt_path, 'r') as f:
            for line_number, line in enumerate(f, start=1):
                parts = line.strip().split('\t')
                if len(parts) != 2 or not parts[0].isdigit() or not parts[1].isalpha() or len(parts[1]) != 1:
                    raise ValueError(f"Invalid format in {site_txt_path} at line {line_number}: '{line.strip()}'")
    except Exception as e:
        raise ValueError(f"Error reading {site_txt_path}: {e}")

def processing_dataset(input_folder):
    for folder in os.listdir(input_folder):
        folder_path = os.path.join(input_folder, folder)
        if os.path.isdir(folder_path) and not folder.startswith('.'): 
            print(f"Processing protein: {folder} ############################")
            protein_pdb = os.path.join(folder_path, "protein.pdb")
            site_pdb = os.path.join(folder_path, "site.pdb")
            site_txt = os.path.join(folder_path, "site.txt")

            processed_pdbs = []
            if os.path.exists(protein_pdb):
                processing(protein_pdb)
                processed_pdbs.append("protein.pdb")

            if os.path.exists(site_pdb):
                processing(site_pdb)
                processed_pdbs.append("site.pdb")

            # If both PDBs exist, skip .mol2 processing
            if "protein.pdb" in processed_pdbs and ("site.pdb" in processed_pdbs or os.path.exists(site_txt)):
                continue
            else:
                for file in os.listdir(folder_path):
                    if file.lower() in ['protein.mol2', 'site.mol2']:
                        pdb_name = file.lower().replace('.mol2', '.pdb')
                        pdb_path = os.path.join(folder_path, pdb_name)
                        if os.path.exists(pdb_path):
                            print(f"{pdb_path} already exists. Skipping conversion.")
                            continue
                        mol_file = os.path.join(folder_path, file)
                        ConvertMol2toPDB(mol_file)
                        mol_found = True
               


        missing = []
        if not os.path.exists(protein_pdb):
            missing.append("protein.pdb")
        if not os.path.exists(site_pdb) and not os.path.exists(site_txt):
            missing.append("site.pdb or site.txt")
        if site_txt and os.path.exists(site_txt):
            try:
                validate_site_txt(site_txt)
            except Exception as e:
                print(f"❌ {e}")
        if missing:
            raise FileNotFoundError(f"❌ Missing in {folder}: {', '.join(missing)}")
