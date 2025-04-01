import glob
import os
import mdtraj as md
from Bio import PDB
import numpy as np



def processing(pdb_file):
    with open(pdb_file, "r") as f:
        lines = f.readlines()[1:]  # Remove first line
        new_header = ["HEADER    Sample PDB file\n", "MODEL        1\n"]
    with open(pdb_file, "w") as file:
        file.writelines(new_header + lines)

# Function to convert .mol2 to .pdb
def convert_mol2_to_pdb(mol2_file):
    pdb_file = f"{os.path.splitext(mol2_file)[0]}.pdb"  # Correct path handling
    try:
        traj = md.load_mol2(mol2_file)
        traj.save_pdb(pdb_file)
        # Modify header only for protein.pdb
        if "protein" in pdb_file:
            processing(pdb_file)
        print(f"Converted {mol2_file} to {pdb_file}")
    except Exception as e:
        print(f"Failed to convert {mol2_file}: {e}")

# Convert all site.mol2 and protein.mol2 files
for mol2_file in glob.glob("./input/*/*.mol2"):
    convert_mol2_to_pdb(mol2_file)

'''

def compute_b_factor_from_nma(protein_pdb):
    structure = parsePDB(protein_pdb)  # Use ProDy to parse PDB

    anm = ANM("Protein NMA")
    anm.buildHessian(structure.select("protein"))  # Select only protein atoms
    print("Hessian built successfully.")
    anm.calcModes(n_modes=10)
    print("Modes calculated successfully.")

    # Compute B-factors from normal modes
    try:
        bfactors = calcSqFlucts(anm[:])  # Compute from all modes
        bfactors = (8 * np.pi**2 / 3) * bfactors  # Convert to B-factors
        print(f"Computed bfactors for {pdb_file}")
    except Exception as e:
        print(f"Failed to compute bfactors {protein_pdb}: {e}")
    # Assign computed B-factors
    try:
        structure.setBetas(bfactors)
        print(f"Set for  {pdb_file}")
    except Exception as e:
        print(f"Failed to set {protein_pdb}: {e}")

    # Save updated PDB file
    writePDB(protein_pdb, structure)
    print(f"Saved {protein_pdb} with computed B-factors.")

for pdb_file in glob.glob("./input/*/protein.pdb"):
    compute_b_factor_from_nma(pdb_file)

'''

