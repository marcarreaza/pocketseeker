import os
import sys
import pandas as pd
from Bio import PDB
from Bio.PDB import DSSP, ShrakeRupley

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(BASE_DIR)

# Define residue mapping
residue_map = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

# Define atom groups
main_chain_atoms = {"N", "CA", "C", "O"}
polar_atoms = {"O", "N", "OH", "NH"}
non_polar_atoms = {"C"}

# Function to compute SASA
def compute_sasa(pdb_file):
    print("Calculating SASA values")
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]
    sa = ShrakeRupley()
    sa.compute(model, level="A")  # Compute SASA at atom level
    sa.compute(model, level = "R")

    sasa_data = []
    for chain in model.get_chains():
        for residue in chain.get_residues():
            res_name = residue.get_resname()
            res_id = residue.get_id()[1]
            # chain_id = chain.id
            # dssp_key = (chain_id, res_id)
            if res_name in residue_map:
                #dssp_values = dssp[dssp_key]  # Extract DSSP data

                # Initialize SASA values
                total_sasa, main_sasa, side_sasa, polar_sasa, non_polar_sasa = 0, 0, 0, 0, 0
                
                for atom in residue.get_atoms():
                    if atom.element == "H":  # Ignore hydrogen atoms
                        continue

                    sasa = atom.sasa
                    total_sasa += sasa

                    if atom.get_name() in main_chain_atoms:
                        main_sasa += sasa
                    else:
                        side_sasa += sasa

                    if atom.element in polar_atoms:
                        polar_sasa += sasa
                    elif atom.element in non_polar_atoms:
                        non_polar_sasa += sasa

                sasa_data.append([res_id, res_name, residue.sasa, total_sasa, main_sasa, side_sasa, polar_sasa, non_polar_sasa])
    
    df = pd.DataFrame(sasa_data, columns=[
        "Position", "Residue", "SASA", "Total_SASA", "Main-chain_SASA", "Side-chain_SASA", "Polar_SASA", "Non-polar_SASA"
            ])
    df = df.sort_values(by="Position")  # Sort by position
    
    return df



### Function to extract input features
def sasa_feature(file):
    try:
        SASA_data = compute_sasa(file)
        return SASA_data
        
    except Exception as e:
        print(f"Error calculating SASA values: {e}")


### For executing as script
if __name__ == "__main__":
    file = sys.argv[1]
    SASA_data = compute_sasa(file)
    output_csv = os.path.join(BASE_DIR, "SASA.csv")
    SASA_data.to_csv(output_csv, index=False)
    print(f"SASA data saved in {output_csv.split('/')[-1]}")