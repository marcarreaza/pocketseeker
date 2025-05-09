import os
import sys
import pandas as pd
from Bio import PDB

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def get_physicochemical_feature (pdb_file):
    print(f"Calculating physicochemical data")
    parser = PDB.PDBParser(QUIET=True)
    model = parser.get_structure("protein", pdb_file)[0]

    dict_physicochemical_characteristics = {
        'CYS': [0, 5.05, 103.1, 4.66, 2.5, 5.5],  # Cysteine: no net charge
        'ASP': [-1, 2.77, 115.1, -4.12, -3.5, 13.0],  # Aspartic Acid: 1 negative charge
        'SER': [0, 5.68, 87.1, -2.84, -0.8, 9.2],  # Serine: no net charge
        'GLN': [0, 5.65, 128.1, -2.76, -3.5, 10.5],  # Glutamine: no net charge
        'LYS': [1, 9.74, 128.2, -4.18, -3.9, 11.3],  # Lysine: 1 positive charge
        'ILE': [0, 6.02, 113.2, 5.58, 4.5, 5.2],  # Isoleucine: no net charge
        'PRO': [0, 6.3, 97.1, -3.03, -1.6, 8.0],  # Proline: no net charge
        'THR': [0, 5.66, 101.1, -1.2, -0.7, 8.6],  # Threonine: no net charge
        'PHE': [0, 5.48, 147.2, 5.27, 2.8, 5.2],  # Phenylalanine: no net charge
        'ASN': [0, 5.41, 114.1, -2.65, -3.5, 11.6],  # Asparagine: no net charge
        'GLY': [0, 5.97, 57.0, -1.62, -0.4, 9.0],  # Glycine: no net charge
        'HIS': [0.1, 7.59, 137.1, 1.28, -3.2, 10.4],  # Histidine: partial positive charge (~0.1)
        'LEU': [0, 5.98, 113.2, 5.01, 3.8, 4.9],  # Leucine: no net charge
        'ARG': [1, 10.76, 156.2, -0.93, -4.5, 10.5],  # Arginine: 1 positive charge
        'TRP': [0, 5.89, 186.2, 5.2, -0.9, 5.4],  # Tryptophan: no net charge
        'ALA': [0, 6.00, 71.1, -0.22, 1.8, 8.1],  # Alanine: no net charge
        'VAL': [0, 5.96, 99.1, 4.45, 4.2, 5.9],  # Valine: no net charge
        'GLU': [-1, 3.22, 129.1, -3.64, -3.5, 12.3],  # Glutamic Acid: 1 negative charge
        'TYR': [0, 5.66, 163.2, 2.15, -1.3, 6.2],  # Tyrosine: no net charge
        'MET': [0, 5.74, 131.2, 3.51, 1.9, 5.7]  # Methionine: no net charge
    }

    residue_data = []
    
    # Iterate over residues
    for residue in model.get_residues():
        res_id = residue.get_id()[1]
        res_name = residue.get_resname()
        na = len(residue)
        com = tuple(residue.center_of_mass())
        if res_name in dict_physicochemical_characteristics:
            residue_data.append([
                res_id,  # Position
                res_name,  # Residue
                na,  # Number of atoms
                *(dict_physicochemical_characteristics[res_name]),  # Number of electrostatic charges,
                    # pI, Mass, Enc, Hydrophobicty (Kyte-Doolittle scale), Polarity (Grantham scale) 
                com[0], 
                com[1],
                com[2]
            ])

    df = pd.DataFrame(residue_data, columns=[
        "Position", "Residue", "Na", "Nec", "pI", "Mass", "Enc", "Hidrofobicity",
        "Polarity", "Center_of_mass_x", "Center_of_mass_y", "Center_of_mass_z"
    ])
    return df.sort_values(by="Position")  # Sort by position




### Function to extract input features
def physicochemical_feature (file):
    try:
        physicochemical_data = get_physicochemical_feature(file)
        return physicochemical_data

    except Exception as e:
        print(f"Error calculating physicochemical data: {e}")


### For executing as script
if __name__ == "__main__":
    file = sys.argv[1]
    physicochemical_data = get_physicochemical_feature(file)
    output_csv = os.path.join(BASE_DIR, "physicochemical.csv")
    physicochemical_data.to_csv(output_csv, index=False)
    print(f"Physicochemical data saved in {output_csv.split('/')[-1]}")