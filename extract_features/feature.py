import os
import pandas as pd
from Bio import PDB
from Bio.PDB import PDBParser

from Bio.SeqUtils.ProtParamData import kd

from Bio.Alphabet import *

from Bio.SeqUtils.ProtParam import ProteinAnalysis



for folder in os.listdir('../input'):
    folder_path = os.path.join('../input', folder)
    if os.path.isdir(folder_path):
        # Construct the path to the protein.pdb inside the folder
        protein_pdb = os.path.join(folder_path, 'protein.pdb')
        # Check if the protein.pdb file exists
        if os.path.exists(protein_pdb):
            print(f"Processing {protein_pdb}")
            # Processing code
            # Load PDB structures
            parser = PDB.PDBParser(QUIET=True)
            protein_structure = parser.get_structure("protein", protein_pdb)
            d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

            dict_hidrofobicity = {aa: kd.get(one_letter_code) for aa, one_letter_code in d.items()}

            dict_physicochemical_characteristics = {
                'CYS': [5.05,103.1,4.66], 'ASP': [2.77,115.1,-4.12], 'SER': [5.68,87.1,-2.84],
                'GLN': [5.65,128.1,-2.76], 'LYS': [9.74,128.2,-4.18], 'ILE': [6.02,113.2,5.58],
                'PRO': [6.3,97.1,-3.03], 'THR': [5.66,101.1,-1.2], 'PHE': [5.48,147.2,5.27],
                'ASN': ,'GLY': , 'HIS': , 'LEU':
            }

            residue_data = []
            for residue in protein_structure.get_residues():
                res_id = residue.get_id()[1]
                res_name = residue.get_resname()
                center_of_mass = tuple(residue.center_of_mass())
                hidrofobicity = dict_hidrofobicity.get(res_name)

                residue_data.append([res_id, res_name, center_of_mass, hidrofobicity])

            df = pd.DataFrame(residue_data, columns=["Position", "Residue", "Center_of_mass", "Hidrofobicity"])
            df = df.sort_values(by="Position")  # Sort by position
                        
            # Save as CSV
            output_file = os.path.join(folder_path, "features.csv")
            df.to_csv(output_file, index=False)
            print(f"Matrix saved to {output_file}")



'''
    #compute andle of all the atoms
    for atom in residue.get_atoms():
        vector = atom.get_vector()
    angle = calc_ang
'''

