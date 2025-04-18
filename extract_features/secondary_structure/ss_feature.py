import sys
import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser, PPBuilder
import os
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Diccionario de conversión de aminoácidos
d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

# ------ Funciones -------
def open_pdb(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    return structure

def SS_extraction(pdb_file):
    structure = open_pdb(pdb_file)
    model = structure[0]
    dssp = PDB.DSSP(model, pdb_file)
    return {key[1][1]: dssp[key][2] for key in dssp.keys()}

def ASA_extraction(pdb_file):
    structure = open_pdb(pdb_file)
    model = structure[0]
    dssp = PDB.DSSP(model, pdb_file)
    return {key[1][1]: dssp[key][3] for key in dssp.keys()}

def phi_psi_extraction(pdb_file):
    structure = open_pdb(pdb_file)
    PHI = {}
    for model in structure:
        ppb = PPBuilder()
        for pp in ppb.build_peptides(model):
            phi_psi_list = pp.get_phi_psi_list()
            for res, (phi, psi) in zip(pp, phi_psi_list):
                PHI[res.get_id()[1]] = (f"{phi:.2f}" if phi else "-", f"{psi:.2f}" if psi else "-")
    return PHI

def angle_between_three_points(p1, p2, p3):
    v1 = p1 - p2
    v2 = p3 - p2
    cosine_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.arccos(cosine_angle)
    return np.degrees(angle)

def theta_tau_extraction(pdb_file):
    structure = open_pdb(pdb_file)
    THETA = {}
    for model in structure:
        for chain in model:
            residues = [res for res in chain.get_residues() if res.get_resname() in d]
            for i in range(len(residues)):
                res_i = residues[i]
                if i in [0, len(residues)-1]:
                    THETA[res_i.get_id()[1]] = ("-", "-")
                    continue
                
                ca_i = residues[i]["CA"].get_coord()
                ca_i_minus_1 = residues[i-1]["CA"].get_coord()
                ca_i_plus_1 = residues[i+1]["CA"].get_coord()
                theta = angle_between_three_points(ca_i_minus_1, ca_i, ca_i_plus_1)
                
                if i not in [1, len(residues)-2]:
                    ca_i_minus_2 = residues[i-2]["CA"].get_coord()
                    ca_i_plus_2 = residues[i+2]["CA"].get_coord()
                    tau = angle_between_three_points(ca_i_minus_2, ca_i, ca_i_plus_2)
                    THETA[res_i.get_id()[1]] = (f"{theta:.2f}", f"{tau:.2f}")
                else:
                    THETA[res_i.get_id()[1]] = (f"{theta:.2f}", "-")
    return THETA



### Function to extract input features
def ss_feature (file):
    try:
        # Intentamos parsear el archivo PDB
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", file)

        SS = SS_extraction(file)
        ASA = ASA_extraction(file)
        PHI_PSI = phi_psi_extraction(file)
        THETA_TAU = theta_tau_extraction(file)

        data = []

        for res_id in SS:
            row = {
                "Res": res_id,
                "SS": SS[res_id],
                "ASA": ASA.get(res_id, "-"),
                "Phi": PHI_PSI.get(res_id, ("-", "-"))[0],
                "Psi": PHI_PSI.get(res_id, ("-", "-"))[1],
                "Theta(i-1=>i+1)": THETA_TAU.get(res_id, ("-", "-"))[0],
                "Tau(i-2=>i+2)": THETA_TAU.get(res_id, ("-", "-"))[1]
            }
            data.append(row)
        print("Calculating secondary structure")

        return pd.DataFrame(data)
        
    except Exception as e:
        print(f"Error calculating secondary structure: {e}")


### For executing as script
if __name__ == "__main__":
    file = sys.argv[1]
    ss_data = ss_feature(file)
    output_csv = os.path.join(BASE_DIR, "ss.csv")
    ss_data.to_csv(output_csv, index=False)
    print(f"Secondary structure data saved in {output_csv.split('/')[-1]}")