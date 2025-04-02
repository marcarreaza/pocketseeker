import csv
import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser, PPBuilder
import os

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

# ------ Procesamiento de múltiples archivos PDB ------
def process_pdb_files(pdb_files, output_file="SS_features.csv"):
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["File", "Res", "SS", "ASA", "Phi", "Psi", "Theta(i-1=>i+1)", "Tau(i-2=>i+2)"])
        
        for pdb_file in pdb_files:
            print(f"Procesando {pdb_file}...")
            file_name = os.path.basename(pdb_file)
            SS = SS_extraction(pdb_file)
            ASA = ASA_extraction(pdb_file)
            PHI_PSI = phi_psi_extraction(pdb_file)
            THETA_TAU = theta_tau_extraction(pdb_file)
            
            for res_id in SS:
                writer.writerow([pdb_file, res_id, SS[res_id], ASA.get(res_id, "-"), PHI_PSI.get(res_id, ("-", "-"))[0],
                                 PHI_PSI.get(res_id, ("-", "-"))[1], THETA_TAU.get(res_id, ("-", "-"))[0],
                                 THETA_TAU.get(res_id, ("-", "-"))[1]])
    print(f"Archivo CSV guardado como {output_file}")

# ------ Ejecutar el script ------
if __name__ == "__main__":
    pdb_files = ["../../input/1a2b_1/protein.pdb", 
                 "../../input/1a2n_1/protein.pdb"
    ] 
    process_pdb_files(pdb_files)
