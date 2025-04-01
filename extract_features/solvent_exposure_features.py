import os
import pandas as pd
from Bio.PDB import PDBParser, HSExposure
from Bio import PDB
from Bio.PDB.ResidueDepth import residue_depth, get_surface

# Computing solvent exposure features
#   - Residue depth. It is the average distance of the atoms of a residue from the solvent accessible surface.
#       - Using command line tool MSMS
#         This module uses Michel Sanner’s MSMS program for the surface calculation. See: http://mgltools.scripps.edu/packages/MSMS
#   - Number of Cα atoms in upper half-sphere (HSE-up)
#   - Number of Cα atoms in lower half-sphere (HSE-down)
#   - Contact number (CN)

# Define residue mapping
residue_map = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def get_HSE(model, RADIUS):
    # Calculate the solvent exposure for CA, CB, and CN atoms with a specified radius
    HSE_up = PDB.HSExposure.HSExposureCA(model, RADIUS)
    HSE_down = PDB.HSExposure.HSExposureCB(model, RADIUS)
    CN = PDB.HSExposure.ExposureCN(model, RADIUS)
    return [HSE_up, HSE_down, CN]

def compute_solvent_exposure(pdb_file):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]
    surface = get_surface(model, MSMS='../programs/msms/msms.x86_64Linux2.2.6.1.staticgcc')
    HSE = get_HSE(model, 15.0) 
    # List to store information
    data = []
    
    # Loop through each chain and residue
    for chain in model.get_chains():
        for residue in chain.get_residues():
            res_name = residue.get_resname()
            res_id = residue.get_id()[1]
            chain_id = residue.get_parent().id
            if res_name in residue_map:
                # Calculate residue depth for each residue
                depth = residue_depth(residue, surface)
                # Try to access the HSE value, if not found, set it to None
                try:
                    hse_CA_value = HSE[0][(chain_id, residue.get_id())]
                    hse_CB_value = HSE[1][(chain_id, residue.get_id())]
                    hse_CN_value = HSE[2][(chain_id, residue.get_id())]
                
                except KeyError:
                    hse_CA_value = None
                    hse_CB_value = None
                    hse_CN_value = None

                # Append residue depth information
                data.append([res_id, res_name, depth, hse_CA_value, hse_CB_value, hse_CN_value])
    
    # Return a DataFrame with the residue depth results
    return pd.DataFrame(data, columns=["Position", "Residue", "Residue_Depth", "HSE-up", "HSE-down", "CN"]).sort_values(by="Position")

# Loop through directories and process each protein
for folder in os.listdir('../input'):
    folder_path = os.path.join('../input', folder)
    if os.path.isdir(folder_path):
        # Path to the unbound (protein) PDB file
        protein_pdb = os.path.join(folder_path, 'protein.pdb')

        if os.path.exists(protein_pdb):
            print(f"Processing {protein_pdb}")

            # Compute residue depth for unbound protein
            residue_depth_df = compute_solvent_exposure(protein_pdb)

            # Save the residue depth results to a CSV file
            output_csv = os.path.join(folder_path, "solvent_exposure_features.csv")
            residue_depth_df.to_csv(output_csv, index=False)
            print(f"Saved results to {output_csv}")