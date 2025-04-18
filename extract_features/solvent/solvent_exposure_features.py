import os
import sys
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

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

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
    try:
        surface = get_surface(model, MSMS=os.path.join(BASE_DIR, '../../programs/msms_mac/msms.x86_64Darwin.2.6.1'))
    except:
        surface = get_surface(model, MSMS=os.path.join(BASE_DIR,'../../programs/msms/msms.x86_64Linux2.2.6.1.staticgcc'))
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
                data.append([res_id, res_name, depth, 
                             #hse_CA_value, hse_CB_value, 
                             hse_CN_value])
    
    # Return a DataFrame with the residue depth results
    return pd.DataFrame(data, columns=["Position", "Residue", "Residue_Depth", 
                                       #"HSE-up", "HSE-down", 
                                       "CN"]).sort_values(by="Position")


### Function to extract input features
def solvent_feature (file):
    try:
        solvent_data = compute_solvent_exposure(file)
        return solvent_data
    
    except Exception as e:
        print(f"Error extracting solvent exposure data: {e}")


### For executing as script
if __name__ == "__main__":
    file = sys.argv[1]
    solvent_data = compute_solvent_exposure(file)
    output_csv = os.path.join(BASE_DIR, "solvent_exposerue.csv")
    solvent_data.to_csv(output_csv, index=False)
    print(f"Solvent exposure data saved in {output_csv.split('/')[-1]}")