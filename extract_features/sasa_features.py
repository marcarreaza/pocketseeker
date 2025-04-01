import os
import pandas as pd
from Bio import PDB
from Bio.PDB import DSSP, ShrakeRupley

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
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]
    sa = ShrakeRupley()
    sa.compute(model, level="A")  # Compute SASA at atom level
    sa.compute(model, level = "R")
    #try:
    #    dssp = DSSP(model, pdb_file)
    #    print(len(dssp))
    #except Exception as e:
    #    raise ValueError(f"DSSP failed on {pdb_file}: {e}")
    #print(dssp.keys())  # Check if any residues exist

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

# Loop through directories and process each protein
for folder in os.listdir('../input'):
    folder_path = os.path.join('../input', folder)
    if os.path.isdir(folder_path):
        # Path to the unbound (protein) PDB file
        protein_pdb = os.path.join(folder_path, 'protein.pdb')

        if os.path.exists(protein_pdb):
            print(f"Processing {protein_pdb}")

            # Compute SASA for unbound protein
            sasa_unbound_df = compute_sasa(protein_pdb)

            # Save the SASA results to a CSV file
            output_csv = os.path.join(folder_path, "sasa_features.csv")
            sasa_unbound_df.to_csv(output_csv, index=False)
            print(f"Saved results to {output_csv}")



'''
residue_data.append([
                            res_id,  # Position
                            res_name,  # Residue
                            tuple(residue.center_of_mass()),  # Center of mass
                            *(dict_physicochemical_characteristics[res_name]),  # pI, Mass, Enc, Hydrophobicty (Kyte-Doolittle scale), Polarity (Grantham scale)
                            round(residue.sasa, 2),  # SASA_first
                            dssp_values[2],  # Secondary Structure
                            round(dssp_values[3], 2),  # Relative ASA
                            dssp_values[4],  # Phi
                            dssp_values[5],  # Psi
                            dssp_values[6],  # NH→O_1_relidx
                            dssp_values[7],  # NH→O_1_energy
                            dssp_values[8],  # O→NH_1_relidx
                            dssp_values[9],  # O→NH_1_energy
                            dssp_values[10], # NH→O_2_relidx
                            dssp_values[11], # NH→O_2_energy
                            dssp_values[12], # O→NH_2_relidx
                            dssp_values[13], # O→NH_2_energy
                            hse_value,  # Half-sphere exposure value
                                # The first number: the count of surrounding atoms (coordination number).
                                # The second number: the number of atoms within the radius.
                                # The third number: the computed HSE value (e.g., np.float64(1.0375910569929647)).
                            hse_CA_value,
                            hse_CB_value,
                            hse_CN_value
'''