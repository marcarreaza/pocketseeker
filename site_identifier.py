import os
import pandas as pd
from Bio import PDB

for folder in os.listdir('./input'):
    folder_path = os.path.join('./input', folder)
    if os.path.isdir(folder_path):
        # Construct the path to the protein.pdb inside the folder
        protein_pdb = os.path.join(folder_path, 'protein.pdb')
        binding_site_pdb = os.path.join(folder_path, 'site.pdb') 
        
        # Check if the protein.pdb file exists
        if os.path.exists(protein_pdb):
            print(f"Processing {protein_pdb}")
            # Processing code
            # Load PDB structures
            parser = PDB.PDBParser(QUIET=True)
            protein_structure = parser.get_structure("protein", protein_pdb)
            binding_site_structure = parser.get_structure("binding_site", binding_site_pdb)

            # Extract residue positions for binding site
            binding_residues = list()
            c = 0
            for residue in binding_site_structure.get_residues():
                res_name = residue.get_resname().strip()
                center_of_mass = residue.center_of_mass()
                binding_residues.append((res_name, tuple(center_of_mass)))
                c +=1

            if len(binding_residues) != c:
                raise ValueError(f"It's not identifying well the residues and the positions for {protein}")
            
            # Create matrix of residues with binding site labels
            residue_data = []
            for residue in protein_structure.get_residues():
                res_id = residue.get_id()[1]  # Residue number
                res_name = residue.get_resname().strip()  # Residue name
                center_of_mass = residue.center_of_mass()
                is_binding_site = "Yes" if (res_name, tuple(center_of_mass)) in binding_residues else "No"
                residue_data.append([res_id, res_name, is_binding_site])
            
            # Convert to DataFrame
            df = pd.DataFrame(residue_data, columns=["Position", "Residue", "Binding_Site"])
            df = df.sort_values(by="Position")  # Sort by position
            
            # Save as CSV
            output_file = os.path.join(folder_path, "residue_matrix.csv")
            df.to_csv(output_file, index=False)
            print(f"Matrix saved to {output_file}")

        else:
            raise ValueError(f"protein.pdb not found in {folder}")