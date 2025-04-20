import os
import pandas as pd
from Bio import PDB


def site_identifier (protein_pdb, binding_site_pdb = None, binding_site_txt = None ):
    if not os.path.exists(protein_pdb):
        raise FileNotFoundError(f"{protein_pdb} not found.")
    print(f"Processing binding sites {protein_pdb}")
    parser = PDB.PDBParser(QUIET=True)
    protein_structure = parser.get_structure("protein", protein_pdb)

    binding_positions = list()

    if binding_site_pdb:
        binding_structure = parser.get_structure("binding_site", binding_site_pdb)
        for residue in binding_structure.get_residues():
            res_name = residue.get_resname().strip()  # Residue number
            center_of_mass = residue.center_of_mass()
            binding_positions.append((res_name, tuple(center_of_mass)))

    elif binding_site_txt:
        if not os.path.exists(binding_site_txt):
            raise FileNotFoundError(f"{binding_site_txt} not found.")

        with open(binding_site_txt, "r") as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        pos, aa = line.split("\t")
                        pos = int(pos)
                        if not aa.isalpha() or len(aa) != 1:
                            raise ValueError
                        binding_positions.append((aa, pos))
                    except ValueError:
                        raise ValueError(f"Invalid line format in {binding_site_txt}: '{line}'")

    # Build residue information
    residue_data = []
    for residue in protein_structure.get_residues():
        res_id = residue.get_id()[1]  # Residue number
        res_name = residue.get_resname().strip()  # Residue name
        center_of_mass = residue.center_of_mass()
        is_binding_site = "Yes" if (res_name, tuple(center_of_mass)) in binding_positions else "No"
        residue_data.append([res_id, res_name, is_binding_site])

    df = pd.DataFrame(residue_data, columns=["Position", "Residue", "Binding_Site"])
    df = df.sort_values(by="Position")
    return df