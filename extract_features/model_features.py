import pandas as pd
import os
import sys
from Bio.PDB import PDBParser, PDBIO

from concativity.concavity_feature import concavity_feature
from core_distance.core_distance_feature import core_distance
from physicochemical.physicochemical_features import physicochemical_feature
from PSSM.PSSM_feature import PSSM_feature
from SASA.sasa_features import sasa_feature
from secondary_structure.ss_feature import ss_feature
from solvent.solvent_exposure_features import solvent_feature
from site_identifier import site_identifier

               
def extract_features (protein_pdb):
    print(protein_pdb)
    # Cargar el archivo PDB
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", protein_pdb)

    # Filter by SASA and Rsidual Depth before extract the others features
    sasa = sasa_feature (protein_pdb)
    filtered_sasa = sasa[(sasa[sasa.columns[sasa.columns.get_loc('SASA'):]]==0).all(axis=1)]
    pos1 = list(filtered_sasa['Position'])

    solvent = solvent_feature (protein_pdb)
    filtered_solvent = solvent[(solvent["Residue_Depth"]>4.0)]
    pos2 = list(filtered_solvent['Position'])

    positions = list(set(pos1) & set(pos2))

    # Crear una nueva estructura eliminando los residuos no deseados
    for model in structure:
        for chain in model:
            residues_a_quitar = []
            for residue in chain:
                res_id = residue.get_id()[1]  # El ID del residuo (número)
                if res_id in positions:
                    residues_a_quitar.append(residue)
            # Eliminar los residuos marcados
            for residue in residues_a_quitar:
                chain.detach_child(residue.id)

    # Guardar el nuevo archivo PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save("filter_pdb.pdb")

    concavity = concavity_feature("filter_pdb.pdb")
    distance_to_core = core_distance("filter_pdb.pdb")
    physicochemical = physicochemical_feature ("filter_pdb.pdb")
    pssm = PSSM_feature ("filter_pdb.pdb")
    ss = ss_feature ("filter_pdb.pdb")
    sasa = sasa_feature("filter_pdb.pdb")
    solvent = solvent_feature ("filter_pdb.pdb")

    df_final = pd.concat([
        concavity,
        distance_to_core["Distance_to_Core"],
        physicochemical.drop(columns=["Position", "Residue"]),
        pssm.drop(columns=["File", "Res"]),
        sasa.drop(columns=["Position", "Residue"]),
        ss.drop(columns=["Res"]),
        solvent.drop(columns=["Position", "Residue"]),
    ], axis=1)

    # Delete the intermediate files 
    files = ['filter_pdb.pdb', 'protein.fasta', 'protein.pssm']
    for file in files:
        if os.path.exists(file):
            os.remove(file)

    return df_final

if __name__ == "__main__":
    dfs = []
    for folder in os.listdir('../input'):
            folder_path = os.path.join('../input', folder)
            if os.path.isdir(folder_path):
                # Path to the unbound (protein) PDB file
                protein_pdb = os.path.join(folder_path, 'protein.pdb')
                binding_site_pdb = os.path.join(folder_path, 'site.pdb')
            
                if os.path.exists(protein_pdb):
                    print(f"Processing {protein_pdb}")

                    df_final = extract_features (protein_pdb)
                    binding_sites = site_identifier(protein_pdb, binding_site_pdb)
                    binding_sites = binding_sites.drop(columns=["Residue"])
                    binding_sites["Binding_Site"] = binding_sites["Binding_Site"].replace({"Yes" : 1, "No" : 0})

                    df_merged = df_final.merge(binding_sites, how='left', left_on='Res', right_on='Position')
                    df_merged = df_merged.drop(columns='Position')
                    df_merged = df_merged.drop(columns=["File"])

                    df_merged.to_csv(os.path.join(folder_path, "features.csv"), index=False)

                    dfs.append(df_merged)

    df_final = pd.concat(dfs, ignore_index=True)
    df_final.to_csv("total_features.csv", index=False)