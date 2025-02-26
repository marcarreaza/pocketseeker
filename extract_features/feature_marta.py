import os
import sys
import pandas as pd
import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.SeqUtils.ProtParamData import kd

# Data needed
d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

# Functions
def extract_sequence(pdb_file):
    if os.path.exists(pdb_file):
        print(f"===> Processing {pdb_file}")
        parser = PDB.PDBParser(QUIET=True)
        protein_structure = parser.get_structure("protein", pdb_file)
        sequence = []
        for model in protein_structure:
            for chain in model:
                for residue in chain:
                    if PDB.is_aa(residue, standard=True):
                        sequence.append(d[residue.get_resname()])
        return "".join(sequence)

seq = extract_sequence("../input/1a2b_1/protein.pdb")
with open("protein.fasta", "w") as f:
    f.write(">protein\n" + seq + "\n")



# Executable
for folder in os.listdir('../input'):
    folder_path = os.path.join('../input', folder)
    if os.path.isdir(folder_path):
        # Construct the path to the protein.pdb inside the folder
        protein_pdb = os.path.join(folder_path, 'protein.pdb')
