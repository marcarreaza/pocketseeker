import glob
import os
import mdtraj as md

for f in glob.glob("./input/*/protein.mol2"):
    print(f"Processing {f}...")
    
    pdb_file = f"{os.path.splitext(f)[0]}.pdb"  # Correct path handling

    try:
        traj = md.load_mol2(f)  # Use 'f' instead of 'mol2_file'
        traj.save_pdb(pdb_file)
        print(f"Converted {f} to {pdb_file}")
    except Exception as e:
        print(f"Failed to convert {f}: {e}")

for f in glob.glob("./input/*/site.mol2"):
    print(f"Processing {f}...")
    
    pdb_file = f"{os.path.splitext(f)[0]}.pdb"  # Correct path handling

    try:
        traj = md.load_mol2(f)  # Use 'f' instead of 'mol2_file'
        traj.save_pdb(pdb_file)
        print(f"Converted {f} to {pdb_file}")
    except Exception as e:
        print(f"Failed to convert {f}: {e}")