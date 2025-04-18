import subprocess
import sys
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def run_chimera(pdb_file):
    print("Executing chimera")
    # Ruta al script de Chimera para resaltar binding sites
    highlight_script = os.path.join(BASE_DIR, 'highlight_binding_sites.py')

    # Ruta completa al ejecutable de Chimera
    chimera_path = '/Applications/Chimera.app/Contents/MacOS/chimera'

    # Ejecutar Chimera con el script y el archivo PDB
    subprocess.run([chimera_path, '--script', highlight_script, pdb_file])

if __name__ == "__main__":
    pdb_file = sys.argv[1]
    run_chimera(pdb_file)