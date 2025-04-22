# PocketSeeker: A ligand binding site predictor for proteins

## 📌 Project Description

PocketSeeker is a bioinformatics tool designed to predict the binding site residues of a protein using machine learning techniques. The core of the application is a Random Forest classifier trained on structural and sequence-based features extracted from protein data. Given a protein structure in PDB format or a PDB ID, PocketSeeker analyzes each residue and determines its likelihood of being part of a binding site.

This tool aims to assist researchers in identifying potential functional regions of proteins, which can be crucial for tasks such as drug design, protein engineering, or functional annotation.


## ⚙️ Setup Instructions

To get started with PocketSeeker, follow the steps below to set up the project on your local machine.

### 1. Clone the repository
Clone the PocketSeeker repository using Git:

`git clone https://github.com/marcarreaza/pocketseeker.git`
`cd pocketseeker`

### 2. Run the setup script
Execute the provided setup script to create and configure the Conda environment:

`bash setup.sh`

This script automatically creates a Conda environment named pocketseeker with all required dependencies. It detects whether the system is running on Linux or macOS and installs the appropriate binaries and tools accordingly.

### 3. Activate the environment
Once the setup is complete, activate the environment manually:

`conda activate pocketseeker`

### 4. (Optional) Add PocketSeeker to your system PATH
To run pocketseeker.py from any location as a command-line tool, you can add the project folder to your system PATH. For example:

`export PATH=$PATH:/path/to/pocketseeker`

Replace /path/to/pocketseeker with the absolute path to the cloned repository. You can also add this line to your .bashrc, .zshrc, or equivalent shell config file to make it permanent.