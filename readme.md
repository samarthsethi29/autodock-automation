# Automated Blind Molecular Docking Pipeline

### Introduction
This Python-based pipeline automates the structural preparation and execution of molecular docking. By providing a Protein Data Bank (PDB) ID and a SMILES string (or ligand name), the script handles structural retrieval, cleaning, protonation, and format conversion required to run **AutoDock Vina** in a headless environment.

---

### What it Does
1.  **Ligand Acquisition**: Fetches Canonical SMILES via PubChem API if a name is provided.
2.  **Structural Retrieval**: Downloads 3D SDF structures from PubChem and PDB files from RCSB.
3.  **Protein & Ligand Prep**: Uses ChimeraX to remove solvents/non-protein atoms and add hydrogens ($H$).
4.  **Format Conversion**: Employs OpenBabel to generate `.pdbqt` files (adding partial charges and defining rotatable bonds).
5.  **Blind Docking Setup**: Automatically calculates the receptor's center of mass to center a $50 \times 50 \times 50$ Å docking grid.
6.  **Simulation**: Executes AutoDock Vina and outputs the top binding poses.

---

### Prerequisites & Step-by-Step Setup

#### 1. System Requirements
Ensure the following are installed and accessible via your terminal:
* **Python 3.x**
* **OpenBabel**: 
    * *macOS*: `brew install open-babel`
    * *Linux*: `sudo apt-get install openbabel`
* **ChimeraX**: [Download UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/download.html)
* **AutoDock Vina**: [Download Vina Executable](https://github.com/ccsb-scripps/AutoDock-Vina/releases)

#### 2. Local Directory Configuration (Crucial)
The script is designed to look for executables in its immediate directory. You must set up the following:

* **AutoDock Vina**: Place the `vina` executable file inside the same folder as the script.
* **ChimeraX Symlink**: Create a symbolic link so the script can invoke the ChimeraX engine.
    * **macOS Command**:
      ```bash
      ln -s /Applications/ChimeraX-1.8.app/Contents/MacOS/ChimeraX ./chimerax
      ```
      *(Note: Update the version number `1.8` to match your installed version).*

#### 3. Necessary Libraries
The script uses standard Python libraries:
* `argparse`: Command-line argument parsing.
* `urllib.request`: API calls to PubChem.
* `subprocess`: Running CLI tools (Vina, ChimeraX, Obabel).
* `pathlib`: Cross-platform file path management.

---

### Usage

Run the script from your terminal using one of the two methods below:

**Method A: Using a Ligand Name**
```bash
python autodock_blind.py --lig "ibuprofen" --pdb_id 1WNY
```

**Method B: Using Smiles String**
```bash
python autodock_blind.py --smiles "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" --pdb_id 1WNY
```