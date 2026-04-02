# Detailed Molecular Docking Pipeline

### 1. Initialization & Dependency Check
* **Environment Setup**: The script uses `pathlib` to resolve absolute paths, ensuring it can find the `vina` binary and `chimerax` symlink regardless of where it is launched.
* **Binary Verification**: It runs a check for `chimerax`, `vina`, and `obabel`. If any are missing from the local folder or the system's `PATH`, the script terminates immediately to prevent partial execution.

### 2. Ligand Data Acquisition (PubChem REST API)
* **SMILES Retrieval**: If a ligand name (e.g., "Caffeine") is provided instead of a SMILES string, the script issues a `GET` request to the **PubChem PUG REST** API to fetch the Canonical SMILES.
* **3D SDF Download**: The script sends a `POST` request to PubChem to download the ligand in **SDF** format. It specifically requests a `3d` record type to ensure the starting coordinates are spatial rather than a flat 2D representation.

### 3. Receptor Preparation (ChimeraX Headless)
* **Automation Script**: The script writes a temporary `.cxc` (ChimeraX command) file.
* **Cleaning Protocol**:
    * `open [PDB_ID]`: Fetches the protein structure directly from the RCSB PDB.
    * `delete solvent`: Removes water molecules that could block potential binding sites.
    * `delete ~protein`: Strips away ions, co-crystallized ligands, and non-amino acid chains.
    * `addh`: Adds necessary hydrogen atoms to the protein to satisfy valency for energy calculations.
* **Export**: Saves the refined protein as `receptor_prepared.pdb`.

### 4. Ligand Preparation (ChimeraX Headless)
* **Protonation**: The downloaded SDF is loaded into ChimeraX.
* **Hydrogen Addition**: Hydrogens are added to the ligand to ensure correct chemical properties for docking.
* **Export**: The ligand is saved as `ligand_prepared.mol2`. The `.mol2` format is used here because it preserves bond-type information (single, double, aromatic) better than standard PDB.

### 5. Format Conversion (OpenBabel)
* **PDBQT Generation**: AutoDock Vina requires the `.pdbqt` format, which includes partial charges ($Q$) and atom types ($T$).
* **Receptor Conversion**: `obabel` converts the PDB to PDBQT using the `-xr` flag, designating it as a rigid receptor.
* **Ligand Conversion**: `obabel` converts the MOL2 to PDBQT using the `-xh` flag. This identifies rotatable bonds and defines the "torsion tree," allowing the ligand to be flexible during docking.

### 6. Search Space Definition (Blind Docking)
* **Center of Mass Calculation**: The script executes a measurement command in ChimeraX to find the geometric center $(x, y, z)$ of the entire protein.
* **Configuration File**: A `conf.txt` is generated with:
    * **Center**: The calculated coordinates from the receptor.
    * **Size**: A large $50 \times 50 \times 50$ Å box. This is the "Blind" aspect, ensuring the box is large enough to cover most protein surfaces.
    * **Exhaustiveness**: Set to `8` (default), determining how thoroughly the algorithm searches for the global minimum energy.

### 7. Vina Execution & Output
* **Simulation**: The script calls the `vina` executable, feeding it the `conf.txt`, the receptor, and the ligand.
* **Global Optimization**: Vina uses a Lamarckian Genetic Algorithm and empirical scoring to find the best-fit poses.
* **Results**: The final binding poses, ranked by affinity ($kcal/mol$), are saved to `results.pdbqt`.