"""
Automated Blind Molecular Docking Pipeline
Usage: python autodock_blind.py --smiles "<SMILES>" --pdb_id <PDB_ID>
"""

import argparse
import glob
import os
import re
import subprocess
import sys
import textwrap
import urllib.parse
import urllib.request
from pathlib import Path


def run(cmd, *, cwd=None):
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        sys.exit(f"[ERROR] Command failed: {' '.join(str(c) for c in cmd)}")


def write_script(path, content):
    path.write_text(textwrap.dedent(content))


def find_local(names):
    script_dir = Path(__file__).resolve().parent
    for name in names:
        p = script_dir / name
        if p.exists():
            return str(p)
    for name in names:
        if subprocess.run(["which", name], capture_output=True).returncode == 0:
            return name
    return None


def chimerax_exec(arg):
    if arg:
        return arg
    result = find_local(["chimerax", "ChimeraX"])
    if result:
        return result
    for pattern in [
        "/Applications/ChimeraX*.app/Contents/MacOS/ChimeraX",
        os.path.expanduser("~/Applications/ChimeraX*.app/Contents/MacOS/ChimeraX"),
    ]:
        matches = sorted(glob.glob(pattern))
        if matches:
            return matches[-1]
    sys.exit(
        "[ERROR] ChimeraX not found. Create a symlink in this folder:\n"
        "  ln -s /Applications/ChimeraX-X.Y.Z.app/Contents/MacOS/ChimeraX ./chimerax"
    )


def vina_exec(arg):
    if arg and arg != "vina":
        return str(Path(arg).resolve())
    result = find_local(["vina", "vina_1.2.7_mac_arm64", "vina_1.2.7_mac_x86_64",
                         "vina_1.2.5_mac_arm64", "vina_1.2.5_mac_x86_64"])
    if result:
        return result
    sys.exit("[ERROR] Vina not found. Place the vina executable in the same folder as this script.")

def get_smiles_from_name(lig_name):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/property/CanonicalSMILES/TXT"
    url = f"{base_url}?name={urllib.parse.quote(lig_name)}"

    try:
        with urllib.request.urlopen(url, timeout=20) as resp:
            smiles = resp.read().decode().strip()
            if smiles:
                return smiles.splitlines()[0]
    except Exception:
        pass

    sys.exit(f"[ERROR] Could not fetch SMILES for ligand name: {lig_name}")
def download_ligand_sdf(smiles, workdir):
    sdf_path = workdir / "ligand.sdf"
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/SDF"
    for record_type in ["3d", None]:
        params = {"smiles": smiles}
        if record_type:
            params["record_type"] = record_type
        try:
            req = urllib.request.Request(
                base_url,
                data=urllib.parse.urlencode(params).encode(),
                headers={"Content-Type": "application/x-www-form-urlencoded"},
            )
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = resp.read()
            if data.strip():
                sdf_path.write_bytes(data)
                return sdf_path
        except Exception:
            continue
    sys.exit("[ERROR] Could not download ligand SDF from PubChem. Check your SMILES string.")


def prepare_ligand(chimerax, workdir):
    mol2_out = workdir / "ligand_prepared.mol2"
    script = workdir / "_prep_ligand.cxc"
    write_script(script, f"""\
        open {workdir / "ligand.sdf"}
        addh
        save {mol2_out} format mol2
        exit
    """)
    run([chimerax, "--nogui", "--script", str(script)], cwd=workdir)
    if not mol2_out.exists():
        sys.exit("[ERROR] ChimeraX did not produce ligand_prepared.mol2")
    return mol2_out


def prepare_receptor(chimerax, pdb_id, workdir):
    pdb_out = workdir / "receptor_prepared.pdb"
    script = workdir / "_prep_receptor.cxc"
    write_script(script, f"""\
        open {pdb_id}
        delete solvent
        delete ~protein
        addh
        save {pdb_out} format pdb
        exit
    """)
    run([chimerax, "--nogui", "--script", str(script)], cwd=workdir)
    if not pdb_out.exists():
        sys.exit("[ERROR] ChimeraX did not produce receptor_prepared.pdb")
    return pdb_out


def convert_to_pdbqt(input_file, output_file, is_receptor):
    flag = ["-xr"] if is_receptor else ["-xh"]
    run(["obabel", str(input_file), "-O", str(output_file)] + flag)
    if not output_file.exists():
        sys.exit(f"[ERROR] OpenBabel did not produce {output_file.name}")


def get_center(chimerax, receptor_pdb, ligand_mol2, workdir):
    out_file = workdir / "_center.txt"
    script = workdir / "_measure_center.cxc"
    write_script(script, f"""\
        open {receptor_pdb}
        open {ligand_mol2}
        select all
        measure center sel
        log save {out_file}
        exit
    """)
    result = subprocess.run(
        [chimerax, "--nogui", "--script", str(script)],
        capture_output=True, text=True, cwd=workdir,
    )
    combined = result.stdout + result.stderr
    if out_file.exists():
        combined += out_file.read_text()
    match = re.search(
        r"Center of mass.*?=\s*\(\s*([-\d.]+),\s*([-\d.]+),\s*([-\d.]+)\s*\)",
        combined,
    )
    if not match:
        sys.exit("[ERROR] Could not parse center of mass from ChimeraX output.")
    return float(match.group(1)), float(match.group(2)), float(match.group(3))


def write_conf(workdir, receptor_pdbqt, ligand_pdbqt, center, exhaustiveness):
    cx, cy, cz = center
    (workdir / "conf.txt").write_text(
        f"receptor = {receptor_pdbqt.name}\n"
        f"ligand = {ligand_pdbqt.name}\n"
        f"\n"
        f"center_x = {cx:.2f}\n"
        f"center_y = {cy:.2f}\n"
        f"center_z = {cz:.2f}\n"
        f"\n"
        f"size_x = 50\n"
        f"size_y = 50\n"
        f"size_z = 50\n"
        f"\n"
        f"exhaustiveness = {exhaustiveness}\n"
    )


def main():
    p = argparse.ArgumentParser(description="Automated blind molecular docking")
    p.add_argument("--smiles", help="Ligand SMILES string")
    p.add_argument("--lig", help="Ligand name (e.g., 'folic acid', 'cocaine')")
    p.add_argument("--pdb_id",         required=True)
    p.add_argument("--workdir",        default="docking_output")
    p.add_argument("--chimerax",       default="")
    p.add_argument("--vina",           default="vina")
    p.add_argument("--exhaustiveness", type=int, default=8)
    args = p.parse_args()

    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    chimerax = chimerax_exec(args.chimerax)
    vina     = vina_exec(args.vina)

    if args.smiles:
        smiles = args.smiles
    elif args.lig:
        print(f"[INFO] Fetching SMILES for ligand: {args.lig}")
        smiles = get_smiles_from_name(args.lig)
    else:
        sys.exit("[ERROR] Provide either --smiles or --lig")
    print(f"Ligand  : {smiles}")
    print(f"Receptor: {args.pdb_id}")
    print(f"Workdir : {workdir}")

    ligand_sdf = download_ligand_sdf(smiles, workdir)
    ligand_mol2    = prepare_ligand(chimerax, workdir)
    receptor_pdb   = prepare_receptor(chimerax, args.pdb_id, workdir)

    ligand_pdbqt   = workdir / "ligand.pdbqt"
    receptor_pdbqt = workdir / "receptor.pdbqt"
    convert_to_pdbqt(ligand_mol2,  ligand_pdbqt,   is_receptor=False)
    convert_to_pdbqt(receptor_pdb, receptor_pdbqt, is_receptor=True)

    center = get_center(chimerax, receptor_pdb, ligand_mol2, workdir)
    print(f"Center  : {center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}")

    write_conf(workdir, receptor_pdbqt, ligand_pdbqt, center, args.exhaustiveness)

    run([vina, "--config", str(workdir / "conf.txt"), "--out", str(workdir / "results.pdbqt")],
        cwd=workdir)

    print(f"Done. Results: {workdir / 'results.pdbqt'}")


if __name__ == "__main__":
    main()