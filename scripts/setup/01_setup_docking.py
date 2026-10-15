#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Phase 1: Molecular Docking -- BMS-626529 into gp120 relay pocket
PDB: 6MEO | Target: Chain G residues 465-472
"""
import subprocess, os, sys, json
from pathlib import Path

WORK = Path.home() / "md-relay-pocket"
INPUT = WORK / "input"
DOCK = WORK / "docking"

RELAY_RESIDUES = list(range(465, 473))  # G:465-472

BMS626529_SMILES = "COc1cnc(C(=O)N2CCN(c3cccc4[nH]nc(C)c34)CC2)c(OC[C@@H]2CC(F)(F)CO2)c1"
BMS626529_NAME = "BMS-626529"

def run(cmd, **kw):
    print(f"  >> {cmd}")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        print(f"  STDERR: {r.stderr[:500]}")
    return r

def compute_pocket_center(pdb_path):
    """Extract center of mass of relay pocket from full PDB."""
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chain = line[21]
            resnum = int(line[22:26].strip())
            if chain == "G" and resnum in RELAY_RESIDUES:
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
    if not xs:
        print("ERROR: No relay pocket atoms found in PDB")
        sys.exit(1)
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)
    print(f"  Relay pocket center: ({cx:.2f}, {cy:.2f}, {cz:.2f}) from {len(xs)} atoms")
    return cx, cy, cz

def prepare_ligand():
    """Convert BMS-626529 SMILES to 3D PDBQT via OpenBabel."""
    smi_file = DOCK / "bms626529.smi"
    mol2_file = DOCK / "bms626529.mol2"
    pdbqt_file = DOCK / "bms626529.pdbqt"

    with open(smi_file, "w") as f:
        f.write(f"{BMS626529_SMILES} {BMS626529_NAME}\n")

    # SMILES -> 3D mol2 with hydrogens and geometry optimization
    run(f"obabel {smi_file} -O {mol2_file} --gen3d --best -h")
    if not mol2_file.exists():
        print("ERROR: mol2 generation failed")
        sys.exit(1)

    # mol2 -> PDBQT for Vina
    run(f"obabel {mol2_file} -O {pdbqt_file} -xr")
    if not pdbqt_file.exists():
        print("ERROR: PDBQT conversion failed")
        sys.exit(1)

    print(f"  Ligand prepared: {pdbqt_file}")
    return pdbqt_file

def prepare_receptor(pdb_path):
    """Convert full PDB to PDBQT receptor for Vina, keeping only chain G near pocket."""
    receptor_pdb = DOCK / "relay_receptor.pdb"
    receptor_pdbqt = DOCK / "relay_receptor.pdbqt"

    # Extract chain G with surrounding shell (20 Å from pocket center)
    cx, cy, cz = compute_pocket_center(pdb_path)

    lines = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            dist = ((x-cx)**2 + (y-cy)**2 + (z-cz)**2) ** 0.5
            if dist < 20.0:
                lines.append(line)
    lines.append("END\n")

    with open(receptor_pdb, "w") as f:
        f.writelines(lines)
    print(f"  Receptor shell: {len(lines)-1} atoms within 20 Å of pocket")

    run(f"obabel {receptor_pdb} -O {receptor_pdbqt} -xr -h")
    if not receptor_pdbqt.exists():
        print("ERROR: Receptor PDBQT failed")
        sys.exit(1)

    return receptor_pdbqt, (cx, cy, cz)

def run_vina(receptor_pdbqt, ligand_pdbqt, center, exhaustiveness=32):
    """Dock BMS-626529 into the relay pocket."""
    cx, cy, cz = center
    config = DOCK / "vina_relay.conf"
    out = DOCK / "bms626529_docked.pdbqt"
    log = DOCK / "vina_relay.log"

    conf_text = f"""receptor = {receptor_pdbqt}
ligand = {ligand_pdbqt}
out = {out}

center_x = {cx:.3f}
center_y = {cy:.3f}
center_z = {cz:.3f}

size_x = 22.0
size_y = 22.0
size_z = 22.0

exhaustiveness = {exhaustiveness}
num_modes = 20
energy_range = 5.0
"""
    with open(config, "w") as f:
        f.write(conf_text)

    print(f"\n=== DOCKING: {BMS626529_NAME} into relay pocket ===")
    print(f"  Center: ({cx:.1f}, {cy:.1f}, {cz:.1f})")
    print(f"  Box: 22×22×22 Å")
    print(f"  Exhaustiveness: {exhaustiveness}")

    r = run(f"vina --config {config} --log {log}")
    if r.returncode != 0:
        print(f"  Vina failed: {r.stderr[:500]}")
        return None

    # Parse results from log
    results = parse_vina_log(log)
    return results

def run_control_docking(pdb_path, ligand_pdbqt):
    """Dock BMS-626529 at known fostemsavir binding site (beta20-beta21) as positive control."""
    # Fostemsavir binds at gp120 β20-β21 hairpin, near residue 375 (HXB2 Phe375)
    # In 6MEO numbering that's approximately G:370-380
    # Use PDB 5U7O coordinates as reference: center ~(125, 170, 142)
    # But we need to compute from 6MEO
    xs, ys, zs = [], [], []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chain = line[21]
            resnum = int(line[22:26].strip())
            if chain == "G" and 370 <= resnum <= 380:
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))

    if not xs:
        # Try wider range around attachment site
        with open(pdb_path) as f:
            for line in f:
                if not line.startswith("ATOM"):
                    continue
                chain = line[21]
                resnum = int(line[22:26].strip())
                if chain == "G" and 415 <= resnum <= 425:
                    xs.append(float(line[30:38]))
                    ys.append(float(line[38:46]))
                    zs.append(float(line[46:54]))

    if not xs:
        print("  WARNING: Could not find fostemsavir control site residues")
        return None

    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)
    print(f"\n=== CONTROL DOCKING: {BMS626529_NAME} at fostemsavir site ===")
    print(f"  Center: ({cx:.1f}, {cy:.1f}, {cz:.1f})")

    receptor_pdb = DOCK / "control_receptor.pdb"
    receptor_pdbqt = DOCK / "control_receptor.pdbqt"

    lines = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            dist = ((x-cx)**2 + (y-cy)**2 + (z-cz)**2) ** 0.5
            if dist < 20.0:
                lines.append(line)
    lines.append("END\n")

    with open(receptor_pdb, "w") as f:
        f.writelines(lines)

    run(f"obabel {receptor_pdb} -O {receptor_pdbqt} -xr -h")

    config = DOCK / "vina_control.conf"
    out = DOCK / "bms626529_control_docked.pdbqt"
    log = DOCK / "vina_control.log"

    conf_text = f"""receptor = {receptor_pdbqt}
ligand = {ligand_pdbqt}
out = {out}

center_x = {cx:.3f}
center_y = {cy:.3f}
center_z = {cz:.3f}

size_x = 22.0
size_y = 22.0
size_z = 22.0

exhaustiveness = 32
num_modes = 20
energy_range = 5.0
"""
    with open(config, "w") as f:
        f.write(conf_text)

    r = run(f"vina --config {config} --log {log}")
    if r.returncode != 0:
        print(f"  Control docking failed: {r.stderr[:500]}")
        return None

    return parse_vina_log(log)

def parse_vina_log(log_path):
    """Extract binding affinities from Vina log."""
    results = []
    in_results = False
    with open(log_path) as f:
        for line in f:
            if "-----+-----" in line:
                in_results = True
                continue
            if in_results and line.strip():
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        mode = int(parts[0])
                        affinity = float(parts[1])
                        rmsd_lb = float(parts[2])
                        rmsd_ub = float(parts[3])
                        results.append({
                            "mode": mode,
                            "affinity_kcal_mol": affinity,
                            "rmsd_lb": rmsd_lb,
                            "rmsd_ub": rmsd_ub
                        })
                    except ValueError:
                        if results:
                            break
    return results

def main():
    pdb_path = INPUT / "hiv_6meo.pdb"
    if not pdb_path.exists():
        print(f"ERROR: PDB not found at {pdb_path}")
        sys.exit(1)

    os.makedirs(DOCK, exist_ok=True)

    print("=" * 60)
    print("PHASE 1: MOLECULAR DOCKING")
    print("  BMS-626529 -> gp120 relay pocket (G:465-472)")
    print("=" * 60)

    # Prepare ligand
    print("\n--- Preparing ligand ---")
    lig_pdbqt = prepare_ligand()

    # Prepare receptor
    print("\n--- Preparing receptor ---")
    rec_pdbqt, center = prepare_receptor(pdb_path)

    # Dock into relay pocket
    relay_results = run_vina(rec_pdbqt, lig_pdbqt, center)

    # Control docking at fostemsavir site
    control_results = run_control_docking(pdb_path, lig_pdbqt)

    # Report
    print("\n" + "=" * 60)
    print("DOCKING RESULTS SUMMARY")
    print("=" * 60)

    if relay_results:
        print(f"\n  RELAY POCKET (G:465-472):")
        for r in relay_results[:5]:
            print(f"    Mode {r['mode']}: {r['affinity_kcal_mol']:.1f} kcal/mol  (RMSD lb={r['rmsd_lb']:.1f} ub={r['rmsd_ub']:.1f})")
        best = relay_results[0]["affinity_kcal_mol"]
        print(f"  BEST AFFINITY: {best:.1f} kcal/mol")
        if best <= -6.0:
            print(f"  -> DRUGGABLE (threshold: ≤ -6.0 kcal/mol)")
        else:
            print(f"  -> Below druggability threshold")
    else:
        print("\n  RELAY POCKET: Docking failed")

    if control_results:
        print(f"\n  FOSTEMSAVIR SITE (control):")
        for r in control_results[:5]:
            print(f"    Mode {r['mode']}: {r['affinity_kcal_mol']:.1f} kcal/mol  (RMSD lb={r['rmsd_lb']:.1f} ub={r['rmsd_ub']:.1f})")
        best_ctrl = control_results[0]["affinity_kcal_mol"]
        print(f"  BEST AFFINITY: {best_ctrl:.1f} kcal/mol")
    else:
        print("\n  FOSTEMSAVIR SITE: Control docking failed")

    # Save results
    output = {
        "ligand": BMS626529_NAME,
        "relay_pocket": relay_results[:10] if relay_results else [],
        "fostemsavir_control": control_results[:10] if control_results else [],
    }
    results_path = DOCK / "docking_results.json"
    with open(results_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\n  Results saved to {results_path}")

if __name__ == "__main__":
    main()
