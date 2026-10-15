#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Dock BMS-626529 into PDB 6MET relay pocket -- cross-structure validation."""
import subprocess, os, sys, json
from pathlib import Path

WORK = Path.home() / "md-relay-pocket"
INPUT = WORK / "input"
DOCK = WORK / "docking_6met"
os.makedirs(DOCK, exist_ok=True)

RELAY_RESIDUES = list(range(465, 473))

def run(cmd):
    print("  >> " + cmd)
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0 and r.stderr:
        print("  STDERR: " + r.stderr[:300])
    return r

def pocket_center(pdb):
    xs, ys, zs = [], [], []
    with open(pdb) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[21] == "G" and int(line[22:26].strip()) in RELAY_RESIDUES:
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
    cx = sum(xs) / len(xs)
    cy = sum(ys) / len(ys)
    cz = sum(zs) / len(zs)
    print("  6MET pocket center: ({:.2f}, {:.2f}, {:.2f}) from {} atoms".format(cx, cy, cz, len(xs)))
    return cx, cy, cz

def extract_shell(pdb, cx, cy, cz, radius=20.0):
    lines = []
    with open(pdb) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            if ((x - cx)**2 + (y - cy)**2 + (z - cz)**2)**0.5 < radius:
                lines.append(line)
    lines.append("END\n")
    return lines

def parse_log(logfile):
    results = []
    in_results = False
    with open(logfile) as f:
        for line in f:
            if "-----+-----" in line:
                in_results = True
                continue
            if in_results and line.strip():
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        results.append({
                            "mode": int(parts[0]),
                            "affinity_kcal_mol": float(parts[1]),
                            "rmsd_lb": float(parts[2]),
                            "rmsd_ub": float(parts[3])
                        })
                    except ValueError:
                        if results:
                            break
    return results

def dock(receptor_pdbqt, ligand_pdbqt, cx, cy, cz, out_pdbqt, log_file, label):
    print("\n=== DOCKING: {} ===".format(label))
    cmd = ("vina --receptor {} --ligand {} --out {} "
           "--center_x {:.3f} --center_y {:.3f} --center_z {:.3f} "
           "--size_x 22 --size_y 22 --size_z 22 "
           "--exhaustiveness 32 --num_modes 20 2>&1 | tee {}").format(
        receptor_pdbqt, ligand_pdbqt, out_pdbqt, cx, cy, cz, log_file)
    r = run(cmd)
    if r.stdout:
        lines = r.stdout.strip().split("\n")
        for l in lines[-8:]:
            print("  " + l)
    return parse_log(str(log_file))

def main():
    pdb = INPUT / "6MET.pdb"
    if not pdb.exists():
        print("ERROR: 6MET.pdb not found at " + str(pdb))
        sys.exit(1)

    # Find existing ligand PDBQT
    lig = WORK / "docking" / "bms626529_lig.pdbqt"
    if not lig.exists():
        lig = WORK / "docking" / "bms626529.pdbqt"
    if not lig.exists():
        print("ERROR: No ligand PDBQT found")
        sys.exit(1)
    print("  Using ligand: " + str(lig))

    # === RELAY POCKET ===
    cx, cy, cz = pocket_center(pdb)

    shell = extract_shell(pdb, cx, cy, cz)
    rec_pdb = DOCK / "6met_relay_receptor.pdb"
    rec_pdbqt = DOCK / "6met_relay_receptor.pdbqt"
    with open(rec_pdb, "w") as f:
        f.writelines(shell)
    print("  Receptor shell: {} atoms".format(len(shell) - 1))

    run("obabel {} -O {} -xr -h".format(rec_pdb, rec_pdbqt))

    relay_results = dock(
        rec_pdbqt, lig, cx, cy, cz,
        DOCK / "bms626529_6met_relay_docked.pdbqt",
        DOCK / "vina_6met_relay.log",
        "BMS-626529 into 6MET relay pocket (G:465-472)")

    # === CONTROL SITE (attachment site G:415-425) ===
    xs2, ys2, zs2 = [], [], []
    with open(pdb) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[21] == "G" and 415 <= int(line[22:26].strip()) <= 425:
                xs2.append(float(line[30:38]))
                ys2.append(float(line[38:46]))
                zs2.append(float(line[46:54]))

    ctrl_results = []
    if xs2:
        cx2 = sum(xs2) / len(xs2)
        cy2 = sum(ys2) / len(ys2)
        cz2 = sum(zs2) / len(zs2)

        shell2 = extract_shell(pdb, cx2, cy2, cz2)
        ctrl_pdb = DOCK / "6met_control_receptor.pdb"
        ctrl_pdbqt = DOCK / "6met_control_receptor.pdbqt"
        with open(ctrl_pdb, "w") as f:
            f.writelines(shell2)

        run("obabel {} -O {} -xr -h".format(ctrl_pdb, ctrl_pdbqt))

        ctrl_results = dock(
            ctrl_pdbqt, lig, cx2, cy2, cz2,
            DOCK / "bms626529_6met_control_docked.pdbqt",
            DOCK / "vina_6met_control.log",
            "BMS-626529 at 6MET attachment site (control)")

    # === RESULTS ===
    output = {
        "experiment": "BMS-626529 cross-structure docking on PDB 6MET",
        "date": "2026-07-16",
        "software": "AutoDock Vina v1.2.7",
        "structure": "PDB 6MET (alternate conformation, same study as 6MEO)",
        "relay_pocket": {
            "target": "gp120 C5 relay pocket (G:465-472)",
            "center_angstrom": [round(cx, 2), round(cy, 2), round(cz, 2)],
            "results": relay_results[:5],
            "best_affinity_kcal_mol": relay_results[0]["affinity_kcal_mol"] if relay_results else None
        },
        "control_site": {
            "target": "gp120 attachment site (G:415-425)",
            "results": ctrl_results[:5],
            "best_affinity_kcal_mol": ctrl_results[0]["affinity_kcal_mol"] if ctrl_results else None
        },
        "cross_structure_comparison": {
            "structure_6meo_relay_kcal_mol": -6.692,
            "structure_6meo_control_kcal_mol": -5.160,
            "structure_6met_relay_kcal_mol": relay_results[0]["affinity_kcal_mol"] if relay_results else None,
            "structure_6met_control_kcal_mol": ctrl_results[0]["affinity_kcal_mol"] if ctrl_results else None
        }
    }

    if relay_results and ctrl_results:
        output["delta_kcal_mol"] = round(
            relay_results[0]["affinity_kcal_mol"] - ctrl_results[0]["affinity_kcal_mol"], 3)

    results_path = DOCK / "6met_docking_results.json"
    with open(results_path, "w") as f:
        json.dump(output, f, indent=2)

    print("\n" + "=" * 60)
    print("6MET CROSS-STRUCTURE DOCKING RESULTS")
    print("=" * 60)
    if relay_results:
        print("  6MET Relay pocket best:  {:.1f} kcal/mol".format(relay_results[0]["affinity_kcal_mol"]))
    if ctrl_results:
        print("  6MET Control site best:  {:.1f} kcal/mol".format(ctrl_results[0]["affinity_kcal_mol"]))
    if relay_results and ctrl_results:
        delta = relay_results[0]["affinity_kcal_mol"] - ctrl_results[0]["affinity_kcal_mol"]
        print("  6MET Delta:              {:.1f} kcal/mol".format(delta))
    print("")
    print("  6MEO Relay pocket best:  -6.7 kcal/mol")
    print("  6MEO Control site best:  -5.2 kcal/mol")
    print("  6MEO Delta:              -1.5 kcal/mol")
    print("")
    if relay_results:
        diff = abs(relay_results[0]["affinity_kcal_mol"] - (-6.692))
        if diff < 1.0:
            print("  CROSS-STRUCTURE CONCORDANCE: {:.1f} kcal/mol difference (GOOD)".format(diff))
        else:
            print("  CROSS-STRUCTURE DIFFERENCE: {:.1f} kcal/mol (investigate)".format(diff))
    print("\n  Results saved to " + str(results_path))

if __name__ == "__main__":
    main()
