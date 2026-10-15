import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Dock BMS-626529 into the ternary complex for IE-MMPBSA preparation.

Workflow:
1. Extract BMS-626529 ligand from 6MU6 PDB (already has the compound)
2. Extract the equilibrated ternary complex receptor from the MD trajectory
3. Prepare PDBQT files for AutoDock Vina
4. Run docking focused on the relay pocket region
5. Select the best pose for MD + IE-MMPBSA
"""
import os
import subprocess
import json
from pathlib import Path
from datetime import datetime

WORK_DIR = Path(os.path.join(DATA_ROOT))
DOCK_DIR = WORK_DIR / "docking"
DOCK_DIR.mkdir(exist_ok=True)
MD_DIR = WORK_DIR / "md"

def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def run(cmd, cwd=None, timeout=600):
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True,
        timeout=timeout, cwd=str(cwd or DOCK_DIR)
    )
    return result.stdout + result.stderr, result.returncode

log("=" * 70)
log("DOCKING: BMS-626529 into Ternary Complex (relay pocket)")
log("=" * 70)

# Step 1: Download 6MU6 PDB (contains BMS-814508/BMS-626529 bound to SOSIP)
log("Step 1: Getting 6MU6 structure for ligand extraction")

pdb_6mu6 = DOCK_DIR / "6mu6.pdb"
if not pdb_6mu6.exists():
    output, rc = run(f"wget -q -O {pdb_6mu6} 'https://files.rcsb.org/download/6MU6.pdb'")
    log(f"  Downloaded 6MU6: rc={rc}")
else:
    log(f"  6MU6 already downloaded")

# Step 2: Extract ligand from 6MU6
# BMS-626529 is listed as BMS or 4XU or similar in the PDB
log("Step 2: Identifying ligand in 6MU6")

# Find HETATM records
with open(pdb_6mu6) as f:
    lines = f.readlines()

hetatms = set()
for line in lines:
    if line.startswith("HETATM"):
        resname = line[17:20].strip()
        if resname not in ["HOH", "NAG", "BMA", "MAN", "FUC", "GAL", "SO4", "GOL", "EDO", "ACT", "PEG", "DMS"]:
            hetatms.add(resname)

log(f"  Non-standard residues in 6MU6: {sorted(hetatms)}")

# The ligand is BMS-814508 (a fostemsavir analog), likely named 9PR or similar
# Let's extract all non-sugar/non-solvent HETATM
ligand_resname = None
for resname in sorted(hetatms):
    count = sum(1 for l in lines if l.startswith("HETATM") and l[17:20].strip() == resname)
    log(f"  {resname}: {count} atoms")
    if count > 20:  # Drug-sized molecule
        ligand_resname = resname

if not ligand_resname:
    log("  WARNING: Could not identify ligand automatically")
    log("  Trying to use first large HETATM group")
    for resname in sorted(hetatms):
        count = sum(1 for l in lines if l.startswith("HETATM") and l[17:20].strip() == resname)
        if count >= 10:
            ligand_resname = resname
            break

log(f"  Using ligand: {ligand_resname}")

# Extract ligand to PDB
ligand_pdb = DOCK_DIR / "bms626529_from_6mu6.pdb"
with open(ligand_pdb, "w") as f:
    for line in lines:
        if line.startswith("HETATM") and line[17:20].strip() == ligand_resname:
            # Take only the first occurrence (first chain)
            chain = line[21]
            break
    
    for line in lines:
        if line.startswith("HETATM") and line[17:20].strip() == ligand_resname and line[21] == chain:
            f.write(line)
    f.write("END\n")

atom_count = sum(1 for l in open(ligand_pdb) if l.startswith("HETATM"))
log(f"  Extracted {atom_count} atoms from chain {chain}")

# Step 3: Convert ligand to PDBQT
log("Step 3: Converting ligand to PDBQT")
ligand_pdbqt = DOCK_DIR / "bms626529.pdbqt"
output, rc = run(f"obabel {ligand_pdb} -O {ligand_pdbqt} --gen3d -p 7.4")
log(f"  obabel conversion: rc={rc}")
if rc != 0:
    # Try without --gen3d since we already have 3D coords
    output, rc = run(f"obabel {ligand_pdb} -O {ligand_pdbqt} -p 7.4")
    log(f"  Retry without gen3d: rc={rc}")

# Step 4: Use pre-extracted receptor from MD trajectory (last frame, equilibrated)
log("Step 4: Using pre-extracted receptor PDB")
receptor_pdb = DOCK_DIR / "receptor_equilibrated.pdb"
if not receptor_pdb.exists():
    log("  ERROR: receptor_equilibrated.pdb not found in docking dir")
    log("  Extract with: echo '1' | gmx trjconv -s md.tpr -f md.xtc -o receptor_equilibrated.pdb -dump 20000 -pbc mol")
    exit(1)
atom_count_rec = sum(1 for l in open(receptor_pdb) if l.startswith("ATOM"))
log(f"  Receptor: {atom_count_rec} atoms")

# Convert receptor to PDBQT
receptor_pdbqt = DOCK_DIR / "receptor.pdbqt"
output, rc = run(f"obabel {receptor_pdb} -O {receptor_pdbqt} -xr")
log(f"  Receptor PDBQT: rc={rc}")

# Step 5: Determine docking box centered on relay pocket (res 372-379)
log("Step 5: Computing docking box from relay pocket coordinates")

# Find relay pocket center from receptor PDB
pocket_coords = []
with open(receptor_pdb) as f:
    for line in f:
        if line.startswith("ATOM"):
            resnum = int(line[22:26].strip())
            atomname = line[12:16].strip()
            if 372 <= resnum <= 379 and atomname == "CA":
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                pocket_coords.append((x, y, z))

if pocket_coords:
    import numpy as np
    coords = np.array(pocket_coords)
    center = coords.mean(axis=0)
    span = coords.max(axis=0) - coords.min(axis=0)
    box_size = span + 15.0  # 15 Å padding
    
    log(f"  Pocket center: ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f})")
    log(f"  Pocket span: ({span[0]:.1f}, {span[1]:.1f}, {span[2]:.1f})")
    log(f"  Box size: ({box_size[0]:.1f}, {box_size[1]:.1f}, {box_size[2]:.1f})")
else:
    log("  WARNING: Could not find pocket CA atoms in receptor")
    log("  Using default box center")
    center = [0, 0, 0]
    box_size = [30, 30, 30]

# Step 6: Run AutoDock Vina
log("Step 6: Running AutoDock Vina docking")

vina_conf = DOCK_DIR / "vina_config.txt"
with open(vina_conf, "w") as f:
    f.write(f"receptor = {receptor_pdbqt}\n")
    f.write(f"ligand = {ligand_pdbqt}\n")
    f.write(f"center_x = {center[0]:.3f}\n")
    f.write(f"center_y = {center[1]:.3f}\n")
    f.write(f"center_z = {center[2]:.3f}\n")
    f.write(f"size_x = {box_size[0]:.1f}\n")
    f.write(f"size_y = {box_size[1]:.1f}\n")
    f.write(f"size_z = {box_size[2]:.1f}\n")
    f.write("exhaustiveness = 32\n")
    f.write("num_modes = 20\n")
    f.write("energy_range = 5\n")

output_pdbqt = DOCK_DIR / "docked_poses.pdbqt"
output, rc = run(f"vina --config {vina_conf} --out {output_pdbqt} 2>&1", timeout=3600)

with open(DOCK_DIR / "vina_output.log", "w") as f:
    f.write(output)

log(f"  Vina returned: rc={rc}")

# Parse docking results
poses = []
for line in output.split("\n"):
    line = line.strip()
    if line and line[0].isdigit() and "affinity" not in line.lower():
        parts = line.split()
        if len(parts) >= 4:
            try:
                mode = int(parts[0])
                affinity = float(parts[1])
                poses.append({"mode": mode, "affinity_kcal": affinity})
                log(f"    Mode {mode}: {affinity:.1f} kcal/mol")
            except ValueError:
                pass

# Save results
results = {
    "ligand_source": "6MU6",
    "ligand_resname": ligand_resname,
    "receptor": "20ns equilibrated ternary complex (last frame)",
    "pocket_center": [float(c) for c in center] if pocket_coords else None,
    "box_size": [float(s) for s in box_size] if pocket_coords else None,
    "exhaustiveness": 32,
    "poses": poses,
}

with open(DOCK_DIR / "docking_results.json", "w") as f:
    json.dump(results, f, indent=2)

log(f"\nResults saved to {DOCK_DIR / 'docking_results.json'}")
log("=" * 70)
log("DONE")
