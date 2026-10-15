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
Prepare ligand-bound ternary complex for short MD + IE-MMPBSA.

Steps:
1. Generate GAFF parameters for JYV ligand using acpype/antechamber
2. Combine with receptor topology
3. Solvate and add ions
4. Energy minimize
5. Equilibrate (short NVT + NPT)
6. Production MD (5-10 ns)
"""
import os
import subprocess
from pathlib import Path
from datetime import datetime

WORK_DIR = Path(os.path.join(DATA_ROOT))
DOCK_DIR = WORK_DIR / "docking"
MD_LIG_DIR = WORK_DIR / "md-ligand-bound"
MD_LIG_DIR.mkdir(exist_ok=True)

def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def run(cmd, cwd=None, timeout=3600):
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True,
        timeout=timeout, cwd=str(cwd or MD_LIG_DIR)
    )
    if result.returncode != 0:
        log(f"  STDERR: {result.stderr[-500:]}")
    return result.stdout + result.stderr, result.returncode

log("=" * 70)
log("PREPARE: Ligand-bound ternary complex for MD + IE-MMPBSA")
log("=" * 70)

# Step 1: Check for acpype/antechamber
log("Step 1: Checking for ligand parameterization tools")
for tool in ["acpype", "antechamber", "obabel", "gmx"]:
    output, rc = run(f"which {tool}")
    status = "FOUND" if rc == 0 else "NOT FOUND"
    log(f"  {tool}: {status} ({output.strip()})")

# Step 2: Generate ligand parameters using obabel + acpype
log("Step 2: Generating ligand GAFF parameters")
ligand_pdb = DOCK_DIR / "bms626529_positioned.pdb"

if not ligand_pdb.exists():
    log(f"  ERROR: Ligand PDB not found: {ligand_pdb}")
    exit(1)

# Convert to mol2 with proper atom types
ligand_mol2 = MD_LIG_DIR / "jyv.mol2"
output, rc = run(f"obabel {ligand_pdb} -O {ligand_mol2} -p 7.4")
log(f"  obabel PDB -> MOL2: rc={rc}")

# Try acpype for GROMACS-ready topology
ligand_acpype_dir = MD_LIG_DIR / "jyv_acpype"
output, rc = run(f"acpype -i {ligand_mol2} -c bcc -n 0 -o gmx 2>&1", timeout=3600)
log(f"  acpype: rc={rc}")

if rc != 0:
    log("  acpype failed. Trying antechamber directly...")
    # Fall back to antechamber + parmchk2
    output, rc = run(f"antechamber -i {ligand_mol2} -fi mol2 -o {MD_LIG_DIR}/jyv_gaff.mol2 -fo mol2 -c bcc -nc 0 -at gaff2 2>&1")
    log(f"  antechamber: rc={rc}")
    if rc == 0:
        output, rc = run(f"parmchk2 -i {MD_LIG_DIR}/jyv_gaff.mol2 -f mol2 -o {MD_LIG_DIR}/jyv.frcmod -s gaff2")
        log(f"  parmchk2: rc={rc}")

# Check what was generated
log("Step 3: Checking generated files")
for pattern in ["jyv*", "JYV*"]:
    output, rc = run(f"ls -la {MD_LIG_DIR}/{pattern} 2>/dev/null; ls -la {MD_LIG_DIR}/*/{pattern} 2>/dev/null")
    if output.strip():
        for line in output.strip().split("\n"):
            if line.strip():
                log(f"  {line.strip()}")

# Try to find the GROMACS topology files from acpype
acpype_outputs = list(MD_LIG_DIR.glob("**/JYV_GMX*.itp")) + list(MD_LIG_DIR.glob("**/jyv_GMX*.itp"))
if acpype_outputs:
    log(f"  Found acpype GROMACS files: {[str(f) for f in acpype_outputs]}")
else:
    log("  No GROMACS topology from acpype found")
    # List the acpype output directory
    output, rc = run(f"find {MD_LIG_DIR} -name '*.itp' -o -name '*.gro' -o -name '*.top' 2>/dev/null")
    if output.strip():
        for line in output.strip().split("\n"):
            log(f"  {line}")

log("")
log("=" * 70)
log("PARAMETERIZATION COMPLETE -- next: combine topologies and solvate")
log("=" * 70)
