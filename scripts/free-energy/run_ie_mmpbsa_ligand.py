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
IE-MMPBSA for BMS-626529 binding to relay pocket in full ternary complex.

Task 1.3: Compute binding free energy of BMS-626529 (UNL) to the
gp120+CD4+CCR5 ternary complex using Interaction Entropy MM-GBSA.

Receptor: full protein (gp120 + CD4 + CCR5, chains A+B+C)
Ligand:   BMS-626529 (UNL, residue name UNL)

CRITICAL: The MD trajectory must be PBC-unwrapped before feeding to
gmx_MMPBSA. Without this, bonded atoms wrapped across periodic boundaries
produce ~160 A distances instead of ~1.2 A, causing sander BOND energy
overflow. This script handles the unwrapping automatically.

Requirements:
  - conda activate mmpbsa (gmx_MMPBSA 1.6.5, AmberTools 26)
  - GROMACS trajectory (md.xtc) + topology (topol.top) + run input (md.tpr)

Usage:
  python3 run_ie_mmpbsa_ligand.py <workdir>

Example:
  python3 run_ie_mmpbsa_ligand.py <path_to_md_workdir>
"""

import sys
import os
import subprocess
import json

workdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA_ROOT, "results", "04-md-simulations")
os.chdir(workdir)
print(f"Working directory: {workdir}")


def run_cmd(cmd_str, desc=""):
    """Run a shell command and return the result."""
    if desc:
        print(f"\n--- {desc} ---")
    print(f"  $ {cmd_str}")
    r = subprocess.run(['bash', '-c', cmd_str], capture_output=True, text=True)
    return r


# Step 1: Create index groups
print("\n=== Step 1: Create index groups ===")
ndx_input = "\nq\n"
result = subprocess.run(
    ["gmx", "make_ndx", "-f", "md.tpr", "-o", "mmpbsa.ndx"],
    input=ndx_input, capture_output=True, text=True
)

lines = result.stdout.split('\n')
protein_grp = None
unl_grp = None
for line in lines:
    line = line.strip()
    if 'Protein' in line and 'Protein-H' not in line and 'non-Protein' not in line:
        parts = line.split()
        if parts and parts[0].isdigit():
            protein_grp = parts[0]
    if 'UNL' in line:
        parts = line.split()
        if parts and parts[0].isdigit():
            unl_grp = parts[0]

print(f"Protein group: {protein_grp}")
print(f"UNL group: {unl_grp}")

if not protein_grp or not unl_grp:
    print("ERROR: Could not find Protein or UNL groups!")
    for line in lines:
        if line.strip() and line.strip()[0].isdigit():
            print(f"  {line.strip()}")
    sys.exit(1)

# Step 2: PBC unwrapping (CRITICAL)
print("\n=== Step 2: Unwrap PBC in trajectory ===")
print("Making molecules whole and centering on protein...")

unwrapped_xtc = "md_pbc_whole.xtc"
if not os.path.exists(unwrapped_xtc):
    run_cmd(
        f'echo -e "{protein_grp}\\n0" | gmx trjconv -f md.xtc -s md.tpr '
        f'-o {unwrapped_xtc} -pbc mol -center',
        "gmx trjconv -pbc mol -center"
    )

if not os.path.exists(unwrapped_xtc):
    print("ERROR: PBC unwrapping failed!")
    sys.exit(1)

size_mb = os.path.getsize(unwrapped_xtc) / 1e6
print(f"  {unwrapped_xtc}: {size_mb:.1f} MB")

# Step 3: Write gmx_MMPBSA input files
print("\n=== Step 3: Write MMPBSA input files ===")
mmpbsa_in = """&general
startframe=1, endframe=200, interval=5,
verbose=2,
/
&gb
igb=5, saltcon=0.15,
/
&decomp
idecomp=1, dec_verbose=0,
print_res="within 6"
/
"""
with open("mmpbsa.in", "w") as f:
    f.write(mmpbsa_in)
print("Wrote mmpbsa.in (igb=5, saltcon=0.15, 40 frames, per-residue decomp within 6 A)")

ie_in = """&general
startframe=1, endframe=1000, interval=5,
verbose=2,
interaction_entropy=1, ie_segment=25,
/
&gb
igb=5, saltcon=0.15,
/
"""
with open("mmpbsa_ie.in", "w") as f:
    f.write(ie_in)
print("Wrote mmpbsa_ie.in (interaction entropy, 200 frames, ie_segment=25)")

# Step 4: Run standard MM-GBSA
print("\n=== Step 4: Launch MM-GBSA (40 frames, per-residue decomp) ===")
cmd = [
    "gmx_MMPBSA", "-O",
    "-i", "mmpbsa.in",
    "-cs", "md.tpr",
    "-ct", unwrapped_xtc,
    "-ci", "mmpbsa.ndx",
    "-cg", protein_grp, unl_grp,
    "-cp", "topol.top",
    "-o", "FINAL_RESULTS_LIGAND.dat",
    "-do", "FINAL_DECOMP_LIGAND.dat",
    "-eo", "FINAL_RESULTS_LIGAND.csv",
    "-nogui",
]

print(f"Command: {' '.join(cmd)}")
result = subprocess.run(cmd, capture_output=True, text=True)
print("STDOUT:", result.stdout[-2000:] if len(result.stdout) > 2000 else result.stdout)
if result.returncode != 0:
    print("STDERR:", result.stderr[-2000:] if len(result.stderr) > 2000 else result.stderr)
    print(f"\nERROR: gmx_MMPBSA failed with return code {result.returncode}")
else:
    print("\nMM-GBSA completed successfully!")

# Step 5: Run IE-MMPBSA (interaction entropy)
print("\n=== Step 5: Launch IE-MMPBSA (interaction entropy) ===")
cmd_ie = [
    "gmx_MMPBSA", "-O",
    "-i", "mmpbsa_ie.in",
    "-cs", "md.tpr",
    "-ct", unwrapped_xtc,
    "-ci", "mmpbsa.ndx",
    "-cg", protein_grp, unl_grp,
    "-cp", "topol.top",
    "-o", "FINAL_RESULTS_IE.dat",
    "-eo", "FINAL_RESULTS_IE.csv",
    "-nogui",
]

print(f"Command: {' '.join(cmd_ie)}")
result_ie = subprocess.run(cmd_ie, capture_output=True, text=True)
print("STDOUT:", result_ie.stdout[-2000:] if len(result_ie.stdout) > 2000 else result_ie.stdout)
if result_ie.returncode != 0:
    print("STDERR:", result_ie.stderr[-2000:] if len(result_ie.stderr) > 2000 else result_ie.stderr)
    print(f"\nERROR: IE-MMPBSA failed with return code {result_ie.returncode}")
else:
    print("\nIE-MMPBSA completed successfully!")

# Step 6: Collect results as JSON
print("\n=== Step 6: Collect results ===")
results = {"method": "IE-MMPBSA", "status": "complete", "files": {}}
for f_name in ["FINAL_RESULTS_LIGAND.dat", "FINAL_DECOMP_LIGAND.dat",
               "FINAL_RESULTS_LIGAND.csv", "FINAL_RESULTS_IE.dat", "FINAL_RESULTS_IE.csv"]:
    if os.path.exists(f_name):
        size = os.path.getsize(f_name)
        results["files"][f_name] = f"{size:,} bytes"
        print(f"  {f_name}: {size:,} bytes")
    else:
        results["files"][f_name] = "NOT FOUND"
        print(f"  {f_name}: NOT FOUND")

with open("ie_mmpbsa_summary.json", "w") as f:
    json.dump(results, f, indent=2)

print("\nDone.")
