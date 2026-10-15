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
Set up and run MM-PBSA per-residue decomposition on the 20ns ternary complex.
Receptor = gp120 (chain A, residues 1-398)
Ligand = CD4 (chain B, residues 1-176, atom-sequential 399-574)
CCR5 context (chain C) is part of the receptor context.

Uses gmx_MMPBSA from the mmpbsa conda environment.
"""
import os
import subprocess
import json
from pathlib import Path
from datetime import datetime

WORK_DIR = Path(os.path.join(DATA_ROOT))
MD_DIR = WORK_DIR / "md"
OUT_DIR = WORK_DIR / "mmpbsa"
OUT_DIR.mkdir(exist_ok=True)

def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def run(cmd, cwd=None):
    """Run a command and return output."""
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True,
        timeout=7200, cwd=str(cwd or MD_DIR)
    )
    if result.returncode != 0:
        log(f"  STDERR: {result.stderr[-500:]}")
    return result.stdout + result.stderr, result.returncode

log("=" * 70)
log("MM-PBSA PER-RESIDUE: 20ns Full Ternary Complex")
log("=" * 70)

# Step 1: Index file already created by make_chain_ndx.py
# Groups: 0=gp120 (6204 atoms), 1=CD4 (2769 atoms), 2=CCR5 (5176 atoms), 3=Receptor (gp120+CCR5)
log("Step 1: Using pre-built chain index file")
ndx_file = str(OUT_DIR / "mmpbsa.ndx")
if not os.path.exists(ndx_file):
    log("  ERROR: Index file not found! Run make_chain_ndx.py first.")
    sys.exit(1)
log(f"  Index file: {ndx_file}")

# Step 2: Create MMPBSA input file
log("Step 2: Writing MMPBSA input")

mmpbsa_input = """
&general
  startframe=1,
  endframe=100,
  interval=1,
  verbose=2,
  keep_files=0,
/
&gb
  igb=5,
  saltcon=0.150,
/
&decomp
  idecomp=2,
  dec_verbose=3,
  print_res="all"
/
"""

input_file = OUT_DIR / "mmpbsa.in"
with open(input_file, "w") as f:
    f.write(mmpbsa_input)

log(f"  Input file: {input_file}")

# Step 3: Extract a subset trajectory (every 10th frame for speed)
log("Step 3: Extracting trajectory subset")
subset_xtc = str(OUT_DIR / "md_subset.xtc")
output, rc = run(f"echo '0' | gmx trjconv -s md.tpr -f md.xtc -o {subset_xtc} -skip 10 -pbc mol")
log(f"  Trajectory subset: {rc}")

# Step 4: Run gmx_MMPBSA
log("Step 4: Running gmx_MMPBSA")
log("  This may take 30-60 minutes...")

# Index groups: 0=gp120, 1=CD4, 2=CCR5, 3=Receptor(gp120+CCR5)
# MM-PBSA: receptor = gp120 (0), ligand = CD4 (1)
# CCR5 is implicit context (part of the simulation but not receptor/ligand)

cmd = (
    f"gmx_MMPBSA -O "
    f"-i {input_file} "
    f"-cs md.tpr "
    f"-ci {ndx_file} "
    f"-cg 0 1 "
    f"-ct {subset_xtc} "
    f"-cp topol.top "
    f"-o {OUT_DIR}/FINAL_RESULTS.dat "
    f"-do {OUT_DIR}/FINAL_DECOMP.dat "
    f"-eo {OUT_DIR}/FINAL_ENERGY.csv "
    f"2>&1"
)

log(f"  Command: {cmd[:100]}...")
output, rc = run(cmd, cwd=MD_DIR)

# Save output
with open(OUT_DIR / "mmpbsa_output.log", "w") as f:
    f.write(output)

log(f"  gmx_MMPBSA returned: {rc}")

if rc == 0:
    log("  MM-PBSA COMPLETE!")
    # Parse results
    results_file = OUT_DIR / "FINAL_RESULTS.dat"
    if results_file.exists():
        with open(results_file) as f:
            content = f.read()
        log(f"  Results file: {len(content)} bytes")
        # Print key lines
        for line in content.split("\n"):
            if "DELTA" in line or "TOTAL" in line or "VDWAALS" in line or "EEL" in line:
                log(f"    {line.strip()}")
else:
    log("  MM-PBSA FAILED -- check mmpbsa_output.log")
    # Print last 50 lines of output
    lines = output.strip().split("\n")
    for line in lines[-50:]:
        log(f"    {line}")

log("")
log("=" * 70)
log("DONE")
log("=" * 70)
