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
Run MM-PBSA with idecomp=3 (pairwise) on the 20ns ternary complex.
idecomp=3 computes pairwise interaction energies between residues,
giving the actual receptor-ligand pairwise binding contributions.
"""
import os
import subprocess
import json
from pathlib import Path
from datetime import datetime

WORK_DIR = Path(os.path.join(DATA_ROOT))
MD_DIR = WORK_DIR / "md"
OUT_DIR = WORK_DIR / "mmpbsa_idecomp3"
OUT_DIR.mkdir(exist_ok=True)

def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def run(cmd, cwd=None, timeout=14400):
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True,
        timeout=timeout, cwd=str(cwd or MD_DIR)
    )
    if result.returncode != 0:
        log(f"  STDERR: {result.stderr[-500:]}")
    return result.stdout + result.stderr, result.returncode

log("=" * 70)
log("MM-PBSA idecomp=3 (PAIRWISE): 20ns Full Ternary Complex")
log("=" * 70)

ndx_file = str(WORK_DIR / "mmpbsa" / "mmpbsa.ndx")
subset_xtc = str(WORK_DIR / "mmpbsa" / "md_subset.xtc")

# Step 1: Write idecomp=3 input (pairwise decomposition)
log("Step 1: Writing idecomp=3 input")
mmpbsa_input = """
&general
  startframe=1,
  endframe=100,
  interval=5,
  verbose=2,
  keep_files=0,
/
&gb
  igb=5,
  saltcon=0.150,
/
&decomp
  idecomp=3,
  dec_verbose=3,
  print_res="all"
/
"""

input_file = OUT_DIR / "mmpbsa_idecomp3.in"
with open(input_file, "w") as f:
    f.write(mmpbsa_input)

# Step 2: Run gmx_MMPBSA
log("Step 2: Running gmx_MMPBSA (idecomp=3 pairwise)")
log("  NOTE: Pairwise decomposition on large systems can be slow (~60-120 min)")

cmd = (
    f"gmx_MMPBSA -O "
    f"-i {input_file} "
    f"-cs md.tpr "
    f"-ci {ndx_file} "
    f"-cg 0 1 "
    f"-ct {subset_xtc} "
    f"-cp topol.top "
    f"-o {OUT_DIR}/FINAL_RESULTS_id3.dat "
    f"-do {OUT_DIR}/FINAL_DECOMP_id3.dat "
    f"-eo {OUT_DIR}/FINAL_ENERGY_id3.csv "
    f"2>&1"
)

log(f"  Command: {cmd[:100]}...")
output, rc = run(cmd, cwd=MD_DIR)

with open(OUT_DIR / "mmpbsa_idecomp3_output.log", "w") as f:
    f.write(output)

if "Timing:" in output:
    log("  CALCULATION COMPLETE (gmx_MMPBSA_ana failure is expected)")
    for line in output.split("\n"):
        if "Total time" in line or "Total GB" in line:
            log(f"    {line.strip()}")
    
    decomp_file = OUT_DIR / "FINAL_DECOMP_id3.dat"
    if decomp_file.exists():
        size = decomp_file.stat().st_size
        log(f"  Pairwise decomp file: {size} bytes")
else:
    log(f"  gmx_MMPBSA returned: {rc}")
    lines = output.strip().split("\n")
    for line in lines[-30:]:
        log(f"    {line}")

log("")
log("=" * 70)
log("DONE")
