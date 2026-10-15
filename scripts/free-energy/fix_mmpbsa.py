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
Fix MM-PBSA atom count mismatch for fragmented system.

Strategy: Use parmed to convert GROMACS topology directly to AMBER format,
avoiding tleap's terminal hydrogen additions. Then use MMPBSA.py directly.
"""
import subprocess
import os
import sys

WORKDIR = os.path.join(DATA_ROOT, "results", "08-binding-energetics")
MDDIR = os.path.join(DATA_ROOT, "topology")
os.chdir(WORKDIR)

# Step 1: Extract protein-only trajectory (no water/ions)
print("Step 1: Extract protein-only trajectory...")
proc = subprocess.run(
    ['gmx', 'trjconv', '-s', f'{MDDIR}/md.tpr', '-f', f'{MDDIR}/md.xtc',
     '-o', 'protein_only.xtc', '-n', f'{MDDIR}/analysis.ndx'],
    input='Protein\n', capture_output=True, text=True
)
if proc.returncode != 0:
    # Try with default Protein group
    proc = subprocess.run(
        ['gmx', 'trjconv', '-s', f'{MDDIR}/md.tpr', '-f', f'{MDDIR}/md.xtc',
         '-o', 'protein_only.xtc'],
        input='1\n', capture_output=True, text=True
    )
print(f"  trjconv exit: {proc.returncode}")

# Step 2: Extract first frame as PDB
print("Step 2: Extract reference PDB...")
proc = subprocess.run(
    ['gmx', 'trjconv', '-s', f'{MDDIR}/md.tpr', '-f', f'{MDDIR}/md.xtc',
     '-o', 'ref_protein.pdb', '-dump', '0'],
    input='1\n', capture_output=True, text=True
)
print(f"  ref PDB exit: {proc.returncode}")

# Step 3: Check atom counts
with open('ref_protein.pdb') as f:
    natom_pdb = sum(1 for line in f if line.startswith('ATOM'))
print(f"  PDB atoms: {natom_pdb}")

# Step 4: Use pdb4amber to clean the PDB (adds missing H, fixes naming)
print("Step 3: pdb4amber processing...")
proc = subprocess.run(
    ['pdb4amber', '-i', 'ref_protein.pdb', '-o', 'amber_ready.pdb',
     '--nohyd', '--dry'],
    capture_output=True, text=True
)
print(f"  pdb4amber exit: {proc.returncode}")
if proc.stderr:
    print(f"  stderr: {proc.stderr[:500]}")

with open('amber_ready.pdb') as f:
    natom_amber = sum(1 for line in f if line.startswith('ATOM'))
print(f"  amber_ready atoms: {natom_amber}")

# Step 5: Build AMBER topology with tleap using --nohyd PDB (no extra H)
print("Step 4: Build AMBER topology with tleap (no extra H)...")
tleap_input = """source leaprc.protein.ff14SB
source leaprc.water.tip3p
com = loadpdb amber_ready.pdb
set default PBRadii mbondi2
saveamberparm com com_noh.prmtop com_noh.inpcrd
quit
"""
with open('tleap_noh.in', 'w') as f:
    f.write(tleap_input)

proc = subprocess.run(['tleap', '-f', 'tleap_noh.in'],
                      capture_output=True, text=True)
print(f"  tleap exit: {proc.returncode}")

# Check topology
import parmed
try:
    t = parmed.load_file('com_noh.prmtop')
    print(f"  prmtop atoms: {len(t.atoms)}")
except Exception as e:
    print(f"  prmtop error: {e}")

print()
print("Done. Check atom counts above for match.")
