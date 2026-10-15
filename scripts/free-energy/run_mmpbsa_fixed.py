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
Fix MM-PBSA atom count mismatch and run calculation.

Problem: tleap adds 37 terminal H at fragment cut points (1715 -> 1752 atoms).
Solution: Strip extra atoms from AMBER topology to match trajectory atom count.
"""
import parmed
import subprocess
import os
import json

WORKDIR = os.path.join(DATA_ROOT, "results", "08-binding-energetics")
MDDIR = os.path.join(DATA_ROOT, "topology")
os.chdir(WORKDIR)

# Step 1: Load the AMBER topology that gmx_MMPBSA already built
print("Step 1: Analyzing atom mismatch...")
parm = parmed.load_file('COM.prmtop', '_GMXMMPBSA_COM.inpcrd')
print(f"  AMBER topology: {len(parm.atoms)} atoms")

# Load the reference PDB (matches trajectory atom count)
ref = parmed.load_file('_GMXMMPBSA_COM.pdb')
print(f"  GROMACS PDB:     {len(ref.atoms)} atoms")
print(f"  Difference:      {len(parm.atoms) - len(ref.atoms)} atoms")

# Step 2: Find which atoms in parm don't exist in ref
# Map by residue number + atom name
ref_atoms = set()
for a in ref.atoms:
    key = (a.residue.number, a.residue.name, a.name)
    ref_atoms.add(key)

extra_indices = []
for i, a in enumerate(parm.atoms):
    key = (a.residue.number, a.residue.name, a.name)
    if key not in ref_atoms:
        extra_indices.append(i + 1)  # 1-indexed for AMBER mask

print(f"\nStep 2: Found {len(extra_indices)} extra atoms to strip")
if len(extra_indices) <= 50:
    for idx in extra_indices:
        a = parm.atoms[idx - 1]
        print(f"  atom {idx}: {a.residue.name}{a.residue.number}.{a.name}")

# Step 3: Strip extra atoms from topology
print(f"\nStep 3: Stripping {len(extra_indices)} atoms from topology...")
mask = ','.join(str(i) for i in extra_indices)
strip_mask = f'@{mask}'
parm.strip(strip_mask)
print(f"  Stripped topology: {len(parm.atoms)} atoms")

# Save stripped topologies
parm.save('COM_fixed.prmtop', overwrite=True)
print("  Saved COM_fixed.prmtop")

# Step 4: Do the same for receptor and ligand topologies
print("\nStep 4: Fixing receptor and ligand topologies...")

rec_parm = parmed.load_file('REC.prmtop', '_GMXMMPBSA_REC.inpcrd')
rec_ref = parmed.load_file('_GMXMMPBSA_REC.pdb')
print(f"  REC prmtop: {len(rec_parm.atoms)} -> target {len(rec_ref.atoms)}")

rec_ref_atoms = set()
for a in rec_ref.atoms:
    rec_ref_atoms.add((a.residue.number, a.residue.name, a.name))

rec_extra = []
for i, a in enumerate(rec_parm.atoms):
    if (a.residue.number, a.residue.name, a.name) not in rec_ref_atoms:
        rec_extra.append(i + 1)

if rec_extra:
    rec_parm.strip(f'@{",".join(str(i) for i in rec_extra)}')
rec_parm.save('REC_fixed.prmtop', overwrite=True)
print(f"  REC fixed: {len(rec_parm.atoms)} atoms")

lig_parm = parmed.load_file('LIG.prmtop', '_GMXMMPBSA_LIG.inpcrd')
lig_ref = parmed.load_file('_GMXMMPBSA_LIG.pdb')
print(f"  LIG prmtop: {len(lig_parm.atoms)} -> target {len(lig_ref.atoms)}")

lig_ref_atoms = set()
for a in lig_ref.atoms:
    lig_ref_atoms.add((a.residue.number, a.residue.name, a.name))

lig_extra = []
for i, a in enumerate(lig_parm.atoms):
    if (a.residue.number, a.residue.name, a.name) not in lig_ref_atoms:
        lig_extra.append(i + 1)

if lig_extra:
    lig_parm.strip(f'@{",".join(str(i) for i in lig_extra)}')
lig_parm.save('LIG_fixed.prmtop', overwrite=True)
print(f"  LIG fixed: {len(lig_parm.atoms)} atoms")

# Step 5: Verify counts
print(f"\nStep 5: Verification")
com_check = parmed.load_file('COM_fixed.prmtop')
rec_check = parmed.load_file('REC_fixed.prmtop')
lig_check = parmed.load_file('LIG_fixed.prmtop')
print(f"  COM: {len(com_check.atoms)} (trajectory: 1715)")
print(f"  REC: {len(rec_check.atoms)} (expected: 1322)")
print(f"  LIG: {len(lig_check.atoms)} (expected: 393)")
print(f"  REC + LIG = {len(rec_check.atoms) + len(lig_check.atoms)}")

if len(com_check.atoms) == 1715:
    print("\n  ATOM COUNTS MATCH! Ready to run MMPBSA.py")
else:
    print(f"\n  WARNING: COM has {len(com_check.atoms)}, expected 1715")
    print("  Will need further adjustment")
