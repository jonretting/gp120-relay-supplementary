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
Fix MM-PBSA atom count mismatch using coordinate-based matching.

tleap renames residues and adds terminal H. Match atoms by 3D coordinates
to identify which AMBER atoms correspond to GROMACS atoms, then strip the rest.
"""
import parmed
import numpy as np
import subprocess
import os

WORKDIR = os.path.join(DATA_ROOT, "results", "08-binding-energetics")
os.chdir(WORKDIR)

print("=== Loading structures ===")
parm = parmed.load_file('COM.prmtop', '_GMXMMPBSA_COM.inpcrd')
ref = parmed.load_file('_GMXMMPBSA_COM.pdb')

print(f"AMBER topology: {len(parm.atoms)} atoms, {len(parm.residues)} residues")
print(f"GROMACS PDB:     {len(ref.atoms)} atoms, {len(ref.residues)} residues")

# Print first few atoms from each to see naming differences
print("\n=== First 10 atoms comparison ===")
print(f"{'AMBER':>40s}  |  {'GROMACS':<40s}")
for i in range(min(20, len(ref.atoms))):
    pa = parm.atoms[i] if i < len(parm.atoms) else None
    ra = ref.atoms[i] if i < len(ref.atoms) else None
    pa_str = f"{pa.residue.name}{pa.residue.number}.{pa.name}" if pa else "---"
    ra_str = f"{ra.residue.name}{ra.residue.number}.{ra.name}" if ra else "---"
    match = "OK" if pa and ra and abs(pa.xx - ra.xx) < 0.01 else "DIFF"
    print(f"  {i:4d}: {pa_str:>35s}  |  {ra_str:<35s}  {match}")

# Find where the extra atoms appear
print("\n=== Finding extra atoms by coordinate matching ===")

# Build coordinate arrays
amber_coords = np.array([(a.xx, a.xy, a.xz) for a in parm.atoms])
gro_coords = np.array([(a.xx, a.xy, a.xz) for a in ref.atoms])

# For each AMBER atom, find nearest GROMACS atom
# Extra atoms added by tleap won't have a close match
matched_amber = []
unmatched_amber = []

gro_used = set()
for i, ac in enumerate(amber_coords):
    dists = np.linalg.norm(gro_coords - ac, axis=1)
    nearest = np.argmin(dists)
    if dists[nearest] < 0.05 and nearest not in gro_used:
        matched_amber.append(i)
        gro_used.add(nearest)
    else:
        unmatched_amber.append(i)

print(f"Matched:   {len(matched_amber)} AMBER atoms have GROMACS counterparts")
print(f"Unmatched: {len(unmatched_amber)} AMBER atoms are tleap additions")

if len(unmatched_amber) <= 50:
    print("\nExtra atoms added by tleap:")
    for idx in unmatched_amber:
        a = parm.atoms[idx]
        print(f"  atom {idx+1}: {a.residue.name}{a.residue.number}.{a.name}")

# Strip unmatched atoms
print(f"\n=== Stripping {len(unmatched_amber)} extra atoms ===")
if unmatched_amber:
    mask = '@' + ','.join(str(i+1) for i in unmatched_amber)
    parm_stripped = parmed.load_file('COM.prmtop', '_GMXMMPBSA_COM.inpcrd')
    parm_stripped.strip(mask)
    parm_stripped.save('COM_fixed.prmtop', overwrite=True)
    print(f"  COM_fixed.prmtop: {len(parm_stripped.atoms)} atoms")

    # Same for receptor
    rec_parm = parmed.load_file('REC.prmtop', '_GMXMMPBSA_REC.inpcrd')
    rec_ref = parmed.load_file('_GMXMMPBSA_REC.pdb')
    rec_amber = np.array([(a.xx, a.xy, a.xz) for a in rec_parm.atoms])
    rec_gro = np.array([(a.xx, a.xy, a.xz) for a in rec_ref.atoms])
    rec_unmatched = []
    rec_used = set()
    for i, ac in enumerate(rec_amber):
        dists = np.linalg.norm(rec_gro - ac, axis=1)
        nearest = np.argmin(dists)
        if dists[nearest] < 0.05 and nearest not in rec_used:
            rec_used.add(nearest)
        else:
            rec_unmatched.append(i)
    if rec_unmatched:
        rec_parm.strip('@' + ','.join(str(i+1) for i in rec_unmatched))
    rec_parm.save('REC_fixed.prmtop', overwrite=True)
    print(f"  REC_fixed.prmtop: {len(rec_parm.atoms)} atoms")

    # Same for ligand
    lig_parm = parmed.load_file('LIG.prmtop', '_GMXMMPBSA_LIG.inpcrd')
    lig_ref = parmed.load_file('_GMXMMPBSA_LIG.pdb')
    lig_amber = np.array([(a.xx, a.xy, a.xz) for a in lig_parm.atoms])
    lig_gro = np.array([(a.xx, a.xy, a.xz) for a in lig_ref.atoms])
    lig_unmatched = []
    lig_used = set()
    for i, ac in enumerate(lig_amber):
        dists = np.linalg.norm(lig_gro - ac, axis=1)
        nearest = np.argmin(dists)
        if dists[nearest] < 0.05 and nearest not in lig_used:
            lig_used.add(nearest)
        else:
            lig_unmatched.append(i)
    if lig_unmatched:
        lig_parm.strip('@' + ','.join(str(i+1) for i in lig_unmatched))
    lig_parm.save('LIG_fixed.prmtop', overwrite=True)
    print(f"  LIG_fixed.prmtop: {len(lig_parm.atoms)} atoms")

    # Verify
    c = parmed.load_file('COM_fixed.prmtop')
    r = parmed.load_file('REC_fixed.prmtop')
    l = parmed.load_file('LIG_fixed.prmtop')
    print(f"\n=== Final verification ===")
    print(f"  COM: {len(c.atoms)} (need 1715)")
    print(f"  REC: {len(r.atoms)} (need 1322)")
    print(f"  LIG: {len(l.atoms)} (need 393)")
    print(f"  REC+LIG: {len(r.atoms)+len(l.atoms)}")

    if len(c.atoms) == 1715 and len(r.atoms) == 1322 and len(l.atoms) == 393:
        print("\n  ALL COUNTS MATCH! Ready for MMPBSA.py")
    else:
        print("\n  Counts don't match exactly - may need further adjustment")
