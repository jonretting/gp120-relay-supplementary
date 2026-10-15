import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Align parameterized mol2 onto positioned PDB using Open Babel.
Then write aligned coordinates as a GRO file for GROMACS.
"""
import os
import subprocess
import numpy as np

WORKDIR = os.path.join(DATA_ROOT, "results", "04-md-simulations")
os.chdir(WORKDIR)

MOL2_FILE = os.path.join(WORKDIR, 'jyv_clean.acpype', 'jyv_clean_bcc_gaff2.mol2')
PDB_FILE = os.path.join(DATA_ROOT, "docking/bms626529_positioned.pdb")

# Use obabel to convert both to SDF (canonical format for alignment)
subprocess.run(['obabel', PDB_FILE, '-O', 'target.sdf', '--gen3d'], check=True,
               capture_output=True)
subprocess.run(['obabel', MOL2_FILE, '-O', 'source.sdf'], check=True,
               capture_output=True)

# Use obfit for structural alignment
result = subprocess.run(
    ['obfit', 'target.sdf', 'source.sdf'],
    capture_output=True, text=True
)
print('obfit result:', result.returncode)
if result.stdout:
    print(result.stdout[:500])
if result.stderr:
    print(result.stderr[:500])

# Alternative approach: use Python openbabel for alignment
from openbabel import openbabel as ob, pybel

# Load both molecules
ref_mol = next(pybel.readfile('pdb', PDB_FILE))
mob_mol = next(pybel.readfile('mol2', MOL2_FILE))

print(f'Reference (PDB): {len(ref_mol.atoms)} atoms')
print(f'Mobile (mol2):   {len(mob_mol.atoms)} atoms')

# Get heavy atom indices
ref_heavy = [(i, a) for i, a in enumerate(ref_mol.atoms) if a.atomicnum > 1]
mob_heavy = [(i, a) for i, a in enumerate(mob_mol.atoms) if a.atomicnum > 1]

print(f'Ref heavy atoms: {len(ref_heavy)}')
print(f'Mob heavy atoms: {len(mob_heavy)}')

# Use Open Babel's align functionality
aligner = ob.OBAlign()
aligner.SetRefMol(ref_mol.OBMol)
aligner.SetTargetMol(mob_mol.OBMol)
aligner.Align()
print(f'OBAlign RMSD: {aligner.GetRMSD():.4f}')

# Apply the rotation to the mobile molecule
aligner.UpdateCoords(mob_mol.OBMol)

# Write aligned mol2
mob_mol.write('mol2', 'aligned_bms626529.mol2', overwrite=True)

# Verify alignment by checking heavy atom distances
ref_coords = np.array([a.coords for _, a in ref_heavy])
aligned_coords = np.array([a.coords for _, a in mob_heavy])

# Per-atom distance
dists = np.sqrt(np.sum((ref_coords - aligned_coords[:len(ref_coords)])**2, axis=1))
print(f'Post-alignment per-atom distances (A):')
print(f'  Mean: {np.mean(dists):.3f}')
print(f'  Max:  {np.max(dists):.3f}')
print(f'  Min:  {np.min(dists):.3f}')

# Write as GRO file (nm units)
print('\nWriting aligned GRO file...')
with open('ligand_aligned.gro', 'w') as f:
    f.write('BMS-626529 aligned to ternary via OBAlign\n')
    f.write(f'{len(mob_mol.atoms)}\n')
    for i, a in enumerate(mob_mol.atoms):
        x_nm = a.coords[0] / 10.0
        y_nm = a.coords[1] / 10.0
        z_nm = a.coords[2] / 10.0
        name = a.OBAtom.GetType()[:5]
        if not name.strip():
            elem = ob.OBElementTable().GetSymbol(a.atomicnum)
            name = f'{elem}{i}'
        f.write(f'{1:>5d}{"UNL":>5s}{name:>5s}{i+1:>5d}{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}\n')
    f.write('   0.00000   0.00000   0.00000\n')

print(f'Wrote ligand_aligned.gro ({len(mob_mol.atoms)} atoms)')

# Also write as PDB for verification
mob_mol.write('pdb', 'aligned_bms626529.pdb', overwrite=True)
print('Wrote aligned_bms626529.pdb for verification')

# Verify some coordinates
print('\nSample aligned coords (A):')
for i in range(min(5, len(mob_mol.atoms))):
    a = mob_mol.atoms[i]
    print(f'  {a.OBAtom.GetType()}: ({a.coords[0]:.3f}, {a.coords[1]:.3f}, {a.coords[2]:.3f})')

print('\nReference PDB coords (A):')
for i in range(min(5, len(ref_mol.atoms))):
    a = ref_mol.atoms[i]
    print(f'  {a.OBAtom.GetType()}: ({a.coords[0]:.3f}, {a.coords[1]:.3f}, {a.coords[2]:.3f})')
