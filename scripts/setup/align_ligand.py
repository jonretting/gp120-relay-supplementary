import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Align the GAFF2-parameterized ligand (83 atoms, with H) to the 
positioned PDB (44 heavy atoms) in the ternary reference frame.

Uses Kabsch superposition on heavy atoms, then applies the same
transformation to hydrogens.
"""
import numpy as np
import os

WORKDIR = os.path.join(DATA_ROOT, "results", "04-md-simulations")
os.chdir(WORKDIR)

POSITIONED_PDB = os.path.join(DATA_ROOT, "docking/bms626529_positioned.pdb")
PARMED_GRO = os.path.join(WORKDIR, 'jyv_clean.acpype', 'jyv_clean_GMX.gro')
PARMED_MOL2 = os.path.join(WORKDIR, 'jyv_clean.acpype', 'jyv_clean_bcc_gaff2.mol2')


def parse_pdb_hetatm(path):
    """Parse HETATM records from PDB, return atom names + coords (Angstroms)."""
    atoms = []
    for line in open(path):
        if line.startswith('HETATM') or line.startswith('ATOM'):
            name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = name[0]
            atoms.append((name, element, np.array([x, y, z])))
    return atoms


def parse_gro(path):
    """Parse GRO file, return atom names + coords (convert nm -> Angstroms)."""
    atoms = []
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[1].strip())
    for i in range(2, 2 + n):
        line = lines[i]
        name = line[10:15].strip()
        x = float(line[20:28]) * 10.0  # nm -> A
        y = float(line[28:36]) * 10.0
        z = float(line[36:44]) * 10.0
        is_h = name.startswith('H') or name.startswith('h')
        atoms.append((name, 'H' if is_h else name[0], np.array([x, y, z])))
    return atoms


def kabsch(P, Q):
    """Compute rotation matrix R and translation to align P onto Q.
    P, Q are Nx3 arrays. Returns R, t such that P_aligned = (P - centroid_P) @ R + centroid_Q.
    """
    centroid_P = P.mean(axis=0)
    centroid_Q = Q.mean(axis=0)
    P_c = P - centroid_P
    Q_c = Q - centroid_Q
    H = P_c.T @ Q_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, d])
    R = Vt.T @ sign_matrix @ U.T
    return R, centroid_P, centroid_Q


# Parse both structures
pdb_atoms = parse_pdb_hetatm(POSITIONED_PDB)
gro_atoms = parse_gro(PARMED_GRO)

print(f'Positioned PDB: {len(pdb_atoms)} atoms (heavy only)')
print(f'Parmed GRO: {len(gro_atoms)} atoms (with H)')

# Extract heavy atom coords from both
pdb_heavy = [(a[0], a[2]) for a in pdb_atoms if a[1] != 'H']
gro_heavy = [(a[0], a[2]) for a in gro_atoms if a[1] != 'H']

print(f'PDB heavy atoms: {len(pdb_heavy)}')
print(f'GRO heavy atoms: {len(gro_heavy)}')

if len(pdb_heavy) != len(gro_heavy):
    print(f'WARNING: Heavy atom count mismatch! PDB={len(pdb_heavy)}, GRO={len(gro_heavy)}')
    n_match = min(len(pdb_heavy), len(gro_heavy))
else:
    n_match = len(pdb_heavy)

P = np.array([a[1] for a in gro_heavy[:n_match]])
Q = np.array([a[1] for a in pdb_heavy[:n_match]])

R, centroid_P, centroid_Q = kabsch(P, Q)

# Apply transformation to ALL atoms (heavy + H)
aligned_coords = []
for name, element, coord in gro_atoms:
    aligned = (coord - centroid_P) @ R + centroid_Q
    aligned_coords.append((name, aligned))

# Compute alignment RMSD on heavy atoms
aligned_heavy = [(coord - centroid_P) @ R + centroid_Q for _, _, coord in gro_atoms if not (_.startswith('H') or _.startswith('h'))]
rmsd = np.sqrt(np.mean(np.sum((np.array(aligned_heavy[:n_match]) - Q)**2, axis=1)))
print(f'Alignment RMSD (heavy atoms): {rmsd:.3f} A')

# Write aligned GRO (in nm)
with open('ligand_aligned.gro', 'w') as f:
    f.write('BMS-626529 aligned to ternary\n')
    f.write(f'{len(aligned_coords)}\n')
    for i, (name, coord) in enumerate(aligned_coords):
        x_nm = coord[0] / 10.0
        y_nm = coord[1] / 10.0
        z_nm = coord[2] / 10.0
        f.write(f'{1:>5d}{"UNL":>5s}{name:>5s}{i+1:>5d}{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}\n')
    f.write('   0.00000   0.00000   0.00000\n')

print(f'Wrote ligand_aligned.gro ({len(aligned_coords)} atoms)')
print('Coords sample (nm):')
for name, coord in aligned_coords[:3]:
    print(f'  {name}: ({coord[0]/10:.3f}, {coord[1]/10:.3f}, {coord[2]/10:.3f})')
