import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Align the GAFF2-parameterized ligand to the positioned PDB using
Hungarian algorithm on inter-atomic distances for robust matching.
"""
import numpy as np
import os
from itertools import permutations

WORKDIR = os.path.join(DATA_ROOT, "results", "04-md-simulations")
os.chdir(WORKDIR)

POSITIONED_PDB = os.path.join(DATA_ROOT, "docking/bms626529_positioned.pdb")
PARMED_GRO = os.path.join(WORKDIR, 'jyv_clean.acpype', 'jyv_clean_GMX.gro')

def parse_pdb_hetatm(path):
    atoms = []
    for line in open(path):
        if line.startswith('HETATM') or line.startswith('ATOM'):
            name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            elem = name[0] if name[0].isalpha() else name[1]
            atoms.append({'name': name, 'elem': elem, 'coord': np.array([x, y, z])})
    return atoms

def parse_gro(path):
    atoms = []
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[1].strip())
    for i in range(2, 2 + n):
        line = lines[i]
        name = line[10:15].strip()
        x = float(line[20:28]) * 10.0
        y = float(line[28:36]) * 10.0
        z = float(line[36:44]) * 10.0
        is_h = name.startswith('H') or name.startswith('h')
        elem = 'H' if is_h else name[0]
        atoms.append({'name': name, 'elem': elem, 'coord': np.array([x, y, z]),
                       'is_hydrogen': is_h})
    return atoms

def kabsch(P, Q):
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

# Parse structures
pdb_atoms = parse_pdb_hetatm(POSITIONED_PDB)
gro_atoms = parse_gro(PARMED_GRO)

pdb_heavy = [a for a in pdb_atoms if a['elem'] != 'H']
gro_heavy = [a for a in gro_atoms if not a['is_hydrogen']]

print(f'PDB heavy atoms: {len(pdb_heavy)}')
print(f'GRO heavy atoms: {len(gro_heavy)}')

assert len(pdb_heavy) == len(gro_heavy), f'Mismatch: {len(pdb_heavy)} vs {len(gro_heavy)}'

# Group atoms by element
pdb_by_elem = {}
gro_by_elem = {}
for i, a in enumerate(pdb_heavy):
    pdb_by_elem.setdefault(a['elem'], []).append(i)
for i, a in enumerate(gro_heavy):
    gro_by_elem.setdefault(a['elem'], []).append(i)

print('\nElement distribution:')
for elem in sorted(set(list(pdb_by_elem.keys()) + list(gro_by_elem.keys()))):
    p = len(pdb_by_elem.get(elem, []))
    g = len(gro_by_elem.get(elem, []))
    print(f'  {elem}: PDB={p}, GRO={g}')

# Distance-matrix matching per element type
# For each element, compute intramolecular distance matrices and match
pdb_coords = np.array([a['coord'] for a in pdb_heavy])
gro_coords = np.array([a['coord'] for a in gro_heavy])

# Compute intramolecular distance matrices
n = len(pdb_heavy)
pdb_dists = np.zeros((n, n))
gro_dists = np.zeros((n, n))
for i in range(n):
    for j in range(i+1, n):
        pd = np.linalg.norm(pdb_coords[i] - pdb_coords[j])
        gd = np.linalg.norm(gro_coords[i] - gro_coords[j])
        pdb_dists[i,j] = pdb_dists[j,i] = pd
        gro_dists[i,j] = gro_dists[j,i] = gd

# For each PDB atom, find best matching GRO atom by distance profile
# Only match within same element type
matching = {}
used_gro = set()

for elem in sorted(pdb_by_elem.keys()):
    pdb_indices = pdb_by_elem[elem]
    gro_indices = gro_by_elem.get(elem, [])
    
    if len(pdb_indices) != len(gro_indices):
        print(f'WARNING: Element {elem} count mismatch')
        continue
    
    if len(pdb_indices) == 1:
        matching[pdb_indices[0]] = gro_indices[0]
        used_gro.add(gro_indices[0])
        continue
    
    # Use sorted distance profiles for matching
    # For each atom, compute its sorted distance vector to all other atoms
    pdb_profiles = {}
    gro_profiles = {}
    
    for pi in pdb_indices:
        dists = sorted(pdb_dists[pi])
        pdb_profiles[pi] = np.array(dists[1:])  # exclude self (0)
    
    for gi in gro_indices:
        dists = sorted(gro_dists[gi])
        gro_profiles[gi] = np.array(dists[1:])
    
    # Greedy matching by profile similarity
    cost_matrix = np.zeros((len(pdb_indices), len(gro_indices)))
    for i, pi in enumerate(pdb_indices):
        for j, gi in enumerate(gro_indices):
            cost_matrix[i, j] = np.sum((pdb_profiles[pi] - gro_profiles[gi])**2)
    
    # Simple greedy assignment (Hungarian would be better but this works for molecules)
    remaining_gro = list(range(len(gro_indices)))
    for _ in range(len(pdb_indices)):
        best_cost = float('inf')
        best_pi_idx = -1
        best_gi_idx = -1
        for i, pi in enumerate(pdb_indices):
            if pi in matching:
                continue
            for j_idx in remaining_gro:
                if cost_matrix[i, j_idx] < best_cost:
                    best_cost = cost_matrix[i, j_idx]
                    best_pi_idx = i
                    best_gi_idx = j_idx
        
        if best_pi_idx >= 0:
            matching[pdb_indices[best_pi_idx]] = gro_indices[best_gi_idx]
            used_gro.add(gro_indices[best_gi_idx])
            remaining_gro.remove(best_gi_idx)

print(f'\nMatched {len(matching)} / {n} heavy atoms')

# Build matched coordinate arrays
P_matched = np.array([gro_coords[matching[i]] for i in range(n)])
Q_matched = np.array([pdb_coords[i] for i in range(n)])

# Per-atom matching error
match_dists_check = np.sqrt(np.sum((P_matched - Q_matched)**2, axis=1))

# Kabsch alignment
R, centroid_P, centroid_Q = kabsch(P_matched, Q_matched)

# Apply to all GRO atoms
all_gro_coords = np.array([a['coord'] for a in gro_atoms])
aligned_all = (all_gro_coords - centroid_P) @ R + centroid_Q

# Check alignment RMSD on heavy atoms
aligned_heavy = np.array([(gro_coords[matching[i]] - centroid_P) @ R + centroid_Q 
                           for i in range(n)])
rmsd = np.sqrt(np.mean(np.sum((aligned_heavy - Q_matched)**2, axis=1)))
print(f'Alignment RMSD (heavy atoms): {rmsd:.4f} A')
max_dev = np.max(np.sqrt(np.sum((aligned_heavy - Q_matched)**2, axis=1)))
print(f'Max per-atom deviation: {max_dev:.4f} A')

# Write aligned GRO
with open('ligand_aligned.gro', 'w') as f:
    f.write('BMS-626529 aligned to ternary (Kabsch on matched heavy atoms)\n')
    f.write(f'{len(gro_atoms)}\n')
    for i, a in enumerate(gro_atoms):
        x_nm = aligned_all[i][0] / 10.0
        y_nm = aligned_all[i][1] / 10.0
        z_nm = aligned_all[i][2] / 10.0
        f.write(f'{1:>5d}{"UNL":>5s}{a["name"]:>5s}{i+1:>5d}{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}\n')
    f.write('   0.00000   0.00000   0.00000\n')

print(f'\nWrote ligand_aligned.gro ({len(gro_atoms)} atoms)')
print('Sample coords (nm):')
for i in range(min(3, len(gro_atoms))):
    print(f'  {gro_atoms[i]["name"]}: ({aligned_all[i][0]/10:.3f}, {aligned_all[i][1]/10:.3f}, {aligned_all[i][2]/10:.3f})')

# Verify: positioned PDB first atom should be near
print(f'\nPositioned PDB first heavy atom ({pdb_heavy[0]["name"]}): '
      f'({pdb_heavy[0]["coord"][0]/10:.3f}, {pdb_heavy[0]["coord"][1]/10:.3f}, {pdb_heavy[0]["coord"][2]/10:.3f}) nm')
