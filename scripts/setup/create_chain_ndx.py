import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Create chain-specific index groups from a multi-chain GROMACS .gro file.

For the gp120-CD4-CCR5 ternary complex:
  Chain A = gp120 (398 residues)
  Chain B = CD4   (176 residues)
  Chain C = CCR5  (313 residues)

Chain boundaries detected by residue numbering reset.
"""

import sys
import os

def parse_gro_chains(gro_path):
    with open(gro_path) as f:
        lines = f.readlines()

    n_atoms = int(lines[1].strip())
    chains = {'A': [], 'B': [], 'C': []}
    chain_labels = ['A', 'B', 'C']
    current_idx = 0
    prev_resid = 0

    for i in range(2, 2 + n_atoms):
        line = lines[i]
        if len(line) < 44:
            continue
        try:
            resid = int(line[:5].strip())
        except ValueError:
            continue

        resname = line[5:10].strip()
        atom_name = line[10:15].strip()
        atom_nr = i - 1  # 1-based atom number

        # Non-protein atoms stop chain detection
        if resname in ('SOL', 'NA', 'CL', 'NA+', 'CL-', 'WAT', 'HOH'):
            break

        if resid < prev_resid and prev_resid > 50:
            current_idx += 1
            if current_idx >= len(chain_labels):
                break
        prev_resid = resid

        chains[chain_labels[current_idx]].append(atom_nr)

    return chains


def write_ndx(ndx_path, groups):
    with open(ndx_path, 'w') as f:
        for name, atoms in groups.items():
            f.write(f'[ {name} ]\n')
            for j, a in enumerate(atoms, 1):
                f.write(f'{a:7d}')
                if j % 15 == 0:
                    f.write('\n')
            if len(atoms) % 15 != 0:
                f.write('\n')


def main():
    gro_path = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser(
        os.path.join(DATA_ROOT, 'md/md.gro'))
    out_ndx = sys.argv[2] if len(sys.argv) > 2 else 'chains.ndx'

    chains = parse_gro_chains(gro_path)

    # Also create a relay pocket group (chain A residues 372-379)
    relay_resids = set(range(372, 380))
    relay_atoms = []
    with open(gro_path) as f:
        lines = f.readlines()
    n_atoms = int(lines[1].strip())
    current_chain = 0
    prev_resid = 0
    for i in range(2, 2 + n_atoms):
        line = lines[i]
        if len(line) < 44:
            continue
        try:
            resid = int(line[:5].strip())
        except ValueError:
            continue
        resname = line[5:10].strip()
        if resname in ('SOL', 'NA', 'CL', 'NA+', 'CL-', 'WAT', 'HOH'):
            break
        if resid < prev_resid and prev_resid > 50:
            current_chain += 1
        prev_resid = resid
        atom_nr = i - 1
        if current_chain == 0 and resid in relay_resids:
            relay_atoms.append(atom_nr)

    groups = {
        'gp120': chains['A'],
        'CD4': chains['B'],
        'CCR5': chains['C'],
        'Relay_Pocket': relay_atoms,
    }

    write_ndx(out_ndx, groups)

    for name, atoms in groups.items():
        print(f'{name}: {len(atoms)} atoms '
              f'(#{atoms[0]}-#{atoms[-1]})' if atoms else f'{name}: 0 atoms')


if __name__ == '__main__':
    main()
