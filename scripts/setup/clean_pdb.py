#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Clean PDB for GROMACS: keep only standard amino acids, fix TYS->TYR."""
import sys
from pathlib import Path

STANDARD_AA = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
               'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}

TYR_ATOMS = {'N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ','OH',
             'H','HA','HB2','HB3','HD1','HD2','HE1','HE2','HH','H1','H2','H3','OXT'}

infile = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('hiv_6meo_fixed.pdb')
outfile = Path(sys.argv[2]) if len(sys.argv) > 2 else Path('hiv_6meo_clean.pdb')

chain_stats = {}
out_lines = []
prev_chain = None

with open(infile) as f:
    for line in f:
        if not line.startswith(('ATOM', 'HETATM')):
            if line.startswith('CRYST'):
                out_lines.append(line)
            continue

        resname = line[17:20].strip()
        atom_name = line[12:16].strip()
        chain = line[21]

        if resname == 'TYS':
            if atom_name not in TYR_ATOMS:
                continue
            line = line[:17] + 'TYR' + line[20:]
            resname = 'TYR'

        if resname not in STANDARD_AA:
            continue

        if prev_chain is not None and chain != prev_chain:
            out_lines.append('TER\n')
        prev_chain = chain

        if chain not in chain_stats:
            chain_stats[chain] = {'atoms': 0, 'residues': set()}
        chain_stats[chain]['atoms'] += 1
        chain_stats[chain]['residues'].add(line[22:26].strip())

        out_lines.append(line)

out_lines.append('TER\n')
out_lines.append('END\n')

with open(outfile, 'w') as f:
    f.writelines(out_lines)

total = sum(v['atoms'] for v in chain_stats.values())
print(f'Clean PDB written to {outfile}: {total} atoms')
for c in sorted(chain_stats):
    info = chain_stats[c]
    print(f'  Chain {c}: {info["atoms"]} atoms, {len(info["residues"])} residues')
