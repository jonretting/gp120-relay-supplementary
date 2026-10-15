#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Inspect SOSIP PDB: chain layout and relay pocket location."""
import sys
from collections import defaultdict

pdb_path = sys.argv[1]
with open(pdb_path) as f:
    lines = f.readlines()

chains = {}
for line in lines:
    if line[:4] == 'ATOM':
        chain = line[21]
        resnum = int(line[22:26].strip())
        if chain not in chains:
            chains[chain] = {'min': resnum, 'max': resnum, 'atoms': 0, 'residues': set()}
        chains[chain]['min'] = min(chains[chain]['min'], resnum)
        chains[chain]['max'] = max(chains[chain]['max'], resnum)
        chains[chain]['atoms'] += 1
        chains[chain]['residues'].add(resnum)

print(f"=== {pdb_path} Chain Summary ===")
for c in sorted(chains):
    n = len(chains[c]['residues'])
    print(f"  Chain {c}: res {chains[c]['min']}-{chains[c]['max']} ({n} unique), {chains[c]['atoms']} atoms")

chain_seq = defaultdict(dict)
for line in lines:
    if line[:4] == 'ATOM' and line[12:16].strip() == 'CA':
        chain = line[21]
        resnum = int(line[22:26].strip())
        resname = line[17:20].strip()
        chain_seq[chain][resnum] = resname

aa3to1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
           'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
           'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M','TYS':'Y'}

target = 'RPGGGDMR'
print(f"\n=== Relay Pocket Search (RPGGGDMR) ===")
relay_chains = []
for c in sorted(chain_seq):
    resnums = sorted(chain_seq[c].keys())
    seq = ''
    for r in resnums:
        seq += aa3to1.get(chain_seq[c][r], 'X')
    idx = seq.find(target)
    if idx >= 0:
        start_res = resnums[idx]
        end_res = resnums[idx + len(target) - 1]
        print(f"  Chain {c}: RPGGGDMR at residues {start_res}-{end_res}")
        relay_chains.append((c, start_res, end_res))
    else:
        print(f"  Chain {c}: no RPGGGDMR")

if relay_chains:
    best = relay_chains[0]
    print(f"\n=== Recommended: Chain {best[0]}, relay at {best[1]}-{best[2]} ===")
    
    # Find center of mass of relay pocket CA atoms
    cx, cy, cz, n = 0, 0, 0, 0
    for line in lines:
        if line[:4] == 'ATOM' and line[21] == best[0] and line[12:16].strip() == 'CA':
            resnum = int(line[22:26].strip())
            if best[1] <= resnum <= best[2]:
                cx += float(line[30:38])
                cy += float(line[38:46])
                cz += float(line[46:54])
                n += 1
    if n:
        cx, cy, cz = cx/n, cy/n, cz/n
        print(f"  Pocket center: ({cx:.1f}, {cy:.1f}, {cz:.1f})")
        
        # Count atoms within 15A of pocket center
        nearby = 0
        nearby_chains = set()
        nearby_residues = set()
        for line in lines:
            if line[:4] == 'ATOM':
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                d = ((x-cx)**2 + (y-cy)**2 + (z-cz)**2)**0.5
                if d <= 15.0:
                    resname = line[17:20].strip()
                    if resname not in ('HOH','WAT','SOL','NA','CL'):
                        nearby += 1
                        nearby_chains.add(line[21])
                        nearby_residues.add((line[21], int(line[22:26].strip())))
        print(f"  Atoms within 15A: {nearby}")
        print(f"  Chains involved: {sorted(nearby_chains)}")
        print(f"  Unique residues within 15A: {len(nearby_residues)}")

# Check for non-standard residues and HETATM ligands
nonstd = set()
hetatm_types = set()
for line in lines:
    if line[:6] == 'HETATM':
        resname = line[17:20].strip()
        if resname not in ('HOH','WAT','SOL','NA','CL','MG','ZN','CA','FE'):
            hetatm_types.add(resname)
    if line[:4] == 'ATOM':
        resname = line[17:20].strip()
        if resname == 'TYS':
            nonstd.add('TYS (sulfotyrosine)')
        elif resname == 'MSE':
            nonstd.add('MSE (selenomethionine)')

if nonstd:
    print(f"\n  Non-standard residues: {', '.join(sorted(nonstd))}")
if hetatm_types:
    print(f"  HETATM ligands: {', '.join(sorted(hetatm_types))}")
