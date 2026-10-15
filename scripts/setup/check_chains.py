import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Check chain composition of the ternary complex GRO file."""
import sys

gro_file = sys.argv[1] if len(sys.argv) > 1 else os.path.join(DATA_ROOT, "topology", "md.gro")

with open(gro_file) as f:
    lines = f.readlines()

print(f"Total lines: {len(lines)}")
print(f"Atoms: {lines[1].strip()}")

protein_lines = []
for line in lines[2:-1]:
    if len(line) > 10:
        resname = line[5:10].strip()
        if resname not in ["SOL", "NA", "CL", "K"]:
            protein_lines.append(line)

print(f"Protein atom lines: {len(protein_lines)}")

resnames = []
prev = None
for l in protein_lines:
    rn = l[0:5].strip()
    rname = l[5:10].strip()
    if rn != prev:
        resnames.append((int(rn), rname))
        prev = rn

print(f"Unique residues: {len(resnames)}")

# Find chain breaks
breaks = [0]
for i in range(1, len(resnames)):
    curr_num = resnames[i][0]
    prev_num = resnames[i-1][0]
    if curr_num <= prev_num:
        breaks.append(i)
        print(f"Chain break at index {i}: res {resnames[i-1]} -> {resnames[i]}")

breaks.append(len(resnames))

for ci in range(len(breaks)-1):
    start = breaks[ci]
    end = breaks[ci+1]
    chain_res = resnames[start:end]
    print(f"\nChain {ci} ({chr(65+ci)}): {len(chain_res)} residues")
    print(f"  Range: {chain_res[0]} to {chain_res[-1]}")
    seq = "".join([r[1][:3] for r in chain_res[:10]])
    print(f"  First 10 residues: {[r for r in chain_res[:10]]}")
    print(f"  Last 5 residues: {[r for r in chain_res[-5:]]}")

    # Count atoms for this chain
    chain_resnums = set(r[0] for r in chain_res)
    chain_atoms = sum(1 for l in protein_lines if int(l[0:5].strip()) in chain_resnums)
    print(f"  Atoms: {chain_atoms}")
