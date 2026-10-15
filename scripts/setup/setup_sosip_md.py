#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Setup reduced-system MD for SOSIP X-ray structures.

Extracts a 15A shell around the relay pocket from chain G,
cleans non-standard residues, writes GROMACS-ready PDB.
"""
import sys
import os
import subprocess
import math

PDB_FILE = sys.argv[1]
RELAY_CHAIN = sys.argv[2] if len(sys.argv) > 2 else 'G'
RELAY_START = int(sys.argv[3]) if len(sys.argv) > 3 else 469
RELAY_END = int(sys.argv[4]) if len(sys.argv) > 4 else 476
SHELL_RADIUS = 15.0
BOX_EDGE = 7.0
BOX_BUFFER = 1.5
SALT_CONC = 0.15

pdb_id = os.path.splitext(os.path.basename(PDB_FILE))[0]
workdir = os.getcwd()

print(f"=== SOSIP MD Setup: {pdb_id} ===")
print(f"  Relay chain: {RELAY_CHAIN}, residues {RELAY_START}-{RELAY_END}")
print(f"  Shell radius: {SHELL_RADIUS} A")

with open(PDB_FILE) as f:
    lines = f.readlines()

# Find relay pocket center of mass (CA atoms)
cx, cy, cz, n = 0, 0, 0, 0
for line in lines:
    if line[:4] == 'ATOM' and line[21] == RELAY_CHAIN and line[12:16].strip() == 'CA':
        resnum = int(line[22:26].strip())
        if RELAY_START <= resnum <= RELAY_END:
            cx += float(line[30:38])
            cy += float(line[38:46])
            cz += float(line[46:54])
            n += 1

if n == 0:
    print("ERROR: No relay pocket CA atoms found!")
    sys.exit(1)

cx, cy, cz = cx/n, cy/n, cz/n
print(f"  Pocket center: ({cx:.1f}, {cy:.1f}, {cz:.1f})")

# Non-standard residue replacements
REPLACE = {'TYS': 'TYR', 'MSE': 'MET', 'SEC': 'CYS', 'CSE': 'CYS', 'HSD': 'HIS', 'HSE': 'HIS'}

# Two-pass: first find which residues have ANY atom within the shell,
# then extract ALL atoms of those residues (complete residues only).
shell_residues = set()
for line in lines:
    if line[:4] != 'ATOM':
        continue
    resname = line[17:20].strip()
    if resname in ('HOH', 'WAT', 'SOL', 'NA', 'CL', 'MG', 'ZN', 'CA', 'FE'):
        continue
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    d = math.sqrt((x-cx)**2 + (y-cy)**2 + (z-cz)**2)
    if d <= SHELL_RADIUS:
        chain = line[21]
        resnum = int(line[22:26].strip())
        shell_residues.add((chain, resnum))

out_lines = []
kept_residues = set()
kept_chains = set()
replaced = {}
prev_chain = None

for line in lines:
    if line[:4] != 'ATOM':
        continue
    resname = line[17:20].strip()
    if resname in ('HOH', 'WAT', 'SOL', 'NA', 'CL', 'MG', 'ZN', 'CA', 'FE'):
        continue
    chain = line[21]
    resnum = int(line[22:26].strip())
    if (chain, resnum) not in shell_residues:
        continue

    if prev_chain is not None and chain != prev_chain:
        out_lines.append("TER\n")
    prev_chain = chain

    if resname in REPLACE:
        new_name = REPLACE[resname]
        replaced[resname] = replaced.get(resname, 0) + 1
        line = line[:17] + f"{new_name:>3s}" + line[20:]

    kept_residues.add((chain, resnum))
    kept_chains.add(chain)
    out_lines.append(line)

out_lines.append("TER\nEND\n")

out_pdb = f"{pdb_id}_reduced.pdb"
with open(out_pdb, 'w') as f:
    f.writelines(out_lines)

n_atoms = sum(1 for l in out_lines if l[:4] == 'ATOM')
print(f"  Extracted: {n_atoms} atoms, {len(kept_residues)} residues, chains {sorted(kept_chains)}")
if replaced:
    for k, v in replaced.items():
        print(f"  Fixed: {k} -> {REPLACE[k]} ({v} atoms)")

# Verify relay pocket is complete
relay_found = 0
for resnum in range(RELAY_START, RELAY_END + 1):
    if (RELAY_CHAIN, resnum) in kept_residues:
        relay_found += 1
print(f"  Relay pocket: {relay_found}/8 residues present")
if relay_found < 8:
    print("  WARNING: Incomplete relay pocket!")

print(f"\n  Output: {out_pdb}")
print(f"\n=== Ready for GROMACS pipeline ===")
print(f"  Next steps:")
print(f"  1. gmx pdb2gmx -f {out_pdb} -o processed.gro -water tip3p -ff amber99sb-ildn -ignh")
print(f"  2. gmx editconf -f processed.gro -o box.gro -c -d {BOX_BUFFER} -bt cubic")
print(f"  3. gmx solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top")
print(f"  4. gmx grompp + gmx genion (0.15 M NaCl)")
print(f"  5. Energy minimization")
print(f"  6. NVT equilibration (100 ps)")
print(f"  7. NPT equilibration (100 ps)")
print(f"  8. Production MD (20 ns)")
