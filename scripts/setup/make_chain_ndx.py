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
Create GROMACS index file with chain-based groups for the ternary complex.
Parses the GRO file to find chain boundaries (where residue number resets).
"""
import os
from pathlib import Path

GRO = Path(os.path.join(DATA_ROOT, "topology", "md.gro"))
NDX = Path(os.path.join(DATA_ROOT, "mmpbsa/mmpbsa.ndx"))

with open(GRO) as f:
    lines = f.readlines()

n_atoms = int(lines[1].strip())
print(f"Total atoms: {n_atoms}")

# Parse atom lines, find chain breaks (where resnum drops)
atoms = []
prev_resnum = -1
chain_id = 0
chain_starts = [0]

for i, line in enumerate(lines[2:2+n_atoms]):
    resnum = int(line[0:5].strip())
    resname = line[5:10].strip()
    atomname = line[10:15].strip()
    atom_num = i + 1  # 1-indexed
    
    if resnum < prev_resnum and resname not in ["SOL", "NA", "CL", "K"]:
        chain_id += 1
        chain_starts.append(i)
    
    atoms.append({
        "num": atom_num,
        "resnum": resnum,
        "resname": resname,
        "atomname": atomname,
        "chain": chain_id if resname not in ["SOL", "NA", "CL", "K"] else -1,
    })
    prev_resnum = resnum

# Identify protein chains
chains = {}
for a in atoms:
    if a["chain"] >= 0:
        if a["chain"] not in chains:
            chains[a["chain"]] = {"atoms": [], "first_resnum": a["resnum"], "first_resname": a["resname"]}
        chains[a["chain"]]["atoms"].append(a["num"])
        chains[a["chain"]]["last_resnum"] = a["resnum"]
        chains[a["chain"]]["last_resname"] = a["resname"]

# Identify chains
chain_names = {}
for cid, info in sorted(chains.items()):
    n_atoms_chain = len(info["atoms"])
    n_res = info["last_resnum"]
    first = info["first_resname"]
    print(f"Chain {cid}: {n_atoms_chain} atoms, residues 1-{n_res}, starts with {first}")
    
    if n_res > 300 and cid == 0:
        chain_names[cid] = "gp120"
    elif n_res < 200 and n_res > 100:
        chain_names[cid] = "CD4"
    elif n_res > 250:
        chain_names[cid] = "CCR5"
    else:
        chain_names[cid] = f"chain_{cid}"

print(f"\nChain mapping: {chain_names}")

# Write index file
with open(NDX, "w") as f:
    # Write each chain as a group
    for cid, info in sorted(chains.items()):
        name = chain_names.get(cid, f"chain_{cid}")
        atom_nums = info["atoms"]
        f.write(f"[ {name} ]\n")
        for i, a in enumerate(atom_nums):
            f.write(f"{a:7d}")
            if (i + 1) % 15 == 0:
                f.write("\n")
        if len(atom_nums) % 15 != 0:
            f.write("\n")
        f.write("\n")
    
    # Also write receptor = gp120 + CCR5
    receptor_atoms = []
    for cid in sorted(chains.keys()):
        if chain_names.get(cid) in ["gp120", "CCR5"]:
            receptor_atoms.extend(chains[cid]["atoms"])
    
    f.write("[ Receptor ]\n")
    for i, a in enumerate(receptor_atoms):
        f.write(f"{a:7d}")
        if (i + 1) % 15 == 0:
            f.write("\n")
    if len(receptor_atoms) % 15 != 0:
        f.write("\n")
    f.write("\n")

print(f"\nIndex file written: {NDX}")
print(f"Groups:")

# Show summary
with open(NDX) as f:
    content = f.read()
import re
groups = re.findall(r'\[ (\w+) \]', content)
for i, g in enumerate(groups):
    count = len(re.findall(r'\d+', content.split(f"[ {g} ]")[1].split("[")[0] if "[" in content.split(f"[ {g} ]")[1] else content.split(f"[ {g} ]")[1]))
    print(f"  {i}: {g} ({count} atoms)")
