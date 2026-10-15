import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Find RPGGGDMR relay pocket in the full ternary complex GRO file."""
import os

aa3to1 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q",
    "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
    "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}

gro = os.path.join(DATA_ROOT, "topology", "md.gro")

with open(gro) as f:
    lines = f.readlines()

n_atoms = int(lines[1].strip())

# Build per-chain sequences
chains = []
current_chain = []
prev_resnum = -1

for line in lines[2:2+n_atoms]:
    resnum = int(line[0:5].strip())
    resname = line[5:10].strip()
    
    if resname in ["SOL", "NA", "CL", "K"]:
        continue
    
    if resnum < prev_resnum and current_chain:
        chains.append(current_chain)
        current_chain = []
    
    if not current_chain or current_chain[-1][0] != resnum:
        current_chain.append((resnum, resname))
    
    prev_resnum = resnum

if current_chain:
    chains.append(current_chain)

for ci, chain in enumerate(chains):
    name = ["gp120", "CD4", "CCR5"][ci] if ci < 3 else f"chain_{ci}"
    seq = "".join(aa3to1.get(r[1], "?") for r in chain)
    print(f"\n{name} ({len(chain)} residues):")
    
    # Search for relay pocket motifs
    for motif in ["RPGGGDMR", "RPAGGDMR", "RPGGGDM", "PGGGDMR", "PGGGD"]:
        idx = seq.find(motif)
        if idx >= 0:
            start_resnum = chain[idx][0]
            end_resnum = chain[idx + len(motif) - 1][0]
            print(f"  Found {motif} at seq position {idx}, residues {start_resnum}-{end_resnum}")
            print(f"  Residues:")
            for j in range(len(motif)):
                rn, rname = chain[idx + j]
                print(f"    {rn:4d} {rname:>3s} ({motif[j]})")
    
    # Also print sequence around residue 69-76 and around the end
    if ci == 0:
        print(f"\n  Sequence around res 69-76: ...{seq[66:79]}...")
        print(f"  Residues 69-76: {[(chain[i][0], chain[i][1]) for i in range(68, min(76, len(chain)))]}")
        
        # Show last 50 residues
        print(f"\n  Last 20 residues:")
        for rn, rname in chain[-20:]:
            print(f"    {rn:4d} {rname:>3s} ({aa3to1.get(rname, '?')})")
