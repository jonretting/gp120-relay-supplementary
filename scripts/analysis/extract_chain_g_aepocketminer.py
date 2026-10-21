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
Extract chain G (gp120 monomer) from multi-chain PDB files for
AE-PocketMiner GNN inference.

AE-PocketMiner (Zhang et al. 2026) is a single-chain GVP-GNN; multi-chain
input is out-of-distribution. This script isolates chain G from each source
PDB and writes standalone PDB files.

Input:  structures/{PDBID}.pdb  (full multi-chain PDB from RCSB)
Output: structures/aepocketminer/{PDBID}_chainG.pdb

Usage:
    python extract_chain_g_aepocketminer.py

Requires: No external dependencies (stdlib only).
"""

import sys

PDBS = ["6MEO", "5FYJ", "5U7O", "6MU6", "4TVP"]

def extract_chain_g(input_path, output_path):
    """Extract chain G ATOM/HETATM records from a PDB file."""
    kept = 0
    with open(input_path, "r") as fin, open(output_path, "w") as fout:
        for line in fin:
            rec = line[:6].strip()
            if rec in ("ATOM", "HETATM"):
                chain = line[21]
                if chain == "G":
                    fout.write(line)
                    kept += 1
            elif rec == "TER":
                if kept > 0:
                    fout.write(line)
            elif rec == "END":
                fout.write(line)
        if not line.startswith("END"):
            fout.write("END\n")
    return kept


def main():
    structures_dir = os.path.join(DATA_ROOT, "structures")
    output_dir = os.path.join(structures_dir, "aepocketminer")
    os.makedirs(output_dir, exist_ok=True)

    for pdb_id in PDBS:
        input_path = os.path.join(structures_dir, f"{pdb_id}.pdb")
        output_path = os.path.join(output_dir, f"{pdb_id}_chainG.pdb")

        if not os.path.isfile(input_path):
            print(f"WARNING: {input_path} not found — skipping {pdb_id}")
            print(f"  Download from: https://files.rcsb.org/download/{pdb_id}.pdb")
            continue

        kept = extract_chain_g(input_path, output_path)
        print(f"{pdb_id}: extracted {kept} ATOM/HETATM lines -> {output_path}")

    print("\nDone. Chain-G PDBs are ready for AE-PocketMiner inference.")


if __name__ == "__main__":
    main()
