#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Compute DCCM for SOSIP 500ns trajectories (5U7O, 6MU6).

Adapts compute_dccm_500ns.py for SOSIP reduced systems (chain G only, ~85 CA).
Relay pocket residues use SOSIP numbering (+4 from HXB2): 469-476.

Usage:
  python3 compute_sosip_dccm_500ns.py <PDB_ID> <covar_ascii.dat> <md.gro>

Requires: gmx covar already run to produce covar_ascii.dat.
  echo 'All_CA All_CA' | gmx covar -s md_500ns.tpr -f md.xtc -n pocket.ndx -ascii covar_ascii.dat -dt 50
"""
import sys
import os
import numpy as np
import json
import math

PDB_ID = sys.argv[1] if len(sys.argv) > 1 else "5U7O"
COVAR_FILE = sys.argv[2] if len(sys.argv) > 2 else "covar_ascii.dat"
GRO_FILE = sys.argv[3] if len(sys.argv) > 3 else "md.gro"

aa3to1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','HIE':'H','HID':'H'
}

print(f"=== {PDB_ID} 500ns SOSIP DCCM ===")

with open(GRO_FILE) as f:
    gro_lines = f.readlines()

resid_map = {}
ca_resnames = {}
ca_idx = 0
for line in gro_lines[2:-1]:
    if len(line) < 44:
        continue
    resid = int(line[:5].strip())
    resname = line[5:10].strip()
    atom_name = line[10:15].strip()
    if atom_name == "CA":
        resid_map[ca_idx] = resid
        ca_resnames[ca_idx] = resname
        ca_idx += 1

N_ATOMS = len(resid_map)
N_COORDS = N_ATOMS * 3
print(f"Total CA atoms: {N_ATOMS}")

seq = ""
for i in range(N_ATOMS):
    seq += aa3to1.get(ca_resnames[i], 'X')

import re
motif_idx = -1
target_motif = None
for variant in ["RPGGGDMR", "RPAGGDMR"]:
    motif_idx = seq.find(variant)
    if motif_idx >= 0:
        target_motif = variant
        break

if motif_idx < 0:
    m = re.search(r'RP[GA]GGDMR', seq)
    if m:
        motif_idx = m.start()
        target_motif = m.group()

if motif_idx < 0:
    print("ERROR: relay motif not found in sequence!")
    sys.exit(1)

relay_indices = list(range(motif_idx, motif_idx + 8))
relay_resids = [resid_map[i] for i in relay_indices]
relay_names = [f"{ca_resnames[i]}{resid_map[i]}" for i in relay_indices]
print(f"Relay motif: {target_motif}")
print(f"Relay CA indices: {relay_indices}")
print(f"Relay residues: {relay_names}")

print(f"\nLoading covariance matrix from {COVAR_FILE}...")
covar = np.zeros((N_COORDS, N_COORDS))
with open(COVAR_FILE) as f:
    row = 0
    col = 0
    for line in f:
        vals = line.strip().split()
        for v in vals:
            covar[row, col] = float(v)
            col += 1
            if col >= N_COORDS:
                col = 0
                row += 1

assert row == N_COORDS and col == 0, f"Parsed {row} rows, {col} leftover (expected {N_COORDS})"
print(f"Covariance matrix: {covar.shape}, Trace: {np.trace(covar):.3f} nm^2")

dccm = np.zeros((N_ATOMS, N_ATOMS))
for i in range(N_ATOMS):
    for j in range(N_ATOMS):
        cov_ij = sum(covar[3*i+a, 3*j+a] for a in range(3))
        cov_ii = sum(covar[3*i+a, 3*i+a] for a in range(3))
        cov_jj = sum(covar[3*j+a, 3*j+a] for a in range(3))
        denom = math.sqrt(cov_ii * cov_jj)
        dccm[i, j] = cov_ij / denom if denom > 0 else 0.0

print(f"DCCM: {dccm.shape}, Range: [{dccm.min():.4f}, {dccm.max():.4f}]")

print("\n=== Intra-Pocket Cross-Correlations ===")
intra_pairs = []
all_positive = True
for a_idx, a in enumerate(relay_indices):
    for b_idx, b in enumerate(relay_indices):
        if b_idx > a_idx:
            c = dccm[a, b]
            pair_name = f"{relay_names[a_idx]}-{relay_names[b_idx]}"
            intra_pairs.append((pair_name, c))
            if c <= 0:
                all_positive = False
            print(f"  {pair_name:30s}: {'+' if c > 0 else '-'}{abs(c):.4f}")

intra_values = [p[1] for p in intra_pairs]
seq_pairs = [dccm[relay_indices[i], relay_indices[i+1]] for i in range(len(relay_indices)-1)]
mean_intra = sum(intra_values) / len(intra_values)
mean_seq = sum(seq_pairs) / len(seq_pairs)

print(f"\n  All 28 positive: {all_positive}")
print(f"  Mean intra-pocket: {mean_intra:+.4f}")
print(f"  Range: {min(intra_values):+.4f} to {max(intra_values):+.4f}")
print(f"  Sequential neighbor mean: {mean_seq:+.4f}")

print("\n--- Comparison to 6MEO 500ns ---")
print(f"  6MEO 500ns: All 28 positive, mean +0.79, range +0.37 to +0.99, seq +0.97")
print(f"  {PDB_ID} 500ns: {'All 28 positive' if all_positive else 'NOT all positive'}, mean {mean_intra:+.4f}, range {min(intra_values):+.4f} to {max(intra_values):+.4f}, seq {mean_seq:+.4f}")

print("\n=== Betweenness Centrality ===")
try:
    import networkx as nx
    G = nx.Graph()
    for i in range(N_ATOMS):
        G.add_node(i)
    for i in range(N_ATOMS):
        for j in range(i+1, N_ATOMS):
            if abs(dccm[i, j]) > 0.2:
                G.add_edge(i, j, weight=abs(dccm[i, j]))

    bc = nx.betweenness_centrality(G, weight='weight')
    bc_sorted = sorted(bc.items(), key=lambda x: x[1], reverse=True)

    print("  Top 10:")
    top10 = []
    for rank, (node, cent) in enumerate(bc_sorted[:10], 1):
        resid = resid_map.get(node, '?')
        rname = ca_resnames.get(node, '?')
        is_relay = '*' if node in relay_indices else ' '
        print(f"    {rank:2d}. {rname}{resid} (idx {node}): BC={cent:.4f} {is_relay}")
        top10.append((f"{rname}{resid}", round(cent, 4)))

    relay_bc = {}
    relay_ranks = {}
    for i in relay_indices:
        resid = resid_map[i]
        rname = ca_resnames[i]
        rank = sum(1 for v in bc.values() if v > bc[i]) + 1
        relay_bc[f"{rname}{resid}"] = round(bc[i], 4)
        relay_ranks[f"{rname}{resid}"] = rank

    relay_in_top10 = sum(1 for i in relay_indices if bc[i] >= bc_sorted[9][1])
    print(f"\n  Relay in top 10: {relay_in_top10}/8")
    print(f"  (6MEO 100ns: 3/8 in top 6)")

    centrality_data = {
        'top_10': top10,
        'relay_residues': relay_bc,
        'relay_ranks': relay_ranks,
        'relay_in_top_10': relay_in_top10,
        'edges_above_0.2': G.number_of_edges(),
        'total_nodes': G.number_of_nodes()
    }
except ImportError:
    print("  NetworkX not available")
    centrality_data = None

results = {
    'pdb_id': PDB_ID,
    'method': 'DCCM from GROMACS covariance (gmx covar), 500ns, dt=50ps',
    'n_ca_atoms': N_ATOMS,
    'relay_motif': target_motif,
    'relay_residues_gromacs': relay_resids,
    'intra_pocket': {
        'all_28_positive': all_positive,
        'mean_correlation': round(mean_intra, 4),
        'min_correlation': round(min(intra_values), 4),
        'max_correlation': round(max(intra_values), 4),
        'sequential_neighbor_mean': round(mean_seq, 4),
        'pairs': {name: round(val, 4) for name, val in intra_pairs}
    },
    'comparison_6meo_500ns': {
        'all_28_positive': True,
        'mean_correlation': 0.79,
        'min_correlation': 0.37,
        'max_correlation': 0.99,
        'sequential_neighbor_mean': 0.97
    }
}

if centrality_data:
    results['betweenness_centrality'] = centrality_data

outfile = f"{PDB_ID}_500ns_dccm.json"
with open(outfile, 'w') as f:
    json.dump(results, f, indent=2)

np.save(f"{PDB_ID}_500ns_dccm.npy", dccm)
print(f"\nSaved: {outfile}, {PDB_ID}_500ns_dccm.npy")
