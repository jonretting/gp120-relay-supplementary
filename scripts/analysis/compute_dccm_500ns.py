#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Compute DCCM (Dynamic Cross-Correlation Matrix) from GROMACS covariance matrix.

The covariance matrix is 3N x 3N (N = number of Cα atoms = 105).
Each element C[3i+a, 3j+b] represents the covariance between
coordinate 'a' of atom 'i' and coordinate 'b' of atom 'j'.

DCCM_ij = Σ_a C[3i+a, 3j+a] / sqrt(Σ_a C[3i+a, 3i+a] * Σ_a C[3j+a, 3j+a])
         = Tr(C_ij) / sqrt(Tr(C_ii) * Tr(C_jj))
"""
import numpy as np
import json
import statistics

COVAR_FILE = "./covar_ascii.dat"
N_ATOMS = 105
N_COORDS = N_ATOMS * 3  # 315

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

assert row == N_COORDS and col == 0, f"Parsed {row} rows, {col} leftover cols"
print(f"Covariance matrix loaded: {covar.shape}")
print(f"Trace: {np.trace(covar):.3f} nm^2")

dccm = np.zeros((N_ATOMS, N_ATOMS))
for i in range(N_ATOMS):
    for j in range(N_ATOMS):
        cov_ij = sum(covar[3*i+a, 3*j+a] for a in range(3))
        cov_ii = sum(covar[3*i+a, 3*i+a] for a in range(3))
        cov_jj = sum(covar[3*j+a, 3*j+a] for a in range(3))
        denom = np.sqrt(cov_ii * cov_jj)
        if denom > 0:
            dccm[i, j] = cov_ij / denom
        else:
            dccm[i, j] = 0.0

print(f"\nDCCM matrix: {dccm.shape}")
print(f"Diagonal check (should all be 1.0): {dccm[0,0]:.4f}, {dccm[50,50]:.4f}, {dccm[104,104]:.4f}")
print(f"Range: [{dccm.min():.4f}, {dccm.max():.4f}]")

# Relay pocket residues: 465-472 in PDB numbering
# Need to find their indices in the 105-residue system
# The residue numbering is preserved from PDB, so we need to identify which
# of the 105 Cα atoms correspond to residues 465-472

# Read residue numbers from the GRO file to build the mapping
resid_map = {}
with open("./md_500ns.gro") as f:
    lines = f.readlines()
    ca_idx = 0
    for line in lines[2:-1]:  # skip header and box line
        resid = int(line[:5].strip())
        atom_name = line[10:15].strip()
        if atom_name == "CA":
            resid_map[ca_idx] = resid
            ca_idx += 1

print(f"\nTotal CA atoms mapped: {len(resid_map)}")

relay_resids = [465, 466, 467, 468, 469, 470, 471, 472]
relay_names = ['ARG465', 'PRO466', 'GLY467', 'GLY468', 'GLY469', 'ASP470', 'MET471', 'ARG472']
relay_indices = [idx for idx, resid in resid_map.items() if resid in relay_resids]
print(f"Relay pocket CA indices: {relay_indices}")
print(f"Relay pocket residues: {[resid_map[i] for i in relay_indices]}")

# Intra-pocket correlations (28 unique pairs)
print("\n=== Intra-Pocket Cross-Correlations (500ns DCCM) ===")
intra_pairs = []
all_positive = True
for a_idx, a in enumerate(relay_indices):
    for b_idx, b in enumerate(relay_indices):
        if b_idx > a_idx:
            c = dccm[a, b]
            pair_name = f"{relay_names[a_idx]}-{relay_names[b_idx]}"
            intra_pairs.append((pair_name, c))
            sign = "+" if c > 0 else "-"
            if c <= 0:
                all_positive = False
            print(f"  {pair_name:25s}: {sign}{abs(c):.4f}")

intra_values = [p[1] for p in intra_pairs]
print(f"\n  All 28 pairs positive: {all_positive}")
print(f"  Mean intra-pocket correlation: {statistics.mean(intra_values):+.4f}")
print(f"  Range: {min(intra_values):+.4f} to {max(intra_values):+.4f}")
print(f"  Sequential neighbor mean: ", end="")
seq_pairs = [dccm[relay_indices[i], relay_indices[i+1]] for i in range(len(relay_indices)-1)]
print(f"{statistics.mean(seq_pairs):+.4f}")

# Compare to 100ns values from paper
print("\n--- Comparison to 100ns DCCM ---")
print(f"  100ns: All 28 intra-pocket pairs positive (C_ij = +0.08 to +0.82)")
print(f"  500ns: {'All 28 positive' if all_positive else 'NOT all positive'} (C_ij = {min(intra_values):+.4f} to {max(intra_values):+.4f})")
print(f"  100ns sequential neighbor mean: >0.75")
print(f"  500ns sequential neighbor mean: {statistics.mean(seq_pairs):+.4f}")

# Betweenness centrality analysis
print("\n=== Betweenness Centrality (500ns) ===")
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
    
    print("  Top 10 most central residues:")
    for rank, (node, cent) in enumerate(bc_sorted[:10], 1):
        resid = resid_map.get(node, '?')
        is_relay = '*' if resid in relay_resids else ' '
        print(f"    {rank:2d}. Residue {resid} (idx {node}): BC = {cent:.4f} {is_relay}")
    
    print("\n  Relay pocket centrality:")
    relay_bc = [(resid_map[i], bc[i]) for i in relay_indices]
    for resid, cent in relay_bc:
        rank = sum(1 for v in bc.values() if v > cent) + 1
        print(f"    {resid}: BC = {cent:.4f} (rank {rank}/{N_ATOMS})")
    
    relay_in_top6 = sum(1 for i in relay_indices if bc[i] >= bc_sorted[5][1])
    print(f"\n  Relay residues in top 6: {relay_in_top6}/8")
    print(f"  100ns comparison: 3 relay glycines in top 6")
except ImportError:
    print("  NetworkX not available, skipping centrality analysis")

# Save results
results = {
    'method': 'DCCM from GROMACS covariance matrix (gmx covar)',
    'trajectory': '500ns reduced system (10,001 frames, dt=50ps)',
    'atoms': f'{N_ATOMS} C-alpha atoms',
    'intra_pocket': {
        'all_28_positive': all_positive,
        'mean_correlation': round(statistics.mean(intra_values), 4),
        'min_correlation': round(min(intra_values), 4),
        'max_correlation': round(max(intra_values), 4),
        'sequential_neighbor_mean': round(statistics.mean(seq_pairs), 4),
        'pairs': {name: round(val, 4) for name, val in intra_pairs}
    },
    'comparison_100ns': {
        'all_28_positive_100ns': True,
        'range_100ns': '+0.08 to +0.82',
        'sequential_neighbor_mean_100ns': '>0.75'
    }
}

# Add centrality if available
try:
    results['betweenness_centrality'] = {
        'top_10': [(resid_map.get(n, n), round(c, 4)) for n, c in bc_sorted[:10]],
        'relay_residues': {str(resid_map[i]): round(bc[i], 4) for i in relay_indices},
        'relay_in_top_6': relay_in_top6,
        'edges_above_0.2': G.number_of_edges(),
        'total_nodes': G.number_of_nodes()
    }
except:
    pass

np.save("./dccm_500ns.npy", dccm)

with open("./dccm_500ns_analysis.json", 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nDCCM matrix saved to dccm_500ns.npy")
print(f"Analysis saved to dccm_500ns_analysis.json")
