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
Allosteric Communication Pathway / Dynamic Network Analysis
Computes residue-residue cross-correlations from the MD trajectory
to map communication pathways through the relay pocket.
"""

import numpy as np
import MDAnalysis as mda
import json, os, time

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

print("=== Allosteric Network / Dynamic Cross-Correlation Analysis ===")
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")

TPR = os.path.join(DATA_ROOT, "topology", "md.tpr")
XTC = os.path.join(DATA_ROOT, "topology", "md.xtc")  # Note: large trajectory not in package
OUT_DIR = os.path.join(DATA_ROOT, "results", "05-elastic-network")
os.makedirs(OUT_DIR, exist_ok=True)

POCKET_RESIDUES = list(range(69, 77))
CONTACT_CUTOFF = 4.5
SKIP_FRAMES = 5

RESNAMES = {69:'ARG', 70:'PRO', 71:'GLY', 72:'GLY', 73:'GLY', 74:'ASP', 75:'MET', 76:'ARG'}
ORIG_MAP = {69:465, 70:466, 71:467, 72:468, 73:469, 74:470, 75:471, 76:472}

def to_py(val):
    """Convert numpy types to Python natives for JSON."""
    if hasattr(val, 'item'):
        return val.item()
    return val

print(f"\nLoading universe: {TPR} + {XTC}")
u = mda.Universe(TPR, XTC)
print(f"  Atoms: {len(u.atoms)}")
print(f"  Residues: {len(u.residues)}")
print(f"  Frames: {len(u.trajectory)}")

ca_atoms = u.select_atoms("name CA")
n_residues = len(ca_atoms)
print(f"  C-alpha atoms: {n_residues}")

print("\n--- Step 1: Computing mean C-alpha positions ---")
positions_all = []
frame_count = 0
for ts in u.trajectory[::SKIP_FRAMES]:
    positions_all.append(ca_atoms.positions.copy())
    frame_count += 1
    if frame_count % 200 == 0:
        print(f"  Frame {frame_count}...")

positions_all = np.array(positions_all)
mean_positions = positions_all.mean(axis=0)
displacements = positions_all - mean_positions
print(f"  Processed {frame_count} frames")

print("\n--- Step 2: Computing DCCM ---")
norms = np.sqrt(np.sum(np.mean(displacements**2, axis=0), axis=1))

dccm = np.zeros((n_residues, n_residues))
for i in range(n_residues):
    if i % 20 == 0:
        print(f"  Row {i}/{n_residues}...")
    for j in range(i, n_residues):
        cij = np.mean(np.sum(displacements[:, i, :] * displacements[:, j, :], axis=1))
        norm = norms[i] * norms[j]
        if norm > 0:
            dccm[i, j] = cij / norm
        dccm[j, i] = dccm[i, j]

print(f"  DCCM shape: {dccm.shape}")

resid_map = [int(r.resid) for r in ca_atoms.residues]

print("\n--- Step 3: Relay pocket cross-correlations ---")
pocket_indices = [i for i, rid in enumerate(resid_map) if rid in POCKET_RESIDUES]
print(f"  Pocket residue indices: {pocket_indices}")
print(f"  Pocket resids: {[resid_map[i] for i in pocket_indices]}")

print("\n  Intra-pocket correlations:")
for i_idx in pocket_indices:
    for j_idx in pocket_indices:
        if i_idx >= j_idx:
            continue
        ri = resid_map[i_idx]
        rj = resid_map[j_idx]
        c = dccm[i_idx, j_idx]
        label = "corr" if c > 0 else "anti"
        if abs(c) > 0.3:
            orig_i = ORIG_MAP.get(ri, ri)
            orig_j = ORIG_MAP.get(rj, rj)
            print(f"    {RESNAMES.get(ri,'?')}{orig_i} <-> {RESNAMES.get(rj,'?')}{orig_j}: {c:+.3f} ({label})")

print("\n  Top correlations: non-pocket -> pocket")
non_pocket_corr = []
for i_idx in range(n_residues):
    ri = resid_map[i_idx]
    if ri in POCKET_RESIDUES:
        continue
    for j_idx in pocket_indices:
        rj = resid_map[j_idx]
        c = float(dccm[i_idx, j_idx])
        non_pocket_corr.append((ri, rj, c))

non_pocket_corr.sort(key=lambda x: -abs(x[2]))
print(f"  {'Residue':>10} {'Pocket Res':>12} {'Correlation':>12} {'Type':>6}")
print("  " + "-" * 44)
for ri, rj, c in non_pocket_corr[:20]:
    label = "+" if c > 0 else "-"
    orig_j = ORIG_MAP.get(rj, rj)
    print(f"  {ri:>10} {RESNAMES.get(rj,'?')}{orig_j:>8} {c:>+12.3f} {label:>6}")

print("\n--- Step 4: Contact network (persistent contacts) ---")
contact_matrix = np.zeros((n_residues, n_residues))
for frame_idx in range(frame_count):
    pos = positions_all[frame_idx]
    for p_idx in pocket_indices:
        dists = np.linalg.norm(pos - pos[p_idx], axis=1)
        close = np.where(dists < CONTACT_CUTOFF)[0]
        for c_idx in close:
            if c_idx != p_idx:
                contact_matrix[p_idx, c_idx] += 1
                contact_matrix[c_idx, p_idx] += 1

contact_freq = contact_matrix / frame_count

print("\n  Persistent contacts (>50% occupancy) to relay pocket:")
persistent = []
for p_idx in pocket_indices:
    rp = resid_map[p_idx]
    for j in range(n_residues):
        rj = resid_map[j]
        if rj in POCKET_RESIDUES:
            continue
        freq = float(contact_freq[p_idx, j])
        if freq > 0.5:
            persistent.append((rp, rj, freq))

persistent.sort(key=lambda x: -x[2])
for rp, rj, freq in persistent[:15]:
    orig_p = ORIG_MAP.get(rp, rp)
    print(f"    {RESNAMES.get(rp,'?')}{orig_p} <-> Res{rj}: {freq:.1%} occupancy")

print("\n--- Step 5: Communication pathway analysis ---")
try:
    import networkx as nx

    G = nx.Graph()
    for i in range(n_residues):
        G.add_node(resid_map[i])

    for i in range(n_residues):
        for j in range(i+1, n_residues):
            c = abs(dccm[i, j])
            if c > 0.2:
                weight = 1.0 / c if c > 0 else 100
                G.add_edge(resid_map[i], resid_map[j], weight=weight, correlation=float(dccm[i, j]))

    print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    bc = nx.betweenness_centrality(G, weight='weight')
    print("\n  Betweenness centrality of pocket residues:")
    for rid in POCKET_RESIDUES:
        if rid in bc:
            orig = ORIG_MAP.get(rid, rid)
            print(f"    {RESNAMES.get(rid,'?')}{orig}: {bc[rid]:.4f}")

    top_bc = sorted(bc.items(), key=lambda x: -x[1])[:10]
    print("\n  Top 10 most central residues in network:")
    for rid, bv in top_bc:
        marker = " <-- POCKET" if rid in POCKET_RESIDUES else ""
        print(f"    Res{rid}: {bv:.4f}{marker}")

    print("\n  Representative communication paths through pocket:")
    for start_rid, end_rid in [(1, 50), (10, 90), (30, 105), (1, 105)]:
        if start_rid in G and end_rid in G:
            try:
                path = nx.shortest_path(G, start_rid, end_rid, weight='weight')
                pocket_in_path = [r for r in path if r in POCKET_RESIDUES]
                path_str = " -> ".join(str(r) for r in path)
                print(f"    {start_rid} -> {end_rid} (len={len(path)}): {path_str}")
                if pocket_in_path:
                    print(f"      Passes through pocket: {[ORIG_MAP.get(r,r) for r in pocket_in_path]}")
            except nx.NetworkXNoPath:
                print(f"    {start_rid} -> {end_rid}: No path found")

    pathway_data = {
        "betweenness_centrality": {str(int(k)): round(float(v), 6) for k, v in bc.items()},
        "top_10_central": [[int(rid), round(float(bv), 6)] for rid, bv in top_bc]
    }
except ImportError:
    print("  NetworkX not available, skipping pathway analysis")
    pathway_data = {}

def convert_numpy(obj):
    """Recursively convert all numpy types to native Python."""
    if isinstance(obj, dict):
        return {str(k): convert_numpy(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [convert_numpy(v) for v in obj]
    if hasattr(obj, 'item'):
        return obj.item()
    return obj

print("\n--- Saving results ---")
output = {
    "analysis": "Allosteric Network / Dynamic Cross-Correlation",
    "frames_analyzed": int(frame_count),
    "skip_frames": int(SKIP_FRAMES),
    "n_residues": int(n_residues),
    "pocket_residues_tpr": [int(r) for r in POCKET_RESIDUES],
    "pocket_residues_original": [465, 466, 467, 468, 469, 470, 471, 472],
    "resid_mapping": {str(int(k)): int(v) for k, v in ORIG_MAP.items()},
    "intra_pocket_correlations": [
        {"res_i": int(ORIG_MAP.get(resid_map[i], resid_map[i])),
         "res_j": int(ORIG_MAP.get(resid_map[j], resid_map[j])),
         "correlation": round(float(dccm[i, j]), 4)}
        for i in pocket_indices for j in pocket_indices if i < j
    ],
    "top_external_correlations": [
        {"external_res": int(ri), "pocket_res": int(ORIG_MAP.get(rj, rj)), "correlation": round(float(c), 4)}
        for ri, rj, c in non_pocket_corr[:30]
    ],
    "persistent_contacts": [
        {"pocket_res": int(ORIG_MAP.get(rp, rp)), "partner_res": int(rj), "occupancy": round(float(freq), 4)}
        for rp, rj, freq in persistent[:20]
    ],
    "pathway_data": pathway_data,
}

output = convert_numpy(output)

outfile = os.path.join(OUT_DIR, "allosteric_network.json")
with open(outfile, 'w') as f:
    json.dump(output, f, indent=2)
print(f"Results saved to: {outfile}")

dccm_pocket = {}
for i_idx in pocket_indices:
    ri = int(ORIG_MAP.get(resid_map[i_idx], resid_map[i_idx]))
    row = {}
    for j_idx in pocket_indices:
        rj = int(ORIG_MAP.get(resid_map[j_idx], resid_map[j_idx]))
        row[str(rj)] = round(float(dccm[i_idx, j_idx]), 4)
    dccm_pocket[str(ri)] = row

dccm_file = os.path.join(OUT_DIR, "dccm_pocket.json")
with open(dccm_file, 'w') as f:
    json.dump(dccm_pocket, f, indent=2)
print(f"Pocket DCCM saved to: {dccm_file}")

print(f"\n=== ALLOSTERIC NETWORK ANALYSIS COMPLETE ===")
print(f"Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}")
