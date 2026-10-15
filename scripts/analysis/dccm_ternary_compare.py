#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
DCCM Analysis on Full Ternary Complex -- Task 6.1 (Decoupling Proof)

Computes Dynamic Cross-Correlation Matrix (DCCM) and betweenness centrality
on the full ternary complex (gp120+CD4+CCR5, ~887 residues).

Designed to run identically on both:
  - Apo ternary trajectory (20ns, ~292K atoms)
  - Drug-bound ternary trajectory (10ns, ~293K atoms)

Then compare cross-domain correlations to test allosteric silencing.

Usage:
  python3 dccm_ternary_compare.py <tpr_file> <xtc_file> <output_prefix> [stride]

Example:
  python3 dccm_ternary_compare.py md.tpr md.xtc apo 5
  python3 dccm_ternary_compare.py md.tpr md.xtc drug_bound 5
"""

import sys
import os
import json
import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import align
except ImportError:
    print("ERROR: MDAnalysis not installed. Run: pip3 install MDAnalysis")
    sys.exit(1)

try:
    import networkx as nx
except ImportError:
    print("ERROR: networkx not installed. Run: pip3 install networkx")
    sys.exit(1)


def compute_dccm(tpr_file, xtc_file, output_prefix, stride=5):
    print(f"Loading trajectory: {xtc_file}")
    print(f"Topology: {tpr_file}")
    
    u = mda.Universe(tpr_file, xtc_file)
    protein = u.select_atoms("protein")
    ca_atoms = u.select_atoms("protein and name CA")
    
    n_ca = len(ca_atoms)
    n_frames = len(u.trajectory)
    n_used = len(range(0, n_frames, stride))
    
    print(f"System: {len(u.atoms)} atoms, {n_ca} CA atoms, {n_frames} frames")
    print(f"Using every {stride}th frame: {n_used} frames")
    
    # Identify chain boundaries for domain analysis
    chains = {}
    for atom in ca_atoms:
        seg = atom.segid
        if seg not in chains:
            chains[seg] = []
        chains[seg].append(atom.index)
    
    print(f"Chains found: {list(chains.keys())}")
    for c, atoms in chains.items():
        print(f"  Chain {c}: {len(atoms)} CA atoms")
    
    # Identify relay pocket residues (465-472 in PDB numbering, or 372-379 in GROMACS ternary)
    # Try both numbering schemes
    relay_indices = []
    relay_resids = []
    for i, atom in enumerate(ca_atoms):
        resid = atom.resid
        resname = atom.resname
        segid = atom.segid
        # Check for RPGGGDMR motif
        if resname in ['ARG', 'PRO', 'GLY', 'ASP', 'MET']:
            pass  # Could be relay
    
    # More robust: find RPGGGDMR by sequence pattern
    ca_resnames = [atom.resname for atom in ca_atoms]
    ca_resids = [atom.resid for atom in ca_atoms]
    ca_segids = [atom.segid for atom in ca_atoms]
    
    target = ['ARG', 'PRO', 'GLY', 'GLY', 'GLY', 'ASP', 'MET', 'ARG']
    for i in range(len(ca_resnames) - 7):
        if ca_resnames[i:i+8] == target:
            relay_indices = list(range(i, i+8))
            relay_resids = ca_resids[i:i+8]
            relay_chain = ca_segids[i]
            print(f"\nRelay pocket found: indices {i}-{i+7}, resids {relay_resids[0]}-{relay_resids[-1]}, chain {relay_chain}")
            break
    
    if not relay_indices:
        print("WARNING: RPGGGDMR motif not found! Searching by resid ranges...")
        for i, (resid, seg) in enumerate(zip(ca_resids, ca_segids)):
            if 370 <= resid <= 380 or 463 <= resid <= 475:
                relay_indices.append(i)
                relay_resids.append(resid)
        if relay_indices:
            print(f"  Found {len(relay_indices)} relay candidates by resid range")
    
    # Collect CA positions across trajectory
    print(f"\nCollecting CA positions across {n_used} frames...")
    positions = np.zeros((n_used, n_ca, 3))
    
    frame_idx = 0
    for ts_idx, ts in enumerate(u.trajectory):
        if ts_idx % stride != 0:
            continue
        positions[frame_idx] = ca_atoms.positions
        frame_idx += 1
        if frame_idx % 50 == 0:
            print(f"  Frame {frame_idx}/{n_used} (time {ts.time:.1f} ps)")
    
    # Compute mean positions
    mean_pos = positions.mean(axis=0)
    
    # Compute fluctuations
    fluct = positions - mean_pos[np.newaxis, :, :]
    
    # Compute DCCM
    print("Computing DCCM...")
    dccm = np.zeros((n_ca, n_ca))
    
    for i in range(n_ca):
        for j in range(i, n_ca):
            cov_ij = np.mean(np.sum(fluct[:, i, :] * fluct[:, j, :], axis=1))
            var_i = np.mean(np.sum(fluct[:, i, :] ** 2, axis=1))
            var_j = np.mean(np.sum(fluct[:, j, :] ** 2, axis=1))
            
            denom = np.sqrt(var_i * var_j)
            if denom > 1e-10:
                dccm[i, j] = cov_ij / denom
                dccm[j, i] = dccm[i, j]
            else:
                dccm[i, j] = 0.0
                dccm[j, i] = 0.0
    
    print(f"DCCM computed: {n_ca}x{n_ca} matrix")
    print(f"  Range: [{dccm.min():.4f}, {dccm.max():.4f}]")
    print(f"  Mean abs correlation: {np.abs(dccm).mean():.4f}")
    
    # Extract relay pocket sub-matrix
    if relay_indices:
        relay_dccm = dccm[np.ix_(relay_indices, relay_indices)]
        print(f"\nRelay pocket DCCM ({len(relay_indices)}x{len(relay_indices)}):")
        n_pos = np.sum(relay_dccm > 0) - len(relay_indices)  # exclude diagonal
        n_neg = np.sum(relay_dccm < 0)
        n_pairs = len(relay_indices) * (len(relay_indices) - 1) // 2
        print(f"  Positive pairs: {n_pos//2}/{n_pairs}")
        print(f"  Intra-pocket mean: {relay_dccm[np.triu_indices(len(relay_indices), k=1)].mean():.4f}")
        print(f"  Intra-pocket range: [{relay_dccm[np.triu_indices(len(relay_indices), k=1)].min():.4f}, {relay_dccm[np.triu_indices(len(relay_indices), k=1)].max():.4f}]")
    
    # Cross-domain correlations
    print("\nCross-domain correlation analysis:")
    
    # Identify domains by chain
    chain_list = list(chains.keys())
    domain_correlations = {}
    
    for ci, chain_i in enumerate(chain_list):
        for cj, chain_j in enumerate(chain_list):
            if ci >= cj:
                continue
            
            idx_i = [k for k, s in enumerate(ca_segids) if s == chain_i]
            idx_j = [k for k, s in enumerate(ca_segids) if s == chain_j]
            
            sub = dccm[np.ix_(idx_i, idx_j)]
            key = f"{chain_i}-{chain_j}"
            domain_correlations[key] = {
                "mean": float(sub.mean()),
                "mean_abs": float(np.abs(sub).mean()),
                "max": float(sub.max()),
                "min": float(sub.min()),
                "n_positive": int(np.sum(sub > 0)),
                "n_negative": int(np.sum(sub < 0)),
                "n_total": sub.size
            }
            print(f"  {key}: mean={sub.mean():.4f}, |mean|={np.abs(sub).mean():.4f}, max={sub.max():.4f}")
    
    # Relay-to-other-chains correlations
    if relay_indices:
        relay_cross = {}
        for chain in chain_list:
            idx_chain = [k for k, s in enumerate(ca_segids) if s == chain]
            non_relay = [k for k in idx_chain if k not in relay_indices]
            if non_relay:
                sub = dccm[np.ix_(relay_indices, non_relay)]
                relay_cross[f"relay-{chain}"] = {
                    "mean": float(sub.mean()),
                    "mean_abs": float(np.abs(sub).mean()),
                    "max": float(sub.max()),
                    "min": float(sub.min()),
                }
                print(f"  relay-{chain}: mean={sub.mean():.4f}, |mean|={np.abs(sub).mean():.4f}")
    
    # Build network and compute centrality
    print("\nBuilding correlation network...")
    threshold = 0.2
    G = nx.Graph()
    
    for i in range(n_ca):
        G.add_node(i, resid=int(ca_resids[i]), resname=ca_resnames[i], chain=ca_segids[i])
    
    for i in range(n_ca):
        for j in range(i+1, n_ca):
            if abs(dccm[i, j]) > threshold:
                G.add_edge(i, j, weight=abs(dccm[i, j]))
    
    print(f"  Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    
    centrality = nx.betweenness_centrality(G, weight='weight')
    
    # Top 20 most central residues
    top_central = sorted(centrality.items(), key=lambda x: x[1], reverse=True)[:20]
    print("\nTop 20 betweenness centrality:")
    centrality_list = []
    for rank, (node, cent) in enumerate(top_central, 1):
        resid = ca_resids[node]
        resname = ca_resnames[node]
        chain = ca_segids[node]
        is_relay = node in relay_indices if relay_indices else False
        marker = " <<<< RELAY" if is_relay else ""
        print(f"  #{rank}: {resname}{resid} (chain {chain}) = {cent:.4f}{marker}")
        centrality_list.append({
            "rank": rank,
            "resid": int(resid),
            "resname": resname,
            "chain": chain,
            "centrality": float(cent),
            "is_relay": is_relay
        })
    
    # Relay pocket centrality
    if relay_indices:
        relay_centrality = {int(ca_resids[i]): float(centrality.get(i, 0)) for i in relay_indices}
        relay_ranks = {}
        all_sorted = sorted(centrality.items(), key=lambda x: x[1], reverse=True)
        for rank, (node, _) in enumerate(all_sorted, 1):
            if node in relay_indices:
                relay_ranks[int(ca_resids[node])] = rank
        
        print(f"\nRelay pocket centrality ranks (out of {n_ca}):")
        for resid in sorted(relay_centrality.keys()):
            rank = relay_ranks.get(resid, -1)
            print(f"  Res {resid}: centrality {relay_centrality[resid]:.4f}, rank #{rank}")
    
    # Save results
    results = {
        "system": output_prefix,
        "trajectory": os.path.basename(xtc_file),
        "topology": os.path.basename(tpr_file),
        "n_atoms": len(u.atoms),
        "n_ca": n_ca,
        "n_frames_total": n_frames,
        "n_frames_used": n_used,
        "stride": stride,
        "chains": {c: len(v) for c, v in chains.items()},
        "dccm_stats": {
            "mean": float(dccm.mean()),
            "mean_abs": float(np.abs(dccm).mean()),
            "range": [float(dccm.min()), float(dccm.max())],
        },
        "relay_pocket": {
            "indices": relay_indices,
            "resids": [int(r) for r in relay_resids] if relay_resids else [],
            "intra_dccm_mean": float(relay_dccm[np.triu_indices(len(relay_indices), k=1)].mean()) if relay_indices else None,
            "intra_dccm_range": [float(relay_dccm[np.triu_indices(len(relay_indices), k=1)].min()),
                                  float(relay_dccm[np.triu_indices(len(relay_indices), k=1)].max())] if relay_indices else None,
            "centrality": relay_centrality if relay_indices else {},
            "centrality_ranks": relay_ranks if relay_indices else {},
        },
        "cross_domain_correlations": domain_correlations,
        "relay_cross_domain": relay_cross if relay_indices else {},
        "top_centrality": centrality_list,
        "network": {
            "nodes": G.number_of_nodes(),
            "edges": G.number_of_edges(),
            "correlation_threshold": threshold,
        }
    }
    
    out_json = f"{output_prefix}_dccm_ternary.json"
    with open(out_json, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {out_json}")
    
    # Save raw DCCM matrix for comparison
    np.save(f"{output_prefix}_dccm_matrix.npy", dccm)
    print(f"Raw DCCM matrix saved to: {output_prefix}_dccm_matrix.npy")
    
    return results


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    
    tpr = sys.argv[1]
    xtc = sys.argv[2]
    prefix = sys.argv[3]
    stride = int(sys.argv[4]) if len(sys.argv) > 4 else 5
    
    compute_dccm(tpr, xtc, prefix, stride)
