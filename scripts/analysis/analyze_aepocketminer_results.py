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
Analyze AE-PocketMiner GNN predictions: per-residue cryptic pocket
probabilities and attention-weight allosteric coupling for the relay
pocket (Arg465--Arg472, RPGGGDMR).

Reads:
  results/11-gnn-coupling/{PDBID}_chainG-preds.npy         (pocket probabilities)
  results/11-gnn-coupling/{PDBID}_chainG-attention_weights.npy  (N x N attention)
  structures/aepocketminer/{PDBID}_chainG.pdb               (for residue numbering)

Writes:
  results/11-gnn-coupling/ev14_aepocketminer_analysis.json  (combined analysis)

Usage:
    python analyze_aepocketminer_results.py

Requires: numpy (pip install numpy)

AE-PocketMiner: Zhang et al. 2026, bioRxiv 10.64898/2026.05.21.726899
Predecessor:    Meller et al. 2023, Nat Commun 14:1177 (PocketMiner)
"""

import numpy as np
import json
import sys

# ---------- Relay pocket residues (PDB numbering, chain G) ----------
RELAY_RESNUMS = [465, 466, 467, 468, 469, 470, 471, 472]
RELAY_LABELS = ["R465", "P466", "G467", "G468", "G469", "D470", "M471", "R472"]

# ---------- Functional region residue lists (HXB2 numbering) ----------
CD4_BINDING = [
    124, 125, 126, 196, 198, 279, 280, 281, 282, 283,
    365, 366, 367, 368, 370, 371, 374, 375, 376,
    427, 428, 429, 430, 455, 456, 457, 458, 459, 460, 461, 474, 475, 476
]
GP41_CONTACT = [
    473, 474, 475, 476, 477, 478, 479, 480, 481, 482,
    483, 484, 485, 486, 487, 488, 489, 490, 491, 492
]
FOSTEMSAVIR = [375, 377, 382, 384, 426, 427, 428, 429, 430, 431, 432, 433, 434, 435]

STRUCTURES = ["6MEO", "5FYJ", "5U7O", "6MU6", "4TVP"]


def get_resnums_from_pdb(pdb_path):
    """Extract residue numbers from PDB backbone atoms (N, CA, C, O)."""
    resnums = []
    seen = set()
    with open(pdb_path, "r") as f:
        for line in f:
            if line[:4] == "ATOM":
                atom_name = line[12:16].strip()
                if atom_name in ("N", "CA", "C", "O"):
                    resnum = int(line[22:26].strip())
                    if resnum not in seen:
                        seen.add(resnum)
                        resnums.append(resnum)
    return np.array(resnums)


def analyze_structure(struct, results_dir, structures_dir):
    """Analyze one structure's AE-PocketMiner output."""
    pdb_path = os.path.join(structures_dir, "aepocketminer", f"{struct}_chainG.pdb")
    preds_path = os.path.join(results_dir, f"{struct}_chainG-preds.npy")
    attn_path = os.path.join(results_dir, f"{struct}_chainG-attention_weights.npy")

    for p in [pdb_path, preds_path, attn_path]:
        if not os.path.isfile(p):
            print(f"ERROR: {p} not found. Run AE-PocketMiner first.", file=sys.stderr)
            return None

    resnums = get_resnums_from_pdb(pdb_path)
    preds = np.load(preds_path).flatten()
    attn_raw = np.load(attn_path)

    # Average over batch and heads -> (N, N)
    if attn_raw.ndim == 4:
        attn = attn_raw[0].mean(axis=0)
    elif attn_raw.ndim == 3:
        attn = attn_raw.mean(axis=0)
    else:
        attn = attn_raw

    n_res = len(resnums)
    preds = preds[:n_res]
    attn = attn[:n_res, :n_res]

    print(f"\n{'='*70}")
    print(f"  {struct} -- chain G ({n_res} residues)")
    print(f"{'='*70}")

    # 1. Relay pocket probabilities
    print(f"\n  RELAY POCKET PROBABILITIES (cryptic pocket score):")
    relay_probs = []
    for resnum, label in zip(RELAY_RESNUMS, RELAY_LABELS):
        idx = np.where(resnums == resnum)[0]
        if len(idx) > 0:
            prob = float(preds[idx[0]])
            relay_probs.append(prob)
            percentile = 100.0 * np.mean(preds <= prob)
            print(f"    {label} ({resnum}): {prob:.4f}  [{percentile:.0f}th percentile]")
        else:
            print(f"    {label} ({resnum}): NOT FOUND in structure")

    if relay_probs:
        mean_prob = float(np.mean(relay_probs))
        gp120_mean = float(np.mean(preds))
        print(f"    ---")
        print(f"    Relay pocket mean: {mean_prob:.4f}")
        print(f"    gp120 overall mean: {gp120_mean:.4f}")
        print(f"    Ratio: {mean_prob/gp120_mean:.2f}x")

    # 2. Allosteric coupling via attention matrix
    print(f"\n  ALLOSTERIC COUPLING (attention from relay pocket -> rest):")
    relay_indices = []
    for resnum in RELAY_RESNUMS:
        idx = np.where(resnums == resnum)[0]
        if len(idx) > 0:
            relay_indices.append(idx[0])

    cd4_overlap = set()
    gp41_overlap = set()
    top_resnums = np.array([])
    top_scores = np.array([])

    if relay_indices:
        relay_attn = attn[relay_indices, :].mean(axis=0)

        # Exclude relay pocket + immediate neighbors
        excluded = set()
        for ri in relay_indices:
            excluded.update([ri - 1, ri, ri + 1])
        valid_mask = np.ones(n_res, dtype=bool)
        for ei in excluded:
            if 0 <= ei < n_res:
                valid_mask[ei] = False

        valid_scores = relay_attn[valid_mask]
        valid_resnums = resnums[valid_mask]

        top_order = np.argsort(valid_scores)[::-1][:40]
        top_resnums = valid_resnums[top_order]
        top_scores = valid_scores[top_order]

        print(f"    Top 20 residues allosterically coupled to relay pocket:")
        for i in range(min(20, len(top_resnums))):
            resnum = int(top_resnums[i])
            score = float(top_scores[i])
            region = ""
            if resnum in CD4_BINDING:
                region = " << CD4 BINDING SITE"
            elif resnum in GP41_CONTACT:
                region = " << gp41 CONTACT REGION"
            elif resnum in FOSTEMSAVIR:
                region = " << FOSTEMSAVIR SITE"
            print(f"      #{i+1:2d}: res {resnum:4d}  score {score:.6f}{region}")

        top40_set = set(int(r) for r in top_resnums)
        cd4_overlap = top40_set & set(CD4_BINDING)
        gp41_overlap = top40_set & set(GP41_CONTACT)
        fost_overlap = top40_set & set(FOSTEMSAVIR)

        print(f"\n    FUNCTIONAL REGION OVERLAP (top 40):")
        print(f"      CD4 binding site: {len(cd4_overlap)} residues  {sorted(cd4_overlap) if cd4_overlap else '(none)'}")
        print(f"      gp41 contact:     {len(gp41_overlap)} residues  {sorted(gp41_overlap) if gp41_overlap else '(none)'}")
        print(f"      Fostemsavir site: {len(fost_overlap)} residues  {sorted(fost_overlap) if fost_overlap else '(none)'}")

    return {
        "structure": struct,
        "n_residues": n_res,
        "relay_pocket_probs": {label: float(p) for label, p in zip(RELAY_LABELS, relay_probs)} if relay_probs else {},
        "relay_pocket_mean_prob": float(np.mean(relay_probs)) if relay_probs else None,
        "gp120_mean_prob": float(np.mean(preds)),
        "top40_coupled_residues": [int(r) for r in top_resnums] if relay_indices else [],
        "top40_scores": [float(s) for s in top_scores] if relay_indices else [],
        "cd4_overlap_count": len(cd4_overlap) if relay_indices else 0,
        "gp41_overlap_count": len(gp41_overlap) if relay_indices else 0,
        "cd4_overlap_residues": sorted(cd4_overlap) if relay_indices and cd4_overlap else [],
        "gp41_overlap_residues": sorted(gp41_overlap) if relay_indices and gp41_overlap else [],
    }


def main():
    results_dir = os.path.join(DATA_ROOT, "results", "11-gnn-coupling")
    structures_dir = os.path.join(DATA_ROOT, "structures")

    all_results = {}
    for struct in STRUCTURES:
        result = analyze_structure(struct, results_dir, structures_dir)
        if result is not None:
            all_results[struct] = result

    if not all_results:
        print("ERROR: No structures analyzed. Check file paths.", file=sys.stderr)
        sys.exit(1)

    out_path = os.path.join(results_dir, "ev14_aepocketminer_analysis.json")
    with open(out_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n\nResults saved to {out_path}")
    print(f"Structures analyzed: {len(all_results)}/5")


if __name__ == "__main__":
    main()
