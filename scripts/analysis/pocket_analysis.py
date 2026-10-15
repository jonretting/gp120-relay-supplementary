#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Pocket-focused MD analysis using internal distances and per-residue metrics.
These metrics are insensitive to PBC artifacts and global structural drift.
"""
import os, sys, json, subprocess
import numpy as np
from pathlib import Path
from datetime import datetime
from itertools import combinations

SIMULATION = sys.argv[1] if len(sys.argv) > 1 else "unknown"
WORK_DIR = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(".")
TPR = sys.argv[3] if len(sys.argv) > 3 else "md.tpr"
XTC = sys.argv[4] if len(sys.argv) > 4 else "md.xtc"

OUTDIR = WORK_DIR / "validation"
OUTDIR.mkdir(exist_ok=True)

def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def run(cmd, input_str=None):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                       input=input_str, timeout=600, cwd=str(WORK_DIR))
    return r.stdout + r.stderr

def parse_xvg(filepath):
    data = []
    with open(filepath) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    data.append([float(x) for x in parts])
                except ValueError:
                    continue
    return data

def extract_ca_coords(pdb_file, target_resnums):
    """Extract CA coordinates for specific residue numbers."""
    coords = {}
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") and " CA " in line:
                resnum = int(line[22:26].strip())
                if resnum in target_resnums:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    resname = line[17:20].strip()
                    chain = line[21]
                    key = (resnum, chain)
                    if key not in coords:
                        coords[key] = {"xyz": np.array([x, y, z]), "resname": resname, "chain": chain}
    return coords

log(f"{'='*70}")
log(f"POCKET-FOCUSED ANALYSIS: {SIMULATION}")
log(f"{'='*70}")

results = {"simulation": SIMULATION, "timestamp": datetime.now().isoformat()}

# Determine which numbering scheme the system uses
# Extract first frame PDB
first_pdb = OUTDIR / "first_frame.pdb"
if not first_pdb.exists():
    clean_xtc = OUTDIR / "clean.xtc"
    if not clean_xtc.exists():
        nj_xtc = str(OUTDIR / "nojump.xtc")
        run(f"echo '0' | gmx trjconv -s {TPR} -f {XTC} -o {nj_xtc} -pbc nojump", "0\n")
        run(f"echo '1\n0' | gmx trjconv -s {TPR} -f {nj_xtc} -o {clean_xtc} -center -pbc mol", "1\n0\n")
    run(f"echo '1' | gmx trjconv -s {TPR} -f {clean_xtc} -o {first_pdb} -dump 0 -pbc mol", "1\n")

# Try each numbering scheme
pocket_numberings = {
    "6MEO": list(range(465, 473)),
    "SOSIP": list(range(469, 477)),
    "reduced": list(range(69, 77)),
}

if first_pdb.exists():
    found_scheme = None
    for scheme, resnums in pocket_numberings.items():
        coords = extract_ca_coords(str(first_pdb), set(resnums))
        # Filter to get unique residues (first occurrence per resnum)
        unique_coords = {}
        for (rn, ch), v in coords.items():
            if rn not in unique_coords:
                unique_coords[rn] = v
        
        if len(unique_coords) >= 6:
            found_scheme = scheme
            log(f"  Found pocket with {scheme} numbering: {len(unique_coords)}/8 residues")
            break
    
    if not found_scheme:
        log("  WARNING: Could not identify pocket residues in first frame!")
        # Try to list all CA residues
        all_ca = {}
        with open(first_pdb) as f:
            for line in f:
                if line.startswith("ATOM") and " CA " in line:
                    rn = int(line[22:26].strip())
                    ch = line[21]
                    rname = line[17:20].strip()
                    all_ca[(rn, ch)] = rname
        
        # Show residue ranges per chain
        chains = {}
        for (rn, ch), rname in sorted(all_ca.items()):
            if ch not in chains:
                chains[ch] = []
            chains[ch].append(rn)
        
        for ch, rns in sorted(chains.items()):
            log(f"    Chain {ch}: residues {min(rns)}-{max(rns)} ({len(rns)} residues)")
        
        results["chain_mapping"] = {ch: {"min": min(rns), "max": max(rns), "count": len(rns)} for ch, rns in chains.items()}
        
        # Try to find RPGGGDMR sequence
        log("  Searching for RPGGGDMR motif...")
        aa3to1 = {"ARG":"R","PRO":"P","GLY":"G","ASP":"D","MET":"M","ALA":"A",
                   "LEU":"L","VAL":"V","ILE":"I","PHE":"F","TRP":"W","TYR":"Y",
                   "HIS":"H","LYS":"K","SER":"S","THR":"T","ASN":"N","GLN":"Q",
                   "CYS":"C","GLU":"E"}
        
        for ch, rns in sorted(chains.items()):
            seq = ""
            rn_list = sorted(rns)
            for rn in rn_list:
                seq += aa3to1.get(all_ca.get((rn, ch), "?"), "?")
            
            # Look for RPGGGDM or RPGGGDMR or RPAGGDMR
            for motif in ["RPGGGDM", "RPAGGDM"]:
                idx = seq.find(motif)
                if idx >= 0:
                    pocket_start = rn_list[idx]
                    pocket_end = rn_list[min(idx + 7, len(rn_list)-1)]
                    pocket_resnums = rn_list[idx:idx+8]
                    log(f"    Found {motif} in chain {ch} at residues {pocket_start}-{pocket_end}")
                    found_scheme = f"chain_{ch}"
                    pocket_numberings[found_scheme] = pocket_resnums
                    break
            if found_scheme:
                break

# Create pocket index for GROMACS
if found_scheme:
    pocket_resnums = pocket_numberings[found_scheme]
    
    # Make custom index with pocket residues
    resid_selection = " | ".join([f"r {r}" for r in pocket_resnums])
    ndx_cmd = f"echo '{resid_selection}\nq' | gmx make_ndx -f {TPR} -o {OUTDIR}/pocket.ndx"
    output = run(ndx_cmd)
    
    # Actually, the above probably won't work well. Let me use a different approach.
    # Extract frames and compute distances in Python
    
    log("")
    log(f"Extracting pocket trajectory (scheme: {found_scheme}, residues: {pocket_resnums})")
    
    # Extract all protein frames as PDB trajectory (every 10th frame for speed)
    # Then parse CA coordinates from each frame
    
    traj_pdb = str(OUTDIR / "protein_traj.pdb")
    clean_xtc = OUTDIR / "clean.xtc"
    if clean_xtc.exists():
        run(f"echo '3' | gmx trjconv -s {TPR} -f {clean_xtc} -o {traj_pdb} -skip 10", "3\n")
    else:
        run(f"echo '3' | gmx trjconv -s {TPR} -f {XTC} -o {traj_pdb} -pbc nojump -skip 10", "3\n")
    
    log("Parsing trajectory PDB...")
    
    if os.path.exists(traj_pdb):
        pocket_set = set(pocket_resnums)
        
        all_frames = []
        current_frame = {}
        frame_num = 0
        
        with open(traj_pdb) as f:
            for line in f:
                if line.startswith("ENDMDL") or (line.startswith("END") and current_frame):
                    if current_frame:
                        all_frames.append(current_frame)
                    current_frame = {}
                    frame_num += 1
                elif line.startswith("ATOM") and " CA " in line:
                    resnum = int(line[22:26].strip())
                    chain = line[21]
                    if resnum in pocket_set:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        key = (resnum, chain)
                        if key not in current_frame:
                            current_frame[key] = np.array([x, y, z])
        
        if current_frame:
            all_frames.append(current_frame)
        
        log(f"  Parsed {len(all_frames)} frames, pocket residues per frame: {len(all_frames[0]) if all_frames else 0}")
        
        if all_frames and len(all_frames[0]) >= 6:
            # For multi-chain systems, pick the chain with most pocket residues
            first = all_frames[0]
            chain_counts = {}
            for (rn, ch) in first:
                chain_counts[ch] = chain_counts.get(ch, 0) + 1
            
            best_chain = max(chain_counts, key=chain_counts.get)
            log(f"  Using chain {best_chain} ({chain_counts[best_chain]} pocket residues)")
            
            # Filter to best chain
            pocket_keys = sorted([k for k in first if k[1] == best_chain])
            n_pocket = len(pocket_keys)
            
            if n_pocket >= 6:
                # Extract positions array: (n_frames, n_residues, 3)
                positions = []
                for frame in all_frames:
                    frame_pos = []
                    valid = True
                    for k in pocket_keys:
                        if k in frame:
                            frame_pos.append(frame[k])
                        else:
                            valid = False
                            break
                    if valid:
                        positions.append(np.array(frame_pos))
                
                positions = np.array(positions)
                n_frames = len(positions)
                log(f"  Valid frames: {n_frames}, pocket atoms: {n_pocket}")
                
                # === Analysis 1: Intra-pocket distances (PBC-independent) ===
                log("\n--- Analysis 1: Intra-pocket CA-CA distances ---")
                pairs = list(combinations(range(n_pocket), 2))
                
                dist_matrix_first = np.zeros((n_pocket, n_pocket))
                for i, j in pairs:
                    d0 = np.linalg.norm(positions[0, i] - positions[0, j])
                    dist_matrix_first[i, j] = d0
                    dist_matrix_first[j, i] = d0
                
                # Track each pair distance over time
                pair_data = {}
                for pi, (i, j) in enumerate(pairs):
                    dists = [np.linalg.norm(positions[t, i] - positions[t, j]) for t in range(n_frames)]
                    rn_i, _ = pocket_keys[i]
                    rn_j, _ = pocket_keys[j]
                    pair_data[f"{rn_i}-{rn_j}"] = {
                        "initial_A": round(float(dists[0]), 2),
                        "mean_A": round(float(np.mean(dists)), 2),
                        "std_A": round(float(np.std(dists)), 2),
                        "final_A": round(float(dists[-1]), 2),
                        "max_deviation_A": round(float(max(abs(d - dists[0]) for d in dists)), 2),
                    }
                
                # Summary stats
                all_means = [v["mean_A"] for v in pair_data.values()]
                all_stds = [v["std_A"] for v in pair_data.values()]
                all_devs = [v["max_deviation_A"] for v in pair_data.values()]
                
                results["intra_pocket_distances"] = {
                    "n_pairs": len(pairs),
                    "avg_mean_dist_A": round(float(np.mean(all_means)), 2),
                    "avg_std_A": round(float(np.mean(all_stds)), 2),
                    "max_deviation_A": round(float(np.max(all_devs)), 2),
                    "pairs": pair_data,
                }
                
                log(f"  {len(pairs)} CA-CA pairs: avg distance = {np.mean(all_means):.2f} A")
                log(f"  Avg fluctuation (std): {np.mean(all_stds):.2f} A")
                log(f"  Max single-pair deviation: {np.max(all_devs):.2f} A")
                
                # === Analysis 2: Pocket internal RMSD (self-aligned) ===
                log("\n--- Analysis 2: Pocket internal RMSD (self-aligned) ---")
                
                ref_pos = positions[0]
                ref_center = ref_pos.mean(axis=0)
                ref_centered = ref_pos - ref_center
                
                pocket_rmsds = []
                for t in range(n_frames):
                    pos = positions[t]
                    pos_center = pos.mean(axis=0)
                    pos_centered = pos - pos_center
                    
                    # Kabsch alignment
                    H = pos_centered.T @ ref_centered
                    U, S, Vt = np.linalg.svd(H)
                    d = np.linalg.det(Vt.T @ U.T)
                    sign_matrix = np.eye(3)
                    sign_matrix[2, 2] = np.sign(d)
                    R = Vt.T @ sign_matrix @ U.T
                    
                    aligned = pos_centered @ R.T
                    rmsd = np.sqrt(np.mean(np.sum((aligned - ref_centered)**2, axis=1)))
                    pocket_rmsds.append(float(rmsd))
                
                q4_pocket = pocket_rmsds[-len(pocket_rmsds)//4:]
                results["pocket_internal_rmsd"] = {
                    "mean_A": round(float(np.mean(pocket_rmsds)), 3),
                    "std_A": round(float(np.std(pocket_rmsds)), 3),
                    "max_A": round(float(np.max(pocket_rmsds)), 3),
                    "final_A": round(float(pocket_rmsds[-1]), 3),
                    "q4_mean_A": round(float(np.mean(q4_pocket)), 3),
                    "q4_std_A": round(float(np.std(q4_pocket)), 3),
                    "n_frames": n_frames,
                }
                
                log(f"  Pocket RMSD (Kabsch-aligned): {np.mean(pocket_rmsds):.3f} +/- {np.std(pocket_rmsds):.3f} A")
                log(f"  Last quarter: {np.mean(q4_pocket):.3f} +/- {np.std(q4_pocket):.3f} A")
                log(f"  Max pocket RMSD: {np.max(pocket_rmsds):.3f} A")
                
                # Time series sample
                indices = np.linspace(0, n_frames-1, 10, dtype=int)
                results["pocket_rmsd_timeseries"] = [
                    {"frame": int(i), "rmsd_A": round(pocket_rmsds[i], 3)} for i in indices
                ]
                for i in indices:
                    log(f"    Frame {i}: {pocket_rmsds[i]:.3f} A")
                
                # === Analysis 3: Pocket compactness (radius of gyration) ===
                log("\n--- Analysis 3: Pocket compactness ---")
                pocket_rg = []
                for t in range(n_frames):
                    pos = positions[t]
                    center = pos.mean(axis=0)
                    rg = np.sqrt(np.mean(np.sum((pos - center)**2, axis=1)))
                    pocket_rg.append(float(rg))
                
                results["pocket_compactness"] = {
                    "rg_mean_A": round(float(np.mean(pocket_rg)), 2),
                    "rg_std_A": round(float(np.std(pocket_rg)), 2),
                    "rg_first_A": round(float(pocket_rg[0]), 2),
                    "rg_last_A": round(float(pocket_rg[-1]), 2),
                    "rg_change_pct": round(float((pocket_rg[-1] - pocket_rg[0]) / pocket_rg[0] * 100), 1),
                }
                
                log(f"  Pocket Rg: {np.mean(pocket_rg):.2f} +/- {np.std(pocket_rg):.2f} A")
                log(f"  First: {pocket_rg[0]:.2f} A -> Last: {pocket_rg[-1]:.2f} A ({results['pocket_compactness']['rg_change_pct']:+.1f}%)")
                
                # === Analysis 4: Per-residue displacement from start ===
                log("\n--- Analysis 4: Per-residue displacement ---")
                displacements = {}
                for ri, key in enumerate(pocket_keys):
                    rn, ch = key
                    d_over_time = [np.linalg.norm(positions[t, ri] - positions[0, ri]) for t in range(n_frames)]
                    displacements[int(rn)] = {
                        "mean_A": round(float(np.mean(d_over_time)), 2),
                        "max_A": round(float(np.max(d_over_time)), 2),
                        "final_A": round(float(d_over_time[-1]), 2),
                    }
                    log(f"  Res {rn}: mean disp = {np.mean(d_over_time):.2f} A, max = {np.max(d_over_time):.2f} A")
                
                results["per_residue_displacement"] = displacements

# Clean up large temp files
traj_pdb = OUTDIR / "protein_traj.pdb"
if traj_pdb.exists():
    size_mb = traj_pdb.stat().st_size / (1024*1024)
    log(f"\nCleaning up trajectory PDB ({size_mb:.0f} MB)")
    os.remove(traj_pdb)

# Save
with open(OUTDIR / "pocket_analysis.json", "w") as f:
    json.dump(results, f, indent=2, default=str)

log("")
log("="*70)
log("POCKET ANALYSIS COMPLETE")
log(f"Results: {OUTDIR / 'pocket_analysis.json'}")
log("="*70)
