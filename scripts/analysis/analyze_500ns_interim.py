import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Interim analysis of 500ns trajectory (~249 ns so far).

Analyzes pocket internal RMSD stability over extended simulation.
Compares to the 20ns ternary result (0.71 A).
"""
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import numpy as np
import json
import os

WORKDIR = os.path.join(DATA_ROOT, "topology")
os.chdir(WORKDIR)

TPR = 'md.tpr'
XTC = 'md.xtc'

POCKET_RESIDS = list(range(69, 77))  # resids 69-76 = 465-472 in PDB numbering

print('=== 500ns Interim Analysis ===')
print(f'Loading trajectory: {TPR} + {XTC}')

u = mda.Universe(TPR, XTC)
print(f'Topology: {len(u.atoms)} atoms, {len(u.residues)} residues')
print(f'Trajectory: {len(u.trajectory)} frames')

last_frame = u.trajectory[-1]
total_time_ns = last_frame.time / 1000.0
print(f'Total time: {total_time_ns:.1f} ns')
print(f'Frame interval: {u.trajectory.dt:.1f} ps')

pocket_sel = ' or '.join([f'resid {r}' for r in POCKET_RESIDS])
pocket = u.select_atoms(f'({pocket_sel}) and name CA')
pocket_all = u.select_atoms(pocket_sel)

print(f'Pocket CA atoms: {len(pocket)}')
print(f'Pocket all atoms: {len(pocket_all)}')
print(f'Pocket residues: {[r.resname + str(r.resid) for r in pocket.residues]}')

# Internal RMSD (aligned on pocket itself)
ref_positions = pocket.positions.copy()
pocket_rmsds = []
times = []

for ts in u.trajectory:
    pocket.translate(-pocket.center_of_mass())
    ref_centered = ref_positions - ref_positions.mean(axis=0)
    
    # Simple RMSD without alignment (use MDAnalysis built-in)
    rmsd_val = rms.rmsd(pocket.positions, ref_positions, superposition=True)
    pocket_rmsds.append(rmsd_val)
    times.append(ts.time / 1000.0)  # ps -> ns

pocket_rmsds = np.array(pocket_rmsds)
times = np.array(times)

# Statistics
mean_rmsd = np.mean(pocket_rmsds)
std_rmsd = np.std(pocket_rmsds)
max_rmsd = np.max(pocket_rmsds)
min_rmsd = np.min(pocket_rmsds)

print(f'\n=== Pocket Internal RMSD ===')
print(f'Mean:   {mean_rmsd:.3f} +/- {std_rmsd:.3f} A')
print(f'Range:  {min_rmsd:.3f} - {max_rmsd:.3f} A')
print(f'Median: {np.median(pocket_rmsds):.3f} A')

# Time-windowed analysis (50ns windows)
window_ns = 50.0
n_windows = int(total_time_ns / window_ns)
print(f'\n=== Time-windowed RMSD ({window_ns:.0f} ns windows) ===')
for w in range(n_windows + 1):
    t_start = w * window_ns
    t_end = min((w + 1) * window_ns, total_time_ns)
    mask = (times >= t_start) & (times < t_end)
    if mask.sum() == 0:
        continue
    w_mean = np.mean(pocket_rmsds[mask])
    w_std = np.std(pocket_rmsds[mask])
    print(f'  {t_start:.0f}-{t_end:.0f} ns: {w_mean:.3f} +/- {w_std:.3f} A ({mask.sum()} frames)')

# Radius of gyration over time
rg_values = []
for ts in u.trajectory:
    rg = pocket_all.radius_of_gyration()
    rg_values.append(rg)

rg_values = np.array(rg_values)
rg_change = (rg_values[-1] - rg_values[0]) / rg_values[0] * 100

print(f'\n=== Pocket Radius of Gyration ===')
print(f'Mean: {np.mean(rg_values):.3f} +/- {np.std(rg_values):.3f} A')
print(f'Start: {rg_values[0]:.3f} A, End: {rg_values[-1]:.3f} A')
print(f'Change: {rg_change:+.1f}%')

# CA-CA pairwise distances stability
print(f'\n=== Intra-pocket CA-CA Distance Stability ===')
n_ca = len(pocket)
pair_stds = []
pair_means = []

u.trajectory[0]
ref_dists = []
for i in range(n_ca):
    for j in range(i+1, n_ca):
        d = np.linalg.norm(pocket[i].position - pocket[j].position)
        ref_dists.append(d)

all_pair_dists = [[] for _ in range(len(ref_dists))]
for ts in u.trajectory:
    pair_idx = 0
    for i in range(n_ca):
        for j in range(i+1, n_ca):
            d = np.linalg.norm(pocket[i].position - pocket[j].position)
            all_pair_dists[pair_idx].append(d)
            pair_idx += 1

for dists in all_pair_dists:
    dists = np.array(dists)
    pair_stds.append(np.std(dists))
    pair_means.append(np.mean(dists))

pair_stds = np.array(pair_stds)
max_pair_dev = np.max([np.max(np.abs(np.array(d) - np.mean(d))) for d in all_pair_dists])

print(f'Number of CA-CA pairs: {len(pair_stds)}')
print(f'Mean pair std: {np.mean(pair_stds):.3f} A')
print(f'Max pair std: {np.max(pair_stds):.3f} A')
print(f'Max single-pair deviation: {max_pair_dev:.3f} A')

# Save results
results = {
    'trajectory_time_ns': float(total_time_ns),
    'n_frames': len(u.trajectory),
    'pocket_rmsd': {
        'mean': float(mean_rmsd),
        'std': float(std_rmsd),
        'max': float(max_rmsd),
        'min': float(min_rmsd),
        'median': float(np.median(pocket_rmsds))
    },
    'pocket_rg': {
        'mean': float(np.mean(rg_values)),
        'std': float(np.std(rg_values)),
        'start': float(rg_values[0]),
        'end': float(rg_values[-1]),
        'change_pct': float(rg_change)
    },
    'pairwise_distances': {
        'n_pairs': len(pair_stds),
        'mean_std': float(np.mean(pair_stds)),
        'max_std': float(np.max(pair_stds)),
        'max_deviation': float(max_pair_dev)
    },
    'comparison_20ns': {
        'pocket_rmsd_20ns': 0.708,
        'note': 'Compare to 20ns ternary (0.708 A) and 100ns fragment (3.040 A)'
    }
}

with open(os.path.join(DATA_ROOT, "500ns_interim_analysis.json"), 'w') as f:
    json.dump(results, f, indent=2)

print(f'\n=== Results saved to 500ns_interim_analysis.json ===')
print(f'\n=== HEADLINE: Pocket RMSD over {total_time_ns:.0f} ns = {mean_rmsd:.3f} +/- {std_rmsd:.3f} A ===')
