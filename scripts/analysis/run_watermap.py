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
WaterMap / Hydration Site Analysis
Analyzes water structure around the relay pocket from the MD trajectory.
Identifies ordered waters, residence times, and displaceable hydration sites.
"""

import numpy as np
import MDAnalysis as mda
import json, os, time

print("=== WaterMap / Hydration Site Analysis ===")
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")

TPR = os.path.join(DATA_ROOT, "topology", "md.tpr")
XTC = os.path.join(DATA_ROOT, "topology", "md.xtc")  # Note: large trajectory not in package
OUT_DIR = os.path.join(DATA_ROOT, "results", "01-pocket-discovery")
os.makedirs(OUT_DIR, exist_ok=True)

POCKET_RESIDUES = list(range(69, 77))
SHELL_RADIUS = 6.0
GRID_SPACING = 1.0
SKIP_FRAMES = 5

print(f"\nLoading universe: {TPR} + {XTC}")
u = mda.Universe(TPR, XTC)
print(f"  Atoms: {len(u.atoms)}")
print(f"  Frames: {len(u.trajectory)}")
print(f"  Analyzing every {SKIP_FRAMES}th frame = {len(u.trajectory)//SKIP_FRAMES} frames")

pocket_sel = " or ".join([f"resid {r}" for r in POCKET_RESIDUES])
pocket = u.select_atoms(f"protein and ({pocket_sel})")
waters = u.select_atoms("name OW")

print(f"  Pocket atoms: {len(pocket)}")
print(f"  Water oxygens: {len(waters)}")

if len(pocket) == 0:
    print("ERROR: No pocket atoms found! Check residue numbering.")
    exit(1)

print("\n--- Step 1: Collecting water positions near pocket ---")
pocket_center_all = []
water_positions = []
frame_count = 0

for ts in u.trajectory[::SKIP_FRAMES]:
    frame_count += 1
    if frame_count % 100 == 0:
        print(f"  Frame {frame_count}...")

    pc = pocket.center_of_mass()
    pocket_center_all.append(pc.copy())

    water_coords = waters.positions
    dists = np.linalg.norm(water_coords - pc, axis=1)
    nearby = water_coords[dists < SHELL_RADIUS]

    for pos in nearby:
        water_positions.append(pos.copy())

water_positions = np.array(water_positions)
pocket_center_mean = np.mean(pocket_center_all, axis=0)
print(f"  Collected {len(water_positions)} water position samples")
print(f"  Mean pocket center: ({pocket_center_mean[0]:.1f}, {pocket_center_mean[1]:.1f}, {pocket_center_mean[2]:.1f})")

print("\n--- Step 2: Grid-based density analysis ---")
grid_extent = SHELL_RADIUS + 1.0
x_range = np.arange(pocket_center_mean[0] - grid_extent, pocket_center_mean[0] + grid_extent, GRID_SPACING)
y_range = np.arange(pocket_center_mean[1] - grid_extent, pocket_center_mean[1] + grid_extent, GRID_SPACING)
z_range = np.arange(pocket_center_mean[2] - grid_extent, pocket_center_mean[2] + grid_extent, GRID_SPACING)

print(f"  Grid dimensions: {len(x_range)} x {len(y_range)} x {len(z_range)} = {len(x_range)*len(y_range)*len(z_range)} voxels")

density = np.zeros((len(x_range), len(y_range), len(z_range)))

for pos in water_positions:
    ix = int((pos[0] - x_range[0]) / GRID_SPACING)
    iy = int((pos[1] - y_range[0]) / GRID_SPACING)
    iz = int((pos[2] - z_range[0]) / GRID_SPACING)
    if 0 <= ix < len(x_range) and 0 <= iy < len(y_range) and 0 <= iz < len(z_range):
        density[ix, iy, iz] += 1

density_per_frame = density / frame_count

print("\n--- Step 3: Identifying hydration sites ---")

nonzero = density_per_frame[density_per_frame > 0]
if len(nonzero) == 0:
    print("ERROR: No water density found at all")
    exit(1)
bulk_mean = nonzero.mean()
threshold = bulk_mean * 1.5
print(f"  Bulk mean density: {bulk_mean:.4f} per frame")
print(f"  Threshold (2x bulk): {threshold:.4f}")

sites = []
visited = set()

for ix in range(1, len(x_range)-1):
    for iy in range(1, len(y_range)-1):
        for iz in range(1, len(z_range)-1):
            val = density_per_frame[ix, iy, iz]
            if val < threshold:
                continue
            if (ix, iy, iz) in visited:
                continue

            is_max = True
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    for dz in [-1, 0, 1]:
                        if dx == 0 and dy == 0 and dz == 0:
                            continue
                        if density_per_frame[ix+dx, iy+dy, iz+dz] > val:
                            is_max = False
                            break
                    if not is_max:
                        break
                if not is_max:
                    break

            if is_max:
                x = float(x_range[ix])
                y = float(y_range[iy])
                z = float(z_range[iz])
                dist_to_center = float(np.linalg.norm([x - pocket_center_mean[0],
                                                       y - pocket_center_mean[1],
                                                       z - pocket_center_mean[2]]))
                sites.append({
                    "id": len(sites) + 1,
                    "x": round(x, 2),
                    "y": round(y, 2),
                    "z": round(z, 2),
                    "density": round(float(val), 4),
                    "density_ratio": round(float(val / bulk_mean), 2),
                    "dist_to_pocket_center": round(dist_to_center, 2),
                })
                for dx in [-1, 0, 1]:
                    for dy in [-1, 0, 1]:
                        for dz in [-1, 0, 1]:
                            visited.add((ix+dx, iy+dy, iz+dz))

sites.sort(key=lambda s: -s["density"])
print(f"  Found {len(sites)} hydration sites above threshold")

print("\n--- Step 4: Classifying hydration sites ---")
structural = []
displaceable = []
peripheral = []

for s in sites:
    if s["density_ratio"] > 2.5 and s["dist_to_pocket_center"] < 5.5:
        s["class"] = "structural"
        structural.append(s)
    elif s["density_ratio"] > 1.5 and s["dist_to_pocket_center"] < 6.5:
        s["class"] = "displaceable"
        displaceable.append(s)
    else:
        s["class"] = "peripheral"
        peripheral.append(s)

print(f"  Structural (keep):      {len(structural)}")
print(f"  Displaceable (targets): {len(displaceable)}")
print(f"  Peripheral:             {len(peripheral)}")

print("\n--- Results ---")
print(f"\n{'ID':>4} {'Position':>24} {'Density':>10} {'Ratio':>8} {'Dist':>8} {'Class':>14}")
print("-" * 72)
for s in sites[:20]:
    pos = f"({s['x']:.1f}, {s['y']:.1f}, {s['z']:.1f})"
    print(f"{s['id']:>4} {pos:>24} {s['density']:>10.4f} {s['density_ratio']:>7.1f}x {s['dist_to_pocket_center']:>7.2f} {s['class']:>14}")

output = {
    "analysis": "WaterMap / Hydration Site Analysis",
    "trajectory": XTC,
    "frames_analyzed": frame_count,
    "skip_frames": SKIP_FRAMES,
    "shell_radius_A": SHELL_RADIUS,
    "grid_spacing_A": GRID_SPACING,
    "pocket_residues_in_tpr": POCKET_RESIDUES,
    "pocket_residues_original": [465, 466, 467, 468, 469, 470, 471, 472],
    "pocket_center": [round(float(x), 2) for x in pocket_center_mean],
    "bulk_mean_density": round(float(bulk_mean), 6),
    "threshold": round(float(threshold), 6),
    "total_sites": len(sites),
    "structural_sites": len(structural),
    "displaceable_sites": len(displaceable),
    "peripheral_sites": len(peripheral),
    "sites": sites[:30],
}

outfile = os.path.join(OUT_DIR, "watermap_results.json")
with open(outfile, 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved to: {outfile}")

print(f"\n=== WaterMap COMPLETE ===")
print(f"Finished: {time.strftime('%Y-%m-%d %H:%M:%S')}")
