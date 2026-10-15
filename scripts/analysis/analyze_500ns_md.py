#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Parse GROMACS XVG output files from 500ns MD analysis."""
import json
import statistics

def parse_xvg(path, skip_header=True):
    """Parse an XVG file, returning list of (x, y1, y2, ...) tuples."""
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                data.append([float(p) for p in parts])
    return data

results = {}

# 1. Pocket RMSD
print("=== Pocket RMSD (CA of relay 465-472) ===")
rmsd_data = parse_xvg(r"D:\md-500ns-analysis\pocket_rmsd.xvg")
rmsd_values = [row[1] for row in rmsd_data]  # nm
rmsd_mean = statistics.mean(rmsd_values)
rmsd_std = statistics.stdev(rmsd_values)
rmsd_max = max(rmsd_values)
rmsd_min = min(rmsd_values)

# Convert nm to Angstrom
rmsd_mean_A = rmsd_mean * 10
rmsd_std_A = rmsd_std * 10
rmsd_max_A = rmsd_max * 10
rmsd_min_A = rmsd_min * 10

print(f"  Mean: {rmsd_mean_A:.2f} +/- {rmsd_std_A:.2f} A")
print(f"  Max:  {rmsd_max_A:.2f} A")
print(f"  Min:  {rmsd_min_A:.2f} A")
print(f"  Frames: {len(rmsd_values)}")

# Time segments for drift analysis
segment_size = len(rmsd_values) // 5
for i in range(5):
    seg = rmsd_values[i*segment_size:(i+1)*segment_size]
    seg_mean = statistics.mean(seg) * 10
    seg_std = statistics.stdev(seg) * 10
    t_start = rmsd_data[i*segment_size][0]
    t_end = rmsd_data[min((i+1)*segment_size-1, len(rmsd_data)-1)][0]
    print(f"  Segment {i+1} ({t_start:.0f}-{t_end:.0f} ns): {seg_mean:.2f} +/- {seg_std:.2f} A")

results['pocket_rmsd'] = {
    'mean_angstrom': round(rmsd_mean_A, 3),
    'std_angstrom': round(rmsd_std_A, 3),
    'max_angstrom': round(rmsd_max_A, 3),
    'min_angstrom': round(rmsd_min_A, 3),
    'frames': len(rmsd_values),
    'unit': 'angstrom',
    'description': 'CA RMSD of relay pocket residues 465-472 vs initial structure, fit to backbone'
}

# 2. RMSF
print("\n=== Per-Residue RMSF (C-alpha) ===")
rmsf_data = parse_xvg(r"D:\md-500ns-analysis\rmsf_ca.xvg")
rmsf_all = {int(row[0]): row[1]*10 for row in rmsf_data}  # nm to A

relay_resids = [465, 466, 467, 468, 469, 470, 471, 472]
relay_names = ['ARG465', 'PRO466', 'GLY467', 'GLY468', 'GLY469', 'ASP470', 'MET471', 'ARG472']

print("  Relay pocket RMSF:")
relay_rmsf = {}
for resid, name in zip(relay_resids, relay_names):
    if resid in rmsf_all:
        v = rmsf_all[resid]
        status = "RIGID" if v < 1.5 else "FLEXIBLE" if v > 3.0 else "MODERATE"
        print(f"    {name}: {v:.3f} A ({status})")
        relay_rmsf[name] = round(v, 3)

all_rmsf_values = list(rmsf_all.values())
relay_rmsf_values = [rmsf_all.get(r, 0) for r in relay_resids if r in rmsf_all]
print(f"\n  System mean RMSF: {statistics.mean(all_rmsf_values):.3f} A")
print(f"  Relay pocket mean RMSF: {statistics.mean(relay_rmsf_values):.3f} A")
print(f"  System median RMSF: {statistics.median(all_rmsf_values):.3f} A")

# Percentile of relay pocket
sorted_all = sorted(all_rmsf_values)
relay_mean_rmsf = statistics.mean(relay_rmsf_values)
below = sum(1 for v in sorted_all if v <= relay_mean_rmsf)
pct = below / len(sorted_all) * 100
print(f"  Relay pocket RMSF percentile: {pct:.1f}% (lower = more rigid)")

results['rmsf'] = {
    'relay_pocket': relay_rmsf,
    'relay_mean_angstrom': round(relay_mean_rmsf, 3),
    'system_mean_angstrom': round(statistics.mean(all_rmsf_values), 3),
    'relay_percentile': round(pct, 1),
    'total_residues': len(all_rmsf_values)
}

# 3. Radius of Gyration
print("\n=== Pocket Radius of Gyration ===")
rg_data = parse_xvg(r"D:\md-500ns-analysis\pocket_rg.xvg")
rg_values = [row[1]*10 for row in rg_data]  # nm to A
rg_mean = statistics.mean(rg_values)
rg_std = statistics.stdev(rg_values)

# Initial and final Rg
rg_initial = rg_values[0]
rg_final = statistics.mean(rg_values[-1000:])
rg_change = (rg_final - rg_initial) / rg_initial * 100

print(f"  Mean: {rg_mean:.2f} +/- {rg_std:.2f} A")
print(f"  Initial: {rg_initial:.2f} A")
print(f"  Final (last 10 ns avg): {rg_final:.2f} A")
print(f"  Change: {rg_change:+.1f}%")

results['pocket_rg'] = {
    'mean_angstrom': round(rg_mean, 2),
    'std_angstrom': round(rg_std, 2),
    'initial_angstrom': round(rg_initial, 2),
    'final_angstrom': round(rg_final, 2),
    'change_percent': round(rg_change, 1)
}

# 4. SASA
print("\n=== Pocket SASA ===")
sasa_data = parse_xvg(r"D:\md-500ns-analysis\pocket_sasa.xvg")
sasa_values = [row[1] for row in sasa_data]  # nm^2
sasa_A2 = [v * 100 for v in sasa_values]  # convert to A^2
sasa_mean = statistics.mean(sasa_A2)
sasa_std = statistics.stdev(sasa_A2)

print(f"  Mean: {sasa_mean:.1f} +/- {sasa_std:.1f} A^2")
print(f"  Frames: {len(sasa_values)}")

# Segments
for i in range(5):
    seg = sasa_A2[i*(len(sasa_A2)//5):(i+1)*(len(sasa_A2)//5)]
    print(f"  Segment {i+1}: {statistics.mean(seg):.1f} +/- {statistics.stdev(seg):.1f} A^2")

results['pocket_sasa'] = {
    'mean_angstrom2': round(sasa_mean, 1),
    'std_angstrom2': round(sasa_std, 1),
    'frames': len(sasa_values)
}

# Summary
print("\n" + "="*60)
print("500ns MD ANALYSIS SUMMARY")
print("="*60)
print(f"System: Reduced relay pocket (105 residues, 33,713 atoms)")
print(f"Duration: 500 ns (50,001 frames)")
print(f"Force field: AMBER99SB-ILDN + TIP3P")
print(f"")
print(f"Pocket RMSD (CA 465-472): {rmsd_mean_A:.2f} +/- {rmsd_std_A:.2f} A (max {rmsd_max_A:.2f})")
print(f"Pocket RMSF (mean):       {relay_mean_rmsf:.3f} A ({pct:.0f}th percentile)")
print(f"Pocket Rg:                {rg_mean:.2f} +/- {rg_std:.2f} A (change: {rg_change:+.1f}%)")
print(f"Pocket SASA:              {sasa_mean:.1f} +/- {sasa_std:.1f} A^2")

# Compare to 100ns
print(f"\nComparison to 100ns reduced system:")
print(f"  100ns pocket RMSD: 2.47 +/- 0.23 A")
print(f"  500ns pocket RMSD: {rmsd_mean_A:.2f} +/- {rmsd_std_A:.2f} A")
print(f"  100ns pocket Rg:   9.70 +/- 0.23 A")
print(f"  500ns pocket Rg:   {rg_mean:.2f} +/- {rg_std:.2f} A")
print(f"  100ns SASA:        1252.6 +/- 43.1 A^2")
print(f"  500ns SASA:        {sasa_mean:.1f} +/- {sasa_std:.1f} A^2")

results['metadata'] = {
    'system': 'Reduced relay pocket (105 residues, 33713 atoms, AMBER99SB-ILDN/TIP3P)',
    'duration_ns': 500,
    'frames': 50001,
    'timestep_ps': 10,
    'cores': 16,
    'completed': 'Aug 6, 2026 22:38 UTC',
    'runtime': '16d 7h 47m',
    'performance_ns_per_day': 24.503,
    'comparison_100ns': {
        'pocket_rmsd_angstrom': '2.47 +/- 0.23',
        'pocket_rg_angstrom': '9.70 +/- 0.23',
        'pocket_sasa_angstrom2': '1252.6 +/- 43.1'
    }
}

with open("./md_500ns_analysis.json", 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to md_500ns_analysis.json")
