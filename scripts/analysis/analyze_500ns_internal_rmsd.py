#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Parse internal pocket RMSD from 500ns simulation."""
import statistics, json

def parse_xvg(path):
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                data.append((float(parts[0]), float(parts[1])))
    return data

# Internal RMSD (CA only, fit to pocket CA)
ca_data = parse_xvg("./pocket_internal_rmsd.xvg")
ca_rmsd = [r[1]*10 for r in ca_data]  # nm to A

# Internal RMSD (all atoms, fit to pocket all)
all_data = parse_xvg("./pocket_all_internal_rmsd.xvg")
all_rmsd = [r[1]*10 for r in all_data]  # nm to A

print("=== INTERNAL Pocket RMSD (fit to pocket, measure pocket) ===")
print(f"\nCA-only (8 atoms):")
print(f"  Mean: {statistics.mean(ca_rmsd):.2f} +/- {statistics.stdev(ca_rmsd):.2f} A")
print(f"  Max:  {max(ca_rmsd):.2f} A")
print(f"  Min:  {min(ca_rmsd):.4f} A")
print(f"  Frames: {len(ca_rmsd)}")

print(f"\nAll heavy atoms (112 atoms):")
print(f"  Mean: {statistics.mean(all_rmsd):.2f} +/- {statistics.stdev(all_rmsd):.2f} A")
print(f"  Max:  {max(all_rmsd):.2f} A")

# Time segments
print("\n--- Time Evolution (CA internal RMSD) ---")
segment_size = len(ca_rmsd) // 5
for i in range(5):
    seg = ca_rmsd[i*segment_size:(i+1)*segment_size]
    t_start = ca_data[i*segment_size][0]
    t_end = ca_data[min((i+1)*segment_size-1, len(ca_data)-1)][0]
    print(f"  {t_start:.0f}-{t_end:.0f} ns: {statistics.mean(seg):.2f} +/- {statistics.stdev(seg):.2f} A (max {max(seg):.2f})")

# Compare to reference values
print("\n--- Comparison ---")
print(f"  100ns reduced system internal RMSD: 2.47 +/- 0.23 A")
print(f"  500ns reduced system internal RMSD: {statistics.mean(ca_rmsd):.2f} +/- {statistics.stdev(ca_rmsd):.2f} A")
print(f"  20ns full ternary internal RMSD:    0.71 +/- 0.15 A")
print(f"  100ns isolated fragment RMSD:       3.04 A")

# Stability assessment
mean_val = statistics.mean(ca_rmsd)
if mean_val < 1.5:
    verdict = "OUTSTANDING - sub-angstrom stability maintained at 500ns"
elif mean_val < 3.0:
    verdict = "GOOD - stable pocket with moderate fluctuations"
elif mean_val < 5.0:
    verdict = "MODERATE - conformational flexibility observed"
else:
    verdict = "SIGNIFICANT MOTION - fragment lacks intermolecular stabilization"

print(f"\n  Verdict: {verdict}")

# Save results
results = {
    'pocket_internal_rmsd_ca': {
        'mean_angstrom': round(statistics.mean(ca_rmsd), 3),
        'std_angstrom': round(statistics.stdev(ca_rmsd), 3),
        'max_angstrom': round(max(ca_rmsd), 3),
        'frames': len(ca_rmsd)
    },
    'pocket_internal_rmsd_all': {
        'mean_angstrom': round(statistics.mean(all_rmsd), 3),
        'std_angstrom': round(statistics.stdev(all_rmsd), 3),
        'max_angstrom': round(max(all_rmsd), 3),
        'frames': len(all_rmsd)
    }
}

# Update the main results file
try:
    with open("./md_500ns_analysis.json") as f:
        main_results = json.load(f)
    main_results['pocket_internal_rmsd_ca'] = results['pocket_internal_rmsd_ca']
    main_results['pocket_internal_rmsd_all'] = results['pocket_internal_rmsd_all']
    with open("./md_500ns_analysis.json", 'w') as f:
        json.dump(main_results, f, indent=2)
    print("\nUpdated md_500ns_analysis.json with internal RMSD data")
except Exception as e:
    print(f"Error updating JSON: {e}")
