import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Parse MMPBSA decomposition results to find top contributing residues."""
import csv
import json
import os

WORKDIR = os.path.join(DATA_ROOT, "results", "08-binding-energetics")
os.chdir(WORKDIR)

# Residue name mapping (AMBER numbering -> original PDB)
# From check_chains.py: residues 0-81 = gp120 chain G, 82-104 = CD4 chain A
ORIG_RESIDS = [
    94,95,100,101,104,253,254,255,256,257,258,259,260,261,
    272,273,274,275,276,277,279,280,281,282,283,284,285,286,287,
    359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,
    381,382,383,385,388,389,420,421,422,423,424,
    444,445,446,447,448,449,450,451,452,
    463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,
    23,24,25,26,27,28,29,35,36,37,38,39,40,41,42,43,44,45,46,47,48,59,62
]

# Parse the decomposition CSV
decomp_file = 'FINAL_DECOMP_MMPBSA.dat'
with open(decomp_file) as f:
    content = f.read()

# Find the "Differences" section with "Total Energy Decomposition"
sections = content.split('DELTAS:')
if len(sections) < 2:
    print("ERROR: Could not find DELTAS section")
    exit(1)

diff_section = sections[1]
parts = diff_section.split('Total Energy Decomposition:')
if len(parts) < 2:
    print("ERROR: Could not find Total Energy Decomposition in DELTAS")
    exit(1)

decomp_data = parts[1].strip()
lines = decomp_data.split('\n')

# Skip header lines
data_lines = []
for line in lines:
    line = line.strip()
    if not line or line.startswith('Residue') or line.startswith(','):
        continue
    if line.startswith('Sidechain') or line.startswith('Backbone'):
        break
    data_lines.append(line)

# Parse each residue
results = []
for i, line in enumerate(data_lines):
    parts = line.split(',')
    if len(parts) < 19:
        continue
    
    res_label = parts[0].strip()
    # Parse: "MET   1" or "R MET   1" etc.
    # The location field follows
    location = parts[1].strip() if len(parts) > 1 else ''
    
    try:
        total_avg = float(parts[-3])
        total_std = float(parts[-2])
        vdw_avg = float(parts[3])
        eel_avg = float(parts[6])
        polar_avg = float(parts[9])
        nonpolar_avg = float(parts[12])
    except (ValueError, IndexError):
        continue
    
    # Map to original residue ID
    orig_resid = ORIG_RESIDS[i] if i < len(ORIG_RESIDS) else -1
    chain = 'G (gp120)' if i < 82 else 'A (CD4)'
    
    results.append({
        'amber_label': res_label,
        'location': location,
        'orig_resid': orig_resid,
        'chain': chain,
        'delta_total': total_avg,
        'delta_total_std': total_std,
        'delta_vdw': vdw_avg,
        'delta_eel': eel_avg,
        'delta_polar_solv': polar_avg,
        'delta_nonpolar_solv': nonpolar_avg,
        'is_relay_pocket': 465 <= orig_resid <= 472 if i < 82 else False
    })

# Sort by absolute contribution
results.sort(key=lambda x: abs(x['delta_total']), reverse=True)

print("=" * 80)
print("MM-GBSA Per-Residue Binding Energy Decomposition")
print("gp120-CD4 Interface (102 frames from 100ns trajectory)")
print("=" * 80)

print(f"\n{'Rank':>4} {'Res':>8} {'Orig#':>6} {'Chain':>12} {'DG_total':>10} {'VdW':>10} {'Elec':>10} {'Polar':>10} {'Relay?':>7}")
print("-" * 80)

for rank, r in enumerate(results[:30], 1):
    relay = " ***" if r['is_relay_pocket'] else ""
    print(f"{rank:4d} {r['amber_label']:>8} {r['orig_resid']:6d} {r['chain']:>12} {r['delta_total']:10.3f} {r['delta_vdw']:10.3f} {r['delta_eel']:10.3f} {r['delta_polar_solv']:10.3f}{relay}")

# Relay pocket residues specifically
print(f"\n{'=' * 80}")
print("RELAY POCKET RESIDUES (465-472)")
print(f"{'=' * 80}")
relay_total = 0
for r in results:
    if r['is_relay_pocket']:
        print(f"  {r['amber_label']:>8} (orig {r['orig_resid']}) DG = {r['delta_total']:8.3f} kcal/mol  (VdW={r['delta_vdw']:.3f}, Elec={r['delta_eel']:.3f})")
        relay_total += r['delta_total']
print(f"\n  Total relay pocket contribution: {relay_total:.3f} kcal/mol")

# Overall binding
print(f"\n{'=' * 80}")
print("OVERALL BINDING FREE ENERGY")
print(f"{'=' * 80}")
total_binding = sum(r['delta_total'] for r in results)
print(f"  Sum of all per-residue contributions: {total_binding:.3f} kcal/mol")

# Save as JSON
json_results = {
    'method': 'MM-GBSA (igb=5, saltcon=0.15M)',
    'frames': 102,
    'trajectory': '100ns reduced system (gp120 + CD4)',
    'interval': 'every 10th frame',
    'overall_binding': {
        'delta_total': -31.372,
        'delta_vdw': -52.545,
        'delta_eel': -364.479,
        'delta_gb': 393.678,
        'delta_esurf': -8.063
    },
    'per_residue': results,
    'relay_pocket_total': relay_total
}

with open('mmpbsa_results.json', 'w') as f:
    json.dump(json_results, f, indent=2)
print(f"\nResults saved to {WORKDIR}/mmpbsa_results.json")
