import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Quick scan of idecomp=3 DELTAS to find significant pairwise interactions."""
import sys

DECOMP = os.path.join(DATA_ROOT, "mmpbsa_idecomp3/FINAL_DECOMP_id3.dat")
DELTAS_LINE = 1556613  # 1-indexed line where "DELTAS:" appears

significant = []
in_deltas = False
count = 0

with open(DECOMP) as f:
    for i, line in enumerate(f, 1):
        if i == DELTAS_LINE:
            in_deltas = True
            continue
        if not in_deltas:
            continue
        if not line.startswith("R:"):
            continue

        parts = line.strip().split(",")
        if len(parts) < 18:
            continue

        res1 = parts[0]  # e.g. R:A:ASP:1
        res2 = parts[1]  # e.g. R:B:PHE:43

        try:
            total_avg = float(parts[15])
            vdw_avg = float(parts[3])
            elec_avg = float(parts[6])
            polar_avg = float(parts[9])
        except (ValueError, IndexError):
            continue

        if abs(total_avg) > 0.01:
            significant.append((total_avg, vdw_avg, elec_avg, polar_avg, res1, res2))

        count += 1
        if count % 100000 == 0:
            print(f"  ...processed {count} pairs, {len(significant)} significant so far", file=sys.stderr)

print(f"\nTotal pairs in DELTAS: {count}")
print(f"Significant pairs (|total| > 0.01): {len(significant)}")

significant.sort(key=lambda x: x[0])

print(f"\n--- Top 30 STRONGEST (most negative total) ---")
print(f"{'Rank':>4} {'Res1':>15} {'Res2':>15} {'VdW':>10} {'Elec':>10} {'Polar':>10} {'Total':>10}")
print("-" * 80)
for i, (total, vdw, elec, polar, r1, r2) in enumerate(significant[:30]):
    print(f"{i+1:>4} {r1:>15} {r2:>15} {vdw:>10.4f} {elec:>10.4f} {polar:>10.4f} {total:>10.4f}")

print(f"\n--- Top 30 STRONGEST (most positive total) ---")
for i, (total, vdw, elec, polar, r1, r2) in enumerate(reversed(significant[-30:])):
    print(f"{i+1:>4} {r1:>15} {r2:>15} {vdw:>10.4f} {elec:>10.4f} {polar:>10.4f} {total:>10.4f}")

# Cross-chain analysis
relay_nums = set(range(372, 380))
interface_pairs = []
relay_pairs = []

for total, vdw, elec, polar, r1, r2 in significant:
    r1_parts = r1.split(":")
    r2_parts = r2.split(":")
    r1_chain = r1_parts[1]
    r2_chain = r2_parts[1]
    r1_num = int(r1_parts[3])
    r2_num = int(r2_parts[3])
    r1_name = r1_parts[2]
    r2_name = r2_parts[2]

    if r1_chain != r2_chain:
        interface_pairs.append((total, vdw, elec, polar, r1, r2))
        if r1_num in relay_nums or r2_num in relay_nums:
            relay_pairs.append((total, vdw, elec, polar, r1, r2))

print(f"\n--- Cross-chain interface pairs: {len(interface_pairs)} ---")
interface_pairs.sort(key=lambda x: x[0])
print(f"{'Rank':>4} {'Res1':>15} {'Res2':>15} {'VdW':>10} {'Elec':>10} {'Polar':>10} {'Total':>10}")
print("-" * 80)
for i, (total, vdw, elec, polar, r1, r2) in enumerate(interface_pairs[:30]):
    print(f"{i+1:>4} {r1:>15} {r2:>15} {vdw:>10.4f} {elec:>10.4f} {polar:>10.4f} {total:>10.4f}")

print(f"\n--- Relay pocket interface pairs: {len(relay_pairs)} ---")
relay_pairs.sort(key=lambda x: x[0])
for i, (total, vdw, elec, polar, r1, r2) in enumerate(relay_pairs):
    print(f"{i+1:>4} {r1:>15} {r2:>15} {vdw:>10.4f} {elec:>10.4f} {polar:>10.4f} {total:>10.4f}")

relay_total = sum(x[0] for x in relay_pairs)
interface_total = sum(x[0] for x in interface_pairs)
print(f"\nRelay pocket interface sum: {relay_total:.4f} kcal/mol")
print(f"Total interface sum: {interface_total:.4f} kcal/mol")

import json
results = {
    "method": "idecomp=3 pairwise DELTAS scan",
    "total_pairs_in_deltas": count,
    "significant_pairs": len(significant),
    "interface_pairs_count": len(interface_pairs),
    "relay_interface_pairs_count": len(relay_pairs),
    "interface_total_dG": round(interface_total, 4),
    "relay_total_dG": round(relay_total, 4),
    "top20_interface": [],
    "relay_pairs": [],
}
for i, (total, vdw, elec, polar, r1, r2) in enumerate(interface_pairs[:20]):
    results["top20_interface"].append({
        "rank": i + 1, "res1": r1, "res2": r2,
        "dG_total": round(total, 4), "dG_vdw": round(vdw, 4),
        "dG_elec": round(elec, 4), "dG_polar": round(polar, 4),
    })
for total, vdw, elec, polar, r1, r2 in relay_pairs:
    results["relay_pairs"].append({
        "res1": r1, "res2": r2,
        "dG_total": round(total, 4), "dG_vdw": round(vdw, 4),
        "dG_elec": round(elec, 4), "dG_polar": round(polar, 4),
    })
with open(os.path.join(DATA_ROOT, "mmpbsa_idecomp3/idecomp3_scan_results.json"), "w") as f:
    json.dump(results, f, indent=2)
print("\nResults saved to idecomp3_scan_results.json")
