import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Extract raw Complex pairwise energies for interface residues from idecomp=3.

Since DELTAS are unreliable on this large system, we extract the raw
Complex-state pairwise interactions between receptor (gp120, R:A:*)
and ligand (CD4, L:B:*) residues. These represent the actual physical
interactions in the bound state.
"""
import json, sys

DECOMP = os.path.join(DATA_ROOT, "mmpbsa_idecomp3/FINAL_DECOMP_id3.dat")
RELAY = set(range(372, 380))

interface = []
count = 0
in_complex = False
past_header = False

with open(DECOMP) as f:
    for line in f:
        if "Complex:" in line:
            in_complex = True
            past_header = False
            continue
        if in_complex and "Resid 1,Resid 2" in line:
            past_header = True
            continue
        if in_complex and past_header and line.startswith(",,"):
            continue
        if "Receptor:" in line or "Ligand:" in line or "DELTAS:" in line:
            break

        if not in_complex or not past_header:
            continue
        if not line.startswith("R:"):
            continue

        parts = line.strip().split(",")
        if len(parts) < 18:
            continue

        r1 = parts[0]
        r2 = parts[1]
        r1_parts = r1.split(":")
        r2_parts = r2.split(":")
        r1_chain = r1_parts[1]
        r2_chain = r2_parts[1]

        if r1_chain == r2_chain:
            count += 1
            continue

        try:
            vdw = float(parts[3])
            elec = float(parts[6])
            polar = float(parts[9])
            nonpol = float(parts[12])
            total = float(parts[15])
        except (ValueError, IndexError):
            count += 1
            continue

        r1_num = int(r1_parts[3])
        r2_num = int(r2_parts[3])
        r1_name = r1_parts[2]
        r2_name = r2_parts[2]

        is_relay = r1_num in RELAY or r2_num in RELAY

        interface.append({
            "r1": r1, "r2": r2,
            "r1_num": r1_num, "r2_num": r2_num,
            "r1_name": r1_name, "r2_name": r2_name,
            "vdw": vdw, "elec": elec, "polar": polar,
            "nonpol": nonpol, "total": total,
            "is_relay": is_relay,
        })
        count += 1

        if count % 200000 == 0:
            print(f"  ...{count} pairs, {len(interface)} interface", file=sys.stderr)

print(f"Total pairs scanned in Complex section: {count}")
print(f"Cross-chain interface pairs: {len(interface)}")

interface.sort(key=lambda x: x["total"])

print(f"\n{'='*85}")
print(f"Top 30 STRONGEST interface pairs (most NEGATIVE total = most favorable)")
print(f"{'='*85}")
print(f"{'Rank':>4} {'gp120 Res':>14} {'CD4 Res':>14} {'VdW':>8} {'Elec':>10} {'Polar':>10} {'Total':>10} {'Relay':>6}")
print("-" * 85)
for i, p in enumerate(interface[:30]):
    if p["r1"].startswith("R:A"):
        gp120 = f"{p['r1_name']}{p['r1_num']}"
        cd4 = f"{p['r2_name']}{p['r2_num']}"
        gnum = p['r1_num']
    else:
        gp120 = f"{p['r2_name']}{p['r2_num']}"
        cd4 = f"{p['r1_name']}{p['r1_num']}"
        gnum = p['r2_num']
    relay = "YES" if p["is_relay"] else ""
    print(f"{i+1:>4} {gp120:>14} {cd4:>14} {p['vdw']:>8.3f} {p['elec']:>10.3f} {p['polar']:>10.3f} {p['total']:>10.3f} {relay:>6}")

relay_iface = [p for p in interface if p["is_relay"]]
relay_iface.sort(key=lambda x: x["total"])

print(f"\n{'='*85}")
print(f"Relay pocket interface pairs ({len(relay_iface)} total)")
print(f"{'='*85}")
print(f"{'gp120 Res':>14} {'CD4 Res':>14} {'VdW':>8} {'Elec':>10} {'Polar':>10} {'Total':>10}")
print("-" * 70)
for p in relay_iface:
    if p["r1"].startswith("R:A"):
        gp120 = f"{p['r1_name']}{p['r1_num']}"
        cd4 = f"{p['r2_name']}{p['r2_num']}"
    else:
        gp120 = f"{p['r2_name']}{p['r2_num']}"
        cd4 = f"{p['r1_name']}{p['r1_num']}"
    print(f"{gp120:>14} {cd4:>14} {p['vdw']:>8.3f} {p['elec']:>10.3f} {p['polar']:>10.3f} {p['total']:>10.3f}")

relay_vdw = sum(p["vdw"] for p in relay_iface)
relay_elec = sum(p["elec"] for p in relay_iface)
relay_polar = sum(p["polar"] for p in relay_iface)
relay_total = sum(p["total"] for p in relay_iface)
iface_total = sum(p["total"] for p in interface)
iface_vdw = sum(p["vdw"] for p in interface)
iface_elec = sum(p["elec"] for p in interface)

print(f"\n--- Summary ---")
print(f"Total interface pairs: {len(interface)}")
print(f"Interface VdW sum: {iface_vdw:.3f}")
print(f"Interface Elec sum: {iface_elec:.3f}")
print(f"Interface Total sum: {iface_total:.3f}")
print(f"\nRelay pocket pairs: {len(relay_iface)}")
print(f"Relay VdW sum: {relay_vdw:.3f}")
print(f"Relay Elec sum: {relay_elec:.3f}")
print(f"Relay Polar sum: {relay_polar:.3f}")
print(f"Relay Total sum: {relay_total:.3f}")
if iface_total != 0:
    print(f"Relay fraction of interface: {relay_total/iface_total*100:.1f}%")

results = {
    "method": "idecomp=3 Complex pairwise (raw, not DELTAS)",
    "note": "DELTAS unreliable on large system; Complex pairwise shows direct physical interactions",
    "total_interface_pairs": len(interface),
    "relay_interface_pairs": len(relay_iface),
    "interface_vdw_sum": round(iface_vdw, 3),
    "interface_elec_sum": round(iface_elec, 3),
    "interface_total_sum": round(iface_total, 3),
    "relay_vdw_sum": round(relay_vdw, 3),
    "relay_elec_sum": round(relay_elec, 3),
    "relay_total_sum": round(relay_total, 3),
    "top30_interface": [],
    "relay_pairs": [],
}

for i, p in enumerate(interface[:30]):
    results["top30_interface"].append({
        "rank": i+1, "res1": p["r1"], "res2": p["r2"],
        "dG_total": round(p["total"], 4),
        "dG_vdw": round(p["vdw"], 4),
        "dG_elec": round(p["elec"], 4),
        "dG_polar": round(p["polar"], 4),
        "is_relay": p["is_relay"],
    })

for p in relay_iface:
    results["relay_pairs"].append({
        "res1": p["r1"], "res2": p["r2"],
        "dG_total": round(p["total"], 4),
        "dG_vdw": round(p["vdw"], 4),
        "dG_elec": round(p["elec"], 4),
        "dG_polar": round(p["polar"], 4),
    })

with open(os.path.join(DATA_ROOT, "mmpbsa_idecomp3/idecomp3_complex_pairwise.json"), "w") as f:
    json.dump(results, f, indent=2)
print(f"\nResults saved to idecomp3_complex_pairwise.json")
