import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Parse the DELTAS section directly from the idecomp=1 DECOMP file."""
import json
import os
from pathlib import Path

DECOMP = Path(os.path.join(DATA_ROOT, "mmpbsa_idecomp1/FINAL_DECOMP_id1.dat"))
OUTDIR = Path(os.path.join(DATA_ROOT, "results", "08-binding-energetics"))

with open(DECOMP) as f:
    content = f.read()

# Find the DELTAS section
delta_start = content.find("DELTAS:")
if delta_start < 0:
    print("ERROR: No DELTAS section found")
    exit(1)

delta_text = content[delta_start:]

# Parse DELTAS Total Energy Decomposition
residues = []
in_total = False
for line in delta_text.split("\n"):
    if "Total Energy Decomposition:" in line:
        in_total = True
        continue
    if "Sidechain Energy Decomposition:" in line or "Backbone Energy Decomposition:" in line:
        in_total = False
        continue
    if not in_total:
        continue
    if line.startswith("Residue,") or line.startswith(",Avg."):
        continue
    if not line.strip():
        continue
    
    parts = line.strip().split(",")
    if len(parts) < 18:
        continue
    
    res_id = parts[0]
    try:
        id_parts = res_id.split(":")
        chain_type = id_parts[0]
        chain = id_parts[1]
        resname = id_parts[2]
        resnum = int(id_parts[3])
        
        data = {
            "id": res_id,
            "chain_type": chain_type,
            "chain": chain,
            "resname": resname,
            "resnum": resnum,
            "d_internal": float(parts[1]),
            "d_vdw": float(parts[3]),
            "d_elec": float(parts[6]),
            "d_polar": float(parts[9]),
            "d_nonpolar": float(parts[12]),
            "d_total": float(parts[15]),
        }
        residues.append(data)
    except (ValueError, IndexError) as e:
        continue

print(f"Parsed {len(residues)} Delta residues")

# Sum check
total_sum = sum(r["d_total"] for r in residues)
vdw_sum = sum(r["d_vdw"] for r in residues)
elec_sum = sum(r["d_elec"] for r in residues)
polar_sum = sum(r["d_polar"] for r in residues)
print(f"\nSum of per-residue ΔG: {total_sum:.4f} kcal/mol (overall: -38.02)")
print(f"Sum ΔVdW: {vdw_sum:.4f}, ΔElec: {elec_sum:.4f}, ΔPolar: {polar_sum:.4f}")

# Sort by |total|
sorted_all = sorted(residues, key=lambda r: abs(r["d_total"]), reverse=True)

print("\n" + "=" * 110)
print("TOP 30 PER-RESIDUE ΔG CONTRIBUTORS -- idecomp=1 DELTAS section")
print("=" * 110)
print(f"{'Rank':>4s} {'Type':>5s} {'Res':>6s} {'#':>4s}  {'ΔG_total':>10s} {'ΔInt':>10s} {'ΔVdW':>10s} {'ΔElec':>10s} {'ΔPolar':>10s}  {'Note':>10s}")
print("-" * 110)
for i, r in enumerate(sorted_all[:30]):
    relay = "**RELAY**" if r["chain_type"] == "R" and r["resnum"] in range(372, 380) else ""
    label = "gp120" if r["chain_type"] == "R" else "CD4"
    print(f"{i+1:4d} {label:>5s} {r['resname']:>6s} {r['resnum']:4d}  {r['d_total']:10.4f} {r['d_internal']:10.4f} {r['d_vdw']:10.4f} {r['d_elec']:10.4f} {r['d_polar']:10.4f}  {relay}")

# Relay pocket
print("\n" + "=" * 110)
print("RELAY POCKET (gp120 res 372-379 = RPGGGDMR = 6MEO 465-472) -- idecomp=1 DELTAS")
print("=" * 110)
pocket_map = {372:"R465", 373:"P466", 374:"G467", 375:"G468",
              376:"G469", 377:"D470", 378:"M471", 379:"R472"}
pocket = [r for r in residues if r["chain_type"] == "R" and r["resnum"] in range(372, 380)]
pocket_total = 0
print(f"{'Res':>6s} {'6MEO':>6s} {'ΔG_total':>10s} {'ΔVdW':>10s} {'ΔElec':>10s} {'ΔPolar':>10s} {'ΔNonPol':>10s}")
print("-" * 70)
for r in sorted(pocket, key=lambda x: x["resnum"]):
    label = pocket_map.get(r["resnum"], "?")
    print(f"{r['resname']:>6s} {label:>6s} {r['d_total']:10.4f} {r['d_vdw']:10.4f} {r['d_elec']:10.4f} {r['d_polar']:10.4f} {r['d_nonpolar']:10.4f}")
    pocket_total += r["d_total"]
print(f"{'TOTAL':>6s} {'':>6s} {pocket_total:10.4f}")

# Save
output = {
    "method": "idecomp=1 DELTAS (direct from gmx_MMPBSA)",
    "overall_dG_kcal": -38.02,
    "sum_per_residue_dG_kcal": round(total_sum, 4),
    "relay_pocket_dG_kcal": round(pocket_total, 4),
    "top30": [
        {
            "rank": i+1, "residue": r["resname"], "resnum_gromacs": r["resnum"],
            "chain": "gp120" if r["chain_type"] == "R" else "CD4",
            "dG_total": round(r["d_total"], 4),
            "dG_vdw": round(r["d_vdw"], 4),
            "dG_elec": round(r["d_elec"], 4),
            "is_relay": r["chain_type"] == "R" and r["resnum"] in range(372, 380),
        }
        for i, r in enumerate(sorted_all[:30])
    ],
    "relay_pocket": [
        {
            "resname": r["resname"],
            "resnum_gromacs": r["resnum"],
            "resnum_6meo": r["resnum"] + 93,
            "dG_total": round(r["d_total"], 4),
            "dG_vdw": round(r["d_vdw"], 4),
            "dG_elec": round(r["d_elec"], 4),
            "dG_polar": round(r["d_polar"], 4),
        }
        for r in sorted(pocket, key=lambda x: x["resnum"])
    ],
}
with open(OUTDIR / "idecomp1_delta_results.json", "w") as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {OUTDIR / 'idecomp1_delta_results.json'}")
