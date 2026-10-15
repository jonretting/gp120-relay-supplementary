import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Parse idecomp=1 MMPBSA decomposition for the ternary complex."""
import json
import os
from pathlib import Path

DECOMP = Path(os.path.join(DATA_ROOT, "mmpbsa_idecomp1/FINAL_DECOMP_id1.dat"))
OUTDIR = Path(os.path.join(DATA_ROOT, "results", "08-binding-energetics"))

with open(DECOMP) as f:
    content = f.read()

def parse_section(text):
    """Parse one section into dict of residue data (Total Energy only)."""
    residues = {}
    in_total = False
    for line in text.split("\n"):
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
            data = {
                "internal": float(parts[1]),
                "vdw": float(parts[3]),
                "elec": float(parts[6]),
                "polar": float(parts[9]),
                "nonpolar": float(parts[12]),
                "total": float(parts[15]),
            }
            residues[res_id] = data
        except (ValueError, IndexError):
            continue
    return residues

# Split into sections
parts = content.split("Receptor:")
complex_text = parts[0]
rest = parts[1] if len(parts) > 1 else ""
parts2 = rest.split("Ligand:")
receptor_text = parts2[0]
ligand_text = parts2[1] if len(parts2) > 1 else ""

complex_res = parse_section(complex_text)
receptor_res = parse_section(receptor_text)
ligand_res = parse_section(ligand_text)

print(f"Complex: {len(complex_res)}, Receptor: {len(receptor_res)}, Ligand: {len(ligand_res)}")

# Compute Delta
deltas = []
for res_id in complex_res:
    comp = complex_res[res_id]
    id_parts = res_id.split(":")
    chain_type = id_parts[0] if len(id_parts) >= 4 else "?"
    chain = id_parts[1] if len(id_parts) >= 4 else "?"
    resname = id_parts[2] if len(id_parts) >= 4 else res_id
    resnum = int(id_parts[3]) if len(id_parts) >= 4 else 0
    
    if chain_type == "R" and res_id in receptor_res:
        ref = receptor_res[res_id]
    elif chain_type == "L" and res_id in ligand_res:
        ref = ligand_res[res_id]
    else:
        continue
    
    delta = {
        "id": res_id,
        "chain_type": chain_type,
        "chain": chain,
        "resname": resname,
        "resnum": resnum,
        "d_internal": round(comp["internal"] - ref["internal"], 4),
        "d_vdw": round(comp["vdw"] - ref["vdw"], 4),
        "d_elec": round(comp["elec"] - ref["elec"], 4),
        "d_polar": round(comp["polar"] - ref["polar"], 4),
        "d_nonpolar": round(comp["nonpolar"] - ref["nonpolar"], 4),
        "d_total": round(comp["total"] - ref["total"], 4),
    }
    deltas.append(delta)

# Compute sum to check if it matches overall ΔG
total_sum = sum(d["d_total"] for d in deltas)
print(f"\nSum of all per-residue ΔG: {total_sum:.2f} kcal/mol (should be ~-38.02)")

# Sort by absolute contribution
sorted_all = sorted(deltas, key=lambda r: abs(r["d_total"]), reverse=True)

print("\n" + "=" * 100)
print(f"TOP 30 PER-RESIDUE ΔG CONTRIBUTORS (|ΔG| ranked) -- idecomp=1")
print("=" * 100)
print(f"{'Rank':>4s} {'Type':>5s} {'Res':>6s} {'#':>4s}  {'ΔG_total':>10s} {'ΔIntrnl':>8s} {'ΔVdW':>8s} {'ΔElec':>10s} {'ΔPolar':>10s}  {'Note':>10s}")
print("-" * 100)
for i, r in enumerate(sorted_all[:30]):
    relay = "**RELAY**" if r["chain_type"] == "R" and r["resnum"] in range(372, 380) else ""
    label = "gp120" if r["chain_type"] == "R" else "CD4"
    print(f"{i+1:4d} {label:>5s} {r['resname']:>6s} {r['resnum']:4d}  {r['d_total']:10.3f} {r['d_internal']:8.3f} {r['d_vdw']:8.3f} {r['d_elec']:10.3f} {r['d_polar']:10.3f}  {relay}")

# Relay pocket 
print("\n" + "=" * 100)
print("RELAY POCKET (gp120 res 372-379 = 6MEO 465-472: RPGGGDMR) -- idecomp=1")
print("=" * 100)
pocket_map = {372:"R(465)", 373:"P(466)", 374:"G(467)", 375:"G(468)", 
              376:"G(469)", 377:"D(470)", 378:"M(471)", 379:"R(472)"}
pocket_residues = [r for r in deltas if r["chain_type"] == "R" and r["resnum"] in range(372, 380)]
pocket_total = 0
print(f"{'Res':>6s} {'6MEO':>8s} {'ΔG_total':>10s} {'ΔIntrn':>8s} {'ΔVdW':>8s} {'ΔElec':>10s} {'ΔPolar':>10s} {'ΔNonP':>8s}")
print("-" * 80)
for r in sorted(pocket_residues, key=lambda x: x["resnum"]):
    label = pocket_map.get(r["resnum"], f"{r['resname']}")
    print(f"{r['resname']:>6s} {label:>8s} {r['d_total']:10.4f} {r['d_internal']:8.4f} {r['d_vdw']:8.4f} {r['d_elec']:10.4f} {r['d_polar']:10.4f} {r['d_nonpolar']:8.4f}")
    pocket_total += r["d_total"]
print(f"{'TOTAL':>6s} {'':>8s} {pocket_total:10.4f}")

# Save
output = {
    "method": "idecomp=1 per-residue (no 1-4 in VdW/EEL)",
    "overall_dG_kcal": -38.02,
    "sum_per_residue_dG_kcal": round(total_sum, 4),
    "relay_pocket_dG_kcal": round(pocket_total, 4),
    "top30": [
        {
            "rank": i+1, "residue": r["resname"], "resnum": r["resnum"],
            "chain": "gp120" if r["chain_type"] == "R" else "CD4",
            "dG_total": r["d_total"], "dG_vdw": r["d_vdw"], "dG_elec": r["d_elec"],
            "is_relay": r["chain_type"] == "R" and r["resnum"] in range(372, 380),
        }
        for i, r in enumerate(sorted_all[:30])
    ],
    "relay_pocket": [
        {
            "resname": r["resname"], "resnum_gromacs": r["resnum"],
            "resnum_6meo": r["resnum"] + 93,
            "dG_total": r["d_total"], "dG_vdw": r["d_vdw"], "dG_elec": r["d_elec"],
            "dG_polar": r["d_polar"], "dG_nonpolar": r["d_nonpolar"],
        }
        for r in sorted(pocket_residues, key=lambda x: x["resnum"])
    ],
}
with open(OUTDIR / "idecomp1_results.json", "w") as f:
    json.dump(output, f, indent=2)
print(f"\nSaved to {OUTDIR / 'idecomp1_results.json'}")
