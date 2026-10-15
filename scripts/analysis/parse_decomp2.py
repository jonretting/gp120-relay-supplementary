import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Parse MMPBSA decomposition: compute Delta = Complex - Receptor - Ligand."""
import json
import os
from pathlib import Path

DECOMP = Path(os.path.join(DATA_ROOT, "mmpbsa/FINAL_DECOMP.dat"))
OUTDIR = Path(os.path.join(DATA_ROOT, "results", "08-binding-energetics"))

with open(DECOMP) as f:
    content = f.read()

def parse_section(text):
    """Parse one section of the decomp file into a dict of residue data."""
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

print(f"Complex residues: {len(complex_res)}")
print(f"Receptor residues: {len(receptor_res)}")
print(f"Ligand residues: {len(ligand_res)}")

# Compute Delta for each residue
# For receptor residues: Delta = Complex(res) - Receptor(res) (ligand contrib = 0)
# For ligand residues: Delta = Complex(res) - Ligand(res) (receptor contrib = 0)
deltas = []

for res_id in complex_res:
    comp = complex_res[res_id]
    
    # Determine if receptor or ligand
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

print(f"Delta residues computed: {len(deltas)}")

# Sort by absolute total contribution
sorted_all = sorted(deltas, key=lambda r: abs(r["d_total"]), reverse=True)

print("\n" + "=" * 100)
print("TOP 20 PER-RESIDUE ΔG CONTRIBUTORS (|ΔG| ranked)")
print("=" * 100)
print(f"{'Rank':>4s} {'Type':>4s} {'Residue':>6s} {'#':>4s}  {'ΔG_total':>10s} {'ΔVdW':>8s} {'ΔElec':>8s} {'ΔPolar':>8s} {'ΔNonPol':>8s}  {'Note':>10s}")
print("-" * 100)
for i, r in enumerate(sorted_all[:20]):
    relay = "**RELAY**" if r["chain_type"] == "R" and r["resnum"] in range(372, 380) else ""
    label = "gp120" if r["chain_type"] == "R" else "CD4"
    print(f"{i+1:4d} {label:>5s} {r['resname']:>6s} {r['resnum']:4d}  {r['d_total']:10.3f} {r['d_vdw']:8.3f} {r['d_elec']:8.3f} {r['d_polar']:8.3f} {r['d_nonpolar']:8.3f}  {relay}")

# Relay pocket residues
print("\n" + "=" * 100)
print("RELAY POCKET RESIDUES (gp120 res 372-379 = 6MEO 465-472: RPGGGDMR)")
print("=" * 100)
pocket_map = {372:"ARG(465)", 373:"PRO(466)", 374:"GLY(467)", 375:"GLY(468)", 
              376:"GLY(469)", 377:"ASP(470)", 378:"MET(471)", 379:"ARG(472)"}

pocket_residues = [r for r in deltas if r["chain_type"] == "R" and r["resnum"] in range(372, 380)]
pocket_total = 0

print(f"{'Residue':>8s} {'6MEO#':>10s} {'ΔG_total':>10s} {'ΔVdW':>8s} {'ΔElec':>8s} {'ΔPolar':>8s} {'ΔNonPol':>8s}")
print("-" * 70)
for r in sorted(pocket_residues, key=lambda x: x["resnum"]):
    label = pocket_map.get(r["resnum"], f"{r['resname']}({r['resnum']+396})")
    print(f"{r['resname']:>8s} {label:>10s} {r['d_total']:10.4f} {r['d_vdw']:8.4f} {r['d_elec']:8.4f} {r['d_polar']:8.4f} {r['d_nonpolar']:8.4f}")
    pocket_total += r["d_total"]
print(f"{'TOTAL':>8s} {'':>10s} {pocket_total:10.4f}")

# Compare with 1.1a results
print("\n" + "=" * 100)
print("COMPARISON: 1.1a (reduced) vs 1.1b (ternary)")
print("=" * 100)
print(f"  1.1a overall ΔG: -31.37 kcal/mol (reduced gp120-CD4, 100ns)")
print(f"  1.1b overall ΔG: -38.02 kcal/mol (ternary gp120-CD4+CCR5, 20ns)")
print(f"  Difference: {-38.02 - (-31.37):.2f} kcal/mol (stronger in ternary context)")
print(f"\n  1.1a relay pocket total: -3.07 kcal/mol")
print(f"  1.1b relay pocket total: {pocket_total:.2f} kcal/mol")

# Save
output = {
    "simulation": "20ns full ternary complex MM-PBSA (1.1b)",
    "comparison_with_1_1a": {
        "overall_dG_1_1a_kcal": -31.37,
        "overall_dG_1_1b_kcal": -38.02,
        "relay_pocket_dG_1_1a_kcal": -3.07,
        "relay_pocket_dG_1_1b_kcal": round(pocket_total, 4),
    },
    "top20_contributors": [
        {
            "rank": i+1,
            "residue": r["resname"],
            "resnum_reduced": r["resnum"],
            "chain": "gp120" if r["chain_type"] == "R" else "CD4",
            "dG_total": r["d_total"],
            "dG_vdw": r["d_vdw"],
            "dG_elec": r["d_elec"],
            "is_relay": r["chain_type"] == "R" and r["resnum"] in range(372, 380),
        }
        for i, r in enumerate(sorted_all[:20])
    ],
    "relay_pocket": {
        "residues": [
            {
                "resname": r["resname"],
                "resnum_reduced": r["resnum"],
                "resnum_6meo": r["resnum"] + 93,
                "dG_total": r["d_total"],
                "dG_vdw": r["d_vdw"],
                "dG_elec": r["d_elec"],
                "dG_polar": r["d_polar"],
                "dG_nonpolar": r["d_nonpolar"],
            }
            for r in sorted(pocket_residues, key=lambda x: x["resnum"])
        ],
        "total_dG": round(pocket_total, 4),
    },
}

with open(OUTDIR / "mmpbsa_ternary_results.json", "w") as f:
    json.dump(output, f, indent=2)

print(f"\nSaved to {OUTDIR / 'mmpbsa_ternary_results.json'}")
