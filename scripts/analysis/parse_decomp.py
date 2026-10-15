import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Parse MMPBSA decomposition results for ternary complex."""
import json
import os
from pathlib import Path

DECOMP = Path(os.path.join(DATA_ROOT, "mmpbsa/FINAL_DECOMP.dat"))
RESULTS = Path(os.path.join(DATA_ROOT, "mmpbsa/FINAL_RESULTS.dat"))
OUTDIR = Path(os.path.join(DATA_ROOT, "results", "08-binding-energetics"))

with open(DECOMP) as f:
    content = f.read()

# Find the Delta section
sections = content.split("Delta:")
if len(sections) < 2:
    print("ERROR: No Delta section found")
    print(f"Sections found: {len(sections)}")
    # Look for alternative markers
    for line in content.split("\n")[:20]:
        print(f"  {line}")
    exit(1)

delta_text = sections[1]

# Parse the Delta Total Energy Decomposition
# Format: Residue,Internal,,,van der Waals,,,Electrostatic,,,Polar Solvation,,,Non-Polar Solv.,,,TOTAL,,
#         ,Avg.,Std. Dev.,SEM,Avg.,Std. Dev.,SEM,...

residues = []
parsing = False
for line in delta_text.split("\n"):
    if "Total Energy Decomposition:" in line:
        parsing = True
        continue
    if "Sidechain Energy Decomposition:" in line:
        break
    if not parsing:
        continue
    if line.startswith("Residue,"):
        continue
    if line.startswith(",Avg."):
        continue
    if not line.strip():
        continue
    
    parts = line.strip().split(",")
    if len(parts) < 18:
        continue
    
    residue_id = parts[0]
    try:
        # Parse: Internal(avg,std,sem), VdW(avg,std,sem), Elec(avg,std,sem), Polar(avg,std,sem), NonPolar(avg,std,sem), Total(avg,std,sem)
        total_avg = float(parts[15])
        total_std = float(parts[16])
        vdw_avg = float(parts[3])
        elec_avg = float(parts[6])
        polar_avg = float(parts[9])
        nonpolar_avg = float(parts[12])
        
        # Parse residue ID: R:A:ASP:1 or L:B:LYS:1
        id_parts = residue_id.split(":")
        if len(id_parts) >= 4:
            chain_type = id_parts[0]  # R=receptor, L=ligand
            chain_letter = id_parts[1]
            resname = id_parts[2]
            resnum = int(id_parts[3])
        else:
            chain_type = "?"
            chain_letter = "?"
            resname = residue_id
            resnum = 0
        
        residues.append({
            "id": residue_id,
            "chain_type": chain_type,
            "chain": chain_letter,
            "resname": resname,
            "resnum": resnum,
            "total_avg": total_avg,
            "total_std": total_std,
            "vdw_avg": vdw_avg,
            "elec_avg": elec_avg,
            "polar_avg": polar_avg,
            "nonpolar_avg": nonpolar_avg,
        })
    except (ValueError, IndexError) as e:
        continue

print(f"Parsed {len(residues)} residues from Delta decomposition")

# Sort by absolute contribution
sorted_all = sorted(residues, key=lambda r: abs(r["total_avg"]), reverse=True)

# Print top 20 overall contributors
print("\n" + "=" * 90)
print("TOP 20 PER-RESIDUE CONTRIBUTORS (|ΔG| ranked)")
print("=" * 90)
print(f"{'Rank':>4s} {'Residue':>10s} {'Chain':>6s} {'#':>4s} {'ΔG_total':>10s} {'VdW':>8s} {'Elec':>8s} {'Polar':>8s} {'NonPol':>8s}")
print("-" * 90)
for i, r in enumerate(sorted_all[:20]):
    relay = " *RELAY*" if r["chain_type"] == "R" and r["resnum"] in range(69, 77) else ""
    print(f"{i+1:4d} {r['resname']:>10s} {r['chain_type']+':'+r['chain']:>6s} {r['resnum']:4d} {r['total_avg']:10.2f} {r['vdw_avg']:8.2f} {r['elec_avg']:8.2f} {r['polar_avg']:8.2f} {r['nonpolar_avg']:8.2f}{relay}")

# Relay pocket residues (69-76 in receptor)
print("\n" + "=" * 90)
print("RELAY POCKET RESIDUES (gp120 res 69-76, maps to 6MEO 465-472)")
print("=" * 90)
pocket_residues = [r for r in residues if r["chain_type"] == "R" and r["resnum"] in range(69, 77)]
pocket_total = 0
print(f"{'Residue':>10s} {'#':>4s} {'ΔG_total':>10s} {'VdW':>8s} {'Elec':>8s} {'Polar':>8s} {'NonPol':>8s}")
print("-" * 70)
for r in sorted(pocket_residues, key=lambda x: x["resnum"]):
    print(f"{r['resname']:>10s} {r['resnum']:4d} {r['total_avg']:10.4f} {r['vdw_avg']:8.4f} {r['elec_avg']:8.4f} {r['polar_avg']:8.4f} {r['nonpolar_avg']:8.4f}")
    pocket_total += r["total_avg"]
print(f"{'TOTAL':>10s} {'':>4s} {pocket_total:10.4f}")

# Save as JSON
output = {
    "simulation": "20ns full ternary complex (gp120-CD4 with CCR5)",
    "overall_binding": {
        "delta_total_kcal": -38.02,
        "delta_vdw_kcal": -86.92,
        "delta_eel_kcal": -157.15,
        "delta_gb_kcal": 218.56,
        "delta_esurf_kcal": -12.51,
    },
    "per_residue_top20": [
        {
            "rank": i+1,
            "residue": r["resname"],
            "resnum": r["resnum"],
            "chain_type": r["chain_type"],
            "chain": r["chain"],
            "dG_total": round(r["total_avg"], 4),
            "dG_vdw": round(r["vdw_avg"], 4),
            "dG_elec": round(r["elec_avg"], 4),
            "is_relay": r["chain_type"] == "R" and r["resnum"] in range(69, 77),
        }
        for i, r in enumerate(sorted_all[:20])
    ],
    "relay_pocket": {
        "residues": [
            {
                "resname": r["resname"],
                "resnum_reduced": r["resnum"],
                "resnum_6meo": r["resnum"] + 396,
                "dG_total": round(r["total_avg"], 4),
                "dG_vdw": round(r["vdw_avg"], 4),
                "dG_elec": round(r["elec_avg"], 4),
                "dG_polar": round(r["polar_avg"], 4),
                "dG_nonpolar": round(r["nonpolar_avg"], 4),
            }
            for r in sorted(pocket_residues, key=lambda x: x["resnum"])
        ],
        "total_dG": round(pocket_total, 4),
    },
}

with open(OUTDIR / "mmpbsa_ternary_results.json", "w") as f:
    json.dump(output, f, indent=2)

print(f"\nResults saved to {OUTDIR / 'mmpbsa_ternary_results.json'}")
