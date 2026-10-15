#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Parse PBC-corrected GROMACS XVG outputs into a validation summary."""
import os, json, sys
import numpy as np
from pathlib import Path

valdir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")

def parse_xvg(filepath):
    data = []
    with open(filepath) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    data.append([float(x) for x in parts])
                except ValueError:
                    continue
    return data

results = {}

# RMSD
rmsd_file = valdir / "rmsd_pbc.xvg"
if rmsd_file.exists():
    rmsd_data = parse_xvg(str(rmsd_file))
    times = [d[0] for d in rmsd_data]
    rmsds = [d[1] for d in rmsd_data]
    rmsds_A = [r * 10 for r in rmsds]
    q4 = rmsds_A[-len(rmsds_A)//4:]
    results["backbone_rmsd"] = {
        "mean_A": round(np.mean(rmsds_A), 2),
        "std_A": round(np.std(rmsds_A), 2),
        "max_A": round(np.max(rmsds_A), 2),
        "final_A": round(rmsds_A[-1], 2),
        "last_quarter_mean_A": round(np.mean(q4), 2),
        "last_quarter_std_A": round(np.std(q4), 2),
        "plateau": np.std(q4) < 0.3,
        "n_points": len(rmsds_A),
    }
    print(f"Backbone RMSD: {np.mean(rmsds_A):.2f} +/- {np.std(rmsds_A):.2f} A  (max {np.max(rmsds_A):.2f}, final {rmsds_A[-1]:.2f})")
    print(f"  Last quarter: {np.mean(q4):.2f} +/- {np.std(q4):.2f} A  plateau={'YES' if np.std(q4) < 0.3 else 'NO'}")

# RMSF
rmsf_file = valdir / "rmsf_pbc.xvg"
if rmsf_file.exists():
    rmsf_data = parse_xvg(str(rmsf_file))
    resnums = [int(d[0]) for d in rmsf_data]
    rmsf_vals = [d[1] * 10 for d in rmsf_data]  # nm -> A
    
    print(f"\nPer-residue RMSF: mean={np.mean(rmsf_vals):.2f} A")
    
    # Check relay pocket across numberings
    for label, pocket_range in [("6MEO(465-472)", range(465, 473)), ("SOSIP(469-476)", range(469, 477)), ("reduced(69-76)", range(69, 77))]:
        pocket = [(r, v) for r, v in zip(resnums, rmsf_vals) if r in pocket_range]
        if pocket:
            pvals = [v for _, v in pocket]
            results[f"rmsf_pocket_{label.split('(')[0]}"] = {
                "residues": [{"resnum": r, "rmsf_A": round(v, 2)} for r, v in pocket],
                "mean_A": round(np.mean(pvals), 2),
            }
            print(f"  Pocket {label}: mean RMSF = {np.mean(pvals):.2f} A")
            for r, v in pocket:
                print(f"    Res {r}: {v:.2f} A")
    
    # Top 5 most flexible
    sorted_rmsf = sorted(zip(resnums, rmsf_vals), key=lambda x: x[1], reverse=True)
    print(f"\n  Top 5 most flexible:")
    for r, v in sorted_rmsf[:5]:
        print(f"    Res {r}: {v:.2f} A")
    results["rmsf_top5"] = [{"resnum": r, "rmsf_A": round(v, 2)} for r, v in sorted_rmsf[:5]]
    results["rmsf_mean_all_A"] = round(np.mean(rmsf_vals), 2)

# Rg
rg_file = valdir / "gyrate_pbc.xvg"
if rg_file.exists():
    rg_data = parse_xvg(str(rg_file))
    rg_vals = [d[1] * 10 for d in rg_data]  # nm -> A
    results["radius_of_gyration"] = {
        "mean_A": round(np.mean(rg_vals), 2),
        "std_A": round(np.std(rg_vals), 2),
    }
    print(f"\nRg: {np.mean(rg_vals):.2f} +/- {np.std(rg_vals):.2f} A")

# Save
with open(valdir / "validation_pbc_corrected.json", "w") as f:
    json.dump(results, f, indent=2)

print(f"\nSaved to {valdir / 'validation_pbc_corrected.json'}")
