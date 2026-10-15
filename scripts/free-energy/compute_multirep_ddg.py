#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Compute multi-rep averaged ddG from all 9 FEP analysis JSONs.

Usage:
  python compute_multirep_ddg.py --indir ./per_rep_results --outdir ./

Expects 9 JSON files in <indir>:
  ternary_rep{1,2,3}_analysis.json
  binary_rep{1,2,3}_analysis.json
  solvent_rep{1,2,3}_analysis.json

These are produced by analyze_fep_rep.py.

Full thermodynamic cycle with solvent correction:
  dG_bind(complex) = dG(solvent) - dG(complex)  [standard: negative = favorable]
  ddG(CD4 effect)  = dG_bind(ternary) - dG_bind(binary)
"""
import argparse
import numpy as np
import json
from pathlib import Path
from datetime import datetime

KJ_TO_KCAL = 0.239006

parser = argparse.ArgumentParser(
    description="Multi-rep FEP thermodynamic cycle: binding FE and CD4 effect")
parser.add_argument("--indir", required=True,
                    help="Directory containing the 9 per-rep analysis JSONs")
parser.add_argument("--outdir", default=".",
                    help="Output directory for fep_multirep_summary.json (default: cwd)")
args = parser.parse_args()

indir = Path(args.indir)

results = {}
for f in sorted(indir.glob("*_rep*_analysis.json")):
    d = json.load(open(f))
    key = f'{d["leg"]}_rep{d["replica"]}'
    results[key] = d

print("=" * 70)
print("COMPLETE FEP INVENTORY (9 legs)")
print("=" * 70)
header = f"{'Leg/Rep':<20} {'TI(dec) kcal/mol':>22} {'MBAR kcal/mol':>22}"
print(header)
print("-" * 70)
for k in sorted(results.keys()):
    r = results[k]
    ti = r["TI_decorrelated"]
    mbar = r.get("MBAR_decorrelated", {})
    mbar_str = f'{mbar["dG_kcal"]:.4f} +/- {mbar["err_kcal"]:.4f}' if mbar else "N/A"
    print(f'{k:<20} {ti["dG_kcal"]:8.4f} +/- {ti["err_kcal"]:.4f}     {mbar_str}')

def leg_stats(leg_name, reps):
    ti_vals = [results[f"{leg_name}_rep{i}"]["TI_decorrelated"]["dG_kcal"] for i in reps]
    mbar_vals = [results[f"{leg_name}_rep{i}"]["MBAR_decorrelated"]["dG_kcal"] for i in reps]
    n = len(reps)
    return {
        "n": n,
        "ti_vals": ti_vals,
        "ti_mean": np.mean(ti_vals),
        "ti_sem": np.std(ti_vals, ddof=1) / np.sqrt(n),
        "ti_sd": np.std(ti_vals, ddof=1),
        "mbar_vals": mbar_vals,
        "mbar_mean": np.mean(mbar_vals),
        "mbar_sem": np.std(mbar_vals, ddof=1) / np.sqrt(n),
        "mbar_sd": np.std(mbar_vals, ddof=1),
    }

tern = leg_stats("ternary", [1, 2, 3])
bina = leg_stats("binary", [1, 2, 3])
solv = leg_stats("solvent", [1, 2, 3])

print()
print("=" * 70)
print("CROSS-REP AVERAGES (3 reps each)")
print("=" * 70)
for name, s in [("Ternary", tern), ("Binary", bina), ("Solvent", solv)]:
    print(f"{name}:")
    print(f"  TI dec:  {s['ti_mean']:.4f} +/- {s['ti_sem']:.4f} SEM  (SD {s['ti_sd']:.4f})  {s['ti_vals']}")
    print(f"  MBAR:    {s['mbar_mean']:.4f} +/- {s['mbar_sem']:.4f} SEM  (SD {s['mbar_sd']:.4f})  {s['mbar_vals']}")

# Decoupling differential (raw, before sign-convention conversion)
print()
print("=" * 70)
print("DECOUPLING COST DIFFERENTIAL")
print("  (higher cost = stronger binding in that environment)")
print("=" * 70)
print()
print(f"  Binary - Solvent:  TI = {bina['ti_mean'] - solv['ti_mean']:+.4f}   MBAR = {bina['mbar_mean'] - solv['mbar_mean']:+.4f} kcal/mol  (drug binding to gp120)")
print(f"  Ternary - Solvent: TI = {tern['ti_mean'] - solv['ti_mean']:+.4f}   MBAR = {tern['mbar_mean'] - solv['mbar_mean']:+.4f} kcal/mol  (drug binding to gp120+CD4)")

# Relative ddG = dG(ternary) - dG(binary)
# Solvent cancels in this difference
ddg_ti = tern["ti_mean"] - bina["ti_mean"]
ddg_ti_err = np.sqrt(tern["ti_sem"]**2 + bina["ti_sem"]**2)
ddg_mbar = tern["mbar_mean"] - bina["mbar_mean"]
ddg_mbar_err = np.sqrt(tern["mbar_sem"]**2 + bina["mbar_sem"]**2)

print()
print("=" * 70)
print("RELATIVE ddG (decoupling convention)")
print("  ternary = gp120 + CD4 + drug")
print("  binary  = gp120 + drug (no CD4)")
print("  ddG_decouple = dG(ternary) - dG(binary)")
print("=" * 70)
print()
print(f"  TI dec:  {ddg_ti:.4f} +/- {ddg_ti_err:.4f} kcal/mol")
print(f"  MBAR:    {ddg_mbar:.4f} +/- {ddg_mbar_err:.4f} kcal/mol")
print()
print("  Interpretation (decoupling: higher dG = stronger binding):")
if ddg_ti < 0:
    cd4_effect = abs(ddg_ti)
    print(f"  >>> ddG < 0: drug is EASIER to remove with CD4 present")
    print(f"  >>> CD4 REDUCES drug binding by {cd4_effect:.2f} kcal/mol")
    print(f"  >>> Drug binds more strongly to gp120 ALONE (binary)")
else:
    print(f"  >>> ddG > 0: drug is HARDER to remove with CD4 present")
    print(f"  >>> CD4 ENHANCES drug binding")

# Standard convention binding free energies
print()
print("=" * 70)
print("BINDING FREE ENERGIES (standard convention: dG_bind = dG_solv - dG_complex)")
print("  Negative = favorable binding")
print("=" * 70)
dg_bind_bin_ti = solv["ti_mean"] - bina["ti_mean"]
dg_bind_bin_ti_err = np.sqrt(solv["ti_sem"]**2 + bina["ti_sem"]**2)
dg_bind_tern_ti = solv["ti_mean"] - tern["ti_mean"]
dg_bind_tern_ti_err = np.sqrt(solv["ti_sem"]**2 + tern["ti_sem"]**2)
dg_bind_bin_mbar = solv["mbar_mean"] - bina["mbar_mean"]
dg_bind_bin_mbar_err = np.sqrt(solv["mbar_sem"]**2 + bina["mbar_sem"]**2)
dg_bind_tern_mbar = solv["mbar_mean"] - tern["mbar_mean"]
dg_bind_tern_mbar_err = np.sqrt(solv["mbar_sem"]**2 + tern["mbar_sem"]**2)
print()
print(f"  gp120 alone (binary):  TI = {dg_bind_bin_ti:.4f} +/- {dg_bind_bin_ti_err:.4f}   MBAR = {dg_bind_bin_mbar:.4f} +/- {dg_bind_bin_mbar_err:.4f} kcal/mol")
print(f"  gp120 + CD4 (ternary): TI = {dg_bind_tern_ti:.4f} +/- {dg_bind_tern_ti_err:.4f}   MBAR = {dg_bind_tern_mbar:.4f} +/- {dg_bind_tern_mbar_err:.4f} kcal/mol")
print()
ddg_bind_ti = dg_bind_tern_ti - dg_bind_bin_ti
ddg_bind_mbar = dg_bind_tern_mbar - dg_bind_bin_mbar
print(f"  ddG_bind (CD4 effect): TI = {ddg_bind_ti:+.4f}   MBAR = {ddg_bind_mbar:+.4f} kcal/mol")
if ddg_bind_ti > 0:
    print(f"  >>> CD4 destabilizes drug binding by +{ddg_bind_ti:.2f} kcal/mol")

print()
print("=" * 70)
print("INTER-REP CONSISTENCY (SD across 3 reps)")
print("=" * 70)
print(f"Ternary  TI SD: {tern['ti_sd']:.4f}   MBAR SD: {tern['mbar_sd']:.4f} kcal/mol")
print(f"Binary   TI SD: {bina['ti_sd']:.4f}   MBAR SD: {bina['mbar_sd']:.4f} kcal/mol")
print(f"Solvent  TI SD: {solv['ti_sd']:.4f}   MBAR SD: {solv['mbar_sd']:.4f} kcal/mol")
print()
for name, s in [("Ternary", tern), ("Binary", bina), ("Solvent", solv)]:
    ti_cv = abs(s["ti_sd"] / s["ti_mean"]) * 100 if s["ti_mean"] != 0 else 0
    mbar_cv = abs(s["mbar_sd"] / s["mbar_mean"]) * 100 if s["mbar_mean"] != 0 else 0
    print(f"{name:10s}  CV(TI): {ti_cv:.2f}%   CV(MBAR): {mbar_cv:.2f}%")

print()
print("=" * 70)
print("TI vs MBAR AGREEMENT")
print("=" * 70)
for name, s in [("Ternary", tern), ("Binary", bina), ("Solvent", solv)]:
    diff = abs(s["ti_mean"] - s["mbar_mean"])
    print(f"{name:10s}  |TI - MBAR| = {diff:.4f} kcal/mol")

# Save comprehensive JSON
timestamp = datetime.now().strftime("%Y-%m-%dT%H:%M")
summary = {
    "timestamp": timestamp,
    "status": "COMPLETE -- all 9 legs analyzed (3x ternary, 3x binary, 3x solvent)",
    "ternary": {
        "n_reps": 3,
        "TI_dec_per_rep_kcal": tern["ti_vals"],
        "TI_dec_mean_kcal": round(tern["ti_mean"], 4),
        "TI_dec_sem_kcal": round(tern["ti_sem"], 4),
        "TI_dec_sd_kcal": round(tern["ti_sd"], 4),
        "MBAR_per_rep_kcal": tern["mbar_vals"],
        "MBAR_mean_kcal": round(tern["mbar_mean"], 4),
        "MBAR_sem_kcal": round(tern["mbar_sem"], 4),
        "MBAR_sd_kcal": round(tern["mbar_sd"], 4),
    },
    "binary": {
        "n_reps": 3,
        "TI_dec_per_rep_kcal": bina["ti_vals"],
        "TI_dec_mean_kcal": round(bina["ti_mean"], 4),
        "TI_dec_sem_kcal": round(bina["ti_sem"], 4),
        "TI_dec_sd_kcal": round(bina["ti_sd"], 4),
        "MBAR_per_rep_kcal": bina["mbar_vals"],
        "MBAR_mean_kcal": round(bina["mbar_mean"], 4),
        "MBAR_sem_kcal": round(bina["mbar_sem"], 4),
        "MBAR_sd_kcal": round(bina["mbar_sd"], 4),
    },
    "solvent": {
        "n_reps": 3,
        "TI_dec_per_rep_kcal": solv["ti_vals"],
        "TI_dec_mean_kcal": round(solv["ti_mean"], 4),
        "TI_dec_sem_kcal": round(solv["ti_sem"], 4),
        "TI_dec_sd_kcal": round(solv["ti_sd"], 4),
        "MBAR_per_rep_kcal": solv["mbar_vals"],
        "MBAR_mean_kcal": round(solv["mbar_mean"], 4),
        "MBAR_sem_kcal": round(solv["mbar_sem"], 4),
        "MBAR_sd_kcal": round(solv["mbar_sd"], 4),
    },
    "system_compositions": {
        "ternary": "gp120 (chain G) + CD4 (chain A) + BMS-626529",
        "binary": "gp120 (chain G) + BMS-626529 (NO CD4)",
        "solvent": "BMS-626529 + water/ions",
        "note": "Verified from topology files. couple-lambda0=vdw-q, couple-lambda1=none (decoupling protocol).",
    },
    "binding_free_energy_standard": {
        "note": "dG_bind = dG_solvent - dG_complex (standard convention, negative = favorable)",
        "binary_gp120_only": {
            "TI_dec_kcal": round(float(solv["ti_mean"] - bina["ti_mean"]), 4),
            "TI_dec_err_kcal": round(float(np.sqrt(solv["ti_sem"]**2 + bina["ti_sem"]**2)), 4),
            "MBAR_kcal": round(float(solv["mbar_mean"] - bina["mbar_mean"]), 4),
            "MBAR_err_kcal": round(float(np.sqrt(solv["mbar_sem"]**2 + bina["mbar_sem"]**2)), 4),
        },
        "ternary_gp120_cd4": {
            "TI_dec_kcal": round(float(solv["ti_mean"] - tern["ti_mean"]), 4),
            "TI_dec_err_kcal": round(float(np.sqrt(solv["ti_sem"]**2 + tern["ti_sem"]**2)), 4),
            "MBAR_kcal": round(float(solv["mbar_mean"] - tern["mbar_mean"]), 4),
            "MBAR_err_kcal": round(float(np.sqrt(solv["mbar_sem"]**2 + tern["mbar_sem"]**2)), 4),
        },
    },
    "cd4_effect_on_drug_binding": {
        "note": "ddG = dG_bind(ternary) - dG_bind(binary). Positive = CD4 reduces drug affinity.",
        "TI_dec_kcal": round(float((solv["ti_mean"] - tern["ti_mean"]) - (solv["ti_mean"] - bina["ti_mean"])), 4),
        "MBAR_kcal": round(float((solv["mbar_mean"] - tern["mbar_mean"]) - (solv["mbar_mean"] - bina["mbar_mean"])), 4),
    },
    "decoupling_ddG_raw": {
        "note": "ddG_decouple = dG(ternary) - dG(binary). Negative = easier to remove drug from ternary.",
        "TI_dec_kcal": round(ddg_ti, 4),
        "TI_dec_err_kcal": round(ddg_ti_err, 4),
        "MBAR_kcal": round(ddg_mbar, 4),
        "MBAR_err_kcal": round(ddg_mbar_err, 4),
    },
    "per_rep_analyses": {k: f"{k}_analysis.json" for k in sorted(results.keys())},
}
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)
out = outdir / "fep_multirep_summary.json"
json.dump(summary, open(out, "w"), indent=2)
print(f"\nSaved to {out}")
