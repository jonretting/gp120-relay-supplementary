#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
FEP Non-Competitiveness Analysis

Analyzes dH/dlambda files from GROMACS FEP simulations using BAR and MBAR
estimators via alchemlyb.

Computes:
  ΔG_drug(ternary):  drug decoupling from gp120*CD4*drug
  ΔG_drug(binary):   drug decoupling from gp120*drug
  ΔΔG = ΔG_drug(ternary) - ΔG_drug(binary)

If |ΔΔG| < 1 kcal/mol -> non-competitive mechanism confirmed.
"""

import os
import sys
import json
import glob
import argparse
import warnings
from pathlib import Path

import numpy as np

DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)

try:
    import alchemlyb
    from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk
    from alchemlyb.estimators import BAR, MBAR
    from alchemlyb.preprocessing import decorrelate_dhdl, decorrelate_u_nk
    HAS_ALCHEMLYB = True
except ImportError:
    HAS_ALCHEMLYB = False
    print("WARNING: alchemlyb not installed. Install with: pip install alchemlyb")

KJ_TO_KCAL = 0.239006


def collect_xvg_files(leg_dir, replica, n_lambda=32, phase="prod"):
    """Collect dhdl.xvg files for all lambda windows of one replica."""
    files = []
    for lam in range(n_lambda):
        win_dir = leg_dir / f"rep{replica}" / f"lambda_{lam:02d}"
        xvg = win_dir / f"{phase}.xvg"
        if not xvg.exists():
            alt_xvg = win_dir / f"dhdl_{phase}.xvg"
            if alt_xvg.exists():
                xvg = alt_xvg
            else:
                xvg_glob = list(win_dir.glob("*.xvg"))
                dhdl_files = [f for f in xvg_glob if "dhdl" in f.name.lower()]
                if dhdl_files:
                    xvg = dhdl_files[0]
                else:
                    print(f"  WARNING: No xvg found for {leg_dir.name}/rep{replica}/lambda_{lam:02d}")
                    continue
        files.append(str(xvg))
    return files


def analyze_leg_bar(xvg_files):
    """Analyze one leg using BAR estimator."""
    dhdl_data = []
    for f in xvg_files:
        try:
            dhdl = extract_dHdl(f, T=300)
            dhdl_data.append(dhdl)
        except Exception as e:
            print(f"  WARNING: Could not parse {f}: {e}")
            continue

    if not dhdl_data:
        return None, None

    import pandas as pd
    combined = pd.concat(dhdl_data)
    combined = decorrelate_dhdl(combined)

    bar = BAR()
    bar.fit(combined)

    delta_g = bar.delta_f_.iloc[0, -1]
    delta_g_err = bar.d_delta_f_.iloc[0, -1]

    return delta_g, delta_g_err


def analyze_leg_mbar(xvg_files):
    """Analyze one leg using MBAR estimator."""
    u_nk_data = []
    for f in xvg_files:
        try:
            u_nk = extract_u_nk(f, T=300)
            u_nk_data.append(u_nk)
        except Exception as e:
            print(f"  WARNING: Could not parse {f}: {e}")
            continue

    if not u_nk_data:
        return None, None

    import pandas as pd
    combined = pd.concat(u_nk_data)
    combined = decorrelate_u_nk(combined)

    mbar = MBAR()
    mbar.fit(combined)

    delta_g = mbar.delta_f_.iloc[0, -1]
    delta_g_err = mbar.d_delta_f_.iloc[0, -1]

    return delta_g, delta_g_err


def analyze_all(base_dir, legs=None, n_replicas=3, n_lambda=32):
    """Run full analysis across all legs and replicas."""
    runs_dir = Path(base_dir) / "runs"
    if legs is None:
        legs = ["ternary", "binary", "solvent"]

    results = {}

    for leg in legs:
        leg_dir = runs_dir / leg
        if not leg_dir.exists():
            print(f"Skipping {leg} -- directory not found")
            continue

        print(f"\n{'='*50}")
        print(f"Analyzing: {leg}")
        print(f"{'='*50}")

        leg_results = {"bar": [], "mbar": []}

        for rep in range(1, n_replicas + 1):
            xvg_files = collect_xvg_files(leg_dir, rep, n_lambda)
            if not xvg_files:
                print(f"  rep{rep}: No data found")
                continue

            print(f"  rep{rep}: {len(xvg_files)} xvg files")

            dg_bar, err_bar = analyze_leg_bar(xvg_files)
            dg_mbar, err_mbar = analyze_leg_mbar(xvg_files)

            if dg_bar is not None:
                leg_results["bar"].append({"dG": float(dg_bar), "err": float(err_bar)})
                print(f"    BAR:  ΔG = {dg_bar:.3f} ± {err_bar:.3f} kJ/mol"
                      f"  ({dg_bar*KJ_TO_KCAL:.3f} ± {err_bar*KJ_TO_KCAL:.3f} kcal/mol)")

            if dg_mbar is not None:
                leg_results["mbar"].append({"dG": float(dg_mbar), "err": float(err_mbar)})
                print(f"    MBAR: ΔG = {dg_mbar:.3f} ± {err_mbar:.3f} kJ/mol"
                      f"  ({dg_mbar*KJ_TO_KCAL:.3f} ± {err_mbar*KJ_TO_KCAL:.3f} kcal/mol)")

        if leg_results["mbar"]:
            dgs = [r["dG"] for r in leg_results["mbar"]]
            mean_dg = np.mean(dgs)
            sem_dg = np.std(dgs) / np.sqrt(len(dgs)) if len(dgs) > 1 else leg_results["mbar"][0]["err"]
            leg_results["mean_mbar_kj"] = float(mean_dg)
            leg_results["sem_mbar_kj"] = float(sem_dg)
            leg_results["mean_mbar_kcal"] = float(mean_dg * KJ_TO_KCAL)
            leg_results["sem_mbar_kcal"] = float(sem_dg * KJ_TO_KCAL)
            print(f"\n  Mean MBAR: ΔG = {mean_dg:.3f} ± {sem_dg:.3f} kJ/mol"
                  f"  ({mean_dg*KJ_TO_KCAL:.3f} ± {sem_dg*KJ_TO_KCAL:.3f} kcal/mol)")

        results[leg] = leg_results

    # Compute ΔΔG
    print(f"\n{'='*60}")
    print("THERMODYNAMIC CYCLE RESULTS")
    print(f"{'='*60}")

    if "ternary" in results and "binary" in results:
        for method in ["mbar", "bar"]:
            key_mean = f"mean_{method}_kcal"
            key_sem = f"sem_{method}_kcal"
            if key_mean in results["ternary"] and key_mean in results["binary"]:
                dg_tern = results["ternary"][key_mean]
                dg_bin = results["binary"][key_mean]
                ddg = dg_tern - dg_bin
                err_tern = results["ternary"][key_sem]
                err_bin = results["binary"][key_sem]
                ddg_err = np.sqrt(err_tern**2 + err_bin**2)

                results[f"ddG_{method}_kcal"] = float(ddg)
                results[f"ddG_{method}_err_kcal"] = float(ddg_err)

                print(f"\n{method.upper()} results:")
                print(f"  ΔG_drug(ternary) = {dg_tern:.3f} ± {err_tern:.3f} kcal/mol")
                print(f"  ΔG_drug(binary)  = {dg_bin:.3f} ± {err_bin:.3f} kcal/mol")
                print(f"  ΔΔG = {ddg:.3f} ± {ddg_err:.3f} kcal/mol")
                print()
                if abs(ddg) < 1.0:
                    print(f"  [OK] |ΔΔG| = {abs(ddg):.3f} < 1.0 kcal/mol -> NON-COMPETITIVE")
                elif abs(ddg) < 2.0:
                    print(f"  ~ |ΔΔG| = {abs(ddg):.3f} -- borderline, may need more sampling")
                else:
                    print(f"  [X] |ΔΔG| = {abs(ddg):.3f} > 2.0 kcal/mol -- significant coupling detected")

    output_file = Path(base_dir) / "fep_results.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_file}")

    return results


def main():
    parser = argparse.ArgumentParser(description="Analyze FEP non-competitiveness results")
    parser.add_argument("--base-dir", default=os.path.join(DATA_ROOT, "results", "06-fep"))
    parser.add_argument("--legs", nargs="+", default=None)
    parser.add_argument("--replicas", type=int, default=3)
    parser.add_argument("--n-lambda", type=int, default=32)
    args = parser.parse_args()

    if not HAS_ALCHEMLYB:
        print("ERROR: alchemlyb is required. Install with:")
        print("  pip install alchemlyb")
        sys.exit(1)

    analyze_all(args.base_dir, legs=args.legs, n_replicas=args.replicas, n_lambda=args.n_lambda)


if __name__ == "__main__":
    main()
