#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Generalized FEP analysis for any leg/rep combination.

Usage:
  python analyze_fep_rep.py ternary 1 --base /path/to/fep/runs
  python analyze_fep_rep.py binary  2 --base /path/to/fep/runs --outdir ./results

Expects FEP data organized as:
  <base>/<system>/rep<N>/lambda_00/prod.xvg
  <base>/<system>/rep<N>/lambda_01/prod.xvg
  ...
  <base>/<system>/rep<N>/lambda_31/prod.xvg

Output: <outdir>/<system>_rep<N>_analysis.json
"""
import warnings
warnings.filterwarnings('ignore')

import argparse
import numpy as np
import pandas as pd
import json
import sys
import io
from pathlib import Path

from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk
from alchemlyb.estimators import TI, BAR, MBAR
from alchemlyb.preprocessing import decorrelate_dhdl, decorrelate_u_nk

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

KJ_TO_KCAL = 0.239006
N_LAMBDA = 32
T = 300

def analyze_leg(base, leg_name, rep):
    xvg_files = [str(base / f'lambda_{i:02d}' / 'prod.xvg') for i in range(N_LAMBDA)]

    missing = [f for f in xvg_files if not Path(f).exists()]
    if missing:
        print(f"ERROR: {len(missing)} missing xvg files")
        for m in missing[:5]:
            print(f"  {m}")
        return None

    print(f"{'='*60}")
    print(f"ANALYZING {leg_name.upper()} REP{rep}")
    print(f"{'='*60}")

    dhdl_raw_list, dhdl_dec_list = [], []
    u_nk_raw_list, u_nk_dec_list = [], []

    for i, f in enumerate(xvg_files):
        d = extract_dHdl(f, T=T)
        dhdl_raw_list.append(d)
        dhdl_dec_list.append(decorrelate_dhdl(d))
        try:
            u = extract_u_nk(f, T=T)
            u_nk_raw_list.append(u)
            u_nk_dec_list.append(decorrelate_u_nk(u))
        except Exception:
            pass
        if i % 8 == 7 or i == N_LAMBDA - 1:
            print(f"  Parsed l_00..l_{i:02d}")

    combined_raw = pd.concat(dhdl_raw_list)
    combined_dec = pd.concat(dhdl_dec_list)
    print(f"  dH/dl: {len(combined_raw)} raw -> {len(combined_dec)} decorrelated")

    ti_raw = TI()
    ti_raw.fit(combined_raw)
    dg_ti_raw = ti_raw.delta_f_.iloc[0, -1]
    err_ti_raw = ti_raw.d_delta_f_.iloc[0, -1]

    ti_dec = TI()
    ti_dec.fit(combined_dec)
    dg_ti_dec = ti_dec.delta_f_.iloc[0, -1]
    err_ti_dec = ti_dec.d_delta_f_.iloc[0, -1]

    dg_bar, err_bar, dg_mbar, err_mbar = None, None, None, None
    n_dec_u = 0

    if u_nk_dec_list:
        combined_u_dec = pd.concat(u_nk_dec_list)
        n_dec_u = len(combined_u_dec)
        print(f"  u_nk: {len(pd.concat(u_nk_raw_list))} raw -> {n_dec_u} decorrelated")

        try:
            bar = BAR()
            bar.fit(combined_u_dec)
            dg_bar = bar.delta_f_.iloc[0, -1]
            err_bar = bar.d_delta_f_.iloc[0, -1]
        except Exception as e:
            print(f"  BAR FAILED: {e}")

        try:
            mbar = MBAR()
            mbar.fit(combined_u_dec)
            dg_mbar = mbar.delta_f_.iloc[0, -1]
            err_mbar = mbar.d_delta_f_.iloc[0, -1]
        except Exception as e:
            print(f"  MBAR FAILED: {e}")

    print(f"\n  {leg_name} rep{rep} RESULTS:")
    print(f"  TI (raw):          {dg_ti_raw:.3f} +/- {err_ti_raw:.3f} kJ/mol  ({dg_ti_raw*KJ_TO_KCAL:.3f} +/- {err_ti_raw*KJ_TO_KCAL:.3f} kcal/mol)")
    print(f"  TI (decorrelated): {dg_ti_dec:.3f} +/- {err_ti_dec:.3f} kJ/mol  ({dg_ti_dec*KJ_TO_KCAL:.3f} +/- {err_ti_dec*KJ_TO_KCAL:.3f} kcal/mol)")
    if dg_bar is not None:
        print(f"  BAR:               {dg_bar:.3f} +/- {err_bar:.3f} kJ/mol  ({dg_bar*KJ_TO_KCAL:.3f} +/- {err_bar*KJ_TO_KCAL:.3f} kcal/mol)")
    if dg_mbar is not None:
        print(f"  MBAR:              {dg_mbar:.3f} +/- {err_mbar:.3f} kJ/mol  ({dg_mbar*KJ_TO_KCAL:.3f} +/- {err_mbar*KJ_TO_KCAL:.3f} kcal/mol)")

    results = {
        "leg": leg_name,
        "replica": rep,
        "n_lambda": N_LAMBDA,
        "temperature_K": T,
        "TI_raw": {"dG_kj": round(float(dg_ti_raw), 4), "err_kj": round(float(err_ti_raw), 4),
                    "dG_kcal": round(float(dg_ti_raw*KJ_TO_KCAL), 4), "err_kcal": round(float(err_ti_raw*KJ_TO_KCAL), 4),
                    "n_samples": int(len(combined_raw))},
        "TI_decorrelated": {"dG_kj": round(float(dg_ti_dec), 4), "err_kj": round(float(err_ti_dec), 4),
                            "dG_kcal": round(float(dg_ti_dec*KJ_TO_KCAL), 4), "err_kcal": round(float(err_ti_dec*KJ_TO_KCAL), 4),
                            "n_samples": int(len(combined_dec))}
    }
    if dg_bar is not None:
        results["BAR_decorrelated"] = {"dG_kj": round(float(dg_bar), 4), "err_kj": round(float(err_bar), 4),
                                        "dG_kcal": round(float(dg_bar*KJ_TO_KCAL), 4), "err_kcal": round(float(err_bar*KJ_TO_KCAL), 4),
                                        "n_samples": int(n_dec_u)}
    if dg_mbar is not None:
        results["MBAR_decorrelated"] = {"dG_kj": round(float(dg_mbar), 4), "err_kj": round(float(err_mbar), 4),
                                         "dG_kcal": round(float(dg_mbar*KJ_TO_KCAL), 4), "err_kcal": round(float(err_mbar*KJ_TO_KCAL), 4),
                                         "n_samples": int(n_dec_u)}

    n_states = ti_raw.delta_f_.shape[1]
    results["per_window_dG_kj"] = [round(float(ti_raw.delta_f_.iloc[0, j]), 4) for j in range(n_states)]

    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="FEP per-leg/rep analysis via alchemlyb (TI, BAR, MBAR)")
    parser.add_argument("system", choices=["ternary", "binary", "solvent"],
                        help="FEP system to analyze")
    parser.add_argument("rep", type=int, choices=[1, 2, 3],
                        help="Replicate number")
    parser.add_argument("--base", required=True,
                        help="Base directory containing <system>/rep<N>/lambda_*/prod.xvg")
    parser.add_argument("--outdir", default=".",
                        help="Output directory for JSON (default: current directory)")
    args = parser.parse_args()

    base = Path(args.base) / args.system / f"rep{args.rep}"
    if not base.exists():
        print(f"ERROR: {base} does not exist")
        sys.exit(1)

    results = analyze_leg(base, args.system, args.rep)
    if results:
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        out_path = outdir / f"{args.system}_rep{args.rep}_analysis.json"
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved to {out_path}")
