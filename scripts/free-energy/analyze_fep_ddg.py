#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Binary rep1 FEP analysis + preliminary ddG calculation.
Same methodology as ternary rep1: TI (raw + decorrelated), BAR, MBAR.
Then computes ddG = dG(ternary) - dG(binary).
"""
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import json
import sys
import os
import io
from pathlib import Path

from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk
from alchemlyb.estimators import TI, BAR, MBAR
from alchemlyb.preprocessing import decorrelate_dhdl, decorrelate_u_nk

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

KJ_TO_KCAL = 0.239006
N_LAMBDA = 32
T = 300

def analyze_leg(base, leg_name):
    xvg_files = [str(base / f'lambda_{i:02d}' / 'prod.xvg') for i in range(N_LAMBDA)]

    print("=" * 60)
    print(f"ANALYZING {leg_name.upper()}")
    print("=" * 60)

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

    # TI raw
    ti_raw = TI()
    ti_raw.fit(combined_raw)
    dg_ti_raw = ti_raw.delta_f_.iloc[0, -1]
    err_ti_raw = ti_raw.d_delta_f_.iloc[0, -1]

    # TI decorrelated
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

    print(f"\n  {leg_name} RESULTS:")
    print(f"  TI (raw):          {dg_ti_raw:.3f} +/- {err_ti_raw:.3f} kJ/mol  ({dg_ti_raw*KJ_TO_KCAL:.3f} +/- {err_ti_raw*KJ_TO_KCAL:.3f} kcal/mol)")
    print(f"  TI (decorrelated): {dg_ti_dec:.3f} +/- {err_ti_dec:.3f} kJ/mol  ({dg_ti_dec*KJ_TO_KCAL:.3f} +/- {err_ti_dec*KJ_TO_KCAL:.3f} kcal/mol)")
    if dg_bar is not None:
        print(f"  BAR:               {dg_bar:.3f} +/- {err_bar} kJ/mol  ({dg_bar*KJ_TO_KCAL:.3f} kcal/mol)")
    if dg_mbar is not None:
        print(f"  MBAR:              {dg_mbar:.3f} +/- {err_mbar:.3f} kJ/mol  ({dg_mbar*KJ_TO_KCAL:.3f} +/- {err_mbar*KJ_TO_KCAL:.3f} kcal/mol)")

    results = {
        "leg": leg_name,
        "replica": 1,
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
        results["BAR_decorrelated"] = {"dG_kj": round(float(dg_bar), 4), "err_kj": float(err_bar) if err_bar is not None else None,
                                        "dG_kcal": round(float(dg_bar*KJ_TO_KCAL), 4), "n_samples": int(n_dec_u)}
    if dg_mbar is not None:
        results["MBAR_decorrelated"] = {"dG_kj": round(float(dg_mbar), 4), "err_kj": round(float(err_mbar), 4),
                                         "dG_kcal": round(float(dg_mbar*KJ_TO_KCAL), 4), "err_kcal": round(float(err_mbar*KJ_TO_KCAL), 4),
                                         "n_samples": int(n_dec_u)}

    n_states = ti_raw.delta_f_.shape[1]
    results["per_window_dG_kj"] = [round(float(ti_raw.delta_f_.iloc[0, j]), 4) for j in range(n_states)]

    return results


def main():
    base_dir = Path(os.environ.get('FEP_BASE', '.'))
    ternary_base = base_dir / 'runs' / 'ternary' / 'rep1'
    binary_base = base_dir / 'runs' / 'binary' / 'rep1'

    ternary_json = base_dir / 'ternary_rep1_analysis.json'
    if ternary_json.exists():
        print("Loading cached ternary rep1 results...")
        with open(ternary_json) as f:
            ternary_results = json.load(f)
        dg_tern_kj = ternary_results["TI_decorrelated"]["dG_kj"]
        err_tern_kj = ternary_results["TI_decorrelated"]["err_kj"]
        dg_tern_mbar_kj = ternary_results.get("MBAR_decorrelated", {}).get("dG_kj")
        print(f"  Ternary TI(dec): {dg_tern_kj:.3f} +/- {err_tern_kj:.3f} kJ/mol")
        if dg_tern_mbar_kj:
            print(f"  Ternary MBAR:    {dg_tern_mbar_kj:.3f} kJ/mol")
    else:
        print("Analyzing ternary rep1...")
        ternary_results = analyze_leg(ternary_base, "ternary")
        dg_tern_kj = ternary_results["TI_decorrelated"]["dG_kj"]
        err_tern_kj = ternary_results["TI_decorrelated"]["err_kj"]
        dg_tern_mbar_kj = ternary_results.get("MBAR_decorrelated", {}).get("dG_kj")

    # Analyze binary
    binary_results = analyze_leg(binary_base, "binary")

    # Save binary results
    out_binary = base_dir / 'binary_rep1_analysis.json'
    with open(out_binary, 'w') as f:
        json.dump(binary_results, f, indent=2)
    print(f"\nBinary results saved to {out_binary}")

    # Compute ddG
    dg_bin_kj = binary_results["TI_decorrelated"]["dG_kj"]
    err_bin_kj = binary_results["TI_decorrelated"]["err_kj"]
    dg_bin_mbar_kj = binary_results.get("MBAR_decorrelated", {}).get("dG_kj")
    err_bin_mbar_kj = binary_results.get("MBAR_decorrelated", {}).get("err_kj")

    print("\n" + "=" * 70)
    print("THERMODYNAMIC CYCLE -- PRELIMINARY ddG (single rep)")
    print("=" * 70)
    print()
    print("ddG = dG(ternary) - dG(binary)")
    print("  If |ddG| < 1 kcal/mol --> NON-COMPETITIVE")
    print("  If |ddG| > 2 kcal/mol --> competitive coupling")
    print()

    # TI
    ddg_ti_kj = dg_tern_kj - dg_bin_kj
    ddg_ti_err_kj = np.sqrt(err_tern_kj**2 + err_bin_kj**2)
    ddg_ti_kcal = ddg_ti_kj * KJ_TO_KCAL
    ddg_ti_err_kcal = ddg_ti_err_kj * KJ_TO_KCAL

    print(f"  TI (decorrelated):")
    print(f"    dG(ternary) = {dg_tern_kj:.3f} +/- {err_tern_kj:.3f} kJ/mol  ({dg_tern_kj*KJ_TO_KCAL:.3f} kcal/mol)")
    print(f"    dG(binary)  = {dg_bin_kj:.3f} +/- {err_bin_kj:.3f} kJ/mol  ({dg_bin_kj*KJ_TO_KCAL:.3f} kcal/mol)")
    print(f"    ddG         = {ddg_ti_kj:.3f} +/- {ddg_ti_err_kj:.3f} kJ/mol  ({ddg_ti_kcal:.3f} +/- {ddg_ti_err_kcal:.3f} kcal/mol)")
    print()

    if abs(ddg_ti_kcal) < 1.0:
        print(f"    >>> |ddG| = {abs(ddg_ti_kcal):.3f} < 1.0 kcal/mol --> NON-COMPETITIVE CONFIRMED")
    elif abs(ddg_ti_kcal) < 2.0:
        print(f"    >>> |ddG| = {abs(ddg_ti_kcal):.3f} -- borderline (1-2 kcal/mol)")
    else:
        print(f"    >>> |ddG| = {abs(ddg_ti_kcal):.3f} > 2.0 kcal/mol --> significant coupling detected")

    # MBAR
    if dg_tern_mbar_kj and dg_bin_mbar_kj:
        ddg_mbar_kj = dg_tern_mbar_kj - dg_bin_mbar_kj
        ddg_mbar_kcal = ddg_mbar_kj * KJ_TO_KCAL
        err_tern_mbar = ternary_results.get("MBAR_decorrelated", {}).get("err_kj", 0)
        ddg_mbar_err_kj = np.sqrt(err_tern_mbar**2 + (err_bin_mbar_kj or 0)**2)
        ddg_mbar_err_kcal = ddg_mbar_err_kj * KJ_TO_KCAL

        print()
        print(f"  MBAR (decorrelated):")
        print(f"    dG(ternary) = {dg_tern_mbar_kj:.3f} kJ/mol  ({dg_tern_mbar_kj*KJ_TO_KCAL:.3f} kcal/mol)")
        print(f"    dG(binary)  = {dg_bin_mbar_kj:.3f} kJ/mol  ({dg_bin_mbar_kj*KJ_TO_KCAL:.3f} kcal/mol)")
        print(f"    ddG         = {ddg_mbar_kj:.3f} +/- {ddg_mbar_err_kj:.3f} kJ/mol  ({ddg_mbar_kcal:.3f} +/- {ddg_mbar_err_kcal:.3f} kcal/mol)")
        print()

        if abs(ddg_mbar_kcal) < 1.0:
            print(f"    >>> |ddG| = {abs(ddg_mbar_kcal):.3f} < 1.0 kcal/mol --> NON-COMPETITIVE CONFIRMED")
        elif abs(ddg_mbar_kcal) < 2.0:
            print(f"    >>> |ddG| = {abs(ddg_mbar_kcal):.3f} -- borderline (1-2 kcal/mol)")
        else:
            print(f"    >>> |ddG| = {abs(ddg_mbar_kcal):.3f} > 2.0 kcal/mol --> significant coupling detected")

    # Save combined results
    cycle = {
        "ternary_rep1": ternary_results,
        "binary_rep1": binary_results,
        "ddG_TI_kj": round(float(ddg_ti_kj), 4),
        "ddG_TI_err_kj": round(float(ddg_ti_err_kj), 4),
        "ddG_TI_kcal": round(float(ddg_ti_kcal), 4),
        "ddG_TI_err_kcal": round(float(ddg_ti_err_kcal), 4),
    }
    if dg_tern_mbar_kj and dg_bin_mbar_kj:
        cycle["ddG_MBAR_kj"] = round(float(ddg_mbar_kj), 4)
        cycle["ddG_MBAR_err_kj"] = round(float(ddg_mbar_err_kj), 4)
        cycle["ddG_MBAR_kcal"] = round(float(ddg_mbar_kcal), 4)
        cycle["ddG_MBAR_err_kcal"] = round(float(ddg_mbar_err_kcal), 4)

    cycle_path = base_dir / 'fep_ddg_preliminary.json'
    with open(cycle_path, 'w') as f:
        json.dump(cycle, f, indent=2)
    print(f"\nFull cycle results saved to {cycle_path}")


if __name__ == '__main__':
    main()
