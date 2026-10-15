#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Comprehensive FEP analysis for ternary rep1.
Runs TI (with per-window decorrelation), BAR, and MBAR on the complete 32-window dataset.
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
BASE = Path(os.environ.get('FEP_BASE', './runs/ternary/rep1'))
N_LAMBDA = 32
T = 300  # K

def run_analysis():
    xvg_files = [str(BASE / f'lambda_{i:02d}' / 'prod.xvg') for i in range(N_LAMBDA)]

    print("=" * 60)
    print("PARSING DATA")
    print("=" * 60)
    dhdl_raw_list = []
    dhdl_dec_list = []
    u_nk_raw_list = []
    u_nk_dec_list = []

    for i, f in enumerate(xvg_files):
        # dH/dl
        d = extract_dHdl(f, T=T)
        dhdl_raw_list.append(d)

        d_dec = decorrelate_dhdl(d)
        dhdl_dec_list.append(d_dec)

        # u_nk
        try:
            u = extract_u_nk(f, T=T)
            u_nk_raw_list.append(u)
            u_dec = decorrelate_u_nk(u)
            u_nk_dec_list.append(u_dec)
            print(f"  l_{i:02d}: {len(d)} raw -> {len(d_dec)} dec dhdl | {len(u)} raw -> {len(u_dec)} dec u_nk")
        except Exception as e:
            print(f"  l_{i:02d}: {len(d)} raw -> {len(d_dec)} dec dhdl | u_nk FAILED: {e}")

    combined_dhdl_raw = pd.concat(dhdl_raw_list)
    combined_dhdl_dec = pd.concat(dhdl_dec_list)
    print(f"\nTotal dH/dl: {len(combined_dhdl_raw)} raw -> {len(combined_dhdl_dec)} decorrelated")

    # -- TI (raw) --
    print("\n" + "=" * 60)
    print("TI (raw, no decorrelation)")
    print("=" * 60)
    ti_raw = TI()
    ti_raw.fit(combined_dhdl_raw)
    dg_ti_raw = ti_raw.delta_f_.iloc[0, -1]
    err_ti_raw = ti_raw.d_delta_f_.iloc[0, -1]
    print(f"  dG = {dg_ti_raw:.3f} +/- {err_ti_raw:.3f} kJ/mol")
    print(f"  dG = {dg_ti_raw*KJ_TO_KCAL:.3f} +/- {err_ti_raw*KJ_TO_KCAL:.3f} kcal/mol")

    # -- TI (decorrelated per-window) --
    print("\n" + "=" * 60)
    print("TI (per-window decorrelated)")
    print("=" * 60)
    ti_dec = TI()
    ti_dec.fit(combined_dhdl_dec)
    dg_ti_dec = ti_dec.delta_f_.iloc[0, -1]
    err_ti_dec = ti_dec.d_delta_f_.iloc[0, -1]
    print(f"  dG = {dg_ti_dec:.3f} +/- {err_ti_dec:.3f} kJ/mol")
    print(f"  dG = {dg_ti_dec*KJ_TO_KCAL:.3f} +/- {err_ti_dec*KJ_TO_KCAL:.3f} kcal/mol")

    dg_bar, err_bar = None, None
    dg_mbar, err_mbar = None, None
    n_dec_u = 0

    if u_nk_dec_list:
        combined_u_nk_dec = pd.concat(u_nk_dec_list)
        n_dec_u = len(combined_u_nk_dec)
        print(f"\nTotal u_nk: {len(pd.concat(u_nk_raw_list))} raw -> {n_dec_u} decorrelated")

        # -- BAR --
        print("\n" + "=" * 60)
        print("BAR (per-window decorrelated)")
        print("=" * 60)
        try:
            bar = BAR()
            bar.fit(combined_u_nk_dec)
            dg_bar = bar.delta_f_.iloc[0, -1]
            err_bar = bar.d_delta_f_.iloc[0, -1]
            print(f"  dG = {dg_bar:.3f} +/- {err_bar:.3f} kJ/mol")
            print(f"  dG = {dg_bar*KJ_TO_KCAL:.3f} +/- {err_bar*KJ_TO_KCAL:.3f} kcal/mol")
        except Exception as e:
            print(f"  BAR FAILED: {e}")

        # -- MBAR --
        print("\n" + "=" * 60)
        print("MBAR (per-window decorrelated)")
        print("=" * 60)
        try:
            mbar = MBAR()
            mbar.fit(combined_u_nk_dec)
            dg_mbar = mbar.delta_f_.iloc[0, -1]
            err_mbar = mbar.d_delta_f_.iloc[0, -1]
            print(f"  dG = {dg_mbar:.3f} +/- {err_mbar:.3f} kJ/mol")
            print(f"  dG = {dg_mbar*KJ_TO_KCAL:.3f} +/- {err_mbar*KJ_TO_KCAL:.3f} kcal/mol")
        except Exception as e:
            print(f"  MBAR FAILED: {e}")

    # -- Summary --
    print("\n" + "=" * 60)
    print("TERNARY REP1 -- CROSS-METHOD COMPARISON")
    print("=" * 60)
    print(f"{'Method':<30} {'dG (kJ/mol)':>20} {'dG (kcal/mol)':>20}")
    print("-" * 70)
    print(f"{'TI (raw)':<30} {dg_ti_raw:8.3f} +/- {err_ti_raw:.3f}    {dg_ti_raw*KJ_TO_KCAL:8.3f} +/- {err_ti_raw*KJ_TO_KCAL:.3f}")
    print(f"{'TI (decorrelated)':<30} {dg_ti_dec:8.3f} +/- {err_ti_dec:.3f}    {dg_ti_dec*KJ_TO_KCAL:8.3f} +/- {err_ti_dec*KJ_TO_KCAL:.3f}")
    if dg_bar is not None:
        print(f"{'BAR (decorrelated)':<30} {dg_bar:8.3f} +/- {err_bar:.3f}    {dg_bar*KJ_TO_KCAL:8.3f} +/- {err_bar*KJ_TO_KCAL:.3f}")
    if dg_mbar is not None:
        print(f"{'MBAR (decorrelated)':<30} {dg_mbar:8.3f} +/- {err_mbar:.3f}    {dg_mbar*KJ_TO_KCAL:8.3f} +/- {err_mbar*KJ_TO_KCAL:.3f}")

    print(f"\n  Binding dG (TI, dec): {-dg_ti_dec*KJ_TO_KCAL:.3f} kcal/mol")
    if dg_mbar is not None:
        print(f"  Binding dG (MBAR):    {-dg_mbar*KJ_TO_KCAL:.3f} kcal/mol")

    # -- Per-window dG profile (from TI raw which works correctly) --
    print("\n" + "=" * 60)
    print("PER-WINDOW CUMULATIVE dG (TI raw)")
    print("=" * 60)
    n_states = ti_raw.delta_f_.shape[1]
    for j in range(n_states):
        val_kj = ti_raw.delta_f_.iloc[0, j]
        err_kj = ti_raw.d_delta_f_.iloc[0, j]
        print(f"  l_{j:02d}: {val_kj:8.3f} +/- {err_kj:5.3f} kJ/mol  ({val_kj*KJ_TO_KCAL:8.3f} +/- {err_kj*KJ_TO_KCAL:5.3f} kcal/mol)")

    # -- Overlap matrix for BAR (quality check) --
    if dg_bar is not None:
        print("\n" + "=" * 60)
        print("BAR PER-PAIR dG")
        print("=" * 60)
        n_bar = bar.delta_f_.shape[0]
        for j in range(n_bar - 1):
            pair_dg = bar.delta_f_.iloc[j, j+1]
            pair_err = bar.d_delta_f_.iloc[j, j+1]
            print(f"  {j:02d}->{j+1:02d}: {pair_dg:8.3f} +/- {pair_err:5.3f} kJ/mol")

    # -- Save JSON --
    results = {
        "leg": "ternary",
        "replica": 1,
        "n_lambda": N_LAMBDA,
        "temperature_K": T,
        "TI_raw": {
            "dG_kj": round(float(dg_ti_raw), 4),
            "err_kj": round(float(err_ti_raw), 4),
            "dG_kcal": round(float(dg_ti_raw * KJ_TO_KCAL), 4),
            "err_kcal": round(float(err_ti_raw * KJ_TO_KCAL), 4),
            "n_samples": int(len(combined_dhdl_raw))
        },
        "TI_decorrelated": {
            "dG_kj": round(float(dg_ti_dec), 4),
            "err_kj": round(float(err_ti_dec), 4),
            "dG_kcal": round(float(dg_ti_dec * KJ_TO_KCAL), 4),
            "err_kcal": round(float(err_ti_dec * KJ_TO_KCAL), 4),
            "n_samples": int(len(combined_dhdl_dec))
        }
    }
    if dg_bar is not None:
        results["BAR_decorrelated"] = {
            "dG_kj": round(float(dg_bar), 4),
            "err_kj": round(float(err_bar), 4),
            "dG_kcal": round(float(dg_bar * KJ_TO_KCAL), 4),
            "err_kcal": round(float(err_bar * KJ_TO_KCAL), 4),
            "n_samples": int(n_dec_u)
        }
    if dg_mbar is not None:
        results["MBAR_decorrelated"] = {
            "dG_kj": round(float(dg_mbar), 4),
            "err_kj": round(float(err_mbar), 4),
            "dG_kcal": round(float(dg_mbar * KJ_TO_KCAL), 4),
            "err_kcal": round(float(err_mbar * KJ_TO_KCAL), 4),
            "n_samples": int(n_dec_u)
        }

    results["per_window_dG_kj"] = [round(float(ti_raw.delta_f_.iloc[0, j]), 4) for j in range(n_states)]
    results["per_window_err_kj"] = [round(float(ti_raw.d_delta_f_.iloc[0, j]), 4) for j in range(n_states)]

    out_path = Path('.') / 'ternary_rep1_analysis.json'
    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    run_analysis()
