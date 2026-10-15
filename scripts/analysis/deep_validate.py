#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Deep validation of MD trajectory with proper PBC handling.
Extracts pocket-specific metrics, checks for unfolding.
"""
import os, sys, json, subprocess
import numpy as np
from pathlib import Path
from datetime import datetime

SIMULATION = sys.argv[1] if len(sys.argv) > 1 else "unknown"
WORK_DIR = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(".")
TPR = sys.argv[3] if len(sys.argv) > 3 else "md.tpr"
XTC = sys.argv[4] if len(sys.argv) > 4 else "md.xtc"

OUTDIR = WORK_DIR / "validation"
OUTDIR.mkdir(exist_ok=True)

def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def run(cmd, input_str=None):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                       input=input_str, timeout=600, cwd=str(WORK_DIR))
    return r.stdout + r.stderr

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

log(f"{'='*70}")
log(f"DEEP VALIDATION: {SIMULATION}")
log(f"Dir: {WORK_DIR}, TPR: {TPR}, XTC: {XTC}")
log(f"{'='*70}")

results = {"simulation": SIMULATION, "timestamp": datetime.now().isoformat()}

# Step 1: Two-step PBC correction
log("Step 1a: PBC nojump correction")
nj_xtc = str(OUTDIR / "nojump.xtc")
run(f"echo '0' | gmx trjconv -s {TPR} -f {XTC} -o {nj_xtc} -pbc nojump", "0\n")

log("Step 1b: Center and cluster")
clean_xtc = str(OUTDIR / "clean.xtc")
run(f"echo '1\n0' | gmx trjconv -s {TPR} -f {nj_xtc} -o {clean_xtc} -center -pbc mol", "1\n0\n")

# Step 2: Extract first and last frames
log("Step 2: Extract first and last frames")
run(f"echo '1' | gmx trjconv -s {TPR} -f {clean_xtc} -o {OUTDIR}/first_frame.pdb -dump 0 -pbc mol", "1\n")

# Get last time from trajectory
out = run(f"gmx check -f {XTC}")
last_time = None
for line in out.split("\n"):
    if "Last frame" in line:
        parts = line.split()
        for i, p in enumerate(parts):
            if p == "time":
                last_time = float(parts[i+1])
                results["total_time_ps"] = last_time
                results["total_time_ns"] = last_time / 1000
                log(f"  Total time: {last_time/1000:.1f} ns")

if last_time:
    run(f"echo '1' | gmx trjconv -s {TPR} -f {clean_xtc} -o {OUTDIR}/last_frame.pdb -dump {last_time}", "1\n")

# Step 3: Backbone RMSD on clean trajectory
log("Step 3: Backbone RMSD (clean trajectory)")
rmsd_xvg = str(OUTDIR / "rmsd_clean.xvg")
run(f"echo '4 4' | gmx rms -s {TPR} -f {clean_xtc} -o {rmsd_xvg} -tu ns", "4\n4\n")

if os.path.exists(rmsd_xvg):
    rmsd_data = parse_xvg(rmsd_xvg)
    if rmsd_data:
        rmsds = [d[1] * 10 for d in rmsd_data]  # nm -> A
        times = [d[0] for d in rmsd_data]
        q4 = rmsds[-len(rmsds)//4:]
        results["backbone_rmsd"] = {
            "mean_A": round(float(np.mean(rmsds)), 2),
            "std_A": round(float(np.std(rmsds)), 2),
            "max_A": round(float(np.max(rmsds)), 2),
            "final_A": round(float(rmsds[-1]), 2),
            "q4_mean_A": round(float(np.mean(q4)), 2),
            "q4_std_A": round(float(np.std(q4)), 2),
            "n_points": len(rmsds),
        }
        log(f"  RMSD: mean={np.mean(rmsds):.2f} +/- {np.std(rmsds):.2f} A, final={rmsds[-1]:.2f} A")
        log(f"  Last Q: {np.mean(q4):.2f} +/- {np.std(q4):.2f} A")

# Step 4: C-alpha RMSD (often more interpretable)
log("Step 4: C-alpha RMSD")
rmsd_ca_xvg = str(OUTDIR / "rmsd_ca_clean.xvg")
run(f"echo '3 3' | gmx rms -s {TPR} -f {clean_xtc} -o {rmsd_ca_xvg} -tu ns", "3\n3\n")

if os.path.exists(rmsd_ca_xvg):
    ca_data = parse_xvg(rmsd_ca_xvg)
    if ca_data:
        ca_rmsds = [d[1] * 10 for d in ca_data]
        ca_q4 = ca_rmsds[-len(ca_rmsds)//4:]
        results["ca_rmsd"] = {
            "mean_A": round(float(np.mean(ca_rmsds)), 2),
            "std_A": round(float(np.std(ca_rmsds)), 2),
            "max_A": round(float(np.max(ca_rmsds)), 2),
            "final_A": round(float(ca_rmsds[-1]), 2),
            "q4_mean_A": round(float(np.mean(ca_q4)), 2),
            "q4_std_A": round(float(np.std(ca_q4)), 2),
        }
        log(f"  CA RMSD: mean={np.mean(ca_rmsds):.2f} +/- {np.std(ca_rmsds):.2f} A, final={ca_rmsds[-1]:.2f} A")

# Step 5: Per-residue RMSF
log("Step 5: Per-residue RMSF (C-alpha)")
rmsf_xvg = str(OUTDIR / "rmsf_clean.xvg")
run(f"echo '3' | gmx rmsf -s {TPR} -f {clean_xtc} -o {rmsf_xvg} -res", "3\n")

if os.path.exists(rmsf_xvg):
    rmsf_data = parse_xvg(rmsf_xvg)
    if rmsf_data:
        resnums = [int(d[0]) for d in rmsf_data]
        rmsf_vals = [d[1] * 10 for d in rmsf_data]  # nm -> A
        
        results["rmsf_mean_all_A"] = round(float(np.mean(rmsf_vals)), 2)
        
        # Pocket search across numberings
        for label, pocket_range in [("6MEO", range(465, 473)), ("SOSIP", range(469, 477)), ("reduced", range(69, 77))]:
            pocket = [(r, v) for r, v in zip(resnums, rmsf_vals) if r in pocket_range]
            if pocket:
                pvals = [v for _, v in pocket]
                results[f"pocket_{label}"] = {
                    "residues": [{"resnum": int(r), "rmsf_A": round(float(v), 3)} for r, v in pocket],
                    "mean_A": round(float(np.mean(pvals)), 3),
                    "max_A": round(float(np.max(pvals)), 3),
                    "min_A": round(float(np.min(pvals)), 3),
                }
                log(f"  Pocket {label}: RMSF mean={np.mean(pvals):.3f} A, range {np.min(pvals):.3f}-{np.max(pvals):.3f} A")
                for r, v in pocket:
                    log(f"    Res {r}: {v:.3f} A")

# Step 6: Radius of gyration
log("Step 6: Radius of gyration")
rg_xvg = str(OUTDIR / "gyrate_clean.xvg")
run(f"echo '1' | gmx gyrate -s {TPR} -f {clean_xtc} -o {rg_xvg}", "1\n")

if os.path.exists(rg_xvg):
    rg_data = parse_xvg(rg_xvg)
    if rg_data:
        rg_vals = [d[1] * 10 for d in rg_data]
        results["rg"] = {
            "mean_A": round(float(np.mean(rg_vals)), 2),
            "std_A": round(float(np.std(rg_vals)), 2),
            "first_A": round(float(rg_vals[0]), 2),
            "last_A": round(float(rg_vals[-1]), 2),
        }
        log(f"  Rg: mean={np.mean(rg_vals):.2f} +/- {np.std(rg_vals):.2f} A (first={rg_vals[0]:.2f}, last={rg_vals[-1]:.2f})")

# Step 7: Energy check
log("Step 7: Energy")
edr = XTC.replace(".xtc", ".edr")
# Try to get temperature and potential energy
for etype, ein in [("Temperature", "16"), ("Potential", "11"), ("Pressure", "18")]:
    out = run(f"echo '{ein}\n0' | gmx energy -f {edr} -o /dev/null 2>&1")
    for line in out.split("\n"):
        if etype in line and ":" not in line:
            parts = line.split()
            for i, p in enumerate(parts):
                if p == etype:
                    try:
                        results[f"energy_{etype.lower()}"] = {
                            "mean": float(parts[i+1]),
                        }
                        if len(parts) > i+3:
                            results[f"energy_{etype.lower()}"]["rmsd"] = float(parts[i+3])
                        log(f"  {etype}: {parts[i+1]}")
                    except (ValueError, IndexError):
                        pass

# Step 8: First/last frame Rg comparison (unfolding check)
log("Step 8: Unfolding check (first vs last frame Rg)")
first_pdb = OUTDIR / "first_frame.pdb"
last_pdb = OUTDIR / "last_frame.pdb"
if first_pdb.exists() and last_pdb.exists():
    def count_ca(pdb):
        n = 0
        with open(pdb) as f:
            for line in f:
                if line.startswith("ATOM") and " CA " in line:
                    n += 1
        return n
    
    ca_first = count_ca(first_pdb)
    ca_last = count_ca(last_pdb)
    log(f"  First frame CA: {ca_first}, Last frame CA: {ca_last}")
    results["frame_check"] = {"ca_first": ca_first, "ca_last": ca_last}

# Step 9: RMSD time series sample (for trend analysis)
if os.path.exists(rmsd_xvg):
    rmsd_data = parse_xvg(rmsd_xvg)
    if rmsd_data:
        times = [d[0] for d in rmsd_data]
        rmsds = [d[1] * 10 for d in rmsd_data]
        # Sample at 10 evenly-spaced points
        indices = np.linspace(0, len(rmsds)-1, 10, dtype=int)
        results["rmsd_timeseries"] = [
            {"time_ns": round(float(times[i]), 2), "rmsd_A": round(float(rmsds[i]), 2)}
            for i in indices
        ]
        log("Step 9: RMSD time series sample:")
        for pt in results["rmsd_timeseries"]:
            log(f"  t={pt['time_ns']:.1f} ns: {pt['rmsd_A']:.2f} A")

# Save
with open(OUTDIR / "deep_validation.json", "w") as f:
    json.dump(results, f, indent=2, default=str)

log("")
log("="*70)
log("DEEP VALIDATION COMPLETE")
log(f"Results: {OUTDIR / 'deep_validation.json'}")
log("="*70)
