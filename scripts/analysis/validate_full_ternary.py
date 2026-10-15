import os

# DATA_ROOT: path to supplementary data package root (hiv1-gp120-cd4-ccr5/)
# Adjust this if running from a different location
DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)
#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Validate 20ns full ternary complex trajectory using GROMACS tools.
Parses XVG outputs into JSON for pipeline integration.
"""

import os
import json
import subprocess
from datetime import datetime
from pathlib import Path

WORK_DIR = Path(os.path.join(DATA_ROOT, "topology"))
OUTDIR = Path(os.path.join(DATA_ROOT, "results", "00-summary"))
OUTDIR.mkdir(exist_ok=True)

def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def run_gmx(cmd, input_str=None):
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True,
        input=input_str, timeout=600, cwd=str(WORK_DIR)
    )
    return result.stdout + result.stderr

def parse_xvg(filepath):
    """Parse a GROMACS XVG file into arrays."""
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

results = {"timestamp": datetime.now().isoformat(), "simulation": "20ns full ternary complex (gp120+CD4+CCR5)"}

log("=" * 70)
log("TRAJECTORY VALIDATION: 20ns Full Ternary Complex")
log("=" * 70)

# 1. Trajectory check
log("Step 1: Trajectory integrity")
output = run_gmx("gmx check -f md.xtc 2>&1")
for line in output.split("\n"):
    if "Last frame" in line:
        log(f"  {line.strip()}")
        parts = line.strip().split()
        for i, p in enumerate(parts):
            if p == "time":
                results["last_time_ps"] = float(parts[i+1])
                results["total_time_ns"] = float(parts[i+1]) / 1000

# 2. System info
log("Step 2: System info")
output = run_gmx("gmx check -s md.tpr 2>&1")
for line in output.split("\n"):
    if "natoms" in line.lower() or "nsteps" in line.lower():
        log(f"  {line.strip()}")

# 3. Temperature
log("Step 3: Temperature")
output = run_gmx("echo '16 0' | gmx energy -f md.edr -o /tmp/temp.xvg 2>&1")
for line in output.split("\n"):
    if "Temperature" in line and "Average" not in line:
        parts = line.split()
        try:
            results["temperature"] = {"mean": float(parts[1]), "rmsd": float(parts[3])}
            log(f"  T = {parts[1]} +/- {parts[3]} K")
        except (IndexError, ValueError):
            log(f"  {line.strip()}")

# 4. Potential energy
log("Step 4: Potential energy")
output = run_gmx("echo '11 0' | gmx energy -f md.edr -o /tmp/pot.xvg 2>&1")
for line in output.split("\n"):
    if "Potential" in line:
        parts = line.split()
        try:
            results["potential_energy"] = {"mean_kJmol": float(parts[1])}
            log(f"  V = {parts[1]} kJ/mol")
        except (IndexError, ValueError):
            log(f"  {line.strip()}")

# 5. Backbone RMSD
log("Step 5: Backbone RMSD")
rmsd_xvg = str(OUTDIR / "rmsd_backbone.xvg")
output = run_gmx(f"echo '4 4' | gmx rms -s md.tpr -f md.xtc -o {rmsd_xvg} -tu ns 2>&1")
if os.path.exists(rmsd_xvg):
    rmsd_data = parse_xvg(rmsd_xvg)
    if rmsd_data:
        import numpy as np
        times = [d[0] for d in rmsd_data]
        rmsds = [d[1] for d in rmsd_data]
        results["backbone_rmsd"] = {
            "mean_nm": float(np.mean(rmsds)),
            "mean_A": float(np.mean(rmsds) * 10),
            "std_A": float(np.std(rmsds) * 10),
            "max_A": float(np.max(rmsds) * 10),
            "final_A": float(rmsds[-1] * 10),
            "n_points": len(rmsds),
            "plateau_reached": float(np.std(rmsds[-len(rmsds)//4:]) * 10) < 0.3
        }
        log(f"  RMSD: mean={np.mean(rmsds)*10:.2f} +/- {np.std(rmsds)*10:.2f} A, max={np.max(rmsds)*10:.2f} A")
        # Last quarter stability
        last_q = rmsds[-len(rmsds)//4:]
        log(f"  Last quarter: mean={np.mean(last_q)*10:.2f} +/- {np.std(last_q)*10:.2f} A (plateau={'YES' if results['backbone_rmsd']['plateau_reached'] else 'NO'})")

# 6. Per-residue RMSF
log("Step 6: Per-residue RMSF")
rmsf_xvg = str(OUTDIR / "rmsf_backbone.xvg")
output = run_gmx(f"echo '4' | gmx rmsf -s md.tpr -f md.xtc -o {rmsf_xvg} -res 2>&1")
if os.path.exists(rmsf_xvg):
    rmsf_data = parse_xvg(rmsf_xvg)
    if rmsf_data:
        import numpy as np
        resnums = [int(d[0]) for d in rmsf_data]
        rmsf_vals = [d[1] for d in rmsf_data]
        
        results["rmsf"] = {
            "mean_all_nm": float(np.mean(rmsf_vals)),
            "mean_all_A": float(np.mean(rmsf_vals) * 10),
        }
        
        # Check for relay pocket residues across possible numberings
        for label, pocket_range in [("6MEO", range(465, 473)), ("SOSIP", range(469, 477)), ("reduced", range(69, 77))]:
            pocket_rmsf = [(r, v) for r, v in zip(resnums, rmsf_vals) if r in pocket_range]
            if pocket_rmsf:
                pocket_vals = [v for _, v in pocket_rmsf]
                results["rmsf"][f"pocket_{label}"] = {
                    "residues": [{"resnum": r, "rmsf_nm": float(v), "rmsf_A": float(v*10)} for r, v in pocket_rmsf],
                    "mean_A": float(np.mean(pocket_vals) * 10),
                }
                log(f"  Pocket ({label}): mean RMSF = {np.mean(pocket_vals)*10:.2f} A")
        
        # Top 10 most flexible
        sorted_rmsf = sorted(zip(resnums, rmsf_vals), key=lambda x: x[1], reverse=True)
        results["rmsf"]["top10_flexible"] = [{"resnum": r, "rmsf_A": float(v*10)} for r, v in sorted_rmsf[:10]]
        log(f"  Top flexible: res {sorted_rmsf[0][0]} ({sorted_rmsf[0][1]*10:.2f} A)")

# 7. Radius of gyration
log("Step 7: Radius of gyration")
rg_xvg = str(OUTDIR / "gyrate.xvg")
output = run_gmx(f"echo '1' | gmx gyrate -s md.tpr -f md.xtc -o {rg_xvg} 2>&1")
if os.path.exists(rg_xvg):
    rg_data = parse_xvg(rg_xvg)
    if rg_data:
        import numpy as np
        rg_vals = [d[1] for d in rg_data]  # column 1 = total Rg
        results["radius_of_gyration"] = {
            "mean_nm": float(np.mean(rg_vals)),
            "mean_A": float(np.mean(rg_vals) * 10),
            "std_A": float(np.std(rg_vals) * 10),
        }
        log(f"  Rg: mean={np.mean(rg_vals)*10:.2f} +/- {np.std(rg_vals)*10:.2f} A")

# 8. Hydrogen bonds
log("Step 8: Hydrogen bonds")
hb_xvg = str(OUTDIR / "hbonds.xvg")
output = run_gmx(f"echo '1 1' | gmx hbond -s md.tpr -f md.xtc -num {hb_xvg} 2>&1")
if os.path.exists(hb_xvg):
    hb_data = parse_xvg(hb_xvg)
    if hb_data:
        import numpy as np
        hb_vals = [d[1] for d in hb_data]
        results["hbonds"] = {
            "mean": float(np.mean(hb_vals)),
            "std": float(np.std(hb_vals)),
        }
        log(f"  H-bonds: mean={np.mean(hb_vals):.0f} +/- {np.std(hb_vals):.0f}")

# Save
with open(OUTDIR / "validation_results.json", "w") as f:
    json.dump(results, f, indent=2)

log("")
log("=" * 70)
log("VALIDATION COMPLETE")
log(f"Results: {OUTDIR / 'validation_results.json'}")
log("=" * 70)
