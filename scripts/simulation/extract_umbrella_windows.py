#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Extract umbrella sampling windows from completed SMD pull trajectory.

Reads pull_pullx.xvg to find timestamps at evenly-spaced COM distance intervals,
then uses `gmx trjconv` to extract snapshot .gro files for each window.

The pull displacement was ~4.8 nm (48 A). We use 0.2 nm spacing, yielding ~24 windows.
Each window will be a starting configuration for a 10 ns umbrella simulation
restrained at that COM distance with k = 1000 kJ/mol/nm^2.

Usage:
    python3 extract_umbrella_windows.py
"""
import subprocess
import os
import sys

DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir)

WORK = os.path.join(DATA_ROOT, "results", "05-smd-v3")
PULL_XTC = os.path.join(WORK, "pull.xtc")
PULL_TPR = os.path.join(WORK, "pull.tpr")
PULLX_XVG = os.path.join(WORK, "pull_pullx.xvg")
INDEX = os.path.join(WORK, "index.ndx")
TOPOL = os.path.join(WORK, "topol.top")
UMBRELLA_MDP = os.path.join(WORK, "umbrella.mdp")
WINDOW_DIR = os.path.join(WORK, "umbrella_windows")

SPACING_NM = 0.2
PULL_RATE = 0.001  # nm/ps from MDP


def read_pullx(path):
    """Read pull coordinate vs time from XVG file."""
    data = []
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                t_ps = float(parts[0])
                pos_nm = float(parts[1])
                data.append((t_ps, pos_nm))
    return data


def find_window_times(data, spacing):
    """Find timestamps closest to evenly-spaced COM distance targets."""
    positions = [d[1] for d in data]
    p_min = min(positions)
    p_max = max(positions)
    
    targets = []
    p = p_min
    while p <= p_max:
        targets.append(p)
        p += spacing
    
    windows = []
    for target in targets:
        best_idx = min(range(len(data)), key=lambda i: abs(data[i][1] - target))
        t_ps, pos = data[best_idx]
        windows.append({
            "target_nm": target,
            "time_ps": t_ps,
            "actual_nm": pos,
            "error_nm": abs(pos - target),
        })
    
    return windows


def run(cmd, input_text=None):
    print(f"  >>> {cmd}")
    r = subprocess.run(cmd, shell=True, input=input_text,
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  STDERR: {r.stderr[-300:]}")
        print(f"  WARNING: exit code {r.returncode}")
    return r


def main():
    os.chdir(WORK)
    os.makedirs(WINDOW_DIR, exist_ok=True)
    
    print("=" * 60)
    print("Umbrella Window Extraction from SMD Pull")
    print("=" * 60)
    
    for f in [PULL_XTC, PULL_TPR, PULLX_XVG, INDEX, TOPOL, UMBRELLA_MDP]:
        if not os.path.exists(f):
            print(f"ERROR: Missing {f}")
            sys.exit(1)
    
    print(f"\nReading pull positions from {PULLX_XVG}...")
    data = read_pullx(PULLX_XVG)
    print(f"  {len(data)} data points")
    print(f"  Time range: {data[0][0]:.1f} - {data[-1][0]:.1f} ps")
    print(f"  Position range: {data[0][1]:.3f} - {data[-1][1]:.3f} nm")
    
    print(f"\nFinding windows at {SPACING_NM} nm spacing...")
    windows = find_window_times(data, SPACING_NM)
    print(f"  {len(windows)} windows identified")
    
    print(f"\n{'Window':>6} {'Target(nm)':>10} {'Actual(nm)':>10} {'Time(ps)':>10} {'Error(nm)':>10}")
    print("-" * 50)
    for i, w in enumerate(windows):
        print(f"{i:6d} {w['target_nm']:10.3f} {w['actual_nm']:10.3f} {w['time_ps']:10.1f} {w['error_nm']:10.4f}")
    
    max_err = max(w["error_nm"] for w in windows)
    print(f"\n  Max position error: {max_err:.4f} nm")
    if max_err > 0.05:
        print("  WARNING: Large position error -- check pull trajectory quality")
    
    print(f"\nExtracting {len(windows)} snapshots from pull.xtc...")
    for i, w in enumerate(windows):
        wdir = os.path.join(WINDOW_DIR, f"window_{i:02d}")
        os.makedirs(wdir, exist_ok=True)
        
        gro_out = os.path.join(wdir, "conf.gro")
        if os.path.exists(gro_out):
            print(f"  Window {i:02d}: conf.gro exists, skipping extraction")
            continue
        
        # echo 0 (System) to trjconv
        run(f"echo 0 | gmx trjconv -s {PULL_TPR} -f {PULL_XTC} -o {gro_out} "
            f"-dump {w['time_ps']:.1f} -n {INDEX}",
            input_text="0\n")
        
        if os.path.exists(gro_out):
            print(f"  Window {i:02d}: extracted at t={w['time_ps']:.1f} ps, "
                  f"pos={w['actual_nm']:.3f} nm")
        else:
            print(f"  Window {i:02d}: EXTRACTION FAILED at t={w['time_ps']:.1f} ps")
    
    print(f"\nGenerating grompp TPR files for umbrella windows...")
    for i, w in enumerate(windows):
        wdir = os.path.join(WINDOW_DIR, f"window_{i:02d}")
        gro_in = os.path.join(wdir, "conf.gro")
        tpr_out = os.path.join(wdir, "umbrella.tpr")
        
        if os.path.exists(tpr_out):
            print(f"  Window {i:02d}: umbrella.tpr exists, skipping grompp")
            continue
        
        if not os.path.exists(gro_in):
            print(f"  Window {i:02d}: missing conf.gro, skipping")
            continue
        
        run(f"gmx grompp -f {UMBRELLA_MDP} -c {gro_in} -r {gro_in} "
            f"-p {TOPOL} -n {INDEX} -o {tpr_out} -maxwarn 3")
        
        if os.path.exists(tpr_out):
            print(f"  Window {i:02d}: grompp OK")
        else:
            print(f"  Window {i:02d}: grompp FAILED")
    
    # Write summary manifest
    manifest_path = os.path.join(WINDOW_DIR, "manifest.txt")
    with open(manifest_path, "w") as mf:
        mf.write(f"# Umbrella sampling windows extracted from SMD pull\n")
        mf.write(f"# Pull trajectory: {PULL_XTC}\n")
        mf.write(f"# Spacing: {SPACING_NM} nm\n")
        mf.write(f"# Date: {subprocess.run('date', capture_output=True, text=True, shell=True).stdout.strip()}\n")
        mf.write(f"#\n")
        mf.write(f"# window  target_nm  actual_nm  time_ps\n")
        for i, w in enumerate(windows):
            mf.write(f"  {i:02d}      {w['target_nm']:.3f}      {w['actual_nm']:.3f}      {w['time_ps']:.1f}\n")
    
    # Write batch launch script
    launch_path = os.path.join(WINDOW_DIR, "run_all_umbrella.sh")
    with open(launch_path, "w") as lf:
        lf.write("#!/bin/bash\n")
        lf.write("# Launch all umbrella sampling windows\n")
        lf.write("# Each window: 10 ns, k=1000 kJ/mol/nm^2\n")
        lf.write("# Run from the umbrella_windows/ directory\n\n")
        lf.write("NTOMP=${1:-12}  # threads per window\n\n")
        for i in range(len(windows)):
            wdir = f"window_{i:02d}"
            lf.write(f"echo '=== Window {i:02d} / {len(windows)-1} ==='\n")
            lf.write(f"cd {wdir}\n")
            lf.write(f"if [ ! -f umbrella.gro ]; then\n")
            lf.write(f"  gmx mdrun -s umbrella.tpr -deffnm umbrella -ntomp $NTOMP -pin on\n")
            lf.write(f"fi\n")
            lf.write(f"cd ..\n\n")
    os.chmod(launch_path, 0o755)
    
    # Write parallel launch script (background all windows)
    par_path = os.path.join(WINDOW_DIR, "run_umbrella_parallel.sh")
    with open(par_path, "w") as pf:
        pf.write("#!/bin/bash\n")
        pf.write("# Launch umbrella windows in parallel batches\n")
        pf.write("# Usage: ./run_umbrella_parallel.sh [threads_per_window] [max_parallel]\n\n")
        pf.write("NTOMP=${1:-4}\n")
        pf.write("MAX_PAR=${2:-4}\n\n")
        pf.write("RUNNING=0\n")
        pf.write("PIDS=()\n\n")
        for i in range(len(windows)):
            wdir = f"window_{i:02d}"
            pf.write(f'if [ ! -f {wdir}/umbrella.gro ]; then\n')
            pf.write(f'  echo "Launching window {i:02d}..."\n')
            pf.write(f'  cd {wdir}\n')
            pf.write(f'  gmx mdrun -s umbrella.tpr -deffnm umbrella -ntomp $NTOMP -pin on > mdrun.log 2>&1 &\n')
            pf.write(f'  PIDS+=($!)\n')
            pf.write(f'  cd ..\n')
            pf.write(f'  RUNNING=$((RUNNING + 1))\n')
            pf.write(f'  if [ $RUNNING -ge $MAX_PAR ]; then\n')
            pf.write(f'    wait ${{PIDS[0]}}\n')
            pf.write(f'    PIDS=("${{PIDS[@]:1}}")\n')
            pf.write(f'    RUNNING=$((RUNNING - 1))\n')
            pf.write(f'  fi\n')
            pf.write(f'fi\n\n')
        pf.write("# Wait for remaining\n")
        pf.write("for pid in ${PIDS[@]}; do wait $pid; done\n")
        pf.write('echo "All umbrella windows complete."\n')
    os.chmod(par_path, 0o755)
    
    print(f"\n{'=' * 60}")
    print(f"EXTRACTION COMPLETE")
    print(f"  Windows: {len(windows)}")
    print(f"  Spacing: {SPACING_NM} nm")
    print(f"  Position range: {windows[0]['target_nm']:.3f} to {windows[-1]['target_nm']:.3f} nm")
    print(f"  Output dir: {WINDOW_DIR}")
    print(f"  Manifest: {manifest_path}")
    print(f"  Sequential launch: {launch_path}")
    print(f"  Parallel launch: {par_path}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
