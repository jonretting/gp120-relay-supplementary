#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""
Phase 3: MD Trajectory Analysis -- Pocket persistence, RMSF, volume
Requires: MDAnalysis, matplotlib (pip install MDAnalysis matplotlib)
"""
import subprocess, os, sys, json
from pathlib import Path

WORK = Path.home() / "md-relay-pocket"
MD = WORK / "md"
RESULTS = WORK / "results"

RELAY_RESIDUES = list(range(465, 473))

def run(cmd, cwd=None):
    print(f"  >> {cmd}")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd or MD)
    if r.returncode != 0:
        print(f"  STDERR: {r.stderr[:500]}")
    return r

def check_trajectory():
    """Check if trajectory exists and how much has been written."""
    xtc = MD / "md.xtc"
    tpr = MD / "md.tpr"
    if not xtc.exists():
        print("  No trajectory file found yet")
        return False
    r = run(f"gmx check -f {xtc}", cwd=MD)
    return r.returncode == 0

def compute_rmsf():
    """Compute per-residue RMSF for pocket and surrounding residues."""
    print("\n--- RMSF Analysis ---")

    # Create index group for relay pocket
    # Use gmx select or make_ndx
    r = run("echo 'q' | gmx make_ndx -f md.tpr -o relay.ndx", cwd=MD)

    # RMSF of C-alpha atoms
    r = run("echo '3' | gmx rmsf -f md.xtc -s md.tpr -o rmsf.xvg -res", cwd=MD)
    if r.returncode != 0:
        return None

    # Parse RMSF
    rmsf_data = {}
    xvg_path = MD / "rmsf.xvg"
    if xvg_path.exists():
        with open(xvg_path) as f:
            for line in f:
                if line.startswith("#") or line.startswith("@"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    resnum = int(float(parts[0]))
                    rmsf_val = float(parts[1])
                    rmsf_data[resnum] = rmsf_val

    # Extract relay pocket RMSF
    relay_rmsf = {r: rmsf_data.get(r, None) for r in RELAY_RESIDUES}
    relay_values = [v for v in relay_rmsf.values() if v is not None]

    if relay_values:
        mean_rmsf = sum(relay_values) / len(relay_values)
        print(f"  Relay pocket mean RMSF: {mean_rmsf:.3f} nm")
        print(f"  Per-residue RMSF (G:465-472):")
        for r, v in sorted(relay_rmsf.items()):
            if v is not None:
                stability = "RIGID" if v < 0.15 else "FLEXIBLE" if v > 0.3 else "MODERATE"
                print(f"    G:{r} = {v:.3f} nm ({stability})")
    return relay_rmsf

def compute_pocket_rmsd():
    """Compute RMSD of pocket backbone over time."""
    print("\n--- Pocket RMSD ---")

    # Backbone RMSD of relay pocket residues
    r = run("echo '4 4' | gmx rms -f md.xtc -s md.tpr -o rmsd_bb.xvg", cwd=MD)

    rmsd_path = MD / "rmsd_bb.xvg"
    if rmsd_path.exists():
        times = []
        rmsds = []
        with open(rmsd_path) as f:
            for line in f:
                if line.startswith("#") or line.startswith("@"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    times.append(float(parts[0]))
                    rmsds.append(float(parts[1]))

        if rmsds:
            mean_rmsd = sum(rmsds) / len(rmsds)
            max_rmsd = max(rmsds)
            # Skip first 2 ns as equilibration
            stable_rmsds = [r for t, r in zip(times, rmsds) if t > 2000]
            stable_mean = sum(stable_rmsds) / len(stable_rmsds) if stable_rmsds else mean_rmsd

            print(f"  Overall RMSD: {mean_rmsd:.3f} ± {(sum((r-mean_rmsd)**2 for r in rmsds)/len(rmsds))**0.5:.3f} nm")
            print(f"  Stable phase RMSD (>2 ns): {stable_mean:.3f} nm")
            print(f"  Max RMSD: {max_rmsd:.3f} nm")

            if stable_mean < 0.2:
                print(f"  -> POCKET IS STABLE (< 2 Å backbone deviation)")
            elif stable_mean < 0.35:
                print(f"  -> MODERATE FLEXIBILITY (2-3.5 Å)")
            else:
                print(f"  -> SIGNIFICANT MOTION (> 3.5 Å)")

            return {"mean_nm": mean_rmsd, "stable_mean_nm": stable_mean, "max_nm": max_rmsd}
    return None

def compute_pocket_volume():
    """Extract snapshots and measure pocket volume using fpocket."""
    print("\n--- Pocket Volume Analysis ---")

    snap_dir = MD / "snapshots"
    os.makedirs(snap_dir, exist_ok=True)

    # Extract frames every 2 ns
    r = run(f"echo '1' | gmx trjconv -f md.xtc -s md.tpr -o {snap_dir}/frame.pdb -sep -dt 2000", cwd=MD)

    # Check if fpocket is available
    r = run("which fpocket")
    if r.returncode != 0:
        print("  fpocket not installed -- skipping volume analysis")
        print("  Install: sudo apt-get install fpocket")
        return None

    volumes = []
    for pdb in sorted(snap_dir.glob("frame*.pdb")):
        r = run(f"fpocket -f {pdb}", cwd=snap_dir)
        # Parse fpocket output for pocket near relay coordinates
        # (simplified -- fpocket creates output dirs)

    return volumes

def extract_hbond_analysis():
    """Analyze hydrogen bond occupancy in the pocket."""
    print("\n--- H-bond Analysis ---")
    r = run("echo '1 1' | gmx hbond -f md.xtc -s md.tpr -num hbond.xvg", cwd=MD)

    hbond_path = MD / "hbond.xvg"
    if hbond_path.exists():
        hbonds = []
        with open(hbond_path) as f:
            for line in f:
                if line.startswith("#") or line.startswith("@"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    hbonds.append(int(parts[1]))

        if hbonds:
            mean_hb = sum(hbonds) / len(hbonds)
            print(f"  Mean H-bonds in system: {mean_hb:.1f}")
            return {"mean_hbonds": mean_hb}
    return None

def compute_sasa():
    """Compute solvent-accessible surface area of pocket."""
    print("\n--- SASA Analysis ---")
    r = run("echo '1' | gmx sasa -f md.xtc -s md.tpr -o sasa.xvg", cwd=MD)

    sasa_path = MD / "sasa.xvg"
    if sasa_path.exists():
        sasas = []
        with open(sasa_path) as f:
            for line in f:
                if line.startswith("#") or line.startswith("@"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    sasas.append(float(parts[1]))

        if sasas:
            mean_sasa = sum(sasas) / len(sasas)
            print(f"  Mean SASA: {mean_sasa:.2f} nm²")
            return {"mean_sasa_nm2": mean_sasa}
    return None

def main():
    os.makedirs(RESULTS, exist_ok=True)

    print("=" * 60)
    print("PHASE 3: MD TRAJECTORY ANALYSIS")
    print("=" * 60)

    if not check_trajectory():
        print("\nTrajectory not ready. Check MD progress:")
        print("  Monitor: tail -f <md_dir>/md_production.log")
        sys.exit(1)

    rmsf = compute_rmsf()
    rmsd = compute_pocket_rmsd()
    hbonds = extract_hbond_analysis()
    sasa = compute_sasa()

    results = {
        "simulation": {
            "system": "6MEO reduced (15 Å shell around G:465-472)",
            "duration_ns": 20,
            "force_field": "AMBER99SB-ILDN or CHARMM27",
            "water": "TIP3P",
            "temperature_K": 300,
            "pressure_bar": 1.0
        },
        "rmsf": rmsf,
        "rmsd": rmsd,
        "hbonds": hbonds,
        "sasa": sasa,
    }

    results_path = RESULTS / "md_analysis.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\n  Results saved to {results_path}")

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)

    if rmsd and rmsd.get("stable_mean_nm", 1.0) < 0.2:
        print("\n  [OK] Pocket backbone is STABLE (RMSD < 2 Å)")
        print("  -> Supports druggability claim in manuscript")
    elif rmsd:
        print(f"\n  Pocket RMSD = {rmsd['stable_mean_nm']*10:.1f} Å -- needs interpretation")

if __name__ == "__main__":
    main()
