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
Validate and analyze 100ns ligand-bound MD trajectory.
Uses MDAnalysis for trajectory analysis + GROMACS for energy checks.
"""

import os
import json
import subprocess
import numpy as np
from datetime import datetime
from pathlib import Path

WORK_DIR = Path(os.path.join(DATA_ROOT, "ligand-bound"))
OUTDIR = WORK_DIR / "validation"
OUTDIR.mkdir(exist_ok=True)

def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)

def run_gmx(cmd, input_str=None):
    """Run a GROMACS command with optional stdin."""
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True,
        input=input_str, timeout=300, cwd=str(WORK_DIR)
    )
    return result.stdout, result.stderr

log("=" * 70)
log("TRAJECTORY VALIDATION: 100ns Ligand-Bound MD")
log("=" * 70)

results = {"timestamp": datetime.now().isoformat(), "simulation": "100ns ligand-bound MD"}

# 1. Trajectory integrity
log("Step 1: Trajectory integrity")
stdout, stderr = run_gmx("gmx check -f md_prod.xtc 2>&1")
output = stdout + stderr
for line in output.split("\n"):
    if "Last frame" in line or "Step" in line or "frame" in line.lower():
        log(f"  {line.strip()}")
results["trajectory_check"] = "PASS" if "Last" in output else "UNKNOWN"

# 2. Energy stability
log("Step 2: Energy extraction")
stdout, stderr = run_gmx("gmx energy -f md_prod.edr -o /dev/null 2>&1", input_str="10 11 12 0\n")
energy_lines = (stdout + stderr).split("\n")
energy_stats = {}
for line in energy_lines:
    for key in ["Potential", "Kinetic", "Total", "Temperature", "Pressure"]:
        if key in line and "Energy" not in line:
            parts = line.split()
            try:
                avg_idx = next(i for i, p in enumerate(parts) if p == key) + 1
                if avg_idx < len(parts):
                    energy_stats[key] = float(parts[avg_idx])
            except (StopIteration, ValueError, IndexError):
                pass

log(f"  Energy stats: {energy_stats}")
results["energy"] = energy_stats

# 3. MDAnalysis-based analysis
log("Step 3: MDAnalysis trajectory analysis")

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align

    tpr = str(WORK_DIR / "md_prod.tpr")
    xtc = str(WORK_DIR / "md_prod.xtc")

    u = mda.Universe(tpr, xtc)
    log(f"  Loaded: {u.atoms.n_atoms} atoms, {u.trajectory.n_frames} frames, dt={u.trajectory.dt} ps")

    results["system"] = {
        "n_atoms": u.atoms.n_atoms,
        "n_frames": u.trajectory.n_frames,
        "dt_ps": float(u.trajectory.dt),
        "total_time_ns": float(u.trajectory.n_frames * u.trajectory.dt / 1000),
    }

    # Identify protein and ligand
    protein = u.select_atoms("protein")
    log(f"  Protein: {protein.n_atoms} atoms, {protein.n_residues} residues")

    # Try to find the ligand (non-protein, non-water, non-ion)
    ligand = u.select_atoms("not protein and not resname SOL NA CL WAT HOH TIP3")
    if ligand.n_atoms > 0:
        lig_resnames = set(ligand.residues.resnames)
        log(f"  Ligand: {ligand.n_atoms} atoms, resnames: {lig_resnames}")
        results["ligand"] = {"n_atoms": ligand.n_atoms, "resnames": list(lig_resnames)}
    else:
        log("  No ligand atoms found by exclusion filter")
        ligand = None

    # 4. Backbone RMSD over time
    log("Step 4: Backbone RMSD")
    backbone = u.select_atoms("backbone")

    ref = u.copy()
    ref.trajectory[0]

    rmsd_data = []
    for ts in u.trajectory:
        align.alignto(u, ref, select="backbone")
        r = rms.rmsd(backbone.positions, ref.select_atoms("backbone").positions)
        rmsd_data.append({"time_ns": float(ts.time / 1000), "rmsd_A": float(r)})

    rmsd_values = [d["rmsd_A"] for d in rmsd_data]
    results["backbone_rmsd"] = {
        "mean_A": float(np.mean(rmsd_values)),
        "std_A": float(np.std(rmsd_values)),
        "max_A": float(np.max(rmsd_values)),
        "final_A": float(rmsd_values[-1]),
        "n_points": len(rmsd_values),
        "stable": float(np.std(rmsd_values[-len(rmsd_values)//4:])) < 0.5
    }
    log(f"  RMSD: mean={np.mean(rmsd_values):.2f} +/- {np.std(rmsd_values):.2f} A, max={np.max(rmsd_values):.2f} A")

    # 5. Per-residue RMSF
    log("Step 5: Per-residue RMSF (CA)")
    ca = u.select_atoms("name CA")

    # Align all frames
    align.AlignTraj(u, ref, select="backbone", in_memory=False).run()

    rmsf_values = []
    positions = np.zeros((u.trajectory.n_frames, ca.n_atoms, 3))
    for i, ts in enumerate(u.trajectory):
        positions[i] = ca.positions

    mean_pos = positions.mean(axis=0)
    rmsf = np.sqrt(np.mean(np.sum((positions - mean_pos)**2, axis=2), axis=0))

    ca_resnums = ca.resids
    ca_resnames = ca.resnames

    # Identify relay pocket residues (check both numbering schemes)
    pocket_6meo = set(range(465, 473))
    pocket_sosip = set(range(469, 477))
    pocket_reduced = set(range(69, 77))

    pocket_rmsf = []
    for i, (resnum, resname, rf) in enumerate(zip(ca_resnums, ca_resnames, rmsf)):
        if resnum in pocket_6meo or resnum in pocket_sosip or resnum in pocket_reduced:
            pocket_rmsf.append({"residue": int(resnum), "resname": str(resname), "rmsf_A": float(rf)})

    results["rmsf"] = {
        "mean_all_A": float(np.mean(rmsf)),
        "pocket_residues": pocket_rmsf,
        "pocket_mean_A": float(np.mean([p["rmsf_A"] for p in pocket_rmsf])) if pocket_rmsf else None,
        "most_flexible_5": sorted(
            [{"residue": int(r), "resname": str(n), "rmsf_A": float(f)} for r, n, f in zip(ca_resnums, ca_resnames, rmsf)],
            key=lambda x: x["rmsf_A"], reverse=True
        )[:5]
    }
    log(f"  RMSF: mean={np.mean(rmsf):.2f} A, pocket mean={results['rmsf']['pocket_mean_A']:.2f} A" if pocket_rmsf else f"  RMSF: mean={np.mean(rmsf):.2f} A")

    # 6. Ligand RMSD if present
    if ligand is not None and ligand.n_atoms > 0:
        log("Step 6: Ligand RMSD")
        u.trajectory[0]
        lig_ref_pos = ligand.positions.copy()

        lig_rmsd = []
        for ts in u.trajectory:
            align.alignto(u, ref, select="backbone")
            r = rms.rmsd(ligand.positions, lig_ref_pos)
            lig_rmsd.append({"time_ns": float(ts.time / 1000), "rmsd_A": float(r)})

        lig_vals = [d["rmsd_A"] for d in lig_rmsd]
        results["ligand_rmsd"] = {
            "mean_A": float(np.mean(lig_vals)),
            "std_A": float(np.std(lig_vals)),
            "max_A": float(np.max(lig_vals)),
            "stable": float(np.std(lig_vals[-len(lig_vals)//4:])) < 1.0
        }
        log(f"  Ligand RMSD: mean={np.mean(lig_vals):.2f} +/- {np.std(lig_vals):.2f} A")

    # 7. Pocket-ligand contacts
    if ligand is not None and pocket_rmsf:
        log("Step 7: Pocket-ligand contacts over time")
        pocket_sel_parts = []
        for numbering in [pocket_6meo, pocket_sosip, pocket_reduced]:
            resid_str = " ".join(str(r) for r in numbering)
            pocket_sel_parts.append(f"(protein and resid {resid_str})")
        
        contacts_over_time = []
        for numbering in [pocket_6meo, pocket_sosip, pocket_reduced]:
            resid_str = " ".join(str(r) for r in numbering)
            try:
                pocket_atoms = u.select_atoms(f"protein and resid {resid_str}")
                if pocket_atoms.n_atoms > 0:
                    contact_counts = []
                    for ts in u.trajectory[::10]:  # every 10th frame
                        dists = np.sqrt(np.sum(
                            (pocket_atoms.positions[:, np.newaxis, :] - ligand.positions[np.newaxis, :, :])**2,
                            axis=2
                        ))
                        n_contacts = np.sum(dists < 4.0)
                        contact_counts.append({"time_ns": float(ts.time / 1000), "n_contacts": int(n_contacts)})
                    
                    if contact_counts:
                        avg_contacts = np.mean([c["n_contacts"] for c in contact_counts])
                        results["pocket_ligand_contacts"] = {
                            "avg_contacts_4A": float(avg_contacts),
                            "pocket_numbering": sorted(list(numbering)),
                            "n_pocket_atoms": pocket_atoms.n_atoms,
                            "n_frames_sampled": len(contact_counts)
                        }
                        log(f"  Avg contacts (4A): {avg_contacts:.1f} (numbering {min(numbering)}-{max(numbering)})")
                        break
            except Exception as e:
                continue

    # 8. Radius of gyration
    log("Step 8: Radius of gyration")
    rgyr = []
    for ts in u.trajectory[::10]:
        rg = protein.radius_of_gyration()
        rgyr.append(float(rg))
    
    results["radius_of_gyration"] = {
        "mean_A": float(np.mean(rgyr)),
        "std_A": float(np.std(rgyr)),
    }
    log(f"  Rg: mean={np.mean(rgyr):.2f} +/- {np.std(rgyr):.2f} A")

except ImportError as e:
    log(f"  MDAnalysis import error: {e}")
    results["mdanalysis_error"] = str(e)
except Exception as e:
    log(f"  Analysis error: {e}")
    import traceback
    traceback.print_exc()
    results["analysis_error"] = str(e)

# Save results
with open(OUTDIR / "validation_results.json", "w") as f:
    json.dump(results, f, indent=2)

log("")
log("=" * 70)
log("VALIDATION COMPLETE")
log(f"Results: {OUTDIR / 'validation_results.json'}")
log("=" * 70)

# Print summary
print("\nSUMMARY:")
print(f"  System: {results.get('system', {}).get('n_atoms', '?')} atoms, {results.get('system', {}).get('n_frames', '?')} frames")
print(f"  Total time: {results.get('system', {}).get('total_time_ns', '?')} ns")
print(f"  Backbone RMSD: {results.get('backbone_rmsd', {}).get('mean_A', '?'):.2f} +/- {results.get('backbone_rmsd', {}).get('std_A', '?'):.2f} A")
if "ligand_rmsd" in results:
    print(f"  Ligand RMSD: {results['ligand_rmsd']['mean_A']:.2f} +/- {results['ligand_rmsd']['std_A']:.2f} A")
print(f"  Rg: {results.get('radius_of_gyration', {}).get('mean_A', '?'):.2f} A")
print(f"  Stable: {results.get('backbone_rmsd', {}).get('stable', '?')}")
