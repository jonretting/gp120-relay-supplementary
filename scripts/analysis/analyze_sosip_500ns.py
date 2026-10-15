#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Comprehensive 500ns SOSIP MD analysis.

Produces a JSON report matching the 5U7O_500ns_analysis.json format:
- Pocket RMSD (internal, CA) with 100ns segments
- Backbone RMSD with 100ns segments
- Per-residue RMSF for relay pocket
- Pocket Rg with 100ns segments
- Pocket SASA with 100ns segments
- Temperature + potential energy
- Intra-pocket distances
- Performance and completion metadata
"""
import sys
import os
import subprocess
import json
import math
import re
from datetime import datetime

PDB_ID = sys.argv[1] if len(sys.argv) > 1 else "6MU6"
TPR = sys.argv[2] if len(sys.argv) > 2 else "md_500ns.tpr"
XTC = sys.argv[3] if len(sys.argv) > 3 else "md.xtc"
EDR = sys.argv[4] if len(sys.argv) > 4 else "md.edr"
GRO = sys.argv[5] if len(sys.argv) > 5 else "md.gro"

aa3to1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','HIE':'H','HID':'H'
}

def run_cmd(cmd, check=True):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if check and r.returncode != 0:
        print(f"  WARNING: Command failed: {cmd[:80]}...")
        print(f"  stderr: {r.stderr[-300:]}")
    return r

def parse_xvg(filename):
    vals = []
    if not os.path.exists(filename):
        return vals
    with open(filename) as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            parts = line.split()
            if len(parts) >= 2:
                vals.append((float(parts[0]), float(parts[1])))
    return vals

def stats(values):
    if not values:
        return {'mean': 0, 'std': 0, 'min': 0, 'max': 0, 'n': 0}
    n = len(values)
    mean = sum(values) / n
    var = sum((v - mean)**2 for v in values) / n if n > 1 else 0
    return {'mean': mean, 'std': math.sqrt(var), 'min': min(values), 'max': max(values), 'n': n}

def segment_stats(data_pairs, nm_to_a=True):
    """Compute 100ns segment statistics from (time_ns, value) pairs."""
    segments = []
    scale = 10.0 if nm_to_a else 1.0
    for seg_start in range(0, 500, 100):
        seg_end = seg_start + 100
        seg_vals = [v * scale for t, v in data_pairs if seg_start <= t < seg_end]
        if seg_vals:
            s = stats(seg_vals)
            segments.append({
                "ns": f"{seg_start}-{seg_end}",
                "mean": round(s['mean'], 3),
                "std": round(s['std'], 3)
            })
    return segments

print(f"=== {PDB_ID} 500ns SOSIP MD Analysis ===")
print(f"TPR: {TPR}, XTC: {XTC}")

# Step 1: Find relay pocket
print("\n--- Step 1: Identify relay pocket ---")
with open(GRO) as f:
    gro_lines = f.readlines()

residues = {}
for line in gro_lines[2:-1]:
    if len(line) < 44:
        continue
    resnum = int(line[0:5].strip())
    resname = line[5:10].strip()
    if resnum not in residues:
        residues[resnum] = resname

sorted_res = sorted(residues.items())
seq = ""
res_list = []
for rnum, rname in sorted_res:
    aa = aa3to1.get(rname, 'X')
    seq += aa
    res_list.append((rnum, rname, aa))

motif_idx = -1
target_motif = None
for variant in ["RPGGGDMR", "RPAGGDMR"]:
    motif_idx = seq.find(variant)
    if motif_idx >= 0:
        target_motif = variant
        break

if motif_idx < 0:
    m = re.search(r'RP[GA]GGDMR', seq)
    if m:
        motif_idx = m.start()
        target_motif = m.group()

if motif_idx < 0:
    print("  ERROR: Relay motif not found!")
    sys.exit(1)

relay_residues = []
relay_resnames = {}
for i in range(8):
    rnum, rname, aa = res_list[motif_idx + i]
    relay_residues.append(rnum)
    relay_resnames[rnum] = rname
    print(f"  {aa} = res {rnum} ({rname})")

print(f"  Motif: {target_motif}")

# Step 2: Create index groups
print("\n--- Step 2: Create index file ---")
pocket_atoms = []
pocket_ca_atoms = []
all_ca_atoms = []

for i, line in enumerate(gro_lines[2:-1], start=1):
    if len(line) < 44:
        continue
    resnum = int(line[0:5].strip())
    atomname = line[10:15].strip()
    if atomname == "CA":
        all_ca_atoms.append(i)
    if resnum in relay_residues:
        pocket_atoms.append(i)
        if atomname == "CA":
            pocket_ca_atoms.append(i)

with open("pocket.ndx", "w") as f:
    f.write("[ Relay_Pocket ]\n")
    for j, a in enumerate(pocket_atoms):
        f.write(f"{a:6d}")
        if (j+1) % 15 == 0: f.write("\n")
    f.write("\n[ Relay_Pocket_CA ]\n")
    for j, a in enumerate(pocket_ca_atoms):
        f.write(f"{a:6d}")
        if (j+1) % 15 == 0: f.write("\n")
    f.write("\n[ All_CA ]\n")
    for j, a in enumerate(all_ca_atoms):
        f.write(f"{a:6d}")
        if (j+1) % 15 == 0: f.write("\n")
    f.write("\n")

print(f"  Pocket atoms: {len(pocket_atoms)}, CA: {len(pocket_ca_atoms)}")
print(f"  Total CA: {len(all_ca_atoms)}, Total residues: {len(sorted_res)}")

# Step 3: Pocket RMSD
print("\n--- Step 3: Pocket RMSD (internal, CA) ---")
run_cmd(f"echo 'Relay_Pocket_CA Relay_Pocket_CA' | gmx rms -s {TPR} -f {XTC} -n pocket.ndx -o pocket_rmsd.xvg -tu ns 2>&1")

pocket_rmsd_data = parse_xvg("pocket_rmsd.xvg")
pocket_rmsd_vals = [v * 10.0 for t, v in pocket_rmsd_data]
pocket_rmsd_st = stats(pocket_rmsd_vals)
print(f"  Mean: {pocket_rmsd_st['mean']:.3f} +/- {pocket_rmsd_st['std']:.3f} A ({pocket_rmsd_st['n']} frames)")

# Step 4: Backbone RMSD
print("\n--- Step 4: Backbone RMSD ---")
run_cmd(f"echo 'Backbone Backbone' | gmx rms -s {TPR} -f {XTC} -o backbone_rmsd.xvg -tu ns 2>&1")

backbone_rmsd_data = parse_xvg("backbone_rmsd.xvg")
backbone_rmsd_vals = [v * 10.0 for t, v in backbone_rmsd_data]
backbone_rmsd_st = stats(backbone_rmsd_vals)
print(f"  Mean: {backbone_rmsd_st['mean']:.3f} +/- {backbone_rmsd_st['std']:.3f} A")

# Step 5: RMSF (per-residue, relay pocket)
print("\n--- Step 5: Per-residue RMSF ---")
run_cmd(f"echo 'Relay_Pocket_CA' | gmx rmsf -s {TPR} -f {XTC} -n pocket.ndx -o pocket_rmsf.xvg -res 2>&1")

rmsf_data = parse_xvg("pocket_rmsf.xvg")
rmsf_per_residue = {}
all_rmsf_vals = []
for resnum_float, rmsf_nm in rmsf_data:
    resnum = int(round(resnum_float))
    rmsf_a = rmsf_nm * 10.0
    if resnum in relay_resnames:
        rname = relay_resnames[resnum]
        rmsf_per_residue[f"{rname}{resnum}"] = round(rmsf_a, 3)
    all_rmsf_vals.append(rmsf_a)

# System-wide RMSF for percentile calculation
run_cmd(f"echo 'All_CA' | gmx rmsf -s {TPR} -f {XTC} -n pocket.ndx -o system_rmsf.xvg -res 2>&1")
system_rmsf_data = parse_xvg("system_rmsf.xvg")
system_rmsf_vals = [v * 10.0 for _, v in system_rmsf_data]

relay_mean_rmsf = sum(all_rmsf_vals) / len(all_rmsf_vals) if all_rmsf_vals else 0
system_mean_rmsf = sum(system_rmsf_vals) / len(system_rmsf_vals) if system_rmsf_vals else 0
if system_rmsf_vals:
    below = sum(1 for v in system_rmsf_vals if v < relay_mean_rmsf)
    relay_percentile = round(100.0 * below / len(system_rmsf_vals), 1)
else:
    relay_percentile = 0

print(f"  Relay mean RMSF: {relay_mean_rmsf:.3f} A (percentile: {relay_percentile})")
for k, v in rmsf_per_residue.items():
    print(f"    {k}: {v:.3f} A")

# Step 6: Pocket Rg
print("\n--- Step 6: Radius of gyration ---")
run_cmd(f"echo 'Relay_Pocket_CA' | gmx gyrate -s {TPR} -f {XTC} -n pocket.ndx -o pocket_rg.xvg -tu ns 2>&1")

pocket_rg_data = parse_xvg("pocket_rg.xvg")
pocket_rg_vals = [v * 10.0 for t, v in pocket_rg_data]
pocket_rg_st = stats(pocket_rg_vals)
print(f"  Mean: {pocket_rg_st['mean']:.2f} +/- {pocket_rg_st['std']:.2f} A")

# Step 7: SASA
print("\n--- Step 7: SASA ---")
run_cmd(f"echo 'Relay_Pocket' | gmx sasa -s {TPR} -f {XTC} -n pocket.ndx -o pocket_sasa.xvg -tu ns 2>&1")

sasa_data = parse_xvg("pocket_sasa.xvg")
sasa_vals = [v * 100.0 for t, v in sasa_data]  # nm^2 -> A^2
sasa_st = stats(sasa_vals)
print(f"  Mean: {sasa_st['mean']:.1f} +/- {sasa_st['std']:.1f} A^2 ({sasa_st['n']} frames)")

# Step 8: Temperature and energy
print("\n--- Step 8: Energy ---")
run_cmd("echo 'Temperature\n' | gmx energy -f " + EDR + " -o temperature.xvg 2>&1")
run_cmd("echo 'Potential\n' | gmx energy -f " + EDR + " -o potential.xvg 2>&1")

temp_data = parse_xvg("temperature.xvg")
temp_vals = [v for _, v in temp_data]
temp_st = stats(temp_vals)

pe_data = parse_xvg("potential.xvg")
pe_vals = [v for _, v in pe_data]
pe_st = stats(pe_vals)
pe_drift = pe_vals[-1] - pe_vals[0] if pe_vals else 0

print(f"  Temperature: {temp_st['mean']:.3f} +/- {temp_st['std']:.3f} K")
print(f"  Potential: {pe_st['mean']:.1f} kJ/mol (drift {pe_drift:.1f})")

# Step 9: Intra-pocket distances
print("\n--- Step 9: Intra-pocket distances ---")
def get_pocket_ca_from_gro(filename, relay_res):
    coords = {}
    try:
        with open(filename) as f:
            lines = f.readlines()
        for line in lines[2:-1]:
            if len(line) < 44: continue
            resnum = int(line[0:5].strip())
            atomname = line[10:15].strip()
            if resnum in relay_res and atomname == "CA":
                x = float(line[20:28]) * 10.0
                y = float(line[28:36]) * 10.0
                z = float(line[36:44]) * 10.0
                coords[resnum] = (x, y, z)
    except FileNotFoundError:
        pass
    return coords

end_coords = get_pocket_ca_from_gro(GRO, relay_residues)

# Extract first frame
run_cmd(f"echo 'System' | gmx trjconv -s {TPR} -f {XTC} -o first_frame.gro -dump 0 2>&1")
start_coords = get_pocket_ca_from_gro("first_frame.gro", relay_residues)

dist_changes = []
if start_coords and end_coords:
    for i, r1 in enumerate(relay_residues):
        for r2 in relay_residues[i+1:]:
            if r1 in start_coords and r2 in start_coords and r1 in end_coords and r2 in end_coords:
                s1, s2 = start_coords[r1], start_coords[r2]
                e1, e2 = end_coords[r1], end_coords[r2]
                d_s = math.sqrt(sum((a-b)**2 for a,b in zip(s1, s2)))
                d_e = math.sqrt(sum((a-b)**2 for a,b in zip(e1, e2)))
                dist_changes.append(abs(d_e - d_s))

dist_st = stats(dist_changes) if dist_changes else {'mean': 0, 'max': 0}
print(f"  Mean distance change: {dist_st['mean']:.2f} A, Max: {dist_st['max']:.2f} A")

# Step 10: Parse performance from md.log
print("\n--- Step 10: Performance ---")
perf = 0.0
completed = ""
with open("md.log") as f:
    for line in f:
        if "Performance:" in line:
            parts = line.split()
            for i, p in enumerate(parts):
                if p == "Performance:":
                    perf = float(parts[i+1])
        if "Finished mdrun" in line:
            completed = line.strip().split("rank 0 ")[-1] if "rank 0" in line else line.strip()

print(f"  Performance: {perf:.2f} ns/day")
print(f"  Completed: {completed}")

# Step 11: Equilibration check
first_100_vals = [v * 10.0 for t, v in pocket_rmsd_data if t < 100]
last_100_vals = [v * 10.0 for t, v in pocket_rmsd_data if t >= 400]
first_mean = sum(first_100_vals)/len(first_100_vals) if first_100_vals else 0
last_mean = sum(last_100_vals)/len(last_100_vals) if last_100_vals else 0
drift = last_mean - first_mean

# Build comprehensive report
print(f"\n{'='*60}")
print(f"  {PDB_ID} 500ns SOSIP -- FINAL REPORT")
print(f"{'='*60}")
print(f"  Pocket RMSD: {pocket_rmsd_st['mean']:.3f} +/- {pocket_rmsd_st['std']:.3f} A")
print(f"  Backbone RMSD: {backbone_rmsd_st['mean']:.3f} +/- {backbone_rmsd_st['std']:.3f} A")
print(f"  Relay RMSF: {relay_mean_rmsf:.3f} A ({relay_percentile}th percentile)")
print(f"  Pocket Rg: {pocket_rg_st['mean']:.2f} +/- {pocket_rg_st['std']:.2f} A")
print(f"  Pocket SASA: {sasa_st['mean']:.1f} +/- {sasa_st['std']:.1f} A^2")
print(f"  Frames: {pocket_rmsd_st['n']}")

report = {
    "pdb_id": PDB_ID,
    "simulation": {
        "duration_ns": 500,
        "timestep_fs": 2,
        "force_field": "amber99sb-ildn",
        "water_model": "tip3p",
        "temperature_K": 300,
        "salt_conc_M": 0.15,
        "performance_ns_day": round(perf, 2),
        "completed": completed
    },
    "relay_pocket": {
        "residues_gromacs": relay_residues,
        "motif": target_motif,
        "n_atoms": len(pocket_atoms),
        "n_ca_atoms": len(pocket_ca_atoms)
    },
    "pocket_rmsd_angstrom": {
        "mean": round(pocket_rmsd_st['mean'], 3),
        "std": round(pocket_rmsd_st['std'], 3),
        "min": round(pocket_rmsd_st['min'], 3),
        "max": round(pocket_rmsd_st['max'], 3),
        "n_frames": pocket_rmsd_st['n'],
        "segments_100ns": segment_stats(pocket_rmsd_data)
    },
    "backbone_rmsd_angstrom": {
        "mean": round(backbone_rmsd_st['mean'], 3),
        "std": round(backbone_rmsd_st['std'], 3),
        "max": round(backbone_rmsd_st['max'], 3),
        "segments_100ns": segment_stats(backbone_rmsd_data)
    },
    "rmsf_angstrom": {
        "relay_pocket": rmsf_per_residue,
        "relay_mean": round(relay_mean_rmsf, 3),
        "system_mean": round(system_mean_rmsf, 3),
        "relay_percentile": relay_percentile,
        "total_residues": len(sorted_res)
    },
    "pocket_rg_angstrom": {
        "mean": round(pocket_rg_st['mean'], 2),
        "std": round(pocket_rg_st['std'], 2),
        "segments_100ns": segment_stats(pocket_rg_data)
    },
    "pocket_sasa_angstrom2": {
        "mean": round(sasa_st['mean'], 1),
        "std": round(sasa_st['std'], 1),
        "n_frames": sasa_st['n'],
        "segments_100ns": segment_stats(sasa_data, nm_to_a=False)
    },
    "temperature_K": {
        "mean": round(temp_st['mean'], 3),
        "std": round(temp_st['std'], 3)
    },
    "potential_energy_kJ_mol": {
        "mean": round(pe_st['mean'], 1),
        "drift": round(pe_drift, 1)
    },
    "equilibration": {
        "equilibrated": abs(drift) < 1.0,
        "first_last_drift_angstrom": round(drift, 2)
    },
    "reference_comparison": {
        "ref_pdb": "6MEO",
        "ref_pocket_rmsd_20ns": 0.708,
        "ref_pocket_rmsd_std": 0.145,
        "self_pocket_rmsd_20ns": round(stats([v*10.0 for t,v in pocket_rmsd_data if t < 20])['mean'], 3) if pocket_rmsd_data else 0,
        "self_pocket_rmsd_20ns_std": round(stats([v*10.0 for t,v in pocket_rmsd_data if t < 20])['std'], 3) if pocket_rmsd_data else 0,
        "ratio_vs_6meo": round(pocket_rmsd_st['mean'] / 0.708, 2) if pocket_rmsd_st['mean'] > 0 else 0
    },
    "intra_pocket_distances": {
        "mean_change_angstrom": round(dist_st['mean'], 2),
        "max_change_angstrom": round(dist_st['max'], 2)
    }
}

# Fix SASA segments -- already in nm^2 from GROMACS, need A^2
sasa_segments = []
for seg_start in range(0, 500, 100):
    seg_end = seg_start + 100
    seg_vals = [v * 100.0 for t, v in sasa_data if seg_start <= t < seg_end]
    if seg_vals:
        s = stats(seg_vals)
        sasa_segments.append({"ns": f"{seg_start}-{seg_end}", "mean": round(s['mean'], 1), "std": round(s['std'], 1)})
report["pocket_sasa_angstrom2"]["segments_100ns"] = sasa_segments

outfile = f"{PDB_ID}_500ns_analysis.json"
with open(outfile, "w") as f:
    json.dump(report, f, indent=2)

print(f"\n  Report saved: {outfile}")
print("=== Analysis complete ===")
