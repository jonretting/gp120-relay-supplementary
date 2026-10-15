#!/usr/bin/env python3
# Author: Jonathan Michael Retting
# Study:  Identification of a Novel Druggable Pocket in the HIV-1 gp120
#         C5 Domain at the Conformational Relay Junction
# License: CC-BY-4.0 (https://creativecommons.org/licenses/by/4.0/)
"""Analyze SOSIP reduced-system MD trajectory.

Reads the .gro topology to find relay pocket residues (RPGGGDMR),
creates GROMACS index groups, runs analysis tools, and produces
a consolidated JSON report.
"""
import sys
import os
import subprocess
import json
import re
import math

PDB_ID = sys.argv[1]
WORKDIR = os.getcwd()

aa3to1 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V','HIE':'H','HID':'H'
}

print(f"=== SOSIP MD Analysis: {PDB_ID} ===")

# Step 1: Find relay pocket residues in the .gro file
print("\n--- Step 1: Identify relay pocket in topology ---")
with open("md.gro") as f:
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

motif_variants = ["RPGGGDMR", "RPAGGDMR"]
target = None
idx = -1
for variant in motif_variants:
    idx = seq.find(variant)
    if idx >= 0:
        target = variant
        break

if idx < 0:
    m = re.search(r'RP[GA]GGDMR', seq)
    if m:
        idx = m.start()
        target = m.group()

if idx < 0:
    print("  ERROR: Relay motif not found in sequence!")
    print(f"  Tried: {motif_variants}")
    print(f"  Sequence: {seq}")
    sys.exit(1)

if target != "RPGGGDMR":
    print(f"  NOTE: Found variant motif {target} (differs from canonical RPGGGDMR)")

relay_residues = []
for i in range(8):
    rnum, rname, aa = res_list[idx + i]
    relay_residues.append(rnum)
    print(f"  {aa} = res {rnum} ({rname})")

relay_start = relay_residues[0]
relay_end = relay_residues[-1]
print(f"  Relay pocket: residues {relay_start}-{relay_end}")

# Step 2: Create index file with relay pocket group
print("\n--- Step 2: Create index groups ---")
with open("md.gro") as f:
    gro_lines = f.readlines()

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
    for j, atom in enumerate(pocket_atoms):
        f.write(f"{atom:6d}")
        if (j + 1) % 15 == 0:
            f.write("\n")
    f.write("\n")
    f.write("[ Relay_Pocket_CA ]\n")
    for j, atom in enumerate(pocket_ca_atoms):
        f.write(f"{atom:6d}")
        if (j + 1) % 15 == 0:
            f.write("\n")
    f.write("\n")
    f.write("[ All_CA ]\n")
    for j, atom in enumerate(all_ca_atoms):
        f.write(f"{atom:6d}")
        if (j + 1) % 15 == 0:
            f.write("\n")
    f.write("\n")

print(f"  Pocket atoms: {len(pocket_atoms)}")
print(f"  Pocket CA atoms: {len(pocket_ca_atoms)}")
print(f"  Total CA atoms: {len(all_ca_atoms)}")

# Step 3: Run GROMACS analyses
print("\n--- Step 3: RMSD analysis ---")

# Pocket RMSD (internal, CA only)
cmd = "echo 'Relay_Pocket_CA Relay_Pocket_CA' | gmx rms -s md.tpr -f md.xtc -n pocket.ndx -o pocket_rmsd.xvg -tu ns 2>&1"
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
if result.returncode != 0:
    print(f"  WARNING: gmx rms failed: {result.stderr[-200:]}")

# Global backbone RMSD
cmd = "echo 'Backbone Backbone' | gmx rms -s md.tpr -f md.xtc -o backbone_rmsd.xvg -tu ns 2>&1"
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

# Pocket Rg
print("\n--- Step 4: Radius of gyration ---")
cmd = "echo 'Relay_Pocket_CA' | gmx gyrate -s md.tpr -f md.xtc -n pocket.ndx -o pocket_rg.xvg -tu ns 2>&1"
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

# Temperature and energy
print("\n--- Step 5: Energy extraction ---")
cmd = "echo 'Temperature\n' | gmx energy -f md.edr -o temperature.xvg 2>&1"
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
for line in result.stderr.split('\n'):
    if 'Temperature' in line and 'Av' not in line:
        continue
    if 'Temperature' in line:
        print(f"  {line.strip()}")

cmd = "echo 'Potential\n' | gmx energy -f md.edr -o potential.xvg 2>&1"
result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
for line in result.stderr.split('\n'):
    if 'Potential' in line and 'Av' not in line:
        continue
    if 'Potential' in line:
        print(f"  {line.strip()}")

# Step 6: Parse XVG files and compute statistics
print("\n--- Step 6: Parse results ---")

def parse_xvg(filename):
    vals = []
    if not os.path.exists(filename):
        return vals
    with open(filename) as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
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
    return {
        'mean': mean,
        'std': math.sqrt(var),
        'min': min(values),
        'max': max(values),
        'n': n
    }

# Pocket RMSD (convert nm to Angstrom)
pocket_rmsd_data = parse_xvg("pocket_rmsd.xvg")
pocket_rmsd_vals = [v[1] * 10.0 for v in pocket_rmsd_data]  # nm -> A
pocket_rmsd_stats = stats(pocket_rmsd_vals)

# Backbone RMSD
backbone_rmsd_data = parse_xvg("backbone_rmsd.xvg")
backbone_rmsd_vals = [v[1] * 10.0 for v in backbone_rmsd_data]
backbone_rmsd_stats = stats(backbone_rmsd_vals)

# Pocket Rg
pocket_rg_data = parse_xvg("pocket_rg.xvg")
pocket_rg_vals = [v[1] * 10.0 for v in pocket_rg_data]
pocket_rg_stats = stats(pocket_rg_vals)

# Temperature
temp_data = parse_xvg("temperature.xvg")
temp_vals = [v[1] for v in temp_data]
temp_stats = stats(temp_vals)

# Potential energy
pe_data = parse_xvg("potential.xvg")
pe_vals = [v[1] for v in pe_data]
pe_stats = stats(pe_vals)
pe_drift = pe_vals[-1] - pe_vals[0] if pe_vals else 0

# Intra-pocket CA-CA distances
print("\n--- Step 7: Intra-pocket distance analysis ---")

# Use first and last frames to check pocket geometry stability
cmd = "echo 'Relay_Pocket_CA' | gmx distance -s md.tpr -f md.xtc -n pocket.ndx -oall pocket_dist.xvg -select 'group \"Relay_Pocket_CA\" permute 1 2 plus group \"Relay_Pocket_CA\" permute 1 3' 2>&1"
# This is complex - let's just compute from the .gro coordinates directly

# Parse the first and last frames for pocket CA coordinates
def get_pocket_ca_from_gro(filename, relay_residues):
    coords = {}
    with open(filename) as f:
        lines = f.readlines()
    for line in lines[2:-1]:
        if len(line) < 44:
            continue
        resnum = int(line[0:5].strip())
        atomname = line[10:15].strip()
        if resnum in relay_residues and atomname == "CA":
            x = float(line[20:28]) * 10.0  # nm -> A
            y = float(line[28:36]) * 10.0
            z = float(line[36:44]) * 10.0
            coords[resnum] = (x, y, z)
    return coords

start_coords = get_pocket_ca_from_gro("npt.gro", relay_residues)
end_coords = get_pocket_ca_from_gro("md.gro", relay_residues)

if start_coords and end_coords:
    dist_changes = []
    for i, r1 in enumerate(relay_residues):
        for r2 in relay_residues[i+1:]:
            if r1 in start_coords and r2 in start_coords and r1 in end_coords and r2 in end_coords:
                s1, s2 = start_coords[r1], start_coords[r2]
                e1, e2 = end_coords[r1], end_coords[r2]
                d_start = math.sqrt(sum((a-b)**2 for a,b in zip(s1, s2)))
                d_end = math.sqrt(sum((a-b)**2 for a,b in zip(e1, e2)))
                dist_changes.append(abs(d_end - d_start))

    dist_change_stats = stats(dist_changes)
    print(f"  Mean CA-CA distance change (start->end): {dist_change_stats['mean']:.2f} A")
    print(f"  Max CA-CA distance change: {dist_change_stats['max']:.2f} A")

# Print summary
print(f"\n{'='*60}")
print(f"  SUMMARY: {PDB_ID} 20ns Reduced-System MD")
print(f"{'='*60}")
print(f"  Pocket Internal RMSD:  {pocket_rmsd_stats['mean']:.3f} +/- {pocket_rmsd_stats['std']:.3f} A (max {pocket_rmsd_stats['max']:.3f} A)")
print(f"  Global Backbone RMSD:  {backbone_rmsd_stats['mean']:.3f} +/- {backbone_rmsd_stats['std']:.3f} A")
print(f"  Pocket Rg:             {pocket_rg_stats['mean']:.2f} +/- {pocket_rg_stats['std']:.2f} A")
if pocket_rg_data:
    rg_start = pocket_rg_data[0][1] * 10.0
    rg_end = pocket_rg_data[-1][1] * 10.0
    rg_change = (rg_end - rg_start) / rg_start * 100
    print(f"  Pocket Rg change:      {rg_start:.2f} -> {rg_end:.2f} A ({rg_change:+.1f}%)")
print(f"  Temperature:           {temp_stats['mean']:.3f} +/- {temp_stats['std']:.3f} K")
print(f"  Potential Energy:      {pe_stats['mean']:.3e} kJ/mol (drift {pe_drift:.1f})")
print(f"  Frames analyzed:       {pocket_rmsd_stats['n']}")

# Compare to 6MEO reference
print(f"\n  --- Comparison to 6MEO reference ---")
print(f"  6MEO pocket RMSD: 0.708 +/- 0.145 A (20ns ternary)")
print(f"  {PDB_ID} pocket RMSD: {pocket_rmsd_stats['mean']:.3f} +/- {pocket_rmsd_stats['std']:.3f} A")
ratio = pocket_rmsd_stats['mean'] / 0.708 if pocket_rmsd_stats['mean'] > 0 else 0
print(f"  Ratio: {ratio:.2f}x")

# Save JSON report
report = {
    "pdb_id": PDB_ID,
    "simulation": {
        "duration_ns": 20,
        "timestep_fs": 2,
        "force_field": "amber99sb-ildn",
        "water_model": "tip3p",
        "temperature_K": 300,
        "salt_conc_M": 0.15
    },
    "relay_pocket": {
        "residues_gromacs": relay_residues,
        "motif": "RPGGGDMR",
        "n_atoms": len(pocket_atoms),
        "n_ca_atoms": len(pocket_ca_atoms)
    },
    "pocket_rmsd_angstrom": {
        "mean": round(pocket_rmsd_stats['mean'], 3),
        "std": round(pocket_rmsd_stats['std'], 3),
        "min": round(pocket_rmsd_stats['min'], 3),
        "max": round(pocket_rmsd_stats['max'], 3),
        "n_frames": pocket_rmsd_stats['n']
    },
    "backbone_rmsd_angstrom": {
        "mean": round(backbone_rmsd_stats['mean'], 3),
        "std": round(backbone_rmsd_stats['std'], 3),
        "max": round(backbone_rmsd_stats['max'], 3)
    },
    "pocket_rg_angstrom": {
        "mean": round(pocket_rg_stats['mean'], 2),
        "std": round(pocket_rg_stats['std'], 2)
    },
    "temperature_K": {
        "mean": round(temp_stats['mean'], 3),
        "std": round(temp_stats['std'], 3)
    },
    "potential_energy_kJ_mol": {
        "mean": round(pe_stats['mean'], 1),
        "drift": round(pe_drift, 1)
    },
    "reference_comparison": {
        "ref_pdb": "6MEO",
        "ref_pocket_rmsd": 0.708,
        "ref_pocket_rmsd_std": 0.145,
        "ratio": round(ratio, 2)
    }
}

if start_coords and end_coords:
    report["intra_pocket_distances"] = {
        "mean_change_angstrom": round(dist_change_stats['mean'], 2),
        "max_change_angstrom": round(dist_change_stats['max'], 2)
    }

with open(f"{PDB_ID}_md_analysis.json", "w") as f:
    json.dump(report, f, indent=2)

print(f"\n  Report saved: {PDB_ID}_md_analysis.json")
