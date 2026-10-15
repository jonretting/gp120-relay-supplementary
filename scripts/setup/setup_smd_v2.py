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
Setup Steered MD v2 -- Elongated rectangular box for CD4 pull from gp120.

Design rationale:
- Rectangular cuboid elongated in Z (pull direction)
- Align dissociation axis along Z
- Asymmetric centering (complex at Z/4)
- Slower pull rate: 0.001 nm/ps (was 0.01)
- geometry = direction (not distance)

Usage: Run from the supplementary data directory.
"""
import subprocess
import os
import sys

WORK = os.path.join(DATA_ROOT, "results", "05-smd-v2")
SRC_GRO = os.path.join(DATA_ROOT, "topology", "md.gro")
SRC_TOP = os.path.join(DATA_ROOT, "topology", "topol.top")
SRC_NDX = os.path.join(DATA_ROOT, "topology", "index.ndx")

os.makedirs(WORK, exist_ok=True)
os.chdir(WORK)

def run(cmd, input_text=None):
    print(f">>> {cmd}")
    r = subprocess.run(cmd, shell=True, input=input_text, capture_output=True, text=True)
    if r.stdout: print(r.stdout[-500:])
    if r.stderr: print(r.stderr[-500:])
    if r.returncode != 0:
        print(f"WARNING: exit code {r.returncode}")
    return r

print("=" * 60)
print("Step 1: Measure current system dimensions")
print("=" * 60)

with open(SRC_GRO) as f:
    lines = f.readlines()
box_line = lines[-1].split()
print(f"Current box: {' '.join(box_line[:3])}")
print(f"Total atoms: {lines[1].strip()}")

xs, ys, zs = [], [], []
for line in lines[2:-1]:
    if len(line) < 44:
        continue
    x = float(line[20:28])
    y = float(line[28:36])
    z = float(line[36:44])
    xs.append(x); ys.append(y); zs.append(z)

print(f"Protein extent: X={min(xs):.1f}-{max(xs):.1f} ({max(xs)-min(xs):.1f} nm)")
print(f"                Y={min(ys):.1f}-{max(ys):.1f} ({max(ys)-min(ys):.1f} nm)")
print(f"                Z={min(zs):.1f}-{max(zs):.1f} ({max(zs)-min(zs):.1f} nm)")

z_extent = max(zs) - min(zs)
xy_extent = max(max(xs)-min(xs), max(ys)-min(ys))

print()
print("=" * 60)
print("Step 2: Calculate new box dimensions")
print("=" * 60)

pull_distance = 5.0  # nm -- we want ~5 nm of additional separation
buffer = 2.5  # nm on each side

xy_box = xy_extent + 2 * 1.5  # 1.5 nm buffer each side in X/Y
z_box = z_extent + pull_distance + 2 * buffer + 5.0  # extra margin for safety

# Round up to nice numbers
import math
xy_box = math.ceil(xy_box)
z_box = math.ceil(z_box)

# Ensure z_box > 2 * (initial_com_dist + pull_distance) / 0.49
# Initial COM ~8 nm + 5 nm pull = 13 nm max COM. Need z > 13/0.49 = 26.5 nm
z_box = max(z_box, 28)

print(f"Target box: {xy_box} x {xy_box} x {z_box} nm (rectangular)")
print(f"Pull distance: {pull_distance} nm")
print(f"Max COM distance: ~{8 + pull_distance} nm")
print(f"0.49 * z_box = {0.49 * z_box:.1f} nm (must be > max COM dist)")

print()
print("=" * 60)
print("Step 3: Extract, align, and re-box")
print("=" * 60)

# Copy source topology
run(f"cp {SRC_TOP} {WORK}/topol.top")
for itp in ["topol_Protein_chain_A.itp", "topol_Protein_chain_B.itp", 
            "topol_Protein_chain_C.itp", "posre_Protein_chain_A.itp",
            "posre_Protein_chain_B.itp", "posre_Protein_chain_C.itp"]:
    src = os.path.join(os.path.dirname(SRC_TOP), itp)
    if os.path.exists(src):
        run(f"cp {src} {WORK}/")

# Copy index file
if os.path.exists(SRC_NDX):
    run(f"cp {SRC_NDX} {WORK}/index.ndx")

# Put complex in new rectangular box, centered asymmetrically (Z/4 from bottom)
z_center = z_box / 4.0
run(f"gmx editconf -f {SRC_GRO} -o boxed.gro -bt cubic "
    f"-box {xy_box} {xy_box} {z_box} "
    f"-center {xy_box/2:.1f} {xy_box/2:.1f} {z_center:.1f} -noc")

print()
print("=" * 60)
print("Step 4: Solvate and add ions")
print("=" * 60)

run("gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top 2>&1")

# Count current charge and add ions
run("gmx grompp -f /dev/null -c solvated.gro -p topol.top -o ions.tpr -maxwarn 10 2>&1 || true")

# Write a minimal MDP for grompp (ions step)
with open("ions.mdp", "w") as f:
    f.write("integrator = steep\nnsteps = 0\n")

run("gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 10 2>&1")
run("echo SOL | gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral -conc 0.15 2>&1")

print()
print("=" * 60)
print("Step 5: Write MDP files")
print("=" * 60)

# Energy minimization
with open("em.mdp", "w") as f:
    f.write("""\
integrator  = steep
emtol       = 500.0
emstep      = 0.01
nsteps      = 50000
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
""")

# NVT equilibration (100 ps, position restraints)
with open("nvt.mdp", "w") as f:
    f.write("""\
integrator  = md
dt          = 0.002
nsteps      = 50000      ; 100 ps
nstxout-compressed = 500
nstenergy   = 500
nstlog      = 500
continuation = no
gen-vel     = yes
gen-temp    = 300
gen-seed    = -1
tcoupl      = V-rescale
tc-grps     = Protein non-Protein
tau-t       = 0.1     0.1
ref-t       = 300     300
pcoupl      = no
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
constraints = h-bonds
constraint-algorithm = LINCS
define      = -DPOSRES
""")

# NPT equilibration (500 ps, position restraints, semi-isotropic for elongated box)
with open("npt.mdp", "w") as f:
    f.write("""\
integrator  = md
dt          = 0.002
nsteps      = 250000     ; 500 ps
nstxout-compressed = 500
nstenergy   = 500
nstlog      = 500
continuation = yes
gen-vel     = no
tcoupl      = V-rescale
tc-grps     = Protein non-Protein
tau-t       = 0.1     0.1
ref-t       = 300     300
pcoupl      = Parrinello-Rahman
pcoupltype  = semiisotropic
tau-p       = 5.0
ref-p       = 1.0 1.0
compressibility = 4.5e-5 4.5e-5
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
constraints = h-bonds
constraint-algorithm = LINCS
define      = -DPOSRES
""")

# SMD pull MDP (5 ns, 0.001 nm/ps along Z)
with open("pull.mdp", "w") as f:
    f.write("""\
integrator  = md
dt          = 0.002
nsteps      = 2500000    ; 5 ns
nstxout-compressed = 500 ; save every 1 ps
nstenergy   = 500
nstlog      = 500
continuation = yes
gen-vel     = no
tcoupl      = V-rescale
tc-grps     = Protein non-Protein
tau-t       = 0.1     0.1
ref-t       = 300     300
pcoupl      = Parrinello-Rahman
pcoupltype  = semiisotropic
tau-p       = 5.0
ref-p       = 1.0 1.0
compressibility = 4.5e-5 4.5e-5
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
constraints = h-bonds
constraint-algorithm = LINCS
nstcomm     = 100
comm-mode   = linear
comm-grps   = System

; Pull code -- direction-based along Z
pull                 = yes
pull-ncoords         = 1
pull-ngroups         = 2
pull-group1-name     = gp120
pull-group2-name     = CD4
pull-coord1-type     = umbrella
pull-coord1-geometry = direction
pull-coord1-vec      = 0 0 1
pull-coord1-dim      = N N Y
pull-coord1-groups   = 1 2
pull-coord1-start    = yes
pull-coord1-rate     = 0.001    ; 0.001 nm/ps = 1 nm/ns
pull-coord1-k        = 1000     ; kJ/mol/nm^2

pull-nstxout         = 50       ; write every 0.1 ps
pull-nstfout         = 50
""")

# Umbrella sampling window MDP template
with open("umbrella.mdp", "w") as f:
    f.write("""\
integrator  = md
dt          = 0.002
nsteps      = 5000000    ; 10 ns per window
nstxout-compressed = 500
nstenergy   = 500
nstlog      = 500
continuation = yes
gen-vel     = no
tcoupl      = V-rescale
tc-grps     = Protein non-Protein
tau-t       = 0.1     0.1
ref-t       = 300     300
pcoupl      = Parrinello-Rahman
pcoupltype  = semiisotropic
tau-p       = 5.0
ref-p       = 1.0 1.0
compressibility = 4.5e-5 4.5e-5
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
constraints = h-bonds
constraint-algorithm = LINCS

; Umbrella restraint
pull                 = yes
pull-ncoords         = 1
pull-ngroups         = 2
pull-group1-name     = gp120
pull-group2-name     = CD4
pull-coord1-type     = umbrella
pull-coord1-geometry = direction
pull-coord1-vec      = 0 0 1
pull-coord1-dim      = N N Y
pull-coord1-groups   = 1 2
pull-coord1-start    = yes
pull-coord1-rate     = 0        ; 0 = restrained, not pulled
pull-coord1-k        = 1000     ; kJ/mol/nm^2

pull-nstxout         = 50
pull-nstfout         = 50
""")

print("MDP files written: em.mdp, nvt.mdp, npt.mdp, pull.mdp, umbrella.mdp")

print()
print("=" * 60)
print("Step 6: Write launch script")
print("=" * 60)

with open("run_smd_v2.sh", "w") as f:
    f.write("""\
#!/bin/bash
set -e

NTOMP=12
cd "$(dirname "$0")"

echo "=== Energy Minimization ==="
gmx grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr -maxwarn 2
gmx mdrun -deffnm em -ntomp $NTOMP -pin on
echo "EM potential: $(grep 'Potential Energy' em.log | tail -1)"

echo "=== NVT Equilibration (100 ps) ==="
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2
gmx mdrun -deffnm nvt -ntomp $NTOMP -pin on

echo "=== NPT Equilibration (500 ps) ==="
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr -maxwarn 2
gmx mdrun -deffnm npt -ntomp $NTOMP -pin on

echo "=== Steered MD (5 ns, 0.001 nm/ps along Z) ==="
gmx grompp -f pull.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -n index.ndx -o pull.tpr -maxwarn 2
nohup gmx mdrun -s pull.tpr -deffnm pull -ntomp $NTOMP -pin on > pull_run.log 2>&1 &
echo "Pull PID: $!"
echo "Monitoring: tail -f pull_run.log"
""")

os.chmod("run_smd_v2.sh", 0o755)
print("Launch script: run_smd_v2.sh")

print()
print("=" * 60)
print("SETUP COMPLETE -- Review and run:")
print(f"  Box: {xy_box} x {xy_box} x {z_box} nm (rectangular)")
print(f"  Pull: 0.001 nm/ps along Z, 5 ns = 5 nm displacement")
print(f"  Windows: extract every 0.2 nm from pull trajectory")
print("=" * 60)
